"""
competitive_hybrid_model.py
---------------------------
Modelagem Híbrida Paralela Competitiva para Lixiviação de Calcina de Zinco.
Baseada na formulação teórica de Sansana et al. (2024):
  "Hybrid modeling for transfer learning in chemical processes"
  Chemical Engineering Science 300 (2024) 120568 (Seção 2.3.2, Eq. 10 e 11).

Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho

Arquitetura Paralela Competitiva:
  1. Base Física (White-Box, f_p):
     - Mantido estritamente inalterado: Balanço Populacional (PBM) com alpha = 0.
     - Predição pura baseada em primeiros princípios: X_wb(t) in [0, 1].
  2. Base Orientada a Dados (Black-Box, f_d):
     - Concorre diretamente na previsão da conversão de zinco: X_bb(t) in [0, 1].
     - Recebe as condições de processo: tempo_min, razao_molar_eta, C_acid_0_mol_L, razao_SL_g_L.
  3. Camada de Fusão / Meta Ponderada (Meta-Model):
     - Formulação de Sansana et al. (2024):
         X_raw = w * X_wb + (1 - w) * X_bb
     - Otimização do peso w in [0, 1] via mínimos quadrados restritos (solução analítica fechada e bounded).
  4. Acoplador Termodinâmico / Restrições Físicas:
     - Conservação de massa: clip(X, 0.0, 1.0)
     - Limite estequiométrico: X <= eta quando eta < 1.0
     - Monotonicidade temporal: dX/dt >= 0 (np.maximum.accumulate por ensaio).
  5. Validação Cruzada Estrita (LOGO-CV):
     - 16 folds (1 por ensaio de bancada).
     - Calibração de f_d e do peso w exclusivamente no conjunto de treino (15 ensaios).
     - Predição cega e cálculo de métricas (R², R²_adj, RMSE, MAE) no ensaio de teste.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Any, Tuple, List, Optional
from scipy.optimize import minimize_scalar

# Garantir saída UTF-8 no Windows
if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

from sklearn.ensemble import (
    GradientBoostingRegressor,
    RandomForestRegressor,
    ExtraTreesRegressor
)
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.pipeline import Pipeline
from sklearn.model_selection import LeaveOneGroupOut, GroupKFold, GroupShuffleSplit
from sklearn.base import clone

# Importações dos módulos locais
from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel


# =============================================================================
# 1. FUNÇÕES DE CÁLCULO DE MÉTRICAS E ESTATÍSTICA
# =============================================================================

def calcular_metricas_estatisticas(y_real: np.ndarray,
                                   y_pred: np.ndarray,
                                   p_parametros: int = 5,
                                   eta: Optional[np.ndarray] = None) -> Dict[str, float]:
    """
    Calcula métricas de desempenho estatístico e consistência física:
      - R² (Coeficiente de Determinação)
      - R²_ajustado = 1 - [(1 - R²)*(n - 1) / (n - p - 1)]
      - RMSE (Raiz do Erro Médio Quadrático)
      - MAE (Erro Médio Absoluto)
      - Violações de Limites Físicos (X < 0 ou X > 1)
      - Violações de Limites Estequiométricos (X > eta quando eta < 1.0)
    """
    y_real = np.asarray(y_real, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    n = len(y_real)
    
    residuos = y_real - y_pred
    ss_res = float(np.sum(residuos ** 2))
    ss_tot = float(np.sum((y_real - np.mean(y_real)) ** 2))
    
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    
    graus_liberdade = n - p_parametros - 1
    if graus_liberdade > 0 and ss_tot > 0:
        r2_ajustado = 1.0 - ((1.0 - r2) * (n - 1) / graus_liberdade)
    else:
        r2_ajustado = r2
        
    rmse = float(np.sqrt(np.mean(residuos ** 2)))
    mae = float(np.mean(np.abs(residuos)))
    
    # Verificação de violações
    viol_limites = np.sum((y_pred < -1e-5) | (y_pred > 1.0 + 1e-5))
    viol_pct = float(viol_limites / n * 100.0)
    
    viol_esteq_pct = 0.0
    if eta is not None:
        eta_arr = np.asarray(eta, dtype=float)
        mascara_deficiente = eta_arr < 1.0
        if np.any(mascara_deficiente):
            viol_esteq = np.sum(y_pred[mascara_deficiente] > (eta_arr[mascara_deficiente] + 1e-4))
            viol_esteq_pct = float(viol_esteq / np.sum(mascara_deficiente) * 100.0)
            
    return {
        "R2": r2,
        "R2_ajustado": r2_ajustado,
        "RMSE": rmse,
        "MAE": mae,
        "Viol_Limites_pct": viol_pct,
        "Viol_Estequiometria_pct": viol_esteq_pct,
        "Violacao_pct": max(viol_pct, viol_esteq_pct)
    }


def verificar_monotonicidade(df: pd.DataFrame, coluna_pred: str) -> float:
    """Calcula a porcentagem de passos temporais onde dX/dt < -1e-4."""
    violacoes = 0
    total_passos = 0
    for _, grp in df.groupby("ensaio_id", sort=False):
        x_vals = grp[coluna_pred].values
        diffs = np.diff(x_vals)
        violacoes += np.sum(diffs < -1e-4)
        total_passos += len(diffs)
    return float(violacoes / total_passos * 100.0) if total_passos > 0 else 0.0


# =============================================================================
# 2. ACOPLADOR TERMODINÂMICO
# =============================================================================

class AcopladorTermodinamico:
    """Aplica conservação de massa, estequiometria e monotonicidade temporal."""
    @staticmethod
    def acoplar(df: pd.DataFrame, x_raw: np.ndarray) -> np.ndarray:
        eta = df['razao_molar_eta'].values
        x_restringido = np.zeros_like(x_raw)
        
        for i in range(len(x_raw)):
            teto = min(1.0, float(eta[i])) if eta[i] < 1.0 else 1.0
            x_restringido[i] = np.clip(x_raw[i], 0.0, teto)
            
        df_temp = df.copy()
        df_temp['X_restringido'] = x_restringido
        
        x_final = []
        if 'ensaio_id' in df_temp.columns:
            for _, grp in df_temp.groupby('ensaio_id', sort=False):
                curva_monotonica = np.maximum.accumulate(grp['X_restringido'].values)
                x_final.extend(curva_monotonica)
        else:
            curva_monotonica = np.maximum.accumulate(df_temp['X_restringido'].values)
            x_final.extend(curva_monotonica)
            
        return np.array(x_final)


# =============================================================================
# 3. META-PONDERADOR COMPETITIVO (Sansana et al., 2024)
# =============================================================================

class MetaPonderadorCompetitivo:
    """
    Implementa a camada de meta-ponderação competitiva:
      f_H(x) = w * f_p(x) + (1 - w) * f_d(x)
    Otimiza w in [0, 1] via mínimos quadrados restritos.
    """
    def __init__(self, limite_inferior: float = 0.0, limite_superior: float = 1.0):
        self.lb = limite_inferior
        self.ub = limite_superior
        self.w_otimo = 0.5

    def ajustar(self, y_real: np.ndarray, y_wb: np.ndarray, y_bb: np.ndarray) -> float:
        """
        Resolve analiticamente o problema de mínimos quadrados restritos:
          min_{w in [0, 1]} sum (y_real - [w * y_wb + (1 - w) * y_bb])^2
        """
        y_real = np.asarray(y_real, dtype=float)
        y_wb = np.asarray(y_wb, dtype=float)
        y_bb = np.asarray(y_bb, dtype=float)
        
        # d_i = y_wb - y_bb
        # e_i = y_real - y_bb
        # Erro residual = e_i - w * d_i
        d = y_wb - y_bb
        e = y_real - y_bb
        
        denom = float(np.sum(d ** 2))
        if denom < 1e-12:
            self.w_otimo = 0.5
        else:
            w_analitico = float(np.sum(d * e) / denom)
            self.w_otimo = float(np.clip(w_analitico, self.lb, self.ub))
            
        return self.w_otimo

    def combinar(self, y_wb: np.ndarray, y_bb: np.ndarray, w: Optional[float] = None) -> np.ndarray:
        """Calcula a predição competitiva bruta: w * y_wb + (1 - w) * y_bb."""
        peso = self.w_otimo if w is None else w
        return peso * np.asarray(y_wb, dtype=float) + (1.0 - peso) * np.asarray(y_bb, dtype=float)


# =============================================================================
# 4. CATÁLOGO DE MODELOS BLACK-BOX PARA PREVISÃO DIRETA DE CONVERSÃO
# =============================================================================

def obter_catalogo_blackbox(random_state: int = 42) -> Dict[str, Any]:
    """
    Retorna o catálogo de algoritmos de Machine Learning concorrentes.
    Todos os modelos recebem x = (tempo_min, razao_molar_eta, C_acid_0_mol_L, razao_SL_g_L)
    e predizem diretamente a conversão de zinco y = X_zn_exp in [0, 1].
    """
    return {
        'Gradient Boosting (GBDT)': GradientBoostingRegressor(
            n_estimators=85,
            max_depth=3,
            learning_rate=0.06,
            random_state=random_state
        ),
        'Extra Trees': ExtraTreesRegressor(
            n_estimators=100,
            max_depth=6,
            min_samples_split=2,
            random_state=random_state
        ),
        'Random Forest': RandomForestRegressor(
            n_estimators=100,
            max_depth=5,
            min_samples_split=2,
            random_state=random_state
        ),
        'Polynomial Ridge': Pipeline([
            ('scaler', StandardScaler()),
            ('poly', PolynomialFeatures(degree=2, include_bias=False)),
            ('ridge', Ridge(alpha=1.0, random_state=random_state))
        ]),
        'SVR (RBF Otimizado)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(C=20.0, epsilon=0.01, gamma='scale'))
        ]),
        'MLP (Tanh + L-BFGS)': Pipeline([
            ('scaler', StandardScaler()),
            ('mlp', MLPRegressor(
                hidden_layer_sizes=(16, 8),
                activation='tanh',
                solver='lbfgs',
                alpha=0.01,
                max_iter=3000,
                random_state=random_state
            ))
        ])
    }


# =============================================================================
# 5. CLASSE PRINCIPAL: MODELO HÍBRIDO PARALELO COMPETITIVO
# =============================================================================

class ModeloHibridoParaleloCompetitivo:
    """
    Implementação completa da arquitetura Híbrida Paralela Competitiva:
      - Mantém o White-Box (PBM com alpha=0) estritamente intacto
      - Treina o Black-Box para conversão direta X_bb
      - Realiza meta-ponderação w * X_wb + (1 - w) * X_bb
      - Aplica acoplamento termodinâmico para 0% violações físicas
    """
    FEATURES_OPERACIONAIS = [
        'tempo_min',
        'razao_molar_eta',
        'C_acid_0_mol_L',
        'razao_SL_g_L'
    ]

    def __init__(self, nome_algoritmo: str = 'Gradient Boosting (GBDT)', random_state: int = 42):
        self.nome_algoritmo = nome_algoritmo
        self.random_state = random_state
        catalogo = obter_catalogo_blackbox(random_state)
        if nome_algoritmo not in catalogo:
            raise ValueError(f"Algoritmo '{nome_algoritmo}' não encontrado no catálogo: {list(catalogo.keys())}")
        self.regressor_bb = catalogo[nome_algoritmo]
        self.pbm = PopulationBalanceModel()
        self.meta_ponderador = MetaPonderadorCompetitivo()
        self.acoplador = AcopladorTermodinamico()
        self.peso_ajustado = 0.5
        self.is_treinado = False

    def simular_whitebox_puro(self, df: pd.DataFrame) -> np.ndarray:
        """Calcula a predição da física pura (PBM com alpha = 0.0) para todos os ensaios."""
        x_wb_lista = []
        if 'ensaio_id' in df.columns:
            for _, grp in df.groupby('ensaio_id', sort=False):
                ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
                eta = float(grp['razao_molar_eta'].iloc[0])
                tempos = grp['tempo_min'].tolist()
                res = self.pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
                x_wb_lista.extend(res['X_pbm'])
        else:
            ca0 = float(df['C_acid_0_mol_L'].iloc[0])
            eta = float(df['razao_molar_eta'].iloc[0])
            tempos = df['tempo_min'].tolist()
            res = self.pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
            x_wb_lista.extend(res['X_pbm'])
        return np.array(x_wb_lista, dtype=float)

    def preparar_features_e_alvo(self, df: pd.DataFrame) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        """Extrai a matriz X operacional, o alvo real y = X_zn_exp (se existir) e os grupos (ensaio_id)."""
        X = df[self.FEATURES_OPERACIONAIS].values
        y = df['X_zn_exp'].values if 'X_zn_exp' in df.columns else None
        grupos = df['ensaio_id'].values if 'ensaio_id' in df.columns else None
        return X, y, grupos

    def treinar(self, df: pd.DataFrame) -> float:
        """
        Treina o Black-Box em 100% dos dados para predição direta de conversão,
        simula o White-Box puro e calibra o peso ótimo w in [0, 1].
        """
        X, y, _ = self.preparar_features_e_alvo(df)
        self.regressor_bb.fit(X, y)
        
        # Simula White-Box puro
        y_wb = self.simular_whitebox_puro(df)
        y_bb = self.regressor_bb.predict(X)
        
        # Calibra o peso ótimo w via Mínimos Quadrados Restritos
        self.peso_ajustado = self.meta_ponderador.ajustar(y, y_wb, y_bb)
        self.is_treinado = True
        return self.peso_ajustado

    def prever(self, df: pd.DataFrame, aplicar_acoplador: bool = True) -> np.ndarray:
        """
        Gera a predição competitiva híbrida combinada:
          X_raw = w * X_wb + (1 - w) * X_bb
        e aplica o acoplador físico.
        """
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de realizar predições.")
        
        X, _, _ = self.preparar_features_e_alvo(df)
        y_wb = self.simular_whitebox_puro(df)
        y_bb = self.regressor_bb.predict(X)
        
        x_raw = self.meta_ponderador.combinar(y_wb, y_bb, self.peso_ajustado)
        
        if aplicar_acoplador:
            return self.acoplador.acoplar(df, x_raw)
        return x_raw

    def validacao_cruzada_logo(self, df: pd.DataFrame, aplicar_acoplador: bool = True) -> Tuple[np.ndarray, Dict[int, float]]:
        """
        Executa validação cruzada Leave-One-Group-Out (LOGO-CV) em 16 folds:
          - Em cada fold k, o regressor Black-Box é treinado nos 15 ensaios restantes.
          - O peso w é calibrado exclusivamente nos 15 ensaios de treino.
          - Prediz-se a conversão no ensaio k cego: w * X_wb + (1 - w) * X_bb.
        Garante ausência total de vazamento de dados (data leakage).
        """
        logo = LeaveOneGroupOut()
        X, y, grupos = self.preparar_features_e_alvo(df)
        y_wb = self.simular_whitebox_puro(df)
        
        y_pred_cv = np.zeros_like(y)
        pesos_cv = {}
        
        for idx_treino, idx_teste in logo.split(X, y, groups=grupos):
            ens_teste_id = int(grupos[idx_teste][0])
            X_tr, y_tr = X[idx_treino], y[idx_treino]
            X_te = X[idx_teste]
            y_wb_tr, y_wb_te = y_wb[idx_treino], y_wb[idx_teste]
            
            # Treina clone do regressor no fold de treino
            reg_fold = clone(self.regressor_bb)
            reg_fold.fit(X_tr, y_tr)
            y_bb_tr = reg_fold.predict(X_tr)
            y_bb_te = reg_fold.predict(X_te)
            
            # Calibra peso w exclusivamente no treino
            w_fold = self.meta_ponderador.ajustar(y_tr, y_wb_tr, y_bb_tr)
            pesos_cv[ens_teste_id] = w_fold
            
            # Predição combinada bruta no ensaio cego de teste
            y_pred_cv[idx_teste] = self.meta_ponderador.combinar(y_wb_te, y_bb_te, w_fold)
            
        if aplicar_acoplador:
            y_pred_cv = self.acoplador.acoplar(df, y_pred_cv)
            
        return y_pred_cv, pesos_cv

    def validacao_cruzada_85_15(self, df: pd.DataFrame, n_splits: int = 6, aplicar_acoplador: bool = True) -> Tuple[np.ndarray, List[float]]:
        """
        Executa validação cruzada agrupada treinando com ~85% dos dados e validando com ~15% (GroupKFold):
          - Em 16 ensaios, k=6 folds particiona em ~13 a 14 ensaios de treino (81.25% a 87.5%)
            e 2 a 3 ensaios de teste cego por fold (12.5% a 18.75%, média exata de ~15%).
          - O Black-Box é treinado nos 85% e o peso w é calibrado nos 85%.
          - Prediz os 15% de ensaios mantidos fora sem qualquer vazamento de dados.
          - Avalia todos os 128 pontos do conjunto de dados de forma estritamente out-of-fold.
        """
        gkf = GroupKFold(n_splits=n_splits)
        X, y, grupos = self.preparar_features_e_alvo(df)
        y_wb = self.simular_whitebox_puro(df)
        
        y_pred_cv = np.zeros_like(y)
        pesos_cv = []
        
        for idx_treino, idx_teste in gkf.split(X, y, groups=grupos):
            X_tr, y_tr = X[idx_treino], y[idx_treino]
            X_te = X[idx_teste]
            y_wb_tr, y_wb_te = y_wb[idx_treino], y_wb[idx_teste]
            
            reg_fold = clone(self.regressor_bb)
            reg_fold.fit(X_tr, y_tr)
            y_bb_tr = reg_fold.predict(X_tr)
            y_bb_te = reg_fold.predict(X_te)
            
            w_fold = self.meta_ponderador.ajustar(y_tr, y_wb_tr, y_bb_tr)
            pesos_cv.append(w_fold)
            
            y_pred_cv[idx_teste] = self.meta_ponderador.combinar(y_wb_te, y_bb_te, w_fold)
            
        if aplicar_acoplador:
            y_pred_cv = self.acoplador.acoplar(df, y_pred_cv)
            
        return y_pred_cv, pesos_cv

    def validacao_holdout_85_15(self, df: pd.DataFrame, test_size: float = 0.15, random_state: int = 42) -> Dict[str, Any]:
        """
        Executa uma divisão de Holdout Agrupado: ~85% de ensaios para treino (13-14) vs ~15% para teste (2-3).
        """
        gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
        X, y, grupos = self.preparar_features_e_alvo(df)
        y_wb = self.simular_whitebox_puro(df)
        
        tr_idx, te_idx = next(gss.split(X, y, groups=grupos))
        
        X_tr, y_tr = X[tr_idx], y[tr_idx]
        X_te, y_te = X[te_idx], y[te_idx]
        y_wb_tr, y_wb_te = y_wb[tr_idx], y_wb[te_idx]
        
        reg = clone(self.regressor_bb)
        reg.fit(X_tr, y_tr)
        y_bb_tr = reg.predict(X_tr)
        y_bb_te = reg.predict(X_te)
        
        w = self.meta_ponderador.ajustar(y_tr, y_wb_tr, y_bb_tr)
        
        pred_tr = self.meta_ponderador.combinar(y_wb_tr, y_bb_tr, w)
        pred_te = self.meta_ponderador.combinar(y_wb_te, y_bb_te, w)
        
        df_tr = df.iloc[tr_idx].copy()
        df_te = df.iloc[te_idx].copy()
        
        pred_tr_acop = self.acoplador.acoplar(df_tr, pred_tr)
        pred_te_acop = self.acoplador.acoplar(df_te, pred_te)
        
        m_tr = calcular_metricas_estatisticas(y_tr, pred_tr_acop, p_parametros=5, eta=df_tr['razao_molar_eta'].values)
        m_te = calcular_metricas_estatisticas(y_te, pred_te_acop, p_parametros=5, eta=df_te['razao_molar_eta'].values)
        
        ens_tr = df_tr['ensaio_id'].unique().tolist()
        ens_te = df_te['ensaio_id'].unique().tolist()
        
        return {
            'ensaios_treino': ens_tr,
            'ensaios_teste': ens_te,
            'pct_treino': len(ens_tr) / len(df['ensaio_id'].unique()) * 100.0,
            'pct_teste': len(ens_te) / len(df['ensaio_id'].unique()) * 100.0,
            'peso_w': w,
            'metricas_treino': m_tr,
            'metricas_teste': m_te,
            'y_treino_real': y_tr,
            'y_treino_pred': pred_tr_acop,
            'y_teste_real': y_te,
            'y_teste_pred': pred_te_acop
        }


# =============================================================================
# 6. PIPELINES DE VALIDAÇÃO (85% TREINO / 15% VALIDAÇÃO E LOGO-CV)
# =============================================================================

def avaliar_catalogo_competitivo_85_15(df: pd.DataFrame,
                                       p_parametros: int = 5,
                                       n_splits: int = 6,
                                       random_state: int = 42) -> pd.DataFrame:
    """
    Executa a validação cruzada agrupada treinando com ~85% e validando com ~15% (GroupKFold k=6)
    para todo o catálogo de regressores Black-Box na arquitetura competitiva.
    """
    catalogo = obter_catalogo_blackbox(random_state)
    gkf = GroupKFold(n_splits=n_splits)
    
    pbm = PopulationBalanceModel()
    y_wb_completo = []
    for _, grp in df.groupby('ensaio_id', sort=False):
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        tempos = grp['tempo_min'].tolist()
        res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        y_wb_completo.extend(res['X_pbm'])
    y_wb_arr = np.array(y_wb_completo, dtype=float)
    df_base = df.copy()
    df_base['X_pbm_puro'] = y_wb_arr
    
    X = df_base[ModeloHibridoParaleloCompetitivo.FEATURES_OPERACIONAIS].values
    y = df_base['X_zn_exp'].values
    grupos = df_base['ensaio_id'].values
    eta = df_base['razao_molar_eta'].values
    
    resultados = []
    ponderador = MetaPonderadorCompetitivo()
    acoplador = AcopladorTermodinamico()
    
    for nome_modelo, regressor_base in catalogo.items():
        # Ajuste em 100% dos dados
        reg_fit = clone(regressor_base)
        reg_fit.fit(X, y)
        y_bb_fit = reg_fit.predict(X)
        w_fit = ponderador.ajustar(y, y_wb_arr, y_bb_fit)
        x_raw_fit = ponderador.combinar(y_wb_arr, y_bb_fit, w_fit)
        x_pred_fit = acoplador.acoplar(df_base, x_raw_fit)
        m_fit = calcular_metricas_estatisticas(y, x_pred_fit, p_parametros=p_parametros, eta=eta)
        
        # Validação 85% treino / 15% validação (GroupKFold k=6)
        y_pred_cv = np.zeros_like(y)
        y_bb_cv = np.zeros_like(y)
        pesos_cv = []
        
        for idx_treino, idx_teste in gkf.split(X, y, groups=grupos):
            X_tr, y_tr = X[idx_treino], y[idx_treino]
            X_te = X[idx_teste]
            y_wb_tr, y_wb_te = y_wb_arr[idx_treino], y_wb_arr[idx_teste]
            
            reg_cv = clone(regressor_base)
            reg_cv.fit(X_tr, y_tr)
            y_bb_tr = reg_cv.predict(X_tr)
            y_bb_te = reg_cv.predict(X_te)
            y_bb_cv[idx_teste] = y_bb_te
            
            w_fold = ponderador.ajustar(y_tr, y_wb_tr, y_bb_tr)
            pesos_cv.append(w_fold)
            
            y_pred_cv[idx_teste] = ponderador.combinar(y_wb_te, y_bb_te, w_fold)
            
        x_pred_cv_acoplado = acoplador.acoplar(df_base, y_pred_cv)
        m_cv = calcular_metricas_estatisticas(y, x_pred_cv_acoplado, p_parametros=p_parametros, eta=eta)
        df_temp_cv = df_base.copy()
        df_temp_cv['pred_cv'] = x_pred_cv_acoplado
        mono_cv = verificar_monotonicidade(df_temp_cv, 'pred_cv')
        
        x_bb_acoplado = acoplador.acoplar(df_base, y_bb_cv)
        m_bb_puro = calcular_metricas_estatisticas(y, x_bb_acoplado, p_parametros=4, eta=eta)
        
        resultados.append({
            'Modelo Black-Box': nome_modelo,
            'w_Ajuste': w_fit,
            'w_medio_85_15': float(np.mean(pesos_cv)),
            'w_min_85_15': float(np.min(pesos_cv)),
            'w_max_85_15': float(np.max(pesos_cv)),
            'R2_Ajuste': m_fit['R2'],
            'R2_adj_Ajuste': m_fit['R2_ajustado'],
            'RMSE_Ajuste': m_fit['RMSE'],
            'R2_Val_85_15': m_cv['R2'],
            'R2_adj_Val_85_15': m_cv['R2_ajustado'],
            'RMSE_Val_85_15': m_cv['RMSE'],
            'MAE_Val_85_15': m_cv['MAE'],
            'Viol_Fisica_pct': m_cv['Viol_Limites_pct'],
            'Viol_Mono_pct': mono_cv,
            'R2_Val_BB_Isolado': m_bb_puro['R2'],
            'RMSE_Val_BB_Isolado': m_bb_puro['RMSE']
        })
        
    df_res = pd.DataFrame(resultados).sort_values(by='R2_Val_85_15', ascending=False).reset_index(drop=True)
    return df_res

def avaliar_catalogo_competitivo_logo_cv(df: pd.DataFrame,
                                        p_parametros: int = 5,
                                        random_state: int = 42) -> pd.DataFrame:
    """
    Executa a validação cruzada Leave-One-Group-Out (LOGO-CV) com 16 folds
    para cada modelo do catálogo Black-Box na arquitetura Paralela Competitiva.
    Garante calibração isenta de data leakage para o regressor e para o peso w.
    """
    catalogo = obter_catalogo_blackbox(random_state)
    logo = LeaveOneGroupOut()
    
    # 1. Simulação prévia do White-Box puro (independente de treino/teste, pois é analítico da física pura)
    pbm = PopulationBalanceModel()
    y_wb_completo = []
    for _, grp in df.groupby('ensaio_id', sort=False):
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        tempos = grp['tempo_min'].tolist()
        res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        y_wb_completo.extend(res['X_pbm'])
    y_wb_arr = np.array(y_wb_completo, dtype=float)
    df_base = df.copy()
    df_base['X_pbm_puro'] = y_wb_arr
    
    X = df_base[ModeloHibridoParaleloCompetitivo.FEATURES_OPERACIONAIS].values
    y = df_base['X_zn_exp'].values
    grupos = df_base['ensaio_id'].values
    eta = df_base['razao_molar_eta'].values
    
    resultados = []
    ponderador = MetaPonderadorCompetitivo()
    acoplador = AcopladorTermodinamico()
    
    for nome_modelo, regressor_base in catalogo.items():
        # --- A. Ajuste em 100% dos dados (Ajuste Global) ---
        reg_fit = clone(regressor_base)
        reg_fit.fit(X, y)
        y_bb_fit = reg_fit.predict(X)
        w_fit = ponderador.ajustar(y, y_wb_arr, y_bb_fit)
        x_raw_fit = ponderador.combinar(y_wb_arr, y_bb_fit, w_fit)
        x_pred_fit = acoplador.acoplar(df_base, x_raw_fit)
        
        m_fit = calcular_metricas_estatisticas(y, x_pred_fit, p_parametros=p_parametros, eta=eta)
        df_temp_fit = df_base.copy()
        df_temp_fit['pred'] = x_pred_fit
        mono_fit = verificar_monotonicidade(df_temp_fit, 'pred')
        
        # --- B. Validação Cruzada Cega (LOGO-CV) ---
        y_pred_cv = np.zeros_like(y)
        y_bb_cv = np.zeros_like(y)
        pesos_cv = []
        
        for idx_treino, idx_teste in logo.split(X, y, groups=grupos):
            X_tr, y_tr = X[idx_treino], y[idx_treino]
            X_te = X[idx_teste]
            y_wb_tr, y_wb_te = y_wb_arr[idx_treino], y_wb_arr[idx_teste]
            
            # Treina regressor nos 15 ensaios
            reg_cv = clone(regressor_base)
            reg_cv.fit(X_tr, y_tr)
            y_bb_tr = reg_cv.predict(X_tr)
            y_bb_te = reg_cv.predict(X_te)
            y_bb_cv[idx_teste] = y_bb_te
            
            # Calibra peso w exclusivamente nos 15 ensaios de treino
            w_fold = ponderador.ajustar(y_tr, y_wb_tr, y_bb_tr)
            pesos_cv.append(w_fold)
            
            # Predição no ensaio de teste
            x_raw_te = ponderador.combinar(y_wb_te, y_bb_te, w_fold)
            y_pred_cv[idx_teste] = x_raw_te
            
        # Aplicação do acoplador físico nas predições agregadas de teste
        x_pred_cv_acoplado = acoplador.acoplar(df_base, y_pred_cv)
        
        m_cv = calcular_metricas_estatisticas(y, x_pred_cv_acoplado, p_parametros=p_parametros, eta=eta)
        df_temp_cv = df_base.copy()
        df_temp_cv['pred_cv'] = x_pred_cv_acoplado
        mono_cv = verificar_monotonicidade(df_temp_cv, 'pred_cv')
        
        # Desempenho do Black-Box puro (isolado) em CV
        x_bb_acoplado = acoplador.acoplar(df_base, y_bb_cv)
        m_bb_puro = calcular_metricas_estatisticas(y, x_bb_acoplado, p_parametros=4, eta=eta)
        
        resultados.append({
            'Modelo Black-Box': nome_modelo,
            'w_Ajuste': w_fit,
            'w_medio_CV': float(np.mean(pesos_cv)),
            'w_min_CV': float(np.min(pesos_cv)),
            'w_max_CV': float(np.max(pesos_cv)),
            'R2_Ajuste': m_fit['R2'],
            'R2_adj_Ajuste': m_fit['R2_ajustado'],
            'RMSE_Ajuste': m_fit['RMSE'],
            'MAE_Ajuste': m_fit['MAE'],
            'R2_CV': m_cv['R2'],
            'R2_adj_CV': m_cv['R2_ajustado'],
            'RMSE_CV': m_cv['RMSE'],
            'MAE_CV': m_cv['MAE'],
            'Viol_Fisica_pct': m_cv['Viol_Limites_pct'],
            'Viol_Mono_pct': mono_cv,
            'R2_CV_BB_Isolado': m_bb_puro['R2'],
            'RMSE_CV_BB_Isolado': m_bb_puro['RMSE']
        })
        
    df_res = pd.DataFrame(resultados).sort_values(by='R2_CV', ascending=False).reset_index(drop=True)
    return df_res


# =============================================================================
# 7. EXECUÇÃO PRINCIPAL
# =============================================================================

def executar_experimento_competitivo():
    print("=" * 110)
    print("MODELAGEM HÍBRIDA PARALELA COMPETITIVA (Sansana et al., 2024)")
    print("Laboratório de Operações e Processos -- Engenharia Química -- UFMG (2026)")
    print("=" * 110)
    
    df = carregar_dados_bancada()
    print(f"-> Dados carregados: {len(df)} pontos experimentais distribuídos em {df['ensaio_id'].nunique()} ensaios.")
    
    print("\n[Etapa 1] Avaliando Catálogo Black-Box na Arquitetura Competitiva via 85% Treino / 15% Validação (GroupKFold k=6)...")
    df_benchmark_85_15 = avaliar_catalogo_competitivo_85_15(df, n_splits=6)
    
    print("\n" + "=" * 115)
    print("RESULTADOS DO BENCHMARK COMPETITIVO: 85% TREINO / 15% VALIDAÇÃO (GroupKFold k=6):")
    print("=" * 115)
    colunas_85_15 = [
        'Modelo Black-Box', 'w_Ajuste', 'w_medio_85_15',
        'R2_Ajuste', 'R2_Val_85_15', 'R2_adj_Val_85_15', 'RMSE_Val_85_15', 'MAE_Val_85_15', 'Viol_Mono_pct'
    ]
    print(df_benchmark_85_15[colunas_85_15].to_string(index=False, justify='center', float_format='{:,.4f}'.format))
    
    melhor_85_15 = df_benchmark_85_15.iloc[0]
    print("\n" + "-" * 115)
    print(f"-> Campeão 85/15: {melhor_85_15['Modelo Black-Box']}")
    print(f"   * Fator de Ponderação Ótimo (w_ajuste): {melhor_85_15['w_Ajuste']:.4f} (White-Box) vs {1.0 - melhor_85_15['w_Ajuste']:.4f} (Black-Box)")
    print(f"   * Peso médio em validação 85/15: {melhor_85_15['w_medio_85_15']:.4f} [min: {melhor_85_15['w_min_85_15']:.4f}, max: {melhor_85_15['w_max_85_15']:.4f}]")
    print(f"   * R² Validação 85/15: {melhor_85_15['R2_Val_85_15']:.4f} | R² Ajustado: {melhor_85_15['R2_adj_Val_85_15']:.4f}")
    print(f"   * RMSE Validação 85/15: {melhor_85_15['RMSE_Val_85_15']:.4f} ({melhor_85_15['RMSE_Val_85_15']*100:.2f}%) | MAE: {melhor_85_15['MAE_Val_85_15']:.4f}")
    print(f"   * R² do Black-Box Isolado: {melhor_85_15['R2_Val_BB_Isolado']:.4f} -> Ganho pela fusão híbrida: +{melhor_85_15['R2_Val_85_15'] - melhor_85_15['R2_Val_BB_Isolado']:.4f}")
    print(f"   * Violações Termodinâmicas e Monotonicidade: 0.0% (100% aderente)")
    print("-" * 115)
    
    print("\n[Etapa 2] Avaliando Holdout Estruturado 85% Treino / 15% Validação...")
    mod_gbdt = ModeloHibridoParaleloCompetitivo('Gradient Boosting (GBDT)')
    res_holdout = mod_gbdt.validacao_holdout_85_15(df, test_size=0.15, random_state=42)
    print(f"   * Ensaios Treino ({res_holdout['pct_treino']:.1f}%): {res_holdout['ensaios_treino']}")
    print(f"   * Ensaios Teste Cego ({res_holdout['pct_teste']:.1f}%): {res_holdout['ensaios_teste']}")
    print(f"   * Peso w calibrado no Treino: {res_holdout['peso_w']:.4f}")
    print(f"   * Treino 85%: R² = {res_holdout['metricas_treino']['R2']:.4f}, RMSE = {res_holdout['metricas_treino']['RMSE']:.4f}")
    print(f"   * Teste Cego 15%: R² = {res_holdout['metricas_teste']['R2']:.4f}, RMSE = {res_holdout['metricas_teste']['RMSE']:.4f}, MAE = {res_holdout['metricas_teste']['MAE']:.4f}")
    
    print("\n[Etapa 3] Comparativo com Split 15/1 (LOGO-CV)...")
    df_benchmark_logo = avaliar_catalogo_competitivo_logo_cv(df)
    melhor_logo = df_benchmark_logo.iloc[0]
    
    print("\n" + "=" * 115)
    print("CONFRONTO DE CONSEQUÊNCIAS: VALIDAÇÃO 85/15 (GroupKFold) vs VALIDAÇÃO 15/1 (LOGO-CV):")
    print("=" * 115)
    print(f"{'Estratégia de Validação':<35} | {'Treino':<12} | {'Validação':<12} | {'R² Validação':<14} | {'RMSE Validação':<14} | {'w Médio':<8}")
    print("-" * 115)
    print(f"{'Validação 85/15 (GroupKFold k=6)':<35} | {'~13 ensaios':<12} | {'~3 ensaios':<12} | "
          f"{melhor_85_15['R2_Val_85_15']:^14.4f} | {melhor_85_15['RMSE_Val_85_15']:^14.4f} | {melhor_85_15['w_medio_85_15']:^8.4f}")
    print(f"{'Validação 15/1 (LOGO-CV 16 folds)':<35} | {'15 ensaios':<12} | {'1 ensaio':<12} | "
          f"{melhor_logo['R2_CV']:^14.4f} | {melhor_logo['RMSE_CV']:^14.4f} | {melhor_logo['w_medio_CV']:^8.4f}")
    print("=" * 115)
    
    # Gerar Gráfico Científico de Paridade (300 DPI)
    print("\n[Etapa 4] Gerando Gráfico Científico de Paridade e Resíduos (300 DPI)...")
    fig_paridade = gerar_grafico_paridade_competitivo(df, mod_gbdt, pasta_saida='scripts/docs')
    
    # Salvar resultados
    os.makedirs('scripts/docs', exist_ok=True)
    caminho_json = 'scripts/docs/benchmark_competitivo_85_15.json'
    df_benchmark_85_15.to_json(caminho_json, orient='records', indent=2)
    print(f"-> Relatório 85/15 exportado para '{caminho_json}'.")
    return df_benchmark_85_15


def gerar_grafico_paridade_competitivo(df: pd.DataFrame,
                                       mod_competitivo: Optional[ModeloHibridoParaleloCompetitivo] = None,
                                       pasta_saida: str = "scripts/docs") -> str:
    """
    Gera gráfico científico de paridade e análise residual em alta resolução (300 DPI)
    para a Modelagem Híbrida Paralela Competitiva (Sansana et al., 2024).
    Painel duplo:
      - Painel A: Paridade com categorização por regime de acidez eta (0.5, 1.0, 1.5, 3.1) e bandas de erro.
      - Painel B: Gráfico de resíduos (X_exp - X_prev) vs X_exp comprovando homocedasticidade e ausência de viés.
    """
    os.makedirs(pasta_saida, exist_ok=True)
    if mod_competitivo is None:
        mod_competitivo = ModeloHibridoParaleloCompetitivo('Gradient Boosting (GBDT)')
        mod_competitivo.treinar(df)
        
    y_exp = df['X_zn_exp'].values
    y_cv, pesos_cv = mod_competitivo.validacao_cruzada_85_15(df, n_splits=6)
    residuos = y_exp - y_cv
    
    m_cv = calcular_metricas_estatisticas(y_exp, y_cv, p_parametros=5, eta=df['razao_molar_eta'].values)
    w_med = float(np.mean(pesos_cv))
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # Configurações de cores e marcadores por regime de eta
    regimes = [
        (0.5, '#d62728', 'o', r'$\eta = 0{,}5$ (Deficiência Ácida)'),
        (1.0, '#ff7f0e', 's', r'$\eta = 1{,}0$ (Estequiométrico)'),
        (1.5, '#2ca02c', '^', r'$\eta = 1{,}5$ (Excesso Moderado)'),
        (3.1, '#1f77b4', 'D', r'$\eta = 3{,}1$ (Alto Excesso)')
    ]
    
    # -------------------------------------------------------------------------
    # PAINEL 1: GRÁFICO DE PARIDADE (X_exp vs X_prev)
    # -------------------------------------------------------------------------
    ax1 = axes[0]
    ax1.plot([0, 1.05], [0, 1.05], 'k--', lw=1.6, label='Ideal ($y=x$)', zorder=2)
    ax1.plot([0, 1.0], [0.05, 1.05], color='#555555', lw=1.0, ls=':', label=r'Margem $\pm 5\%$', zorder=2)
    ax1.plot([0.05, 1.05], [0, 1.0], color='#555555', lw=1.0, ls=':', zorder=2)
    ax1.fill_between([0, 1.05], [0.05, 1.10], [-0.05, 1.0], color='gray', alpha=0.08, label=r'Faixa $\pm 5\%$')
    
    for eta_val, cor, marker, label_reg in regimes:
        mask = np.isclose(df['razao_molar_eta'].values, eta_val, atol=0.05)
        ax1.scatter(
            y_exp[mask], y_cv[mask],
            color=cor, marker=marker, s=55, alpha=0.85,
            edgecolors='black', lw=0.6, label=label_reg, zorder=4
        )
        
    texto_metricas = (
        f"Validação 85/15 (GroupKFold)\n"
        f"$R^2 = {m_cv['R2']:.4f}$\n"
        f"$R^2_{{ajustado}} = {m_cv['R2_ajustado']:.4f}$\n"
        f"$RMSE = {m_cv['RMSE']*100:.2f}\\%$\n"
        f"$MAE = {m_cv['MAE']*100:.2f}\\%$\n"
        f"Peso White-Box $w = {w_med:.4f}$\n"
        f"Violações Físicas: 0.0%"
    )
    ax1.text(
        0.04, 0.96, texto_metricas,
        transform=ax1.transAxes, fontsize=9.0, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', edgecolor='#cccccc', alpha=0.92)
    )
    
    ax1.set_title("A) Gráfico de Paridade: Modelo Híbrido Paralelo Competitivo\n"
                  r"Formulação Sansana et al. (2024): $X_H = w X_{wb} + (1-w) X_{bb}$",
                  fontsize=10.5, fontweight='bold')
    ax1.set_xlabel('Conversão Experimental de Zinco ($X_{Zn}^{exp}$)', fontsize=10)
    ax1.set_ylabel('Conversão Prevista pelo Modelo Híbrido ($X_{Zn}^{prev}$)', fontsize=10)
    ax1.set_xlim(-0.02, 1.05)
    ax1.set_ylim(-0.02, 1.05)
    ax1.grid(True, linestyle='--', alpha=0.35)
    ax1.legend(loc='lower right', fontsize=8.0, frameon=True, framealpha=0.92)

    # -------------------------------------------------------------------------
    # PAINEL 2: ANÁLISE DE RESÍDUOS (X_exp - X_prev) vs X_exp
    # -------------------------------------------------------------------------
    ax2 = axes[1]
    ax2.axhline(0.0, color='black', lw=1.4, linestyle='-', zorder=2)
    ax2.axhline(0.03, color='#e65100', lw=1.0, linestyle='--', label=r'Tolerância $\pm 3\%$', zorder=2)
    ax2.axhline(-0.03, color='#e65100', lw=1.0, linestyle='--', zorder=2)
    ax2.axhline(0.05, color='#b71c1c', lw=0.9, linestyle=':', label=r'Tolerância $\pm 5\%$', zorder=2)
    ax2.axhline(-0.05, color='#b71c1c', lw=0.9, linestyle=':', zorder=2)
    ax2.fill_between([-0.02, 1.05], 0.03, -0.03, color='#2ca02c', alpha=0.10, label=r'Faixa de Alta Precisão ($\pm 3\%$)')
    
    for eta_val, cor, marker, label_reg in regimes:
        mask = np.isclose(df['razao_molar_eta'].values, eta_val, atol=0.05)
        ax2.scatter(
            y_exp[mask], residuos[mask],
            color=cor, marker=marker, s=55, alpha=0.85,
            edgecolors='black', lw=0.6, label=label_reg, zorder=4
        )
        
    texto_residuos = (
        f"Distribuição dos Resíduos:\n"
        f"Resíduo Médio: {np.mean(residuos):+.4f}\n"
        f"Desvio Padrão: {np.std(residuos):.4f}\n"
        f"Resíduo Máximo: {np.max(np.abs(residuos)):.4f}\n"
        f"Pontos dentro de $\\pm 3\\%$: {np.sum(np.abs(residuos) <= 0.03) / len(residuos)*100:.1f}%"
    )
    ax2.text(
        0.04, 0.96, texto_residuos,
        transform=ax2.transAxes, fontsize=9.0, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f9fa', edgecolor='#cccccc', alpha=0.92)
    )
    
    ax2.set_title("B) Análise de Resíduos Experimentais vs Conversão\n"
                  r"Homocedasticidade: $\Delta X = X_{Zn}^{exp} - X_{Zn}^{prev}$",
                  fontsize=10.5, fontweight='bold')
    ax2.set_xlabel('Conversão Experimental de Zinco ($X_{Zn}^{exp}$)', fontsize=10)
    ax2.set_ylabel(r'Erro Residual Instantâneo ($\Delta X$)', fontsize=10)
    ax2.set_xlim(-0.02, 1.05)
    ax2.set_ylim(-0.08, 0.08)
    ax2.grid(True, linestyle='--', alpha=0.35)
    ax2.legend(loc='lower right', fontsize=8.0, frameon=True, framealpha=0.92)

    plt.tight_layout()
    caminho_fig = os.path.join(pasta_saida, "paridade_modelo_competitivo.png")
    plt.savefig(caminho_fig)
    plt.close()
    print(f"-> Gráfico de paridade do modelo competitivo salvo em: {caminho_fig}")
    return caminho_fig


if __name__ == '__main__':
    executar_experimento_competitivo()
