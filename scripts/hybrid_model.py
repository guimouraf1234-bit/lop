"""
hybrid_model.py
---------------
Modelo Híbrido Cinza (Grey-Box) para Lixiviação de Calcina de Zinco em Reator Batelada.
Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho

Estrutura da Modelagem Híbrida:
  1. Base Física (White-Box): PBM via Método das Características com distribuição RRB.
  2. Machine Learning Residual (Black-Box): Aprende o desvio Delta_X = X_exp - X_pbm.
     Catálogo com múltiplos regressores otimizados (Gradient Boosting, Extra Trees,
     Random Forest, Gaussian Process / Kriging, SVR Otimizado, MLP Tanh + L-BFGS).
  3. Acoplador Termodinâmico: Garante dX/dt >= 0, conservação de massa e limites estequiométricos.
  4. Validação Cruzada Estrita (LOGO-CV): Avalia a generalização cega (Leave-One-Group-Out).
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Any, Tuple, List

# Algoritmos de Machine Learning e Validação
from sklearn.ensemble import (
    RandomForestRegressor,
    ExtraTreesRegressor,
    GradientBoostingRegressor
)
from sklearn.neural_network import MLPRegressor
from sklearn.svm import SVR
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import LeaveOneGroupOut, GroupShuffleSplit, GroupKFold
from sklearn.base import clone

# Importações dos módulos locais do projeto
from data_loader import carregar_dados_bancada, carregar_dados_com_residuo
from population_balance import PopulationBalanceModel


# Parâmetros grey-box (6 atributos termodinâmicos e operacionais)
FEATURES_FISICAS = [
    'tempo_min',           # Tempo de reação (cinética temporal)
    'razao_molar_eta',     # Razão estequiométrica ácido/zincita
    'C_acid_0_mol_L',      # Concentração inicial de H2SO4
    'razao_SL_g_L',        # Razão sólido-líquido
    'X_pbm',               # Conversão teórica do PBM
    'C_acid_pbm_mol_L'     # Concentração teórica de ácido remanescente
]


def preparar_features_e_alvo(df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extrai a matriz de atributos (X), o alvo residual (y = delta_X) 
    e os identificadores de ensaio (grupos) para validação cruzada.
    """
    if 'C_acid_pbm_mol_L' not in df.columns:
        df['C_acid_pbm_mol_L'] = df['C_acid_0_mol_L'] * np.maximum(
            0.0, 1.0 - (df['X_pbm'] / df['razao_molar_eta'])
        )

    y = None
    if 'delta_X' in df.columns:
        y = df['delta_X'].values
    elif 'X_zn_exp' in df.columns and 'X_pbm' in df.columns:
        df['delta_X'] = df['X_zn_exp'] - df['X_pbm']
        y = df['delta_X'].values

    X = df[FEATURES_FISICAS].values
    grupos = df['ensaio_id'].values if 'ensaio_id' in df.columns else None

    return X, y, grupos


# Definição e configuração dos modelos de machine learning
def obter_catalogo_modelos(random_state: int = 42) -> Dict[str, Any]:
    """
    Retorna o catálogo de regressores residuais avaliados no TCC:
      - Gradient Boosting (GBDT): Campeão geral com minimização sequencial de resíduos (R²_CV = 0.9934).
      - Extra Trees: Árvores extremamente aleatorizadas com excelente suavização de variância (R²_CV = 0.9932).
      - Random Forest: Robusto contra ruídos e não-linearidades (R²_CV = 0.9907).
      - Gaussian Process (Kriging): Processo Gaussiano com kernel Matérn e ruído WhiteKernel (R²_CV = 0.9878).
      - SVR (Otimizado): Kernel RBF com regularização balanceada C=5.0 e eps=0.005 (R²_CV = 0.9870).
      - Rede Neural (MLP Tanh + L-BFGS): Ativação suave tanh e Quasi-Newton de segunda ordem (R²_CV = 0.9688).
      - SVR (Baseline): Kernel RBF com C=1.0 e eps=0.01 (R²_CV = 0.9833).
      - Rede Neural (MLP ReLU + Adam): Baseline linear por partes com Adam estocástico (R²_CV = 0.8013).
    """
    return {
        'Gradient Boosting (GBDT)': GradientBoostingRegressor(
            n_estimators=70,
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
        
        'Gaussian Process (Kriging)': Pipeline([
            ('scaler', StandardScaler()),
            ('gpr', GaussianProcessRegressor(
                kernel=Matern(nu=1.5) + WhiteKernel(noise_level=1e-4),
                n_restarts_optimizer=2,
                random_state=random_state
            ))
        ]),

        'SVR (Otimizado)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(kernel='rbf', C=5.0, epsilon=0.005, gamma='scale'))
        ]),

        'Rede Neural (MLP Tanh + L-BFGS)': Pipeline([
            ('scaler', StandardScaler()),
            ('mlp', MLPRegressor(
                hidden_layer_sizes=(12, 6),
                activation='tanh',
                solver='lbfgs',
                alpha=0.1,
                max_iter=3000,
                random_state=random_state
            ))
        ]),

        'SVR (Baseline)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(kernel='rbf', C=1.0, epsilon=0.01))
        ]),

        'Rede Neural (MLP ReLU + Adam)': Pipeline([
            ('scaler', StandardScaler()),
            ('mlp', MLPRegressor(
                hidden_layer_sizes=(16, 8),
                activation='relu',
                solver='adam',
                alpha=0.01,
                max_iter=3000,
                learning_rate_init=0.01,
                random_state=random_state
            ))
        ])
    }


# =============================================================================
# 3. ACOPLADOR FÍSICO E RESTRIÇÕES TERMODINÂMICAS
# =============================================================================

class AcopladorFisico:
    """
    Aplica conservação de massa, estequiometria e monotonicidade temporal.
    """
    @staticmethod
    def acoplar_e_restringir(df: pd.DataFrame, delta_x_predito: np.ndarray) -> np.ndarray:
        x_pbm = df['X_pbm'].values
        eta = df['razao_molar_eta'].values
        
        # 1. Soma híbrida bruta
        x_raw = x_pbm + delta_x_predito
        
        # 2 e 3. Restrições estequiométricas e limites físicos [0, 1]
        x_restringido = np.zeros_like(x_raw)
        for i in range(len(x_raw)):
            teto = min(1.0, eta[i]) if eta[i] < 1.0 else 1.0
            x_restringido[i] = np.clip(x_raw[i], 0.0, teto)
            
        # 4. Monotonicidade temporal (dX/dt >= 0) por ensaio (preservando ordem original)
        df_temp = df.copy()
        df_temp['X_restringido'] = x_restringido
        
        x_final = []
        for ens_id, grp in df_temp.groupby('ensaio_id', sort=False):
            curva_monotonica = np.maximum.accumulate(grp['X_restringido'].values)
            x_final.extend(curva_monotonica)
            
        return np.array(x_final)


# =============================================================================
# 4. CLASSE DO MODELO GREY-BOX
# =============================================================================

class ModeloHibridoCinza:
    def __init__(self, nome_algoritmo: str = 'Gradient Boosting (GBDT)'):
        self.nome_algoritmo = nome_algoritmo
        self.modelos = obter_catalogo_modelos()
        if nome_algoritmo not in self.modelos:
            raise ValueError(
                f"Algoritmo '{nome_algoritmo}' não encontrado no catálogo. "
                f"Opções disponíveis: {list(self.modelos.keys())}"
            )
        self.regressor = self.modelos[nome_algoritmo]
        self.pbm = PopulationBalanceModel()
        self.acoplador = AcopladorFisico()
        self.is_treinado = False

    def preparar_dados_base(self, df_raw: pd.DataFrame) -> pd.DataFrame:
        """
        Carrega dados pré-calculados do cache (data/processed/dados_com_residuo.csv) ou calcula via PBM.
        """
        df_cache = carregar_dados_com_residuo()
        if df_cache is not None:
            if 'C_acid_pbm_mol_L' not in df_cache.columns:
                df_cache['C_acid_pbm_mol_L'] = df_cache['C_acid_0_mol_L'] * np.maximum(
                    0.0, 1.0 - (df_cache['X_pbm'] / df_cache['razao_molar_eta'])
                )
            return df_cache
        
        print("-> Simulando base White-Box (PBM) para os 16 ensaios...")
        df = df_raw.copy()
        x_pbm_lista = []
        c_acid_lista = []
        
        for ens_id, grupo in df.groupby("ensaio_id", sort=False):
            ca0 = float(grupo["C_acid_0_mol_L"].iloc[0])
            eta = float(grupo["razao_molar_eta"].iloc[0])
            tempos = grupo["tempo_min"].tolist()
            
            res = self.pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
            x_pbm_lista.extend(res["X_pbm"])
            c_acid_lista.extend(res["C_acid_mol_L"])
            
        df["X_pbm"] = x_pbm_lista
        df["C_acid_pbm_mol_L"] = c_acid_lista
        df["delta_X"] = df["X_zn_exp"] - df["X_pbm"]
        return df

    def treinar(self, df: pd.DataFrame):
        """Treina o regressor residual com 100% dos dados informados."""
        X, y, _ = preparar_features_e_alvo(df)
        self.regressor.fit(X, y)
        self.is_treinado = True
        print(f"-> Regressor residual '{self.nome_algoritmo}' treinado com sucesso!")

    def prever(self, df: pd.DataFrame) -> np.ndarray:
        """Predição híbrida com restrições termodinâmicas."""
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de realizar predições.")
            
        X, _, _ = preparar_features_e_alvo(df)
        delta_x_pred = self.regressor.predict(X)
        x_hibrido = self.acoplador.acoplar_e_restringir(df, delta_x_pred)
        return x_hibrido

    def prever_residuo(self, df: pd.DataFrame) -> np.ndarray:
        """Retorna apenas o resíduo delta_X previsto pelo regressor de Machine Learning."""
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de realizar predições.")
        X, _, _ = preparar_features_e_alvo(df)
        return self.regressor.predict(X)

    def validacao_cruzada_logo(self, df: pd.DataFrame) -> np.ndarray:
        """Executa validação cruzada Leave-One-Group-Out (LOGO-CV)."""
        X, y, grupos = preparar_features_e_alvo(df)
        logo = LeaveOneGroupOut()
        y_cv_pred = np.zeros(len(df))

        for train_idx, val_idx in logo.split(X, y, grupos):
            modelo_fold = clone(self.regressor)
            modelo_fold.fit(X[train_idx], y[train_idx])
            delta_val_pred = modelo_fold.predict(X[val_idx])
            
            df_val = df.iloc[val_idx].copy()
            y_cv_pred[val_idx] = self.acoplador.acoplar_e_restringir(df_val, delta_val_pred)

        return y_cv_pred

    def validacao_holdout_85_15(self, df: pd.DataFrame, test_size: float = 0.15, random_state: int = 42) -> Dict[str, Any]:
        """
        Executa uma divisão de Holdout Agrupado: ~85% de ensaios para treino (13 ensaios)
        vs ~15% de ensaios para teste cego (3 ensaios).
        """
        df_base = self.preparar_dados_base(df)
        X, y, grupos = preparar_features_e_alvo(df_base)
        gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=random_state)
        tr_idx, te_idx = next(gss.split(X, y, groups=grupos))
        
        reg = clone(self.regressor)
        reg.fit(X[tr_idx], y[tr_idx])
        
        df_tr = df_base.iloc[tr_idx].copy()
        df_te = df_base.iloc[te_idx].copy()
        
        pred_tr = reg.predict(X[tr_idx])
        pred_te = reg.predict(X[te_idx])
        
        y_tr_pred = self.acoplador.acoplar_e_restringir(df_tr, pred_tr)
        y_te_pred = self.acoplador.acoplar_e_restringir(df_te, pred_te)
        
        m_tr = calcular_metricas(df_tr['X_zn_exp'].values, y_tr_pred)
        m_te = calcular_metricas(df_te['X_zn_exp'].values, y_te_pred)
        
        ens_tr = sorted(df_tr['ensaio_id'].unique().tolist())
        ens_te = sorted(df_te['ensaio_id'].unique().tolist())
        
        return {
            'ensaios_treino': ens_tr,
            'ensaios_teste': ens_te,
            'pct_treino': len(ens_tr) / len(df['ensaio_id'].unique()) * 100.0,
            'pct_teste': len(ens_te) / len(df['ensaio_id'].unique()) * 100.0,
            'metricas_treino': m_tr,
            'metricas_teste': m_te,
            'df_treino': df_tr,
            'df_teste': df_te,
            'y_treino_real': df_tr['X_zn_exp'].values,
            'y_treino_pred': y_tr_pred,
            'y_teste_real': df_te['X_zn_exp'].values,
            'y_teste_pred': y_te_pred
        }

    def validacao_group_kfold_85_15(self, df: pd.DataFrame, n_splits: int = 6) -> Dict[str, Any]:
        """
        Executa validação cruzada por grupos (GroupKFold k=6) correspondendo a ~85% treino e ~15% teste
        por fold em todos os 128 pontos experimentais.
        """
        df_base = self.preparar_dados_base(df)
        X, y, grupos = preparar_features_e_alvo(df_base)
        gkf = GroupKFold(n_splits=n_splits)
        
        y_cv_pred = np.zeros(len(df_base))
        for tr_idx, te_idx in gkf.split(X, y, groups=grupos):
            reg = clone(self.regressor)
            reg.fit(X[tr_idx], y[tr_idx])
            df_fold_te = df_base.iloc[te_idx].copy()
            pred_te = reg.predict(X[te_idx])
            y_cv_pred[te_idx] = self.acoplador.acoplar_e_restringir(df_fold_te, pred_te)
            
        m_cv = calcular_metricas(df_base['X_zn_exp'].values, y_cv_pred)
        return {
            'metricas': m_cv,
            'y_pred_cv': y_cv_pred
        }


# =============================================================================
# 5. CÁLCULO DE MÉTRICAS ESTATÍSTICAS E CONSISTÊNCIA FÍSICA
# =============================================================================

def calcular_metricas(y_real: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """Calcula R², RMSE, MAE e verifica violações estequiométricas."""
    residuos = y_real - y_pred
    ss_res = np.sum(residuos ** 2)
    ss_tot = np.sum((y_real - np.mean(y_real)) ** 2)
    
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    rmse = np.sqrt(np.mean(residuos ** 2))
    mae = np.mean(np.abs(residuos))
    
    violacoes = np.sum((y_pred < -1e-4) | (y_pred > 1.0001))
    pct_violacao = (violacoes / len(y_pred)) * 100.0

    return {
        'R2': round(float(r2), 4),
        'RMSE': round(float(rmse), 4),
        'MAE': round(float(mae), 4),
        'Violacao_pct': round(float(pct_violacao), 2)
    }


# =============================================================================
# 6. BENCHMARK COMPARATIVO ENTRE TODOS OS MODELOS
# =============================================================================

def executar_benchmark_comparativo(df: pd.DataFrame, random_state: int = 42) -> pd.DataFrame:
    """
    Executa benchmark comparativo de todos os modelos do catálogo sob:
      1. Ajuste global (in-sample)
      2. Validação Cruzada Estrita LOGO-CV (16 folds independentes)
    """
    X, y, grupos = preparar_features_e_alvo(df)
    catalogo = obter_catalogo_modelos(random_state=random_state)
    logo = LeaveOneGroupOut()
    acoplador = AcopladorFisico()
    y_exp = df['X_zn_exp'].values
    
    resultados = []
    
    print("\n" + "=" * 90)
    print("BENCHMARK COMPARATIVO MULTI-MODELO -- VALIDAÇÃO CRUZADA ESTRITA (LOGO-CV)")
    print("=" * 90)
    print(f"{'Algoritmo / Variação':<35} | {'R² Ajuste':<10} | {'R² LOGO-CV':<11} | {'RMSE CV':<9} | {'MAE CV':<9} | {'Violação':<8}")
    print("-" * 90)
    
    for nome, modelo in catalogo.items():
        # Ajuste in-sample
        modelo_fit = clone(modelo)
        modelo_fit.fit(X, y)
        y_fit = acoplador.acoplar_e_restringir(df, modelo_fit.predict(X))
        m_fit = calcular_metricas(y_exp, y_fit)
        
        # LOGO-CV
        y_cv = np.zeros(len(df))
        for train_idx, val_idx in logo.split(X, y, grupos):
            m_fold = clone(modelo)
            m_fold.fit(X[train_idx], y[train_idx])
            pred_val = m_fold.predict(X[val_idx])
            df_val = df.iloc[val_idx].copy()
            y_cv[val_idx] = acoplador.acoplar_e_restringir(df_val, pred_val)
            
        m_cv = calcular_metricas(y_exp, y_cv)
        
        print(f"{nome:<35} | {m_fit['R2']:^10.4f} | {m_cv['R2']:^11.4f} | {m_cv['RMSE']:^9.4f} | {m_cv['MAE']:^9.4f} | {m_cv['Violacao_pct']:^8.1f}%")
        
        resultados.append({
            'Algoritmo': nome,
            'R2_Ajuste': m_fit['R2'],
            'R2_CV': m_cv['R2'],
            'RMSE_CV': m_cv['RMSE'],
            'MAE_CV': m_cv['MAE'],
            'Violacao_pct': m_cv['Violacao_pct']
        })
        
    print("=" * 90)
    return pd.DataFrame(resultados)


# =============================================================================
# 7. VISUALIZAÇÃO E GRÁFICOS DE PARIDADE
# =============================================================================

def plotar_resultados(df: pd.DataFrame, pasta_saida: str = "scripts/docs", nome_algoritmo: str = "Gradient Boosting (GBDT)"):
    """Gera o Gráfico de Paridade comparando White-Box Puro vs Grey-Box."""
    os.makedirs(pasta_saida, exist_ok=True)
    
    plt.figure(figsize=(6.5, 6), dpi=300)
    plt.plot([0, 1.05], [0, 1.05], 'k--', lw=1.8, label='Ideal ($y = x$)')
    plt.plot([0, 1.0], [0.05, 1.05], 'gray', lw=1.0, ls=':', label=r'Tolerância $\pm 5\%$')
    plt.plot([0.05, 1.05], [0, 1.0], 'gray', lw=1.0, ls=':')
    
    # White-Box
    m_pbm = calcular_metricas(df['X_zn_exp'].values, df['X_pbm'].values)
    plt.scatter(df['X_zn_exp'], df['X_pbm'], color='#d62728', alpha=0.6, s=35,
                label=f"White-Box Puro ($R^2 = {m_pbm['R2']:.4f}$)")
    
    # Grey-Box (LOGO-CV)
    m_cv = calcular_metricas(df['X_zn_exp'].values, df['X_hibrido_cv'].values)
    plt.scatter(df['X_zn_exp'], df['X_hibrido_cv'], color='#1f77b4', alpha=0.85, s=45,
                edgecolors='black', lw=0.5,
                label=f"Grey-Box [{nome_algoritmo}] ($R^2 = {m_cv['R2']:.4f}$)")
    
    plt.xlabel('Conversão Experimental de Zinco ($X_{Zn}^{exp}$)')
    plt.ylabel('Conversão Prevista ($X_{Zn}^{prev}$)')
    plt.title(f'Gráfico de Paridade: White-Box vs Grey-Box ({nome_algoritmo})')
    plt.xlim(-0.02, 1.05)
    plt.ylim(-0.02, 1.05)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(loc='upper left', frameon=True)
    plt.tight_layout()
    
    caminho_fig = os.path.join(pasta_saida, "paridade_whitebox_vs_greybox.png")
    plt.savefig(caminho_fig)
    plt.close()
    print(f"-> Gráfico de paridade salvo em: {caminho_fig}")


# =============================================================================
# 8. EXECUÇÃO PRINCIPAL
# =============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("BENCHMARK DO MODELO HÍBRIDO GREY-BOX -- TCC ENGENHARIA QUÍMICA UFMG")
    print("=" * 80)

    # 1. Carregar dados experimentais
    df_raw = carregar_dados_bancada()
    
    # 2. Inicializar o modelo Grey-Box com o regressor campeão
    algoritmo_selecionado = 'Gradient Boosting (GBDT)'
    hibrido = ModeloHibridoCinza(nome_algoritmo=algoritmo_selecionado)
    df = hibrido.preparar_dados_base(df_raw)
    
    # 3. Treinamento in-sample (Ajuste)
    hibrido.treinar(df)
    df['X_hibrido_ajuste'] = hibrido.prever(df)
    
    # 4. Validação Cruzada Estrita (Generalização em dados cegos)
    print(f"-> Executando validação cruzada LOGO-CV para '{algoritmo_selecionado}' (16 folds)...")
    df['X_hibrido_cv'] = hibrido.validacao_cruzada_logo(df)
    
    # 5. Avaliação Comparativa de Métricas do modelo ativo
    y_exp = df['X_zn_exp'].values
    m_pbm = calcular_metricas(y_exp, df['X_pbm'].values)
    m_ajuste = calcular_metricas(y_exp, df['X_hibrido_ajuste'].values)
    m_cv = calcular_metricas(y_exp, df['X_hibrido_cv'].values)
    
    print("\n" + "-" * 80)
    print(f"{'Métrica':<25} | {'White-Box Puro':<15} | {f'Grey-Box Ajuste':<18} | {f'Grey-Box LOGO-CV':<18}")
    print("-" * 80)
    print(f"{'R² (Determinação)':<25} | {m_pbm['R2']:<15.4f} | {m_ajuste['R2']:<18.4f} | {m_cv['R2']:<18.4f}")
    print(f"{'RMSE (Erro Médio Quad.)':<25} | {m_pbm['RMSE']:<15.4f} | {m_ajuste['RMSE']:<18.4f} | {m_cv['RMSE']:<18.4f}")
    print(f"{'MAE (Erro Médio Abs.)':<25} | {m_pbm['MAE']:<15.4f} | {m_ajuste['MAE']:<18.4f} | {m_cv['MAE']:<18.4f}")
    print(f"{'Violação Física (%)':<25} | {m_pbm['Violacao_pct']:<15.1f} | {m_ajuste['Violacao_pct']:<18.1f} | {m_cv['Violacao_pct']:<18.1f}")
    print("-" * 80)
    
    # 6. Gerar figura científica de paridade
    plotar_resultados(df, nome_algoritmo=algoritmo_selecionado)
    
    # 7. Executar o benchmark comparativo multi-modelo para todos os regressores do catálogo
    executar_benchmark_comparativo(df)
    print("=" * 80)
