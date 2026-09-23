"""
serial_hybrid_model.py
----------------------
Modelo Híbrido Serial (Grey-Box Cinético) para Lixiviação de Calcina de Zinco.
Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho

Arquitetura Híbrida Serial:
  1. Variáveis de Entrada Operacionais:
     - Razão molar estequiométrica (eta = H2SO4 / ZnO)
     - Concentração inicial de ácido (C_acid_0_mol_L)
  2. Modelo Black-Box (Regressor de Machine Learning):
     - Estima o parâmetro de amortecimento cinético da tese: alpha = f_ML(eta, C_A0)
     - Garante consistência física de amortecimento: alpha >= 0
  3. Modelo White-Box Fundamental (Balanço Populacional - PBM):
     - Recebe o parâmetro alpha predito pelo Blackbox
     - Resolve a taxa de retração superficial:
       v(t) = (2 / rho_s) * max(ks * C_Af(t) - alpha * [C_A0 - C_Af(t)], 0)
     - Integra o PBM com distribuição Granulométrica RRB analítica
     - Produz X(t) e C_Af(t) com 100% de consistência física e monotonicidade estrita.
  4. Métricas de Avaliação:
     - Coeficiente de Determinação (R²)
     - Coeficiente de Determinação Ajustado (R²_ajustado = 1 - (1 - R²)*(n - 1)/(n - p - 1))
     - Erro Médio Quadrático (RMSE) e Erro Médio Absoluto (MAE)
     - Verificação estrita de violações termodinâmicas e estequiométricas.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from typing import Dict, Any, Tuple, List, Optional
from scipy.optimize import minimize_scalar

# Garantir compatibilidade de saída UTF-8 no Windows
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
from sklearn.compose import TransformedTargetRegressor
from sklearn.pipeline import Pipeline
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.base import clone

# Importações dos módulos físicos locais
from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel


# =============================================================================
# 1. FUNÇÕES DE CÁLCULO DE MÉTRICAS E ESTATÍSTICA
# =============================================================================

def calcular_metricas_estatisticas(y_real: np.ndarray,
                                   y_pred: np.ndarray,
                                   p_parametros: int = 2,
                                   eta: Optional[np.ndarray] = None) -> Dict[str, float]:
    """
    Calcula métricas de desempenho estatístico e consistência física:
      - R² (Coeficiente de Determinação)
      - R²_ajustado = 1 - [(1 - R²)*(n - 1) / (n - p - 1)]
      - RMSE (Raiz do Erro Médio Quadrático)
      - MAE (Erro Médio Absoluto)
      - Violações de Limites Físicos (X < 0 ou X > 1)
      - Violações de Limites Estequiométricos (X > eta quando eta < 1.0)
    
    Parâmetros:
      y_real: valores experimentais (X_zn_exp)
      y_pred: valores preditos pelo modelo
      p_parametros: número de graus de liberdade / parâmetros explicativos ajustados
      eta: array opcional com a razão estequiométrica de cada ponto
    """
    y_real = np.asarray(y_real, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    n = len(y_real)
    
    residuos = y_real - y_pred
    ss_res = float(np.sum(residuos ** 2))
    ss_tot = float(np.sum((y_real - np.mean(y_real)) ** 2))
    
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    
    # R² Ajustado
    graus_liberdade = n - p_parametros - 1
    if graus_liberdade > 0 and ss_tot > 0:
        r2_ajustado = 1.0 - ((1.0 - r2) * (n - 1) / graus_liberdade)
    else:
        r2_ajustado = r2
        
    rmse = float(np.sqrt(np.mean(residuos ** 2)))
    mae = float(np.mean(np.abs(residuos)))
    
    violacoes_mascara = (y_pred < -1e-4) | (y_pred > 1.0001)
    if eta is not None:
        eta_arr = np.asarray(eta, dtype=float)
        violacoes_mascara = violacoes_mascara | ((eta_arr < 1.0) & (y_pred > eta_arr + 1e-4))
        
    violacoes_totais = int(np.sum(violacoes_mascara))
    pct_violacao = float((violacoes_totais / n) * 100.0)
    
    return {
        'R2': round(float(r2), 4),
        'R2_ajustado': round(float(r2_ajustado), 4),
        'RMSE': round(float(rmse), 4),
        'MAE': round(float(mae), 4),
        'Violacao_pct': round(float(pct_violacao), 2),
        'SS_res': round(float(ss_res), 6),
        'SS_tot': round(float(ss_tot), 6),
        'N': n,
        'p': p_parametros
    }



def verificar_monotonicidade(df: pd.DataFrame, coluna_pred: str = 'X_pred') -> float:
    """Calcula o percentual de passos temporais onde dX/dt < 0 (violação da irreversibilidade)."""
    violacoes = 0
    total_passos = 0
    for _, grp in df.groupby('ensaio_id', sort=False):
        vals = grp[coluna_pred].values
        diffs = np.diff(vals)
        violacoes += np.sum(diffs < -1e-4)
        total_passos += len(diffs)
    return float((violacoes / total_passos) * 100.0) if total_passos > 0 else 0.0


# =============================================================================
# 2. CATÁLOGO DE REGRESSORES BLACKBOX PARA ALPHA
# =============================================================================

def obter_catalogo_regressores_alpha(random_state: int = 42) -> Dict[str, Any]:
    """
    Retorna o catálogo de regressores Blackbox para estimar alpha a partir de (eta, C_A0):
      1. Gradient Boosting (GBDT) - Campeão com árvores aditivas regularizadas
      2. Polynomial Ridge - Regressão quadrática de 2ª ordem com penalidade L2 suave
      3. Random Forest - Floresta aleatória com agregação bootstrap
      4. Extra Trees - Árvores extremamente aleatorizadas
      5. Support Vector Regression (SVR) - Kernel RBF com margem insensível
      6. Rede Neural (MLP) - Perceptron multicamadas com ativação não-linear suave
    """
    return {
        'Gradient Boosting (GBDT)': GradientBoostingRegressor(
            n_estimators=45,
            max_depth=3,
            learning_rate=0.08,
            subsample=0.9,
            random_state=random_state
        ),
        
        'Polynomial Ridge': Pipeline([
            ('poly', PolynomialFeatures(degree=2, include_bias=False)),
            ('scaler', StandardScaler()),
            ('ridge', Ridge(alpha=1.5))
        ]),
        
        'Random Forest': RandomForestRegressor(
            n_estimators=60,
            max_depth=4,
            min_samples_split=2,
            random_state=random_state
        ),
        
        'Extra Trees': ExtraTreesRegressor(
            n_estimators=60,
            max_depth=4,
            min_samples_split=2,
            random_state=random_state
        ),
        
        'SVR (RBF)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(kernel='rbf', C=10000.0, epsilon=400.0, gamma='scale'))
        ]),
        
        'Rede Neural (MLP)': TransformedTargetRegressor(
            regressor=Pipeline([
                ('scaler', StandardScaler()),
                ('mlp', MLPRegressor(
                    hidden_layer_sizes=(8, 4),
                    activation='tanh',
                    solver='lbfgs',
                    alpha=0.05,
                    max_iter=3000,
                    random_state=random_state
                ))
            ]),
            transformer=StandardScaler()
        )
    }


# =============================================================================
# 3. CLASSE DO MODELO HÍBRIDO SERIAL (ALPHA-BLACKBOX -> WHITEBOX PBM)
# =============================================================================

class ModeloHibridoAlphaSerial:
    """
    Implementação da Arquitetura Híbrida Serial:
      Etapa 1 (Black-Box):
        Entrada: Condições operacionais de alimentação (eta, C_acid_0_mol_L).
        Saída: Parâmetro cinético de desaceleração/amortecimento empírico alpha [um/min].
        Restrição física: alpha >= 0 (passivação/amortecimento positivo).
      
      Etapa 2 (White-Box):
        Entrada: alpha predito pela Blackbox + tempos operacionais.
        Núcleo: Balanço Populacional (PBM) acoplado ao Shrinking Core Model (SCM).
        Saída: Trajetória completa de conversão X(t) e concentração de ácido livre C_Af(t).
    """
    def __init__(self,
                 nome_algoritmo: str = 'Gradient Boosting (GBDT)',
                 random_state: int = 42):
        self.nome_algoritmo = nome_algoritmo
        self.random_state = random_state
        self.catalogo = obter_catalogo_regressores_alpha(random_state=random_state)
        
        if nome_algoritmo not in self.catalogo:
            raise ValueError(
                f"Algoritmo '{nome_algoritmo}' não encontrado no catálogo. "
                f"Opções disponíveis: {list(self.catalogo.keys())}"
            )
            
        self.regressor = self.catalogo[nome_algoritmo]
        self.pbm = PopulationBalanceModel()
        self.alfas_otimos_calibrados: Dict[int, float] = {}
        self.tabela_ensaios_treino: Optional[pd.DataFrame] = None
        self.is_treinado = False

    def calibrar_alfas_otimos_ensaios(self,
                                      df_dados: pd.DataFrame,
                                      verbose: bool = False) -> Dict[int, float]:
        """
        Calcula o valor ótimo de alpha para cada ensaio que minimiza a soma
        dos erros quadráticos (SSE) entre o PBM e os dados experimentais.
        """
        if isinstance(df_dados, tuple):
            df_dados = df_dados[0]
        alfas = {}
        if verbose:
            print("-> Calibrando alfas ótimos ensaio a ensaio via minimização de SSE...")
            
        for ens_id, grp in df_dados.groupby('ensaio_id', sort=False):
            ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
            eta = float(grp['razao_molar_eta'].iloc[0])
            tempos = grp['tempo_min'].tolist()
            y_exp = grp['X_zn_exp'].values
            
            def loss_alpha(a_val):
                res = self.pbm.simular_ensaio(
                    tempos_min=tempos,
                    C_acid_0_mol_L=ca0,
                    eta=eta,
                    alpha=float(a_val)
                )
                return float(np.sum((y_exp - res['X_pbm']) ** 2))
            
            # Limite físico: amortecimento alpha >= 0 até 35000 um/min
            opt = minimize_scalar(loss_alpha, bounds=(0.0, 35000.0), method='bounded')
            alpha_opt = max(0.0, float(opt.x))
            alfas[int(ens_id)] = alpha_opt
            
            if verbose:
                print(f"   Ensaio {ens_id:2d} (eta={eta:3.1f}, CA0={ca0:3.1f} M) -> alpha* = {alpha_opt:8.1f} um/min")
                
        self.alfas_otimos_calibrados = alfas
        return alfas

    def extrair_atributos_operacionais(self, df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
        """
        Extrai as variáveis de entrada operacionais únicas por ensaio (eta, CA0)
        e o alvo de alpha calibrado.
        """
        resumo_ensaios = []
        for ens_id, grp in df.groupby('ensaio_id', sort=False):
            ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
            eta = float(grp['razao_molar_eta'].iloc[0])
            sl = float(grp['razao_SL_g_L'].iloc[0]) if 'razao_SL_g_L' in grp.columns else 0.0
            
            alpha_alvo = self.alfas_otimos_calibrados.get(int(ens_id), 0.0)
            resumo_ensaios.append({
                'ensaio_id': int(ens_id),
                'razao_molar_eta': eta,
                'C_acid_0_mol_L': ca0,
                'razao_SL_g_L': sl,
                'alpha_alvo': alpha_alvo
            })
            
        df_ens = pd.DataFrame(resumo_ensaios)
        X_op = df_ens[['razao_molar_eta', 'C_acid_0_mol_L']].values
        y_alpha = df_ens['alpha_alvo'].values
        return X_op, y_alpha, df_ens

    def treinar(self, df_dados: pd.DataFrame, verbose: bool = True):
        """
        Treina o modelo híbrido serial:
          1. Calibra os alfas ótimos para todos os ensaios da base.
          2. Treina o regressor Blackbox para mapear (eta, C_A0) -> alpha.
        """
        if isinstance(df_dados, tuple):
            df_dados = df_dados[0]
        if not self.alfas_otimos_calibrados:
            self.calibrar_alfas_otimos_ensaios(df_dados, verbose=verbose)
            
        X_op, y_alpha, df_ens = self.extrair_atributos_operacionais(df_dados)
        self.tabela_ensaios_treino = df_ens
        
        self.regressor.fit(X_op, y_alpha)
        self.is_treinado = True
        
        if verbose:
            print(f"-> Regressor Blackbox '{self.nome_algoritmo}' treinado com sucesso!")

    def prever_alpha(self, eta: float, C_acid_0_mol_L: float) -> float:
        """Estima o parâmetro alpha via Blackbox para uma condição de entrada."""
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de prever alpha.")
        X_in = np.array([[float(eta), float(C_acid_0_mol_L)]])
        alpha_pred = float(self.regressor.predict(X_in)[0])
        # Restrição termodinâmica: amortecimento não pode ser negativo
        return float(max(0.0, alpha_pred))

    def prever(self, df_dados: pd.DataFrame) -> Tuple[np.ndarray, Dict[int, float]]:
        """
        Executa a predição completa do modelo Híbrido Serial para os ensaios no DataFrame.
        Retorna:
          - Array de conversões preditas X(t) para cada registro de df_dados
          - Dicionário com o alpha predito para cada ensaio
        """
        if isinstance(df_dados, tuple):
            df_dados = df_dados[0]
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de realizar predições.")
            
        y_pred = []
        alfas_preditos = {}
        
        for ens_id, grp in df_dados.groupby('ensaio_id', sort=False):
            ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
            eta = float(grp['razao_molar_eta'].iloc[0])
            tempos = grp['tempo_min'].tolist()
            
            alpha_hat = self.prever_alpha(eta=eta, C_acid_0_mol_L=ca0)
            alfas_preditos[int(ens_id)] = alpha_hat
            
            res_pbm = self.pbm.simular_ensaio(
                tempos_min=tempos,
                C_acid_0_mol_L=ca0,
                eta=eta,
                alpha=alpha_hat
            )
            y_pred.extend(res_pbm['X_pbm'])
            
        return np.array(y_pred, dtype=float), alfas_preditos

    def validacao_cruzada_logo(self, df_dados: pd.DataFrame) -> Tuple[np.ndarray, Dict[int, float]]:
        """
        Executa validação cruzada Leave-One-Group-Out (LOGO-CV) em 16 folds:
          - Em cada fold k, o Blackbox é treinado nos 15 ensaios restantes
          - O Blackbox estima alpha_hat para o ensaio k deixado de fora
          - O Whitebox PBM simula o ensaio k com o alpha_hat cego.
        Garante avaliação estrita sem vazamento de dados (*data leakage*).
        """
        if isinstance(df_dados, tuple):
            df_dados = df_dados[0]
        if not self.alfas_otimos_calibrados:
            self.calibrar_alfas_otimos_ensaios(df_dados, verbose=False)
            
        X_op, y_alpha, df_ens = self.extrair_atributos_operacionais(df_dados)
        logo = LeaveOneGroupOut()
        grupos = df_ens['ensaio_id'].values
        
        alfas_cv = {}
        y_cv_dict = {}
        
        for train_idx, val_idx in logo.split(X_op, y_alpha, grupos):
            ens_val_id = int(grupos[val_idx[0]])
            
            # Treina regressor apenas com os outros 15 ensaios
            mod_fold = clone(self.regressor)
            mod_fold.fit(X_op[train_idx], y_alpha[train_idx])
            
            alpha_hat_val = float(max(0.0, mod_fold.predict(X_op[val_idx])[0]))
            alfas_cv[ens_val_id] = alpha_hat_val
            
            # Simula o ensaio de teste cego no PBM
            grp_val = df_dados[df_dados['ensaio_id'] == ens_val_id]
            ca0 = float(grp_val['C_acid_0_mol_L'].iloc[0])
            eta = float(grp_val['razao_molar_eta'].iloc[0])
            tempos = grp_val['tempo_min'].tolist()
            
            res_val = self.pbm.simular_ensaio(
                tempos_min=tempos,
                C_acid_0_mol_L=ca0,
                eta=eta,
                alpha=alpha_hat_val
            )
            y_cv_dict[ens_val_id] = res_val['X_pbm']
            
        # Recompor na ordem original do DataFrame
        y_cv_total = []
        for ens_id, grp in df_dados.groupby('ensaio_id', sort=False):
            y_cv_total.extend(y_cv_dict[int(ens_id)])
            
        return np.array(y_cv_total, dtype=float), alfas_cv


# =============================================================================
# 4. BENCHMARK DOS REGRESSORES BLACKBOX PARA ALPHA
# =============================================================================

def executar_benchmark_regressores_alpha(df_dados: pd.DataFrame,
                                         random_state: int = 42) -> pd.DataFrame:
    """
    Compara todos os algoritmos Blackbox candidatos para a estimação de alpha
    sob Ajuste Global (in-sample) e sob LOGO-CV (teste cego).
    """
    catalogo = obter_catalogo_regressores_alpha(random_state=random_state)
    y_exp = df_dados['X_zn_exp'].values
    
    # 1. Instância base para calibrar os alfas ótimos uma única vez
    modelo_base = ModeloHibridoAlphaSerial(nome_algoritmo='Gradient Boosting (GBDT)', random_state=random_state)
    modelo_base.calibrar_alfas_otimos_ensaios(df_dados, verbose=False)
    
    print("\n" + "=" * 105)
    print("BENCHMARK MULTI-MODELO DA ARQUITETURA HÍBRIDA SERIAL (ESTIMAÇÃO DE ALPHA VIA BLACKBOX -> WHITEBOX)")
    print("=" * 105)
    header = f"{'Regresssor Black-Box':<28} | {'R² Fit':<8} | {'R²_adj Fit':<11} | {'R² CV':<8} | {'R²_adj CV':<11} | {'RMSE CV':<8} | {'MAE CV':<8}"
    print(header)
    print("-" * 105)
    
    resultados = []
    
    for nome in catalogo.keys():
        hibrido_serial = ModeloHibridoAlphaSerial(nome_algoritmo=nome, random_state=random_state)
        # Compartilha a calibração de alfas ótimos para otimizar velocidade
        hibrido_serial.alfas_otimos_calibrados = modelo_base.alfas_otimos_calibrados
        
        # 1. Ajuste In-Sample
        hibrido_serial.treinar(df_dados, verbose=False)
        y_fit, _ = hibrido_serial.prever(df_dados)
        m_fit = calcular_metricas_estatisticas(y_exp, y_fit, p_parametros=2)
        
        # 2. LOGO-CV
        y_cv, _ = hibrido_serial.validacao_cruzada_logo(df_dados)
        m_cv = calcular_metricas_estatisticas(y_exp, y_cv, p_parametros=2)
        
        linha = (f"{nome:<28} | {m_fit['R2']:^8.4f} | {m_fit['R2_ajustado']:^11.4f} | "
                 f"{m_cv['R2']:^8.4f} | {m_cv['R2_ajustado']:^11.4f} | {m_cv['RMSE']:^8.4f} | {m_cv['MAE']:^8.4f}")
        print(linha)
        
        resultados.append({
            'Algoritmo_BlackBox': nome,
            'R2_Fit': m_fit['R2'],
            'R2_adj_Fit': m_fit['R2_ajustado'],
            'RMSE_Fit': m_fit['RMSE'],
            'MAE_Fit': m_fit['MAE'],
            'R2_CV': m_cv['R2'],
            'R2_adj_CV': m_cv['R2_ajustado'],
            'RMSE_CV': m_cv['RMSE'],
            'MAE_CV': m_cv['MAE'],
            'Violacao_CV_pct': m_cv['Violacao_pct']
        })
        
    print("=" * 105)
    return pd.DataFrame(resultados)


# =============================================================================
# 5. EXECUÇÃO PRINCIPAL E DEMONSTRAÇÃO
# =============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("EXECUÇÃO DO MODELO HÍBRIDO SERIAL (BLACKBOX -> ALPHA -> WHITEBOX PBM)")
    print("TCC Engenharia Química UFMG (2026)")
    print("=" * 80)
    
    # Carregar dados experimentais de bancada (128 pontos)
    df_bancada = carregar_dados_bancada()
    print(f"-> Base de dados carregada: {len(df_bancada)} pontos (16 ensaios de bancada).")
    
    # Inicializar e calibrar modelo campeão
    campeao = 'Gradient Boosting (GBDT)'
    modelo_serial = ModeloHibridoAlphaSerial(nome_algoritmo=campeao)
    modelo_serial.calibrar_alfas_otimos_ensaios(df_bancada, verbose=True)
    modelo_serial.treinar(df_bancada, verbose=True)
    
    # Ajuste In-Sample
    y_fit, alfas_pred_fit = modelo_serial.prever(df_bancada)
    df_bancada['X_serial_fit'] = y_fit
    m_fit = calcular_metricas_estatisticas(df_bancada['X_zn_exp'].values, y_fit, p_parametros=2)
    
    # Validação Cruzada Estrita (LOGO-CV)
    print("\n-> Executando validação cruzada Leave-One-Group-Out (LOGO-CV - 16 folds)...")
    y_cv, alfas_pred_cv = modelo_serial.validacao_cruzada_logo(df_bancada)
    df_bancada['X_serial_cv'] = y_cv
    m_cv = calcular_metricas_estatisticas(df_bancada['X_zn_exp'].values, y_cv, p_parametros=2)
    
    # Tabela de Alfas Comparativos
    print("\n" + "-" * 75)
    print(f"{'Ensaio':<8} | {'eta':<6} | {'CA0 (M)':<9} | {'alpha Otimo (um/min)':<20} | {'alpha Fit':<12} | {'alpha CV':<12}")
    print("-" * 75)
    for ens_id in range(1, 17):

        grp = df_bancada[df_bancada['ensaio_id'] == ens_id]
        eta = float(grp['razao_molar_eta'].iloc[0])
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        a_opt = modelo_serial.alfas_otimos_calibrados[ens_id]
        a_fit = alfas_pred_fit[ens_id]
        a_cv = alfas_pred_cv[ens_id]
        print(f"{ens_id:^8d} | {eta:^6.1f} | {ca0:^9.1f} | {a_opt:^18.1f} | {a_fit:^15.1f} | {a_cv:^15.1f}")
    print("-" * 75)
    
    # Executar benchmark multi-modelo
    df_bench = executar_benchmark_regressores_alpha(df_bancada)
