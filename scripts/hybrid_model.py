"""
hybrid_model.py
---------------
Modelo Híbrido Cinza (Grey-Box) para Lixiviação de Calcina de Zinco em Reator Batelada.
Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho

Arquitetura: PARALELA COOPERATIVA (modelagem do resíduo; Sansana et al., 2024)
  1. Base física (White-Box): PBM (Método das Características, RRB), configurado por `kinetics.ConfigWB`.
  2. Machine Learning residual: aprende ΔX = X_exp - X_pbm a partir de atributos operacionais e de estado.
  3. Acoplador físico: X = restringir(X_pbm + ΔX_ML) com teto estequiométrico η_real/ν, X(0) = 0 e
     monotonia (regressão isotônica) -- as restrições são aplicadas SOBRE a saída bruta do ML, e a saída
     bruta é auditada (`AcopladorFisico.diagnosticar`).
  4. Validação: protocolos em `validacao.py` (GroupKFold repetido, LOGO, deixa-um-η-fora, ...).

Conjuntos de atributos do ML (R5/I6): 'A' = (t, η_real, C_A0); 'B' = A + X_pbm; 'C' = B + C_Af,pbm + [H+]_pbm.
O conjunto 'legado' (6 atributos com η nominal e S/L) existe só para reproduzir os números anteriores à revisão.

Domínio de aplicabilidade (I8): t <= 15 min, η_real em [0,53; 3,30], C_A0 em [0,1; 1,5] mol/L, T = 40 °C.
Fora dele o modelo avisa; o ML não tem informação ali e o híbrido degenera para o White-Box.
"""

import warnings
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, RegressorMixin, clone
from sklearn.ensemble import ExtraTreesRegressor, GradientBoostingRegressor, RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from sklearn.svm import SVR

from kinetics import ConfigWB
from population_balance import PopulationBalanceModel, MM_ZNO_G_MOL, TEOR_ZNO


# =============================================================================
# 1. ATRIBUTOS DO ML
# =============================================================================

CONJUNTOS_ENTRADA: Dict[str, List[str]] = {
    # Anterior à revisão: η nominal, S/L (função exata de C_A0/η) e atributos do PBM
    'legado': ['tempo_min', 'razao_molar_eta', 'C_acid_0_mol_L', 'razao_SL_g_L', 'X_pbm', 'C_acid_pbm_mol_L'],
    # Só as três variáveis independentes do planejamento (sem informação do PBM)
    'A': ['tempo_min', 'razao_molar_eta_real', 'C_acid_0_mol_L'],
    # A + conversão do White-Box (o ML "vê" onde a teoria está operando)
    'B': ['tempo_min', 'razao_molar_eta_real', 'C_acid_0_mol_L', 'X_pbm'],
    # B + estado químico do licor segundo o White-Box (exige força motriz [H+]: indisponível no White-Box de laboratório)
    'C': ['tempo_min', 'razao_molar_eta_real', 'C_acid_0_mol_L', 'X_pbm', 'C_acid_pbm_mol_L', 'H_pbm_mol_L'],
}
FEATURES_FISICAS = CONJUNTOS_ENTRADA['legado']

# Conjunto principal: definido pela ablação (experimentos_batelada.py, etapa "ablacao") e fixado aqui.
ENTRADAS_PRINCIPAL = 'B'

# White-Box principal: de laboratório (CAf livre, nu = 1, alpha = 0, ks(T) de Balarini); sem dados de conversão ou pH.
# A especiação por pH (V0, V1, V2) e nu a priori = 1,06 ficam como análises comparativas (etapa "whitebox", "nu").
CONFIG_WB_PRINCIPAL = ConfigWB.laboratorio()

# Domínio de aplicabilidade
DOMINIO_TEMPO_MAX_MIN = 15.0
DOMINIO_ETA_REAL = (0.53, 3.30)
DOMINIO_CA0_MOL_L = (0.1, 1.5)


def dados_de_entrada(df: pd.DataFrame, entradas: Any) -> np.ndarray:
    """Matriz de atributos na ordem do conjunto escolhido (nome em CONJUNTOS_ENTRADA ou lista de colunas)."""
    colunas = CONJUNTOS_ENTRADA[entradas] if isinstance(entradas, str) else list(entradas)
    faltam = [c for c in colunas if c not in df.columns]
    if faltam:
        raise KeyError(f"Colunas ausentes para o conjunto de atributos '{entradas}': {faltam}. "
                       f"Gere o quadro com construir_base()/construir_base_condicao().")
    X = df[colunas].to_numpy(dtype=float)
    if not np.all(np.isfinite(X)):
        raise ValueError(f"Atributos não finitos no conjunto '{entradas}' (a força motriz do White-Box é H2SO4?).")
    return X


def preparar_features_e_alvo(df: pd.DataFrame, entradas: Any = ENTRADAS_PRINCIPAL
                             ) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
    """Matriz de atributos X, alvo residual y = ΔX (se houver) e grupos (ensaio_id). Não altera `df`."""
    if (isinstance(entradas, str) and entradas == 'legado' and 'C_acid_pbm_mol_L' not in df.columns
            and 'razao_molar_eta_real' in df.columns):
        df = df.assign(C_acid_pbm_mol_L=df['C_acid_0_mol_L'] * np.maximum(
            0.0, 1.0 - df['X_pbm'] / df['razao_molar_eta_real']))
    X = dados_de_entrada(df, entradas)
    if 'delta_X' in df.columns:
        y = df['delta_X'].to_numpy(dtype=float)
    elif 'X_zn_exp' in df.columns and 'X_pbm' in df.columns:
        y = (df['X_zn_exp'] - df['X_pbm']).to_numpy(dtype=float)
    else:
        y = None
    grupos = df['ensaio_id'].to_numpy() if 'ensaio_id' in df.columns else None
    return X, y, grupos


# =============================================================================
# 2. BASE FÍSICA: QUADRO COM A SIMULAÇÃO DO WHITE-BOX
# =============================================================================

_CACHE_PBM: Dict[Tuple, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
_CACHE_PBM_MAX = 256


def _assinatura_dados(df: pd.DataFrame) -> int:
    """Identifica as ENTRADAS DO PBM (ensaio, tempo, CA0, eta real). Os dados medidos não entram na chave."""
    cols = ['ensaio_id', 'tempo_min', 'C_acid_0_mol_L', 'razao_molar_eta_real']
    return hash(df[cols].to_numpy(dtype=float).tobytes())


def construir_base(df_raw: pd.DataFrame, config: ConfigWB, rapido: bool = True,
                   usar_cache: bool = True) -> pd.DataFrame:
    """
    Simula o White-Box (PBM com `config`) em cada ensaio de `df_raw` e devolve uma CÓPIA com as colunas
    X_pbm, C_acid_pbm_mol_L, H_pbm_mol_L e (se houver X_zn_exp) delta_X = X_zn_exp - X_pbm.
    Sempre recalcula a física (sem cache em disco). O cache em memória guarda só a saída do PBM, por
    (config, entradas do PBM); delta_X é sempre recalculado dos dados recebidos.
    """
    chave = (config, _assinatura_dados(df_raw), rapido)
    df = df_raw.copy()
    if usar_cache and chave in _CACHE_PBM:
        x, c, h = (a.copy() for a in _CACHE_PBM[chave])
    else:
        n = len(df)
        x, c, h = np.zeros(n), np.zeros(n), np.full(n, np.nan)
        pbm = PopulationBalanceModel.para_config(config)
        simular = pbm.simular_ensaio_rapido if rapido else pbm.simular_ensaio
        ids = df['ensaio_id'].to_numpy()
        for g in pd.unique(ids):
            pos = np.where(ids == g)[0]
            sub = df.iloc[pos]
            r = simular(sub['tempo_min'].tolist(), float(sub['C_acid_0_mol_L'].iloc[0]),
                        float(sub['razao_molar_eta_real'].iloc[0]), config=config)
            x[pos], c[pos], h[pos] = r['X_pbm'], r['C_acid_mol_L'], r['C_H_mol_L']
        if usar_cache:
            if len(_CACHE_PBM) >= _CACHE_PBM_MAX:
                _CACHE_PBM.pop(next(iter(_CACHE_PBM)))
            _CACHE_PBM[chave] = (x.copy(), c.copy(), h.copy())
    df['X_pbm'] = x
    df['C_acid_pbm_mol_L'] = c
    df['H_pbm_mol_L'] = h
    if 'X_zn_exp' in df.columns:
        df['delta_X'] = df['X_zn_exp'] - df['X_pbm']
    return df


def construir_base_condicao(eta_real: float, C_acid_0_mol_L: float, tempos_min, config: ConfigWB,
                            ensaio_id: int = 0, eta_nominal: Optional[float] = None) -> pd.DataFrame:
    """Quadro (já com o White-Box simulado) para uma condição arbitrária, p.ex. curvas contínuas e mapas em η."""
    tempos = np.asarray(tempos_min, dtype=float)
    pbm = PopulationBalanceModel.para_config(config)
    r = pbm.simular_ensaio_rapido(tempos.tolist(), C_acid_0_mol_L, eta_real, config=config)
    sl = C_acid_0_mol_L * MM_ZNO_G_MOL / (eta_real * TEOR_ZNO)   # g/L = m_B0 / V
    return pd.DataFrame({
        'ensaio_id': ensaio_id,
        'tempo_min': tempos,
        'razao_molar_eta': eta_real if eta_nominal is None else eta_nominal,
        'razao_molar_eta_real': eta_real,
        'C_acid_0_mol_L': C_acid_0_mol_L,
        'razao_SL_g_L': sl,
        'X_pbm': r['X_pbm'],
        'C_acid_pbm_mol_L': r['C_acid_mol_L'],
        'H_pbm_mol_L': r['C_H_mol_L'],
    })


# =============================================================================
# 3. CATÁLOGO DE REGRESSORES
# =============================================================================

def selecionar_entradas_gp(X: np.ndarray) -> np.ndarray:
    """(LEGADO) Entradas do GP a partir de FEATURES_FISICAS: [log(1 + t), eta, CA0]."""
    return np.column_stack([np.log1p(X[:, 0]), X[:, 1], X[:, 2]])


def obter_catalogo_modelos(random_state: int = 42) -> Dict[str, Any]:
    """
    (LEGADO) Catálogo de regressores avaliados no TCC antes da revisão, para 'legado' (6 atributos).
    Mantido para reproduzir os números anteriores e para `continuous_reactor_model.py`.
    O benchmark atual usa `obter_especificacoes()`.
    """
    return {
        'Gradient Boosting (GBDT)': GradientBoostingRegressor(
            n_estimators=70, max_depth=3, learning_rate=0.06, random_state=random_state),
        'Extra Trees': ExtraTreesRegressor(
            n_estimators=100, max_depth=6, min_samples_split=2, random_state=random_state),
        'Random Forest': RandomForestRegressor(
            n_estimators=100, max_depth=5, min_samples_split=2, random_state=random_state),
        # Entradas independentes do planejamento fatorial (log t, eta, CA0): S/L, X_pbm e
        # C_acid_pbm são funções delas. Amplitude ajustável e escala própria por entrada (ARD).
        'Gaussian Process (Kriging)': Pipeline([
            ('entradas', FunctionTransformer(selecionar_entradas_gp)),
            ('scaler', StandardScaler()),
            ('gpr', GaussianProcessRegressor(
                kernel=ConstantKernel(1.0) * Matern(length_scale=[1.0, 1.0, 1.0], nu=1.5)
                       + WhiteKernel(noise_level=1e-3),
                normalize_y=True, n_restarts_optimizer=5, random_state=random_state))
        ]),
        'SVR (Otimizado)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(kernel='rbf', C=5.0, epsilon=0.005, gamma='scale'))]),
        'Rede Neural (MLP Tanh)': Pipeline([
            ('scaler', StandardScaler()),
            ('mlp', MLPRegressor(hidden_layer_sizes=(12, 6), activation='tanh', solver='lbfgs',
                                 alpha=0.1, max_iter=3000, random_state=random_state))]),
        'SVR (GridSearch)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(kernel='rbf', C=2.0, epsilon=0.01, gamma='scale'))]),
        'SVR (Baseline)': Pipeline([
            ('scaler', StandardScaler()),
            ('svr', SVR(kernel='rbf', C=1.0, epsilon=0.01))]),
        'Rede Neural (MLP ReLU + Adam)': Pipeline([
            ('scaler', StandardScaler()),
            ('mlp', MLPRegressor(hidden_layer_sizes=(16, 8), activation='relu', solver='adam',
                                 alpha=0.01, max_iter=3000, learning_rate_init=0.01, random_state=random_state))]),
    }


def _log1p_primeira_coluna(X):
    """Transformação do tempo (primeira coluna) para o GP: a cinética é rápida no início."""
    X = np.array(X, dtype=float, copy=True)
    X[:, 0] = np.log1p(X[:, 0])
    return X


class ComAlvoPadronizado(BaseEstimator, RegressorMixin):
    """
    Padroniza o alvo (média 0, desvio 1) antes de ajustar e desfaz na predição. Necessário para a caixa-preta
    (alvo = X em [0, 1], com média ~0,7): sem isso, SVR/MLP revertem para um intercepto arbitrário fora da região
    de treino, o que seria uma baseline injusta. Repassa `groups` quando o regressor interno o exige.
    """
    def __init__(self, regressor=None):
        self.regressor = regressor

    @property
    def precisa_grupos(self):
        return getattr(self.regressor, 'precisa_grupos', False)

    def fit(self, X, y, groups=None):
        y = np.asarray(y, dtype=float)
        self.media_, self.desvio_ = float(y.mean()), float(y.std()) or 1.0
        self.regressor_ = clone(self.regressor)
        alvo = (y - self.media_) / self.desvio_
        if self.precisa_grupos:
            self.regressor_.fit(X, alvo, groups=groups)
        else:
            self.regressor_.fit(X, alvo)
        return self

    def predict(self, X):
        return self.regressor_.predict(X) * self.desvio_ + self.media_


class AjusteAninhado(BaseEstimator, RegressorMixin):
    """
    Escolhe, por validação cruzada INTERNA agrupada por ensaio, o melhor de `candidatos` usando só os dados
    de treino recebidos, e reajusta o escolhido em todo o treino. Evita escolher hiperparâmetros com os
    mesmos ensaios usados para reportar o desempenho (CV aninhada).
    """
    precisa_grupos = True

    def __init__(self, candidatos=(), n_splits_max: int = 7):
        self.candidatos = candidatos
        self.n_splits_max = n_splits_max

    def fit(self, X, y, groups=None):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        if groups is None:
            raise ValueError("AjusteAninhado exige `groups` (ensaio_id) para a CV interna")
        grupos = np.asarray(groups)
        k = min(self.n_splits_max, len(np.unique(grupos)))
        if len(self.candidatos) == 1 or k < 2:
            melhor, scores = 0, [np.nan] * len(self.candidatos)
        else:
            splits = list(GroupKFold(n_splits=k).split(X, y, grupos))
            scores = []
            for cand in self.candidatos:
                sse, n = 0.0, 0
                for tr, va in splits:
                    m = clone(cand).fit(X[tr], y[tr])
                    e = y[va] - m.predict(X[va])
                    sse += float(e @ e)
                    n += len(va)
                scores.append(sse / n)
            melhor = int(np.argmin(scores))
        self.indice_melhor_ = melhor
        self.scores_internos_ = scores
        self.melhor_ = clone(self.candidatos[melhor]).fit(X, y)
        return self

    def predict(self, X):
        return self.melhor_.predict(np.asarray(X, dtype=float))


def fabricar_gp(n_entradas: int, random_state: int = 42):
    """Processo gaussiano: C * Matérn(ν = 3/2, ARD) + ruído branco, alvo normalizado, log(1+t) no tempo."""
    return Pipeline([
        ('entradas', FunctionTransformer(_log1p_primeira_coluna)),
        ('scaler', StandardScaler()),
        ('gpr', GaussianProcessRegressor(
            kernel=ConstantKernel(1.0) * Matern(length_scale=[1.0] * n_entradas, nu=1.5)
                   + WhiteKernel(noise_level=1e-3),
            normalize_y=True, n_restarts_optimizer=5, random_state=random_state))
    ])


def fabricar_extra_trees(n_entradas: int, random_state: int = 42):
    return ExtraTreesRegressor(n_estimators=100, max_depth=6, min_samples_split=2, random_state=random_state)


def fabricar_random_forest(n_entradas: int, random_state: int = 42):
    return RandomForestRegressor(n_estimators=100, max_depth=5, min_samples_split=2, random_state=random_state)


def fabricar_gbdt(n_entradas: int, random_state: int = 42):
    return GradientBoostingRegressor(n_estimators=70, max_depth=3, learning_rate=0.06, random_state=random_state)


def fabricar_svr_tunado(n_entradas: int, random_state: int = 42):
    """
    SVR RBF com alvo padronizado; C, γ e ε escolhidos por CV aninhada (grade 3 x 5 x 3 = 45; ε em unidades do alvo
    padronizado). Mesma grade para o híbrido (alvo ΔX) e para a caixa-preta (alvo X). Tempo em log(1 + t).
    """
    candidatos = [Pipeline([('entradas', FunctionTransformer(_log1p_primeira_coluna)), ('scaler', StandardScaler()),
                            ('svr', SVR(kernel='rbf', C=c, gamma=g, epsilon=e))])
                  for c in (1.0, 10.0, 100.0) for g in (0.1, 0.3, 1.0, 3.0, 10.0) for e in (0.01, 0.05, 0.2)]
    return ComAlvoPadronizado(AjusteAninhado(candidatos=candidatos))


def fabricar_mlp_tunada(n_entradas: int, random_state: int = 42):
    """MLP tanh + L-BFGS com alvo padronizado; arquitetura e regularização por CV aninhada (grade 2 x 2); log(1 + t)."""
    candidatos = [Pipeline([('entradas', FunctionTransformer(_log1p_primeira_coluna)), ('scaler', StandardScaler()),
                            ('mlp', MLPRegressor(hidden_layer_sizes=h, activation='tanh', solver='lbfgs',
                                                 alpha=a, max_iter=3000, random_state=random_state))])
                  for h in ((12, 6), (8,)) for a in (0.1, 1.0)]
    return ComAlvoPadronizado(AjusteAninhado(candidatos=candidatos))


@dataclass(frozen=True)
class ModeloSpec:
    """Nome, rótulos e fábrica de um regressor: ÚNICO lugar com essas informações (R8)."""
    nome: str
    slug: str
    titulo: str
    familia: str
    fabrica: Callable[[int, int], Any]      # (n_entradas, random_state) -> estimador

    def criar(self, n_entradas: int, random_state: int = 42):
        return self.fabrica(n_entradas, random_state)


def obter_especificacoes() -> Dict[str, ModeloSpec]:
    """Regressores do benchmark atual (chave = slug)."""
    specs = [
        ModeloSpec('Gaussian Process (Kriging)', 'gp', 'Processo Gaussiano (Matérn 3/2, ARD)',
                   'Processo Gaussiano', fabricar_gp),
        ModeloSpec('Extra Trees', 'extra_trees', 'Extra Trees', 'Árvores / Bagging', fabricar_extra_trees),
        ModeloSpec('Gradient Boosting (GBDT)', 'gbdt', 'Gradient Boosting', 'Árvores / Boosting', fabricar_gbdt),
        ModeloSpec('Random Forest', 'random_forest', 'Random Forest', 'Árvores / Bagging', fabricar_random_forest),
        ModeloSpec('SVR (tunado)', 'svr', 'SVR RBF (C, ε por CV aninhada)', 'Support Vector', fabricar_svr_tunado),
        ModeloSpec('Rede Neural (MLP Tanh)', 'mlp_tanh', 'MLP tanh + L-BFGS (CV aninhada)',
                   'Rede Neural', fabricar_mlp_tunada),
    ]
    return {s.slug: s for s in specs}


def especificacao_por_nome(nome: str) -> ModeloSpec:
    specs = obter_especificacoes()
    for s in specs.values():
        if nome in (s.nome, s.slug):
            return s
    raise ValueError(f"Algoritmo '{nome}' não encontrado. Opções: {[s.nome for s in specs.values()]}")


def ajustar_regressor(regressor, X, y, grupos=None):
    """Ajusta repassando `groups` apenas aos estimadores que o exigem (AjusteAninhado)."""
    if getattr(regressor, 'precisa_grupos', False):
        return regressor.fit(X, y, groups=grupos)
    return regressor.fit(X, y)


# =============================================================================
# 4. ACOPLADOR FÍSICO E RESTRIÇÕES
# =============================================================================

class AcopladorFisico:
    """
    Impõe à saída bruta do modelo (X_bruto = X_pbm + ΔX_ML, ou X direto na caixa-preta):
      1. 0 <= X <= teto, com teto = min(1, η_real / ν): o ácido é o limitante e consome ν mol por mol de ZnO
         (mesma estequiometria do PBM; o teto com η nominal da versão anterior era uma correção oculta);
      2. X(t = 0) = 0 (condição inicial conhecida);
      3. monotonia dX/dt >= 0 por ensaio, pela projeção L2 (regressão isotônica) -- não pelo máximo acumulado,
         que é viesado para cima e depende da ordem das linhas.
    modo='legado' reproduz a versão anterior (teto com η nominal, sem X(0) = 0, máximo acumulado em ordem de linha).
    """
    def __init__(self, nu: float = 1.0, modo: str = 'fisico'):
        if modo not in ('fisico', 'legado'):
            raise ValueError("modo deve ser 'fisico' ou 'legado'")
        self.nu = float(nu)
        self.modo = modo

    def teto(self, df: pd.DataFrame) -> np.ndarray:
        if self.modo == 'legado':
            eta = df['razao_molar_eta'].to_numpy(dtype=float)
            return np.where(eta < 1.0, np.minimum(1.0, eta), 1.0)
        if 'razao_molar_eta_real' not in df.columns:
            raise KeyError("O acoplador físico exige a coluna 'razao_molar_eta_real'")
        return np.minimum(1.0, df['razao_molar_eta_real'].to_numpy(dtype=float) / self.nu)

    def restringir(self, df: pd.DataFrame, x_bruto: np.ndarray) -> np.ndarray:
        x_bruto = np.asarray(x_bruto, dtype=float)
        teto = self.teto(df)
        x = np.clip(x_bruto, 0.0, teto)
        ids = df['ensaio_id'].to_numpy()

        if self.modo == 'legado':
            saida = np.empty_like(x)
            for g in pd.unique(ids):
                pos = np.where(ids == g)[0]
                saida[pos] = np.maximum.accumulate(x[pos])
            return saida

        t = df['tempo_min'].to_numpy(dtype=float)
        saida = np.empty_like(x)
        for g in pd.unique(ids):
            pos = np.where(ids == g)[0]
            ordem = pos[np.argsort(t[pos], kind='stable')]
            y = x[ordem].copy()
            y[t[ordem] <= 0.0] = 0.0
            if len(ordem) > 1:
                iso = IsotonicRegression(y_min=0.0, y_max=float(teto[ordem].max()), increasing=True)
                y = iso.fit_transform(t[ordem], y)
            saida[ordem] = y
        saida[t <= 0.0] = 0.0
        return saida

    def acoplar_e_restringir(self, df: pd.DataFrame, delta_x_predito: np.ndarray) -> np.ndarray:
        return self.restringir(df, df['X_pbm'].to_numpy(dtype=float) + np.asarray(delta_x_predito, dtype=float))

    def diagnosticar(self, df: pd.DataFrame, x_bruto: np.ndarray, tol: float = 1e-3) -> Dict[str, float]:
        """Auditoria da saída BRUTA (antes das restrições): onde o ML viola a física e quanto."""
        x_bruto = np.asarray(x_bruto, dtype=float)
        teto = self.teto(df)
        t = df['tempo_min'].to_numpy(dtype=float)
        ids = df['ensaio_id'].to_numpy()
        pos_t = t > 0.0
        n = int(pos_t.sum())
        acima = int(np.sum((x_bruto > teto + tol) & pos_t))
        abaixo = int(np.sum((x_bruto < -tol) & pos_t))
        quedas, passos = 0, 0
        for g in pd.unique(ids):
            pos = np.where(ids == g)[0]
            ordem = pos[np.argsort(t[pos], kind='stable')]
            d = np.diff(x_bruto[ordem])
            quedas += int(np.sum(d < -tol))
            passos += len(d)
        x0 = np.abs(x_bruto[~pos_t]) if np.any(~pos_t) else np.array([0.0])
        return {
            'pct_acima_teto': 100.0 * acima / max(n, 1),
            'pct_abaixo_zero': 100.0 * abaixo / max(n, 1),
            'pct_quedas': 100.0 * quedas / max(passos, 1),
            'x0_max_abs': float(x0.max()),
        }


# =============================================================================
# 5. MODELO GREY-BOX
# =============================================================================

class ModeloHibridoCinza:
    """
    Híbrido paralelo cooperativo: X = restringir(X_pbm + ΔX_ML). `entradas='legado'` + `ConfigWB.legado()`
    reproduzem o modelo anterior à revisão (acoplador no modo legado).
    """
    def __init__(self, nome_algoritmo: str = 'Gaussian Process (Kriging)', entradas: str = ENTRADAS_PRINCIPAL,
                 config_wb: Optional[ConfigWB] = None, acoplador: Optional[AcopladorFisico] = None,
                 random_state: int = 42):
        self.nome_algoritmo = nome_algoritmo
        self.entradas = entradas
        self.config_wb = config_wb if config_wb is not None else (
            ConfigWB.legado() if entradas == 'legado' else CONFIG_WB_PRINCIPAL)
        self.random_state = random_state
        self.acoplador = acoplador if acoplador is not None else AcopladorFisico(
            nu=self.config_wb.nu, modo='legado' if entradas == 'legado' else 'fisico')
        self.n_entradas = len(CONJUNTOS_ENTRADA[entradas])
        if entradas == 'legado':
            catalogo = obter_catalogo_modelos(random_state)
            if nome_algoritmo not in catalogo:
                raise ValueError(f"Algoritmo '{nome_algoritmo}' não encontrado no catálogo legado. "
                                 f"Opções: {list(catalogo.keys())}")
            self.regressor = catalogo[nome_algoritmo]
        else:
            self.regressor = especificacao_por_nome(nome_algoritmo).criar(self.n_entradas, random_state)
        self.is_treinado = False
        self.dominio_ = None

    def preparar_dados_base(self, df_raw: pd.DataFrame) -> pd.DataFrame:
        """Simula o White-Box configurado e devolve o quadro com X_pbm, C_acid_pbm_mol_L, H_pbm_mol_L, delta_X."""
        return construir_base(df_raw, self.config_wb)

    def treinar(self, df: pd.DataFrame):
        """Treina o regressor residual com 100% dos dados informados."""
        X, y, grupos = preparar_features_e_alvo(df, self.entradas)
        ajustar_regressor(self.regressor, X, y, grupos)
        self.is_treinado = True
        if {'razao_molar_eta_real', 'C_acid_0_mol_L', 'tempo_min'} <= set(df.columns):
            self.dominio_ = {
                'tempo': (float(df['tempo_min'].min()), float(df['tempo_min'].max())),
                'eta_real': (float(df['razao_molar_eta_real'].min()), float(df['razao_molar_eta_real'].max())),
                'ca0': (float(df['C_acid_0_mol_L'].min()), float(df['C_acid_0_mol_L'].max())),
            }
        return self

    def fora_do_dominio(self, df: pd.DataFrame) -> np.ndarray:
        """True onde a condição está fora do domínio visto no treino (extrapolação)."""
        if self.dominio_ is None:
            return np.zeros(len(df), dtype=bool)
        t = df['tempo_min'].to_numpy(dtype=float)
        eta = df['razao_molar_eta_real'].to_numpy(dtype=float)
        ca0 = df['C_acid_0_mol_L'].to_numpy(dtype=float)
        d = self.dominio_
        tol = 1e-9
        return ((t < d['tempo'][0] - tol) | (t > d['tempo'][1] + tol) |
                (eta < d['eta_real'][0] - 0.02) | (eta > d['eta_real'][1] + 0.1) |
                (ca0 < d['ca0'][0] - tol) | (ca0 > d['ca0'][1] + tol))

    def prever_residuo(self, df: pd.DataFrame) -> np.ndarray:
        """Apenas o resíduo ΔX previsto pelo regressor de Machine Learning."""
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de realizar predições.")
        X, _, _ = preparar_features_e_alvo(df, self.entradas)
        return self.regressor.predict(X)

    def prever(self, df: pd.DataFrame, avisar_dominio: bool = True) -> np.ndarray:
        """Predição híbrida com as restrições físicas."""
        if avisar_dominio:
            fora = self.fora_do_dominio(df)
            if fora.any():
                warnings.warn(f"{int(fora.sum())} de {len(df)} pontos fora do domínio de treino "
                              f"({self.dominio_}); o ML não tem informação ali e o resultado é extrapolação.",
                              UserWarning, stacklevel=2)
        return self.acoplador.acoplar_e_restringir(df, self.prever_residuo(df))

    def prever_com_incerteza(self, df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """(X_híbrido, desvio-padrão preditivo do ΔX). Só para o processo gaussiano."""
        if not self.is_treinado:
            raise RuntimeError("O modelo precisa ser treinado antes de realizar predições.")
        if not (isinstance(self.regressor, Pipeline) and 'gpr' in self.regressor.named_steps):
            raise NotImplementedError("Incerteza analítica disponível apenas para o Processo Gaussiano")
        X, _, _ = preparar_features_e_alvo(df, self.entradas)
        Xt = self.regressor[:-1].transform(X)
        media, desvio = self.regressor.named_steps['gpr'].predict(Xt, return_std=True)
        return self.acoplador.acoplar_e_restringir(df, media), desvio

    def validacao_cruzada_group_kfold(self, df: pd.DataFrame, n_splits: int = 8) -> np.ndarray:
        """GroupKFold sobre um quadro já simulado (White-Box fixo; para o protocolo completo use validacao.py)."""
        X, y, grupos = preparar_features_e_alvo(df, self.entradas)
        return self._oof(df, X, y, grupos, GroupKFold(n_splits=n_splits).split(X, y, grupos))

    def validacao_cruzada_logo(self, df: pd.DataFrame) -> np.ndarray:
        """Leave-One-Group-Out (15 ensaios de treino, 1 de teste), sobre um quadro já simulado."""
        X, y, grupos = preparar_features_e_alvo(df, self.entradas)
        return self._oof(df, X, y, grupos, LeaveOneGroupOut().split(X, y, grupos))

    def _oof(self, df, X, y, grupos, splits) -> np.ndarray:
        saida = np.zeros(len(df))
        for tr, va in splits:
            m = clone(self.regressor)
            ajustar_regressor(m, X[tr], y[tr], grupos[tr])
            saida[va] = self.acoplador.acoplar_e_restringir(df.iloc[va], m.predict(X[va]))
        return saida


# =============================================================================
# 6. MÉTRICAS (LEGADO)
# =============================================================================

def calcular_metricas(y_real: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """(LEGADO) R², RMSE, MAE e violações no valor já restringido. Para a avaliação atual use validacao.py."""
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


if __name__ == "__main__":
    print("Este módulo é uma biblioteca. Para reproduzir todos os resultados do TCC (batelada) execute:\n"
          "    python scripts/experimentos_batelada.py   # gera lop/resultados/*.csv\n"
          "    python scripts/figuras_batelada.py        # gera lop/resultados/figuras/*.png")
