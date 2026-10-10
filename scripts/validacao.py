"""
validacao.py
------------
Protocolos de validação, baselines e métricas do modelo híbrido da lixiviação em batelada.

Princípios (ver lop/docs/plano_validade.md):
  - Todas as dobras são feitas por ENSAIO INTEIRO (o PBM é recalculado a cada dobra, de modo que os
    parâmetros do White-Box que dependem dos dados, como nu, só usam ensaios de treino).
  - Quatro protocolos: GroupKFold(8) repetido, LOGO, deixa-um-nível-de-η-fora (extrapolação em η),
    deixa-um-nível-de-C_A0-fora e extrapolação temporal (treina t <= 5 min, prevê 15 min).
  - Métricas só em t > 0 (em t = 0 o resíduo é trivial); ganho sobre o White-Box e seu teto pelo piso de ruído.
  - Baselines obrigatórios: White-Box sozinho, tabela média por (η, t), Ridge e caixa-preta (mesmos regressores
    prevendo X direto, sem física). O híbrido só é "válido" onde supera todos eles.
  - Incerteza por bootstrap de ensaios (IC95%) e comparações pareadas; empates são declarados.
"""

import math
import os
import sys
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import norm

_DIR_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
if _DIR_SCRIPTS not in sys.path:
    sys.path.insert(0, _DIR_SCRIPTS)
# os processos filhos do joblib (loky) herdam o PYTHONPATH e conseguem importar estes módulos
os.environ['PYTHONPATH'] = _DIR_SCRIPTS + os.pathsep + os.environ.get('PYTHONPATH', '')

from kinetics import ConfigWB, NU_A_PRIORI, especiar_escalar, ph_de_h  # noqa: E402
from hybrid_model import (AcopladorFisico, ModeloSpec, ajustar_regressor, construir_base,  # noqa: E402
                          dados_de_entrada, obter_especificacoes, preparar_features_e_alvo)

NIVEIS_ETA = (0.5, 1.0, 1.5, 3.1)

# Instante (min) a partir do qual a curva de cada nível nominal de η já está no patamar (usado no piso de ruído)
PATAMAR_A_PARTIR_DE_MIN = {0.5: 1.0, 1.0: 2.0, 1.5: 2.0, 3.1: 3.0}


# =============================================================================
# 1. PISO DE RUÍDO
# =============================================================================

def piso_ruido(df: pd.DataFrame) -> Dict[str, float]:
    """
    Estima o desvio-padrão do erro experimental (σ) nos patamares das curvas, onde a conversão não muda
    mais e a dispersão é só ruído de medida (os dados são médias de triplicatas).

    sigma           : desvio-padrão combinado em torno da média de cada ensaio no patamar (graus de liberdade corrigidos)
    sigma_diferencas: estimador de von Neumann, sqrt(média(Δ²)/2) das diferenças entre amostras consecutivas
                      no patamar (insensível a uma deriva lenta)
    """
    ss, dof, d2, nd = 0.0, 0, 0.0, 0
    por_eta = {}
    for eta, t_min in PATAMAR_A_PARTIR_DE_MIN.items():
        sub = df[(df['razao_molar_eta'] == eta) & (df['tempo_min'] >= t_min)]
        dev = sub['X_zn_exp'] - sub.groupby('ensaio_id')['X_zn_exp'].transform('mean')
        k = sub['ensaio_id'].nunique()
        n = len(sub)
        ss_eta, dof_eta = float(np.sum(dev ** 2)), n - k
        por_eta[eta] = math.sqrt(ss_eta / dof_eta)
        ss += ss_eta
        dof += dof_eta
        for _, g in sub.sort_values(['ensaio_id', 'tempo_min']).groupby('ensaio_id'):
            d = np.diff(g['X_zn_exp'].to_numpy())
            d2 += float(np.sum(d ** 2))
            nd += len(d)
    return {'sigma': math.sqrt(ss / dof), 'sigma_diferencas': math.sqrt(d2 / (2 * nd)),
            'graus_liberdade': dof, 'por_eta': por_eta}


# =============================================================================
# 2. ESQUEMAS DE VALIDAÇÃO (sempre por ensaio inteiro)
# =============================================================================

@dataclass(frozen=True)
class Dobra:
    treino: Tuple[int, ...]            # ensaio_id usados no treino
    teste: Tuple[int, ...]             # ensaio_id avaliados
    rep: int = 0                       # repetição (partições aleatórias) ou sorteio (curva de aprendizado)
    rotulo: str = ''
    t_corte: Optional[float] = None    # extrapolação temporal: treina em t <= t_corte, testa em t > t_corte


def _ids(df: pd.DataFrame) -> np.ndarray:
    return np.array(sorted(df['ensaio_id'].unique()))


def esquema_gkf8_sklearn(df: pd.DataFrame, n_splits: int = 8) -> List[Dobra]:
    """GroupKFold determinístico do sklearn (a partição usada no TCC): 14 ensaios de treino, 2 de teste."""
    ids = df['ensaio_id'].to_numpy()
    dobras = []
    for k, (tr, va) in enumerate(GroupKFold(n_splits=n_splits).split(np.zeros((len(df), 1)), groups=ids)):
        dobras.append(Dobra(tuple(sorted(set(ids[tr]))), tuple(sorted(set(ids[va]))), 0, f'f{k}'))
    return dobras


def esquema_gkf8_repetido(df: pd.DataFrame, n_rep: int = 20, seed0: int = 0, n_splits: int = 8) -> List[Dobra]:
    """`n_rep` partições aleatórias dos 16 ensaios em 8 dobras de 2 ensaios (semente seed0 + rep)."""
    ids = _ids(df)
    dobras = []
    for rep in range(n_rep):
        perm = np.random.default_rng(seed0 + rep).permutation(ids)
        for k, teste in enumerate(np.array_split(perm, n_splits)):
            treino = np.setdiff1d(ids, teste)
            dobras.append(Dobra(tuple(int(i) for i in treino), tuple(int(i) for i in sorted(teste)), rep, f'r{rep}f{k}'))
    return dobras


def esquema_logo(df: pd.DataFrame) -> List[Dobra]:
    ids = _ids(df)
    return [Dobra(tuple(int(j) for j in ids if j != i), (int(i),), 0, f'ens{i}') for i in ids]


def esquema_deixa_nivel(df: pd.DataFrame, coluna: str) -> List[Dobra]:
    """Deixa fora, de uma vez, todos os ensaios de um nível de `coluna` (razao_molar_eta ou C_acid_0_mol_L)."""
    nivel_por_ensaio = df.groupby('ensaio_id')[coluna].first()
    ids = _ids(df)
    dobras = []
    for nivel in sorted(nivel_por_ensaio.unique()):
        teste = [int(i) for i in ids if nivel_por_ensaio[i] == nivel]
        treino = [int(i) for i in ids if nivel_por_ensaio[i] != nivel]
        dobras.append(Dobra(tuple(treino), tuple(teste), 0, f'{coluna}={nivel}'))
    return dobras


def esquema_extrapolacao_temporal(df: pd.DataFrame, t_corte: float = 5.0) -> List[Dobra]:
    ids = tuple(int(i) for i in _ids(df))
    return [Dobra(ids, ids, 0, f't<={t_corte:g}->t>{t_corte:g}', t_corte)]


def dobras_curva_aprendizado(df: pd.DataFrame, ns: Sequence[int], n_sorteios: int = 30, seed0: int = 1000,
                             n_teste: int = 2) -> List[Dobra]:
    """Em cada sorteio, os `n_teste` últimos ensaios da permutação são o teste; o treino usa os n primeiros."""
    ids = _ids(df)
    dobras = []
    for s in range(n_sorteios):
        perm = np.random.default_rng(seed0 + s).permutation(ids)
        teste = tuple(int(i) for i in sorted(perm[-n_teste:]))
        for n in ns:
            dobras.append(Dobra(tuple(int(i) for i in sorted(perm[:n])), teste, s, f'n={n}'))
    return dobras


# =============================================================================
# 3. nu (FATOR DE FECHAMENTO ESTEQUIOMÉTRICO) DENTRO DA DOBRA
# =============================================================================

def estimar_nu(resp: pd.DataFrame, ids_treino: Sequence[int]) -> float:
    """
    nu = mediana, sobre os ensaios de TREINO com η nominal <= 1,5, de (1 - CAf/CA0) * η_real / X_final.
    Atenção: CAf final da Tabela A1.1 é calculado a partir do pH (especiação ideal de ácido puro), não é
    titulação; nos ensaios com η = 0,5 o ácido é totalmente consumido e nu = η_real / X_final não depende dele.
    Sem ensaios elegíveis no treino, devolve o valor a priori pela composição da calcina.
    """
    r = resp[resp['ensaio_id'].isin(list(ids_treino)) & (resp['razao_molar_eta'] <= 1.5) & (resp['X_zn_final'] > 0.2)]
    if r.empty:
        return NU_A_PRIORI
    nu_i = (1.0 - r['C_acid_final_mol_L'] / r['C_acid_0_mol_L']) * r['razao_molar_eta_real'] / r['X_zn_final']
    return float(np.median(nu_i))


def escolher_nu(resp: pd.DataFrame, ids_treino: Sequence[int], modo: str, nu_fixo: Optional[float] = None) -> float:
    if modo == 'dobra':
        return estimar_nu(resp, ids_treino)
    if modo == 'a_priori':
        return NU_A_PRIORI
    if modo == 'fixo':
        if nu_fixo is None:
            raise ValueError("modo 'fixo' exige nu_fixo")
        return float(nu_fixo)
    raise ValueError("nu_modo deve ser 'dobra', 'a_priori' ou 'fixo'")


# =============================================================================
# 4. ESTIMADORES (baselines, híbrido paralelo, caixa-preta)
# =============================================================================

def fabricar_ridge(n_entradas: int, random_state: int = 42):
    return Pipeline([('scaler', StandardScaler()), ('ridge', Ridge(alpha=1.0))])


SPEC_RIDGE = ModeloSpec('Ridge', 'ridge', 'Regressão linear (Ridge)', 'Linear', fabricar_ridge)


class WhiteBoxPuro:
    """Só a física: X = X_pbm (depois das restrições do acoplador)."""
    tipo = 'white-box'

    def __init__(self, **_):
        self.nome = 'White-Box'

    def ajustar(self, base: pd.DataFrame):
        return self

    def prever_bruto(self, base: pd.DataFrame) -> np.ndarray:
        return base['X_pbm'].to_numpy(dtype=float)


class TabelaEtaT:
    """Sem ML: soma ao White-Box o resíduo médio dos ensaios de treino com o mesmo nível de η e o mesmo t."""
    tipo = 'baseline'

    def __init__(self, **_):
        self.nome = 'Tabela (η, t)'

    def ajustar(self, base: pd.DataFrame):
        self.tab_ = base.groupby(['razao_molar_eta', 'tempo_min'])['delta_X'].mean()
        self.tab_t_ = base.groupby('tempo_min')['delta_X'].mean()
        return self

    def prever_bruto(self, base: pd.DataFrame) -> np.ndarray:
        d = [self.tab_.get((e, t), self.tab_t_.get(t, 0.0)) for e, t in zip(base['razao_molar_eta'], base['tempo_min'])]
        return base['X_pbm'].to_numpy(dtype=float) + np.array(d, dtype=float)


class ResidualML:
    """Híbrido paralelo cooperativo: X = X_pbm + ΔX_ML(atributos). Com spec Ridge é só um baseline linear."""

    def __init__(self, spec: ModeloSpec, entradas: str, random_state: int = 42, rotulo: Optional[str] = None):
        self.spec, self.entradas, self.random_state = spec, entradas, random_state
        self.tipo = 'baseline' if spec.slug == 'ridge' else 'hibrido'
        self.nome = rotulo or f'Híbrido {spec.slug} [{entradas}]'

    def ajustar(self, base: pd.DataFrame):
        X, y, g = preparar_features_e_alvo(base, self.entradas)
        self.reg_ = self.spec.criar(X.shape[1], self.random_state)
        ajustar_regressor(self.reg_, X, y, g)
        return self

    def prever_bruto(self, base: pd.DataFrame) -> np.ndarray:
        X, _, _ = preparar_features_e_alvo(base, self.entradas)
        return base['X_pbm'].to_numpy(dtype=float) + self.reg_.predict(X)

    def prever_desvio(self, base: pd.DataFrame) -> np.ndarray:
        """Desvio-padrão preditivo (analítico), só para o processo gaussiano."""
        if not (isinstance(self.reg_, Pipeline) and 'gpr' in self.reg_.named_steps):
            raise NotImplementedError
        X, _, _ = preparar_features_e_alvo(base, self.entradas)
        _, sd = self.reg_.named_steps['gpr'].predict(self.reg_[:-1].transform(X), return_std=True)
        return sd


class CaixaPretaML:
    """Sem física: o mesmo regressor prevê X direto de (t, η_real, C_A0); as restrições do acoplador valem igual."""
    tipo = 'caixa-preta'

    def __init__(self, spec: ModeloSpec, random_state: int = 42):
        self.spec, self.random_state = spec, random_state
        self.nome = f'Caixa-preta {spec.slug}'

    def ajustar(self, base: pd.DataFrame):
        X = dados_de_entrada(base, 'A')
        self.reg_ = self.spec.criar(X.shape[1], self.random_state)
        ajustar_regressor(self.reg_, X, base['X_zn_exp'].to_numpy(dtype=float), base['ensaio_id'].to_numpy())
        return self

    def prever_bruto(self, base: pd.DataFrame) -> np.ndarray:
        return self.reg_.predict(dados_de_entrada(base, 'A'))


def criar_estimador(codigo: Tuple, random_state: int = 42):
    """
    codigo: ('wb',) | ('tabela',) | ('ridge', entradas) | ('hibrido', slug, entradas) | ('caixa', slug).
    São tuplas simples para poderem ser enviadas aos processos do joblib.
    """
    tipo = codigo[0]
    specs = obter_especificacoes()
    if tipo == 'wb':
        return WhiteBoxPuro()
    if tipo == 'tabela':
        return TabelaEtaT()
    if tipo == 'ridge':
        return ResidualML(SPEC_RIDGE, codigo[1], random_state, rotulo='Ridge (resíduo)')
    if tipo == 'hibrido':
        return ResidualML(specs[codigo[1]], codigo[2], random_state)
    if tipo == 'caixa':
        return CaixaPretaML(specs[codigo[1]], random_state)
    raise ValueError(f'código de estimador desconhecido: {codigo}')


# =============================================================================
# 5. EXECUÇÃO DE UM ESQUEMA (previsões fora-da-dobra, em formato longo)
# =============================================================================

def _rodar_dobra(df_raw: pd.DataFrame, resp: pd.DataFrame, dobra: Dobra, codigos: Sequence[Tuple],
                 cfg_base: ConfigWB, nu_modo: str, nu_fixo: Optional[float], esquema: str, random_state: int
                 ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    warnings.filterwarnings('ignore')
    nu = escolher_nu(resp, dobra.treino, nu_modo, nu_fixo)
    cfg = cfg_base.com_nu(nu)
    base = construir_base(df_raw, cfg)
    ids = base['ensaio_id'].to_numpy()
    t = base['tempo_min'].to_numpy(dtype=float)
    m_tr = np.isin(ids, dobra.treino)
    m_te = np.isin(ids, dobra.teste)
    if dobra.t_corte is not None:
        m_tr &= t <= dobra.t_corte
        m_te &= t > dobra.t_corte
    tr = base[m_tr].reset_index(drop=True)
    te = base[m_te].reset_index(drop=True)
    acop = AcopladorFisico(nu=cfg.nu)
    teto = acop.teto(te)

    linhas, diags = [], []
    for cod in codigos:
        est = criar_estimador(cod, random_state).ajustar(tr)
        xb = est.prever_bruto(te)
        xh = acop.restringir(te, xb)
        linhas.append(pd.DataFrame({
            'esquema': esquema, 'rep': dobra.rep, 'dobra': dobra.rotulo, 'modelo': est.nome, 'tipo': est.tipo,
            'ensaio_id': te['ensaio_id'].to_numpy(), 'tempo_min': te['tempo_min'].to_numpy(),
            'eta_nominal': te['razao_molar_eta'].to_numpy(), 'eta_real': te['razao_molar_eta_real'].to_numpy(),
            'ca0': te['C_acid_0_mol_L'].to_numpy(), 'X_exp': te['X_zn_exp'].to_numpy(), 'X_pbm': te['X_pbm'].to_numpy(),
            'teto': teto, 'x_bruto': xb, 'x_hat': xh, 'nu': nu}))
        d = acop.diagnosticar(te, xb)
        d.update({'esquema': esquema, 'rep': dobra.rep, 'dobra': dobra.rotulo, 'modelo': est.nome,
                  'n_pts': int((te['tempo_min'] > 0).sum())})
        diags.append(d)
    return pd.concat(linhas, ignore_index=True), pd.DataFrame(diags)


def executar_esquema(df_raw: pd.DataFrame, resp: pd.DataFrame, nome_esquema: str, dobras: Sequence[Dobra],
                     codigos: Sequence[Tuple], cfg_base: ConfigWB, nu_modo: str = 'dobra',
                     nu_fixo: Optional[float] = None, n_jobs: int = -1, random_state: int = 42
                     ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Devolve (previsões fora-da-dobra em formato longo, auditoria da saída bruta por dobra e modelo)."""
    args = [(df_raw, resp, d, list(codigos), cfg_base, nu_modo, nu_fixo, nome_esquema, random_state) for d in dobras]
    if n_jobs == 1 or len(args) == 1:
        res = [_rodar_dobra(*a) for a in args]
    else:
        with parallel_config(backend='loky', n_jobs=n_jobs, inner_max_num_threads=1):
            res = Parallel()(delayed(_rodar_dobra)(*a) for a in args)
    pred = pd.concat([r[0] for r in res], ignore_index=True)
    diag = pd.concat([r[1] for r in res], ignore_index=True)
    return pred, diag


# =============================================================================
# 6. MÉTRICAS (somente t > 0), GANHO, BOOTSTRAP DE ENSAIOS E COMPARAÇÕES PAREADAS
# =============================================================================

REF_WB = 'White-Box'


def _metricas(g: pd.DataFrame) -> Dict[str, float]:
    e = g['X_exp'].to_numpy() - g['x_hat'].to_numpy()
    sst = float(np.sum((g['X_exp'] - g['X_exp'].mean()) ** 2))
    return {'MSE': float(np.mean(e ** 2)), 'RMSE': float(np.sqrt(np.mean(e ** 2))), 'MAE': float(np.mean(np.abs(e))),
            'vies': float(-np.mean(e)), 'R2': float(1.0 - np.sum(e ** 2) / sst) if sst > 0 else float('nan'),
            'n': int(len(g))}


def metricas_por_repeticao(pred: pd.DataFrame, sigma: float) -> pd.DataFrame:
    """Uma linha por (esquema, rep, modelo): RMSE/MAE/viés/R² em t > 0, ganho sobre o White-Box e RMSE por nível de η."""
    p = pred[pred['tempo_min'] > 0]
    linhas = []
    for (esq, rep, mod), g in p.groupby(['esquema', 'rep', 'modelo'], sort=False):
        m = _metricas(g)
        m.update({'esquema': esq, 'rep': rep, 'modelo': mod, 'tipo': g['tipo'].iloc[0]})
        for eta in NIVEIS_ETA:
            ge = g[g['eta_nominal'] == eta]
            m[f'RMSE_eta{eta:g}'] = float(np.sqrt(np.mean((ge['X_exp'] - ge['x_hat']) ** 2))) if len(ge) else float('nan')
        linhas.append(m)
    t = pd.DataFrame(linhas)
    wb = t[t['modelo'] == REF_WB].set_index(['esquema', 'rep'])['MSE']
    t['MSE_WB'] = [wb.get((e, r), float('nan')) for e, r in zip(t['esquema'], t['rep'])]
    t['ganho'] = 1.0 - t['MSE'] / t['MSE_WB']
    t['teto_ganho'] = 1.0 - sigma ** 2 / t['MSE_WB']
    return t


def _sse_por_ensaio(pred: pd.DataFrame) -> pd.DataFrame:
    """SSE e N por (esquema, modelo, ensaio), média sobre as repetições (cada ensaio é testado 1x por repetição)."""
    p = pred[pred['tempo_min'] > 0].copy()
    p['se'] = (p['X_exp'] - p['x_hat']) ** 2
    por_rep = p.groupby(['esquema', 'modelo', 'rep', 'ensaio_id'])['se'].agg(['sum', 'count']).reset_index()
    return por_rep.groupby(['esquema', 'modelo', 'ensaio_id'])[['sum', 'count']].mean().reset_index()


def _indices_bootstrap(n_ensaios: int, n_boot: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).integers(0, n_ensaios, size=(n_boot, n_ensaios))


def resumo_benchmark(pred: pd.DataFrame, diag: pd.DataFrame, sigma: float, n_boot: int = 2000,
                     seed: int = 0) -> pd.DataFrame:
    """
    Uma linha por (esquema, modelo): média e desvio entre repetições, IC95% por bootstrap de ensaios,
    ganho/teto, RMSE por nível de η e a auditoria da saída bruta (violações ANTES do acoplador).
    """
    por_rep = metricas_por_repeticao(pred, sigma)
    cols_media = ['RMSE', 'MAE', 'vies', 'R2', 'ganho', 'teto_ganho'] + [f'RMSE_eta{e:g}' for e in NIVEIS_ETA]
    agg = por_rep.groupby(['esquema', 'modelo', 'tipo'], sort=False)
    res = agg[cols_media].mean()
    res['RMSE_dp_reps'] = agg['RMSE'].std(ddof=1)
    res['n_rep'] = agg['RMSE'].count()
    res['n_pontos_por_rep'] = agg['n'].mean()
    res = res.reset_index()

    # IC95% (bootstrap de ensaios) do RMSE
    sse = _sse_por_ensaio(pred)
    lo, hi = {}, {}
    for (esq, mod), g in sse.groupby(['esquema', 'modelo']):
        s, n = g['sum'].to_numpy(), g['count'].to_numpy()
        idx = _indices_bootstrap(len(g), n_boot, seed)
        rmse_b = np.sqrt(s[idx].sum(axis=1) / n[idx].sum(axis=1))
        lo[(esq, mod)], hi[(esq, mod)] = np.percentile(rmse_b, [2.5, 97.5])
    res['RMSE_ic95_inf'] = [lo[(e, m)] for e, m in zip(res['esquema'], res['modelo'])]
    res['RMSE_ic95_sup'] = [hi[(e, m)] for e, m in zip(res['esquema'], res['modelo'])]

    # auditoria da saída bruta
    d = diag.copy()
    d['peso'] = d['n_pts']
    def _pond(g, c): return float(np.average(g[c], weights=g['peso'])) if g['peso'].sum() > 0 else float('nan')
    aud = d.groupby(['esquema', 'modelo']).apply(lambda g: pd.Series({
        'pct_acima_teto': _pond(g, 'pct_acima_teto'), 'pct_abaixo_zero': _pond(g, 'pct_abaixo_zero'),
        'pct_quedas': _pond(g, 'pct_quedas'), 'x0_max_abs': float(g['x0_max_abs'].max())}), include_groups=False).reset_index()
    return res.merge(aud, on=['esquema', 'modelo'], how='left')


def comparacoes_pareadas(pred: pd.DataFrame, pares: Sequence[Tuple[str, str]], n_boot: int = 2000,
                         seed: int = 0) -> pd.DataFrame:
    """
    Para cada (A, B) e esquema: Δ = RMSE(A) - RMSE(B) com IC95% por bootstrap de ensaios (mesmas reamostragens
    para os dois modelos). `empate` = o IC contém 0. Para esquemas repetidos, informa também o desvio-padrão
    da diferença pareada entre repetições (regra do plano: empate se |Δ| < dp).
    """
    sse = _sse_por_ensaio(pred)
    por_rep = metricas_por_repeticao(pred, sigma=0.0)
    linhas = []
    for esq in sse['esquema'].unique():
        se = sse[sse['esquema'] == esq]
        for a, b in pares:
            ga, gb = se[se['modelo'] == a].set_index('ensaio_id'), se[se['modelo'] == b].set_index('ensaio_id')
            if ga.empty or gb.empty:
                continue
            ids = ga.index.intersection(gb.index)
            sa, na = ga.loc[ids, 'sum'].to_numpy(), ga.loc[ids, 'count'].to_numpy()
            sb, nb = gb.loc[ids, 'sum'].to_numpy(), gb.loc[ids, 'count'].to_numpy()
            idx = _indices_bootstrap(len(ids), n_boot, seed)
            d_b = np.sqrt(sa[idx].sum(1) / na[idx].sum(1)) - np.sqrt(sb[idx].sum(1) / nb[idx].sum(1))
            delta = math.sqrt(sa.sum() / na.sum()) - math.sqrt(sb.sum() / nb.sum())
            ic = np.percentile(d_b, [2.5, 97.5])
            ra = por_rep[(por_rep['esquema'] == esq) & (por_rep['modelo'] == a)].set_index('rep')['RMSE']
            rb = por_rep[(por_rep['esquema'] == esq) & (por_rep['modelo'] == b)].set_index('rep')['RMSE']
            dif = (ra - rb).dropna()
            dp = float(dif.std(ddof=1)) if len(dif) > 1 else float('nan')
            linhas.append({'esquema': esq, 'A': a, 'B': b, 'delta_RMSE': delta, 'ic95_inf': ic[0], 'ic95_sup': ic[1],
                           'empate_ic': bool(ic[0] <= 0.0 <= ic[1]), 'dp_dif_pareada_reps': dp,
                           'empate_dp': bool(abs(float(dif.mean())) < dp) if len(dif) > 1 else float('nan')})
    return pd.DataFrame(linhas)


def resumo_curva_aprendizado(pred: pd.DataFrame) -> pd.DataFrame:
    """RMSE (t > 0) nos ensaios fora do treino em função do número de ensaios de treino (rótulo 'n=..')."""
    p = pred[pred['tempo_min'] > 0].copy()
    p['n_treino'] = p['dobra'].str.extract(r'n=(\d+)')[0].astype(int)
    p['se'] = (p['X_exp'] - p['x_hat']) ** 2
    por_sorteio = p.groupby(['modelo', 'n_treino', 'rep'])['se'].mean().reset_index()
    por_sorteio['rmse'] = np.sqrt(por_sorteio['se'])
    out = por_sorteio.groupby(['modelo', 'n_treino']).agg(
        RMSE=('se', lambda s: float(np.sqrt(s.mean()))), RMSE_dp_sorteios=('rmse', 'std'),
        n_sorteios=('rmse', 'count')).reset_index()
    out['erro_padrao'] = out['RMSE_dp_sorteios'] / np.sqrt(out['n_sorteios'])
    return out


# =============================================================================
# 7. INTERVALOS DE PREDIÇÃO E COBERTURA EMPÍRICA
# =============================================================================

def _q_conforme(tr: pd.DataFrame, codigo: Tuple, acop: AcopladorFisico, nivel: float, random_state: int) -> float:
    """Quantil conforme (split/CV+): erro absoluto fora-do-ensaio no treino (LOGO interno), nível `nivel`."""
    erros = []
    ids = np.unique(tr['ensaio_id'].to_numpy())
    for g in ids:
        a = tr[tr['ensaio_id'] != g].reset_index(drop=True)
        b = tr[tr['ensaio_id'] == g].reset_index(drop=True)
        est = criar_estimador(codigo, random_state).ajustar(a)
        xh = acop.restringir(b, est.prever_bruto(b))
        erros.append(np.abs(b['X_zn_exp'].to_numpy() - xh)[b['tempo_min'].to_numpy() > 0])
    e = np.sort(np.concatenate(erros))
    k = min(int(math.ceil((len(e) + 1) * nivel)), len(e))
    return float(e[k - 1])


def _rodar_dobra_intervalos(df_raw, resp, dobra, slugs, entradas, cfg_base, nu_modo, nivel, esquema, random_state):
    warnings.filterwarnings('ignore')
    nu = escolher_nu(resp, dobra.treino, nu_modo)
    cfg = cfg_base.com_nu(nu)
    base = construir_base(df_raw, cfg)
    ids = base['ensaio_id'].to_numpy()
    tr = base[np.isin(ids, dobra.treino)].reset_index(drop=True)
    te = base[np.isin(ids, dobra.teste)].reset_index(drop=True)
    acop = AcopladorFisico(nu=cfg.nu)
    teto = acop.teto(te)
    z = float(norm.ppf(0.5 + nivel / 2.0))
    linhas = []
    for slug in slugs:
        cod = ('hibrido', slug, entradas)
        est = criar_estimador(cod, random_state).ajustar(tr)
        xh = acop.restringir(te, est.prever_bruto(te))
        metodos = []
        if slug == 'gp':
            sd = est.prever_desvio(te)
            metodos.append(('GP (analítico)', np.clip(xh - z * sd, 0.0, teto), np.clip(xh + z * sd, 0.0, teto)))
        q = _q_conforme(tr, cod, acop, nivel, random_state)
        metodos.append(('Conforme (LOGO interno)', np.clip(xh - q, 0.0, teto), np.clip(xh + q, 0.0, teto)))
        for metodo, lo, hi in metodos:
            linhas.append(pd.DataFrame({
                'esquema': esquema, 'modelo': est.nome, 'metodo': metodo, 'ensaio_id': te['ensaio_id'],
                'tempo_min': te['tempo_min'], 'eta_nominal': te['razao_molar_eta'], 'X_exp': te['X_zn_exp'],
                'x_hat': xh, 'inf': lo, 'sup': hi}))
    return pd.concat(linhas, ignore_index=True)


def cobertura_intervalos(df_raw, resp, nome_esquema, dobras, slugs, entradas, cfg_base, nivel=0.95,
                         nu_modo='dobra', n_jobs=-1, random_state=42) -> pd.DataFrame:
    """Pontos com intervalo [inf, sup] fora-da-dobra; a cobertura é a fração de X_exp dentro do intervalo (t > 0)."""
    args = [(df_raw, resp, d, list(slugs), entradas, cfg_base, nu_modo, nivel, nome_esquema, random_state) for d in dobras]
    if n_jobs == 1 or len(args) == 1:
        res = [_rodar_dobra_intervalos(*a) for a in args]
    else:
        with parallel_config(backend='loky', n_jobs=n_jobs, inner_max_num_threads=1):
            res = Parallel()(delayed(_rodar_dobra_intervalos)(*a) for a in args)
    return pd.concat(res, ignore_index=True)


def resumo_cobertura(pts: pd.DataFrame, nivel: float = 0.95) -> pd.DataFrame:
    p = pts[pts['tempo_min'] > 0].copy()
    p['cobre'] = (p['X_exp'] >= p['inf'] - 1e-12) & (p['X_exp'] <= p['sup'] + 1e-12)
    p['largura'] = p['sup'] - p['inf']
    r = p.groupby(['esquema', 'modelo', 'metodo']).agg(cobertura=('cobre', 'mean'), largura_media=('largura', 'mean'),
                                                     n=('cobre', 'count')).reset_index()
    r['nivel_nominal'] = nivel
    return r


# =============================================================================
# 8. VALIDAÇÃO TERMODINÂMICA DO WHITE-BOX (pH medido: única medida independente do estado ácido)
# =============================================================================

def tabela_termodinamica(resp: pd.DataFrame, config: ConfigWB, t_final_min: float = 15.0) -> pd.DataFrame:
    """
    Por ensaio: pH_0 previsto (solução lixiviante, sem Zn) e pH_f previsto de duas formas,
      - 'espec'  : especiação aplicada ao estado final calculado com X_final MEDIDO (análise de Zn) e nu da
                   config -> testa só a química (independe da cinética);
      - 'ponta'  : estado final simulado pelo PBM em t = 15 min -> testa o modelo completo.
    O pH é -log10[H+] (sem coeficientes de atividade). A medição foi feita na temperatura de reação (40 °C).
    """
    from population_balance import PopulationBalanceModel
    pbm = PopulationBalanceModel.para_config(config)
    k2, kz = config.k2, config.kz
    linhas = []
    for _, r in resp.iterrows():
        ca0, eta = float(r['C_acid_0_mol_L']), float(r['razao_molar_eta_real'])
        h0 = especiar_escalar(ca0, 2.0 * ca0, k2)['h']
        xf = float(r['X_zn_final'])
        caf = max(ca0 * (1.0 - config.nu * xf / eta), 0.0)
        hf = especiar_escalar(ca0, 2.0 * caf, k2, xf * ca0 / eta if kz else 0.0, kz)['h'] if caf > 1e-9 else float('nan')
        sim = pbm.simular_ensaio_rapido([0.0, t_final_min], ca0, eta, config=config)
        linhas.append({'ensaio_id': int(r['ensaio_id']), 'eta_nominal': r['razao_molar_eta'], 'ca0': ca0,
                       'pH0_med': r['pH_0'], 'pH0_prev': ph_de_h(h0),
                       'pHf_med': r['pH_final'], 'pHf_espec': ph_de_h(hf) if hf == hf else float('nan'),
                       'pHf_ponta': float(sim['pH_pbm'][-1]), 'CAf_pbm_15min': float(sim['C_acid_mol_L'][-1])})
    t = pd.DataFrame(linhas)
    # pH_f só é comparável quando sobra ácido (η >= 1); com η = 0,5 o pH final (5-6) é de hidrólise/tamponamento
    t['pHf_comparavel'] = t['eta_nominal'] >= 1.0
    return t


def erro_termodinamico(t: pd.DataFrame) -> Dict[str, float]:
    e0 = (t['pH0_prev'] - t['pH0_med']).to_numpy()
    m = t['pHf_comparavel'].to_numpy()
    ef_e = (t['pHf_espec'] - t['pHf_med']).to_numpy()[m]
    ef_p = (t['pHf_ponta'] - t['pHf_med']).to_numpy()[m]
    def rms(v): return float(np.sqrt(np.mean(np.asarray(v) ** 2)))
    return {'RMSE_pH0': rms(e0), 'RMSE_pHf_espec': rms(ef_e), 'vies_pHf_espec': float(np.mean(ef_e)),
            'RMSE_pHf_ponta': rms(ef_p), 'vies_pHf_ponta': float(np.mean(ef_p)),
            'RMSE_pH_combinado_espec': rms(np.concatenate([e0, ef_e])),
            'RMSE_pH_combinado_ponta': rms(np.concatenate([e0, ef_p]))}


def metricas_x_whitebox(df_raw: pd.DataFrame, config: ConfigWB) -> Dict[str, float]:
    """RMSE/MAE/R² (t > 0) e R² global do White-Box puro (sem acoplador) com a configuração dada."""
    b = construir_base(df_raw, config)
    y, x = b['X_zn_exp'].to_numpy(), b['X_pbm'].to_numpy()
    m = b['tempo_min'].to_numpy() > 0
    def r2(a, p): return float(1 - np.sum((a - p) ** 2) / np.sum((a - a.mean()) ** 2))
    return {'RMSE_t>0': float(np.sqrt(np.mean((y[m] - x[m]) ** 2))), 'MAE_t>0': float(np.mean(np.abs(y[m] - x[m]))),
            'R2_t>0': r2(y[m], x[m]), 'R2_todos': r2(y, x)}


def escada_whitebox(df_raw: pd.DataFrame, resp: pd.DataFrame, nu: float = NU_A_PRIORI) -> pd.DataFrame:
    """
    Degraus do White-Box, da formulação da tese até o modelo com especiação e nu. Métricas de X e, quando a
    força motriz é [H+], também os erros de pH. As linhas finais variam a especiação (V0, V1, V2).
    """
    degraus = [
        ('Tese (α = 0)', ConfigWB.tese(0.0)),
        ('Tese (α = 5500)', ConfigWB.tese(5500.0)),
        ('+ RRB em base mássica', ConfigWB(rrb='massica', ks_arrhenius=False, forca_motriz='H2SO4')),
        ('+ ks(T) de Arrhenius', ConfigWB(rrb='massica', ks_arrhenius=True, forca_motriz='H2SO4')),
        ('+ força motriz [H+] (V0: K2 = 0,0099, TCC)', ConfigWB(especiacao='V0')),
        ('V1: K2(T) por van\'t Hoff', ConfigWB(especiacao='V1')),
        ('V2: V1 + par iônico ZnSO4', ConfigWB(especiacao='V2')),
        ('V0 + nu', ConfigWB(especiacao='V0', nu=nu)),
        ('V1 + nu', ConfigWB(especiacao='V1', nu=nu)),
        ('V2 + nu', ConfigWB(especiacao='V2', nu=nu)),
    ]
    linhas = []
    for nome, cfg in degraus:
        r = {'degrau': nome, 'config': cfg.rotulo()}
        r.update(metricas_x_whitebox(df_raw, cfg))
        if cfg.forca_motriz == 'H+':
            r.update(erro_termodinamico(tabela_termodinamica(resp, cfg)))
        linhas.append(r)
    return pd.DataFrame(linhas)


def selecionar_whitebox(df_raw: pd.DataFrame, resp: pd.DataFrame, nu: float = NU_A_PRIORI) -> Tuple[ConfigWB, pd.DataFrame]:
    """
    Regra (D2): entre V0, V1 e V2 (todas com nu a priori), escolhe a de menor erro combinado de pH (pH_0 dos
    16 ensaios e pH_f dos 12 com η >= 1) pela especiação sobre o estado final MEDIDO. O ajuste de X só desempata.
    """
    linhas = []
    for esp in ('V0', 'V1', 'V2'):
        cfg = ConfigWB(especiacao=esp, nu=nu)
        r = {'especiacao': esp}
        r.update(erro_termodinamico(tabela_termodinamica(resp, cfg)))
        r.update(metricas_x_whitebox(df_raw, cfg))
        linhas.append(r)
    tab = pd.DataFrame(linhas).sort_values(['RMSE_pH_combinado_espec', 'RMSE_t>0']).reset_index(drop=True)
    return ConfigWB(especiacao=tab.loc[0, 'especiacao'], nu=nu), tab
