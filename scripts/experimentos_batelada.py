"""
experimentos_batelada.py
------------------------
PIPELINE ÚNICO dos resultados do TCC (batelada). Todo número do texto sai daqui:

    python scripts/experimentos_batelada.py                 # tudo (protocolo completo)
    python scripts/experimentos_batelada.py --rapido        # versão reduzida (desenvolvimento e testes)
    python scripts/experimentos_batelada.py --etapas whitebox,benchmark

Etapas (e arquivos gerados em lop/resultados/):
  dados          piso_ruido.csv
  whitebox       escada_whitebox.csv, selecao_whitebox.csv, termodinamica_por_ensaio.csv   (Fase 2, D2)
  ablacao        ablacao_entradas.csv                                                      (I6)
  benchmark      benchmark_batelada.csv, comparacoes_pareadas.csv, predicoes_oof.csv.gz    (I1, I2)
  nu             sensibilidade_nu.csv                                                      (proveniência de nu)
  aprendizado    curva_aprendizado.csv                                                     (I3)
  intervalos     cobertura_intervalos.csv, cobertura_intervalos_pontos.csv.gz              (I7)
  residuo        residuo_aprendido.csv, gp_escalas.csv, residuo_patamar_por_ensaio.csv     (I9)
  dominio        dominio_incerteza_eta.csv                                                 (I8)
  escada_hibrido escada_hibrido.csv                                                       (Fase 3)
  curvas         curvas_logo.csv (curvas contínuas fora-da-dobra, GP e ET, para as figuras)
Sempre grava config.json (versões, sementes, White-Box e atributos principais).

Reprodutibilidade: sementes fixas, floats gravados com 8 algarismos significativos, ordem estável.
"""

import argparse
import json
import os
import platform
import sys
import time
import warnings

import numpy as np
import pandas as pd

DIR_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
if DIR_SCRIPTS not in sys.path:
    sys.path.insert(0, DIR_SCRIPTS)

import validacao as V  # noqa: E402
from data_loader import carregar_dados_bancada, carregar_respostas_finais  # noqa: E402
from hybrid_model import (CONFIG_WB_PRINCIPAL, ENTRADAS_PRINCIPAL, CONJUNTOS_ENTRADA, DOMINIO_CA0_MOL_L,  # noqa: E402
                          DOMINIO_ETA_REAL, DOMINIO_TEMPO_MAX_MIN, ModeloHibridoCinza, construir_base,
                          construir_base_condicao, obter_especificacoes, preparar_features_e_alvo)
from kinetics import ConfigWB, NU_A_PRIORI  # noqa: E402

DIR_RESULTADOS = os.path.abspath(os.path.join(DIR_SCRIPTS, '..', 'resultados'))
SLUGS_PRINCIPAIS = ('gp', 'extra_trees', 'gbdt', 'random_forest', 'svr', 'mlp_tanh')
SLUGS_CAIXA_PRETA = ('gp', 'extra_trees', 'gbdt', 'svr')
ETAPAS = ('dados', 'whitebox', 'ablacao', 'benchmark', 'nu', 'aprendizado', 'intervalos', 'residuo',
          'dominio', 'escada_hibrido', 'curvas')


class Parametros:
    """Tamanho dos experimentos (completo ou reduzido)."""
    def __init__(self, rapido: bool):
        self.rapido = rapido
        self.n_rep_gkf = 3 if rapido else 20
        self.n_rep_ablacao = 3 if rapido else 10
        self.n_sorteios_curva = 4 if rapido else 30
        self.ns_curva = (4, 8, 14) if rapido else (4, 6, 8, 10, 12, 14)
        self.n_rep_escada = 2 if rapido else 10
        self.n_sorteios_escada = 3 if rapido else 15
        self.n_boot = 300 if rapido else 2000
        self.semente = 0


def salvar(df: pd.DataFrame, nome: str, saida: str):
    caminho = os.path.join(saida, nome)
    if nome.endswith('.gz'):
        df.to_csv(caminho, index=False, float_format='%.8g', encoding='utf-8',
                  compression={'method': 'gzip', 'mtime': 0})
    else:
        df.to_csv(caminho, index=False, float_format='%.8g', encoding='utf-8')
    print(f'   -> {os.path.relpath(caminho, os.path.dirname(DIR_RESULTADOS))} ({len(df)} linhas)', flush=True)


def codigos_benchmark(entradas: str):
    """Baselines + híbridos + caixas-pretas do benchmark principal."""
    cods = [('wb',), ('tabela',), ('ridge', entradas)]
    cods += [('hibrido', s, entradas) for s in SLUGS_PRINCIPAIS]
    cods += [('caixa', s) for s in SLUGS_CAIXA_PRETA]
    return cods


# -----------------------------------------------------------------------------
# Esquemas do benchmark
# -----------------------------------------------------------------------------

def esquemas_benchmark(df: pd.DataFrame, p: Parametros):
    """(nome, dobras, nu_modo). A extrapolação temporal usa nu a priori (nu da dobra usaria o estado final de 15 min)."""
    return [
        ('GKF8 (partição do TCC)', V.esquema_gkf8_sklearn(df), 'dobra'),
        (f'GKF8 repetido ({p.n_rep_gkf} partições)', V.esquema_gkf8_repetido(df, p.n_rep_gkf, p.semente), 'dobra'),
        ('LOGO (16 ensaios)', V.esquema_logo(df), 'dobra'),
        ('Deixa-um-η-fora (extrapolação em η)', V.esquema_deixa_nivel(df, 'razao_molar_eta'), 'dobra'),
        ('Deixa-um-C_A0-fora', V.esquema_deixa_nivel(df, 'C_acid_0_mol_L'), 'dobra'),
        ('Extrapolação temporal (t<=5 -> 15 min)', V.esquema_extrapolacao_temporal(df, 5.0), 'a_priori'),
    ]


# -----------------------------------------------------------------------------
# Regra de escolha dos atributos (ablação)
# -----------------------------------------------------------------------------

def escolher_entradas(resumo: pd.DataFrame, por_rep: pd.DataFrame, ordem=('A', 'B', 'C')) -> dict:
    """
    resumo: uma linha por (esquema, entradas, modelo) com RMSE médio. por_rep: RMSE por (rep, entradas) no GKF8
    repetido (média sobre os modelos). Escolhe o conjunto de menor RMSE médio; se algum conjunto MAIS SIMPLES
    (anterior em `ordem`) tiver diferença pareada média menor que o desvio-padrão da diferença pareada entre
    repetições, escolhe o mais simples (empate declarado).
    """
    score = resumo.groupby('entradas')['RMSE'].mean()
    melhor = score.idxmin()
    escolhido = melhor
    empates = []
    for cand in ordem:
        if cand == melhor or ordem.index(cand) >= ordem.index(melhor):
            continue
        d = (por_rep[cand] - por_rep[melhor]).dropna()
        if len(d) > 1 and abs(d.mean()) < d.std(ddof=1):
            empates.append(cand)
    if empates:
        escolhido = min(empates, key=ordem.index)
    return {'melhor_por_rmse': melhor, 'escolhido': escolhido, 'empates_mais_simples': empates,
            'rmse_medio': {k: float(v) for k, v in score.items()}}


# -----------------------------------------------------------------------------
# Etapas
# -----------------------------------------------------------------------------

def etapa_dados(ctx):
    pr = V.piso_ruido(ctx['df'])
    linhas = [{'grandeza': 'sigma (patamares)', 'valor': pr['sigma']},
              {'grandeza': 'sigma (diferenças, von Neumann)', 'valor': pr['sigma_diferencas']},
              {'grandeza': 'graus de liberdade', 'valor': pr['graus_liberdade']}]
    linhas += [{'grandeza': f'sigma eta={k:g}', 'valor': v} for k, v in pr['por_eta'].items()]
    salvar(pd.DataFrame(linhas), 'piso_ruido.csv', ctx['saida'])
    ctx['sigma'] = pr['sigma']


def etapa_whitebox(ctx):
    df, resp = ctx['df'], ctx['resp']
    salvar(V.escada_whitebox(df, resp), 'escada_whitebox.csv', ctx['saida'])
    cfg, tab = V.selecionar_whitebox(df, resp)
    salvar(tab, 'selecao_whitebox.csv', ctx['saida'])
    por_ensaio = []
    for esp in ('V0', 'V1', 'V2'):
        t = V.tabela_termodinamica(resp, ConfigWB(especiacao=esp, nu=NU_A_PRIORI))
        t.insert(0, 'especiacao', esp)
        por_ensaio.append(t)
    salvar(pd.concat(por_ensaio, ignore_index=True), 'termodinamica_por_ensaio.csv', ctx['saida'])
    ctx['wb_selecionado'] = cfg
    consistente = cfg.especiacao == CONFIG_WB_PRINCIPAL.especiacao
    ctx['consistencia']['whitebox'] = bool(consistente)
    print(f"   White-Box escolhido pela regra de pH: {cfg.rotulo()} "
          f"({'consistente' if consistente else 'DIVERGE'} com hybrid_model.CONFIG_WB_PRINCIPAL)", flush=True)
    if not consistente:
        warnings.warn('A regra de seleção do White-Box diverge de hybrid_model.CONFIG_WB_PRINCIPAL; atualize a constante.')


def etapa_ablacao(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    linhas, por_rep = [], {}
    esquemas = [(f'GKF8 repetido ({p.n_rep_ablacao})', V.esquema_gkf8_repetido(df, p.n_rep_ablacao, p.semente)),
                ('Deixa-um-η-fora', V.esquema_deixa_nivel(df, 'razao_molar_eta'))]
    for ent in ('A', 'B', 'C'):
        cods = [('wb',)] + [('hibrido', s, ent) for s in ('gp', 'extra_trees', 'gbdt', 'svr')]
        for nome, dobras in esquemas:
            pred, diag = V.executar_esquema(df, resp, nome, dobras, cods, CONFIG_WB_PRINCIPAL, 'dobra', None, ctx['n_jobs'])
            r = V.resumo_benchmark(pred, diag, ctx['sigma'], n_boot=p.n_boot)
            r = r[r['modelo'] != V.REF_WB].copy()
            r['entradas'] = ent
            r['esquema_ablacao'] = nome
            linhas.append(r)
            if nome.startswith('GKF8'):
                m = V.metricas_por_repeticao(pred, ctx['sigma'])
                m = m[m['modelo'] != V.REF_WB]
                por_rep[ent] = m.groupby('rep')['RMSE'].mean()
        print(f'   ablacao: conjunto {ent} concluido', flush=True)
    ab = pd.concat(linhas, ignore_index=True)
    ab['atributos'] = ab['entradas'].map(lambda e: ','.join(CONJUNTOS_ENTRADA[e]))
    resumo = ab[['entradas', 'modelo', 'esquema_ablacao', 'RMSE']]
    sel = escolher_entradas(resumo, pd.DataFrame(por_rep))
    ab['escolhido'] = ab['entradas'] == sel['escolhido']
    salvar(ab[['esquema_ablacao', 'entradas', 'atributos', 'modelo', 'RMSE', 'RMSE_dp_reps', 'RMSE_ic95_inf',
               'RMSE_ic95_sup', 'ganho', 'teto_ganho', 'RMSE_eta0.5', 'RMSE_eta1', 'RMSE_eta1.5', 'RMSE_eta3.1',
               'escolhido']], 'ablacao_entradas.csv', ctx['saida'])
    ctx['selecao_entradas'] = sel
    consistente = sel['escolhido'] == ENTRADAS_PRINCIPAL
    ctx['consistencia']['entradas'] = bool(consistente)
    print(f"   atributos: melhor por RMSE = {sel['melhor_por_rmse']}; escolhido = {sel['escolhido']} "
          f"(empates mais simples: {sel['empates_mais_simples']}) "
          f"({'consistente' if consistente else 'DIVERGE'} com hybrid_model.ENTRADAS_PRINCIPAL)", flush=True)
    if not consistente:
        warnings.warn('A regra de seleção dos atributos diverge de hybrid_model.ENTRADAS_PRINCIPAL; atualize a constante.')


def etapa_benchmark(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    cods = codigos_benchmark(ENTRADAS_PRINCIPAL)
    preds, diags, comps = [], [], []
    for nome, dobras, nu_modo in esquemas_benchmark(df, p):
        t0 = time.time()
        pred, diag = V.executar_esquema(df, resp, nome, dobras, cods, CONFIG_WB_PRINCIPAL, nu_modo, None, ctx['n_jobs'])
        preds.append(pred)
        diags.append(diag)
        resumo = V.resumo_benchmark(pred, diag, ctx['sigma'], n_boot=p.n_boot)
        hib = resumo[resumo['tipo'] == 'hibrido'].sort_values('RMSE')
        melhor = hib.iloc[0]['modelo']
        pares = [(melhor, m) for m in resumo['modelo'] if m != melhor]
        pares += [(f'Híbrido {s} [{ENTRADAS_PRINCIPAL}]', f'Caixa-preta {s}') for s in SLUGS_CAIXA_PRETA]
        c = V.comparacoes_pareadas(pred, pares, n_boot=p.n_boot)
        c['melhor_hibrido'] = melhor
        comps.append(c)
        print(f'   {nome}: melhor hibrido = {melhor} (RMSE {hib.iloc[0]["RMSE"]:.4f}) [{time.time() - t0:.0f}s]', flush=True)
    pred, diag = pd.concat(preds, ignore_index=True), pd.concat(diags, ignore_index=True)
    resumo = V.resumo_benchmark(pred, diag, ctx['sigma'], n_boot=p.n_boot)
    salvar(resumo, 'benchmark_batelada.csv', ctx['saida'])
    salvar(pd.concat(comps, ignore_index=True).drop_duplicates(['esquema', 'A', 'B']), 'comparacoes_pareadas.csv', ctx['saida'])
    salvar(pred, 'predicoes_oof.csv.gz', ctx['saida'])
    ctx['benchmark'] = resumo


def etapa_nu(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    cods = [('wb',), ('hibrido', 'gp', ENTRADAS_PRINCIPAL), ('hibrido', 'extra_trees', ENTRADAS_PRINCIPAL)]
    linhas = []
    esq = [('GKF8 (partição do TCC)', V.esquema_gkf8_sklearn(df)),
           ('Deixa-um-η-fora', V.esquema_deixa_nivel(df, 'razao_molar_eta'))]
    for nu_modo, nu_fixo, rotulo in (('dobra', None, 'nu estimado na dobra'), ('a_priori', None, f'nu a priori ({NU_A_PRIORI})'),
                                     ('fixo', 1.0, 'nu = 1 (sem consumo extra de ácido)')):
        for nome, dobras in esq:
            pred, diag = V.executar_esquema(df, resp, nome, dobras, cods, CONFIG_WB_PRINCIPAL, nu_modo, nu_fixo, ctx['n_jobs'])
            r = V.resumo_benchmark(pred, diag, ctx['sigma'], n_boot=p.n_boot)
            r.insert(0, 'nu_modo', rotulo)
            nus = pred.groupby('dobra')['nu'].first()
            r['nu_medio'] = float(nus.mean())
            r['nu_min'] = float(nus.min())
            r['nu_max'] = float(nus.max())
            linhas.append(r)
    salvar(pd.concat(linhas, ignore_index=True), 'sensibilidade_nu.csv', ctx['saida'])


def etapa_aprendizado(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    e = ENTRADAS_PRINCIPAL
    cods = [('wb',), ('tabela',)] + [('hibrido', s, e) for s in ('gp', 'extra_trees', 'gbdt')] + \
           [('caixa', s) for s in ('gp', 'extra_trees', 'gbdt')]
    dobras = V.dobras_curva_aprendizado(df, p.ns_curva, p.n_sorteios_curva)
    pred, _ = V.executar_esquema(df, resp, 'Curva de aprendizado', dobras, cods, CONFIG_WB_PRINCIPAL, 'dobra', None, ctx['n_jobs'])
    salvar(V.resumo_curva_aprendizado(pred), 'curva_aprendizado.csv', ctx['saida'])


def etapa_intervalos(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    esq = [('Deixa-um-η-fora', V.esquema_deixa_nivel(df, 'razao_molar_eta'))]
    if not p.rapido:
        esq.insert(0, ('LOGO (16 ensaios)', V.esquema_logo(df)))
    pts = []
    for nome, dobras in esq:
        pts.append(V.cobertura_intervalos(df, resp, nome, dobras, ('gp', 'extra_trees', 'gbdt'), ENTRADAS_PRINCIPAL,
                                          CONFIG_WB_PRINCIPAL, 0.95, 'dobra', ctx['n_jobs']))
    pts = pd.concat(pts, ignore_index=True)
    salvar(V.resumo_cobertura(pts), 'cobertura_intervalos.csv', ctx['saida'])
    salvar(pts, 'cobertura_intervalos_pontos.csv.gz', ctx['saida'])


def _treinar_com_tudo(ctx, slug, entradas=ENTRADAS_PRINCIPAL):
    cfg = CONFIG_WB_PRINCIPAL
    m = ModeloHibridoCinza(obter_especificacoes()[slug].nome, entradas=entradas, config_wb=cfg)
    base = m.preparar_dados_base(ctx['df'])
    m.treinar(base)
    return m, base


def etapa_residuo(ctx):
    """Dependência do resíduo aprendido nas variáveis de operação e escalas de comprimento do GP (I9)."""
    gp, base = _treinar_com_tudo(ctx, 'gp')
    et, _ = _treinar_com_tudo(ctx, 'extra_trees')
    cfg = CONFIG_WB_PRINCIPAL
    tempos = np.array([0.5, 1.0, 2.0, 5.0, 15.0])
    linhas = []
    for ca0 in (0.1, 0.5, 1.0, 1.5):
        for eta in np.round(np.arange(0.55, 3.31, 0.05), 3):
            g = construir_base_condicao(float(eta), ca0, tempos, cfg)
            dgp, det = gp.prever_residuo(g), et.prever_residuo(g)
            xgp, xet = gp.acoplador.acoplar_e_restringir(g, dgp), et.acoplador.acoplar_e_restringir(g, det)
            for i, t in enumerate(tempos):
                linhas.append({'ca0': ca0, 'eta_real': eta, 'tempo_min': t, 'X_pbm': g['X_pbm'].iloc[i],
                               'delta_gp': dgp[i], 'delta_et': det[i], 'X_hat_gp': xgp[i], 'X_hat_et': xet[i]})
    salvar(pd.DataFrame(linhas), 'residuo_aprendido.csv', ctx['saida'])
    # escalas do GP (em unidades padronizadas das entradas) para os conjuntos A, B e C
    esc = []
    for ent in ('A', 'B', 'C'):
        m, _ = _treinar_com_tudo(ctx, 'gp', ent)
        k = m.regressor.named_steps['gpr'].kernel_
        nomes = CONJUNTOS_ENTRADA[ent]
        for nome, ls in zip(nomes, k.k1.k2.length_scale):
            esc.append({'entradas': ent, 'atributo': nome, 'escala_de_comprimento': float(ls),
                        'amplitude': float(np.sqrt(k.k1.k1.constant_value)), 'ruido_branco': float(k.k2.noise_level)})
    salvar(pd.DataFrame(esc), 'gp_escalas.csv', ctx['saida'])
    # resíduo no patamar (15 min) por ensaio: o que o ML precisa explicar
    b = construir_base(ctx['df'], cfg)
    pat = b[b['tempo_min'] == 15.0][['ensaio_id', 'razao_molar_eta', 'razao_molar_eta_real', 'C_acid_0_mol_L',
                                    'X_zn_exp', 'X_pbm', 'delta_X']].reset_index(drop=True)
    salvar(pat, 'residuo_patamar_por_ensaio.csv', ctx['saida'])


def etapa_dominio(ctx):
    """Concordância entre modelos e incerteza do GP entre os níveis de η ensaiados, em C_A0 não amostrados (I8)."""
    cfg = CONFIG_WB_PRINCIPAL
    modelos = {s: _treinar_com_tudo(ctx, s)[0] for s in ('gp', 'extra_trees', 'gbdt', 'svr')}
    tempos = np.array([0.5, 2.0, 15.0])
    linhas = []
    for ca0 in (0.3, 0.75, 1.25):
        for eta in np.round(np.arange(0.55, 3.31, 0.05), 3):
            g = construir_base_condicao(float(eta), ca0, tempos, cfg)
            xs = {s: m.acoplador.acoplar_e_restringir(g, m.prever_residuo(g)) for s, m in modelos.items()}
            _, sd = modelos['gp'].prever_com_incerteza(g)
            mat = np.vstack(list(xs.values()))
            for i, t in enumerate(tempos):
                linhas.append({'ca0': ca0, 'eta_real': eta, 'tempo_min': t, 'X_pbm': g['X_pbm'].iloc[i],
                               'X_gp': xs['gp'][i], 'X_extra_trees': xs['extra_trees'][i], 'X_gbdt': xs['gbdt'][i],
                               'X_svr': xs['svr'][i], 'discordancia_max_min': float(mat[:, i].max() - mat[:, i].min()),
                               'sd_gp': float(sd[i])})
    salvar(pd.DataFrame(linhas), 'dominio_incerteza_eta.csv', ctx['saida'])


def etapa_escada_hibrido(ctx):
    """O híbrido sobre degraus do White-Box: mais física melhora a extrapolação em η e o desempenho com poucos dados?"""
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    e = ENTRADAS_PRINCIPAL
    cods = [('wb',), ('hibrido', 'gp', e), ('hibrido', 'extra_trees', e), ('caixa', 'gp'), ('caixa', 'extra_trees')]
    degraus = [
        ('1. Tese (α = 0)', ConfigWB.tese(0.0), 'fixo', 1.0),
        ('2. White-Box do TCC (V0, nu = 1)', ConfigWB(especiacao='V0'), 'fixo', 1.0),
        ('3. V0 + nu', ConfigWB(especiacao='V0'), 'dobra', None),
        ('4. V2 + nu (principal)', CONFIG_WB_PRINCIPAL, 'dobra', None),
    ]
    esq = [(f'GKF8 repetido ({p.n_rep_escada})', V.esquema_gkf8_repetido(df, p.n_rep_escada, p.semente)),
           ('Deixa-um-η-fora', V.esquema_deixa_nivel(df, 'razao_molar_eta')),
           ('Poucos dados (n = 4 ensaios de treino)', [d for d in V.dobras_curva_aprendizado(df, (4,), p.n_sorteios_escada)]),
           ('Poucos dados (n = 8 ensaios de treino)', [d for d in V.dobras_curva_aprendizado(df, (8,), p.n_sorteios_escada)])]
    linhas = []
    for nome_d, cfg, nu_modo, nu_fixo in degraus:
        for nome_e, dobras in esq:
            pred, diag = V.executar_esquema(df, resp, nome_e, dobras, cods, cfg, nu_modo, nu_fixo, ctx['n_jobs'])
            r = V.resumo_benchmark(pred, diag, ctx['sigma'], n_boot=p.n_boot)
            r.insert(0, 'degrau_whitebox', nome_d)
            linhas.append(r)
        print(f'   escada hibrida: {nome_d} concluido', flush=True)
    salvar(pd.concat(linhas, ignore_index=True), 'escada_hibrido.csv', ctx['saida'])


def etapa_curvas(ctx):
    """Curvas contínuas fora-da-dobra (LOGO) para o White-Box e os híbridos GP e ET, para as figuras."""
    df, resp = ctx['df'], ctx['resp']
    t_fino = np.round(np.linspace(0.0, 15.0, 151), 4)
    linhas = []
    for dobra in V.esquema_logo(df):
        nu = V.escolher_nu(resp, dobra.treino, 'dobra')
        cfg = CONFIG_WB_PRINCIPAL.com_nu(nu)
        base = construir_base(df, cfg)
        tr = base[base['ensaio_id'].isin(dobra.treino)].reset_index(drop=True)
        ens = dobra.teste[0]
        r = df[df['ensaio_id'] == ens].iloc[0]
        g = construir_base_condicao(float(r['razao_molar_eta_real']), float(r['C_acid_0_mol_L']), t_fino, cfg,
                                    ensaio_id=ens, eta_nominal=float(r['razao_molar_eta']))
        acop = V.AcopladorFisico(nu=cfg.nu)
        linha = pd.DataFrame({'ensaio_id': ens, 'tempo_min': t_fino, 'X_pbm': g['X_pbm'].to_numpy()})
        for slug in ('gp', 'extra_trees'):
            est = V.ResidualML(obter_especificacoes()[slug], ENTRADAS_PRINCIPAL).ajustar(tr)
            linha[f'X_{slug}'] = acop.restringir(g, est.prever_bruto(g))
            if slug == 'gp':
                sd = est.prever_desvio(g)
                linha['sd_gp'] = sd
        linhas.append(linha)
    salvar(pd.concat(linhas, ignore_index=True), 'curvas_logo.csv', ctx['saida'])


FUNCOES = {'dados': etapa_dados, 'whitebox': etapa_whitebox, 'ablacao': etapa_ablacao, 'benchmark': etapa_benchmark,
           'nu': etapa_nu, 'aprendizado': etapa_aprendizado, 'intervalos': etapa_intervalos,
           'residuo': etapa_residuo, 'dominio': etapa_dominio, 'escada_hibrido': etapa_escada_hibrido,
           'curvas': etapa_curvas}


def gravar_config(ctx):
    import sklearn, scipy
    cfg = {
        'versoes': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                    'scikit-learn': sklearn.__version__, 'scipy': scipy.__version__},
        'modo': 'rapido' if ctx['p'].rapido else 'completo',
        'sementes': {'modelos': 42, 'particoes': ctx['p'].semente, 'curva_aprendizado': 1000, 'bootstrap': 0},
        'white_box_principal': CONFIG_WB_PRINCIPAL.para_dict(),
        'nu_a_priori': NU_A_PRIORI,
        'atributos_principais': {'conjunto': ENTRADAS_PRINCIPAL, 'colunas': CONJUNTOS_ENTRADA[ENTRADAS_PRINCIPAL]},
        'dominio_de_aplicabilidade': {'tempo_max_min': DOMINIO_TEMPO_MAX_MIN, 'eta_real': DOMINIO_ETA_REAL,
                                      'ca0_mol_L': DOMINIO_CA0_MOL_L, 'T_celsius': 40.0},
        'parametros': {k: v for k, v in vars(ctx['p']).items()},
        'sigma_piso_ruido': ctx.get('sigma'),
        'consistencia_com_constantes': ctx['consistencia'],
    }
    with open(os.path.join(ctx['saida'], 'config.json'), 'w', encoding='utf-8') as f:
        json.dump(cfg, f, indent=2, ensure_ascii=False, sort_keys=True)


def main(argv=None):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    warnings.filterwarnings('ignore')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rapido', action='store_true', help='tamanho reduzido (desenvolvimento/testes)')
    ap.add_argument('--saida', default=DIR_RESULTADOS, help='pasta de saída (padrão: lop/resultados)')
    ap.add_argument('--n-jobs', type=int, default=-1, help='processos paralelos (-1 = todos)')
    ap.add_argument('--etapas', default=','.join(ETAPAS), help=f'lista separada por vírgulas de {ETAPAS}')
    a = ap.parse_args(argv)
    etapas = [e.strip() for e in a.etapas.split(',') if e.strip()]
    for e in etapas:
        if e not in FUNCOES:
            ap.error(f'etapa desconhecida: {e}')
    os.makedirs(a.saida, exist_ok=True)

    df, resp = carregar_dados_bancada(), carregar_respostas_finais()
    ctx = {'df': df, 'resp': resp, 'p': Parametros(a.rapido), 'saida': a.saida, 'n_jobs': a.n_jobs,
           'consistencia': {}, 'sigma': V.piso_ruido(df)['sigma']}
    print(f"PIPELINE BATELADA ({'rapido' if a.rapido else 'completo'}) | White-Box: {CONFIG_WB_PRINCIPAL.rotulo()} "
          f"| atributos: {ENTRADAS_PRINCIPAL}", flush=True)
    t_ini = time.time()
    for e in ETAPAS:                      # ordem fixa, independente da ordem pedida
        if e in etapas:
            t0 = time.time()
            print(f'[{e}]', flush=True)
            FUNCOES[e](ctx)
            print(f'[{e}] {time.time() - t0:.0f}s', flush=True)
    gravar_config(ctx)
    print(f'Concluído em {time.time() - t_ini:.0f}s. Resultados em {a.saida}', flush=True)


if __name__ == '__main__':
    main()
