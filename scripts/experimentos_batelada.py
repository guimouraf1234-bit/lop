"""
experimentos_batelada.py
------------------------
PIPELINE ÚNICO dos resultados do TCC (batelada). Todo número do texto sai daqui:

    python scripts/experimentos_batelada.py                 # tudo (protocolo completo)
    python scripts/experimentos_batelada.py --rapido        # versão reduzida (desenvolvimento e testes)
    python scripts/experimentos_batelada.py --etapas whitebox,benchmark

Etapas (e arquivos gerados em lop/resultados/):
  dados          piso_ruido.csv
  whitebox       escada_whitebox.csv, selecao_whitebox.csv, termodinamica_por_ensaio.csv   (diagnóstico de pH; comparativo)
  ablacao        ablacao_entradas.csv                                                      (I6)
  benchmark      benchmark_batelada.csv, comparacoes_pareadas.csv, predicoes_oof.csv.gz    (I1, I2)
  nu             sensibilidade_nu.csv                                                      (proveniência de nu)
  aprendizado    curva_aprendizado.csv                                                     (I3)
  intervalos     cobertura_intervalos.csv, cobertura_intervalos_pontos.csv.gz              (I7)
  residuo        residuo_aprendido.csv, gp_escalas.csv, residuo_patamar_por_ensaio.csv     (I9)
  dominio        dominio_incerteza_eta.csv                                                 (I8)
  escada_hibrido escada_hibrido.csv                                                       (Fase 3)
  curvas         curvas_logo.csv (curvas contínuas fora-da-dobra, GP e ET, para as figuras)
  particao_comum particao_comum_*.csv: teste cego 7, 8, 14 + CV nas 4 dobras da Tabela B.3 (mesma partição e
                 métricas dos modelos serial e black-box do relatório, para comparação direta)
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
          'dominio', 'escada_hibrido', 'curvas', 'particao_comum')
# O White-Box principal tem nu fixo (sem estimativa nos dados de treino): nenhum dado de conversão entra no fator de fechamento
NU_MODO_PRINCIPAL = 'fixo'
NU_FIXO_PRINCIPAL = CONFIG_WB_PRINCIPAL.nu
# 'C' exige [H+] do White-Box; com força motriz CAf (principal) não há [H+] e o conjunto não existe
CONJUNTOS_DISPONIVEIS = tuple(e for e in ('A', 'B', 'C') if not (e == 'C' and CONFIG_WB_PRINCIPAL.forca_motriz == 'H2SO4'))


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
    """(nome, dobras). Todos com nu fixo do White-Box principal: nenhuma dobra estima nu nos dados."""
    return [
        ('GKF8 (partição do TCC)', V.esquema_gkf8_sklearn(df)),
        (f'GKF8 repetido ({p.n_rep_gkf} partições)', V.esquema_gkf8_repetido(df, p.n_rep_gkf, p.semente)),
        ('LOGO (16 ensaios)', V.esquema_logo(df)),
        ('Deixa-um-η-fora (extrapolação em η)', V.esquema_deixa_nivel(df, 'razao_molar_eta')),
        ('Deixa-um-C_A0-fora', V.esquema_deixa_nivel(df, 'C_acid_0_mol_L')),
        ('Extrapolação temporal (t<=5 -> 15 min)', V.esquema_extrapolacao_temporal(df, 5.0)),
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
    print(f"   diagnóstico: especiação de menor erro de pH = {cfg.rotulo()} (comparativo; "
          f"White-Box principal = {CONFIG_WB_PRINCIPAL.rotulo()})", flush=True)


def etapa_ablacao(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    linhas, por_rep = [], {}
    esquemas = [(f'GKF8 repetido ({p.n_rep_ablacao})', V.esquema_gkf8_repetido(df, p.n_rep_ablacao, p.semente)),
                ('Deixa-um-η-fora', V.esquema_deixa_nivel(df, 'razao_molar_eta'))]
    for ent in CONJUNTOS_DISPONIVEIS:
        cods = [('wb',)] + [('hibrido', s, ent) for s in ('gp', 'extra_trees', 'gbdt', 'svr')]
        for nome, dobras in esquemas:
            pred, diag = V.executar_esquema(df, resp, nome, dobras, cods, CONFIG_WB_PRINCIPAL,
                                            NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL, ctx['n_jobs'])
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
    for nome, dobras in esquemas_benchmark(df, p):
        t0 = time.time()
        pred, diag = V.executar_esquema(df, resp, nome, dobras, cods, CONFIG_WB_PRINCIPAL,
                                        NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL, ctx['n_jobs'])
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
    # Análise comparativa: o White-Box principal usa nu fixo (NU_FIXO_PRINCIPAL); aqui se mede o efeito de estimar nu
    for nu_modo, nu_fixo, rotulo in (('dobra', None, 'nu estimado na dobra (comparativo)'),
                                     ('a_priori', None, f'nu a priori ({NU_A_PRIORI}, comparativo)'),
                                     (NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL, 'nu = 1 (White-Box principal)')):
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
    pred, _ = V.executar_esquema(df, resp, 'Curva de aprendizado', dobras, cods, CONFIG_WB_PRINCIPAL,
                                 NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL, ctx['n_jobs'])
    salvar(V.resumo_curva_aprendizado(pred), 'curva_aprendizado.csv', ctx['saida'])


def etapa_intervalos(ctx):
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    esq = [('Deixa-um-η-fora', V.esquema_deixa_nivel(df, 'razao_molar_eta'))]
    if not p.rapido:
        esq.insert(0, ('LOGO (16 ensaios)', V.esquema_logo(df)))
    pts = []
    for nome, dobras in esq:
        pts.append(V.cobertura_intervalos(df, resp, nome, dobras, ('gp', 'extra_trees', 'gbdt'), ENTRADAS_PRINCIPAL,
                                          CONFIG_WB_PRINCIPAL, 0.95, NU_MODO_PRINCIPAL, ctx['n_jobs'],
                                          nu_fixo=NU_FIXO_PRINCIPAL))
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
    for ent in CONJUNTOS_DISPONIVEIS:
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
        ('3. V0 + nu (estimado na dobra, comparativo)', ConfigWB(especiacao='V0'), 'dobra', None),
        ('4. CAf, nu = 1 (principal)', CONFIG_WB_PRINCIPAL, NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL),
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
    cfg = CONFIG_WB_PRINCIPAL
    linhas = []
    for dobra in V.esquema_logo(df):
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


# -----------------------------------------------------------------------------
# Partição comum às abordagens do relatório (teste cego 7, 8, 14 + 4 dobras da Tabela B.3)
# -----------------------------------------------------------------------------

SEMENTES_COMITE = tuple(range(42, 52))
SLUGS_PARTICAO_COMUM = SLUGS_PRINCIPAIS + ('knn', 'hist_gbdt', 'adaboost', 'kernel_ridge', 'poly_ridge', 'gp_rbf')


def codigos_particao_comum(entradas: str):
    """Baselines + cada regressor como ML do híbrido e como caixa-preta (mesmo regressor, sem física)."""
    cods = [('wb',), ('tabela',), ('ridge', entradas)]
    cods += [('hibrido', s, entradas) for s in SLUGS_PARTICAO_COMUM]
    cods += [('caixa', s) for s in SLUGS_PARTICAO_COMUM]
    return cods      # as mesmas do modelo serial e do black-box (42 a 51)
PASSO_MALHA_MIN = 0.25                      # malha da inconsistência física (como no modelo serial)
TOL_INCONSISTENCIA = 1e-3                   # mesma tolerância de AcopladorFisico.diagnosticar


def _metricas_pp(y, yhat) -> dict:
    """Métricas do relatório, em pontos percentuais de conversão (R² adimensional)."""
    y, yhat = np.asarray(y, float), np.asarray(yhat, float)
    e = y - yhat
    sst = float(np.sum((y - y.mean()) ** 2))
    return {'pontos': int(len(y)), 'R2': float(1 - np.sum(e ** 2) / sst) if sst > 0 else float('nan'),
            'RMSE_pp': 100 * float(np.sqrt(np.mean(e ** 2))), 'MAE_pp': 100 * float(np.mean(np.abs(e))),
            'erro_max_pp': 100 * float(np.max(np.abs(e)))}


def _pct_inconsistente(malha: pd.DataFrame, x: np.ndarray, teto: np.ndarray, tol: float = TOL_INCONSISTENCIA) -> float:
    """% dos nós (t > 0) com conversão decrescente em relação ao nó anterior, negativa ou acima do teto estequiométrico."""
    d = pd.DataFrame({'e': malha['ensaio_id'].to_numpy(), 't': malha['tempo_min'].to_numpy(dtype=float),
                      'x': np.asarray(x, float), 'teto': teto})
    d = d[d['t'] > 0].sort_values(['e', 't'])
    queda = d.groupby('e')['x'].diff() < -tol
    ruim = queda | (d['x'] > d['teto'] + tol) | (d['x'] < -tol)
    return 100.0 * float(ruim.mean())


def _ajustar_e_prever(cod, semente, tr, alvos):
    """Ajusta um estimador nos ensaios de treino e devolve a saída BRUTA em cada quadro de `alvos`."""
    warnings.filterwarnings('ignore')
    est = V.criar_estimador(cod, semente).ajustar(tr)
    return [est.prever_bruto(a) for a in alvos]


def etapa_particao_comum(ctx):
    """
    Mesma partição e mesma avaliação dos modelos serial e black-box do relatório (Apêndice B, Seção 5.4):
      1. seleção do regressor só pela CV nos 13 ensaios de treino (4 dobras da Tabela B.3, 10 sementes);
      2. modelo final = comitê das 10 sementes (média da saída bruta, depois o acoplador), treinado nos 13 ensaios;
      3. teste cego (Ensaios 7, 8 e 14) consultado uma única vez; métricas em p.p., só em t > 0;
      4. inconsistência física na malha de 0,25 min, antes e depois do acoplador.
    Todos os modelos do benchmark passam pelo mesmo protocolo (baselines e caixas-pretas incluídos).
    """
    df, resp, p = ctx['df'], ctx['resp'], ctx['p']
    sementes = SEMENTES_COMITE[:2] if p.rapido else SEMENTES_COMITE
    cods = codigos_particao_comum(ENTRADAS_PRINCIPAL)
    nome_por_cod = {c: V.criar_estimador(c).nome for c in cods}
    treino, teste = V.ensaios_treino_comum(df), V.ENSAIOS_TESTE_CEGO
    dobras = V.esquema_cv_tabela_b3(df)
    rotulo_dobra = {d.rotulo: d.teste for d in dobras}

    # 1. Validação cruzada (somente treino) ---------------------------------
    preds, diags = [], []
    for s in sementes:
        pr, dg = V.executar_esquema(df, resp, 'CV 4 dobras (Tabela B.3)', dobras, cods, CONFIG_WB_PRINCIPAL,
                                    NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL, ctx['n_jobs'], random_state=s)
        pr['rep'], dg['rep'] = s, s
        preds.append(pr)
        diags.append(dg)
    pcv, dcv = pd.concat(preds, ignore_index=True), pd.concat(diags, ignore_index=True)
    cv = V.resumo_benchmark(pcv, dcv, ctx['sigma'], n_boot=p.n_boot)
    hib = cv[cv['tipo'] == 'hibrido'].sort_values('RMSE')
    escolhido = hib.iloc[0]['modelo']
    cv['selecionado'] = cv['modelo'] == escolhido
    salvar(cv, 'particao_comum_cv.csv', ctx['saida'])

    q = pcv[pcv['tempo_min'] > 0].assign(se=lambda d: (d['X_exp'] - d['x_hat']) ** 2)
    rm = lambda g: 100 * np.sqrt(g['se'].mean())
    por_dobra = (q.groupby(['modelo', 'dobra', 'rep']).apply(rm, include_groups=False).rename('RMSE_pp').reset_index()
                 .groupby(['modelo', 'dobra'])['RMSE_pp'].agg(RMSE_pp='mean', RMSE_dp_sementes_pp='std').reset_index())
    por_dobra['ensaios_validacao'] = por_dobra['dobra'].map(lambda r: ', '.join(map(str, rotulo_dobra[r])))
    salvar(por_dobra, 'particao_comum_cv_dobras.csv', ctx['saida'])
    por_ens = (q.groupby(['modelo', 'ensaio_id', 'dobra', 'rep']).apply(rm, include_groups=False).rename('RMSE_pp')
               .reset_index().groupby(['modelo', 'ensaio_id', 'dobra'])['RMSE_pp'].mean().reset_index())
    salvar(por_ens, 'particao_comum_cv_ensaios.csv', ctx['saida'])
    r = hib.iloc[0]
    print(f"   CV (13 ensaios, 4 dobras, {len(sementes)} sementes): selecionado = {escolhido} "
          f"(RMSE fora da dobra {100 * r['RMSE']:.2f} ± {100 * r['RMSE_dp_reps']:.2f} p.p.)", flush=True)

    # 2. Modelo final: comitê de sementes treinado nos 13 ensaios -----------
    cfg = CONFIG_WB_PRINCIPAL.com_nu(V.escolher_nu(resp, treino, NU_MODO_PRINCIPAL, NU_FIXO_PRINCIPAL))
    acop = V.AcopladorFisico(nu=cfg.nu)
    base = construir_base(df, cfg).sort_values(['ensaio_id', 'tempo_min']).reset_index(drop=True)
    tr = base[base['ensaio_id'].isin(treino)].reset_index(drop=True)
    t_malha = np.round(np.arange(0.0, 15.0 + 1e-9, PASSO_MALHA_MIN), 4)
    malha = []
    for e, g in base.groupby('ensaio_id'):
        malha.append(construir_base_condicao(float(g['razao_molar_eta_real'].iloc[0]), float(g['C_acid_0_mol_L'].iloc[0]),
                                             t_malha, cfg, ensaio_id=int(e), eta_nominal=float(g['razao_molar_eta'].iloc[0])))
    malha = pd.concat(malha, ignore_index=True)
    teto_pts, teto_malha = acop.teto(base), acop.teto(malha)

    tarefas = [(c, s) for c in cods for s in sementes]
    from joblib import Parallel, delayed, parallel_config
    with parallel_config(backend='loky', n_jobs=ctx['n_jobs'], inner_max_num_threads=1):
        saidas = Parallel()(delayed(_ajustar_e_prever)(c, s, tr, [base, malha]) for c, s in tarefas)
    brutos = {}
    for (c, s), (xb_pts, xb_malha) in zip(tarefas, saidas):
        brutos.setdefault(c, {})[s] = (xb_pts, xb_malha)

    # 3. Métricas: treino, teste cego (consultado uma vez) e global ---------
    pos = base['tempo_min'].to_numpy() > 0
    ids = base['ensaio_id'].to_numpy()
    ids_m = malha['ensaio_id'].to_numpy()
    conjuntos = [('Treino', treino), ('Teste cego', teste), ('Global', tuple(sorted(set(treino) | set(teste))))]
    linhas, preds_fin, curvas = [], [], []
    for c in cods:
        nome = nome_por_cod[c]
        xb_pts = np.mean([v[0] for v in brutos[c].values()], axis=0)
        xb_malha = np.mean([v[1] for v in brutos[c].values()], axis=0)
        xh_pts, xh_malha = acop.restringir(base, xb_pts), acop.restringir(malha, xb_malha)
        tipo = V.criar_estimador(c).tipo
        for conj, ens in conjuntos:
            m = pos & np.isin(ids, ens)
            mm = np.isin(ids_m, ens)
            met = _metricas_pp(base.loc[m, 'X_zn_exp'], xh_pts[m])
            err = np.abs(base.loc[m, 'X_zn_exp'].to_numpy() - xh_pts[m])
            k = np.argmax(err)
            linha = {'modelo': nome, 'tipo': tipo, 'conjunto': conj, 'ensaios': len(ens), **met,
                     'erro_max_ensaio': int(ids[m][k]), 'erro_max_t_min': float(base.loc[m, 'tempo_min'].to_numpy()[k]),
                     'inconsistencia_bruta_pct': _pct_inconsistente(malha[mm], xb_malha[mm], teto_malha[mm]),
                     'inconsistencia_final_pct': _pct_inconsistente(malha[mm], xh_malha[mm], teto_malha[mm])}
            if conj == 'Teste cego':
                for e in teste:
                    me = pos & (ids == e)
                    linha[f'RMSE_E{e}_pp'] = 100 * float(np.sqrt(np.mean((base.loc[me, 'X_zn_exp'] - xh_pts[me]) ** 2)))
            linhas.append(linha)
        preds_fin.append(pd.DataFrame({'modelo': nome, 'tipo': tipo, 'ensaio_id': ids,
                                       'conjunto': np.where(np.isin(ids, teste), 'Teste cego', 'Treino'),
                                       'tempo_min': base['tempo_min'], 'X_exp': base['X_zn_exp'],
                                       'X_pbm': base['X_pbm'], 'x_bruto': xb_pts, 'x_hat': xh_pts}))
        if tipo in ('hibrido', 'caixa-preta', 'white-box'):
            curvas.append(pd.DataFrame({'modelo': nome, 'ensaio_id': ids_m, 'tempo_min': malha['tempo_min'],
                                        'x_bruto': xb_malha, 'x_hat': xh_malha}))
    des = pd.DataFrame(linhas)
    des['selecionado'] = des['modelo'] == escolhido
    salvar(des, 'particao_comum_desempenho.csv', ctx['saida'])
    salvar(pd.concat(preds_fin, ignore_index=True), 'particao_comum_predicoes.csv', ctx['saida'])
    salvar(pd.concat(curvas, ignore_index=True), 'particao_comum_curvas.csv', ctx['saida'])

    # Variabilidade entre sementes do modelo selecionado
    c_sel = next(c for c in cods if nome_por_cod[c] == escolhido)
    mt = pos & np.isin(ids, teste)
    sem = []
    for s, (xb, _) in brutos[c_sel].items():
        xh = acop.restringir(base, xb)
        sem.append({'modelo': escolhido, 'semente': s, **{f'{k}_teste': v for k, v in
                                                          _metricas_pp(base.loc[mt, 'X_zn_exp'], xh[mt]).items()}})
    xs = np.array([acop.restringir(base, v[0])[pos] for v in brutos[c_sel].values()])
    faixa = 100 * (xs.max(0) - xs.min(0))
    sem = pd.DataFrame(sem)
    sem['faixa_media_pp_112'] = float(faixa.mean())
    sem['faixa_max_pp_112'] = float(faixa.max())
    salvar(sem, 'particao_comum_sementes.csv', ctx['saida'])

    t = des[(des['modelo'] == escolhido) & (des['conjunto'] == 'Teste cego')].iloc[0]
    print(f"   teste cego ({escolhido}): RMSE {t['RMSE_pp']:.2f} p.p., R² {t['R2']:.3f}, MAE {t['MAE_pp']:.2f} p.p., "
          f"erro máx. {t['erro_max_pp']:.2f} p.p. (E{t['erro_max_ensaio']}, {t['erro_max_t_min']:g} min)", flush=True)
    ctx['particao_comum'] = {'ensaios_teste_cego': list(teste), 'dobras_tabela_b3': [list(d) for d in V.DOBRAS_TABELA_B3],
                             'sementes': list(sementes), 'selecionado_pela_cv': escolhido}


FUNCOES = {'dados': etapa_dados, 'whitebox': etapa_whitebox, 'ablacao': etapa_ablacao, 'benchmark': etapa_benchmark,
           'nu': etapa_nu, 'aprendizado': etapa_aprendizado, 'intervalos': etapa_intervalos,
           'residuo': etapa_residuo, 'dominio': etapa_dominio, 'escada_hibrido': etapa_escada_hibrido,
           'curvas': etapa_curvas, 'particao_comum': etapa_particao_comum}


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
        'particao_comum': ctx.get('particao_comum'),
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
