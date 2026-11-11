"""
figuras_batelada.py
-------------------
Figuras do TCC (batelada). Lê APENAS os CSV gerados por experimentos_batelada.py (nada é recalculado) e grava
PNG em lop/resultados/figuras/. Cada figura tem o seu "gêmeo em tabela" nos CSV (valores exatos).

    python scripts/figuras_batelada.py [--resultados DIR] [--saida DIR] [--figuras f01,f02,...]

Convenções visuais (método de visualização do projeto; paleta validada com validate_palette.js):
  - cor segue a ENTIDADE em todas as figuras: híbrido = azul, caixa-preta = laranja, White-Box = tinta,
    baselines (tabela, Ridge) = cinza; o modelo (GP, ET, ...) é dado pelo MARCADOR (identidade nunca só por cor,
    para impressão em tons de cinza);
  - atributos A < B < C e C_A0 (ordinais) usam a rampa azul, do claro ao escuro;
  - traços finos, grade em hairline sólida e recessiva, texto em tons de tinta, rótulos diretos seletivos,
    anel da cor da superfície nos marcadores;
  - TÍTULOS QUANTITATIVOS são montados a partir dos dados (nenhum número fixo no código).
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib import ticker                        # noqa: E402
from matplotlib.lines import Line2D                  # noqa: E402

DIR_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
if DIR_SCRIPTS not in sys.path:
    sys.path.insert(0, DIR_SCRIPTS)
from data_loader import carregar_dados_bancada       # noqa: E402
from hybrid_model import CONFIG_WB_PRINCIPAL         # noqa: E402

DIR_RESULTADOS = os.path.abspath(os.path.join(DIR_SCRIPTS, '..', 'resultados'))

# ---- tokens (superfície clara) ----
SUP, TINTA, TINTA2, MUDO = '#fcfcfb', '#0b0b0b', '#52514e', '#898781'
GRADE, EIXO = '#e1e0d9', '#c3c2b7'
AZUL, LARANJA = '#2a78d6', '#eb6834'
RAMPA = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']      # ordinal (passos 250, 400, 550, 700)

MARC = {'gp': 'o', 'extra_trees': 's', 'gbdt': '^', 'random_forest': 'v', 'svr': 'D', 'mlp_tanh': 'P',
        'knn': 'h', 'hist_gbdt': '>', 'adaboost': '<', 'kernel_ridge': 'd', 'poly_ridge': '*', 'gp_rbf': 'p'}
ROT = {'gp': 'GP', 'extra_trees': 'Extra Trees', 'gbdt': 'GBDT', 'random_forest': 'Random Forest',
       'svr': 'SVR', 'mlp_tanh': 'MLP', 'knn': 'k-NN', 'hist_gbdt': 'HistGBDT', 'adaboost': 'AdaBoost',
       'kernel_ridge': 'Kernel Ridge', 'poly_ridge': 'Polinomial', 'gp_rbf': 'GP (RBF)'}
ORDEM_MODELOS = ['White-Box', 'Tabela (η, t)', 'Ridge (resíduo)'] + \
                [f'Híbrido {s} [B]' for s in ('gp', 'extra_trees', 'gbdt', 'random_forest', 'svr', 'mlp_tanh')] + \
                [f'Caixa-preta {s}' for s in ('gp', 'extra_trees', 'gbdt', 'svr')]


class FormatadorVirgula(ticker.ScalarFormatter):
    """Rótulos numéricos com vírgula decimal (texto em português)."""
    def __call__(self, x, pos=None):
        return super().__call__(x, pos).replace('.', ',')


def configurar():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Segoe UI', 'DejaVu Sans'], 'font.size': 9,
        'axes.titlesize': 10, 'axes.titleweight': 'semibold', 'axes.labelsize': 9,
        'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
        'figure.facecolor': SUP, 'axes.facecolor': SUP, 'savefig.facecolor': SUP,
        'text.color': TINTA, 'axes.labelcolor': TINTA2, 'xtick.color': MUDO, 'ytick.color': MUDO,
        'xtick.labelcolor': TINTA2, 'ytick.labelcolor': TINTA2, 'axes.edgecolor': EIXO, 'axes.linewidth': 0.8,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRADE,
        'grid.linewidth': 0.7, 'grid.linestyle': '-', 'axes.axisbelow': True, 'legend.frameon': False,
        'legend.fontsize': 8.5, 'lines.linewidth': 1.6, 'lines.markersize': 6, 'savefig.dpi': 200,
        'axes.unicode_minus': False,
    })


def info_modelo(nome: str) -> dict:
    """Cor (entidade), marcador (modelo) e rótulo curto."""
    if nome == 'White-Box':
        return {'cor': TINTA, 'm': 'o', 'rot': 'White-Box', 'grupo': 'White-Box'}
    if nome.startswith('Tabela'):
        return {'cor': MUDO, 'm': 's', 'rot': 'Tabela (η, t)', 'grupo': 'Baseline'}
    if nome.startswith('Ridge'):
        return {'cor': MUDO, 'm': 'X', 'rot': 'Ridge (resíduo)', 'grupo': 'Baseline'}
    h = re.match(r'Híbrido (\w+) \[(\w)\]', nome)
    if h:
        return {'cor': AZUL, 'm': MARC[h.group(1)], 'rot': f'Híbrido {ROT[h.group(1)]}', 'grupo': 'Híbrido'}
    c = re.match(r'Caixa-preta (\w+)', nome)
    if c:
        return {'cor': LARANJA, 'm': MARC[c.group(1)], 'rot': f'Caixa-preta {ROT[c.group(1)]}', 'grupo': 'Caixa-preta'}
    raise ValueError(nome)


def rotulo_esquema(s: str) -> str:
    """Nome de protocolo para exibição (símbolos corretos, sem contagem de partições)."""
    s = re.sub(r' \(\d+( partições)?\)', '', s)
    return s.replace('<=', '≤').replace('->', '→').replace('C_A0', 'C$_{A0}$')


def pt(ax, x, y, cor, m, ms=6.5, anel=1.3, **kw):
    """Marcador com anel da cor da superfície (legível onde se sobrepõe)."""
    return ax.plot(x, y, ls='none', marker=m, color=cor, ms=ms, markeredgecolor=SUP, markeredgewidth=anel,
                   zorder=kw.pop('zorder', 4), **kw)


def fmt(v, casas=3):
    return f'{v:.{casas}f}'.replace('.', ',')


def ler(R, nome):
    return pd.read_csv(os.path.join(R, nome))


def sigma(R):
    return float(ler(R, 'piso_ruido.csv').set_index('grandeza').loc['sigma (patamares)', 'valor'])


def virgulas(fig):
    for ax in fig.axes:
        for eixo in (ax.xaxis, ax.yaxis):
            if type(eixo.get_major_formatter()) is ticker.ScalarFormatter:
                eixo.set_major_formatter(FormatadorVirgula())


def salvar_fig(fig, saida, nome):
    virgulas(fig)
    caminho = os.path.join(saida, nome)
    fig.savefig(caminho, bbox_inches='tight')
    plt.close(fig)
    print(f'   -> {os.path.relpath(caminho, os.path.dirname(os.path.dirname(saida)))}', flush=True)


def legenda_grupos(fig, y=-0.01, baselines=True, extras=()):
    hs = [Line2D([], [], color=TINTA, marker='o', ls='-', ms=5.5, label='White-Box')]
    if baselines:
        hs.append(Line2D([], [], color=MUDO, marker='s', ls='-', ms=5.5, label='Baselines (tabela, Ridge)'))
    hs += [Line2D([], [], color=AZUL, marker='o', ls='-', ms=5.5, label='Híbrido (White-Box + ML)'),
           Line2D([], [], color=LARANJA, marker='o', ls='-', ms=5.5, label='Caixa-preta (só ML)')] + list(extras)
    fig.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, y), ncol=len(hs), columnspacing=1.6,
               handletextpad=0.5)


# =============================================================================
# F01 -- 16 ensaios, fora-da-dobra (LOGO)
# =============================================================================

def f01_cinetica_logo(R, S):
    dados = carregar_dados_bancada()
    curvas = ler(R, 'curvas_logo.csv')
    pred = pd.read_csv(os.path.join(R, 'predicoes_oof.csv.gz'))
    pl = pred[pred['esquema'].str.startswith('LOGO') & (pred['tempo_min'] > 0)].copy()
    pl['e2'] = (pl['X_exp'] - pl['x_hat']) ** 2
    rm = pl.groupby(['ensaio_id', 'modelo'])['e2'].mean().pow(0.5)
    etas, cas = [0.5, 1.0, 1.5, 3.1], [0.1, 0.5, 1.0, 1.5]
    fig, axes = plt.subplots(4, 4, figsize=(11, 9.2), sharex=True, sharey=True)
    for i, eta in enumerate(etas):
        for j, ca in enumerate(cas):
            ax = axes[i, j]
            ens = int(dados[(dados['razao_molar_eta'] == eta) & (dados['C_acid_0_mol_L'] == ca)]['ensaio_id'].iloc[0])
            c = curvas[curvas['ensaio_id'] == ens]
            e = dados[dados['ensaio_id'] == ens]
            teto = min(1.0, float(e['razao_molar_eta_real'].iloc[0]) / CONFIG_WB_PRINCIPAL.nu)
            ax.plot(c['tempo_min'], c['X_pbm'], color=TINTA, lw=1.4, zorder=2)
            lo = np.clip(c['X_gp'] - 1.96 * c['sd_gp'], 0, teto)
            hi = np.clip(c['X_gp'] + 1.96 * c['sd_gp'], 0, teto)
            ax.fill_between(c['tempo_min'], lo, hi, color=AZUL, alpha=0.14, lw=0, zorder=1)
            ax.plot(c['tempo_min'], c['X_gp'], color=AZUL, lw=1.8, zorder=3)
            ax.plot(e['tempo_min'], e['X_zn_exp'], ls='none', marker='o', ms=5, mfc=SUP, mec=TINTA, mew=1.1, zorder=5)
            r_wb = rm.get((ens, 'White-Box'), np.nan)
            r_gp = rm.get((ens, 'Híbrido gp [B]'), np.nan)
            ax.text(0.97, 0.06, f'RMSE  WB {fmt(r_wb)}\nGP {fmt(r_gp)}', transform=ax.transAxes, ha='right',
                    va='bottom', fontsize=7.4, color=TINTA2, linespacing=1.25)
            ax.set_title(f'η = {eta:g}  ·  $C_{{A0}}$ = {ca:g} mol/L'.replace('.', ','), fontsize=8.8, color=TINTA2,
                         fontweight='normal', pad=3)
            ax.set_xlim(-0.3, 15.5)
            ax.set_ylim(-0.02, 1.06)
            ax.set_xticks([0, 5, 10, 15])
            if j == 0:
                ax.set_ylabel('Conversão $X_{Zn}$')
            if i == 3:
                ax.set_xlabel('Tempo (min)')
    hs = [Line2D([], [], ls='none', marker='o', ms=5, mfc=SUP, mec=TINTA, mew=1.1, label='Experimental'),
          Line2D([], [], color=TINTA, lw=1.4, label='White-Box de laboratório (CAf, ν = 1)'),
          Line2D([], [], color=AZUL, lw=1.8, label='Híbrido GP, treinado SEM o ensaio (LOGO); faixa de 95%')]
    fig.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, 0.975), ncol=3, columnspacing=1.8)
    fig.suptitle('Cinética dos 16 ensaios: cada curva vem de um modelo que não viu aquele ensaio', y=1.0,
                 fontsize=11, fontweight='semibold')
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    salvar_fig(fig, S, 'f01_cinetica_16_ensaios_logo.png')


# =============================================================================
# F02 -- benchmark por protocolo de validação
# =============================================================================

def f02_benchmark(R, S):
    b = ler(R, 'benchmark_batelada.csv')
    sig = sigma(R)
    esquemas = list(dict.fromkeys(b['esquema']))
    fig, axes = plt.subplots(2, 3, figsize=(12.4, 8.2), sharey=True)
    ordem = [m for m in ORDEM_MODELOS if m in set(b['modelo'])]
    y = {m: len(ordem) - 1 - k for k, m in enumerate(ordem)}
    for ax, esq in zip(axes.ravel(), esquemas):
        g = b[b['esquema'] == esq].set_index('modelo')
        wb = float(g.loc['White-Box', 'RMSE'])
        xmax = round(1.7 * wb, 3)
        for m in ordem:
            if m not in g.index:
                continue
            r = g.loc[m]
            inf = r['RMSE_ic95_inf'] if np.isfinite(r['RMSE_ic95_inf']) else r['RMSE']
            sup = r['RMSE_ic95_sup'] if np.isfinite(r['RMSE_ic95_sup']) else r['RMSE']
            st = info_modelo(m)
            ax.plot([inf, min(sup, xmax)], [y[m]] * 2, color=st['cor'], lw=1.3, alpha=0.9, solid_capstyle='round', zorder=2)
            if r['RMSE'] <= xmax - 0.004:
                pt(ax, r['RMSE'], y[m], st['cor'], st['m'])
            else:
                ax.annotate(f"{fmt(r['RMSE'], 2)} →", (xmax, y[m]), ha='right', va='center', fontsize=7.8,
                            color=st['cor'] if st['cor'] != MUDO else TINTA2, xytext=(-2, 0), textcoords='offset points')
        ax.axvline(wb, color=TINTA, lw=0.9, zorder=1)
        ax.axvline(sig, color=MUDO, lw=0.9, zorder=1)
        ax.text(wb, len(ordem) - 0.35, ' White-Box', fontsize=7.6, color=TINTA, ha='left', va='bottom')
        ax.text(sig, len(ordem) - 0.35, 'σ ', fontsize=7.6, color=TINTA2, ha='right', va='bottom')
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.7, len(ordem) - 0.1)
        ax.set_title(rotulo_esquema(esq), fontsize=9.2, pad=12)
        ax.set_yticks([y[m] for m in ordem])
        ax.set_yticklabels([info_modelo(m)['rot'] for m in ordem])
        ax.grid(axis='y', visible=False)
        ax.set_xlabel('RMSE (t > 0), com IC95% por bootstrap de ensaios')
    for ax in axes.ravel()[len(esquemas):]:
        ax.axis('off')
    fig.suptitle('Erro fora-da-dobra por protocolo de validação (σ = piso de ruído experimental)',
                 y=1.0, fontsize=11, fontweight='semibold')
    legenda_grupos(fig, y=-0.025)
    fig.tight_layout(rect=[0, 0.01, 1, 0.975])
    salvar_fig(fig, S, 'f02_benchmark_por_protocolo.png')


# =============================================================================
# F03 -- curva de aprendizado
# =============================================================================

def f03_curva_aprendizado(R, S):
    c = ler(R, 'curva_aprendizado.csv')
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    for m in c['modelo'].unique():
        g = c[c['modelo'] == m].sort_values('n_treino')
        st = info_modelo(m)
        ax.errorbar(g['n_treino'], g['RMSE'], yerr=g['erro_padrao'], color=st['cor'], lw=1.6, elinewidth=0.7,
                    capsize=0, alpha=0.9, zorder=3)
        pt(ax, g['n_treino'], g['RMSE'], st['cor'], st['m'])
    ax.set_yscale('log')
    ax.set_xticks(sorted(c['n_treino'].unique()))
    ax.set_xlabel('Ensaios usados no treino (de 16; 2 ensaios fora para teste)')
    ax.set_ylabel('RMSE nos ensaios fora do treino (t > 0), escala log')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: fmt(v, 2)))
    ax.yaxis.set_minor_formatter(ticker.FuncFormatter(lambda v, _: fmt(v, 2) if round(v, 3) in (0.03, 0.05, 0.07, 0.15) else ''))
    n0, n1 = c['n_treino'].min(), c['n_treino'].max()
    def melhor(prefixo, n):
        g = c[(c['n_treino'] == n) & c['modelo'].str.startswith(prefixo)]
        return g['RMSE'].min()
    ax.set_title(f"Com {n0} ensaios: melhor híbrido {fmt(melhor('Híbrido', n0))} × melhor caixa-preta {fmt(melhor('Caixa', n0))};"
                 f" com {n1}: {fmt(melhor('Híbrido', n1))} × {fmt(melhor('Caixa', n1))}", fontsize=9.6)
    hs = [Line2D([], [], color=TINTA, marker='o', label='White-Box'), Line2D([], [], color=MUDO, marker='s', label='Tabela (η, t)'),
          Line2D([], [], color=AZUL, marker='o', label='Híbrido'), Line2D([], [], color=LARANJA, marker='o', label='Caixa-preta')]
    hs += [Line2D([], [], color=TINTA2, marker=MARC[s], ls='none', label=ROT[s]) for s in ('gp', 'extra_trees', 'gbdt')]
    ax.legend(handles=hs, ncol=2, loc='upper right', columnspacing=1.5)
    fig.tight_layout()
    salvar_fig(fig, S, 'f03_curva_de_aprendizado.png')


# =============================================================================
# F04 / F05 -- escada e validação termodinâmica do White-Box
# =============================================================================

def f04_escada_whitebox(R, S):
    e = ler(R, 'escada_whitebox.csv')
    rot = ['Tese (α = 0)', 'Tese (α = 5500)', '+ RRB mássica', '+ k_s(T)', '+ [H⁺] V0 (K₂ = 0,0099; TCC)',
           "V1: K₂(T) van't Hoff", 'V2: V1 + par iônico ZnSO₄', 'V0 + ν', 'V1 + ν', 'V2 + ν',
           'Laboratório: CAf, ν = 1 (principal)']
    n = len(e)
    yy = np.arange(n)[::-1]
    sig = sigma(R)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 5.2), sharey=True, gridspec_kw={'width_ratios': [1.15, 1]})
    for k in range(n):
        pt(a1, e.loc[k, 'RMSE_t>0'], yy[k], TINTA, 'o')
        a1.text(e.loc[k, 'RMSE_t>0'] + 0.002, yy[k], fmt(e.loc[k, 'RMSE_t>0']), va='center', fontsize=8, color=TINTA2)
    a1.axvline(sig, color=MUDO, lw=0.9)
    a1.text(sig + 0.001, -0.45, 'piso de ruído σ', fontsize=7.8, color=TINTA2, va='bottom')
    a1.set_xlim(0, 0.115)
    a1.set_ylim(-0.7, n - 0.4)
    a1.set_xlabel('RMSE da conversão (t > 0), 112 pontos')
    a1.set_title('Ajuste de $X$ (112 pontos)', fontsize=9.6)
    a1.set_yticks(yy)
    a1.set_yticklabels(rot)
    for lab in a1.get_yticklabels():
        if 'principal' in lab.get_text():
            lab.set_fontweight('semibold')
            lab.set_color(TINTA)
    a1.grid(axis='y', visible=False)
    ok = e['RMSE_pH0'].notna()
    for k in range(n):
        if not ok[k]:
            a2.text(0.03, yy[k], 'sem [H⁺] (força motriz = H₂SO₄)', va='center', fontsize=7.8, color=MUDO)
            continue
        pt(a2, e.loc[k, 'RMSE_pH0'], yy[k], MUDO, 'o')
        pt(a2, e.loc[k, 'RMSE_pHf_espec'], yy[k], TINTA, 'D')
        a2.text(e.loc[k, 'RMSE_pHf_espec'] + 0.04, yy[k], fmt(e.loc[k, 'RMSE_pHf_espec'], 2), va='center', fontsize=8, color=TINTA2)
    a2.set_xlim(0, 1.35)
    a2.set_xlabel('RMSE do pH previsto (unidades de pH)')
    a2.set_title('pH medido (única medida independente do estado ácido)', fontsize=9.6)
    a2.grid(axis='y', visible=False)
    hs = [Line2D([], [], color=MUDO, marker='o', ls='none', label='pH₀ (16 ensaios)'),
          Line2D([], [], color=TINTA, marker='D', ls='none', label='pH final, η ≥ 1 (12 ensaios), estado final medido')]
    a2.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=1)
    fig.tight_layout()
    salvar_fig(fig, S, 'f04_escada_whitebox.png')


def f05_ph_final(R, S):
    t = ler(R, 'termodinamica_por_ensaio.csv')
    t = t[t['pHf_comparavel']]
    fig, ax = plt.subplots(figsize=(5.8, 5.4))
    lim = (-0.6, 3.1)
    ax.plot(lim, lim, color=EIXO, lw=1.0, zorder=1)
    cfg = {'V0': (MUDO, 'o', 'V0 (K₂ = 0,0099)'), 'V1': (TINTA2, 's', 'V1 (K₂(T))'), 'V2': (TINTA, 'D', 'V2 (+ ZnSO₄⁰)')}
    fator = {}
    for esp, (cor, m, rot) in cfg.items():
        g = t[t['especiacao'] == esp]
        pt(ax, g['pHf_med'], g['pHf_espec'], cor, m, label=rot)
        fator[esp] = float((10.0 ** (g['pHf_espec'] - g['pHf_med'])).max())
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect('equal')
    ax.set_xlabel('pH final medido')
    ax.set_ylabel('pH final previsto (especiação sobre o estado final medido)')
    ax.set_title(f"[H⁺] previsto até {fator['V0']:.0f}× menor que o medido sem o par iônico (V0);\n"
                 f"até {fator['V2']:.1f}× com ele (V2)".replace('.', ','), fontsize=9.2)
    ax.legend(loc='upper left')
    fig.tight_layout()
    salvar_fig(fig, S, 'f05_ph_final_previsto_vs_medido.png')


# =============================================================================
# F06 -- ablação dos atributos
# =============================================================================

def f06_ablacao(R, S):
    a = ler(R, 'ablacao_entradas.csv')
    a['slug'] = a['modelo'].str.extract(r'Híbrido (\w+) \[')[0]
    esq = list(dict.fromkeys(a['esquema_ablacao']))
    slugs = ['gp', 'extra_trees', 'gbdt', 'svr']
    cores_ent = RAMPA[:3]
    fig, axes = plt.subplots(1, len(esq), figsize=(6.0 * len(esq), 4.6), sharey=False)
    for ax, e in zip(np.atleast_1d(axes), esq):
        g = a[a['esquema_ablacao'] == e]
        for i, ent in enumerate(['A', 'B', 'C']):
            for k, s in enumerate(slugs):
                r = g[(g['entradas'] == ent) & (g['slug'] == s)]
                if r.empty:
                    continue
                r = r.iloc[0]
                x = k + (i - 1) * 0.22
                ax.plot([x, x], [r['RMSE_ic95_inf'], r['RMSE_ic95_sup']], color=cores_ent[i], lw=1.2)
                pt(ax, x, r['RMSE'], cores_ent[i], MARC[s])
        ax.set_xticks(range(len(slugs)))
        ax.set_xticklabels([ROT[s] for s in slugs])
        ax.set_ylabel('RMSE (t > 0), com IC95%')
        ax.set_title(rotulo_esquema(e), fontsize=9.8)
        ax.grid(axis='x', visible=False)
    hs = [Line2D([], [], color=RAMPA[0], marker='o', ls='-', label='A: t, η, $C_{A0}$ (sem física)'),
          Line2D([], [], color=RAMPA[1], marker='o', ls='-', label='B: A + $X_{PBM}$'),
          Line2D([], [], color=RAMPA[2], marker='o', ls='-', label='C: B + $C_{Af}$, [H⁺] do PBM')]
    fig.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.04), ncol=3)
    fig.suptitle('Atributos do ML: efeito de informar a conversão ($X_{PBM}$) e o estado químico do White-Box',
                 y=1.02, fontsize=10.6, fontweight='semibold')
    fig.tight_layout()
    salvar_fig(fig, S, 'f06_ablacao_atributos.png')


# =============================================================================
# F07 -- extrapolação por nível de η deixado de fora
# =============================================================================

def f07_extrapolacao_eta(R, S):
    b = ler(R, 'benchmark_batelada.csv')
    g = b[b['esquema'].str.startswith('Deixa-um-η')].set_index('modelo')
    sig = sigma(R)
    modelos = ['White-Box', 'Tabela (η, t)', 'Híbrido extra_trees [B]', 'Híbrido gp [B]', 'Caixa-preta extra_trees', 'Caixa-preta gp']
    cols = {0.5: 'RMSE_eta0.5', 1.0: 'RMSE_eta1', 1.5: 'RMSE_eta1.5', 3.1: 'RMSE_eta3.1'}
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    ymax = 0.16
    off = np.linspace(-0.3, 0.3, len(modelos))
    for k, m in enumerate(modelos):
        st = info_modelo(m)
        for i, (eta, col) in enumerate(cols.items()):
            v = g.loc[m, col]
            x = i + off[k]
            if v <= ymax:
                pt(ax, x, v, st['cor'], st['m'])
            else:
                ax.annotate(f'↑ {fmt(v, 2)}', (x, ymax), ha='center', va='top', fontsize=7.4,
                            color=st['cor'] if st['cor'] != MUDO else TINTA2, xytext=(0, -2), textcoords='offset points')
    n_h = sum(min(g.loc[m, c] for m in ('Híbrido extra_trees [B]', 'Híbrido gp [B]')) < g.loc['White-Box', c] for c in cols.values())
    n_c = sum(min(g.loc[m, c] for m in ('Caixa-preta extra_trees', 'Caixa-preta gp')) < g.loc['White-Box', c] for c in cols.values())
    ax.axhline(sig, color=MUDO, lw=0.9)
    ax.text(3.55, sig, 'piso de ruído σ', fontsize=7.8, color=TINTA2, va='bottom', ha='right')
    ax.set_xticks(range(4))
    ax.set_xticklabels([f'η = {e:g}'.replace('.', ',') for e in cols])
    ax.set_xlim(-0.55, 3.55)
    ax.set_ylim(0, ymax)
    ax.set_xlabel('Nível de η deixado FORA do treino (os seus 4 ensaios são o teste)')
    ax.set_ylabel('RMSE no nível ausente (t > 0)')
    ax.grid(axis='x', visible=False)
    ax.set_title(f'Extrapolar em η: o melhor híbrido supera o White-Box em {n_h} de 4 níveis;\n'
                 f'a melhor caixa-preta, em {n_c} de 4', fontsize=9.8)
    legenda_grupos(fig, y=-0.06)
    hs = [Line2D([], [], color=TINTA2, marker=MARC[s], ls='none', label=ROT[s]) for s in ('gp', 'extra_trees')]
    ax.legend(handles=hs, loc='upper left', ncol=2)
    fig.tight_layout()
    salvar_fig(fig, S, 'f07_extrapolacao_por_nivel_de_eta.png')


# =============================================================================
# F08 -- cobertura dos intervalos
# =============================================================================

def rotulo_intervalo(modelo: str, metodo: str) -> str:
    slug = re.match(r'Híbrido (\w+) \[', modelo).group(1)
    return ROT[slug] + '\n' + ('analítico' if 'analít' in metodo else 'conforme')


def f08_cobertura(R, S):
    c = ler(R, 'cobertura_intervalos.csv')
    esq = list(dict.fromkeys(c['esquema']))
    fig, axes = plt.subplots(1, len(esq), figsize=(5.6 * len(esq), 4.6), sharey=True)
    for ax, e in zip(np.atleast_1d(axes), esq):
        g = c[c['esquema'] == e].reset_index(drop=True)
        for k, r in g.iterrows():
            slug = re.match(r'Híbrido (\w+) \[', r['modelo']).group(1)
            ax.plot([k, k], [0, r['cobertura']], color=AZUL, lw=1.4)
            pt(ax, k, r['cobertura'], AZUL, MARC[slug])
            ax.text(k, r['cobertura'] - 0.035, f"{r['cobertura'] * 100:.0f}%", ha='center', va='top', fontsize=8.4, color=TINTA2)
        ax.axhline(0.95, color=TINTA, lw=0.9)
        ax.text(len(g) - 0.5, 0.955, 'nominal 95%', ha='right', va='bottom', fontsize=7.8, color=TINTA)
        ax.set_xticks(range(len(g)))
        ax.set_xticklabels([rotulo_intervalo(m, me) for m, me in zip(g['modelo'], g['metodo'])], fontsize=8)
        ax.set_ylim(0, 1.08)
        ax.set_title(f"{rotulo_esquema(e)}: cobertura de {g['cobertura'].min() * 100:.0f} a {g['cobertura'].max() * 100:.0f}%",
                     fontsize=9.6)
        ax.grid(axis='x', visible=False)
    np.atleast_1d(axes)[0].set_ylabel('Cobertura empírica do intervalo de 95% (t > 0)')
    fig.suptitle('Intervalos de 95% avaliados em pontos fora-da-dobra', y=1.02, fontsize=10.6, fontweight='semibold')
    fig.tight_layout()
    salvar_fig(fig, S, 'f08_cobertura_dos_intervalos.png')


# =============================================================================
# F09 -- o resíduo que o ML precisa explicar
# =============================================================================

def eixo_eta_log(ax):
    ax.set_xscale('log')
    ax.set_xticks([0.5, 1, 2, 3])
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:g}'.replace('.', ',')))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())


def f09_residuo(R, S):
    p = ler(R, 'residuo_patamar_por_ensaio.csv')
    l = ler(R, 'residuo_aprendido.csv')
    niveis = sorted(p['C_acid_0_mol_L'].unique())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.2, 4.6), sharey=True)
    for i, ca in enumerate(niveis):
        g = p[p['C_acid_0_mol_L'] == ca]
        pt(a1, g['razao_molar_eta_real'], g['delta_X'], RAMPA[i], 'o', ms=7.5)
    for ax in (a1, a2):
        ax.axhline(0, color=EIXO, lw=0.9)
        eixo_eta_log(ax)
        ax.set_xlabel('η real (Eq. 4.3), escala log')
    a1.set_ylabel('$X_{exp} - X_{PBM}$ em 15 min')
    pior = p.loc[p['delta_X'].idxmin()]
    a1.set_title(f"O que o White-Box de laboratório não explica: pior em η real = {pior['razao_molar_eta_real']:.2f}".replace('.', ',') +
                 f" ({pior['delta_X']:.2f})".replace('.', ','), fontsize=9.4)
    hs = [Line2D([], [], color=RAMPA[i], marker='o', ls='none', label=f'$C_{{A0}}$ = {ca:g} mol/L'.replace('.', ',')) for i, ca in enumerate(niveis)]
    a1.legend(handles=hs, loc='lower right')
    g = l[(l['ca0'] == 1.0) & (l['tempo_min'] == 15.0)]
    a2.plot(g['eta_real'], g['delta_gp'], color=AZUL, lw=1.8, label='Híbrido GP')
    a2.plot(g['eta_real'], g['delta_et'], color=AZUL, lw=1.2, marker='s', ms=3.5, markevery=3, label='Híbrido Extra Trees')
    pe = p[p['C_acid_0_mol_L'] == 1.0]
    a2.plot(pe['razao_molar_eta_real'], pe['delta_X'], ls='none', marker='o', ms=7.5, mfc=SUP, mec=TINTA, mew=1.2,
            label='Ensaios com $C_{A0}$ = 1,0')
    a2.set_title('O que o ML aprendeu ($C_{A0}$ = 1,0; $t$ = 15 min): informação só nos 4 níveis de η ensaiados', fontsize=9.4)
    a2.legend(loc='lower right')
    fig.tight_layout()
    salvar_fig(fig, S, 'f09_residuo_do_whitebox_e_aprendido.png')


# =============================================================================
# F10 -- domínio de aplicabilidade e incerteza entre níveis de η
# =============================================================================

def f10_dominio(R, S):
    d = ler(R, 'dominio_incerteza_eta.csv')
    p = ler(R, 'residuo_patamar_por_ensaio.csv')
    niveis_eta = sorted(p.groupby('razao_molar_eta')['razao_molar_eta_real'].mean().round(3))
    cas = sorted(d['ca0'].unique())
    tempos = [2.0, 15.0]
    cor = RAMPA[:3]
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.0), sharex=True)
    for i, t in enumerate(tempos):
        for j, (col, tit) in enumerate([('discordancia_max_min', 'Discordância entre GP, ET, GBDT e SVR\n(máx − mín da conversão prevista)'),
                                        ('sd_gp', 'Desvio-padrão preditivo do GP\n(do resíduo)')]):
            ax = axes[i, j]
            for k, ca in enumerate(cas):
                g = d[(d['ca0'] == ca) & (d['tempo_min'] == t)]
                ax.plot(g['eta_real'], g[col], color=cor[k], lw=1.8, label=f'$C_{{A0}}$ = {ca:g}'.replace('.', ','))
            for e in niveis_eta:
                ax.axvline(e, color=EIXO, lw=0.9, zorder=0)
            eixo_eta_log(ax)
            if i == 0:
                ax.set_title(tit, fontsize=9.4)
            if j == 0:
                ax.set_ylabel(f't = {t:g} min')
            if i == 1:
                ax.set_xlabel('η real, escala log (linhas cinza: níveis ensaiados)')
    axes[0, 0].legend(title='$C_{A0}$ não ensaiado (mol/L)', loc='upper right', title_fontsize=8)
    mx = d.loc[d['tempo_min'] == 15.0, 'discordancia_max_min'].max()
    fig.suptitle(f'Entre os níveis ensaiados de η os modelos discordam em até {fmt(mx, 2)} de conversão (t = 15 min)',
                 y=1.0, fontsize=11, fontweight='semibold')
    fig.tight_layout()
    salvar_fig(fig, S, 'f10_dominio_de_aplicabilidade.png')


# =============================================================================
# F11 -- o híbrido sobre degraus do White-Box
# =============================================================================

def f11_escada_hibrido(R, S):
    e = ler(R, 'escada_hibrido.csv')
    degraus = list(dict.fromkeys(e['degrau_whitebox']))
    rot = ['1\nTese\n(α = 0)', '2\nTCC\n(V0, ν = 1)', '3\nV0 + ν', '4\nCAf, ν = 1\n(principal)']
    esq = list(dict.fromkeys(e['esquema']))
    fig, axes = plt.subplots(1, len(esq), figsize=(3.9 * len(esq), 4.6), sharey=False)
    modelos = ['White-Box', 'Híbrido gp [B]', 'Híbrido extra_trees [B]', 'Caixa-preta gp', 'Caixa-preta extra_trees']
    for ax, q in zip(np.atleast_1d(axes), esq):
        g = e[e['esquema'] == q]
        for m in modelos:
            st = info_modelo(m)
            v = [float(g[(g['degrau_whitebox'] == d) & (g['modelo'] == m)]['RMSE'].iloc[0]) for d in degraus]
            ax.plot(range(len(degraus)), v, color=st['cor'], lw=1.6, alpha=0.95)
            pt(ax, range(len(degraus)), v, st['cor'], st['m'])
        ax.set_xticks(range(len(degraus)))
        ax.set_xticklabels(rot, fontsize=7.8)
        ax.set_title(rotulo_esquema(q), fontsize=9)
        ax.grid(axis='x', visible=False)
    np.atleast_1d(axes)[0].set_ylabel('RMSE (t > 0)')
    legenda_grupos(fig, y=-0.07, baselines=False,
                   extras=[Line2D([], [], color=TINTA2, marker=MARC[s], ls='none', label=ROT[s]) for s in ('gp', 'extra_trees')])
    fig.suptitle('O mesmo ML sobre quatro versões do White-Box (a caixa-preta não depende do White-Box)', y=1.03,
                 fontsize=10.6, fontweight='semibold')
    fig.tight_layout()
    salvar_fig(fig, S, 'f11_hibrido_sobre_degraus_do_whitebox.png')


# =============================================================================
# F12 -- paridade fora-da-dobra
# =============================================================================

def f12_paridade(R, S):
    p = pd.read_csv(os.path.join(R, 'predicoes_oof.csv.gz'))
    sig = sigma(R)
    esq = [e for e in p['esquema'].unique() if e.startswith('GKF8 (partição')][0]
    p = p[(p['esquema'] == esq) & (p['tempo_min'] > 0)]
    modelos = ['White-Box', 'Tabela (η, t)', 'Híbrido gp [B]', 'Híbrido extra_trees [B]', 'Caixa-preta gp', 'Caixa-preta extra_trees']
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 7.2), sharex=True, sharey=True)
    xs = np.linspace(0.3, 1.05, 50)
    for ax, m in zip(axes.ravel(), modelos):
        g = p[p['modelo'] == m]
        st = info_modelo(m)
        ax.fill_between(xs, xs - 2 * sig, xs + 2 * sig, color=EIXO, alpha=0.35, lw=0)
        ax.plot(xs, xs, color=EIXO, lw=1.0)
        pt(ax, g['X_exp'], g['x_hat'], st['cor'], st['m'], ms=4.4, anel=0.7, alpha=0.9)
        r = float(np.sqrt(np.mean((g['X_exp'] - g['x_hat']) ** 2)))
        ax.set_title(f"{st['rot']}   RMSE {fmt(r)}", fontsize=9.2)
        ax.set_xlim(0.3, 1.05)
        ax.set_ylim(0.3, 1.05)
        ax.set_aspect('equal')
    for ax in axes[1]:
        ax.set_xlabel('Conversão experimental')
    for ax in axes[:, 0]:
        ax.set_ylabel('Conversão prevista (fora-da-dobra)')
    fig.suptitle('Paridade (GroupKFold do TCC, t > 0); faixa cinza = ±2σ do erro experimental', y=1.0,
                 fontsize=10.6, fontweight='semibold')
    fig.tight_layout()
    salvar_fig(fig, S, 'f12_paridade_fora_da_dobra.png')


# =============================================================================
# F13 -- partição comum do relatório: Ensaio 1 (treino) e teste cego 7, 8, 14
# =============================================================================

def _curvas_teste_cego(R, S, hib, arquivo, pastas=''):
    """Ensaio 1 (treino) e teste cego: híbrido `hib`, a caixa-preta com o mesmo regressor e o White-Box."""
    des = ler(R, 'particao_comum_desempenho.csv')
    cur = ler(R, 'particao_comum_curvas.csv')
    pred = ler(R, 'particao_comum_predicoes.csv')
    caixa = 'Caixa-preta ' + re.match(r'Híbrido (\w+) ', hib).group(1)
    teste = sorted(pred.loc[pred['conjunto'] == 'Teste cego', 'ensaio_id'].unique())
    ensaios = [1] + [int(e) for e in teste]
    df = carregar_dados_bancada().groupby('ensaio_id').first()
    ordem = [m for m in ('White-Box', caixa, hib) if m in set(cur['modelo'])]
    fig, axes = plt.subplots(1, len(ensaios), figsize=(2.75 * len(ensaios), 3.1), sharey=True)
    for ax, e in zip(axes, ensaios):
        for m in ordem:
            g = cur[(cur['modelo'] == m) & (cur['ensaio_id'] == e)].sort_values('tempo_min')
            st = info_modelo(m)
            ax.plot(g['tempo_min'], g['x_hat'], color=st['cor'], lw=1.8 if m == hib else 1.2,
                    ls='--' if m == 'White-Box' else '-', zorder=3 if m == hib else 2)
        g = pred[(pred['modelo'] == hib) & (pred['ensaio_id'] == e)]
        ax.plot(g['tempo_min'], g['X_exp'], ls='none', marker='o', ms=5, color=TINTA, markerfacecolor=SUP,
                markeredgewidth=1.2, zorder=5)
        r = df.loc[e]
        ax.set_title(f"Ensaio {e}{'*' if e in teste else ''}  ($C_{{A0}}$ = {fmt(r['C_acid_0_mol_L'], 2)}; "
                     f"η = {fmt(r['razao_molar_eta'], 1)})", fontsize=9)
        ax.set_xlabel('t (min)')
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel('$X_{Zn}$')
    t = des[(des['modelo'] == hib) & (des['conjunto'] == 'Teste cego')].iloc[0]
    tc = des[(des['modelo'] == caixa) & (des['conjunto'] == 'Teste cego')]
    hs = [Line2D([], [], color=info_modelo(m)['cor'], ls='--' if m == 'White-Box' else '-', label=info_modelo(m)['rot'])
          for m in ordem] + [Line2D([], [], ls='none', marker='o', color=TINTA, markerfacecolor=SUP, label='Experimental')]
    fig.legend(handles=hs, loc='lower center', bbox_to_anchor=(0.5, -0.07), ncol=len(hs))
    extra = f"; caixa-preta {fmt(tc['RMSE_pp'].iloc[0], 2)} p.p." if len(tc) else ''
    fig.suptitle(f"Teste cego (Ensaios {', '.join(map(str, teste))}; * = teste): {info_modelo(hib)['rot']} com "
                 f"RMSE {fmt(t['RMSE_pp'], 2)} p.p. e R² {fmt(t['R2'])}{extra}", y=1.02, fontsize=10.4,
                 fontweight='semibold')
    fig.tight_layout()
    destino = os.path.join(S, pastas) if pastas else S
    os.makedirs(destino, exist_ok=True)
    salvar_fig(fig, destino, arquivo)


def f13_teste_cego(R, S):
    """Regressor selecionado pela CV (figura principal)."""
    des = ler(R, 'particao_comum_desempenho.csv')
    _curvas_teste_cego(R, S, des.loc[des['selecionado'], 'modelo'].iloc[0], 'f13_teste_cego_particao_comum.png')


def f13_por_regressor(R, S):
    """Uma figura igual à f13 para CADA regressor usado como ML do híbrido (pasta figuras/particao_comum/)."""
    des = ler(R, 'particao_comum_desempenho.csv')
    cv = ler(R, 'particao_comum_cv.csv').set_index('modelo')['RMSE']
    hibs = sorted(des.loc[des['tipo'] == 'hibrido', 'modelo'].unique(), key=lambda m: cv[m])
    for k, h in enumerate(hibs, 1):
        slug = re.match(r'Híbrido (\w+) ', h).group(1)
        _curvas_teste_cego(R, S, h, f'f13_{k:02d}_{slug}.png', 'particao_comum')


def f14_resumo_particao_comum(R, S):
    """RMSE da CV (só treino) e do teste cego para cada regressor, como híbrido e como caixa-preta."""
    des = ler(R, 'particao_comum_desempenho.csv')
    cv = ler(R, 'particao_comum_cv.csv').set_index('modelo')
    t = des[des['conjunto'] == 'Teste cego'].set_index('modelo')
    hibs = sorted(t.index[t['tipo'] == 'hibrido'], key=lambda m: cv.loc[m, 'RMSE'])
    slugs = [re.match(r'Híbrido (\w+) ', h).group(1) for h in hibs]
    xmax = 12.0
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 0.42 * len(slugs) + 1.6), sharey=True)
    y = np.arange(len(slugs))[::-1]
    for ax, (tit, col) in zip(axes, [('Validação cruzada (13 ensaios, 4 dobras)', 'cv'),
                                     ('Teste cego (Ensaios 7, 8, 14)', 'teste')]):
        for nome, dy in (('Híbrido', 0.17), ('Caixa-preta', -0.17)):
            vals = []
            for s_ in slugs:
                m = f'Híbrido {s_} [B]' if nome == 'Híbrido' else f'Caixa-preta {s_}'
                if m not in t.index:
                    vals.append(np.nan)
                    continue
                vals.append(100 * cv.loc[m, 'RMSE'] if col == 'cv' else t.loc[m, 'RMSE_pp'])
            cor = AZUL if nome == 'Híbrido' else LARANJA
            ax.barh(y + dy, np.minimum(vals, xmax), height=0.32, color=cor, label=nome)
            for yy, v in zip(y + dy, vals):
                if np.isfinite(v):
                    ax.text(min(v, xmax) + 0.15, yy, fmt(v, 2) + ('  (fora da escala)' if v > xmax else ''),
                            va='center', fontsize=7.2, color=TINTA2)
        wb = 100 * cv.loc['White-Box', 'RMSE'] if col == 'cv' else t.loc['White-Box', 'RMSE_pp']
        ax.axvline(wb, color=TINTA, ls='--', lw=1.0)
        ax.text(wb - 0.1, y.max() + 0.55, f'White-Box {fmt(wb, 1)}', fontsize=7.5, ha='right', color=TINTA2)
        ax.set_xlim(0, xmax + 3)
        ax.set_title(tit, fontsize=9.5)
        ax.set_xlabel('RMSE (p.p.)')
        ax.grid(axis='y', visible=False)
    axes[0].set_yticks(y, [ROT.get(s_, s_) for s_ in slugs])
    axes[0].legend(loc='lower right')
    fig.suptitle('Cada regressor como ML do híbrido (azul) e sozinho (laranja); ordenado pela CV', y=1.0,
                 fontsize=10.4, fontweight='semibold')
    fig.tight_layout()
    salvar_fig(fig, S, 'f14_resumo_regressores_particao_comum.png')


FIGURAS = {'f01': f01_cinetica_logo, 'f02': f02_benchmark, 'f03': f03_curva_aprendizado, 'f04': f04_escada_whitebox,
           'f05': f05_ph_final, 'f06': f06_ablacao, 'f07': f07_extrapolacao_eta, 'f08': f08_cobertura,
           'f09': f09_residuo, 'f10': f10_dominio, 'f11': f11_escada_hibrido, 'f12': f12_paridade,
           'f13': f13_teste_cego, 'f13r': f13_por_regressor, 'f14': f14_resumo_particao_comum}


def main(argv=None):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--resultados', default=DIR_RESULTADOS)
    ap.add_argument('--saida', default=None, help='padrão: <resultados>/figuras')
    ap.add_argument('--figuras', default=','.join(FIGURAS))
    a = ap.parse_args(argv)
    saida = a.saida or os.path.join(a.resultados, 'figuras')
    os.makedirs(saida, exist_ok=True)
    configurar()
    for f in [x.strip() for x in a.figuras.split(',') if x.strip()]:
        FIGURAS[f](a.resultados, saida)


if __name__ == '__main__':
    main()
