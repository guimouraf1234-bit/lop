"""Gráficos do híbrido paralelo (GP + acoplador) na paleta do relatório, com malha temporal contendo t = 1e-6 min.

Rodar na raiz de um checkout da main (usa scripts/ e entregas/resultados_gp/modelos_gp_13_ensaios.joblib).
Saídas: entregas/graficos_paralelo/ (figuras, tabelas 5.4.1 e 5.4.2 em CSV, valores da sensibilidade).
"""
import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, 'scripts')
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from joblib import load
import gerar_resultados_gp_relatorio as G
from data_loader import carregar_dados_bancada
from hybrid_model import construir_base, construir_base_condicao, CONFIG_WB_PRINCIPAL as CFG
import validacao as V

OUT = Path('entregas/graficos_paralelo')
MODELOS = load(G.OUT / 'modelos_gp_13_ensaios.joblib')       # 10 GPs (sementes 42-51), 13 ensaios de treino
TEMPOS = np.r_[0.0, 1e-6, np.arange(.25, 15.001, .25)]       # malha com o ponto t = 1e-6 min
COR_WB, COR_GP, COR_LIM = G.COR_WB, G.COR_GP, '#555555'
plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 13, 'axes.linewidth': 1.1})
virg = lambda v, _=None: f'{v:.1f}'.replace('.', ',')


def hibrido(malha):
    malha = malha.copy()
    malha['x_bruto'] = np.mean([m.prever_bruto(malha) for m in MODELOS], axis=0)
    malha['x_hat'] = G.ACOP.restringir(malha, malha.x_bruto)
    return malha


def metricas(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float)
    e = y - p
    return {'R2': 1 - (e ** 2).sum() / ((y - y.mean()) ** 2).sum(), 'RMSE': np.sqrt((e ** 2).mean()),
            'MAE': np.abs(e).mean(), 'Erro_max': np.abs(e).max()}


# ---------------------------------------------------------------- 16 ensaios originais (Tabela 5.4.1)
b = construir_base(carregar_dados_bancada(), CFG).sort_values(['ensaio_id', 'tempo_min']).reset_index(drop=True)
malha = pd.concat([construir_base_condicao(float(g.razao_molar_eta_real.iloc[0]), float(g.C_acid_0_mol_L.iloc[0]),
                                           TEMPOS, CFG, ensaio_id=int(e), eta_nominal=float(g.razao_molar_eta.iloc[0]))
                   for e, g in b.groupby('ensaio_id')], ignore_index=True)
curvas = hibrido(malha)
b['x_hat'] = G.pontos_da_malha(b, curvas, curvas.x_hat)
treino = V.ensaios_treino_comum(b)
t541 = []
for nome, ids in [('Treino', treino), ('Teste cego', V.ENSAIOS_TESTE_CEGO), ('Global', range(1, 17))]:
    p = b[b.ensaio_id.isin(ids) & (b.tempo_min > 0)]
    t541.append({'Conjunto': nome, 'Ensaios': len(list(ids)), 'Pontos': len(p), **metricas(p.X_zn_exp, p.x_hat)})
t541 = pd.DataFrame(t541); t541.to_csv(OUT / 'tabela_5_4_1.csv', index=False)

# ---------------------------------------------------------------- ensaios novos E1-E4 (Tabela 5.4.2 e figura)
med = pd.read_csv(G.OUT / 'predicoes_ensaios_novos_gp.csv')
cond = med.groupby('ensaio_id').first()
mn = pd.concat([construir_base_condicao(float(r.eta_ZnO), float(r.C_acid_0_mol_L), TEMPOS, CFG,
                                        ensaio_id=e, eta_nominal=float(r.razao_molar_eta)) for e, r in cond.iterrows()],
               ignore_index=True)
cn = hibrido(mn)
med['x_hat'] = G.pontos_da_malha(med, cn, cn.x_hat)
med['X_pbm'] = G.pontos_da_malha(med, cn, cn.X_pbm)
t542 = []
for e, g in med.groupby('ensaio_id'):
    t542.append({'Ensaio': e, 'C_A0': g.C_acid_0_mol_L.iloc[0], 'eta': g.razao_molar_eta.iloc[0],
                 **metricas(g.X_relatorio, g.x_hat), 'RMSE_white_box': metricas(g.X_relatorio, g.X_pbm)['RMSE']})
t542 = pd.DataFrame(t542); t542.to_csv(OUT / 'tabela_5_4_2.csv', index=False)
cn.to_csv(OUT / 'curvas_ensaios_novos.csv', index=False)

TRACO_GP = (0, (6, 3))
LEG = [Line2D([], [], color=COR_WB, ls=':', lw=2.4, marker='^', ms=8, label='White-box (novo)'),
       Line2D([], [], color=COR_GP, ls=TRACO_GP, lw=2.8, marker='P', ms=9, label='Híbrido paralelo (GP)'),
       Line2D([], [], color=COR_LIM, ls=':', lw=1.2, label=r'Limite estequiométrico, min(1, $\eta$)'),
       Line2D([], [], ls='none', marker='o', mfc='white', mec='black', mew=1.3, ms=7, label='Experimental')]


def painel(ax, g, tp, yp, titulo, rmse_wb, rmse_gp, rotulo_y=True):
    """g: curva (tempo_min, X_pbm, x_hat); tp, yp: medições incluindo t = 0."""
    eta, ca = g.razao_molar_eta.iloc[0], g.C_acid_0_mol_L.iloc[0]
    ax.plot(g.tempo_min, g.X_pbm, ':', color=COR_WB, lw=2.4)
    ax.plot(g.tempo_min, g.x_hat, ls=TRACO_GP, color=COR_GP, lw=2.8)
    ax.axhline(min(1, eta), color=COR_LIM, ls=':', lw=1.2)
    ax.scatter(tp, yp, s=50, fc='white', ec='black', lw=1.3, zorder=6, clip_on=False)
    ax.text(.03, .96, titulo, transform=ax.transAxes, fontweight='bold', va='top', fontsize=15)
    txt = (f'$C_{{A0}}$ = {ca:.2f} mol/L, $\\eta$ = {eta:.2f}\nRMSE: white-box {rmse_wb:.3f}; GP {rmse_gp:.3f}').replace('.', ',')
    ax.text(.97, .05, txt, transform=ax.transAxes, ha='right', va='bottom', fontsize=12, linespacing=1.6)
    ax.set(xlim=(-.3, 15.5), ylim=(0, 1.12), xticks=[0, 5, 10, 15], yticks=np.arange(0, 1.01, .2))
    ax.yaxis.set_major_formatter(virg); ax.grid(color='#bfbfbf', lw=.8, ls=':')
    ax.tick_params(direction='in', top=True, right=True)


def rmse(y, p):
    return float(np.sqrt(np.mean((np.asarray(y) - np.asarray(p)) ** 2)))


# ensaios novos sem o E1, em pirâmide (E2 em cima; E3 e E4 embaixo)
fig = plt.figure(figsize=(10, 7.8))
gs = fig.add_gridspec(2, 4)
eixos = {'E2': fig.add_subplot(gs[0, 1:3]), 'E3': fig.add_subplot(gs[1, 0:2]), 'E4': fig.add_subplot(gs[1, 2:4])}
for e, ax in eixos.items():
    g = cn[cn.ensaio_id == e]
    m = med[med.ensaio_id == e]
    painel(ax, g, np.r_[0, m.tempo_min], np.r_[0, m.X_relatorio], f'Ensaio {e}',
           rmse(m.X_relatorio, m.X_pbm), rmse(m.X_relatorio, m.x_hat))
fig.supxlabel('Tempo de reação, $t$ (min)', fontsize=15); fig.supylabel(r'Conversão de zinco, $X_{Zn}$ (–)', fontsize=15)
fig.legend(handles=LEG, loc='upper center', ncol=4, frameon=False, fontsize=12.5, bbox_to_anchor=(.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, .95))
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'ensaios_novos_paralelo.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)

# teste cego: Ensaios 7, 8 e 14 (7 em cima, centralizado; 8 e 14 embaixo)
fig = plt.figure(figsize=(10, 7.8))
gs = fig.add_gridspec(2, 4)
eixos = {7: fig.add_subplot(gs[0, 1:3]), 8: fig.add_subplot(gs[1, 0:2]), 14: fig.add_subplot(gs[1, 2:4])}
for e, ax in eixos.items():
    g = curvas[curvas.ensaio_id == e]
    m = b[(b.ensaio_id == e)]
    mp = m[m.tempo_min > 0]
    painel(ax, g, m.tempo_min, m.X_zn_exp, f'Ensaio {e}*', rmse(mp.X_zn_exp, mp.X_pbm), rmse(mp.X_zn_exp, mp.x_hat))
fig.supxlabel('Tempo de reação, $t$ (min)', fontsize=15); fig.supylabel(r'Conversão de zinco, $X_{Zn}$ (–)', fontsize=15)
fig.legend(handles=LEG, loc='upper center', ncol=4, frameon=False, fontsize=12.5, bbox_to_anchor=(.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, .95))
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'teste_cego_paralelo.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)

# ---------------------------------------------------------------- sensibilidade a η (t = 15 e 0,5 min)
FATOR_REAL = 1.069382        # η real / η nominal nos ensaios originais
etas = np.round(np.arange(.1, 5.0001, .1), 1)
cas = [.1, .5, 1., 1.5]
ms = pd.concat([construir_base_condicao(eta * FATOR_REAL, ca, TEMPOS, CFG, ensaio_id=f'{ca}_{eta}', eta_nominal=eta)
                for ca in cas for eta in etas], ignore_index=True)
cs = hibrido(ms)
ESTILO = {.1: dict(ls=':', marker='v', alpha=.40), .5: dict(ls='--', marker='s', alpha=.60),
          1.: dict(ls='-.', marker='D', alpha=.80), 1.5: dict(ls='-', marker='o', alpha=1.)}
sens = []
for t_alvo, rot in [(15., '15'), (.5, '0,5')]:
    s = cs[np.isclose(cs.tempo_min, t_alvo)][['C_acid_0_mol_L', 'razao_molar_eta', 'X_pbm', 'x_hat']]
    s = pd.concat([pd.DataFrame({'C_acid_0_mol_L': cas, 'razao_molar_eta': 0., 'X_pbm': 0., 'x_hat': 0.}), s])
    s = s.sort_values(['C_acid_0_mol_L', 'razao_molar_eta']); s.insert(0, 'tempo_min', t_alvo); sens.append(s)
    exp = b[np.isclose(b.tempo_min, t_alvo)]
    fig, ax = plt.subplots(figsize=(10, 5.2))
    ax.axvspan(.5, 3.1, color='#ececec', zorder=0)
    ax.plot([0, 1, 5], [0, 1, 1], ':', color=COR_LIM, lw=1.2, zorder=1)
    for ca in cas:
        g, st = s[s.C_acid_0_mol_L == ca], ESTILO[ca]
        ax.plot(g.razao_molar_eta, g.x_hat, color=COR_GP, lw=2.4, ls=st['ls'], alpha=st['alpha'], zorder=3, clip_on=False)
        gm = g[np.isclose((g.razao_molar_eta * 2) % 1, 0) | np.isclose((g.razao_molar_eta * 2) % 1, 1)]
        ax.plot(gm.razao_molar_eta, gm.x_hat, ls='none', marker=st['marker'], ms=6, color=COR_GP, alpha=st['alpha'], zorder=6, clip_on=False)
        q = exp[np.isclose(exp.C_acid_0_mol_L, ca)]
        ax.scatter(q.razao_molar_eta, q.X_zn_exp, s=70, marker=st['marker'], fc='white', ec='black', lw=1.5, zorder=5)
    ax.set(xlim=(0, 5), ylim=(-.03, 1.05), xticks=np.arange(0, 5.01, .5), yticks=np.arange(0, 1.01, .2))
    ax.xaxis.set_major_formatter(virg); ax.yaxis.set_major_formatter(virg)
    ax.set_xlabel('η (–)'); ax.set_ylabel(fr'$X_{{\mathrm{{Zn}}}}$ em $t$ = {rot} min (–)')
    ax.grid(color='#bfbfbf', lw=.8, ls=':'); ax.tick_params(direction='in', top=True, right=True)
    h = []
    for ca in cas:
        st = ESTILO[ca]
        h.append((Line2D([], [], color=COR_GP, ls=st['ls'], marker=st['marker'], alpha=st['alpha'], lw=2.4),
                  Line2D([], [], ls='none', marker=st['marker'], mfc='white', mec='black', mew=1.5, ms=8)))
    lab = [fr'$C_{{\mathrm{{A0}}}}$ = {virg(c)} mol/L' for c in cas]
    from matplotlib.legend_handler import HandlerTuple
    vazios = tuple(Line2D([], [], ls='none', marker=ESTILO[c]['marker'], mfc='white', mec='black', mew=1.5, ms=8) for c in cas)
    extra = [vazios, Line2D([], [], color=COR_LIM, ls=':', lw=1.2), Patch(color='#ececec')]
    ordem_h = [h[0], h[1], h[2], h[3], extra[0], extra[1], extra[2]]
    ordem_l = lab + ['Medições (símbolos vazios)', r'min(1, η)', 'Faixa de η do treino']
    fig.tight_layout(rect=(0, 0, 1, .86))
    fig.legend(ordem_h, ordem_l, handler_map={tuple: HandlerTuple(ndivide=None, pad=.6)},
               loc='upper center', ncol=4, frameon=False, fontsize=12, bbox_to_anchor=(.5, 1.0), columnspacing=1.4)
    nome = 'sensibilidade_eta_15min' if t_alvo == 15 else 'sensibilidade_eta_0,5min'
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{nome}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
pd.concat(sens).to_csv(OUT / 'sensibilidade_eta.csv', index=False)

pd.set_option('display.float_format', lambda v: f'{v:.3f}')
print(t541.to_string(index=False)); print(); print(t542.to_string(index=False))
