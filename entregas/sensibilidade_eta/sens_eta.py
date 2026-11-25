"""Sensibilidade a η: X em t = 15 min do híbrido paralelo (GP + acoplador), η de 0 a 5 em passos de 0,1."""
import warnings; warnings.filterwarnings('ignore')
import sys; sys.path.insert(0, 'scripts')  # rodar a partir da raiz, na main
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from joblib import load
import gerar_resultados_gp_relatorio as G
from data_loader import carregar_dados_bancada
from hybrid_model import construir_base, construir_base_condicao, CONFIG_WB_PRINCIPAL as CFG

modelos = load(G.OUT / 'modelos_gp_13_ensaios.joblib')
FATOR_REAL = 1.069382            # η real / η nominal dos ensaios originais
tempos = np.arange(0, 15.001, .25)
etas = np.round(np.arange(0, 5.0001, .1), 1)
cas = [0.1, 0.5, 1.0, 1.5]
linhas = []
for ca0 in cas:
    malhas = []
    for k, eta in enumerate(etas):
        if eta == 0:
            continue
        malhas.append(construir_base_condicao(eta * FATOR_REAL, ca0, tempos, CFG,
                                              ensaio_id=f'{ca0}_{k}', eta_nominal=eta))
    M = pd.concat(malhas, ignore_index=True)
    M['x_bruto'] = np.mean([m.prever_bruto(M) for m in modelos], axis=0)
    M['x_hat'] = G.ACOP.restringir(M, M.x_bruto)
    f = M[np.isclose(M.tempo_min, 15)]
    linhas.append(pd.DataFrame({'C_A0': ca0, 'eta': 0.0, 'X_pbm_15': 0.0, 'X_hibrido_15': 0.0}, index=[0]))
    linhas.append(pd.DataFrame({'C_A0': ca0, 'eta': f.razao_molar_eta.values,
                                'X_pbm_15': f.X_pbm.values, 'X_hibrido_15': f.x_hat.values}))
R = pd.concat(linhas, ignore_index=True)
R.to_csv('entregas/sensibilidade_eta/sensibilidade_eta_hibrido_paralelo.csv', index=False)

base = construir_base(carregar_dados_bancada(), CFG)
exp = base[np.isclose(base.tempo_min, 15)]

plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 13,
                     'axes.linewidth': 1.2})
COR = G.COR_GP
estilos = {0.1: dict(ls=':', marker='v', alpha=.45), 0.5: dict(ls='--', marker='s', alpha=.6),
           1.0: dict(ls='-.', marker='D', alpha=.8), 1.5: dict(ls='-', marker='o', alpha=1.0)}
fig, ax = plt.subplots(figsize=(10, 5.2))
ax.axvspan(0.5, 3.1, color='#ececec', zorder=0)
ax.plot([0, 1, 5], [0, 1, 1], ':', color='#555555', lw=1.2, zorder=1)
for ca0 in cas:
    g = R[R.C_A0 == ca0]
    e = estilos[ca0]
    ax.plot(g.eta, g.X_hibrido_15, color=COR, lw=2.4, ls=e['ls'], alpha=e['alpha'], zorder=3)
    gm = g[np.isclose((g.eta * 2) % 1, 0)]       # marcadores a cada 0,5
    ax.plot(gm.eta, gm.X_hibrido_15, ls='none', marker=e['marker'], ms=6, color=COR, alpha=e['alpha'], zorder=3)
mk = {0.1: 'v', 0.5: 's', 1.0: 'D', 1.5: 'o'}
for ca0 in cas:
    q = exp[np.isclose(exp.C_acid_0_mol_L, ca0)]
    ax.scatter(q.razao_molar_eta, q.X_zn_exp, s=60, marker='o' if ca0 == 1.5 else mk[ca0],
               facecolors='white', edgecolors='black', linewidths=1.4, zorder=5)
ax.set_xlim(0, 5); ax.set_ylim(0, 1.05)
ax.set_xticks(np.arange(0, 5.01, .5)); ax.set_yticks(np.arange(0, 1.01, .2))
fmt = lambda v, _: f'{v:.1f}'.replace('.', ',')
ax.xaxis.set_major_formatter(fmt); ax.yaxis.set_major_formatter(fmt)
ax.set_xlabel(r'$\eta$ (–)'); ax.set_ylabel(r'$X_{\mathrm{Zn}}$ em $t$ = 15 min (–)')
ax.set_title('Híbrido paralelo')
ax.grid(color='#bfbfbf', lw=.8, ls=':')
ax.tick_params(direction='in', top=True, right=True)
h = [Line2D([], [], color=COR, ls=estilos[c]['ls'], marker=estilos[c]['marker'], alpha=estilos[c]['alpha'],
            lw=2.4, label=fr'$C_{{\mathrm{{A0}}}}$ = {str(c).replace(".", ",")} mol/L') for c in cas]
h += [Line2D([], [], ls='none', marker='o', mfc='white', mec='black', mew=1.4, ms=7, label='Experimental (16 ensaios)'),
      Line2D([], [], color='#555555', ls=':', lw=1.2, label=r'min(1, $\eta$)'),
      Patch(color='#ececec', label=r'Faixa de $\eta$ do treino')]
fig.legend(handles=h, loc='upper center', ncol=4, frameon=False, fontsize=11, bbox_to_anchor=(.5, 1.07))
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'entregas/sensibilidade_eta/sensibilidade_eta_hibrido_paralelo.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
print(R.pivot(index='eta', columns='C_A0', values='X_hibrido_15').iloc[::5].round(3).to_string())
