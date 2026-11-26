"""Malha do planejamento (C_A0 x η): 16 ensaios de Bortot Coelho (2017) e ensaios novos E1-E4."""
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# (ensaio, C_A0 mol/L, η nominal) - Bortot Coelho (2017)
originais = [(1, .1, .5), (2, .1, 3.1), (3, .5, .5), (4, 1., .5), (5, 1.5, .5), (6, .1, 1.), (7, .5, 3.1),
             (8, .5, 1.), (9, 1., 1.), (10, 1.5, 1.), (11, .1, 1.5), (12, 1., 3.1), (13, .5, 1.5),
             (14, 1., 1.5), (15, 1.5, 1.5), (16, 1.5, 3.1)]
teste = {7, 8, 14}
novos = [('E1', .25, .75), ('E2', .75, 1.25), ('E3', .75, 2.25), ('E4', .25, .40)]
COR_TREINO, COR_NOVO = '#0072B2', '#56B4E9'

plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 13, 'axes.linewidth': 1.2})
fig, ax = plt.subplots(figsize=(8, 5.6))
ax.add_patch(plt.Rectangle((.1, .5), 1.4, 2.6, fill=False, ls='--', lw=1.1, ec='#888888', zorder=1))
for n, ca, eta in originais:
    if n in teste:
        ax.scatter(ca, eta, s=150, marker='s', facecolors='white', edgecolors='black', linewidths=1.8, zorder=4)
    else:
        ax.scatter(ca, eta, s=120, marker='o', color=COR_TREINO, edgecolors='black', linewidths=.8, zorder=4)
    ax.annotate(str(n), (ca, eta), xytext=(8, 6), textcoords='offset points', fontsize=11)
for e, ca, eta in novos:
    ax.scatter(ca, eta, s=190, marker='*', color=COR_NOVO, edgecolors='black', linewidths=.8, zorder=5)
    ax.annotate(e, (ca, eta), xytext=(9, -4), textcoords='offset points', fontsize=12, fontweight='bold')
ax.set_xlim(0, 1.65); ax.set_ylim(0, 3.4)
ax.set_xticks([0, .1, .25, .5, .75, 1., 1.5]); ax.set_yticks([0, .4, .5, .75, 1., 1.25, 1.5, 2.25, 3.1])
fmt = lambda v, _: (f'{v:g}').replace('.', ',')
ax.xaxis.set_major_formatter(fmt); ax.yaxis.set_major_formatter(fmt)
ax.set_xlabel(r'Concentração inicial de ácido, $C_{\mathrm{A0}}$ (mol/L)')
ax.set_ylabel(r'Razão molar, $\eta$ (–)')
ax.grid(color='#bfbfbf', lw=.8, ls=':'); ax.tick_params(direction='in', top=True, right=True)
h = [Line2D([], [], ls='none', marker='o', ms=10, mfc=COR_TREINO, mec='black', label='Treino (13 ensaios)'),
     Line2D([], [], ls='none', marker='s', ms=10, mfc='white', mec='black', mew=1.8, label='Teste cego (7, 8 e 14)'),
     Line2D([], [], ls='none', marker='*', ms=14, mfc=COR_NOVO, mec='black', label='Ensaios novos (E1 a E4)'),
     Line2D([], [], ls='--', color='#888888', label='Domínio do planejamento original')]
fig.legend(handles=h, loc='upper center', ncol=2, frameon=False, fontsize=11.5, bbox_to_anchor=(.53, 1.06))
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'entregas/malha_ensaios/malha_ensaios.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
