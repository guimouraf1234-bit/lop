"""Malha do planejamento: quadrados encostados (um por condição), ensaios numerados e E1-E4 com borda vermelha."""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

originais = [(1, .1, .5), (2, .1, 3.1), (3, .5, .5), (4, 1., .5), (5, 1.5, .5), (6, .1, 1.), (7, .5, 3.1),
             (8, .5, 1.), (9, 1., 1.), (10, 1.5, 1.), (11, .1, 1.5), (12, 1., 3.1), (13, .5, 1.5),
             (14, 1., 1.5), (15, 1.5, 1.5), (16, 1.5, 3.1)]
teste = {7, 8, 14}
novos = [('E1', .25, .75), ('E2', .75, 1.25), ('E3', .75, 2.25), ('E4', .25, .40)]
etas = [.4, .5, .75, 1., 1.25, 1.5, 2.25, 3.1]      # colunas
cas = [.1, .25, .5, .75, 1., 1.5]                   # linhas
AZUL, AZUL_CLARO, FUNDO, VERMELHO = '#0072B2', '#56B4E9', '#F4F8FB', '#D62728'
virg = lambda v: (f'{v:.2f}'.rstrip('0').rstrip('.') if v not in (1., 3.1) else f'{v:.1f}').replace('.', ',')

plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 14})
fig, ax = plt.subplots(figsize=(10.5, 6.6))
# malha: quadrados encostados uns nos outros
for i in range(len(etas)):
    for j in range(len(cas)):
        ax.add_patch(Rectangle((i, j), 1, 1, fc=FUNDO, ec='#9DB7CC', lw=1.0, zorder=1))
def celula(eta, ca, fc, txt, cor_txt):
    i, j = etas.index(eta), cas.index(ca)
    ax.add_patch(Rectangle((i, j), 1, 1, fc=fc, ec='#9DB7CC', lw=1.0, zorder=2))
    ax.text(i + .5, j + .5, txt, ha='center', va='center', fontsize=15, fontweight='bold', color=cor_txt, zorder=4)
    return i, j
for n, ca, eta in originais:
    i, j = celula(eta, ca, 'white' if n in teste else AZUL, str(n), 'black' if n in teste else 'white')
    if n in teste:
        ax.add_patch(Rectangle((i + .06, j + .06), .88, .88, fill=False, ec='black', lw=1.6, zorder=3))
for e, ca, eta in novos:
    i, j = celula(eta, ca, AZUL_CLARO, e, 'black')
    ax.add_patch(Rectangle((i + .04, j + .04), .92, .92, fill=False, ec=VERMELHO, lw=3.0, zorder=5))
# contorno do domínio do planejamento original (η 0,5-3,1; C_A0 0,1-1,5)
ax.add_patch(Rectangle((1.02, .02), 6.96, 5.96, fill=False, ec='#1F1F1F', lw=2.0, ls='--', zorder=6))
ax.set_xlim(-.02, 8.02); ax.set_ylim(-.02, 6.02); ax.set_aspect('equal')
ax.set_xticks([i + .5 for i in range(len(etas))]); ax.set_xticklabels([virg(e) for e in etas])
ax.set_yticks([j + .5 for j in range(len(cas))]); ax.set_yticklabels([virg(c) for c in cas])
ax.tick_params(length=0, pad=6)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel(r'Razão molar, $\eta$ (–)'); ax.set_ylabel(r'$C_{\mathrm{A0}}$ (mol/L)')
h = [Line2D([], [], ls='none', marker='s', ms=14, mfc=AZUL, mec='#9DB7CC', label='Treino'),
     Line2D([], [], ls='none', marker='s', ms=14, mfc='white', mec='black', label='Teste cego'),
     Line2D([], [], ls='none', marker='s', ms=14, mfc=AZUL_CLARO, mec='#9DB7CC', label='Ensaios novos'),
     Line2D([], [], ls='none', marker='s', ms=15, mfc='none', mec=VERMELHO, mew=2.4, label='E1 a E4'),
     Line2D([], [], color='#1F1F1F', ls='--', lw=2.0, label='Domínio original')]
fig.legend(handles=h, loc='upper center', ncol=5, frameon=False, fontsize=12.5, bbox_to_anchor=(.53, 1.0))
fig.tight_layout(rect=(0, 0, 1, .94))
for ext in ('png', 'pdf'):
    fig.savefig(f'entregas/malha_ensaios/malha_ensaios.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
