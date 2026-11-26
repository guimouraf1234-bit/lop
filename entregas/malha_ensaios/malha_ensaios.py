"""Malha do planejamento: cada quadrado é um experimento (16 de Bortot Coelho, 2017, e E1-E4)."""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

originais = [(1, .1, .5), (2, .1, 3.1), (3, .5, .5), (4, 1., .5), (5, 1.5, .5), (6, .1, 1.), (7, .5, 3.1),
             (8, .5, 1.), (9, 1., 1.), (10, 1.5, 1.), (11, .1, 1.5), (12, 1., 3.1), (13, .5, 1.5),
             (14, 1., 1.5), (15, 1.5, 1.5), (16, 1.5, 3.1)]
teste = {7, 8, 14}
novos = [('E1', .25, .75), ('E2', .75, 1.25), ('E3', .75, 2.25), ('E4', .25, .40)]
cas = [.1, .25, .5, .75, 1., 1.5]
etas = [.4, .5, .75, 1., 1.25, 1.5, 2.25, 3.1]
COR_TREINO, COR_TESTE, COR_NOVO, COR_VAZIO = '#0072B2', '#FFFFFF', '#56B4E9', '#F2F2F2'
virg = lambda v: f'{v:.2f}'.rstrip('0').rstrip('.').replace('.', ',') if v not in (1., 3.1) else f'{v:.1f}'.replace('.', ',')

plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 13})
fig, ax = plt.subplots(figsize=(10, 6.4))
for i in range(len(etas)):
    for j in range(len(cas)):
        ax.add_patch(Rectangle((i, j), 1, 1, fc=COR_VAZIO, ec='white', lw=3))
def quadrado(ca, eta, cor, rotulo, cor_txt, borda='black', lw=1.6):
    i, j = etas.index(eta), cas.index(ca)
    ax.add_patch(Rectangle((i + .04, j + .04), .92, .92, fc=cor, ec=borda, lw=lw))
    ax.text(i + .5, j + .5, rotulo, ha='center', va='center', fontsize=15, fontweight='bold', color=cor_txt)
for n, ca, eta in originais:
    if n in teste:
        quadrado(ca, eta, COR_TESTE, str(n), 'black')
    else:
        quadrado(ca, eta, COR_TREINO, str(n), 'white')
for e, ca, eta in novos:
    quadrado(ca, eta, COR_NOVO, e, 'black', borda='#D62728', lw=3.2)
ax.set_xlim(0, len(etas)); ax.set_ylim(0, len(cas))
ax.set_xticks([i + .5 for i in range(len(etas))]); ax.set_xticklabels([virg(e) for e in etas])
ax.set_yticks([j + .5 for j in range(len(cas))]); ax.set_yticklabels([virg(c) for c in cas])
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlabel(r'Razão molar, $\eta$ (–)')
ax.set_ylabel(r'$C_{\mathrm{A0}}$ (mol/L)')
ax.set_aspect('equal')
h = [Patch(fc=COR_TREINO, ec='black', label='Treino (13 ensaios)'),
     Patch(fc=COR_TESTE, ec='black', label='Teste cego (7, 8 e 14)'),
     Patch(fc=COR_NOVO, ec='#D62728', lw=2.5, label='Ensaios novos (E1 a E4)'),
     Patch(fc=COR_VAZIO, ec='#cccccc', label='Condição não ensaiada')]
fig.legend(handles=h, loc='upper center', ncol=4, frameon=False, fontsize=12, bbox_to_anchor=(.5, 1.0))
fig.tight_layout(rect=(0, 0, 1, .95))
for ext in ('png', 'pdf'):
    fig.savefig(f'entregas/malha_ensaios/malha_ensaios.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
