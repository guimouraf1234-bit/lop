"""Planejamento experimental: domínio original (retângulo), um quadrado por ensaio e E1-E4 circulados."""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

originais = [(1, .1, .5), (2, .1, 3.1), (3, .5, .5), (4, 1., .5), (5, 1.5, .5), (6, .1, 1.), (7, .5, 3.1),
             (8, .5, 1.), (9, 1., 1.), (10, 1.5, 1.), (11, .1, 1.5), (12, 1., 3.1), (13, .5, 1.5),
             (14, 1., 1.5), (15, 1.5, 1.5), (16, 1.5, 3.1)]
teste = {7, 8, 14}
novos = [('E1', .25, .75), ('E2', .75, 1.25), ('E3', .75, 2.25), ('E4', .25, .40)]
AZUL, AZUL_CLARO, VERMELHO = '#0072B2', '#56B4E9', '#D62728'

plt.rcParams.update({'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 14})
fig, ax = plt.subplots(figsize=(11, 5.2))
# domínio do planejamento original: retângulo horizontal (η no eixo x)
ax.add_patch(Rectangle((.5, .1), 2.6, 1.4, fc='#EAF4FB', ec=AZUL, lw=1.6, ls='--', zorder=0))
ax.text(1.88, 1.3, 'Domínio do planejamento original', ha='center', va='center', fontsize=12, color=AZUL, style='italic')
S = 700        # área do quadrado em pt² (tamanho fixo na tela, independe da escala dos eixos)
def quadrado(eta, ca, fc, ec, txt, cor_txt):
    ax.scatter(eta, ca, s=S, marker='s', fc=fc, ec=ec, lw=1.4, zorder=3)
    ax.text(eta, ca, txt, ha='center', va='center', fontsize=12.5, fontweight='bold', color=cor_txt, zorder=4)
for n, ca, eta in originais:
    if n in teste:
        quadrado(eta, ca, 'white', 'black', str(n), 'black')
    else:
        quadrado(eta, ca, AZUL, AZUL, str(n), 'white')
for e, ca, eta in novos:
    quadrado(eta, ca, AZUL_CLARO, 'black', e, 'black')
    ax.scatter(eta, ca, s=S * 2.6, marker='o', fc='none', ec=VERMELHO, lw=2.2, zorder=5)
ax.set_xlim(.2, 3.35); ax.set_ylim(-.05, 1.72)
ax.set_xticks([.4, .5, .75, 1., 1.25, 1.5, 2.25, 3.1]); ax.set_yticks([.1, .25, .5, .75, 1., 1.5])
fmt = lambda v, _: (f'{v:.2f}'.rstrip('0').rstrip('.') if v not in (1., 3.1) else f'{v:.1f}').replace('.', ',')
ax.xaxis.set_major_formatter(fmt); ax.yaxis.set_major_formatter(fmt)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.tick_params(direction='out', labelsize=12.5)
ax.set_xlabel(r'Razão molar, $\eta$ (–)'); ax.set_ylabel(r'$C_{\mathrm{A0}}$ (mol/L)')
h = [Line2D([], [], ls='none', marker='s', ms=13, mfc=AZUL, mec=AZUL, label='Treino'),
     Line2D([], [], ls='none', marker='s', ms=13, mfc='white', mec='black', label='Teste cego'),
     Line2D([], [], ls='none', marker='s', ms=13, mfc=AZUL_CLARO, mec='black', label='Ensaios novos'),
     Line2D([], [], ls='none', marker='o', ms=17, mfc='none', mec=VERMELHO, mew=2, label='E1 a E4')]
fig.legend(handles=h, loc='upper center', ncol=4, frameon=False, fontsize=12.5, bbox_to_anchor=(.5, 1.03))
fig.tight_layout(rect=(0, 0, 1, .94))
for ext in ('png', 'pdf'):
    fig.savefig(f'entregas/malha_ensaios/malha_ensaios.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
