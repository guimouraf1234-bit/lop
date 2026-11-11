"""Plano de ensaios em blocos: grade categórica (C_A0 x η) com a malha atual e os 6 ensaios propostos."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from data_loader import carregar_respostas_finais

raiz = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
r = carregar_respostas_finais()

etas = [0.3, 0.53, 0.8, 1.07, 1.3, 1.6, 2.4, 3.2]          # colunas
ca0s = [0.1, 0.3, 0.5, 1.0, 1.25, 1.5, 2.0]               # linhas (de baixo para cima)
nivel_eta = {0.5: 0.53, 1.0: 1.07, 1.5: 1.6, 3.1: 3.2}    # nominal -> coluna
dentro = {(0.8, 1.25): "D1", (1.3, 0.3): "D2", (2.4, 1.25): "D3"}
fora = {(0.3, 0.5): "F1", (1.07, 2.0): "F2"}
f3 = (1.07, 0.5)  # mesma posição do ensaio 8, a 60 °C

AZUL, VERDE, VERMELHO, VAZIO = "#1f4e79", "#2a9d4a", "#d62728", "#eef0f3"
fig, ax = plt.subplots(figsize=(10, 7), dpi=150)
atuais = {(nivel_eta[e], c): int(i) for e, c, i in zip(r.razao_molar_eta, r.C_acid_0_mol_L, r.ensaio_id)}

def bloco(x, y, w, cor, texto, tam=12, h=1):
    ax.add_patch(Rectangle((x + 0.05, y + 0.05), w - 0.1, h - 0.1, fc=cor, ec="white", lw=2))
    ax.text(x + w / 2, y + h / 2, texto, ha="center", va="center", color="white" if cor != VAZIO else "#999",
            fontsize=tam, weight="bold")

for j, e in enumerate(etas):
    for i, c in enumerate(ca0s):
        if (e, c) == f3:
            bloco(j, i + 0.5, 1, AZUL, str(atuais[(e, c)]), 10, h=0.5)
            bloco(j, i, 1, VERMELHO, "F3 60 °C", 7.5, h=0.5)
        elif (e, c) in atuais:
            bloco(j, i, 1, AZUL, str(atuais[(e, c)]))
        elif (e, c) in dentro:
            bloco(j, i, 1, VERDE, dentro[(e, c)], 13)
        elif (e, c) in fora:
            bloco(j, i, 1, VERMELHO, fora[(e, c)], 13)
        else:
            bloco(j, i, 1, VAZIO, "")

# contorno do domínio da malha atual (η 0,53 a 3,2; C_A0 0,1 a 1,5)
ax.add_patch(Rectangle((1, 0), 7, 6, fill=False, ec="black", lw=2.2, ls="--", clip_on=False))

ax.set_xticks([j + 0.5 for j in range(len(etas))], [f"{e:g}".replace(".", ",") for e in etas])
ax.set_yticks([i + 0.5 for i in range(len(ca0s))], [f"{c:g}".replace(".", ",") for c in ca0s])
ax.set_xlabel("η real (H₂SO₄ / Zn)"); ax.set_ylabel("C_A0 (mol/L)")
ax.set_xlim(0, len(etas)); ax.set_ylim(0, len(ca0s))
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("Malha atual e 6 ensaios propostos (eixos fora de escala)", pad=12)
ax.legend(handles=[Patch(fc=AZUL, label="Malha atual (nº do ensaio)"),
                   Patch(fc=VERDE, label="Refinar o modelo (dentro)"),
                   Patch(fc=VERMELHO, label="Validar o modelo (fora)"),
                   Line2D([], [], color="black", lw=2, ls="--", label="Domínio atual (40 °C)")],
          loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, frameon=False, fontsize=10)
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(os.path.join(raiz, "resultados", "figuras", "f13_plano_blocos.png"))
