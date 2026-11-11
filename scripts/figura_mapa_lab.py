"""Figura no estilo do mapa de discordância: grade η nominal x C_A0, 16 ensaios, domínio e 6 ensaios propostos."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from data_loader import carregar_respostas_finais

raiz = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
d = pd.read_csv(os.path.join(raiz, "resultados", "mapa_discordancia_lab.csv"))
etas = sorted(d.eta_nominal.unique()); ca0s = sorted(d.C_A0.unique())
M = d.pivot(index="C_A0", columns="eta_nominal", values="discordancia_media").loc[ca0s, etas].values
massa = d.pivot(index="C_A0", columns="eta_nominal", values="massa_g").loc[ca0s, etas].values
r = carregar_respostas_finais()
atuais = {(e, c) for e, c in zip(r.razao_molar_eta, r.C_acid_0_mol_L)}
etas_lbl = [3.1 if e == 3.1 else e for e in etas]

# (η nominal, C_A0)
refinar = {"D1": (0.75, 1.25), "D2": (1.25, 0.25), "D3": (2.5, 1.25)}
validar = {"F1": (0.35, 1.0), "F2": (0.75, 2.0), "F3": (4.0, 0.5)}

fig, ax = plt.subplots(figsize=(10, 6.8), dpi=150)
for i, c in enumerate(ca0s):
    for j, e in enumerate(etas):
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc="#eef2f7", ec="white", lw=2))
        if (e, c) in atuais:
            ax.scatter(j, i, s=200, c="#2a9d4a", ec="k", zorder=4)
ax.set_xlim(-0.5, len(etas) - 0.5); ax.set_ylim(-0.5, len(ca0s) - 0.5)

j0, j1 = etas.index(0.5), etas.index(3.1); i0, i1 = ca0s.index(0.1), ca0s.index(1.5)
ax.add_patch(Rectangle((j0 - 0.5, i0 - 0.5), j1 - j0 + 1, i1 - i0 + 1, fill=False, ec="black", ls="--", lw=1.8, zorder=5))

def caixa(nome, pos, cor):
    j, i = etas.index(pos[0]), ca0s.index(pos[1])
    ax.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False, ec=cor, lw=3, zorder=6))
    ax.text(j - 0.42, i + 0.42, nome, ha="left", va="top", color=cor, fontsize=10, weight="bold", zorder=7,
            bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))


ax.set_xticks(range(len(etas)), [f"{e:g}" for e in etas])
ax.set_yticks(range(len(ca0s)), [f"{c:g}" for c in ca0s])
ax.set_xlabel("Razão molar η (nominal)"); ax.set_ylabel("$C_{A0}$ (mol/L)")
ax.set_title("Malha de ensaios", weight="bold")
ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="#2a9d4a", mec="k", ms=11, label="16 ensaios da tese (40 °C)"),
                   Line2D([], [], color="black", ls="--", lw=1.8, label="Domínio da malha atual")],
          loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, frameon=False, fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(raiz, "resultados", "figuras", "f14_mapa_ensaios_lab.png"))

tab = []
for grupo, dic in (("refinar", refinar), ("validar", validar)):
    for n, (e, c) in dic.items():
        lin = d[(d.eta_nominal == e) & (d.C_A0 == c)].iloc[0]
        tab.append({"id": n, "grupo": grupo, "eta_nominal": e, "eta_real": lin.eta_real, "C_A0_mol_L": c,
                    "massa_solido_g": lin.massa_g, "volume_L": 0.4, "T_C": 40,
                    "discordancia_media": round(lin.discordancia_media, 3),
                    "X_hibrido_15min": round(lin.X_hibrido_medio_15min, 2)})
t = pd.DataFrame(tab)
t.to_csv(os.path.join(raiz, "resultados", "plano_6_ensaios_lab.csv"), index=False)
print(t.to_string(index=False))
