import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from data_loader import carregar_respostas_finais

raiz = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
r = carregar_respostas_finais()
dentro = [("D1", 0.80, 1.25), ("D2", 1.30, 0.30), ("D3", 2.40, 1.25)]
fora = [("F1", 0.30, 0.50, ""), ("F2", 1.07, 2.00, ""), ("F3", 1.07, 0.50, "60 °C")]

fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150)
ax.add_patch(plt.Rectangle((0.535, 0.1), 3.29 - 0.535, 1.4, fill=False, ls="--", lw=1.3, ec="#666"))
ax.scatter(r.razao_molar_eta_real, r.C_acid_0_mol_L, s=60, c="#1f4e79", ec="w", label="Malha atual (16 ensaios, 40 °C)", zorder=3)
ax.scatter([p[1] for p in dentro], [p[2] for p in dentro], s=130, marker="D", c="#2a9d4a", ec="k", label="Refinar o modelo (dentro)", zorder=4)
ax.scatter([p[1] for p in fora], [p[2] for p in fora], s=170, marker="*", c="#d62728", ec="k", label="Validar o modelo (fora)", zorder=4)
for n, x, y in dentro:
    ax.annotate(n, (x, y), xytext=(8, 6), textcoords="offset points", fontsize=10, weight="bold")
for n, x, y, t in fora:
    ax.annotate(n + (f" ({t})" if t else ""), (x, y), xytext=(8, 6), textcoords="offset points", fontsize=10, weight="bold")
ax.set_xlim(0, 3.8); ax.set_ylim(0, 2.3)
ax.set_xlabel("η real (H₂SO₄ / Zn)"); ax.set_ylabel("C_A0 (mol/L)")
ax.set_title("Malha atual e 6 ensaios propostos")
ax.legend(loc="upper right", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(raiz, "resultados", "figuras", "f12_plano_simples.png"))
