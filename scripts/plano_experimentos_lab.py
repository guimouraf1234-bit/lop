"""Plano de ensaios de laboratório: pontos novos dentro (refino) e fora (validação) da malha (CA0, eta) da tese.
Massa de sólido pela relação medida nos 16 ensaios: eta = 42,78 * CA0 / m  (V = 0,4 L)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from data_loader import carregar_respostas_finais

K = 42.78  # g . (mol/L)^-1 ; eta = K * CA0 / m
raiz = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# (id, tipo, CA0, eta_real, T_C, motivo)
plano = [
    ("D1", "dentro", 0.30, 0.80, 40, "pico de discordância entre 0,53 e 1,07 (início do patamar)"),
    ("D2", "dentro", 1.25, 0.80, 40, "idem, com CA0 alto"),
    ("D3", "dentro", 0.30, 1.30, 40, "patamar salta de 0,86 para 0,97 entre η 1,07 e 1,6"),
    ("D4", "dentro", 1.25, 1.30, 40, "idem, com CA0 alto"),
    ("D5", "dentro", 0.30, 2.40, 40, "maior vão da malha (1,6 a 3,2), onde o GP tem sd máxima"),
    ("D6", "dentro", 1.25, 2.40, 40, "idem, com CA0 alto"),
    ("D7", "dentro", 1.00, 1.07, 40, "réplica do ensaio 9 (piso de ruído; patamar 0,85 vs 1,0 do PBM)"),
    ("F1", "fora", 0.50, 0.30, 40, "η abaixo do domínio: testa o teto X = η/ν"),
    ("F2", "fora", 0.50, 5.00, 40, "η acima do domínio: saturação X = 1 e Fe"),
    ("F3", "fora", 2.00, 1.07, 40, "CA0 acima de 1,5: força iônica e associação de ZnSO4"),
    ("F4", "fora", 0.05, 1.07, 40, "CA0 abaixo de 0,1: [H+] muito baixa"),
    ("F5", "fora", 2.00, 0.75, 40, "canto: CA0 alto e η na região de maior incerteza"),
    ("F6", "fora", 0.50, 1.07, 60, "T acima de 40 °C: Arrhenius (ks de Balarini)"),
    ("F7", "fora", 0.50, 1.07, 25, "T abaixo de 40 °C: idem"),
]
df = pd.DataFrame(plano, columns=["id", "tipo", "CA0_mol_L", "eta_real", "T_C", "motivo"])
df["massa_solido_g"] = (K * df.CA0_mol_L / df.eta_real).round(2)
df["volume_L"] = 0.4
df.to_csv(os.path.join(raiz, "resultados", "plano_experimentos_lab.csv"), index=False, encoding="utf-8")
print(df.drop(columns="motivo").to_string(index=False))

# --- figura ---
r = carregar_respostas_finais()
inc = pd.read_csv(os.path.join(raiz, "resultados", "dominio_incerteza_eta.csv"))
inc = inc[inc.tempo_min == 15].groupby("eta_real").discordancia_max_min.mean()

fig, ax = plt.subplots(figsize=(10, 6.2), dpi=150)
xmax, ymax = 5.6, 2.2
# faixa de incerteza (discordância entre modelos) ao fundo, só depende de eta
e = np.linspace(inc.index.min(), inc.index.max(), 200)
u = np.interp(e, inc.index, inc.values)
im = ax.imshow(u[None, :], extent=[e[0], e[-1], 0, ymax], aspect="auto", cmap="Oranges", vmin=0, vmax=0.08, alpha=0.6, origin="lower", zorder=0)
cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.035)
cb.set_label("Discordância entre modelos (X, t = 15 min)")
# domínio da malha
ax.add_patch(plt.Rectangle((0.535, 0.1), 3.29 - 0.535, 1.4, fill=False, ls="--", lw=1.5, ec="#444", zorder=2))
ax.text(0.58, 1.53, "domínio da malha atual", fontsize=9, color="#444")
ax.scatter(r.razao_molar_eta_real, r.C_acid_0_mol_L, s=70, c="#1f4e79", edgecolor="w", zorder=4, label="Malha atual (16 ensaios, 40 °C)")
ax.axvline(1 / 1.0, color="grey", lw=0.8, ls=":")
ax.text(1.02, 2.12, "η = ν = 1 (estequiométrico)", fontsize=8, color="grey")
d = df[df.tipo == "dentro"]; f = df[(df.tipo == "fora") & (df.T_C == 40)]; t = df[df.T_C != 40]
ax.scatter(d.eta_real, d.CA0_mol_L, s=90, marker="D", c="#2a9d4a", edgecolor="k", zorder=5, label="Dentro: refinar o modelo (7)")
ax.scatter(f.eta_real, f.CA0_mol_L, s=110, marker="*", c="#d62728", edgecolor="k", zorder=5, label="Fora: validar (5, 40 °C)")
ax.scatter(t.eta_real, t.CA0_mol_L, s=110, marker="s", c="#9467bd", edgecolor="k", zorder=5, label="Fora: temperatura 25 e 60 °C (2)")
for _, p in df.iterrows():
    off = {"F6": (8, 8), "F7": (8, -14)}.get(p.id, (6, 6))
    ax.annotate(p.id + (f" ({p.T_C:.0f} °C)" if p.T_C != 40 else ""), (p.eta_real, p.CA0_mol_L), xytext=off, textcoords="offset points", fontsize=8)
ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)
ax.set_xlabel("η real (razão molar H₂SO₄ / Zn)"); ax.set_ylabel("C_A0 (mol/L)")
ax.set_title("Malha atual e ensaios propostos para o laboratório")
ax.legend(loc="upper right", fontsize=8, framealpha=0.95)
ax.grid(alpha=0.25)
fig.tight_layout()
saida = os.path.join(raiz, "resultados", "figuras", "f11_plano_experimentos_lab.png")
fig.savefig(saida)
print(saida)
