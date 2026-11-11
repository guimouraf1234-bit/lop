"""Curvas previstas para os ensaios de validação F1-F3 (fora da malha), com os modelos treinados nos 16 ensaios.
Script avulso. Saída: resultados/figuras/f15_simulacao_F1_F3.png e resultados/simulacao_F1_F3.csv."""
import os
import sys
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figuras_batelada as F
from data_loader import carregar_dados_bancada
from hybrid_model import (CONFIG_WB_PRINCIPAL, ENTRADAS_PRINCIPAL, ModeloHibridoCinza, construir_base_condicao,
                          obter_especificacoes)
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

R = F.DIR_RESULTADOS
cfg = CONFIG_WB_PRINCIPAL
plano = pd.read_csv(os.path.join(R, "plano_6_ensaios_lab.csv"))
plano = plano[plano.grupo == "validar"].reset_index(drop=True)

df = carregar_dados_bancada()
modelos = {}
for s in ("gp", "extra_trees", "mlp_tanh"):
    m = ModeloHibridoCinza(obter_especificacoes()[s].nome, entradas=ENTRADAS_PRINCIPAL, config_wb=cfg)
    m.treinar(m.preparar_dados_base(df))
    modelos[s] = m

t = np.round(np.linspace(0.0, 15.0, 151), 4)
linhas = []
for _, p in plano.iterrows():
    g = construir_base_condicao(float(p.eta_real), float(p.C_A0_mol_L), t, cfg)
    xgp, sd = modelos["gp"].prever_com_incerteza(g)
    d = pd.DataFrame({"id": p.id, "tempo_min": t, "X_pbm": g["X_pbm"].to_numpy(), "X_gp": xgp, "sd_gp": sd})
    for s in ("extra_trees", "mlp_tanh"):
        d[f"X_{s}"] = modelos[s].acoplador.acoplar_e_restringir(g, modelos[s].prever_residuo(g))
    linhas.append(d)
sim = pd.concat(linhas, ignore_index=True)
sim.to_csv(os.path.join(R, "simulacao_F1_F3.csv"), index=False)

F.configurar()
COR_ET, COR_MLP = "#1c5cab", "#0d366b"
fig, axes = plt.subplots(1, 3, figsize=(11, 3.9), sharey=True)
for ax, (_, p) in zip(axes, plano.iterrows()):
    c = sim[sim.id == p.id]
    teto = min(1.0, p.eta_real / cfg.nu)
    lo = np.clip(c.X_gp - 1.96 * c.sd_gp, 0, teto); hi = np.clip(c.X_gp + 1.96 * c.sd_gp, 0, teto)
    ax.fill_between(c.tempo_min, lo, hi, color=F.AZUL, alpha=0.14, lw=0, zorder=1)
    ax.plot(c.tempo_min, c.X_pbm, color=F.TINTA, lw=1.4, zorder=2)
    ax.plot(c.tempo_min, c.X_gp, color=F.AZUL, lw=1.8, zorder=3)
    ax.plot(c.tempo_min, c.X_extra_trees, color=COR_ET, lw=1.5, ls="--", zorder=3)
    ax.plot(c.tempo_min, c.X_mlp_tanh, color=COR_MLP, lw=1.5, ls=":", zorder=3)
    for tt in (0.5, 1, 2, 3, 4, 5, 15):
        ax.axvline(tt, color=F.EIXO, lw=0.6, ls=(0, (2, 3)), zorder=0)
    f = c[c.tempo_min == 15.0].iloc[0]
    ax.text(0.97, 0.05, f"X(15 min)\nWB {F.fmt(f.X_pbm, 2)}  ·  GP {F.fmt(f.X_gp, 2)} ± {F.fmt(1.96 * f.sd_gp, 2)}\n"
            f"ET {F.fmt(f.X_extra_trees, 2)}  ·  MLP {F.fmt(f.X_mlp_tanh, 2)}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.6, color=F.TINTA2, linespacing=1.3)
    ax.set_title(f"{p.id}:  η = {p.eta_nominal:g}  ·  $C_{{A0}}$ = {p.C_A0_mol_L:g} mol/L  ·  {p.massa_solido_g:.1f} g"
                 .replace(".", ","), fontsize=9, color=F.TINTA2, fontweight="normal")
    ax.set_xlim(-0.3, 15.5); ax.set_ylim(-0.02, 1.06); ax.set_xticks([0, 5, 10, 15])
    ax.set_xlabel("Tempo (min)")
axes[0].set_ylabel("Conversão $X_{Zn}$")
hs = [Line2D([], [], color=F.TINTA, lw=1.4, label="White-Box (CAf, ν = 1)"),
      Line2D([], [], color=F.AZUL, lw=1.8, label="Híbrido GP; faixa de 95%"),
      Line2D([], [], color=COR_ET, lw=1.5, ls="--", label="Híbrido Extra Trees"),
      Line2D([], [], color=COR_MLP, lw=1.5, ls=":", label="Híbrido MLP tanh"),
      Line2D([], [], color=F.EIXO, lw=0.8, ls=(0, (2, 3)), label="Tempos de amostragem")]
fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.93), ncol=5, columnspacing=1.4)
fig.suptitle("Ensaios de validação (fora da malha): previsão dos modelos treinados com os 16 ensaios, 40 °C",
             y=1.0, fontsize=11, fontweight="semibold")
fig.tight_layout(rect=[0, 0, 1, 0.86])
F.salvar_fig(fig, os.path.join(R, "figuras"), "f15_simulacao_F1_F3.png")
print(sim[sim.tempo_min.isin([0.5, 1, 2, 5, 15])].round(3).to_string(index=False))
