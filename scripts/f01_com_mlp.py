"""Variante da f01 com o híbrido MLP tanh (LOGO) sobreposto ao GP. Script avulso: não altera o pipeline principal.
Saída: resultados/figuras/f01b_cinetica_16_ensaios_logo_mlp.png (+ resultados/curvas_logo_mlp.csv)."""
import os
import sys
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figuras_batelada as F
import validacao as V
from data_loader import carregar_dados_bancada
from hybrid_model import (CONFIG_WB_PRINCIPAL, ENTRADAS_PRINCIPAL, construir_base, construir_base_condicao,
                          obter_especificacoes)
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

R, S = F.DIR_RESULTADOS, os.path.join(F.DIR_RESULTADOS, "figuras")
df = carregar_dados_bancada()
cfg = CONFIG_WB_PRINCIPAL
caminho_csv = os.path.join(R, "curvas_logo_mlp.csv")

if os.path.exists(caminho_csv):
    mlp = pd.read_csv(caminho_csv)
else:
    t_fino = np.round(np.linspace(0.0, 15.0, 151), 4)
    base = construir_base(df, cfg)
    linhas = []
    for dobra in V.esquema_logo(df):
        tr = base[base["ensaio_id"].isin(dobra.treino)].reset_index(drop=True)
        ens = dobra.teste[0]
        r = df[df["ensaio_id"] == ens].iloc[0]
        g = construir_base_condicao(float(r["razao_molar_eta_real"]), float(r["C_acid_0_mol_L"]), t_fino, cfg,
                                    ensaio_id=ens, eta_nominal=float(r["razao_molar_eta"]))
        est = V.ResidualML(obter_especificacoes()["mlp_tanh"], ENTRADAS_PRINCIPAL).ajustar(tr)
        x = V.AcopladorFisico(nu=cfg.nu).restringir(g, est.prever_bruto(g))
        linhas.append(pd.DataFrame({"ensaio_id": ens, "tempo_min": t_fino, "X_mlp": x}))
        print("LOGO ensaio", ens, flush=True)
    mlp = pd.concat(linhas, ignore_index=True)
    mlp.to_csv(caminho_csv, index=False)

F.configurar()
curvas = F.ler(R, "curvas_logo.csv")
pred = pd.read_csv(os.path.join(R, "predicoes_oof.csv.gz"))
pl = pred[pred["esquema"].str.startswith("LOGO") & (pred["tempo_min"] > 0)].copy()
pl["e2"] = (pl["X_exp"] - pl["x_hat"]) ** 2
rm = pl.groupby(["ensaio_id", "modelo"])["e2"].mean().pow(0.5)

COR_MLP = "#eb6834"  # laranja só para destacar a nova curva nesta prévia
etas, cas = [0.5, 1.0, 1.5, 3.1], [0.1, 0.5, 1.0, 1.5]
fig, axes = plt.subplots(4, 4, figsize=(11, 9.2), sharex=True, sharey=True)
for i, eta in enumerate(etas):
    for j, ca in enumerate(cas):
        ax = axes[i, j]
        ens = int(df[(df["razao_molar_eta"] == eta) & (df["C_acid_0_mol_L"] == ca)]["ensaio_id"].iloc[0])
        c = curvas[curvas["ensaio_id"] == ens]
        m = mlp[mlp["ensaio_id"] == ens]
        e = df[df["ensaio_id"] == ens]
        teto = min(1.0, float(e["razao_molar_eta_real"].iloc[0]) / cfg.nu)
        ax.plot(c["tempo_min"], c["X_pbm"], color=F.TINTA, lw=1.4, zorder=2)
        lo = np.clip(c["X_gp"] - 1.96 * c["sd_gp"], 0, teto)
        hi = np.clip(c["X_gp"] + 1.96 * c["sd_gp"], 0, teto)
        ax.fill_between(c["tempo_min"], lo, hi, color=F.AZUL, alpha=0.14, lw=0, zorder=1)
        ax.plot(c["tempo_min"], c["X_gp"], color=F.AZUL, lw=1.8, zorder=3)
        ax.plot(m["tempo_min"], m["X_mlp"], color=COR_MLP, lw=1.8, ls="--", zorder=4)
        ax.plot(e["tempo_min"], e["X_zn_exp"], ls="none", marker="o", ms=5, mfc=F.SUP, mec=F.TINTA, mew=1.1, zorder=5)
        r = [rm.get((ens, k), np.nan) for k in ("White-Box", "Híbrido gp [B]", "Híbrido mlp_tanh [B]")]
        ax.text(0.97, 0.06, f"RMSE  WB {F.fmt(r[0])}\nGP {F.fmt(r[1])}\nMLP {F.fmt(r[2])}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7.4, color=F.TINTA2, linespacing=1.25)
        ax.set_title(f"η = {eta:g}  ·  $C_{{A0}}$ = {ca:g} mol/L".replace(".", ","), fontsize=8.8, color=F.TINTA2,
                     fontweight="normal", pad=3)
        ax.set_xlim(-0.3, 15.5); ax.set_ylim(-0.02, 1.06); ax.set_xticks([0, 5, 10, 15])
        if j == 0:
            ax.set_ylabel("Conversão $X_{Zn}$")
        if i == 3:
            ax.set_xlabel("Tempo (min)")
hs = [Line2D([], [], ls="none", marker="o", ms=5, mfc=F.SUP, mec=F.TINTA, mew=1.1, label="Experimental"),
      Line2D([], [], color=F.TINTA, lw=1.4, label="White-Box de laboratório (CAf, ν = 1)"),
      Line2D([], [], color=F.AZUL, lw=1.8, label="Híbrido GP (LOGO); faixa de 95%"),
      Line2D([], [], color=COR_MLP, lw=1.8, ls="--", label="Híbrido MLP tanh (LOGO)")]
fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.975), ncol=4, columnspacing=1.5)
fig.suptitle("Cinética dos 16 ensaios: cada curva vem de um modelo que não viu aquele ensaio", y=1.0,
             fontsize=11, fontweight="semibold")
fig.tight_layout(rect=[0, 0, 1, 0.945])
F.salvar_fig(fig, S, "f01b_cinetica_16_ensaios_logo_mlp.png")
print("RMSE LOGO médio:", pl.groupby("modelo")["e2"].mean().pow(0.5).loc[["White-Box", "Híbrido gp [B]", "Híbrido mlp_tanh [B]"]].round(4).to_dict())
