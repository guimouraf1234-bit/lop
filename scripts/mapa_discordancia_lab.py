"""Mapa de discordância entre híbridos (GP, ET, GBDT, SVR) numa grade η nominal x C_A0, para escolher ensaios novos.
Treina os 4 híbridos com os 16 ensaios (White-Box principal) e salva resultados/mapa_discordancia_lab.csv."""
import os
import sys
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
from data_loader import carregar_dados_bancada
from hybrid_model import CONFIG_WB_PRINCIPAL, ENTRADAS_PRINCIPAL, ModeloHibridoCinza, construir_base_condicao, obter_especificacoes

raiz = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ETAS = [0.35, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.1, 4.0]      # nominal
CA0S = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
REAL_POR_NOMINAL = 1.0694   # η real / η nominal nos ensaios da tese (Eq. 4.3)
K_MASSA = 42.78             # η real = K * C_A0 / m, com V = 0,4 L
TEMPOS = np.array([0.5, 1.0, 2.0, 5.0, 15.0])

df = carregar_dados_bancada()
modelos = {}
for s in ("gp", "extra_trees", "gbdt", "svr"):
    m = ModeloHibridoCinza(obter_especificacoes()[s].nome, entradas=ENTRADAS_PRINCIPAL, config_wb=CONFIG_WB_PRINCIPAL)
    m.treinar(m.preparar_dados_base(df))
    modelos[s] = m

linhas = []
for en in ETAS:
    er = en * REAL_POR_NOMINAL
    for c in CA0S:
        g = construir_base_condicao(float(er), float(c), TEMPOS, CONFIG_WB_PRINCIPAL)
        mat = np.vstack([m.acoplador.acoplar_e_restringir(g, m.prever_residuo(g)) for m in modelos.values()])
        disc = mat.max(axis=0) - mat.min(axis=0)
        linhas.append({"eta_nominal": en, "eta_real": round(er, 3), "C_A0": c, "massa_g": round(K_MASSA * c / er, 2),
                       "discordancia_media": float(disc.mean()), "discordancia_15min": float(disc[-1]),
                       "X_pbm_15min": float(g["X_pbm"].iloc[-1]), "X_hibrido_medio_15min": float(mat[:, -1].mean())})
d = pd.DataFrame(linhas)
d.to_csv(os.path.join(raiz, "resultados", "mapa_discordancia_lab.csv"), index=False)
print(d.pivot(index="C_A0", columns="eta_nominal", values="discordancia_media").iloc[::-1].round(3).to_string())
print(d.pivot(index="C_A0", columns="eta_nominal", values="massa_g").iloc[::-1].round(1).to_string())
