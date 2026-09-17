"""
Passo 4: Diagnóstico do Resíduo da Física Pura (White-Box)
Calcula X_pbm para todos os 128 pontos e quantifica Delta_X = X_exp - X_pbm
"""

import pandas as pd
from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel


def calcular_residuos():
    df = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    
    x_pbm_lista = []
    
    for _, grupo in df.groupby("ensaio_id"):
        ca0 = grupo["C_acid_0_mol_L"].iloc[0]
        eta = grupo["razao_molar_eta"].iloc[0]
        tempos = grupo["tempo_min"].tolist()
        
        res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        x_pbm_lista.extend(res["X_pbm"])
        
    df["X_pbm"] = x_pbm_lista
    df["delta_X"] = df["X_zn_exp"] - df["X_pbm"]
    
    # Salva o arquivo com os resíduos prontos
    df.to_csv("scripts/dados_com_residuo.csv", index=False)
    return df


if __name__ == "__main__":
    df = calcular_residuos()
    print("Resíduos calculados e salvos com sucesso!")
