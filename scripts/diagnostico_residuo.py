import os
import pandas as pd
from data_loader import carregar_dados_bancada, obter_caminhos
from population_balance import PopulationBalanceModel


def calcular_residuos():
    df = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    
    x_pbm_lista = []
    c_acid_lista = []
    
    for _, grupo in df.groupby("ensaio_id"):
        ca0 = float(grupo["C_acid_0_mol_L"].iloc[0])
        eta = float(grupo["razao_molar_eta"].iloc[0])
        tempos = grupo["tempo_min"].tolist()
        
        res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        x_pbm_lista.extend(res["X_pbm"])
        c_acid_lista.extend(res["C_acid_mol_L"])
        
    df["X_pbm"] = x_pbm_lista
    df["C_acid_pbm_mol_L"] = c_acid_lista
    df["delta_X"] = df["X_zn_exp"] - df["X_pbm"]
    
    caminhos = obter_caminhos()
    caminho_saida = caminhos["dados_com_residuo"]
    os.makedirs(os.path.dirname(caminho_saida), exist_ok=True)
    df.to_csv(caminho_saida, index=False)
    print(f"Resíduos calculados e salvos com sucesso em: {caminho_saida}")
    return df


if __name__ == "__main__":
    calcular_residuos()
