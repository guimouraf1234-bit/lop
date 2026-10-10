"""
Módulo de Data Augmentation Temporal com Preservação de Princípios Físicos
Lixiviação de Zinco (Calcina de Três Marias)

Metodologia:
1. Interpolação PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) 
   sobre os pontos experimentais com garantia estrita de monotonicidade temporal (dX/dt >= 0).
2. Respeito rigoroso ao limite estequiométrico máximo (teto = min(1.0, eta)).
3. Resolução analítica/numérica do modelo mecanicista PBM (Population Balance Model)
   em todos os instantes temporais aumentados para cálculo de X_pbm e C_acid_pbm_mol_L.
4. Cálculo exato do resíduo delta_X = X_aug - X_pbm.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# Adiciona o diretório de scripts do projeto ao path
diretorio_atual = os.path.dirname(os.path.abspath(__file__))
diretorio_lop = os.path.dirname(diretorio_atual)
diretorio_scripts = os.path.join(diretorio_lop, "scripts")
if diretorio_scripts not in sys.path:
    sys.path.insert(0, diretorio_scripts)

from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel

FEATURES_FISICAS = [
    'tempo_min',
    'razao_molar_eta',
    'C_acid_0_mol_L',
    'razao_SL_g_L',
    'X_pbm',
    'C_acid_pbm_mol_L'
]


def criar_grade_temporal_densa() -> np.ndarray:
    """
    Gera uma grade temporal de 40 pontos por ensaio:
    - Alta densidade na fase rápida de dissolução (0 a 5 min)
    - Amostragem regular no platô assintótico (5 a 15 min)
    - Preserva estritamente os 8 instantes experimentais originais
    """
    t_rapido = np.linspace(0.0, 5.0, 25)
    t_lento = np.linspace(6.0, 15.0, 10)
    t_orig = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0])
    
    t_total = np.unique(np.concatenate([t_rapido, t_lento, t_orig]))
    t_total.sort()
    return t_total


def aumentar_ensaio(grp: pd.DataFrame, t_grid: np.ndarray, pbm: PopulationBalanceModel) -> pd.DataFrame:
    """
    Realiza a interpolação PCHIP e a simulação PBM para um único ensaio.
    """
    ens_id = int(grp['ensaio_id'].iloc[0])
    eta = float(grp['razao_molar_eta'].iloc[0])
    ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
    razao_sl = float(grp['razao_SL_g_L'].iloc[0])
    t_exp = grp['tempo_min'].values
    x_exp = grp['X_zn_exp'].values
    
    # 1. Teto físico estequiométrico
    teto = min(1.0, eta) if eta < 1.0 else 1.0
    
    # 2. Monotonicidade termodinâmica antes do PCHIP
    x_mono = np.clip(np.maximum.accumulate(x_exp), 0.0, teto)
    
    # 3. Interpolação PCHIP que preserva formato e monotonicidade
    pchip = PchipInterpolator(t_exp, x_mono)
    x_aug = np.clip(pchip(t_grid), 0.0, teto)
    
    # 4. Solução mecanicista PBM para a grade temporal densa
    res_pbm = pbm.simular_ensaio(
        tempos_min=t_grid.tolist(),
        C_acid_0_mol_L=ca0,
        eta=eta,
        alpha=0.0
    )
    x_pbm = np.array(res_pbm['X_pbm'])
    c_acid_pbm = np.array(res_pbm['C_acid_mol_L'])
    
    # 5. Cálculo do resíduo exato
    delta_x = x_aug - x_pbm
    
    return pd.DataFrame({
        'ensaio_id': ens_id,
        'tempo_min': t_grid,
        'razao_molar_eta': eta,
        'C_acid_0_mol_L': ca0,
        'razao_SL_g_L': razao_sl,
        'X_pbm': x_pbm,
        'C_acid_pbm_mol_L': c_acid_pbm,
        'X_zn_exp': x_aug,
        'delta_X': delta_x
    })


def gerar_base_aumentada(salvar_csv: bool = True) -> pd.DataFrame:
    """
    Gera o dataset completo aumentado (16 ensaios x 40 pontos = 640 linhas)
    e opcionalmente salva em CSV na pasta teste_augmentation.
    """
    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    t_grid = criar_grade_temporal_densa()
    
    dfs_aumentados = []
    for _, grp in df_exp.groupby('ensaio_id', sort=True):
        df_ens = aumentar_ensaio(grp, t_grid, pbm)
        dfs_aumentados.append(df_ens)
        
    df_final = pd.concat(dfs_aumentados, ignore_index=True)
    
    if salvar_csv:
        caminho_csv = os.path.join(diretorio_atual, "dados_aumentados_640pts.csv")
        df_final.to_csv(caminho_csv, index=False)
        print(f"-> Base aumentada salva com sucesso: {caminho_csv} ({len(df_final)} linhas)")
        
    return df_final


if __name__ == '__main__':
    df_aug = gerar_base_aumentada(salvar_csv=True)
    print("Resumo estatístico da base aumentada:")
    print(df_aug[['tempo_min', 'X_pbm', 'X_zn_exp', 'delta_X']].describe())
