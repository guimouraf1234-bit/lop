import numpy as np

# Constante universal dos gases [J/(mol*K)]
R_GAS = 8.31446

# Parâmetros Cinéticos da Dissolução da Zincita (Balarini et al., 2025)
KS_REF_UM_MIN = 18000.0        # [um/min] a 40 °C
EA_KJ_MOL = 13.45              # Energia de ativação [kJ/mol]
EA_J_MOL = EA_KJ_MOL * 1000.0  # [J/mol]
T_REF_K = 313.15               # 40 °C em Kelvin

# Densidade molar da fase sólida (Bortot Coelho, 2017)
RHO_S_MOL_L = 69.2             # [mol/L de sólido]


def calcular_ks(temperatura_celsius: float) -> float:
    """
    Calcula a constante cinética superficial ks [um/min] ajustada por Arrhenius.
    """
    t_kelvin = temperatura_celsius + 273.15
    fator_arrhenius = np.exp(-(EA_J_MOL / R_GAS) * (1.0 / t_kelvin - 1.0 / T_REF_K))
    return float(KS_REF_UM_MIN * fator_arrhenius)


def concentracao_acido_livre(C_acid_0: float, X_zn: float, eta: float) -> float:
    """
    Calcula a concentração instantânea de H2SO4 livre no licor:
    CAf(t) = CA0 * (1 - X_zn / eta)
    """
    if eta <= 0:
        return 0.0
    c_livre = C_acid_0 * (1.0 - (X_zn / eta))
    # O ácido nunca pode ser negativo fisicamente
    return float(max(c_livre, 0.0))


def taxa_retricao_linear_v(C_acid_mol_L: float,
                           ks_um_min: float = KS_REF_UM_MIN,
                           alpha: float = 0.0,
                           C_acid_0_mol_L: float = 0.0) -> float:
    """
    Calcula a velocidade linear de diminuição do diâmetro da partícula:
    v(t) = dD/dt = (2 / rho_s) * [ ks * CAf - alpha * (CA0 - CAf) ]  [um/min]
    
    - alpha = 0.0 -> Modelo fundamental mecanicista puro
    - alpha = 5500.0 -> Modelo com amortecimento empírico da tese de 2017
    """
    c_f = max(C_acid_mol_L, 0.0)
    c_0 = max(C_acid_0_mol_L, c_f)
    
    # Força motriz da reação química superficial
    forca_motriz = max(ks_um_min * c_f - alpha * (c_0 - c_f), 0.0)
    
    taxa = (2.0 / RHO_S_MOL_L) * forca_motriz
    return float(taxa)
