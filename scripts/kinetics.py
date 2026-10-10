
import numpy as np

# Constante universal dos gases [J/(mol*K)]
R_GAS = 8.31446

# Parâmetros cinéticos (Balarini, 2009; Tabela 5.9 de Bortot Coelho, 2017)
# ks = 3,00e-2 cm/s = 1,8e4 um/min, estimado a 30 °C com C_A = [H+] (Balarini et al., 2025, Fig. 8).
# A tese o aplica diretamente a 40 °C sobre a molaridade de H2SO4.
KS_REF_UM_MIN = 18000.0        # [um/min]
ALPHA_TESE_UM_MIN = 5500.0     # Parâmetro ajustável alpha da tese [um/min] (item 5.2.4, Tabela 5.9)
# O White-Box deste trabalho usa alpha = 0 (Eq. 3.42): o resíduo que alpha corrigia fica a cargo do ML
EA_KJ_MOL = 13.45              # Energia de ativação [kJ/mol] (resumo de Balarini et al., 2025 e tese).
# A Tabela 4 do artigo traz 13,75 kJ/mol para a reação química; o efeito em ks(40 °C) é de 0,4% (1,186 contra 1,190).
EA_J_MOL = EA_KJ_MOL * 1000.0  # [J/mol]
T_REF_K = 303.15               # 30 °C em Kelvin (temperatura de estimação do ks)

# Segunda dissociação do ácido sulfúrico, HSO4- <-> H+ + SO4--. Valor que reproduz o
# [H+] = 0,0689 mol/L usado por Balarini para 6 g/L de H2SO4 (0,0612 mol/L).
K2_HSO4_MOL_L = 0.0099

# Densidade molar da zincita (Tabela 5.9; 6,92e4 mol/m3)
RHO_S_MOL_L = 69.2             # [mol/L de sólido]


def calcular_ks(temperatura_celsius: float) -> float:
    """
    Calcula a constante cinética superficial ks [um/min] em base [H+],
    corrigida por Arrhenius a partir da temperatura de estimação (30 °C).
    """
    t_kelvin = temperatura_celsius + 273.15
    fator_arrhenius = np.exp(-(EA_J_MOL / R_GAS) * (1.0 / t_kelvin - 1.0 / T_REF_K))
    return float(KS_REF_UM_MIN * fator_arrhenius)


def concentracao_acido_livre(C_acid_0: float, X_zn: float, eta: float) -> float:
    """
    Calcula a concentração instantânea de H2SO4 livre no licor (Eq. 3.51):
    CAf(t) = CA0 * (1 - X_zn / eta)
    X_zn é a conversão da zincita e eta = mols H2SO4 / mols ZnO (Eq. 4.3).
    """
    if eta <= 0:
        return 0.0
    c_livre = C_acid_0 * (1.0 - (X_zn / eta))
    # O ácido nunca pode ser negativo fisicamente
    return float(max(c_livre, 0.0))


def concentracao_h_livre(C_acid_0: float, C_acid_livre: float, K2: float = K2_HSO4_MOL_L) -> float:
    """
    Calcula [H+] no licor pela especiação HSO4- <-> H+ + SO4-- [mol/L].
    Prótons totais = 2*CAf (a dissolução da zincita consome 2 H+ por mol) e
    sulfato total = CA0 (o ZnSO4 formado mantém o sulfato em solução):
        h + [HSO4-] = 2*CAf,   [HSO4-] = h*(CA0 - [HSO4-]) / K2
    O sulfato do ZnSO4 tampona o H+ ao longo da lixiviação.
    """
    H_total = 2.0 * max(C_acid_livre, 0.0)
    S_total = max(C_acid_0, 0.0)
    a = S_total + H_total + K2
    hso4 = (a - np.sqrt(max(a * a - 4.0 * S_total * H_total, 0.0))) / 2.0
    return float(max(H_total - hso4, 0.0))


def taxa_retricao_linear_v(C_acid_mol_L: float,
                           ks_um_min: float,
                           alpha: float = 0.0,
                           C_acid_0_mol_L: float = 0.0,
                           C_motriz_mol_L: float = None) -> float:
    """
    Calcula a velocidade linear de diminuição do diâmetro da partícula (Eq. 4.1):
    v(t) = -dD/dt = (2 / rho_s) * [ ks * C - alpha * (CA0 - CAf) ]  [um/min]

    C é a concentração da força motriz: [H+] (C_motriz_mol_L) na formulação corrigida,
    ou CAf (H2SO4 livre) na formulação da tese, quando C_motriz_mol_L não é informado.

    - alpha = 0    -> Eq. 3.42 de Balarini (2009), White-Box deste trabalho (padrão)
    - alpha = 5500 -> Eq. 4.1 completa da tese (referência)
    """
    c_f = max(C_acid_mol_L, 0.0)
    c_0 = max(C_acid_0_mol_L, c_f)
    c_motriz = c_f if C_motriz_mol_L is None else max(C_motriz_mol_L, 0.0)

    # Força motriz da reação química superficial
    forca_motriz = max(ks_um_min * c_motriz - alpha * (c_0 - c_f), 0.0)

    taxa = (2.0 / RHO_S_MOL_L) * forca_motriz
    return float(taxa)
