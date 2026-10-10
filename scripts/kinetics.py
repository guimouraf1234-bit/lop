"""
kinetics.py
-----------
Cinética de dissolução da zincita, estequiometria do ácido e especiação do H2SO4.

Todo parâmetro do White-Box tem proveniência declarada (ver `ConfigWB`):
  - k_s, E_a:   Balarini (2009) e Balarini et al. (2025); k_s referido a 30 °C.
  - K2(T):      Apêndice A4 da tese (pKa2 = 1,96 a 25 °C, ΔH° = -5,2 kcal/mol, van't Hoff).
  - K_ZnSO4(T): IUPAC (Powell et al., Pure Appl. Chem. 85, 2249, 2013, Parte 5).
  - nu:         fator de fechamento estequiométrico do ácido (ver `nu_a_priori`).
Nenhum deles é ajustado aos 128 pontos de conversão.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Dict, Optional

import numpy as np
from scipy.optimize import brentq

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

# Densidade molar da zincita (Tabela 5.9; 6,92e4 mol/m3)
RHO_S_MOL_L = 69.2             # [mol/L de sólido]

# ---------------------------------------------------------------------------
# Especiação do H2SO4
# ---------------------------------------------------------------------------
# V0 (legado): K2 constante = 0,0099 mol/L, valor que reproduz o [H+] = 0,0689 mol/L usado por Balarini
#   para 6 g/L de H2SO4 (0,0612 mol/L) a 30 °C. Não tem dependência com a temperatura.
K2_HSO4_MOL_L = 0.0099

# V1: K2(T) por van't Hoff a partir dos dados do Apêndice A4 da tese (HSO4- <-> H+ + SO4--).
#   pKa2 = 1,96 a 25 °C (A4.1) e ΔH° = -5,2 kcal/mol (A4.1). Observação: a mesma tese lista Ka2 = 12,383e-3
#   (pKa2 = 1,91) na página seguinte; adota-se o pKa2 do texto (1,96), que reproduz o [H+] de Balarini a 30 °C.
PKA2_25C = 1.96
DH_DISSOC_HSO4_J_MOL = -5.2 * 4184.0
T_PADRAO_K = 298.15

# V2: V1 + par iônico ZnSO4(aq), Zn2+ + SO4(2-) <-> ZnSO4(aq).
#   IUPAC Technical Report, Powell et al., Pure Appl. Chem. 85 (2013) 2249 (DOI 10.1351/PAC-REP-13-06-03),
#   Parte 5 (sistemas Zn2+ + OH-, Cl-, CO3(2-), SO4(2-), PO4(3-)): log10 K1° = 2,304 +/- 0,042 (Im = 0, 25 °C)
#   e ΔrH° = +6,0 +/- 0,5 kJ/mol. ATENÇÃO: valores obtidos de resumo secundário (o texto primário não pôde ser
#   aberto); conferir no relatório antes de citar no TCC.
#   A constante é de força iônica zero e é aplicada SEM correção de atividade: superestima o pareamento em
#   licor concentrado (I > 1 mol/L). V2 deve ser lida como cota superior do efeito do par iônico.
LOGK_ZNSO4_25C = 2.304
DH_ZNSO4_J_MOL = 6000.0

ESPECIACOES = ('V0', 'V1', 'V2')

# Valor a priori de nu (ver `nu_a_priori`)
NU_A_PRIORI = 1.06


def _van_t_hoff(k_ref: float, delta_h_j_mol: float, t_celsius: float, t_ref_k: float = T_PADRAO_K) -> float:
    """K(T) = K(Tref) * exp[-ΔH/R * (1/T - 1/Tref)] (ΔH independente da temperatura)."""
    t_k = t_celsius + 273.15
    return float(k_ref * math.exp(-(delta_h_j_mol / R_GAS) * (1.0 / t_k - 1.0 / t_ref_k)))


def k2_hso4(t_celsius: float, especiacao: str = 'V1') -> float:
    """K2 da dissociação do HSO4- [mol/L] na especiação escolhida (V0: constante; V1/V2: van't Hoff)."""
    if especiacao == 'V0':
        return K2_HSO4_MOL_L
    return _van_t_hoff(10.0 ** (-PKA2_25C), DH_DISSOC_HSO4_J_MOL, t_celsius)


def k_znso4(t_celsius: float) -> float:
    """Constante de associação Zn2+ + SO4(2-) -> ZnSO4(aq) [L/mol], em I = 0, por van't Hoff."""
    return _van_t_hoff(10.0 ** LOGK_ZNSO4_25C, DH_ZNSO4_J_MOL, t_celsius)


def nu_a_priori() -> float:
    """
    Fator de fechamento estequiométrico do ácido, nu = mol de H2SO4 consumidos / mol de ZnO dissolvidos,
    estimado só pela composição da calcina (sem usar nenhum dado de conversão).

    Contas (por 100 g de calcina; composição elementar de ensaios_bancada/parametros_processo.json):
      ZnO: 76,1 g / 81,38 g/mol = 0,935 mol.
      Óxidos que também consomem 1 mol de ácido por mol: CaO (1,4% Ca = 0,035 mol), CuO (0,4% Cu = 0,0063 mol),
      CdO (0,4% Cd = 0,0036 mol), MnO (0,2% Mn = 0,0036 mol)  -> 0,0485 mol (se dissolverem por completo).
      Ferrita: ~6% do Fe lixiviado (X_Fe medido) -> 0,0025 mol de ZnFe2O4 x 4 H2SO4 = 0,0101 mol. Como a
      conversão (Eq. 4.5) desconta o Zn da ferrita, todo esse ácido é extra em relação ao ZnO dissolvido.
      nu = 1 + (0,0484 + 0,0101) / 0,935 = 1,063 (1,03 se o Ca já estiver como CaSO4 e não consumir ácido).
      Adota-se 1,06 com incerteza de +/- 0,02. Concorda com nu = eta/X_final = 1,07 dos ensaios com eta = 0,5
      (ácido totalmente consumido), a única estimativa direta que não depende do pH.
    """
    return NU_A_PRIORI


@dataclass(frozen=True)
class ConfigWB:
    """
    Configuração do modelo White-Box (PBM). Congelada e gravada junto com todo resultado.

    rrb          'massica'  : RRB integrada em base mássica (correta; este trabalho)
                 'numerica' : Eq. 3.49 da tese (RRB mássica tratada como densidade numérica, ponderada por D^3)
    ks_arrhenius True       : ks referido a 30 °C e corrigido por Arrhenius até T_celsius
                 False      : ks = 1,8e4 um/min a qualquer T (tese)
    forca_motriz 'H+'       : força motriz = [H+] pela especiação (base em que ks foi estimado)
                 'H2SO4'    : força motriz = molaridade do H2SO4 livre (tese)
    especiacao   V0 / V1 / V2 (só vale com forca_motriz = 'H+'; ver constantes acima)
    nu           mol de H2SO4 consumidos por mol de ZnO dissolvido: C_Af = C_A0 (1 - nu X / eta)
    """
    rrb: str = 'massica'
    ks_arrhenius: bool = True
    forca_motriz: str = 'H+'
    especiacao: str = 'V0'
    nu: float = 1.0
    T_celsius: float = 40.0
    alpha: float = 0.0

    def __post_init__(self):
        if self.rrb not in ('massica', 'numerica'):
            raise ValueError("rrb deve ser 'massica' ou 'numerica'")
        if self.forca_motriz not in ('H+', 'H2SO4'):
            raise ValueError("forca_motriz deve ser 'H+' ou 'H2SO4'")
        if self.especiacao not in ESPECIACOES:
            raise ValueError(f"especiacao deve ser uma de {ESPECIACOES}")
        if not self.nu > 0:
            raise ValueError("nu deve ser positivo")

    # ---- fábricas ----
    @classmethod
    def legado(cls) -> 'ConfigWB':
        """White-Box do TCC antes da revisão: RRB mássica, ks(T), [H+] com K2 = 0,0099, nu = 1."""
        return cls()

    @classmethod
    def tese(cls, alpha: float = 0.0) -> 'ConfigWB':
        """Formulação original da dissertação (Apêndice A2.1): RRB numérica, ks fixo, molaridade de H2SO4."""
        return cls(rrb='numerica', ks_arrhenius=False, forca_motriz='H2SO4', alpha=alpha)

    def com_nu(self, nu: float) -> 'ConfigWB':
        return replace(self, nu=float(nu))

    # ---- constantes derivadas ----
    @property
    def ks_um_min(self) -> float:
        return calcular_ks(self.T_celsius) if self.ks_arrhenius else KS_REF_UM_MIN

    @property
    def k2(self) -> float:
        return k2_hso4(self.T_celsius, self.especiacao)

    @property
    def kz(self) -> Optional[float]:
        return k_znso4(self.T_celsius) if self.especiacao == 'V2' else None

    def rotulo(self) -> str:
        if self.forca_motriz == 'H2SO4':
            base = f"tese[rrb={self.rrb},ks_T={'s' if self.ks_arrhenius else 'n'}]"
        else:
            base = f"{self.especiacao}[rrb={self.rrb}]"
        return f"{base} nu={self.nu:.3f} alpha={self.alpha:g} T={self.T_celsius:g}C"

    def para_dict(self) -> Dict:
        return {'rrb': self.rrb, 'ks_arrhenius': self.ks_arrhenius, 'forca_motriz': self.forca_motriz,
                'especiacao': self.especiacao, 'nu': self.nu, 'T_celsius': self.T_celsius, 'alpha': self.alpha}


def calcular_ks(temperatura_celsius: float) -> float:
    """
    Calcula a constante cinética superficial ks [um/min] em base [H+],
    corrigida por Arrhenius a partir da temperatura de estimação (30 °C).
    """
    t_kelvin = temperatura_celsius + 273.15
    fator_arrhenius = np.exp(-(EA_J_MOL / R_GAS) * (1.0 / t_kelvin - 1.0 / T_REF_K))
    return float(KS_REF_UM_MIN * fator_arrhenius)


def concentracao_acido_livre(C_acid_0: float, X_zn: float, eta: float, nu: float = 1.0) -> float:
    """
    Calcula a concentração instantânea de H2SO4 livre no licor (Eq. 3.51, generalizada por nu):
    CAf(t) = CA0 * (1 - nu * X_zn / eta)
    X_zn é a conversão da zincita, eta = mols H2SO4 / mols ZnO (Eq. 4.3) e nu = mols de ácido consumidos por
    mol de ZnO dissolvido (nu = 1 recupera a Eq. 3.51 da tese).
    """
    if eta <= 0:
        return 0.0
    c_livre = C_acid_0 * (1.0 - (nu * X_zn / eta))
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


def especiar_escalar(S_total: float, H_total: float, K2: float,
                     zn_total: float = 0.0, kz: Optional[float] = None) -> Dict[str, float]:
    """
    Especiação do sistema H+/HSO4-/SO4(2-) (e ZnSO4(aq), se kz for dado), em concentrações (sem atividades).

        h + b = H_total                      (prótons totais = 2 CAf)
        b + s + z = S_total                  (sulfato total = CA0)
        K2 = h s / b                         (2ª dissociação)
        z = kz s (zn_total - z)              (par iônico; Zn2+ livre = zn_total - z)

    Sem par iônico (kz None ou zn_total = 0) a solução é a raiz da quadrática (a mesma de
    `concentracao_h_livre`). Com par iônico, resolve-se a equação monótona f(s) = K2 b - (H - b) s por brentq.
    Retorna h, hso4, so4, znso4 [mol/L].
    """
    S = max(S_total, 0.0)
    H = max(H_total, 0.0)
    if S <= 0.0 or H <= 0.0:
        return {'h': 0.0, 'hso4': 0.0, 'so4': S, 'znso4': 0.0}

    if kz is None or zn_total <= 0.0:
        a = S + H + K2
        b = (a - math.sqrt(max(a * a - 4.0 * S * H, 0.0))) / 2.0
        return {'h': max(H - b, 0.0), 'hso4': b, 'so4': S - b, 'znso4': 0.0}

    zn = zn_total

    def b_de_s(s: float) -> float:
        z = kz * s * zn / (1.0 + kz * s)
        return S - s - z

    def f(s: float) -> float:
        b = b_de_s(s)
        return K2 * b - (H - b) * s

    # f(0+) = K2 S > 0 e f(S) < 0 (b(S) <= 0); a raiz física é única e está em (0, S)
    s = brentq(f, 1e-15, S, xtol=1e-15, rtol=1e-12, maxiter=200)
    b = b_de_s(s)
    z = kz * s * zn / (1.0 + kz * s)
    return {'h': max(H - b, 0.0), 'hso4': b, 'so4': s, 'znso4': z}


def especiar(C_acid_0: float, C_acid_livre: float, config: ConfigWB, zn_dissolvido: float = 0.0) -> Dict[str, float]:
    """Especiação segundo a configuração do White-Box. zn_dissolvido = X * CA0 / eta [mol/L] (só V2)."""
    return especiar_escalar(C_acid_0, 2.0 * C_acid_livre, config.k2, zn_dissolvido, config.kz)


def ph_de_h(h_mol_l: float) -> float:
    """pH = -log10([H+]) (sem coeficiente de atividade); NaN se [H+] = 0."""
    return float(-math.log10(h_mol_l)) if h_mol_l > 0 else float('nan')


def estado_quimico(C_acid_0: float, X: float, eta: float, config: ConfigWB):
    """
    Estado químico do licor dada a conversão X: (CAf, força motriz [mol/L], [H+]).
    Compartilhado entre a integração de referência e a rápida do PBM.
    """
    c_acid = concentracao_acido_livre(C_acid_0, X, eta, config.nu)
    if config.forca_motriz == 'H2SO4':
        return c_acid, c_acid, float('nan')
    zn = X * C_acid_0 / eta if config.especiacao == 'V2' else 0.0
    h = especiar_escalar(C_acid_0, 2.0 * c_acid, config.k2, zn, config.kz)['h']
    return c_acid, h, h


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
