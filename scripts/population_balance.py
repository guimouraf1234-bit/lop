"""
Solucionador do Balanço Populacional (PBM) para Lixiviação em Batelada
Modelo da Dissertação de Bortot Coelho (2017): item 4.3, Tabela 5.9 e Apêndice A2.1,
com alpha = 0 (Eq. 3.42) como White-Box do modelo híbrido.

Formulações disponíveis:
  - 'corrigida' (padrão): RRB em base mássica, ks referido a 30 °C (Arrhenius até a
    temperatura do ensaio) e força motriz [H+] pela especiação HSO4-/SO4--.
  - 'tese': como no Apêndice A2.1 (RRB em base numérica com peso D^3, ks = 1,8e4 um/min
    sobre a molaridade de H2SO4).
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from typing import Dict, Any, List
from kinetics import (calcular_ks, taxa_retricao_linear_v, concentracao_acido_livre,
                      concentracao_h_livre, KS_REF_UM_MIN)


# Parâmetros granulométricos RRB do concentrado ustulado (Tabela 5.9)
# A RRB foi ajustada a frações mássicas passantes (Tabela 5.4 da tese).
RRB_M = 1.022          # Coeficiente de uniformidade
RRB_D632 = 41.65       # Diâmetro característico [um]
D_MIN_UM = 0.1         # Limite inferior da distribuição [um]
D_MAX_UM = 297.0       # Limite superior da distribuição [um] (peneira 50#)

# Cálculo da razão molar eta (Eq. 4.3)
TEOR_ZNO = 0.761       # Fração mássica de zincita no concentrado ustulado
MM_ZNO_G_MOL = 81.38   # Massa molar da zincita [g/mol]

FORMULACOES = ('corrigida', 'tese')


def calcular_eta(volume_L: float, C_acid_0_mol_L: float, massa_solido_g: float) -> float:
    """Razão molar eta = V*CA0 / (mB0*T_ZnO/MM_ZnO) (Eq. 4.3)."""
    return float(volume_L * C_acid_0_mol_L / (massa_solido_g * TEOR_ZNO / MM_ZNO_G_MOL))


class PopulationBalanceModel:
    """
    Implementa a solução do Balanço Populacional em batelada pelo Método das Características.

    formulacao='corrigida': f(D,0) é a densidade mássica q3 da RRB; cada classe retém a
        fração (1 - delta_D/D)^3 da sua massa.
    formulacao='tese': f(D,0) é tratada como densidade numérica (Eq. 3.45) e a conversão
        vem do momento de ordem 3 (Eq. 3.49), como no Apêndice A2.1.
    """
    def __init__(self, m: float = RRB_M, d632: float = RRB_D632,
                 d_min: float = D_MIN_UM, d_max: float = D_MAX_UM,
                 formulacao: str = 'corrigida'):
        if formulacao not in FORMULACOES:
            raise ValueError(f"formulacao deve ser uma de {FORMULACOES}")
        self.m = m
        self.d632 = d632
        self.d_min = d_min
        self.d_max = d_max
        self.formulacao = formulacao
        if formulacao == 'tese':
            # Momento de ordem 3 da população inicial (denominador da Eq. 3.49)
            self.I0, _ = quad(lambda D: D ** 3.0 * self.rrb_pdf(D), d_min, d_max, limit=200)
        else:
            # Massa total da distribuição truncada em [D_min, D_max]
            self.I0, _ = quad(self.rrb_pdf, d_min, d_max, limit=200)

    def rrb_pdf(self, D: float) -> float:
        """Densidade da distribuição RRB f(D,0) [1/um]."""
        if D <= 0:
            return 0.0
        term1 = (self.m / self.d632) * ((D / self.d632) ** (self.m - 1.0))
        term2 = np.exp(- ((D / self.d632) ** self.m))
        return float(term1 * term2)

    def calcular_conversao_zincita(self, delta_D: float) -> float:
        """
        Calcula a conversão da zincita a partir do encolhimento acumulado delta_D [um].
        corrigida: 1 - X = integral[max(delta_D, D_min), D_max] (1 - delta_D/D)^3 q3(D) dD / integral q3 dD
        tese (Eq. 3.49): X = 1 - integral[D_min + delta_D, D_max] (D - delta_D)^3 f dD / integral D^3 f dD
        """
        if delta_D <= 0:
            return 0.0
        if delta_D >= self.d_max - self.d_min:
            return 1.0

        if self.formulacao == 'tese':
            integrand = lambda D: ((D - delta_D) ** 3.0) * self.rrb_pdf(D)
            restante, _ = quad(integrand, self.d_min + delta_D, self.d_max, limit=200, epsrel=1e-6)
        else:
            integrand = lambda D: ((1.0 - delta_D / D) ** 3.0) * self.rrb_pdf(D)
            restante, _ = quad(integrand, max(delta_D, self.d_min), self.d_max, limit=200, epsrel=1e-6)

        fracao_restante = np.clip(restante / self.I0, 0.0, 1.0)
        return float(1.0 - fracao_restante)

    def simular_ensaio(self,
                       tempos_min: List[float],
                       C_acid_0_mol_L: float,
                       eta: float,
                       T_celsius: float = 40.0,
                       alpha: float = 0.0) -> Dict[str, np.ndarray]:
        """
        Simula a cinética de lixiviação para os tempos especificados.

        Parâmetros:
            tempos_min: lista de instantes de tempo [min]
            C_acid_0_mol_L: concentração inicial de H2SO4 [mol/L]
            eta: razão molar H2SO4 / ZnO calculada pela Eq. 4.3 (ver calcular_eta)
            T_celsius: temperatura da reação [°C] (padrão 40 °C)
            alpha: parâmetro ajustável da Eq. 4.1 [um/min] (0 = Eq. 3.42, padrão; 5500 = tese)

        Retorna dicionário com:
            - tempo_min: tempos avaliados
            - X_pbm: conversão da zincita (0 a 1), comparável a X_zn_exp (Eq. 4.5)
            - C_acid_mol_L: concentração de ácido livre (H2SO4) ao longo do tempo
            - C_H_mol_L: concentração de H+ livre (formulação corrigida; NaN na da tese)
            - delta_D_um: retração acumulada do diâmetro
        """
        corrigida = self.formulacao == 'corrigida'
        ks = calcular_ks(T_celsius) if corrigida else KS_REF_UM_MIN
        t_eval = np.array(tempos_min, dtype=float)
        t_max = float(np.max(t_eval))

        def ode_delta_d(t, y):
            delta_D = max(float(y[0]), 0.0)
            x_zno = self.calcular_conversao_zincita(delta_D)
            c_acid = concentracao_acido_livre(C_acid_0_mol_L, x_zno, eta)
            c_h = concentracao_h_livre(C_acid_0_mol_L, c_acid) if corrigida else None
            v = taxa_retricao_linear_v(
                C_acid_mol_L=c_acid,
                ks_um_min=ks,
                alpha=alpha,
                C_acid_0_mol_L=C_acid_0_mol_L,
                C_motriz_mol_L=c_h
            )
            return [v]

        if t_max <= 0.0:
            delta_D_vec = np.zeros(len(t_eval))
        else:
            sol = solve_ivp(
                fun=ode_delta_d,
                t_span=(0.0, t_max),
                y0=[0.0],
                t_eval=t_eval,
                method='RK45',
                rtol=1e-6,
                atol=1e-8
            )
            delta_D_vec = sol.y[0]
        x_zno_vec = np.array([self.calcular_conversao_zincita(d) for d in delta_D_vec])
        c_acid_vec = np.array([concentracao_acido_livre(C_acid_0_mol_L, x, eta) for x in x_zno_vec])
        if corrigida:
            c_h_vec = np.array([concentracao_h_livre(C_acid_0_mol_L, c) for c in c_acid_vec])
        else:
            c_h_vec = np.full(len(t_eval), np.nan)

        return {
            'tempo_min': t_eval,
            'delta_D_um': delta_D_vec,
            'X_pbm': x_zno_vec,
            'C_acid_mol_L': c_acid_vec,
            'C_H_mol_L': c_h_vec
        }
