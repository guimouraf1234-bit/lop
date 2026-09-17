"""
Solucionador do Balanço Populacional (PBM) para Lixiviação em Batelada
Equações 3.12 a 3.16 da Entrega Parcial II e Dissertação de Bortot Coelho (2017)
"""

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import gamma
from typing import Dict, Any, List
from kinetics import calcular_ks, taxa_retricao_linear_v, concentracao_acido_livre


# Parâmetros Granulométricos RRB da Calcina de Três Marias (Bortot Coelho, 2017)
RRB_M = 1.022          # Coeficiente de uniformidade
RRB_D632 = 41.65       # Diâmetro característico [um]
FRACAO_ZINCITA = 0.87  # Fração de zinco contido na zincita solúvel (Tabela 3.3)


class PopulationBalanceModel:
    """
    Implementa a solução do Balanço Populacional em batelada pelo Método das Características.
    """
    def __init__(self, m: float = RRB_M, d632: float = RRB_D632, f_zno: float = FRACAO_ZINCITA):
        self.m = m
        self.d632 = d632
        self.f_zno = f_zno
        
        # Terceiro momento volumétrico inicial analítico: I0 = integral(D^3 * f(D,0) dD)
        # Para a distribuição de Weibull/RRB, a integral exata é (d632^3) * Gamma(1 + 3/m)
        self.I0 = (self.d632 ** 3.0) * float(gamma(1.0 + 3.0 / self.m))

    def rrb_pdf(self, D: float) -> float:
        """Função Densidade de Probabilidade da distribuição RRB."""
        if D <= 0:
            return 0.0
        term1 = (self.m / self.d632) * ((D / self.d632) ** (self.m - 1.0))
        term2 = np.exp(- ((D / self.d632) ** self.m))
        return float(term1 * term2)

    def calcular_conversao_zincita(self, delta_D: float) -> float:
        """
        Calcula a conversão da zincita a partir do encolhimento acumulado delta_D [um].
        1 - X = integral[delta_D a inf] (D - delta_D)^3 * f(D,0) dD / I0
        """
        if delta_D <= 0:
            return 0.0
        
        D_upper = max(self.d632 * 20.0, delta_D + 600.0)
        integrand = lambda D: ((D - delta_D) ** 3.0) * self.rrb_pdf(D)
        I_restante, _ = quad(integrand, delta_D, D_upper, limit=100, epsrel=1e-5)
        
        fracao_restante = np.clip(I_restante / self.I0, 0.0, 1.0)
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
            eta: razão molar estequiométrica H2SO4 / ZnO
            T_celsius: temperatura da reação [°C] (padrão 40 °C)
            alpha: fator de amortecimento empírico [um/min] (0.0 = mecanicista puro; 5500.0 = dissertação 2017)
            
        Retorna dicionário com:
            - tempo_min: tempos avaliados
            - X_pbm: conversão fracional da zincita (0 a 1), correspondente a X_Zn da Equação (4.5) da tese
            - X_zno: sinônimo de X_pbm
            - X_total_rocha: conversão ponderada pelo teor de zincita na rocha (f_ZnO * X_zno)
            - C_acid_mol_L: concentração de ácido livre ao longo do tempo
            - delta_D_um: retração acumulada do diâmetro
        """
        ks = calcular_ks(T_celsius)
        t_eval = np.array(tempos_min, dtype=float)
        t_max = float(np.max(t_eval))
        
        def ode_delta_d(t, y):
            delta_D = max(float(y[0]), 0.0)
            x_zno = self.calcular_conversao_zincita(delta_D)
            c_acid = concentracao_acido_livre(C_acid_0_mol_L, x_zno, eta)
            v = taxa_retricao_linear_v(
                C_acid_mol_L=c_acid,
                ks_um_min=ks,
                alpha=alpha,
                C_acid_0_mol_L=C_acid_0_mol_L
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
                rtol=1e-5,
                atol=1e-7
            )
            delta_D_vec = sol.y[0]
        x_zno_vec = np.array([self.calcular_conversao_zincita(d) for d in delta_D_vec])
        c_acid_vec = np.array([concentracao_acido_livre(C_acid_0_mol_L, x, eta) for x in x_zno_vec])
        
        # Recuperação mineral ponderada pelo teor de zincita na calcina bruta (87%)
        x_total_rocha = self.f_zno * x_zno_vec

        return {
            'tempo_min': t_eval,
            'delta_D_um': delta_D_vec,
            'X_pbm': x_zno_vec,
            'X_zno': x_zno_vec,
            'X_total_rocha': x_total_rocha,
            'C_acid_mol_L': c_acid_vec
        }

