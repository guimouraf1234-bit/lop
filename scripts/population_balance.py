"""
Solucionador do Balanço Populacional (PBM) para Lixiviação em Batelada
Modelo da Dissertação de Bortot Coelho (2017): item 4.3, Tabela 5.9 e Apêndice A2.1,
com alpha = 0 (Eq. 3.42) como White-Box do modelo híbrido.

A configuração física (base da RRB, k_s(T), força motriz, especiação do H2SO4, nu, T, alpha) é um
`kinetics.ConfigWB` congelado. Sem `config`, `simular_ensaio` reproduz exatamente o comportamento legado:
  - formulacao='corrigida' (padrão): RRB em base mássica, ks referido a 30 °C (Arrhenius até a
    temperatura do ensaio) e força motriz [H+] pela especiação HSO4-/SO4-- (K2 = 0,0099), nu = 1.
  - formulacao='tese': como no Apêndice A2.1 (RRB em base numérica com peso D^3, ks = 1,8e4 um/min
    sobre a molaridade de H2SO4).

Dois integradores equivalentes (diferença < 1e-4 em X, verificada em tests/test_pbm.py):
  - `simular_ensaio`        : referência, integra a quadratura a cada passo da EDO (lento, ~0,25 s/ensaio);
  - `simular_ensaio_rapido` : usa a tabela X(ΔD), interpolada uma vez por instância (~10 ms/ensaio).
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp, quad

from kinetics import taxa_retricao_linear_v, ConfigWB, estado_quimico, ph_de_h


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

# Cache das tabelas X(ΔD) por (formulação, m, d632, d_min, d_max)
_TABELAS_CONVERSAO: Dict[Tuple, Tuple[np.ndarray, np.ndarray]] = {}
# Cache de instâncias por base da RRB (usado por `para_config`)
_INSTANCIAS: Dict[str, 'PopulationBalanceModel'] = {}


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

    @classmethod
    def para_config(cls, config: ConfigWB) -> 'PopulationBalanceModel':
        """Instância (compartilhada) cuja base da RRB corresponde a `config.rrb`."""
        formulacao = 'corrigida' if config.rrb == 'massica' else 'tese'
        if formulacao not in _INSTANCIAS:
            _INSTANCIAS[formulacao] = cls(formulacao=formulacao)
        return _INSTANCIAS[formulacao]

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

    # ------------------------------------------------------------------
    # Configuração e tabela de conversão
    # ------------------------------------------------------------------
    def _resolver_config(self, T_celsius: float, alpha: float, config: Optional[ConfigWB]) -> ConfigWB:
        base_inst = 'massica' if self.formulacao == 'corrigida' else 'numerica'
        if config is None:
            if self.formulacao == 'tese':
                return ConfigWB.tese(alpha=alpha)
            return ConfigWB(alpha=alpha, T_celsius=T_celsius)
        if config.rrb != base_inst:
            raise ValueError(f"config.rrb='{config.rrb}' é incompatível com esta instância "
                             f"(formulacao='{self.formulacao}'); use PopulationBalanceModel.para_config(config)")
        return config

    def tabela_conversao(self) -> Tuple[np.ndarray, np.ndarray]:
        """Tabela X(ΔD) usada pelo integrador rápido (calculada uma vez por configuração granulométrica)."""
        chave = (self.formulacao, self.m, self.d632, self.d_min, self.d_max)
        if chave not in _TABELAS_CONVERSAO:
            # malha fina nos primeiros 5 um, onde a conversão sobe mais rápido
            dd = np.concatenate([np.linspace(0.0, 5.0, 501),
                                 np.linspace(5.01, self.d_max - self.d_min, 3000)])
            xx = np.array([self.calcular_conversao_zincita(d) for d in dd])
            _TABELAS_CONVERSAO[chave] = (dd, xx)
        return _TABELAS_CONVERSAO[chave]

    # ------------------------------------------------------------------
    # Integradores
    # ------------------------------------------------------------------
    def simular_ensaio(self,
                       tempos_min: List[float],
                       C_acid_0_mol_L: float,
                       eta: float,
                       T_celsius: float = 40.0,
                       alpha: float = 0.0,
                       config: Optional[ConfigWB] = None) -> Dict[str, np.ndarray]:
        """
        Integração de REFERÊNCIA (quadratura recalculada a cada passo da EDO).

        Parâmetros:
            tempos_min: lista de instantes de tempo [min]
            C_acid_0_mol_L: concentração inicial de H2SO4 [mol/L]
            eta: razão molar H2SO4 / ZnO calculada pela Eq. 4.3 (ver calcular_eta)
            T_celsius, alpha: usados só quando `config` não é dado (compatibilidade com o código legado)
            config: ConfigWB; se dado, substitui T_celsius/alpha e define nu, especiação, etc.

        Retorna dicionário com:
            - tempo_min: tempos avaliados
            - X_pbm: conversão da zincita (0 a 1), comparável a X_zn_exp (Eq. 4.5)
            - C_acid_mol_L: concentração de ácido livre (H2SO4) ao longo do tempo
            - C_H_mol_L: concentração de H+ livre (NaN na formulação da tese)
            - pH_pbm: -log10([H+]) (NaN na formulação da tese)
            - delta_D_um: retração acumulada do diâmetro
        """
        cfg = self._resolver_config(T_celsius, alpha, config)
        ks = cfg.ks_um_min
        t_eval = np.array(tempos_min, dtype=float)
        t_max = float(np.max(t_eval))

        def ode_delta_d(t, y):
            delta_D = max(float(y[0]), 0.0)
            x_zno = self.calcular_conversao_zincita(delta_D)
            c_acid, c_motriz, _ = estado_quimico(C_acid_0_mol_L, x_zno, eta, cfg)
            return [taxa_retricao_linear_v(C_acid_mol_L=c_acid, ks_um_min=ks, alpha=cfg.alpha,
                                           C_acid_0_mol_L=C_acid_0_mol_L, C_motriz_mol_L=c_motriz)]

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
        return self._montar_saida(t_eval, delta_D_vec, x_zno_vec, C_acid_0_mol_L, eta, cfg)

    def simular_ensaio_rapido(self,
                              tempos_min: List[float],
                              C_acid_0_mol_L: float,
                              eta: float,
                              T_celsius: float = 40.0,
                              alpha: float = 0.0,
                              config: Optional[ConfigWB] = None) -> Dict[str, np.ndarray]:
        """Mesma interface e saída de `simular_ensaio`, usando a tabela X(ΔD) (≈ 25x mais rápido)."""
        cfg = self._resolver_config(T_celsius, alpha, config)
        ks = cfg.ks_um_min
        dd, xx = self.tabela_conversao()
        t_in = np.asarray(tempos_min, dtype=float)
        ordem = np.argsort(t_in, kind='stable')
        t_ord = t_in[ordem]
        t_max = float(t_ord[-1])

        def ode_delta_d(t, y):
            x = float(np.interp(y[0] if y[0] > 0.0 else 0.0, dd, xx, right=1.0))
            c_acid, c_motriz, _ = estado_quimico(C_acid_0_mol_L, x, eta, cfg)
            return [taxa_retricao_linear_v(C_acid_mol_L=c_acid, ks_um_min=ks, alpha=cfg.alpha,
                                           C_acid_0_mol_L=C_acid_0_mol_L, C_motriz_mol_L=c_motriz)]

        if t_max <= 0.0:
            delta_ord = np.zeros(len(t_ord))
        else:
            sol = solve_ivp(ode_delta_d, (0.0, t_max), [0.0], t_eval=t_ord,
                            method='LSODA', rtol=1e-7, atol=1e-9)
            if not sol.success:
                raise RuntimeError(f"Falha na integração do PBM: {sol.message}")
            delta_ord = sol.y[0]
        x_ord = np.interp(np.maximum(delta_ord, 0.0), dd, xx, right=1.0)

        # devolve na ordem original dos tempos
        delta_D_vec = np.empty_like(delta_ord)
        x_zno_vec = np.empty_like(x_ord)
        delta_D_vec[ordem] = delta_ord
        x_zno_vec[ordem] = x_ord
        return self._montar_saida(t_in, delta_D_vec, x_zno_vec, C_acid_0_mol_L, eta, cfg)

    @staticmethod
    def _montar_saida(t_eval, delta_D_vec, x_zno_vec, C_acid_0, eta, cfg: ConfigWB) -> Dict[str, np.ndarray]:
        c_acid_vec = np.empty(len(t_eval))
        c_h_vec = np.empty(len(t_eval))
        for i, x in enumerate(x_zno_vec):
            c_acid_vec[i], _, c_h_vec[i] = estado_quimico(C_acid_0, float(x), eta, cfg)
        ph_vec = np.array([ph_de_h(h) if not math.isnan(h) else float('nan') for h in c_h_vec])
        return {
            'tempo_min': np.asarray(t_eval, dtype=float),
            'delta_D_um': delta_D_vec,
            'X_pbm': x_zno_vec,
            'C_acid_mol_L': c_acid_vec,
            'C_H_mol_L': c_h_vec,
            'pH_pbm': ph_vec,
        }
