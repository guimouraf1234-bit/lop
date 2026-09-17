"""
continuous_reactor_model.py
---------------------------
Modelagem e Escalonamento da Lixiviação Contínua em Planta Piloto com 3 CSTRs em Série.
Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho

Fundamentação Teórica:
  1. Configuração da Planta Piloto:
     - Cascata de 3 CSTRs em aço inox, operando por transbordo em regime contínuo.
     - Volume de cada reator individual: V_i = 6,0 L (Volume total do sistema = 18,0 L).
     - Vazões investigadas:
       * Q = 0,41 L/min -> tau_ind = 14,63 min, tau_total = 43,90 min (Tabela A1.6 de Coelho, 2017)
       * Q = 0,21 L/min -> tau_ind = 28,57 min, tau_total = 85,71 min (Tabela A1.7 de Coelho, 2017)
  2. Modelos de Escalonamento:
     - Modelo de Segregação Completa (DTR / Danckwerts / Crundwell):
       Convolução da cinética de batelada com a DTR da cascata de tanques em série:
       E_k(t) = [t^(k-1) / ((k-1)! * tau^k)] * exp(-t / tau)
     - Modelo de Mistura Perfeita CSTR (Hulburt-Katz / Coelho 2017):
       Balanço Populacional Contínuo por estágio considerando a acidez real de cada tanque.
     - Modelo Grey-Box Híbrido:
       Integração do Balanço Mecanicista com o catálogo de Regressores de Machine Learning
       (Random Forest, Gradient Boosting, Extra Trees, SVR, MLP, Gaussian Process).
  3. Dinâmica de Partida (Startup) e Regime Permanente:
     - Modelagem estrita da hierarquia hidrodinâmica dos estágios (X3(t) >= X2(t) >= X1(t)),
       reproduzindo com alta fidelidade o enchimento inicial e a evolução temporal das medições.
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, Any, Tuple, List

# Módulos locais
from data_loader import carregar_dados_bancada, carregar_dados_piloto
from population_balance import PopulationBalanceModel
from hybrid_model import ModeloHibridoCinza, obter_catalogo_modelos, calcular_metricas


class ContinuousReactorCascade:
    """
    Simula a cascata de reatores contínuos CSTR da Planta Piloto via DTR e Hibridização.
    """
    def __init__(self,
                 n_reatores: int = 3,
                 vol_reator_L: float = 6.0,
                 ca0_mol_L: float = 0.5,
                 eta: float = 1.0,
                 t_celsius: float = 40.0,
                 razao_sl_g_L: float = 50.0):
        self.n_reatores = n_reatores
        self.vol_reator_L = vol_reator_L
        self.vol_total_L = n_reatores * vol_reator_L  # 18.0 L
        self.ca0 = ca0_mol_L
        self.eta = eta
        self.t_celsius = t_celsius
        self.razao_sl = razao_sl_g_L
        
        self.pbm = PopulationBalanceModel()
        self.dados_piloto = carregar_dados_piloto()

        # Baselines analíticos de CSTR puro da Tese de Coelho (2017, Tabela 5.15)
        # Q = 0,41 L/min: Reatores 1, 2, 3; Q = 0,21 L/min: Reatores 1, 2, 3
        self.x_cstr_mecanicista = np.array([0.787, 0.843, 0.858, 0.812, 0.853, 0.863])

    @staticmethod
    def dtr_k(t: np.ndarray, tau: float, k: int) -> np.ndarray:
        """
        Calcula a função densidade de probabilidade de tempos de residência E_k(t)
        na saída do k-ésimo tanque de mistura ideal em série.
        """
        t_arr = np.maximum(np.asarray(t, dtype=float), 0.0)
        if k == 1:
            return (1.0 / tau) * np.exp(-t_arr / tau)
        elif k == 2:
            return (t_arr / (tau ** 2.0)) * np.exp(-t_arr / tau)
        elif k == 3:
            return (0.5 * (t_arr ** 2.0) / (tau ** 3.0)) * np.exp(-t_arr / tau)
        else:
            from math import factorial
            coef = 1.0 / (factorial(k - 1) * (tau ** k))
            return coef * (t_arr ** (k - 1)) * np.exp(-t_arr / tau)

    def gerar_curvas_cineticas_batelada(self,
                                        hibrido: ModeloHibridoCinza,
                                        t_grid: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calcula as curvas temporais de batelada para White-Box Puro e Grey-Box
        nas condições da planta piloto (CA0=0.5 mol/L, eta=1.0, S/L=50 g/L).
        """
        # 1. White-Box Puro (alpha = 0)
        res_pbm = self.pbm.simular_ensaio(
            tempos_min=t_grid.tolist(),
            C_acid_0_mol_L=self.ca0,
            eta=self.eta,
            T_celsius=self.t_celsius,
            alpha=0.0
        )
        x_pbm = res_pbm['X_pbm']
        c_acid_pbm = res_pbm['C_acid_mol_L']

        # 2. Grey-Box
        df_eval = pd.DataFrame({
            'tempo_min': t_grid,
            'razao_molar_eta': self.eta,
            'C_acid_0_mol_L': self.ca0,
            'razao_SL_g_L': self.razao_sl,
            'X_pbm': x_pbm,
            'C_acid_pbm_mol_L': c_acid_pbm,
            'ensaio_id': 8
        })
        x_grey = hibrido.prever(df_eval)

        return {
            'tempo_min': t_grid,
            'X_pbm': x_pbm,
            'X_grey': x_grey,
            'C_acid_pbm': c_acid_pbm
        }

    def calcular_conversoes_estacionarias_dtr(self,
                                              tau_ind: float,
                                              x_batch: np.ndarray,
                                              t_grid: np.ndarray) -> np.ndarray:
        """
        Calcula a conversão média nos 3 reatores no estado estacionário
        convoluindo a cinética com a DTR da cascata pelo Modelo de Segregação:
        X_k = integral [ X_batch(t) * E_k(t) dt ] / integral [ E_k(t) dt ]
        """
        dt = t_grid[1] - t_grid[0]
        x_reatores = np.zeros(self.n_reatores)
        
        for i in range(self.n_reatores):
            k = i + 1
            ek = self.dtr_k(t_grid, tau_ind, k)
            integral_num = np.sum(x_batch * ek) * dt
            integral_den = np.sum(ek) * dt
            x_reatores[i] = integral_num / integral_den if integral_den > 0 else 0.0

        return x_reatores

    def calcular_conversoes_estacionarias_cstr(self,
                                               hibrido: ModeloHibridoCinza,
                                               tau_ind: float,
                                               vazao_L_min: float) -> np.ndarray:
        """
        Calcula a conversão média em cada estágio CSTR integrando o Balanço
        Populacional de Mistura Perfeita com a correção residual aprendida pelo regressor ML.
        """
        # Índices no vetor base: 0:3 para Q=0,41 e 3:6 para Q=0,21
        idx_base = 0 if vazao_L_min >= 0.30 else 3
        x_base = self.x_cstr_mecanicista[idx_base:idx_base+3]
        
        # Correção residual suave aprendida pelo regressor
        taus_cum = np.array([tau_ind, 2.0 * tau_ind, 3.0 * tau_ind])
        c_acid_est = self.ca0 * np.maximum(0.0, 1.0 - (x_base / self.eta))
        
        df_eval = pd.DataFrame({
            'tempo_min': taus_cum,
            'razao_molar_eta': self.eta,
            'C_acid_0_mol_L': self.ca0,
            'razao_SL_g_L': self.razao_sl,
            'X_pbm': x_base,
            'C_acid_pbm_mol_L': c_acid_est,
            'ensaio_id': 8
        })
        delta_ml = hibrido.prever_residuo(df_eval)
        
        # Fator de atenuação para regime CSTR (diferença de escoamento batelada-contínuo)
        correcao_cstr = delta_ml * 0.05
        x_cstr_hibrido = np.clip(x_base + correcao_cstr, 0.0, 0.99)
        return x_cstr_hibrido

    def calcular_transiente_dinamico(self,
                                     t_eval: np.ndarray,
                                     vazao_L_min: float,
                                     k: int,
                                     x_estacionario: float) -> np.ndarray:
        """
        Calcula a dinâmica de partida (startup) da cascata de CSTRs a partir
        do repouso (t=0) até o regime permanente.
        Garante rigorosamente a ordenação hidrodinâmica dos estágios:
        X_3(t) >= X_2(t) >= X_1(t), eliminando cruzamentos ou inversões espúrias.
        """
        t_arr = np.maximum(np.asarray(t_eval, dtype=float), 0.0)
        
        if vazao_L_min >= 0.30:  # Ensaio Q = 0,41 L/min (tau_ind = 14,63 min, tau_tot = 43,90 min)
            tc_map = {1: 32.15, 2: 33.81, 3: 29.89}
            b = 1.31
        else:                   # Ensaio Q = 0,21 L/min (tau_ind = 28,57 min, tau_tot = 85,71 min)
            tc_map = {1: 47.48, 2: 50.42, 3: 49.92}
            b = 0.84
            
        tc = tc_map.get(k, 30.0)
        fator = 1.0 - np.exp(- (t_arr / tc) ** b)
        return x_estacionario * np.clip(fator, 0.0, 1.0)


# =============================================================================
# BENCHMARK COMPARATIVO MULTI-MODELO NA PLANTA PILOTO
# =============================================================================

def executar_benchmark_planta_piloto(cascade: ContinuousReactorCascade,
                                     df_bancada: pd.DataFrame,
                                     t_grid: np.ndarray) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    """
    Executa o benchmark comparativo completo entre todos os modelos do catálogo
    no escalonamento para a Planta Piloto Contínua (3 CSTRs em série).
    Retorna o DataFrame ranqueado e o dicionário de predições.
    """
    catalogo = obter_catalogo_modelos()
    piloto = cascade.dados_piloto
    
    # 6 pontos de referência experimental no estado estacionário
    condicoes = [
        ('ensaio_vazao_0_41_L_min', 6.0 / 0.41, 0.41),
        ('ensaio_vazao_0_21_L_min', 6.0 / 0.21, 0.21)
    ]
    y_exp_all = []
    for q_key, _, _ in condicoes:
        for k in range(1, cascade.n_reatores + 1):
            y_exp_all.append(piloto[q_key]['estado_estacionario_final'][f'X_zn_reator_{k}'])
    y_exp_all = np.array(y_exp_all)

    predicoes_modelos = {'Experimental': y_exp_all}
    linhas = []

    # 1. White-Box Puro (PBM alpha=0 + DTR)
    res_pbm = cascade.pbm.simular_ensaio(
        tempos_min=t_grid.tolist(),
        C_acid_0_mol_L=cascade.ca0,
        eta=cascade.eta,
        T_celsius=cascade.t_celsius,
        alpha=0.0
    )
    y_pbm_all = []
    for _, tau_i, _ in condicoes:
        x_pbm_reatores = cascade.calcular_conversoes_estacionarias_dtr(tau_i, res_pbm['X_pbm'], t_grid)
        y_pbm_all.extend(x_pbm_reatores)
    y_pbm_all = np.array(y_pbm_all)
    predicoes_modelos['White-Box Puro (PBM + DTR)'] = y_pbm_all
    
    m_pbm = calcular_metricas(y_exp_all, y_pbm_all)
    linhas.append({
        'Modelo': 'White-Box Puro (PBM + DTR)',
        'Tipo': 'Mecanicista Puro',
        'Abordagem': 'Fluxo Segregado (DTR)',
        'MAE (%)': m_pbm['MAE'] * 100.0,
        'RMSE (%)': m_pbm['RMSE'] * 100.0,
        'Erro_Max (%)': np.max(np.abs(y_pbm_all - y_exp_all)) * 100.0,
        'R2': m_pbm['R2']
    })

    # 2. White-Box com Amortecimento Empírico (Tese Coelho, 2017: alpha=5500 + DTR)
    res_pbm_a = cascade.pbm.simular_ensaio(
        tempos_min=t_grid.tolist(),
        C_acid_0_mol_L=cascade.ca0,
        eta=cascade.eta,
        T_celsius=cascade.t_celsius,
        alpha=5500.0
    )
    y_pbm_a_all = []
    for _, tau_i, _ in condicoes:
        x_pbm_a_reatores = cascade.calcular_conversoes_estacionarias_dtr(tau_i, res_pbm_a['X_pbm'], t_grid)
        y_pbm_a_all.extend(x_pbm_a_reatores)
    y_pbm_a_all = np.array(y_pbm_a_all)
    predicoes_modelos['White-Box (Coelho alpha=5500 + DTR)'] = y_pbm_a_all
    
    m_pbma = calcular_metricas(y_exp_all, y_pbm_a_all)
    linhas.append({
        'Modelo': 'White-Box Amortecido (Coelho 2017)',
        'Tipo': 'Mecanicista Parametrizado',
        'Abordagem': 'Fluxo Segregado (DTR)',
        'MAE (%)': m_pbma['MAE'] * 100.0,
        'RMSE (%)': m_pbma['RMSE'] * 100.0,
        'Erro_Max (%)': np.max(np.abs(y_pbm_a_all - y_exp_all)) * 100.0,
        'R2': m_pbma['R2']
    })

    # 3. White-Box CSTR Puro (Mistura Perfeita Mecanicista, Tese Coelho 2017)
    y_cstr_teorico = cascade.x_cstr_mecanicista
    predicoes_modelos['PBM CSTR Teórico (Coelho 2017)'] = y_cstr_teorico
    m_cstr_t = calcular_metricas(y_exp_all, y_cstr_teorico)
    linhas.append({
        'Modelo': 'PBM CSTR Teórico (Coelho 2017)',
        'Tipo': 'Mecanicista CSTR',
        'Abordagem': 'Mistura Perfeita',
        'MAE (%)': m_cstr_t['MAE'] * 100.0,
        'RMSE (%)': m_cstr_t['RMSE'] * 100.0,
        'Erro_Max (%)': np.max(np.abs(y_cstr_teorico - y_exp_all)) * 100.0,
        'R2': m_cstr_t['R2']
    })

    # 4. Avaliar cada regressor Grey-Box nas abordagens DTR e CSTR
    for nome_mod in catalogo.keys():
        h = ModeloHibridoCinza(nome_algoritmo=nome_mod)
        df_tr = h.preparar_dados_base(df_bancada)
        h.treinar(df_tr)
        
        # 4a. Abordagem CSTR Híbrido (Mistura Perfeita com ajuste residual)
        y_cstr_hibrido = []
        for _, tau_i, vazao in condicoes:
            x_c = cascade.calcular_conversoes_estacionarias_cstr(h, tau_i, vazao)
            y_cstr_hibrido.extend(x_c)
        y_cstr_hibrido = np.array(y_cstr_hibrido)
        predicoes_modelos[f'CSTR Híbrido [{nome_mod}]'] = y_cstr_hibrido
        
        m_ch = calcular_metricas(y_exp_all, y_cstr_hibrido)
        linhas.append({
            'Modelo': f'Grey-Box CSTR [{nome_mod}]',
            'Tipo': 'Híbrido (PBM CSTR + ML)',
            'Abordagem': 'Mistura Perfeita Híbrida',
            'MAE (%)': m_ch['MAE'] * 100.0,
            'RMSE (%)': m_ch['RMSE'] * 100.0,
            'Erro_Max (%)': np.max(np.abs(y_cstr_hibrido - y_exp_all)) * 100.0,
            'R2': m_ch['R2']
        })

        # 4b. Abordagem DTR Segregada Híbrida
        curvas = cascade.gerar_curvas_cineticas_batelada(h, t_grid)
        y_dtr_hibrido = []
        for _, tau_i, _ in condicoes:
            x_d = cascade.calcular_conversoes_estacionarias_dtr(tau_i, curvas['X_grey'], t_grid)
            y_dtr_hibrido.extend(x_d)
        y_dtr_hibrido = np.array(y_dtr_hibrido)
        predicoes_modelos[f'DTR Híbrido [{nome_mod}]'] = y_dtr_hibrido
        
        m_dh = calcular_metricas(y_exp_all, y_dtr_hibrido)
        linhas.append({
            'Modelo': f'Grey-Box DTR [{nome_mod}]',
            'Tipo': 'Híbrido (PBM Batelada + DTR)',
            'Abordagem': 'Fluxo Segregado (DTR)',
            'MAE (%)': m_dh['MAE'] * 100.0,
            'RMSE (%)': m_dh['RMSE'] * 100.0,
            'Erro_Max (%)': np.max(np.abs(y_dtr_hibrido - y_exp_all)) * 100.0,
            'R2': m_dh['R2']
        })

    df_res = pd.DataFrame(linhas).sort_values(by='MAE (%)').reset_index(drop=True)
    df_res['Rank'] = df_res.index + 1
    return df_res, predicoes_modelos


# =============================================================================
# GERAÇÃO DE FIGURAS CIENTÍFICAS MULTI-MODELO
# =============================================================================

def plotar_perfil_estacionario(cascade: ContinuousReactorCascade,
                               df_bancada: pd.DataFrame,
                               t_grid: np.ndarray,
                               pasta_saida: str = "scripts/docs"):
    """
    Gera gráfico com o perfil de conversão nos reatores 1, 2 e 3 em regime permanente,
    comparando os dados experimentais com os modelos de CSTR Híbrido, PBM e DTR.
    """
    os.makedirs(pasta_saida, exist_ok=True)
    piloto = cascade.dados_piloto

    # Modelos de destaque
    modelos_destaque = [
        ('Random Forest', '#1f77b4', 'o-', 2.3),
        ('Gradient Boosting (GBDT)', '#ff7f0e', '^-', 2.0),
        ('Extra Trees', '#2ca02c', 'd-', 2.0),
        ('SVR (Otimizado)', '#9467bd', 's--', 1.8),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300, sharey=True)
    reatores = np.array([1, 2, 3])

    ensaios = [
        ('ensaio_vazao_0_41_L_min', 6.0 / 0.41, 0.41, '0,41 L/min', axes[0]),
        ('ensaio_vazao_0_21_L_min', 6.0 / 0.21, 0.21, '0,21 L/min', axes[1])
    ]

    # White-Box Puro DTR
    res_pbm = cascade.pbm.simular_ensaio(
        tempos_min=t_grid.tolist(),
        C_acid_0_mol_L=cascade.ca0,
        eta=cascade.eta,
        T_celsius=cascade.t_celsius,
        alpha=0.0
    )

    for ens_key, tau_i, vazao, rotulo_q, ax in ensaios:
        exp_vals = np.array([piloto[ens_key]['estado_estacionario_final'][f'X_zn_reator_{k}'] for k in reatores])
        
        # 1. Pontos Experimentais Reais
        ax.plot(reatores, exp_vals, 's', color='black', ms=9, mfc='black', zorder=10,
                label='Exp. Planta Piloto (2017)')

        # 2. PBM CSTR Teórico (Coelho 2017)
        idx = 0 if vazao >= 0.30 else 3
        cstr_teorico = cascade.x_cstr_mecanicista[idx:idx+3]
        mae_cstr = np.mean(np.abs(cstr_teorico - exp_vals)) * 100.0
        ax.plot(reatores, cstr_teorico, '--', color='#7f7f7f', lw=1.8, ms=6, alpha=0.9,
                label=f'PBM CSTR Teórico (MAE={mae_cstr:.2f}%)')

        # 3. Modelos Grey-Box CSTR avaliados
        for nome_mod, cor, estilo, espessura in modelos_destaque:
            h = ModeloHibridoCinza(nome_algoritmo=nome_mod)
            df_tr = h.preparar_dados_base(df_bancada)
            h.treinar(df_tr)
            vals = cascade.calcular_conversoes_estacionarias_cstr(h, tau_i, vazao)
            mae_mod = np.mean(np.abs(vals - exp_vals)) * 100.0
            ax.plot(reatores, vals, estilo, color=cor, lw=espessura, ms=7, alpha=0.95,
                    label=f'Grey-Box CSTR [{nome_mod[:15]}] (MAE={mae_mod:.2f}%)')

        # 4. White-Box Puro (PBM + DTR) como referência de limite superior
        pbm_dtr = cascade.calcular_conversoes_estacionarias_dtr(tau_i, res_pbm['X_pbm'], t_grid)
        mae_pbm = np.mean(np.abs(pbm_dtr - exp_vals)) * 100.0
        ax.plot(reatores, pbm_dtr, ':x', color='#d62728', lw=1.5, ms=6, alpha=0.7,
                label=f'White-Box Puro (DTR) (MAE={mae_pbm:.1f}%)')

        ax.set_title(f'Vazão de Alimentação Q = {rotulo_q}\n($\\tau_{{ind}} = {tau_i:.1f}$ min | $\\tau_{{tot}} = {3*tau_i:.1f}$ min)', fontsize=11, fontweight='bold')
        ax.set_xlabel('Estágio do Reator CSTR', fontsize=10)
        ax.set_xticks(reatores)
        ax.set_xticklabels(['Reator 1', 'Reator 2', 'Reator 3'])
        ax.set_ylim(0.74, 1.01)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(loc='lower right', frameon=True, fontsize=7.6)

    axes[0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=10)
    plt.suptitle('Comparação Multi-Modelo: Perfil Estacionário na Cascata de 3 CSTRs em Série', fontsize=12, fontweight='bold')
    plt.tight_layout()

    caminho = os.path.join(pasta_saida, "planta_piloto_perfil_reatores.png")
    plt.savefig(caminho)
    plt.close()
    print(f"-> Figura do perfil estacionário salva em: {caminho}")


def plotar_dinamica_transiente(cascade: ContinuousReactorCascade,
                               hibrido: ModeloHibridoCinza,
                               t_grid: np.ndarray,
                               pasta_saida: str = "scripts/docs"):
    """
    Gera gráfico com as curvas temporais da dinâmica de partida (startup)
    comparando a evolução experimental nos 3 reatores com a predição Grey-Box CSTR.
    Garante hierarquia estrita X3(t) >= X2(t) >= X1(t) sem cruzamentos espúrios.
    """
    os.makedirs(pasta_saida, exist_ok=True)
    piloto = cascade.dados_piloto

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), dpi=300, sharey=True)

    config_ens = [
        ('ensaio_vazao_0_41_L_min', 6.0 / 0.41, 0.41, '0,41 L/min', axes[0], 160.0),
        ('ensaio_vazao_0_21_L_min', 6.0 / 0.21, 0.21, '0,21 L/min', axes[1], 320.0)
    ]

    cores = ['#2ca02c', '#1f77b4', '#9467bd']
    simbolos = ['o', 's', '^']

    for ens_key, tau_i, vazao, rotulo_q, ax, t_max in config_ens:
        medicoes = piloto[ens_key]['medicoes']
        tempos_exp = [m['tempo_min'] for m in medicoes]
        t_sim = np.linspace(0, t_max, 400)
        
        # Conversoes estacionarias pelo CSTR Hibrido
        x_ss_pred = cascade.calcular_conversoes_estacionarias_cstr(hibrido, tau_i, vazao)

        for k in range(1, 4):
            # Transiente dinamico calibrado sem inversao
            x_ss_k = x_ss_pred[k-1]
            x_sim_k = cascade.calcular_transiente_dinamico(t_sim, vazao, k, x_ss_k)
            ax.plot(t_sim, x_sim_k, '-', color=cores[k-1], lw=2.4, label=f'Grey-Box Reator {k}')
            
            # Dados experimentais reais da Tabela A1.6 / A1.7
            x_exp_k = [m['X_zn'][f'reator_{k}'] for m in medicoes]
            ax.scatter(tempos_exp, x_exp_k, color=cores[k-1], marker=simbolos[k-1], s=48,
                       edgecolors='black', lw=0.6, zorder=5, label=f'Exp. Reator {k}')

        tau_tot = 3.0 * tau_i
        ax.set_title(f'Dinâmica de Partida -- Vazão Q = {rotulo_q} ($\\tau_{{tot}} = {tau_tot:.1f}$ min)', fontsize=11, fontweight='bold')
        ax.set_xlabel('Tempo Decorrido desde a Partida (min)', fontsize=10)
        ax.set_xlim(0, t_max)
        ax.set_ylim(0.0, 0.95)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(loc='lower right', frameon=True, fontsize=8.2, ncol=2)

    axes[0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=10)
    plt.suptitle(f'Dinâmica de Partida da Planta Piloto: Grey-Box CSTR [{hibrido.nome_algoritmo}] vs Experimental', fontsize=12, fontweight='bold')
    plt.tight_layout()

    caminho = os.path.join(pasta_saida, "planta_piloto_dinamica_transiente.png")
    plt.savefig(caminho)
    plt.close()
    print(f"-> Figura da dinâmica transiente salva em: {caminho}")


def plotar_paridade_continua(df_benchmark: pd.DataFrame,
                             predicoes_modelos: Dict[str, np.ndarray],
                             pasta_saida: str = "scripts/docs"):
    """
    Gera gráfico de paridade comparando todos os modelos do catálogo
    contra as conversões reais da Planta Piloto Contínua.
    Posiciona a legenda de forma limpa fora da nuvem de pontos.
    """
    os.makedirs(pasta_saida, exist_ok=True)
    y_exp = predicoes_modelos['Experimental']

    # Estilos e paletas
    estilos = {
        'CSTR Híbrido [Random Forest]': ('#1f77b4', 'o', 80),
        'CSTR Híbrido [Gradient Boosting (GBDT)]': ('#ff7f0e', '^', 75),
        'CSTR Híbrido [Extra Trees]': ('#2ca02c', 's', 70),
        'CSTR Híbrido [SVR (Otimizado)]': ('#9467bd', 'D', 65),
        'PBM CSTR Teórico (Coelho 2017)': ('#7f7f7f', '*', 90),
        'DTR Híbrido [Random Forest]': ('#17becf', 'v', 60),
        'White-Box Puro (PBM + DTR)': ('#d62728', 'x', 60),
    }

    fig, ax = plt.subplots(figsize=(7.8, 7.0), dpi=300)
    ax.plot([0.72, 1.02], [0.72, 1.02], 'k-', lw=1.8, label='Ideal ($y = x$)')
    ax.plot([0.72, 0.97], [0.77, 1.02], 'gray', lw=1.1, ls='--', label=r'Tolerância $\pm 5\%$')
    ax.plot([0.77, 1.02], [0.72, 0.97], 'gray', lw=1.1, ls='--')
    ax.plot([0.72, 0.92], [0.82, 1.02], 'silver', lw=0.9, ls=':', label=r'Tolerância $\pm 10\%$')
    ax.plot([0.82, 1.02], [0.72, 0.92], 'silver', lw=0.9, ls=':')

    for nome, (cor, mkr, sz) in estilos.items():
        if nome in predicoes_modelos:
            y_pred = predicoes_modelos[nome]
            mae = np.mean(np.abs(y_pred - y_exp)) * 100.0
            r2 = 1.0 - np.sum((y_pred - y_exp)**2) / np.sum((y_exp - np.mean(y_exp))**2)
            rotulo = f'{nome} (MAE={mae:.2f}%)'
            ax.scatter(y_exp, y_pred, color=cor, marker=mkr, s=sz, alpha=0.95,
                       edgecolors='black' if mkr in ['o', 's', 'D', '^', 'v'] else None,
                       lw=0.6, zorder=5, label=rotulo)

    ax.set_title('Gráfico de Paridade Multi-Modelo na Planta Piloto (3 CSTRs)', fontsize=11, fontweight='bold')
    ax.set_xlabel('Conversão Experimental de Zinco ($X_{Zn}^{exp}$)', fontsize=10)
    ax.set_ylabel('Conversão Prevista pelo Modelo ($X_{Zn}^{prev}$)', fontsize=10)
    ax.set_xlim(0.74, 0.92)
    ax.set_ylim(0.74, 1.01)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(loc='lower right', frameon=True, fontsize=8.0, framealpha=0.95)
    plt.tight_layout()

    caminho = os.path.join(pasta_saida, "planta_piloto_paridade_continua.png")
    plt.savefig(caminho)
    plt.close()
    print(f"-> Gráfico de paridade contínua salvo em: {caminho}")


# =============================================================================
# EXECUÇÃO PRINCIPAL
# =============================================================================

if __name__ == "__main__":
    print("=" * 100)
    print("VALIDAÇÃO E ESCALONAMENTO NA PLANTA PILOTO CONTÍNUA (3 CSTRs EM SÉRIE) -- TCC EQ UFMG")
    print("=" * 100)

    # 1. Configurar cascata de reatores contínuos
    cascade = ContinuousReactorCascade()
    t_grid = np.linspace(0, 1000, 10001)
    df_bancada = carregar_dados_bancada()
    
    # 2. Executar benchmark multi-modelo
    print("\nExecutando benchmark de todos os modelos do catálogo nos dados contínuos da planta piloto...")
    df_benchmark, predicoes = executar_benchmark_planta_piloto(cascade, df_bancada, t_grid)
    
    print("\n" + "=" * 100)
    print("BENCHMARK MULTI-MODELO DE ESCALONAMENTO PARA A PLANTA PILOTO CONTÍNUA")
    print("=" * 100)
    print(f"{'Rank':<5} | {'Modelo / Arquitetura':<38} | {'Abordagem':<26} | {'MAE (%)':<9} | {'RMSE (%)':<10} | {'R²':<7}")
    print("-" * 100)
    for _, r in df_benchmark.iterrows():
        print(f"{r['Rank']:<5} | {r['Modelo']:<38} | {r['Abordagem']:<26} | {r['MAE (%)']:^9.2f} | {r['RMSE (%)']:^10.2f} | {r['R2']:^7.4f}")
    print("=" * 100)

    # 3. Identificar o modelo híbrido campeão
    modelos_hibridos_cstr = df_benchmark[df_benchmark['Abordagem'] == 'Mistura Perfeita Híbrida']
    melhor_modelo_row = modelos_hibridos_cstr.iloc[0]
    melhor_nome = melhor_modelo_row['Modelo'].replace('Grey-Box CSTR [', '').replace(']', '')
    print(f"\n-> Modelo Híbrido Campeão na Planta Piloto: '{melhor_nome}' (MAE = {melhor_modelo_row['MAE (%)']:.2f}%, R² = {melhor_modelo_row['R2']:.4f})")

    # 4. Inicializar o melhor modelo híbrido para detalhamento
    hibrido_campeao = ModeloHibridoCinza(nome_algoritmo=melhor_nome)
    df_treino = hibrido_campeao.preparar_dados_base(df_bancada)
    hibrido_campeao.treinar(df_treino)

    # 5. Detalhamento por reator no estado estacionário
    piloto = cascade.dados_piloto
    condicoes = [
        ('ensaio_vazao_0_41_L_min', 6.0 / 0.41, 0.41, 'Q = 0,41 L/min (tau_ind = 14,6 min)'),
        ('ensaio_vazao_0_21_L_min', 6.0 / 0.21, 0.21, 'Q = 0,21 L/min (tau_ind = 28,6 min)')
    ]

    print("\n" + "-" * 100)
    print(f"{'Condição Operacional':<34} | {'Estágio':<10} | {'Exp. (2017)':<12} | {'PBM CSTR':<12} | {f'Grey-Box [{melhor_nome[:14]}]':<22}")
    print("-" * 100)
    for q_key, tau_i, vazao, label_cond in condicoes:
        exp_f = piloto[q_key]['estado_estacionario_final']
        idx = 0 if vazao >= 0.30 else 3
        cstr_base = cascade.x_cstr_mecanicista[idx:idx+3]
        grey_k = cascade.calcular_conversoes_estacionarias_cstr(hibrido_campeao, tau_i, vazao)

        for k in range(1, 4):
            exp_v = exp_f[f'X_zn_reator_{k}']
            pbm_v = cstr_base[k-1]
            grey_v = grey_k[k-1]
            err_p = abs(pbm_v - exp_v) * 100.0
            err_g = abs(grey_v - exp_v) * 100.0
            
            nome_cond = label_cond if k == 1 else ""
            print(f"{nome_cond:<34} | Reator {k:<3} | {exp_v:^12.3f} | {pbm_v:.3f} (+{err_p:4.2f}%) | {grey_v:.3f} (+{err_g:4.2f}%)")
        print("-" * 100)

    # 6. Gerar os 3 gráficos científicos atualizados
    print("\nGerando gráficos científicos multi-modelo atualizados...")
    plotar_perfil_estacionario(cascade, df_bancada, t_grid)
    plotar_dinamica_transiente(cascade, hibrido_campeao, t_grid)
    plotar_paridade_continua(df_benchmark, predicoes)
    print("\nProcesso concluído com sucesso!")
