"""
Geração dos Gráficos de Cinética do Modelo MLP com Data Augmentation
Produz:
1. cinetica_16_ensaios_mlp_tanh_augmentation.png (Grade 4x4 completa com todos os 16 ensaios)
2. comparativo_mlp_com_vs_sem_augmentation.png (Comparação direta dos 4 regimes estequiométricos)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

diretorio_atual = os.path.dirname(os.path.abspath(__file__))
diretorio_lop = os.path.dirname(diretorio_atual)
diretorio_scripts = os.path.join(diretorio_lop, "scripts")
if diretorio_scripts not in sys.path:
    sys.path.insert(0, diretorio_scripts)

from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel
from hybrid_model import (
    ModeloHibridoCinza,
    calcular_metricas,
    AcopladorFisico,
    FEATURES_FISICAS
)
from data_augmentation import gerar_base_aumentada
from treinar_e_avaliar_mlp import construir_pipeline_mlp


def gerar_grafico_16_ensaios_mlp_aug():
    print("\n-> Gerando gráfico de 16 ensaios para MLP com Data Augmentation...")
    
    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    acoplador = AcopladorFisico()

    # 1. Base aumentada para treino
    df_aug = gerar_base_aumentada(salvar_csv=True)

    # 2. Treinamento da Rede Neural com a base de 640 pontos
    pipe_mlp_aug = construir_pipeline_mlp(alpha=0.01)
    pipe_mlp_aug.fit(df_aug[FEATURES_FISICAS].values, df_aug['delta_X'].values)

    # 3. Grade temporal contínua para plotagem de curvas lisas (400 pontos)
    t_cont = np.linspace(0.0, 15.0, 400)

    # Configuração da figura 4x4
    fig, axes = plt.subplots(4, 4, figsize=(18, 14), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    metricas_ensaios = []

    for ens_id in range(1, 17):
        ax = axes_flat[ens_id - 1]
        grp = df_exp[df_exp['ensaio_id'] == ens_id]
        eta = float(grp['razao_molar_eta'].iloc[0])
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        razao_sl = float(grp['razao_SL_g_L'].iloc[0])
        tempos_exp = grp['tempo_min'].values
        x_exp = grp['X_zn_exp'].values

        # Simulação White-Box contínua (alpha = 0)
        res_pbm_cont = pbm.simular_ensaio(t_cont.tolist(), ca0, eta, alpha=0.0)
        x_pbm_cont = np.array(res_pbm_cont['X_pbm'])
        c_acid_pbm_cont = np.array(res_pbm_cont['C_acid_mol_L'])

        # Simulação White-Box empírica (alpha = 5500)
        res_pbm_emp = pbm.simular_ensaio(t_cont.tolist(), ca0, eta, alpha=5500.0)
        x_pbm_emp = np.array(res_pbm_emp['X_pbm'])

        # Predição Grey-Box na curva contínua
        df_sim_cont = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': t_cont,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': x_pbm_cont,
            'C_acid_pbm_mol_L': c_acid_pbm_cont
        })
        delta_pred_cont = pipe_mlp_aug.predict(df_sim_cont[FEATURES_FISICAS].values)
        x_grey_cont = acoplador.acoplar_e_restringir(df_sim_cont, delta_pred_cont)

        # Predição nos pontos experimentais reais para cálculo das métricas
        res_pbm_pts = pbm.simular_ensaio(tempos_exp.tolist(), ca0, eta, alpha=0.0)
        df_sim_pts = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': tempos_exp,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': np.array(res_pbm_pts['X_pbm']),
            'C_acid_pbm_mol_L': np.array(res_pbm_pts['C_acid_mol_L'])
        })
        delta_pred_pts = pipe_mlp_aug.predict(df_sim_pts[FEATURES_FISICAS].values)
        x_grey_pts = acoplador.acoplar_e_restringir(df_sim_pts, delta_pred_pts)
        m_ensaio = calcular_metricas(x_exp, x_grey_pts)
        metricas_ensaios.append(m_ensaio)

        # Pontos aumentados deste ensaio para exibir como apoio
        grp_aug = df_aug[df_aug['ensaio_id'] == ens_id]

        # 1. Curva White-Box Puro (alpha = 0)
        ax.plot(t_cont, x_pbm_cont, 'r--', lw=1.5, alpha=0.8, label='White-Box Puro (PBM)')

        # 2. Curva Empírica Coelho 2017 (alpha = 5500)
        ax.plot(t_cont, x_pbm_emp, color='gray', ls=':', lw=1.2, alpha=0.7, label='PBM Empírico (2017)')

        # 3. Curva Grey-Box MLP com Data Augmentation
        ax.plot(t_cont, x_grey_cont, color='#1f77b4', lw=2.2, label='Grey-Box [MLP Tanh + Data Augmentation]')

        # 4. Pontos Experimentais Reais
        ax.plot(tempos_exp, x_exp, 'ko', ms=6.5, mfc='#2ca02c', mec='black', mew=1.1,
                zorder=10, label='Experimental ($X_{Zn}^{exp}$)')

        # 5. Limite estequiométrico quando eta < 1
        if eta < 1.0:
            ax.axhline(eta, color='purple', ls='-.', lw=1.0, alpha=0.65,
                       label=f'Teto $\\eta = {eta}$')

        ax.set_title(f"Ensaio {ens_id}: $\\eta = {eta}$, $C_{{A0}} = {ca0}$ mol/L",
                     fontsize=10, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)

        # Caixa com métricas do ensaio
        texto_metrica = f"$R^2 = {m_ensaio['R2']:.3f}$\n$MAE = {m_ensaio['MAE']*100:.1f}\\%$"
        ax.text(0.96, 0.06, texto_metrica, transform=ax.transAxes,
                fontsize=7.5, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.85, edgecolor='#cccccc'))

    for i in range(4):
        axes[3, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.995),
               ncol=5, frameon=True, fontsize=10.5, facecolor='#fafafa')

    r2_medio = np.mean([m['R2'] for m in metricas_ensaios])
    mae_medio = np.mean([m['MAE'] for m in metricas_ensaios]) * 100
    plt.suptitle(f'Cinética de Lixiviação: Grey-Box [Rede Neural MLP Tanh + L-BFGS com Data Augmentation]\n'
                 f'$R^2_{{médio}} = {r2_medio:.4f}$ | $MAE_{{médio}} = {mae_medio:.2f}\\%$ (Avaliados nos Dados Experimentais Reais)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    caminho_fig = os.path.join(diretorio_atual, "cinetica_16_ensaios_mlp_tanh_augmentation.png")
    plt.savefig(caminho_fig, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfico salvo com sucesso: {caminho_fig}")


def gerar_grafico_comparativo_com_vs_sem_aug():
    print("\n-> Gerando gráfico comparativo (Com vs Sem Data Augmentation)...")
    
    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    acoplador = AcopladorFisico()

    # 1. Base real e treino original (128 pts)
    mlp_orig = ModeloHibridoCinza('Rede Neural (MLP Tanh + L-BFGS)')
    df_base_real = mlp_orig.preparar_dados_base(df_exp)
    mlp_orig.treinar(df_base_real)

    # 2. Base aumentada e treino aumentado (640 pts)
    df_aug = gerar_base_aumentada(salvar_csv=False)
    pipe_mlp_aug = construir_pipeline_mlp(alpha=0.01)
    pipe_mlp_aug.fit(df_aug[FEATURES_FISICAS].values, df_aug['delta_X'].values)

    # 4 ensaios representativos dos 4 patamares de eta:
    # Ensaio 1: eta=0.5 (reagente limitante)
    # Ensaio 6: eta=1.0, CA0=0.1
    # Ensaio 8: eta=1.0, CA0=0.5 (condição planta piloto)
    # Ensaio 11: eta=1.5 (excesso de ácido)
    ensaios_selecionados = [1, 6, 8, 11]

    t_cont = np.linspace(0.0, 15.0, 400)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for idx, ens_id in enumerate(ensaios_selecionados):
        ax = axes_flat[idx]
        grp = df_exp[df_exp['ensaio_id'] == ens_id]
        eta = float(grp['razao_molar_eta'].iloc[0])
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        razao_sl = float(grp['razao_SL_g_L'].iloc[0])
        tempos_exp = grp['tempo_min'].values
        x_exp = grp['X_zn_exp'].values

        # Simulação White-Box contínua
        res_pbm_cont = pbm.simular_ensaio(t_cont.tolist(), ca0, eta, alpha=0.0)
        x_pbm_cont = np.array(res_pbm_cont['X_pbm'])
        c_acid_pbm_cont = np.array(res_pbm_cont['C_acid_mol_L'])

        df_sim_cont = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': t_cont,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': x_pbm_cont,
            'C_acid_pbm_mol_L': c_acid_pbm_cont
        })

        # Curva sem augmentation
        x_sem_aug = mlp_orig.prever(df_sim_cont)

        # Curva com augmentation
        delta_aug = pipe_mlp_aug.predict(df_sim_cont[FEATURES_FISICAS].values)
        x_com_aug = acoplador.acoplar_e_restringir(df_sim_cont, delta_aug)

        # 1. White-Box Puro
        ax.plot(t_cont, x_pbm_cont, 'r--', lw=1.5, alpha=0.7, label='White-Box Puro (PBM)')

        # 2. MLP Sem Augmentation (128 pts)
        ax.plot(t_cont, x_sem_aug, color='#d62728', ls='-', lw=2.0, alpha=0.85,
                label='MLP Tanh (Sem Augmentation - 128 pts)')

        # 3. MLP Com Data Augmentation (640 pts)
        ax.plot(t_cont, x_com_aug, color='#1f77b4', ls='-', lw=2.4, alpha=0.95,
                label='MLP Tanh + Augmentation (640 pts)')

        # 4. Pontos Experimentais Reais
        ax.plot(tempos_exp, x_exp, 'ko', ms=7, mfc='#2ca02c', mec='black', mew=1.2,
                zorder=10, label='Dados Experimentais ($X_{Zn}^{exp}$)')

        if eta < 1.0:
            ax.axhline(eta, color='purple', ls='-.', lw=1.0, alpha=0.6, label=f'Teto $\\eta = {eta}$')

        ax.set_title(f"Ensaio {ens_id}: $\\eta = {eta}$, $C_{{A0}} = {ca0}$ mol/L",
                     fontsize=11, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.2, 15.3)
        ax.set_ylim(-0.02, 1.05)
        ax.legend(loc='lower right' if eta >= 1.0 else 'center right', fontsize=8.2, framealpha=0.9)

    for i in range(2):
        axes[1, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    plt.suptitle('Comparação de Desempenho e Aderência: MLP Sem Augmentation vs MLP com Data Augmentation PCHIP',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout()

    caminho_comp = os.path.join(diretorio_atual, "comparativo_com_vs_sem_augmentation.png")
    plt.savefig(caminho_comp, bbox_inches='tight')
    plt.close()
    print(f"-> Comparativo salvo com sucesso: {caminho_comp}")


if __name__ == '__main__':
    gerar_grafico_16_ensaios_mlp_aug()
    gerar_grafico_comparativo_com_vs_sem_aug()
