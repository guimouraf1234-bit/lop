"""
Geração dos Gráficos de Cinética do Modelo SVR com Data Augmentation
Produz:
1. cinetica_16_ensaios_svr_augmentation.png (Grade 4x4 completa com todos os 16 ensaios)
2. comparativo_svr_com_vs_sem_augmentation.png (Comparativo direto dos 4 regimes estequiométricos)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

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


def construir_pipeline_svr(C: float = 5.0, epsilon: float = 0.005) -> Pipeline:
    return Pipeline([
        ('scaler', StandardScaler()),
        ('svr', SVR(kernel='rbf', C=C, epsilon=epsilon, gamma='scale'))
    ])


def gerar_grafico_16_ensaios_svr_aug():
    print("\n-> Gerando gráfico de 16 ensaios para SVR com Data Augmentation...")

    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    acoplador = AcopladorFisico()

    caminho_csv = os.path.join(diretorio_atual, "dados_aumentados_640pts.csv")
    if os.path.exists(caminho_csv):
        df_aug = pd.read_csv(caminho_csv)
    else:
        df_aug = gerar_base_aumentada(salvar_csv=True)

    # Treinamento do SVR com a base aumentada de 640 pontos
    pipe_svr_aug = construir_pipeline_svr(C=5.0, epsilon=0.005)
    pipe_svr_aug.fit(df_aug[FEATURES_FISICAS].values, df_aug['delta_X'].values)

    t_cont = np.linspace(0.0, 15.0, 400)
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

        # Simulação White-Box contínua
        res_pbm_cont = pbm.simular_ensaio(t_cont.tolist(), ca0, eta, alpha=0.0)
        x_pbm_cont = np.array(res_pbm_cont['X_pbm'])
        c_acid_pbm_cont = np.array(res_pbm_cont['C_acid_mol_L'])

        # Simulação White-Box empírica (2017)
        res_pbm_emp = pbm.simular_ensaio(t_cont.tolist(), ca0, eta, alpha=5500.0)
        x_pbm_emp = np.array(res_pbm_emp['X_pbm'])

        # Predição contínua Grey-Box SVR com Augmentation
        df_sim_cont = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': t_cont,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': x_pbm_cont,
            'C_acid_pbm_mol_L': c_acid_pbm_cont
        })
        delta_pred_cont = pipe_svr_aug.predict(df_sim_cont[FEATURES_FISICAS].values)
        x_grey_cont = acoplador.acoplar_e_restringir(df_sim_cont, delta_pred_cont)

        # Predição nos pontos discretos para métricas
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
        delta_pred_pts = pipe_svr_aug.predict(df_sim_pts[FEATURES_FISICAS].values)
        x_grey_pts = acoplador.acoplar_e_restringir(df_sim_pts, delta_pred_pts)
        m_ensaio = calcular_metricas(x_exp, x_grey_pts)
        metricas_ensaios.append(m_ensaio)

        # 1. White-Box Puro
        ax.plot(t_cont, x_pbm_cont, 'r--', lw=1.5, alpha=0.8, label='White-Box Puro (PBM)')

        # 2. PBM Empírico 2017
        ax.plot(t_cont, x_pbm_emp, color='gray', ls=':', lw=1.2, alpha=0.7, label='PBM Empírico (2017)')

        # 3. Curva Grey-Box SVR com Data Augmentation
        ax.plot(t_cont, x_grey_cont, color='#9467bd', lw=2.2, label='Grey-Box [SVR Otimizado + Data Augmentation]')

        # 4. Pontos Experimentais Reais
        ax.plot(tempos_exp, x_exp, 'ko', ms=6.5, mfc='#2ca02c', mec='black', mew=1.1,
                zorder=10, label='Experimental ($X_{Zn}^{exp}$)')

        # 5. Teto estequiométrico
        if eta < 1.0:
            ax.axhline(eta, color='purple', ls='-.', lw=1.0, alpha=0.65, label=f'Teto $\\eta = {eta}$')

        ax.set_title(f"Ensaio {ens_id}: $\\eta = {eta}$, $C_{{A0}} = {ca0}$ mol/L",
                     fontsize=10, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)

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
    plt.suptitle(f'Cinética de Lixiviação: Grey-Box [SVR Otimizado com Data Augmentation]\n'
                 f'$R^2_{{médio}} = {r2_medio:.4f}$ | $MAE_{{médio}} = {mae_medio:.2f}\\%$ (Avaliados nos Dados Experimentais Reais)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    caminho_fig = os.path.join(diretorio_atual, "cinetica_16_ensaios_svr_augmentation.png")
    plt.savefig(caminho_fig, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfico salvo com sucesso: {caminho_fig}")


def gerar_grafico_comparativo_svr_com_vs_sem_aug():
    print("\n-> Gerando comparativo SVR (Sem Augmentation vs Com Augmentation)...")

    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    acoplador = AcopladorFisico()

    # 1. SVR Original (128 pts)
    h_svr_orig = ModeloHibridoCinza('SVR (Otimizado)')
    df_base_real = h_svr_orig.preparar_dados_base(df_exp)
    h_svr_orig.treinar(df_base_real)

    # 2. SVR Augmentation (640 pts)
    caminho_csv = os.path.join(diretorio_atual, "dados_aumentados_640pts.csv")
    df_aug = pd.read_csv(caminho_csv)
    pipe_svr_aug = construir_pipeline_svr(C=5.0, epsilon=0.005)
    pipe_svr_aug.fit(df_aug[FEATURES_FISICAS].values, df_aug['delta_X'].values)

    ensaios_alvo = [1, 6, 8, 11]
    t_cont = np.linspace(0.0, 15.0, 400)

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for idx, ens_id in enumerate(ensaios_alvo):
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

        # Predição Sem Augmentation
        x_svr_sem = h_svr_orig.prever(df_sim_cont)

        # Predição Com Augmentation
        delta_svr_aug = pipe_svr_aug.predict(df_sim_cont[FEATURES_FISICAS].values)
        x_svr_com = acoplador.acoplar_e_restringir(df_sim_cont, delta_svr_aug)

        # 1. White-Box Puro
        ax.plot(t_cont, x_pbm_cont, 'r--', lw=1.5, alpha=0.7, label='White-Box Puro (PBM)')

        # 2. SVR Sem Augmentation (128 pts)
        ax.plot(t_cont, x_svr_sem, color='#1f77b4', ls='-', lw=2.2, alpha=0.9,
                label='SVR Otimizado (Sem Augmentation - 128 pts)')

        # 3. SVR Com Augmentation (640 pts)
        ax.plot(t_cont, x_svr_com, color='#9467bd', ls='-', lw=2.2, alpha=0.9,
                label='SVR Otimizado + Augmentation (640 pts)')

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

    plt.suptitle('Comparação de Suavidade Cinética: SVR Sem Augmentation vs SVR com Data Augmentation PCHIP',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout()

    caminho_comp = os.path.join(diretorio_atual, "comparativo_svr_com_vs_sem_augmentation.png")
    plt.savefig(caminho_comp, bbox_inches='tight')
    plt.close()
    print(f"-> Comparativo SVR salvo com sucesso: {caminho_comp}")


if __name__ == '__main__':
    gerar_grafico_16_ensaios_svr_aug()
    gerar_grafico_comparativo_svr_com_vs_sem_aug()
