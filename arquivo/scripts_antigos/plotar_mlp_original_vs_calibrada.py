"""
plotar_mlp_original_vs_calibrada.py
------------------------------------
Gera o gráfico comparativo da Grade 4x4 (16 ensaios de bancada) contendo EXCLUSIVAMENTE:
  1. White-Box Puro (PBM)
  2. Grey-Box MLP Original (12, 6, alpha=0.10)
  3. Grey-Box MLP Calibrada (12, 6, alpha=0.05)
  4. Pontos Experimentais Reais (X_zn_exp)

Avaliação em TESTE CEGO estrito (Leave-One-Group-Out: 15 ensaios treino / 1 ensaio teste).
"""

import os
import sys
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import clone

diretorio_scripts = os.path.dirname(os.path.abspath(__file__))
if diretorio_scripts not in sys.path:
    sys.path.append(diretorio_scripts)

from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel
from hybrid_model import (
    ModeloHibridoCinza,
    preparar_features_e_alvo,
    calcular_metricas,
    AcopladorFisico
)


def gerar_grafico_original_vs_calibrada():
    print("=" * 90)
    print("GERANDO COMPARATIVO EXCLUSIVO: MLP ORIGINAL vs MLP CALIBRADA (16 ENSAIOS)")
    print("=" * 90)

    # 1. Carrega dados
    df_exp = carregar_dados_bancada()
    hibrido_aux = ModeloHibridoCinza('Rede Neural (MLP Tanh)')
    df_base = hibrido_aux.preparar_dados_base(df_exp)
    X, y, grupos = preparar_features_e_alvo(df_base)

    pbm = PopulationBalanceModel()
    acoplador = AcopladorFisico()
    logo = LeaveOneGroupOut()
    t_cont = np.linspace(0.0, 15.0, 301)

    # 2. Pré-computa simulação PBM White-Box
    print("-> Pré-computando curvas PBM White-Box...")
    sim_base = {}
    for ens_id, grp in df_exp.groupby('ensaio_id', sort=True):
        eta = float(grp['razao_molar_eta'].iloc[0])
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        razao_sl = float(grp['razao_SL_g_L'].iloc[0])
        tempos_exp = grp['tempo_min'].values
        x_exp = grp['X_zn_exp'].values

        res_pbm = pbm.simular_ensaio(t_cont.tolist(), ca0, eta, alpha=0.0)
        sim_base[ens_id] = {
            'eta': eta,
            'ca0': ca0,
            'razao_sl': razao_sl,
            'tempos_exp': tempos_exp,
            'x_exp': x_exp,
            'x_pbm_cont': np.array(res_pbm['X_pbm']),
            'c_acid_pbm_cont': np.array(res_pbm['C_acid_mol_L'])
        }

    # 3. Definição dos 2 modelos a comparar
    modelos_comparados = {
        'original': {
            'rotulo': 'Grey-Box [MLP Original: (12, 6), $\\alpha=0.10$]',
            'pipeline': Pipeline([
                ('scaler', StandardScaler()),
                ('mlp', MLPRegressor(
                    hidden_layer_sizes=(12, 6),
                    activation='tanh',
                    solver='lbfgs',
                    alpha=0.10,
                    max_iter=3000,
                    random_state=42
                ))
            ]),
            'cor': '#1f77b4',       # Azul clássico
            'estilo': '-',
            'largura': 2.3
        },
        'calibrada': {
            'rotulo': 'Grey-Box [MLP Calibrada: (12, 6), $\\alpha=0.05$]',
            'pipeline': Pipeline([
                ('scaler', StandardScaler()),
                ('mlp', MLPRegressor(
                    hidden_layer_sizes=(12, 6),
                    activation='tanh',
                    solver='lbfgs',
                    alpha=0.05,
                    max_iter=3000,
                    random_state=42
                ))
            ]),
            'cor': '#9467bd',       # Roxo elegante
            'estilo': '-',
            'largura': 2.3
        }
    }

    # 4. Simulação em teste cego Leave-One-Group-Out
    resultados = {}
    for chave, mod in modelos_comparados.items():
        print(f"-> Simulando {mod['rotulo']} sob LOGO-CV...")
        sim_ensaios = {}
        y_real_all, y_pred_all = [], []

        for fold_idx, (train_idx, val_idx) in enumerate(logo.split(X, y, grupos)):
            ens_id = int(np.unique(grupos[val_idx])[0])
            b = sim_base[ens_id]

            m_fold = clone(mod['pipeline'])
            m_fold.fit(X[train_idx], y[train_idx])

            # Predição contínua
            df_cont = pd.DataFrame({
                'ensaio_id': ens_id,
                'tempo_min': t_cont,
                'razao_molar_eta': b['eta'],
                'C_acid_0_mol_L': b['ca0'],
                'razao_SL_g_L': b['razao_sl'],
                'X_pbm': b['x_pbm_cont'],
                'C_acid_pbm_mol_L': b['c_acid_pbm_cont']
            })
            Xc, _, _ = preparar_features_e_alvo(df_cont)
            delta_c = m_fold.predict(Xc)
            x_cont = acoplador.acoplar_e_restringir(df_cont, delta_c)

            # Predição discreta nos pontos experimentais
            df_pts = df_base[df_base['ensaio_id'] == ens_id].copy()
            Xp, _, _ = preparar_features_e_alvo(df_pts)
            delta_p = m_fold.predict(Xp)
            x_pts = acoplador.acoplar_e_restringir(df_pts, delta_p)

            m_local = calcular_metricas(b['x_exp'], x_pts)

            y_real_all.extend(b['x_exp'])
            y_pred_all.extend(x_pts)

            sim_ensaios[ens_id] = {
                'x_cont': x_cont,
                'x_pts': x_pts,
                'm_local': m_local
            }

        m_glob = calcular_metricas(np.array(y_real_all), np.array(y_pred_all))
        print(f"   Métricas LOGO: R² = {m_glob['R2']:.4f} | RMSE = {m_glob['RMSE']:.4f} | MAE = {m_glob['MAE']*100:.2f}%")

        resultados[chave] = {
            'meta': mod,
            'metricas': m_glob,
            'sim_ensaios': sim_ensaios
        }

    # 5. Renderização da Grade 4x4
    print("\n-> Renderizando Grade 4x4 com White-Box, Original, Calibrada e Pontos Experimentais...")
    fig, axes = plt.subplots(4, 4, figsize=(18, 14), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ens_id in range(1, 17):
        ax = axes_flat[ens_id - 1]
        b = sim_base[ens_id]

        # 1. White-Box Puro
        ax.plot(t_cont, b['x_pbm_cont'], 'r--', lw=1.5, alpha=0.80, label='White-Box Puro (PBM)')

        # 2. MLP Original
        res_orig = resultados['original']
        ax.plot(t_cont, res_orig['sim_ensaios'][ens_id]['x_cont'],
                color=res_orig['meta']['cor'], ls=res_orig['meta']['estilo'],
                lw=res_orig['meta']['largura'], alpha=0.90,
                label=res_orig['meta']['rotulo'])

        # 3. MLP Calibrada
        res_calib = resultados['calibrada']
        ax.plot(t_cont, res_calib['sim_ensaios'][ens_id]['x_cont'],
                color=res_calib['meta']['cor'], ls=res_calib['meta']['estilo'],
                lw=res_calib['meta']['largura'], alpha=0.90,
                label=res_calib['meta']['rotulo'])

        # 4. Pontos Experimentais Reais
        ax.plot(b['tempos_exp'], b['x_exp'], 'ko', ms=6, mfc='#2ca02c', mec='black', mew=1.0,
                zorder=10, label='Experimental ($X_{Zn}^{exp}$)')

        # Limite estequiométrico quando eta < 1
        if b['eta'] < 1.0:
            ax.axhline(b['eta'], color='purple', ls='-.', lw=1.0, alpha=0.5,
                       label=f'Teto $\\eta = {b["eta"]}$')

        ax.set_title(f"Ensaio {ens_id}: $\\eta = {b['eta']}$, $C_{{A0}} = {b['ca0']}$ mol/L",
                     fontsize=9.5, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)

        # Caixa com métricas comparativas locais de cada ensaio
        m_o = res_orig['sim_ensaios'][ens_id]['m_local']
        m_c = res_calib['sim_ensaios'][ens_id]['m_local']
        texto_box = (
            f"Orig.: $R^2={m_o['R2']:.3f}$ | $MAE={m_o['MAE']*100:.1f}\\%$\n"
            f"Calib.: $R^2={m_c['R2']:.3f}$ | $MAE={m_c['MAE']*100:.1f}\\%$"
        )
        ax.text(0.96, 0.06, texto_box, transform=ax.transAxes,
                fontsize=7.2, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='#ffffff', alpha=0.88, edgecolor='#cccccc'))

    for i in range(4):
        axes[3, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    # Legenda global única no topo
    handles, labels = axes[0, 0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(0.5, 0.995),
               ncol=5, frameon=True, fontsize=10.0, facecolor='#fafafa')

    m_o_g = resultados['original']['metricas']
    m_c_g = resultados['calibrada']['metricas']
    subtitulo_global = (
        f"Cinética de Lixiviação de Zinco em Reator Batelada: White-Box vs MLP Original vs MLP Calibrada\n"
        f"AVALIAÇÃO EM TESTE CEGO (LOGO-CV: 15 Treino / 1 Teste) — "
        f"MLP Original: R²={m_o_g['R2']:.4f}, MAE={m_o_g['MAE']*100:.2f}% | "
        f"MLP Calibrada: R²={m_c_g['R2']:.4f}, MAE={m_c_g['MAE']*100:.2f}%"
    )
    plt.suptitle(subtitulo_global, fontsize=12.0, fontweight='bold', y=1.025)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    nome_arquivo = "comparativo_mlp_original_vs_calibrada.png"
    pastas_destino = [
        os.path.abspath(os.path.join(diretorio_scripts, "..", "graficos")),
        os.path.abspath(os.path.join(diretorio_scripts, "..", "graficos2")),
        os.path.join(diretorio_scripts, "graficos"),
        os.path.join(diretorio_scripts, "docs")
    ]

    for p in pastas_destino:
        os.makedirs(p, exist_ok=True)

    caminho_principal = os.path.join(pastas_destino[0], nome_arquivo)
    plt.savefig(caminho_principal, bbox_inches='tight')
    plt.close()
    print(f"\n-> Gráfico salvo com sucesso em: {caminho_principal}")

    for p in pastas_destino[1:]:
        dest = os.path.join(p, nome_arquivo)
        shutil.copyfile(caminho_principal, dest)
        print(f"   -> Copiado para: {dest}")

    print("\nExecução concluída com total sucesso!")
    return {
        'arquivo': caminho_principal,
        'metricas_original': m_o_g,
        'metricas_calibrada': m_c_g
    }


if __name__ == "__main__":
    gerar_grafico_original_vs_calibrada()
