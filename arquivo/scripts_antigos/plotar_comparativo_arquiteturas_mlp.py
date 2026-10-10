"""
plotar_comparativo_arquiteturas_mlp.py
---------------------------------------
Gera um painel comparativo 4x4 (16 ensaios de bancada) confrontando 4 abordagens
arquiteturais distintas da Rede Neural Multilayer Perceptron com ativação Tangente Hiperbólica
e otimizador L-BFGS, todas avaliadas sob validação em TESTE CEGO (Leave-One-Group-Out):

  1. Abordagem A (Azul): Original TCC -> (12, 6), alpha = 0.10 (Referência clássica suave)
  2. Abordagem B (Verde): Rasa 1-Camada -> (12,), alpha = 0.03 (Aproximação universal sem quinas)
  3. Abordagem C (Laranja): Compacta 2-Camadas -> (10, 5), alpha = 0.05 (Menos parâmetros livres)
  4. Abordagem D (Roxo): Calibrada 2-Camadas -> (12, 6), alpha = 0.05 (Mesma topologia, alpha ajustado)
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

ARQUITETURAS = [
    {
        'id': 'orig',
        'nome': 'Original TCC: (12, 6), alpha=0.10',
        'rotulo_legenda': 'MLP Original: (12, 6), $\\alpha=0.10$',
        'layers': (12, 6),
        'alpha': 0.10,
        'cor': '#1f77b4',       # Azul clássico
        'estilo': '-',
        'largura': 2.2
    },
    {
        'id': 'rasa_12',
        'nome': 'Rasa 1-Camada: (12,), alpha=0.03',
        'rotulo_legenda': 'MLP Rasa 1-Camada: (12,), $\\alpha=0.03$',
        'layers': (12,),
        'alpha': 0.03,
        'cor': '#2ca02c',       # Verde esmeralda
        'estilo': '-',
        'largura': 2.0
    },
    {
        'id': 'comp_10_5',
        'nome': 'Compacta: (10, 5), alpha=0.05',
        'rotulo_legenda': 'MLP Compacta: (10, 5), $\\alpha=0.05$',
        'layers': (10, 5),
        'alpha': 0.05,
        'cor': '#ff7f0e',       # Laranja
        'estilo': '--',
        'largura': 2.0
    },
    {
        'id': 'ajust_12_6',
        'nome': 'Calibrada: (12, 6), alpha=0.05',
        'rotulo_legenda': 'MLP Calibrada: (12, 6), $\\alpha=0.05$',
        'layers': (12, 6),
        'alpha': 0.05,
        'cor': '#9467bd',       # Roxo / Púrpura
        'estilo': ':',
        'largura': 2.3
    }
]


def gerar_comparativo_arquiteturas():
    print("=" * 90)
    print("GERADOR DE COMPARATIVO VISUAL DE ARQUITETURAS MLP TANH (16 ENSAIOS)")
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

    # 3. Treina cada arquitetura sob Leave-One-Group-Out
    resultados_arqs = {}
    for arq in ARQUITETURAS:
        print(f"\n-> Simulando arquitetura: {arq['nome']}...")
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('mlp', MLPRegressor(
                hidden_layer_sizes=arq['layers'],
                activation='tanh',
                solver='lbfgs',
                alpha=arq['alpha'],
                max_iter=3000,
                random_state=42
            ))
        ])

        sim_ensaios = {}
        y_real_all, y_pred_all = [], []

        for fold_idx, (train_idx, val_idx) in enumerate(logo.split(X, y, grupos)):
            ens_id = int(np.unique(grupos[val_idx])[0])
            b = sim_base[ens_id]

            m_fold = clone(pipe)
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

            # Predição discreta
            df_pts = df_base[df_base['ensaio_id'] == ens_id].copy()
            Xp, _, _ = preparar_features_e_alvo(df_pts)
            delta_p = m_fold.predict(Xp)
            x_pts = acoplador.acoplar_e_restringir(df_pts, delta_p)

            y_real_all.extend(b['x_exp'])
            y_pred_all.extend(x_pts)

            sim_ensaios[ens_id] = {
                'x_cont': x_cont,
                'x_pts': x_pts
            }

        m_glob = calcular_metricas(np.array(y_real_all), np.array(y_pred_all))
        print(f"   Métricas LOGO: R² = {m_glob['R2']:.4f} | RMSE = {m_glob['RMSE']:.4f} | MAE = {m_glob['MAE']*100:.2f}%")

        resultados_arqs[arq['id']] = {
            'meta': arq,
            'metricas': m_glob,
            'sim_ensaios': sim_ensaios
        }

    # 4. Renderização da Grade 4x4 com as 4 curvas sobrepostas
    print("\n-> Renderizando Grade 4x4 com as 4 arquiteturas sobrepostas...")
    fig, axes = plt.subplots(4, 4, figsize=(19, 14.5), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ens_id in range(1, 17):
        ax = axes_flat[ens_id - 1]
        b = sim_base[ens_id]

        # White-Box
        ax.plot(t_cont, b['x_pbm_cont'], color='#cc0000', ls=':', lw=1.3, alpha=0.65, label='White-Box (PBM)')

        # Plotar cada uma das 4 arquiteturas
        for arq in ARQUITETURAS:
            r = resultados_arqs[arq['id']]
            x_c = r['sim_ensaios'][ens_id]['x_cont']
            ax.plot(t_cont, x_c, color=arq['cor'], ls=arq['estilo'], lw=arq['largura'], alpha=0.9,
                    label=arq['rotulo_legenda'])

        # Pontos experimentais reais
        ax.plot(b['tempos_exp'], b['x_exp'], 'ko', ms=6.5, mfc='#ffff66', mec='black', mew=1.2,
                zorder=10, label='Experimental ($X_{Zn}^{exp}$)')

        # Teto estequiométrico se eta < 1
        if b['eta'] < 1.0:
            ax.axhline(b['eta'], color='#666666', ls='-.', lw=0.9, alpha=0.5)

        ax.set_title(f"Ensaio {ens_id}: $\\eta = {b['eta']}$, $C_{{A0}} = {b['ca0']}$ mol/L",
                     fontsize=9.5, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)

    for i in range(4):
        axes[3, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    # Monta legenda superior completa
    handles, labels = axes[0, 0].get_legend_handles_labels()
    # Remove duplicatas mantendo ordem
    by_label = dict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), loc='upper center', bbox_to_anchor=(0.5, 0.995),
               ncol=3, frameon=True, fontsize=10.0, facecolor='#fafafa')

    # Subtítulo com resumo estatístico comparativo
    subtitulo_stats = (
        "Comparativo de Curvatura e Desempenho Visual entre Arquiteturas MLP (Tanh + L-BFGS) em Teste Cego (LOGO-CV)\n"
        f"Original (12, 6, α=0.10): R²={resultados_arqs['orig']['metricas']['R2']:.4f}, MAE={resultados_arqs['orig']['metricas']['MAE']*100:.2f}%  |  "
        f"Rasa (12,, α=0.03): R²={resultados_arqs['rasa_12']['metricas']['R2']:.4f}, MAE={resultados_arqs['rasa_12']['metricas']['MAE']*100:.2f}%  |  "
        f"Compacta (10, 5, α=0.05): R²={resultados_arqs['comp_10_5']['metricas']['R2']:.4f}, MAE={resultados_arqs['comp_10_5']['metricas']['MAE']*100:.2f}%  |  "
        f"Calibrada (12, 6, α=0.05): R²={resultados_arqs['ajust_12_6']['metricas']['R2']:.4f}, MAE={resultados_arqs['ajust_12_6']['metricas']['MAE']*100:.2f}%"
    )
    plt.suptitle(subtitulo_stats, fontsize=11.5, fontweight='bold', y=1.035)
    plt.tight_layout(rect=[0, 0, 1, 0.94])

    nome_arquivo = "comparativo_4_arquiteturas_mlp_16_ensaios.png"
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

    print("\n" + "=" * 90)
    print("CONCLUÍDO! RESUMO COMPARATIVO DAS 4 ARQUITETURAS (LOGO-CV):")
    print("=" * 90)
    for arq in ARQUITETURAS:
        m = resultados_arqs[arq['id']]['metricas']
        print(f"{arq['nome']:<40} | R² = {m['R2']:.4f} | RMSE = {m['RMSE']:.4f} | MAE = {m['MAE']*100:.2f}%")
    print("=" * 90)


if __name__ == "__main__":
    gerar_comparativo_arquiteturas()
