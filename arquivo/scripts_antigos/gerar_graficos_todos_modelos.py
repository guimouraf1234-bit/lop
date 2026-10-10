"""
gerar_graficos_todos_modelos.py
--------------------------------
Gera o conjunto completo de visualizações cinéticas científicas de alta resolução
para TODOS os 8 algoritmos de Machine Learning catalogados em hybrid_model.py,
adotando ESTRITAMENTE a validação cruzada agrupada em TESTE CEGO (GroupKFold com 8 dobras / 85-15).

Para cada modelo, produz:
  1. Grade 4x4 Completa (cinetica_16_ensaios_{slug}.png):
     - Curva White-Box (PBM da tese de Bortot Coelho 2017, alpha=5500)
     - Curva Grey-Box gerada 100% em Teste Cego Out-of-Fold (treino em 14 ensaios, teste no ensaio plotado)
     - Pontos Experimentais Reais (X_zn_exp)
     - R² e MAE de teste cego em cada quadrante.

  2. Comparação Ponto a Ponto Sequencial dos 128 Pontos (comparacao_ponto_a_ponto_128_{slug}.png):
     - Painel Superior: Trajetória discreta ponto a ponto dos 128 pontos ordenados por ensaio (1 a 16)
     - Painel Inferior: Diagrama de resíduos Delta_X em teste cego vs faixa de tolerância de +-5%.

Modelos avaliados:
  1. Extra Trees Regressor
  2. Gradient Boosting (GBDT)
  3. Random Forest Regressor
  4. Gaussian Process (Kriging)
  5. SVR (Otimizado) [Kernel RBF]
  6. SVR (Baseline)
  7. Rede Neural (MLP Tanh + L-BFGS)
  8. Rede Neural (MLP ReLU + Adam)
"""

import os
import sys
import glob
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import GroupKFold
from sklearn.base import clone

# Adiciona o diretório dos scripts ao path
diretorio_scripts = os.path.dirname(os.path.abspath(__file__))
if diretorio_scripts not in sys.path:
    sys.path.append(diretorio_scripts)

from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel
from kinetics import ALPHA_TESE_UM_MIN
from hybrid_model import (
    ModeloHibridoCinza,
    preparar_features_e_alvo,
    obter_catalogo_modelos,
    calcular_metricas,
    AcopladorFisico
)

MODELOS_AVALIADOS = [
    {
        'nome': 'Extra Trees',
        'slug': 'extra_trees',
        'titulo': 'Extra Trees Regressor'
    },
    {
        'nome': 'Gradient Boosting (GBDT)',
        'slug': 'gradient_boosting',
        'slug_alt': 'gbdt',
        'titulo': 'Gradient Boosting (GBDT)'
    },
    {
        'nome': 'Random Forest',
        'slug': 'random_forest',
        'titulo': 'Random Forest Regressor'
    },
    {
        'nome': 'Gaussian Process (Kriging)',
        'slug': 'gaussian_process',
        'slug_alt': 'kriging',
        'titulo': 'Gaussian Process Regressor (Kriging)'
    },
    {
        'nome': 'SVR (Otimizado)',
        'slug': 'svr_otimizado',
        'slug_alt': 'svr',
        'titulo': 'SVR (Otimizado - Kernel RBF)'
    },
    {
        'nome': 'SVR (Baseline)',
        'slug': 'svr_baseline',
        'titulo': 'SVR (Baseline - C=1.0, eps=0.01)'
    },
    {
        'nome': 'Rede Neural (MLP Tanh)',
        'slug': 'mlp_tanh',
        'titulo': 'Rede Neural (MLP Tanh + L-BFGS)'
    },
    {
        'nome': 'Rede Neural (MLP ReLU + Adam)',
        'slug': 'mlp_relu',
        'titulo': 'Rede Neural (MLP ReLU + Adam)'
    }
]


def limpar_diretorios_graficos(pastas: list):
    """Remove gráficos PNG anteriores para garantir que todos os arquivos sejam gerados em teste cego."""
    print("-> Limpando pastas de gráficos anteriores...")
    total_removidos = 0
    for pasta in pastas:
        if os.path.exists(pasta):
            arquivos = glob.glob(os.path.join(pasta, "*.png"))
            for f in arquivos:
                try:
                    os.remove(f)
                    total_removidos += 1
                except Exception as e:
                    print(f"   Aviso ao remover {f}: {e}")
    print(f"   Total de {total_removidos} arquivos antigos removidos com sucesso.")


def precomputar_simulacoes_base(df_exp: pd.DataFrame, pbm: PopulationBalanceModel, t_cont: np.ndarray):
    """
    Pré-computa as simulações mecanicistas White-Box (modelo da tese, alpha=5500) para os 16 ensaios,
    otimizando a resolução numérica do PBM para todos os modelos.
    """
    dados_base = {}
    print("\n-> Pré-computando curvas White-Box do PBM para os 16 ensaios de bancada...")
    
    for ens_id, grp in df_exp.groupby('ensaio_id', sort=True):
        eta = float(grp['razao_molar_eta'].iloc[0])
        eta_real = float(grp['razao_molar_eta_real'].iloc[0])  # Eq. 4.3
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        razao_sl = float(grp['razao_SL_g_L'].iloc[0])
        tempos_exp = grp['tempo_min'].values
        x_exp = grp['X_zn_exp'].values

        res_pbm_cont = pbm.simular_ensaio(
            tempos_min=t_cont.tolist(),
            C_acid_0_mol_L=ca0,
            eta=eta_real
        )
        x_pbm_cont = np.array(res_pbm_cont['X_pbm'])
        c_acid_pbm_cont = np.array(res_pbm_cont['C_acid_mol_L'])

        res_pbm_emp = PopulationBalanceModel(formulacao='tese').simular_ensaio(tempos_min=t_cont.tolist(), C_acid_0_mol_L=ca0, eta=eta_real,
                                       alpha=ALPHA_TESE_UM_MIN)
        x_pbm_emp = np.array(res_pbm_emp['X_pbm'])  # modelo completo da tese (Eq. 4.1)

        dados_base[ens_id] = {
            'eta': eta,
            'ca0': ca0,
            'razao_sl': razao_sl,
            'tempos_exp': tempos_exp,
            'x_exp': x_exp,
            't_cont': t_cont,
            'x_pbm_cont': x_pbm_cont,
            'c_acid_pbm_cont': c_acid_pbm_cont,
            'x_pbm_emp': x_pbm_emp
        }

    return dados_base


def gerar_graficos_modelo_teste_cego(
    nome_modelo: str,
    titulo_modelo: str,
    slug_modelo: str,
    df_exp: pd.DataFrame,
    df_base: pd.DataFrame,
    X: np.ndarray,
    y: np.ndarray,
    grupos: np.ndarray,
    dados_base: dict,
    catalogo: dict,
    pasta_destino: str,
    slug_alt: str = None
):
    """
    Executa GroupKFold (8 dobras: 85/15) e plota:
      1. cinetica_16_ensaios_{slug}.png
      2. comparacao_ponto_a_ponto_128_{slug}.png
    """
    print(f"\n[{nome_modelo}] Executando GroupKFold (8 dobras) e gerando predições cegas...")
    modelo_base = catalogo[nome_modelo]
    acoplador = AcopladorFisico()
    gkf = GroupKFold(n_splits=8)

    simulacoes_oof = {}
    y_real_global = []
    y_pred_global = []

    for fold_idx, (train_idx, val_idx) in enumerate(gkf.split(X, y, grupos)):
        ens_teste = np.unique(grupos[val_idx])

        # Treina estritamente com os 14 ensaios da partição de treino
        modelo_fold = clone(modelo_base)
        modelo_fold.fit(X[train_idx], y[train_idx])

        for ens_id in ens_teste:
            b = dados_base[ens_id]

            # 1. Predição discreta nos pontos experimentais reais (teste cego)
            df_pts = df_base[df_base['ensaio_id'] == ens_id].copy()
            X_pts, _, _ = preparar_features_e_alvo(df_pts)
            delta_pts = modelo_fold.predict(X_pts)
            x_grey_pts = acoplador.acoplar_e_restringir(df_pts, delta_pts)

            m_ens = calcular_metricas(b['x_exp'], x_grey_pts)
            y_real_global.extend(b['x_exp'])
            y_pred_global.extend(x_grey_pts)

            # 2. Predição contínua na grade fina de tempo (trajetória suave de teste cego)
            df_sim_cont = pd.DataFrame({
                'ensaio_id': ens_id,
                'tempo_min': b['t_cont'],
                'razao_molar_eta': b['eta'],
                'C_acid_0_mol_L': b['ca0'],
                'razao_SL_g_L': b['razao_sl'],
                'X_pbm': b['x_pbm_cont'],
                'C_acid_pbm_mol_L': b['c_acid_pbm_cont']
            })
            X_cont, _, _ = preparar_features_e_alvo(df_sim_cont)
            delta_cont = modelo_fold.predict(X_cont)
            x_grey_cont = acoplador.acoplar_e_restringir(df_sim_cont, delta_cont)

            simulacoes_oof[ens_id] = {
                'eta': b['eta'],
                'ca0': b['ca0'],
                'razao_sl': b['razao_sl'],
                'tempos_exp': b['tempos_exp'],
                'x_exp': b['x_exp'],
                't_cont': b['t_cont'],
                'x_pbm_cont': b['x_pbm_cont'],
                'x_pbm_emp': b['x_pbm_emp'],
                'x_grey_cont': x_grey_cont,
                'x_grey_pts': x_grey_pts,
                'm_grey': m_ens,
                'fold': fold_idx + 1
            }

    m_global = calcular_metricas(np.array(y_real_global), np.array(y_pred_global))
    print(f"  -> R²_CV: {m_global['R2']:.4f} | RMSE_CV: {m_global['RMSE']:.4f} | MAE_CV: {m_global['MAE']*100:.2f}% | Violação: {m_global['Violacao_pct']:.1f}%")

    # =========================================================================
    # FIGURA 1: Grade 4x4 dos 16 Ensaios em Teste Cego
    # =========================================================================
    fig1, axes1 = plt.subplots(4, 4, figsize=(18, 14), dpi=300, sharex=True, sharey=True)
    axes1_flat = axes1.flatten()

    for ens_id in range(1, 17):
        ax = axes1_flat[ens_id - 1]
        s = simulacoes_oof[ens_id]

        # White-Box (alpha = 0)
        ax.plot(s['t_cont'], s['x_pbm_cont'], 'r--', lw=1.6, alpha=0.85, label='White-Box (PBM, $\\alpha = 0$)')
        ax.plot(s['t_cont'], s['x_pbm_emp'], color='gray', ls=':', lw=1.2, alpha=0.7, label='Tese 2017 ($\\alpha = 5500$)')
        # Grey-Box Teste Cego
        ax.plot(s['t_cont'], s['x_grey_cont'], color='#0052cc', lw=2.2, label=f'Grey-Box [{nome_modelo}] (Teste Cego)')
        # Experimental
        ax.plot(s['tempos_exp'], s['x_exp'], 'ko', ms=6, mfc='#2ca02c', mec='black', mew=1.0, zorder=10, label='Experimental ($X_{Zn}^{exp}$)')

        if s['eta'] < 1.0:
            ax.axhline(s['eta'], color='purple', ls='-.', lw=1.0, alpha=0.6, label=f'Teto $\\eta = {s["eta"]}$')

        ax.set_title(f"Ensaio {ens_id}: $\\eta = {s['eta']}$, $C_{{A0}} = {s['ca0']}$ mol/L [Fold {s['fold']}]",
                     fontsize=9.5, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)

        texto_metrica = f"$R^2_{{teste}} = {s['m_grey']['R2']:.3f}$\n$MAE_{{teste}} = {s['m_grey']['MAE']*100:.1f}\\%$"
        ax.text(0.96, 0.06, texto_metrica, transform=ax.transAxes,
                fontsize=7.5, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='#f0f4f8', alpha=0.90, edgecolor='#8baac9'))

    for i in range(4):
        axes1[3, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes1[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    handles1, labels1 = axes1[0, 0].get_legend_handles_labels()
    fig1.legend(handles1, labels1, loc='upper center', bbox_to_anchor=(0.5, 0.995),
                ncol=5, frameon=True, fontsize=10.5, facecolor='#fafafa')

    subtitulo_global = (
        f"Cinética de Lixiviação de Zinco em Reator Batelada: White-Box vs Grey-Box [{titulo_modelo}]\n"
        f"AVALIAÇÃO EM TESTE CEGO (GroupKFold 8 Dobras / 85-15) — Métricas Globais: "
        f"$R^2_{{CV}} = {m_global['R2']:.4f}$ | $RMSE_{{CV}} = {m_global['RMSE']:.4f}$ | $MAE_{{CV}} = {m_global['MAE']*100:.2f}\\%$"
    )
    plt.suptitle(subtitulo_global, fontsize=12.5, fontweight='bold', y=1.025)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    arq_cinetica = f"cinetica_16_ensaios_{slug_modelo}.png"
    p_cinetica = os.path.join(pasta_destino, arq_cinetica)
    plt.savefig(p_cinetica, bbox_inches='tight')
    plt.close()

    if slug_alt:
        shutil.copyfile(p_cinetica, os.path.join(pasta_destino, f"cinetica_16_ensaios_{slug_alt}.png"))

    # =========================================================================
    # FIGURA 2: Comparação Ponto a Ponto Sequencial dos 128 Pontos Experimentais
    # =========================================================================
    df_exp_sort = df_exp.sort_values(by=['ensaio_id', 'tempo_min']).reset_index(drop=True)
    X_exp_pts = df_exp_sort['X_zn_exp'].values
    indices = np.arange(1, len(df_exp_sort) + 1)

    # Coleta predições ordenadas correspondentes aos 128 pontos
    X_pbm_pts_list = []
    X_grey_pts_list = []
    for ens_id in range(1, 17):
        s = simulacoes_oof[ens_id]
        X_pbm_pts_list.extend(dados_base[ens_id]['x_pbm_cont'][np.isin(np.round(dados_base[ens_id]['t_cont'], 2), np.round(s['tempos_exp'], 2))][:len(s['tempos_exp'])])
        X_grey_pts_list.extend(s['x_grey_pts'])

    # Garante consistência dimensional
    X_pbm_pts = df_base.sort_values(by=['ensaio_id', 'tempo_min'])['X_pbm'].values
    X_grey_pts = np.array(X_grey_pts_list)

    fig2, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(16, 8.5), dpi=300,
                                          sharex=True, gridspec_kw={'height_ratios': [2.2, 1.2]})

    # Painel Superior: Conversão Ponto a Ponto
    ax_top.plot(indices, X_exp_pts, 'ko-', lw=1.2, ms=5, mfc='#2ca02c', mec='black', mew=1.0,
                label='Experimental ($X_{Zn}^{exp}$)')
    ax_top.plot(indices, X_pbm_pts, 'r.--', lw=1.0, ms=4, alpha=0.75,
                label='White-Box Puro ($X_{PBM}$)')
    ax_top.plot(indices, X_grey_pts, color='#0052cc', ls='-', marker='s', lw=1.4, ms=4.5, alpha=0.9,
                label=f'Grey-Box [{titulo_modelo}] (Teste Cego)')

    for e in range(1, 16):
        x_div = e * 8 + 0.5
        ax_top.axvline(x_div, color='#dddddd', ls=':', lw=1.0)
        ax_bot.axvline(x_div, color='#dddddd', ls=':', lw=1.0)

    ax_top.set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')
    ax_top.set_title(f'Comparação Ponto a Ponto nos 128 Pontos Experimentais em Teste Cego (GroupKFold 85/15): Grey-Box [{titulo_modelo}]',
                     fontsize=12.5, fontweight='bold', pad=8)
    ax_top.set_ylim(-0.05, 1.08)
    ax_top.grid(True, linestyle='--', alpha=0.4)
    ax_top.legend(loc='upper left', frameon=True, fontsize=10, facecolor='#fafafa')

    for e in range(1, 17):
        x_centro = (e - 1) * 8 + 4.5
        ax_top.text(x_centro, 1.03, f"E{e}", fontsize=8.5, ha='center', fontweight='bold', color='#444444')

    # Painel Inferior: Resíduos ponto a ponto
    res_pbm = X_exp_pts - X_pbm_pts
    res_grey = X_exp_pts - X_grey_pts

    ax_bot.axhline(0.0, color='black', lw=1.2)
    ax_bot.axhspan(-0.05, 0.05, color='#d4edda', alpha=0.6, label='Faixa de Tolerância $\\pm 5\\%$')
    ax_bot.plot(indices, res_pbm, 'r.--', lw=0.9, ms=4, alpha=0.7,
                label=f'Resíduo White-Box (MAE = {np.mean(np.abs(res_pbm))*100:.1f}%)')
    ax_bot.plot(indices, res_grey, color='#0052cc', ls='-', marker='s', lw=1.2, ms=4.5, alpha=0.9,
                label=f'Resíduo Grey-Box [{titulo_modelo}] (MAE = {m_global["MAE"]*100:.2f}%)')

    ax_bot.set_xlabel('Índice Sequencial do Ponto Experimental (1 a 128)', fontsize=11, fontweight='medium')
    ax_bot.set_ylabel('Erro $\\Delta X = X_{exp} - X_{prev}$', fontsize=10, fontweight='medium')
    ax_bot.set_ylim(-0.25, 0.45)
    ax_bot.set_xlim(0, 129)
    ax_bot.grid(True, linestyle='--', alpha=0.4)
    ax_bot.legend(loc='upper right', frameon=True, fontsize=9.0, facecolor='#fafafa')

    plt.tight_layout()
    arq_ponto_a_ponto = f"comparacao_ponto_a_ponto_128_{slug_modelo}.png"
    p_ponto_a_ponto = os.path.join(pasta_destino, arq_ponto_a_ponto)
    plt.savefig(p_ponto_a_ponto, bbox_inches='tight')
    plt.close()

    if slug_alt:
        shutil.copyfile(p_ponto_a_ponto, os.path.join(pasta_destino, f"comparacao_ponto_a_ponto_128_{slug_alt}.png"))

    return {
        'nome': nome_modelo,
        'slug': slug_modelo,
        'r2_cv': m_global['R2'],
        'rmse_cv': m_global['RMSE'],
        'mae_cv': m_global['MAE'] * 100.0,
        'violacao': m_global['Violacao_pct'],
        'arq_cinetica': arq_cinetica,
        'arq_ponto_a_ponto': arq_ponto_a_ponto
    }


def executar():
    print("=" * 95)
    print("GERADOR CIENTÍFICO MULTI-MODELO -- VALIDAÇÃO CRUZADA ESTRITA (GROUP-KFOLD 85/15)")
    print("Todos os gráficos cinéticos serão gerados em TESTE CEGO OUT-OF-FOLD")
    print("=" * 95)

    pasta_scripts_graficos = os.path.join(diretorio_scripts, "graficos")
    pasta_scripts_graficos_acent = os.path.join(diretorio_scripts, "gráficos")
    pasta_lop_graficos = os.path.abspath(os.path.join(diretorio_scripts, "..", "graficos"))
    pasta_lop_graficos_acent = os.path.abspath(os.path.join(diretorio_scripts, "..", "gráficos"))
    pasta_docs = os.path.join(diretorio_scripts, "docs")

    pastas_alvo = [
        pasta_scripts_graficos,
        pasta_scripts_graficos_acent,
        pasta_lop_graficos,
        pasta_lop_graficos_acent,
        pasta_docs
    ]

    for p in pastas_alvo:
        os.makedirs(p, exist_ok=True)

    # 1. Limpeza de gráficos antigos
    limpar_diretorios_graficos(pastas_alvo)

    # 2. Carregamento de dados e pré-computação mecanicista
    df_exp = carregar_dados_bancada()
    hibrido_aux = ModeloHibridoCinza('Extra Trees')
    df_base = hibrido_aux.preparar_dados_base(df_exp)
    X, y, grupos = preparar_features_e_alvo(df_base)

    pbm = PopulationBalanceModel()
    t_cont = np.linspace(0.0, 15.0, 301)
    dados_base = precomputar_simulacoes_base(df_exp, pbm, t_cont)
    catalogo = obter_catalogo_modelos()

    # 3. Geração para cada modelo
    resumos = []
    for item in MODELOS_AVALIADOS:
        res = gerar_graficos_modelo_teste_cego(
            nome_modelo=item['nome'],
            titulo_modelo=item['titulo'],
            slug_modelo=item['slug'],
            df_exp=df_exp,
            df_base=df_base,
            X=X,
            y=y,
            grupos=grupos,
            dados_base=dados_base,
            catalogo=catalogo,
            pasta_destino=pasta_scripts_graficos,
            slug_alt=item.get('slug_alt')
        )
        resumos.append(res)

    # 4. Sincronização para todas as pastas de gráficos
    print("\n-> Sincronizando novos gráficos em teste cego para todas as pastas do repositório...")
    arquivos_gerados = os.listdir(pasta_scripts_graficos)
    for arq in arquivos_gerados:
        origem = os.path.join(pasta_scripts_graficos, arq)
        if os.path.isfile(origem):
            for pasta_dest in [pasta_scripts_graficos_acent, pasta_lop_graficos, pasta_lop_graficos_acent, pasta_docs]:
                destino = os.path.join(pasta_dest, arq)
                shutil.copyfile(origem, destino)

    print("\n" + "=" * 95)
    print("CONCLUÍDO COM SUCESSO! TABELA CONSOLIDADA DE TESTE CEGO (GROUP-KFOLD 85/15):")
    print("=" * 95)
    print(f"{'Modelo / Algoritmo':<35} | {'R² CV (85/15)':<14} | {'RMSE CV':<10} | {'MAE CV (%)':<12} | {'Violação':<8}")
    print("-" * 95)
    for r in resumos:
        print(f"{r['nome']:<35} | {r['r2_cv']:^14.4f} | {r['rmse_cv']:^10.4f} | {r['mae_cv']:^12.2f}% | {r['violacao']:^8.1f}%")
    print("=" * 95)
    print(f"\nTodos os gráficos ponto a ponto e grades 4x4 estão salvos em:")
    print(f" -> {pasta_lop_graficos}")
    print(f" -> {pasta_scripts_graficos}")


if __name__ == "__main__":
    executar()
