"""
plotar_cinetica.py
------------------
Plota as curvas cinéticas de conversão de zinco (X_Zn) em função do tempo (t),
comparando a conversão gerada pelos modelos (White-Box Puro e Grey-Box Híbrido)
com os pontos experimentais das Tabelas A1.1 e A1.4 de Bortot Coelho (2017).

Gera visualizações científicas de alta qualidade:
  1. Grade 4x4 completa: todos os 16 ensaios individuais de bancada com pontos experimentais.
  2. Painel agrupado 2x2 por razão estequiométrica (eta = 0.5, 1.0, 1.5, 3.1) variando CA0.
  3. Gráfico ponto a ponto sequencial dos 128 pontos experimentais com diagrama de resíduos.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Adiciona o diretório scripts ao path
diretorio_scripts = os.path.dirname(os.path.abspath(__file__))
if diretorio_scripts not in sys.path:
    sys.path.append(diretorio_scripts)

from data_loader import carregar_dados_bancada, obter_caminhos
from population_balance import PopulationBalanceModel
from kinetics import ALPHA_TESE_UM_MIN
from hybrid_model import ModeloHibridoCinza, AcopladorFisico, calcular_metricas


def obter_slug_modelo(nome: str) -> str:
    mapeamento = {
        'Gradient Boosting (GBDT)': 'gbdt',
        'Extra Trees': 'extra_trees',
        'Random Forest': 'random_forest',
        'SVR (Otimizado)': 'svr',
        'Gaussian Process (Kriging)': 'gaussian_process',
        'Rede Neural (MLP Tanh + L-BFGS)': 'mlp_tanh',
        'SVR (Baseline)': 'svr_baseline',
        'Rede Neural (MLP ReLU + Adam)': 'mlp_relu'
    }
    return mapeamento.get(nome, 'modelo')


def gerar_graficos_cinetica(pasta_saida: str = None, algoritmo_hibrido: str = "SVR (Otimizado)"):
    if pasta_saida is None:
        pasta_saida = os.path.join(diretorio_scripts, "docs")
    os.makedirs(pasta_saida, exist_ok=True)
    slug = obter_slug_modelo(algoritmo_hibrido)

    print("=" * 80)
    print("GERADOR DE GRÁFICOS CINÉTICOS: CONVERSÃO GERADA VS PONTOS EXPERIMENTAIS")
    print(f"Modelo Híbrido selecionado: {algoritmo_hibrido} (Slug: {slug})")
    print("=" * 80)

    # 1. Carrega dados de bancada
    df_exp = carregar_dados_bancada()
    
    # 2. Inicializa e treina o modelo híbrido
    hibrido = ModeloHibridoCinza(nome_algoritmo=algoritmo_hibrido)
    df_treino = hibrido.preparar_dados_base(df_exp)
    hibrido.treinar(df_treino)
    
    pbm = PopulationBalanceModel()
    t_continuo = np.linspace(0.0, 15.0, 301)  # Malha fina de tempo de 0 a 15 min

    # Dicionário para armazenar simulações de cada ensaio
    simulacoes = {}
    
    print("\n-> Simulando cinética contínua para os 16 ensaios de bancada...")
    for ens_id, grp in df_exp.groupby('ensaio_id', sort=True):
        eta = float(grp['razao_molar_eta'].iloc[0])
        eta_real = float(grp['razao_molar_eta_real'].iloc[0])  # Eq. 4.3
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        razao_sl = float(grp['razao_SL_g_L'].iloc[0])
        tempos_exp = grp['tempo_min'].values
        x_exp = grp['X_zn_exp'].values
        
        # Simulação White-Box no grid contínuo (alpha = 0)
        res_pbm_cont = pbm.simular_ensaio(
            tempos_min=t_continuo.tolist(),
            C_acid_0_mol_L=ca0,
            eta=eta_real
        )
        x_pbm_cont = res_pbm_cont['X_pbm']
        c_acid_pbm_cont = res_pbm_cont['C_acid_mol_L']
        res_pbm_emp = PopulationBalanceModel(formulacao='tese').simular_ensaio(tempos_min=t_continuo.tolist(), C_acid_0_mol_L=ca0, eta=eta_real,
                                       alpha=ALPHA_TESE_UM_MIN)
        x_pbm_emp = np.array(res_pbm_emp['X_pbm'])  # modelo completo da tese (Eq. 4.1)
        
        # Predição Grey-Box no grid contínuo
        df_sim_cont = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': t_continuo,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': x_pbm_cont,
            'C_acid_pbm_mol_L': c_acid_pbm_cont
        })
        x_grey_cont = hibrido.prever(df_sim_cont)
        
        # Predição nos pontos experimentais exatos para cálculo de métricas
        df_exp_pts = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': tempos_exp,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': grp['X_pbm'].values if 'X_pbm' in grp.columns else pbm.simular_ensaio(tempos_exp.tolist(), ca0, eta_real)['X_pbm'],
            'C_acid_pbm_mol_L': ca0 * np.maximum(0.0, 1.0 - ((grp['X_pbm'].values if 'X_pbm' in grp.columns else pbm.simular_ensaio(tempos_exp.tolist(), ca0, eta_real)['X_pbm']) / eta_real))
        })
        x_grey_pts = hibrido.prever(df_exp_pts)
        
        m_pbm = calcular_metricas(x_exp, df_exp_pts['X_pbm'].values)
        m_grey = calcular_metricas(x_exp, x_grey_pts)
        
        simulacoes[ens_id] = {
            'eta': eta,
            'ca0': ca0,
            'razao_sl': razao_sl,
            'tempos_exp': tempos_exp,
            'x_exp': x_exp,
            't_cont': t_continuo,
            'x_pbm_cont': x_pbm_cont,
            'x_pbm_emp': x_pbm_emp,
            'x_grey_cont': x_grey_cont,
            'x_grey_pts': x_grey_pts,
            'm_pbm': m_pbm,
            'm_grey': m_grey
        }

    # =========================================================================
    # FIGURA 1: Grade 4x4 Completa (16 Ensaios Individuais)
    # =========================================================================
    print("-> Gerando Figura 1: Grade 4x4 dos 16 ensaios de bancada...")
    fig, axes = plt.subplots(4, 4, figsize=(18, 14), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for ens_id in range(1, 17):
        ax = axes_flat[ens_id - 1]
        s = simulacoes[ens_id]
        
        # 1. Curva White-Box (alpha = 0) e referência da tese (alpha = 5500)
        ax.plot(s['t_cont'], s['x_pbm_cont'], 'r--', lw=1.6, alpha=0.85, 
                label='White-Box (PBM, $\\alpha = 0$)')
        ax.plot(s['t_cont'], s['x_pbm_emp'], color='gray', ls=':', lw=1.2, alpha=0.7, label='Tese 2017 ($\\alpha = 5500$)')
        
        # 3. Curva Grey-Box Híbrida
        ax.plot(s['t_cont'], s['x_grey_cont'], color='#1f77b4', lw=2.2,
                label=f'Grey-Box [{algoritmo_hibrido}]')
        
        # 4. Pontos Experimentais Reais
        ax.plot(s['tempos_exp'], s['x_exp'], 'ko', ms=6, mfc='#2ca02c', mec='black', mew=1.0,
                zorder=10, label='Experimental ($X_{Zn}^{exp}$)')

        # Limite estequiométrico quando eta < 1
        if s['eta'] < 1.0:
            ax.axhline(s['eta'], color='purple', ls='-.', lw=1.0, alpha=0.6,
                       label=f'Teto $\\eta = {s["eta"]}$')
            
        ax.set_title(f"Ensaio {ens_id}: $\\eta = {s['eta']}$, $C_{{A0}} = {s['ca0']}$ mol/L",
                     fontsize=10, fontweight='bold', pad=4)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)
        
        # Caixa de texto com métricas do ensaio
        texto_metrica = f"$R^2_{{GB}} = {s['m_grey']['R2']:.3f}$\n$MAE_{{GB}} = {s['m_grey']['MAE']*100:.1f}\\%$"
        ax.text(0.96, 0.06, texto_metrica, transform=ax.transAxes,
                fontsize=7.5, verticalalignment='bottom', horizontalalignment='right',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.85, edgecolor='#cccccc'))

    # Rótulos dos eixos externos
    for i in range(4):
        axes[3, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    # Legenda global única no topo
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.995),
               ncol=5, frameon=True, fontsize=10.5, facecolor='#fafafa')

    plt.suptitle(f'Cinética de Lixiviação de Zinco em Reator Batelada: White-Box vs Grey-Box [{algoritmo_hibrido}]',
                 fontsize=14, fontweight='bold', y=1.015)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    caminho_fig1 = os.path.join(pasta_saida, f"cinetica_16_ensaios_{slug}.png")
    plt.savefig(caminho_fig1, bbox_inches='tight')
    plt.close()
    print(f"-> Figura 1 salva em: {caminho_fig1}")

    # =========================================================================
    # FIGURA 2: Curvas Agrupadas por Razão Estequiométrica (Painel 2x2 por eta)
    # =========================================================================
    print("-> Gerando Figura 2: Painel 2x2 agrupado por razão estequiométrica eta...")
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 11), dpi=300, sharex=True, sharey=True)
    etas_unicos = [0.5, 1.0, 1.5, 3.1]
    
    cores_ca0 = {
        0.1: '#1f77b4',   # Azul
        0.5: '#ff7f0e',   # Laranja
        1.0: '#2ca02c',   # Verde
        1.5: '#d62728'    # Vermelho
    }
    simbolos_ca0 = {
        0.1: 'o',
        0.5: 's',
        1.0: '^',
        1.5: 'D'
    }

    subtitulos = {
        0.5: r"(A) Deficiência Ácida ($\eta = 0{,}5$) -- Estagnação Estequiométrica",
        1.0: r"(B) Equivalência ($\eta = 1{,}0$) -- Equilíbrio Ácido/Base",
        1.5: r"(C) Leve Excesso ($\eta = 1{,}5$) -- Cinética Rápida",
        3.1: r"(D) Grande Excesso ($\eta = 3{,}1$) -- Dissolução Quase Total"
    }

    for idx_eta, eta_val in enumerate(etas_unicos):
        r = idx_eta // 2
        c = idx_eta % 2
        ax = axes2[r, c]
        
        # Filtra os 4 ensaios com este eta
        for ens_id, s in simulacoes.items():
            if abs(s['eta'] - eta_val) < 1e-4:
                ca0 = s['ca0']
                cor = cores_ca0[ca0]
                mkr = simbolos_ca0[ca0]
                
                # Curva contínua do modelo híbrido
                ax.plot(s['t_cont'], s['x_grey_cont'], color=cor, lw=2.2, alpha=0.9,
                        label=f'$C_{{A0}} = {ca0}$ mol/L (Grey-Box)')
                
                # Curva White-Box tracejada fina
                ax.plot(s['t_cont'], s['x_pbm_cont'], color=cor, ls='--', lw=1.2, alpha=0.5)
                
                # Pontos experimentais
                ax.plot(s['tempos_exp'], s['x_exp'], marker=mkr, color=cor, ls='', ms=7.5,
                        mfc='white', mew=2.0, zorder=8,
                        label=f'$C_{{A0}} = {ca0}$ mol/L (Exp.)')
        
        # Linha de teto quando eta < 1
        if eta_val < 1.0:
            ax.axhline(eta_val, color='black', ls='-.', lw=1.3, alpha=0.7,
                       label=f'Limite Teórico $\\eta = {eta_val}$')

        ax.set_title(subtitulos[eta_val], fontsize=11, fontweight='bold', pad=6)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(-0.2, 15.3)
        ax.set_ylim(-0.02, 1.05)
        ax.legend(loc='lower right' if eta_val >= 1.0 else 'upper left', 
                  fontsize=8.0, framealpha=0.9, ncol=2)

    for i in range(2):
        axes2[1, i].set_xlabel('Tempo de Lixiviação $t$ (min)', fontsize=11, fontweight='medium')
        axes2[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    plt.suptitle(r'Influência da Concentração Inicial de Ácido ($C_{A0}$) e Razão Estequiométrica ($\eta$):' + '\n' +
                 f'Comparação Ponto a Ponto entre Dados Experimentais e Modelo Grey-Box [{algoritmo_hibrido}]',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout()

    caminho_fig2 = os.path.join(pasta_saida, f"cinetica_agrupada_por_eta_{slug}.png")
    plt.savefig(caminho_fig2, bbox_inches='tight')
    plt.close()
    print(f"-> Figura 2 salva em: {caminho_fig2}")

    # =========================================================================
    # FIGURA 3: Comparação Ponto a Ponto Sequencial dos 128 Pontos Experimentais
    # =========================================================================
    print("-> Gerando Figura 3: Gráfico Ponto a Ponto sequencial dos 128 pontos com Resíduos...")
    
    # Monta dataframe ordenado por ensaio e tempo
    df_exp_sort = df_exp.sort_values(by=['ensaio_id', 'tempo_min']).reset_index(drop=True)
    X_exp_pts = df_exp_sort['X_zn_exp'].values
    
    # Calcula previsões exatas
    df_pred = df_treino.sort_values(by=['ensaio_id', 'tempo_min']).reset_index(drop=True)
    X_pbm_pts = df_pred['X_pbm'].values
    X_grey_pts = hibrido.prever(df_pred)
    
    indices = np.arange(1, len(df_exp_sort) + 1)
    
    fig3, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(16, 8.5), dpi=300, 
                                          sharex=True, gridspec_kw={'height_ratios': [2.2, 1.2]})
    
    # Painel Superior: Conversão Ponto a Ponto
    ax_top.plot(indices, X_exp_pts, 'ko-', lw=1.2, ms=5, mfc='black', 
                label='Experimental ($X_{Zn}^{exp}$)')
    ax_top.plot(indices, X_pbm_pts, 'r.--', lw=1.0, ms=4, alpha=0.8,
                label='White-Box Puro ($X_{PBM}$)')
    ax_top.plot(indices, X_grey_pts, 'b.-', lw=1.4, ms=4.5, alpha=0.9,
                label=f'Grey-Box [{algoritmo_hibrido}]')
    
    # Linhas divisórias entre ensaios
    for e in range(1, 16):
        x_div = e * 8 + 0.5
        ax_top.axvline(x_div, color='#dddddd', ls=':', lw=1.0)
        ax_bot.axvline(x_div, color='#dddddd', ls=':', lw=1.0)
        
    ax_top.set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')
    ax_top.set_title(f'Comparação Ponto a Ponto nos 128 Pontos Experimentais: Grey-Box [{algoritmo_hibrido}]',
                     fontsize=13, fontweight='bold', pad=8)
    ax_top.set_ylim(-0.05, 1.08)
    ax_top.grid(True, linestyle='--', alpha=0.4)
    ax_top.legend(loc='upper left', frameon=True, fontsize=10, facecolor='#fafafa')

    # Adiciona rótulos dos ensaios no topo
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
    ax_bot.plot(indices, res_grey, 'b.-', lw=1.2, ms=4.5, alpha=0.9,
                label=f'Resíduo Grey-Box [{algoritmo_hibrido}] (MAE = {np.mean(np.abs(res_grey))*100:.1f}%)')
    
    ax_bot.set_xlabel('Índice Sequencial do Ponto Experimental (1 a 128)', fontsize=11, fontweight='medium')
    ax_bot.set_ylabel('Erro $\\Delta X = X_{exp} - X_{prev}$', fontsize=10, fontweight='medium')
    ax_bot.set_ylim(-0.25, 0.45)
    ax_bot.set_xlim(0, 129)
    ax_bot.grid(True, linestyle='--', alpha=0.4)
    ax_bot.legend(loc='upper right', frameon=True, fontsize=9.0, facecolor='#fafafa')

    plt.tight_layout()
    caminho_fig3 = os.path.join(pasta_saida, f"comparacao_ponto_a_ponto_128_{slug}.png")
    plt.savefig(caminho_fig3, bbox_inches='tight')
    plt.close()
    print(f"-> Figura 3 salva em: {caminho_fig3}")

    print("\nTodos os gráficos foram gerados e salvos com sucesso!")
    return {
        'figura_16_ensaios': caminho_fig1,
        'figura_agrupada_eta': caminho_fig2,
        'figura_ponto_a_ponto': caminho_fig3
    }


def plotar_comparacao_suavidade(pasta_saida: str = None):
    """
    Gera um gráfico comparativo direto da suavidade das curvas:
    Árvores (Gradient Boosting) vs Métodos de Kernel Suave (SVR Otimizado com Kernel RBF)
    evidenciando a ausência de degraus na cinética contínua.
    """
    if pasta_saida is None:
        pasta_saida = os.path.join(diretorio_scripts, "docs")
    os.makedirs(pasta_saida, exist_ok=True)
    
    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    t_cont = np.linspace(0.0, 15.0, 400)
    
    # Modelos a comparar
    h_gbdt = ModeloHibridoCinza('Gradient Boosting (GBDT)')
    dft_gbdt = h_gbdt.preparar_dados_base(df_exp)
    h_gbdt.treinar(dft_gbdt)
    
    h_svr = ModeloHibridoCinza('SVR (Otimizado)')
    dft_svr = h_svr.preparar_dados_base(df_exp)
    h_svr.treinar(dft_svr)

    h_gp = ModeloHibridoCinza('Gaussian Process (Kriging)')
    dft_gp = h_gp.preparar_dados_base(df_exp)
    h_gp.treinar(dft_gp)
    
    # Ensaios representativos dos 4 regimes:
    # Ensaio 1: eta=0.5, CA0=0.1
    # Ensaio 6: eta=1.0, CA0=0.1
    # Ensaio 8: eta=1.0, CA0=0.5 (planta piloto)
    # Ensaio 11: eta=1.5, CA0=0.1
    ensaios_alvo = [1, 6, 8, 11]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300, sharex=True, sharey=True)
    axes_flat = axes.flatten()
    
    for idx, ens_id in enumerate(ensaios_alvo):
        ax = axes_flat[idx]
        grp = df_exp[df_exp['ensaio_id'] == ens_id]
        eta = float(grp['razao_molar_eta'].iloc[0])
        eta_real = float(grp['razao_molar_eta_real'].iloc[0])  # Eq. 4.3
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        razao_sl = float(grp['razao_SL_g_L'].iloc[0])
        tempos_exp = grp['tempo_min'].values
        x_exp = grp['X_zn_exp'].values
        
        # Simulação White-Box contínua
        res_pbm = pbm.simular_ensaio(t_cont.tolist(), ca0, eta_real)
        
        df_sim = pd.DataFrame({
            'ensaio_id': ens_id,
            'tempo_min': t_cont,
            'razao_molar_eta': eta,
            'C_acid_0_mol_L': ca0,
            'razao_SL_g_L': razao_sl,
            'X_pbm': res_pbm['X_pbm'],
            'C_acid_pbm_mol_L': res_pbm['C_acid_mol_L']
        })
        
        y_gbdt = h_gbdt.prever(df_sim)
        y_svr = h_svr.prever(df_sim)
        y_gp = h_gp.prever(df_sim)
        
        # 1. White-Box Puro
        ax.plot(t_cont, res_pbm['X_pbm'], 'r--', lw=1.5, alpha=0.7, label='White-Box (PBM, $\\alpha = 0$)')
        
        # 2. Gradient Boosting (Degraus de Árvore)
        ax.plot(t_cont, y_gbdt, color='#ff7f0e', ls='-', lw=2.0, alpha=0.85, 
                label='Grey-Box [GBDT] (Árvores - degraus)')
        
        # 3. SVR Otimizado (RBF Kernel - Suave C_infinito)
        ax.plot(t_cont, y_svr, color='#1f77b4', ls='-', lw=2.4, alpha=0.95, 
                label='Grey-Box [SVR Otimizado] (Kernel RBF Suave)')
        
        # 4. Gaussian Process
        ax.plot(t_cont, y_gp, color='#2ca02c', ls=':', lw=2.0, alpha=0.9, 
                label='Grey-Box [Processo Gaussiano]')
        
        # 5. Pontos Experimentais
        ax.plot(tempos_exp, x_exp, 'ko', ms=7, mfc='black', zorder=10, label='Dados Experimentais')
        
        if eta < 1.0:
            ax.axhline(eta, color='purple', ls='-.', lw=1.0, alpha=0.6, label=f'Teto $\\eta = {eta}$')
            
        ax.set_title(f"Ensaio {ens_id}: $\\eta = {eta}$, $C_{{A0}} = {ca0}$ mol/L",
                     fontsize=11, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.set_xlim(-0.2, 15.3)
        ax.set_ylim(-0.02, 1.05)
        ax.legend(loc='lower right' if eta >= 1.0 else 'upper left', fontsize=8.0, framealpha=0.9)

    for i in range(2):
        axes[1, i].set_xlabel('Tempo de Lixiviação (min)', fontsize=11, fontweight='medium')
        axes[i, 0].set_ylabel('Conversão de Zinco ($X_{Zn}$)', fontsize=11, fontweight='medium')

    plt.suptitle('Comparação de Suavidade Cinética: Modelos de Árvore (GBDT) vs Modelos de Kernel Contínuo (SVR e Kriging)',
                 fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    caminho_comp = os.path.join(pasta_saida, "comparacao_suavidade_arvores_vs_kernel.png")
    plt.savefig(caminho_comp, bbox_inches='tight')
    plt.close()
    print(f"-> Gráfico comparativo de suavidade salvo em: {caminho_comp}")
    return caminho_comp


if __name__ == "__main__":
    import sys
    algoritmo = sys.argv[1] if len(sys.argv) > 1 else "SVR (Otimizado)"
    gerar_graficos_cinetica(algoritmo_hibrido=algoritmo)
    print("\n-> Gerando comparação visual de suavidade entre GBDT, SVR e Kriging...")
    plotar_comparacao_suavidade()
