"""
gerar_curvas_cineticas_serial.py
---------------------------------
Gera gráficos científicos em alta resolução (300 DPI) para o
Modelo Híbrido Serial Cinético (alpha-GBDT -> PBM com distribuição RRB):
  1. scripts/docs/curvas_cineticas_16_ensaios_serial_85_15.png:
     Matriz completa 4x4 (16 ensaios) com pontos experimentais, curva contínua
     do Híbrido Serial, linhas de resíduo ponto a ponto e destaque para o Teste Cego (15%).
  2. scripts/docs/curvas_cineticas_ensaios_criticos_serial_85_15.png:
     Painel 2x2 focado nos 4 regimes estequiométricos (eta = 0.5, 1.0, 1.5, 3.1).
  3. scripts/docs/paridade_e_residuos_ponto_a_ponto_serial_85_15.png:
     Painel com Paridade Ponto a Ponto e Resíduos Individuais no Treino (85%) e Teste (15%).

Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Garantir saída UTF-8 no Windows
if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

# Importações dos módulos locais
from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel
from serial_hybrid_model import ModeloHibridoAlphaSerial, calcular_metricas_estatisticas


def configurar_estilo():
    """Configura o estilo matplotlib com padrão de publicação internacional."""
    plt.rcParams.update({
        'font.size': 10.0,
        'axes.labelsize': 11.0,
        'axes.titlesize': 11.5,
        'xtick.labelsize': 9.5,
        'ytick.labelsize': 9.5,
        'legend.fontsize': 9.0,
        'figure.titlesize': 13.5,
        'font.family': 'sans-serif',
        'mathtext.fontset': 'dejavusans',
        'axes.linewidth': 1.1,
        'grid.linewidth': 0.7,
        'grid.alpha': 0.4
    })


# =============================================================================
# 1. MATRIZ 4x4 COM TODOS OS 16 ENSAIOS DO MODELO SERIAL (85/15)
# =============================================================================

def gerar_matriz_16_ensaios_serial(df: pd.DataFrame,
                                   res_holdout: dict,
                                   pbm: PopulationBalanceModel,
                                   pasta_saida: str = "scripts/docs") -> str:
    """Gera matriz 4x4 (16 subplots) com curvas cinéticas ponto a ponto do modelo serial."""
    os.makedirs(pasta_saida, exist_ok=True)
    fig, axes = plt.subplots(4, 4, figsize=(18.5, 15.0), dpi=300, sharex=True, sharey=True)
    axes = axes.flatten()
    
    t_fino = np.linspace(0.0, 15.0, 150)
    ensaios_unicos = sorted(df['ensaio_id'].unique())
    ensaios_teste = set(res_holdout['ensaios_teste'])  # {1, 2, 6}
    alfas_preditos = res_holdout['alfas_preditos']
    
    for i, ens_id in enumerate(ensaios_unicos):
        ax = axes[i]
        grp = df[df['ensaio_id'] == ens_id]
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        t_exp = grp['tempo_min'].values
        y_exp = grp['X_zn_exp'].values
        
        eh_teste = ens_id in ensaios_teste
        status_tag = "TESTE CEGO (15%)" if eh_teste else "TREINO (85%)"
        cor_curva = '#d62728' if eh_teste else '#2ca02c'
        
        # 1. Simulação White-Box puro (alpha = 0)
        res_wb = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        
        # 2. Simulação Tese Coelho (2017) (alpha = 5500)
        res_coelho = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=5500.0)
        
        # 3. Predição Contínua Híbrida Serial (alpha predito pelo Blackbox)
        alpha_hat = alfas_preditos[ens_id]
        res_serial_cont = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_hat)
        x_serial_cont = res_serial_cont['X_pbm']
        
        # 4. Predição nos pontos experimentais
        res_serial_pts = pbm.simular_ensaio(tempos_min=t_exp.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_hat)
        y_serial_pts = np.array(res_serial_pts['X_pbm'])
        
        # Conexão ponto a ponto de resíduo (barras verticais)
        for tp, yr, yp in zip(t_exp, y_exp, y_serial_pts):
            ax.plot([tp, tp], [yr, yp], color='#d62728' if eh_teste else '#e65100', lw=1.1, alpha=0.75, zorder=3)
            
        # Plot experimental
        ax.scatter(t_exp, y_exp, color='black', s=45, zorder=6, edgecolors='white', lw=0.6,
                   label='Experimental ($X_{\\mathrm{Zn}}^{\\mathrm{exp}}$)')
        
        # Curvas previstas
        ax.plot(t_fino, x_serial_cont, color=cor_curva, lw=2.3, zorder=5,
                label=f'Híbrido Serial [{status_tag}]')
        ax.plot(t_fino, res_wb['X_pbm'], '--', color='#7f7f7f', lw=1.2, zorder=2,
                label=r'White-Box Puro ($\alpha=0$)')
        ax.plot(t_fino, res_coelho['X_pbm'], ':', color='#ff7f0e', lw=1.2, zorder=1,
                label=r'Tese Coelho ($\alpha=5500$)')
        
        # Teto estequiométrico quando eta < 1
        if eta < 1.0:
            ax.axhline(eta, color='gray', linestyle='-.', lw=0.8, alpha=0.6)
            
        # Métricas do ensaio
        m_ens = calcular_metricas_estatisticas(y_exp, y_serial_pts, p_parametros=2, eta=np.array([eta]*len(y_exp)))
        
        # Destaque de borda para ensaios de teste cego (15%)
        if eh_teste:
            for spine in ax.spines.values():
                spine.set_edgecolor('#d62728')
                spine.set_linewidth(2.2)
            ax.set_facecolor('#fff5f5')
        else:
            ax.set_facecolor('#ffffff')
            
        # Título
        ax.set_title(rf"Ensaio {ens_id} [{status_tag}]: $\eta={eta:.1f}$, $C_{{A0}}={ca0:.1f}$ M",
                     fontsize=9.8, fontweight='bold' if eh_teste else 'normal',
                     color='#b30000' if eh_teste else 'black')
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, linestyle='--', alpha=0.35)
        
        # Rótulos de eixos apenas nas bordas externas
        if i >= 12:
            ax.set_xlabel('Tempo de Lixiviação (min)', fontsize=10.0)
        if i % 4 == 0:
            ax.set_ylabel(r'Conversão $X_{\mathrm{Zn}}$', fontsize=10.0)
            
        # Caixa de métricas discreta
        caixa_texto = (f"$\\hat{{\\alpha}} = {alpha_hat:.0f}\\ \\mu\\mathrm{{m/min}}$\n"
                       f"$R^2 = {m_ens['R2']:.3f}$\n"
                       f"$RMSE = {m_ens['RMSE']*100:.1f}\\%$")
        ax.text(0.48, 0.12, caixa_texto, transform=ax.transAxes, fontsize=8.0,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#ffffff', edgecolor='#bbbbbb', alpha=0.9))
                
    # Legenda única superior
    handles, labels = axes[0].get_legend_handles_labels()
    # Adicionar legenda de teste e treino
    fig.legend(handles[:4], labels[:4], loc='upper center', bbox_to_anchor=(0.5, 0.995),
               ncol=4, fontsize=10.0, frameon=True, framealpha=0.95)
    
    fig.suptitle("Validação Cinética Ponto a Ponto: Modelo Híbrido Serial Cinético (alpha-GBDT -> PBM)\n"
                 "Confronto em Todos os 16 Ensaios de Bancada (Particionamento 85% Treino / 15% Teste Cego em Bordas Vermelhas)",
                 fontsize=13.0, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    caminho_fig = os.path.join(pasta_saida, "curvas_cineticas_16_ensaios_serial_85_15.png")
    plt.savefig(caminho_fig, bbox_inches='tight')
    plt.close()
    print(f"-> Painel 4x4 do modelo serial salvo em: {caminho_fig}")
    return caminho_fig


# =============================================================================
# 2. PAINEL 2x2 NOS 4 REGIMES ESTEQUIOMÉTRICOS CRÍTICOS
# =============================================================================

def gerar_painel_ensaios_criticos_serial(df: pd.DataFrame,
                                         res_holdout: dict,
                                         pbm: PopulationBalanceModel,
                                         pasta_saida: str = "scripts/docs") -> str:
    """Gera painel 2x2 focado nos 4 regimes de acidez eta (0.5, 1.0, 1.5, 3.1) em grande formato com resíduos."""
    os.makedirs(pasta_saida, exist_ok=True)
    ensaios_alvo = [1, 6, 14, 2]  # Regime 1: eta=0.5 (Teste); Regime 2: eta=1.0 (Teste); Regime 3: eta=1.5 (Treino); Regime 4: eta=3.1 (Teste)
    alfas_preditos = res_holdout['alfas_preditos']
    ensaios_teste = set(res_holdout['ensaios_teste'])
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), dpi=300)
    axes = axes.flatten()
    t_fino = np.linspace(0.0, 15.0, 150)
    
    for i, ens_id in enumerate(ensaios_alvo):
        ax = axes[i]
        grp = df[df['ensaio_id'] == ens_id]
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        t_exp = grp['tempo_min'].values
        y_exp = grp['X_zn_exp'].values
        
        eh_teste = ens_id in ensaios_teste
        status_tag = "Teste Cego (15%)" if eh_teste else "Treino (85%)"
        cor_curva = '#d62728' if eh_teste else '#2ca02c'
        
        alpha_hat = alfas_preditos[ens_id]
        res_wb = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        res_coelho = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=5500.0)
        res_serial = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_hat)
        
        # Predição nos pontos
        res_pts = pbm.simular_ensaio(tempos_min=t_exp.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_hat)
        y_pts = np.array(res_pts['X_pbm'])
        
        # Linhas de resíduo ponto a ponto
        for tp, yr, yp in zip(t_exp, y_exp, y_pts):
            ax.plot([tp, tp], [yr, yp], color='#d62728' if eh_teste else '#e65100', lw=1.3, alpha=0.8, zorder=3)
            
        ax.plot(t_fino, res_serial['X_pbm'], color=cor_curva, lw=2.5, zorder=5,
                label=f'Híbrido Serial Cinético ({status_tag})')
        ax.plot(t_fino, res_wb['X_pbm'], '--', color='#7f7f7f', lw=1.4, zorder=2,
                label=r'White-Box Fundamental ($\alpha=0$)')
        ax.plot(t_fino, res_coelho['X_pbm'], ':', color='#ff7f0e', lw=1.3, zorder=1,
                label=r'Tese Bortot Coelho ($\alpha=5500\ \mu\mathrm{m/min}$)')
        ax.scatter(t_exp, y_exp, color='black', s=55, zorder=6, edgecolors='white', lw=0.8,
                   label='Experimental ($X_{\\mathrm{Zn}}^{\\mathrm{exp}}$)')
        
        if eta < 1.0:
            ax.axhline(eta, color='gray', linestyle='-.', lw=1.0, alpha=0.7, label=f'Teto Estequiométrico ($\\eta={eta:.1f}$)')
            
        m_ens = calcular_metricas_estatisticas(y_exp, y_pts, p_parametros=2, eta=np.array([eta]*len(y_exp)))
        
        titulo_str = (f"Ensaio {ens_id} [{status_tag}] — $\\eta = {eta:.1f}$, $C_{{A0}} = {ca0:.1f}$ M\n"
                      f"$\\hat{{\\alpha}}_{{\\mathrm{{pred}}}} = {alpha_hat:.0f}\\ \\mu\\mathrm{{m/min}}$ | "
                      f"$R^2 = {m_ens['R2']:.4f}$ | $RMSE = {m_ens['RMSE']*100:.2f}\\%$")
        ax.set_title(titulo_str, fontsize=11.0, fontweight='bold', pad=8)
        ax.set_xlim(-0.3, 15.5)
        ax.set_ylim(-0.02, 1.05)
        ax.set_xlabel('Tempo de Lixiviação (min)', fontsize=10.5)
        ax.set_ylabel('Conversão de Zinco ($X_{\\mathrm{Zn}}$)', fontsize=10.5)
        ax.grid(True, linestyle='--', alpha=0.45)
        ax.legend(loc='lower right', fontsize=9.2, framealpha=0.92)
        
    fig.suptitle('Perfis Cinéticos Detalhados do Modelo Híbrido Serial nos Regimes Críticos de Acidez (85/15)',
                 fontsize=14.0, fontweight='bold', y=0.99)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    caminho_fig = os.path.join(pasta_saida, "curvas_cineticas_ensaios_criticos_serial_85_15.png")
    plt.savefig(caminho_fig)
    plt.close()
    print(f"-> Painel de ensaios críticos do modelo serial salvo em: {caminho_fig}")
    return caminho_fig


# =============================================================================
# 3. PARIDADE E DISPERSÃO DE RESÍDUOS PONTO A PONTO
# =============================================================================

def gerar_paridade_e_residuos_ponto_a_ponto(res_holdout: dict,
                                            pasta_saida: str = "scripts/docs") -> str:
    """Gera painel com gráfico de paridade ponto a ponto e gráfico de resíduo temporal de todos os 128 pontos."""
    os.makedirs(pasta_saida, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), dpi=300)
    
    y_tr_real = res_holdout['y_treino_real']
    y_tr_pred = res_holdout['y_treino_pred']
    y_te_real = res_holdout['y_teste_real']
    y_te_pred = res_holdout['y_teste_pred']
    
    df_tr = res_holdout['df_treino']
    df_te = res_holdout['df_teste']
    t_tr = df_tr['tempo_min'].values
    t_te = df_te['tempo_min'].values
    
    m_tr = res_holdout['metricas_treino']
    m_te = res_holdout['metricas_teste']
    
    # --- Painel 1: Paridade Ponto a Ponto ---
    ax1 = axes[0]
    t_lin = np.linspace(-0.02, 1.05, 200)
    ax1.plot(t_lin, t_lin, 'k-', lw=1.8, label='Paridade Ideal ($y = x$)', zorder=3)
    ax1.fill_between(t_lin, t_lin - 0.05, t_lin + 0.05, color='#2ca02c', alpha=0.12,
                     label=r'Margem $\pm 5\%$', zorder=1)
    ax1.plot(t_lin, t_lin + 0.05, color='#2ca02c', ls='--', lw=0.9, alpha=0.6)
    ax1.plot(t_lin, t_lin - 0.05, color='#2ca02c', ls='--', lw=0.9, alpha=0.6)
    
    ax1.scatter(y_tr_real, y_tr_pred, color='#2ca02c', s=42, alpha=0.65,
                label=f'Treino 85% (13 ensaios, N=104)\n$R^2 = {m_tr["R2"]:.4f}$ | $RMSE = {m_tr["RMSE"]*100:.2f}\\%$',
                zorder=4)
    ax1.scatter(y_te_real, y_te_pred, color='#d62728', s=65, marker='^', alpha=0.92, edgecolors='black', lw=0.9,
                label=f'Teste Cego 15% (Ensaios 1, 2, 6, N=24)\n$R^2 = {m_te["R2"]:.4f}$ | $RMSE = {m_te["RMSE"]*100:.2f}\\%$ | $MAE = {m_te["MAE"]*100:.2f}\\%$',
                zorder=5)
    
    ax1.set_xlim(-0.02, 1.05)
    ax1.set_ylim(-0.02, 1.05)
    ax1.set_xlabel('Conversão Experimental de Zinco ($X_{\\mathrm{Zn}}^{\\mathrm{exp}}$)')
    ax1.set_ylabel('Conversão Prevista pelo Modelo Serial ($X_{\\mathrm{Zn}}^{\\mathrm{prev}}$)')
    ax1.set_title('(A) Gráfico de Paridade Ponto a Ponto (85/15)', fontweight='bold')
    ax1.grid(True, linestyle='--', alpha=0.45)
    ax1.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, framealpha=0.92)
    
    # --- Painel 2: Resíduo Temporal Ponto a Ponto ---
    ax2 = axes[1]
    e_tr = y_tr_real - y_tr_pred
    e_te = y_te_real - y_te_pred
    
    ax2.axhline(0.0, color='black', ls='-', lw=1.5, zorder=2)
    ax2.axhline(0.05, color='#2ca02c', ls='--', lw=1.0, alpha=0.7, label=r'Tolerância $\pm 0{,}05$')
    ax2.axhline(-0.05, color='#2ca02c', ls='--', lw=1.0, alpha=0.7)
    
    ax2.scatter(t_tr - 0.1, e_tr, color='#2ca02c', s=42, alpha=0.65, label='Resíduos Treino (85%)', zorder=3)
    ax2.scatter(t_te + 0.1, e_te, color='#d62728', s=65, marker='^', alpha=0.92, edgecolors='black', lw=0.9,
                label='Resíduos Teste Cego (15%)', zorder=4)
    
    # Linhas de zero a resíduo para teste
    for tp, res_val in zip(t_te, e_te):
        ax2.plot([tp + 0.1, tp + 0.1], [0.0, res_val], color='#d62728', lw=1.0, alpha=0.6, zorder=2)
        
    ax2.set_xlim(-0.5, 15.5)
    ax2.set_ylim(-0.25, 0.25)
    ax2.set_xlabel('Tempo de Lixiviação (min)')
    ax2.set_ylabel('Resíduo Pontual $\\Delta X = X_{\\mathrm{Zn}}^{\\mathrm{exp}} - X_{\\mathrm{Zn}}^{\\mathrm{prev}}$')
    ax2.set_title('(B) Dispersão de Resíduos Ponto a Ponto no Tempo (85/15)', fontweight='bold')
    ax2.grid(True, linestyle='--', alpha=0.45)
    ax2.legend(loc='lower right', frameon=True, fancybox=True, shadow=True, framealpha=0.92)
    
    fig.suptitle('Diagnóstico Ponto a Ponto Completo do Modelo Híbrido Serial Cinético (Divisão 85/15)',
                 fontsize=13.5, fontweight='bold', y=0.99)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    caminho_fig = os.path.join(pasta_saida, "paridade_e_residuos_ponto_a_ponto_serial_85_15.png")
    plt.savefig(caminho_fig)
    plt.close()
    print(f"-> Diagnóstico ponto a ponto salvo em: {caminho_fig}")
    return caminho_fig


# =============================================================================
# EXECUÇÃO PRINCIPAL
# =============================================================================

def main():
    configurar_estilo()
    pasta_saida = "scripts/docs"
    os.makedirs(pasta_saida, exist_ok=True)
    
    print("=" * 100)
    print("GERAÇÃO DE GRÁFICOS PONTO A PONTO DO MODELO HÍBRIDO SERIAL CINÉTICO (SPLIT 85/15)")
    print("=" * 100)
    
    df = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    
    print("-> Inicializando e calibrando Modelo Híbrido Serial...")
    mod_serial = ModeloHibridoAlphaSerial(nome_algoritmo='Gradient Boosting (GBDT)', random_state=42)
    mod_serial.calibrar_alfas_otimos_ensaios(df, verbose=False)
    
    print("-> Executando particionamento estruturado 85% Treino / 15% Teste Cego...")
    res_holdout = mod_serial.validacao_holdout_85_15(df, test_size=0.15, random_state=42)
    
    print(f"   * Ensaios Treino ({res_holdout['pct_treino']:.1f}%): {res_holdout['ensaios_treino']}")
    print(f"   * Ensaios Teste Cego ({res_holdout['pct_teste']:.1f}%): {res_holdout['ensaios_teste']}")
    print(f"   * Métricas Treino: R² = {res_holdout['metricas_treino']['R2']:.4f}, RMSE = {res_holdout['metricas_treino']['RMSE']*100:.2f}%")
    print(f"   * Métricas Teste:  R² = {res_holdout['metricas_teste']['R2']:.4f}, RMSE = {res_holdout['metricas_teste']['RMSE']*100:.2f}%")
    
    print("\n-> Gerando Gráficos Científicos em 300 DPI...")
    f1 = gerar_matriz_16_ensaios_serial(df, res_holdout, pbm, pasta_saida)
    f2 = gerar_painel_ensaios_criticos_serial(df, res_holdout, pbm, pasta_saida)
    f3 = gerar_paridade_e_residuos_ponto_a_ponto(res_holdout, pasta_saida)
    
    print("\n" + "=" * 100)
    print("GRÁFICOS PONTO A PONTO DO HÍBRIDO SERIAL GERADOS COM SUCESSO!")
    print(f"  1. Matriz 16 ensaios: {f1}")
    print(f"  2. Ensaios críticos: {f2}")
    print(f"  3. Paridade e resíduos: {f3}")
    print("=" * 100)


if __name__ == "__main__":
    main()
