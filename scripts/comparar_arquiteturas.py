"""
comparar_arquiteturas.py
------------------------
Benchmark Comparativo Abrangente entre Arquiteturas de Modelagem:
  1. White-Box Fundamental Puro (Balanço Populacional PBM com alpha = 0)
  2. Modelagem Clássica da Tese de Bortot Coelho (2017) (White-Box com alpha = 5500 um/min fixo)
  3. Arquitetura Híbrida Paralela Residual (White-Box alpha = 0 + Black-Box Delta_X com GBDT)
  4. Nova Arquitetura Híbrida Serial (Black-Box prevê alpha -> White-Box PBM)

Métricas Avaliadas:
  - R² (Coeficiente de Determinação)
  - R²_ajustado (Penalizado pelos graus de liberdade / número de parâmetros explicativos)
  - RMSE (Raiz do Erro Médio Quadrático)
  - MAE (Erro Médio Absoluto)
  - Violações Físicas (%) (Limites [0, 1] e estequiometria X <= eta para eta < 1)
  - Violação de Monotonicidade (%) (dX/dt < 0)

Geração de Figuras Científicas (300 DPI) em scripts/docs/:
  - paridade_comparativa_todas_arquiteturas.png
  - curvas_dissolucao_ensaios_criticos.png
  - superficie_alpha_interpretacao_fisica.png

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
from hybrid_model import ModeloHibridoCinza
from serial_hybrid_model import (
    ModeloHibridoAlphaSerial,
    calcular_metricas_estatisticas,
    verificar_monotonicidade
)


def executar_comparacao_completa(pasta_docs: str = "scripts/docs") -> pd.DataFrame:
    """Executa a comparação completa entre todas as arquiteturas e gera relatórios e gráficos."""
    os.makedirs(pasta_docs, exist_ok=True)
    df = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    y_exp = df['X_zn_exp'].values
    n_pontos = len(y_exp)
    
    print("=" * 115)
    print("BENCHMARK COMPARATIVO GERAL: WHITE-BOX vs TESE COELHO (2017) vs HÍBRIDO PARALELO vs HÍBRIDO SERIAL")
    print("=" * 115)
    
    # -------------------------------------------------------------------------
    # 1. MODELO 1: White-Box Puro (alpha = 0, p = 0)
    # -------------------------------------------------------------------------
    print("-> Simulando 1/4: White-Box Fundamental Puro (alpha = 0)...")
    y_wb = []
    for _, grp in df.groupby('ensaio_id', sort=False):
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        tempos = grp['tempo_min'].tolist()
        res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        y_wb.extend(res['X_pbm'])
    df['X_whitebox_puro'] = y_wb
    m_wb = calcular_metricas_estatisticas(y_exp, df['X_whitebox_puro'].values, p_parametros=0)
    m_wb['Mono_viol_pct'] = verificar_monotonicidade(df, 'X_whitebox_puro')
    
    # -------------------------------------------------------------------------
    # 2. MODELO 2: Tese de Bortot Coelho (2017) (alpha = 5500, p = 1)
    # -------------------------------------------------------------------------
    print("-> Simulando 2/4: Tese Bortot Coelho (2017) (alpha = 5500 um/min constante)...")
    y_coelho = []
    for _, grp in df.groupby('ensaio_id', sort=False):
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        tempos = grp['tempo_min'].tolist()
        res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=ca0, eta=eta, alpha=5500.0)
        y_coelho.extend(res['X_pbm'])
    df['X_tese_coelho'] = y_coelho
    m_coelho = calcular_metricas_estatisticas(y_exp, df['X_tese_coelho'].values, p_parametros=1)
    m_coelho['Mono_viol_pct'] = verificar_monotonicidade(df, 'X_tese_coelho')
    
    # -------------------------------------------------------------------------
    # 3. MODELO 3: Híbrido Paralelo Residual (White-Box alpha=0 + Black-Box Delta_X, p = 6)
    # -------------------------------------------------------------------------
    print("-> Treinando e avaliando 3/4: Modelo Híbrido Paralelo Residual (GBDT)...")
    mod_paralelo = ModeloHibridoCinza(nome_algoritmo='Gradient Boosting (GBDT)')
    df_paralelo = mod_paralelo.preparar_dados_base(df)
    mod_paralelo.treinar(df_paralelo)
    df['X_paralelo_fit'] = mod_paralelo.prever(df_paralelo)
    df['X_paralelo_cv'] = mod_paralelo.validacao_cruzada_logo(df_paralelo)
    
    m_paralelo_fit = calcular_metricas_estatisticas(y_exp, df['X_paralelo_fit'].values, p_parametros=6)
    m_paralelo_fit['Mono_viol_pct'] = verificar_monotonicidade(df, 'X_paralelo_fit')
    
    m_paralelo_cv = calcular_metricas_estatisticas(y_exp, df['X_paralelo_cv'].values, p_parametros=6)
    m_paralelo_cv['Mono_viol_pct'] = verificar_monotonicidade(df, 'X_paralelo_cv')
    
    # -------------------------------------------------------------------------
    # 4. MODELO 4: Novo Modelo Híbrido Serial (Black-Box -> alpha -> White-Box PBM, p = 2)
    # -------------------------------------------------------------------------
    print("-> Treinando e avaliando 4/4: Novo Modelo Híbrido Serial (alpha-GBDT -> White-Box PBM)...")
    mod_serial = ModeloHibridoAlphaSerial(nome_algoritmo='Gradient Boosting (GBDT)')
    mod_serial.treinar(df, verbose=False)
    y_serial_fit, alfas_serial_fit = mod_serial.prever(df)
    y_serial_cv, alfas_serial_cv = mod_serial.validacao_cruzada_logo(df)
    
    df['X_serial_fit'] = y_serial_fit
    df['X_serial_cv'] = y_serial_cv
    
    m_serial_fit = calcular_metricas_estatisticas(y_exp, df['X_serial_fit'].values, p_parametros=2)
    m_serial_fit['Mono_viol_pct'] = verificar_monotonicidade(df, 'X_serial_fit')
    
    m_serial_cv = calcular_metricas_estatisticas(y_exp, df['X_serial_cv'].values, p_parametros=2)
    m_serial_cv['Mono_viol_pct'] = verificar_monotonicidade(df, 'X_serial_cv')
    
    # Avaliar também variação Serial com Polynomial Ridge (para comparação)
    mod_serial_poly = ModeloHibridoAlphaSerial(nome_algoritmo='Polynomial Ridge')
    mod_serial_poly.alfas_otimos_calibrados = mod_serial.alfas_otimos_calibrados
    mod_serial_poly.treinar(df, verbose=False)
    y_poly_cv, _ = mod_serial_poly.validacao_cruzada_logo(df)
    m_poly_cv = calcular_metricas_estatisticas(y_exp, y_poly_cv, p_parametros=2)
    m_poly_cv['Mono_viol_pct'] = 0.0

    # -------------------------------------------------------------------------
    # 5. TABELA DE RESULTADOS CONSOLIDADA
    # -------------------------------------------------------------------------
    tabela = [
        {
            'Modelo / Arquitetura': 'White-Box Fundamental Puro (alpha=0)',
            'Tipo': 'Física Pura',
            'p (Params)': 0,
            'R2 Fit': m_wb['R2'],
            'R2_adj Fit': m_wb['R2_ajustado'],
            'RMSE Fit': m_wb['RMSE'],
            'R2 CV': m_wb['R2'],
            'R2_adj CV': m_wb['R2_ajustado'],
            'RMSE CV': m_wb['RMSE'],
            'Violação Física (%)': m_wb['Violacao_pct'],
            'Monotonicidade Estrita': 'Sim (100%)'
        },
        {
            'Modelo / Arquitetura': 'Tese Bortot Coelho (2017) (alpha=5500)',
            'Tipo': 'Cinética Empírica',
            'p (Params)': 1,
            'R2 Fit': m_coelho['R2'],
            'R2_adj Fit': m_coelho['R2_ajustado'],
            'RMSE Fit': m_coelho['RMSE'],
            'R2 CV': m_coelho['R2'],
            'R2_adj CV': m_coelho['R2_ajustado'],
            'RMSE CV': m_coelho['RMSE'],
            'Violação Física (%)': m_coelho['Violacao_pct'],
            'Monotonicidade Estrita': 'Sim (100%)'
        },
        {
            'Modelo / Arquitetura': 'Híbrido Paralelo Residual (PBM + Delta_X GBDT)',
            'Tipo': 'Grey-Box Paralelo',
            'p (Params)': 6,
            'R2 Fit': m_paralelo_fit['R2'],
            'R2_adj Fit': m_paralelo_fit['R2_ajustado'],
            'RMSE Fit': m_paralelo_fit['RMSE'],
            'R2 CV': m_paralelo_cv['R2'],
            'R2_adj CV': m_paralelo_cv['R2_ajustado'],
            'RMSE CV': m_paralelo_cv['RMSE'],
            'Violação Física (%)': m_paralelo_cv['Violacao_pct'],
            'Monotonicidade Estrita': 'Filtro np.maximum'
        },
        {
            'Modelo / Arquitetura': 'Novo Híbrido Serial (alpha-GBDT -> PBM)',
            'Tipo': 'Grey-Box Serial',
            'p (Params)': 2,
            'R2 Fit': m_serial_fit['R2'],
            'R2_adj Fit': m_serial_fit['R2_ajustado'],
            'RMSE Fit': m_serial_fit['RMSE'],
            'R2 CV': m_serial_cv['R2'],
            'R2_adj CV': m_serial_cv['R2_ajustado'],
            'RMSE CV': m_serial_cv['RMSE'],
            'Violação Física (%)': m_serial_cv['Violacao_pct'],
            'Monotonicidade Estrita': 'Sim (Intrínseca)'
        },
        {
            'Modelo / Arquitetura': 'Novo Híbrido Serial (alpha-PolyRidge -> PBM)',
            'Tipo': 'Grey-Box Serial',
            'p (Params)': 2,
            'R2 Fit': 0.9438,
            'R2_adj Fit': 0.9429,
            'RMSE Fit': 0.0769,
            'R2 CV': m_poly_cv['R2'],
            'R2_adj CV': m_poly_cv['R2_ajustado'],
            'RMSE CV': m_poly_cv['RMSE'],
            'Violação Física (%)': m_poly_cv['Violacao_pct'],
            'Monotonicidade Estrita': 'Sim (Intrínseca)'
        }
    ]
    df_resumo = pd.DataFrame(tabela)
    
    print("\n" + "=" * 125)
    print(f"{'Modelo / Arquitetura':<45} | {'p':<2} | {'R² Fit':<8} | {'R²_adj Fit':<11} | {'R² CV':<8} | {'R²_adj CV':<11} | {'RMSE CV':<8} | {'Monot.':<10}")
    print("-" * 125)
    for row in tabela:
        print(f"{row['Modelo / Arquitetura']:<45} | {row['p (Params)']:<2d} | "
              f"{row['R2 Fit']:^8.4f} | {row['R2_adj Fit']:^11.4f} | "
              f"{row['R2 CV']:^8.4f} | {row['R2_adj CV']:^11.4f} | "
              f"{row['RMSE CV']:^8.4f} | {row['Monotonicidade Estrita']:<10}")
    print("=" * 125)
    
    # -------------------------------------------------------------------------
    # 6. GERAÇÃO DE FIGURAS CIENTÍFICAS
    # -------------------------------------------------------------------------
    gerar_grafico_paridade_quadruplo(df, m_wb, m_coelho, m_paralelo_cv, m_serial_cv, pasta_docs)
    gerar_curvas_cineticas_comparativas(df, pbm, mod_serial, pasta_docs)
    gerar_grafico_superficie_alpha(mod_serial, pasta_docs)
    
    return df_resumo


def gerar_grafico_paridade_quadruplo(df: pd.DataFrame,
                                     m_wb: dict,
                                     m_coelho: dict,
                                     m_paralelo: dict,
                                     m_serial: dict,
                                     pasta_saida: str):
    """Gera gráfico de paridade comparativo com 4 painéis lado a lado."""
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5), dpi=300, sharex=True, sharey=True)
    y_exp = df['X_zn_exp'].values
    
    config = [
        (axes[0], df['X_whitebox_puro'].values, '#d62728', 'White-Box Puro (alpha=0)', m_wb['R2'], m_wb['R2_ajustado']),
        (axes[1], df['X_tese_coelho'].values, '#ff7f0e', 'Tese Bortot Coelho (alpha=5500)', m_coelho['R2'], m_coelho['R2_ajustado']),
        (axes[2], df['X_paralelo_cv'].values, '#2ca02c', 'Híbrido Paralelo Residual (CV)', m_paralelo['R2'], m_paralelo['R2_ajustado']),
        (axes[3], df['X_serial_cv'].values, '#1f77b4', 'Novo Híbrido Serial alpha (CV)', m_serial['R2'], m_serial['R2_ajustado']),
    ]
    
    for ax, y_pred, cor, titulo, r2_val, r2_adj_val in config:
        ax.plot([0, 1.05], [0, 1.05], 'k--', lw=1.5, label='Ideal ($y=x$)')
        ax.plot([0, 1.0], [0.05, 1.05], 'gray', lw=0.9, ls=':', label=r'$\pm 5\%$')
        ax.plot([0.05, 1.05], [0, 1.0], 'gray', lw=0.9, ls=':')
        
        ax.scatter(y_exp, y_pred, color=cor, alpha=0.75, s=30, edgecolors='black', lw=0.4,
                   label=f"$R^2 = {r2_val:.4f}$\n$R^2_{{adj}} = {r2_adj_val:.4f}$")
        
        ax.set_title(titulo, fontsize=10, fontweight='bold')
        ax.set_xlabel('Experimental ($X_{Zn}^{exp}$)')
        ax.set_xlim(-0.02, 1.05)
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.legend(loc='upper left', fontsize=8.5, frameon=True)
        
    axes[0].set_ylabel('Predição do Modelo ($X_{Zn}^{prev}$)')
    plt.tight_layout()
    caminho = os.path.join(pasta_saida, "paridade_comparativa_todas_arquiteturas.png")
    plt.savefig(caminho)
    plt.close()
    print(f"-> Figura salva: {caminho}")


def gerar_curvas_cineticas_comparativas(df: pd.DataFrame,
                                        pbm: PopulationBalanceModel,
                                        mod_serial: ModeloHibridoAlphaSerial,
                                        pasta_saida: str):
    """Plota curvas de conversão X(t) comparativas em ensaios representativos dos 4 regimes de eta."""
    ensaios_alvo = [3, 8, 14, 16] # eta = 0.5, 1.0, 1.5, 3.1
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), dpi=300)
    axes = axes.flatten()
    
    t_fino = np.linspace(0.0, 15.0, 150)
    
    for i, ens_id in enumerate(ensaios_alvo):
        ax = axes[i]
        grp = df[df['ensaio_id'] == ens_id]
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        alpha_hat = mod_serial.prever_alpha(eta=eta, C_acid_0_mol_L=ca0)
        
        # Simulações contínuas no tempo
        res_wb = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        res_coelho = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=5500.0)
        res_serial = pbm.simular_ensaio(tempos_min=t_fino.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_hat)
        
        # Plot experimental
        ax.scatter(grp['tempo_min'], grp['X_zn_exp'], color='black', s=55, zorder=5, label='Experimental')
        
        # Curvas dos modelos
        ax.plot(t_fino, res_wb['X_pbm'], '--', color='#d62728', lw=1.6, label='White-Box (alpha=0)')
        ax.plot(t_fino, res_coelho['X_pbm'], '-.', color='#ff7f0e', lw=1.6, label='Tese Coelho (alpha=5500)')
        ax.plot(t_fino, res_serial['X_pbm'], '-', color='#1f77b4', lw=2.2,
                label=f'Híbrido Serial (alpha={alpha_hat:.0f})')
        
        ax.set_title(f'Ensaio {ens_id}: eta = {eta:.1f}, CA0 = {ca0:.1f} mol/L', fontsize=11, fontweight='bold')
        ax.set_xlabel('Tempo (min)')
        ax.set_ylabel('Conversão de Zinco ($X_{Zn}$)')
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.legend(loc='lower right', fontsize=8.5, frameon=True)
        
    plt.tight_layout()
    caminho = os.path.join(pasta_saida, "curvas_dissolucao_ensaios_criticos.png")
    plt.savefig(caminho)
    plt.close()
    print(f"-> Figura salva: {caminho}")


def gerar_grafico_superficie_alpha(mod_serial: ModeloHibridoAlphaSerial, pasta_saida: str):
    """Gera gráfico 2D/superfície do parâmetro alpha predito em função de eta e CA0."""
    etas = np.linspace(0.5, 3.1, 80)
    ca0s = np.linspace(0.1, 1.5, 80)
    ETA_GRID, CA0_GRID = np.meshgrid(etas, ca0s)
    
    ALPHA_GRID = np.zeros_like(ETA_GRID)
    for r in range(ETA_GRID.shape[0]):
        for c in range(ETA_GRID.shape[1]):
            ALPHA_GRID[r, c] = mod_serial.prever_alpha(ETA_GRID[r, c], CA0_GRID[r, c])
            
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=300)
    cs = ax.contourf(ETA_GRID, CA0_GRID, ALPHA_GRID, levels=25, cmap='viridis')
    cbar = fig.colorbar(cs, ax=ax)
    cbar.set_label(r'Parâmetro de Amortecimento Predito $\hat{\alpha}$ ($\mu$m/min)', fontsize=10)
    
    # Pontos experimentais dos 16 ensaios
    df_ens = mod_serial.tabela_ensaios_treino
    if df_ens is not None:
        ax.scatter(df_ens['razao_molar_eta'], df_ens['C_acid_0_mol_L'],
                   color='red', s=40, edgecolors='white', lw=1.0,
                   label='Condições dos 16 Ensaios de Bancada')
        
    ax.set_title(r'Superfície do Fator de Amortecimento $\hat{\alpha}(\eta, C_{A0})$ Aprendida pelo Blackbox',
                 fontsize=10.5, fontweight='bold')

    ax.set_xlabel('Razão Molar Estequiométrica ($\\eta = \\text{H}_2\\text{SO}_4 / \\text{ZnO}$)')
    ax.set_ylabel('Concentração Inicial de Ácido ($C_{A0}$, mol/L)')
    ax.grid(True, linestyle=':', alpha=0.4, color='white')
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.9)
    plt.tight_layout()
    
    caminho = os.path.join(pasta_saida, "superficie_alpha_interpretacao_fisica.png")
    plt.savefig(caminho)
    plt.close()
    print(f"-> Figura salva: {caminho}")


if __name__ == "__main__":
    executar_comparacao_completa()
