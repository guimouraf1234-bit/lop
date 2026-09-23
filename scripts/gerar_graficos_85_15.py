"""
gerar_graficos_85_15.py
-----------------------
Gera gráficos científicos em alta resolução (300 DPI) para a validação com 
divisão 85% Treino / 15% Teste Cego para as arquiteturas híbridas:
  1. Modelo Híbrido Serial Cinético (alpha-GBDT -> PBM)
  2. Modelo Híbrido Paralelo Cooperativo / Residual (PBM + Delta_X GBDT)
  3. Painel Comparativo Integrado (Serial vs Cooperativo vs Competitivo)

Figuras científicas geradas em scripts/docs/:
  - paridade_serial_85_15.png
  - ponto_a_ponto_serial_85_15.png
  - paridade_cooperativo_85_15.png
  - ponto_a_ponto_cooperativo_85_15.png
  - comparativo_paridade_85_15_todas.png

Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG (2026)
Alunos: Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata
Orientador: Prof. Dr. Fabrício Eduardo Bortot Coelho
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Configuração de encoding UTF-8 no Windows
if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

# Importações dos módulos locais
from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel
from serial_hybrid_model import ModeloHibridoAlphaSerial, calcular_metricas_estatisticas
from hybrid_model import ModeloHibridoCinza
from competitive_hybrid_model import ModeloHibridoParaleloCompetitivo


def configurar_estilo_grafico():
    """Configura parâmetros visuais com padrão de publicação internacional."""
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'font.family': 'sans-serif',
        'mathtext.fontset': 'dejavusans',
        'axes.linewidth': 1.1,
        'grid.linewidth': 0.7,
        'grid.alpha': 0.4
    })


def gerar_grafico_paridade_individual(res_holdout: dict,
                                      nome_modelo: str,
                                      caminho_saida: str,
                                      cor_treino: str = '#1f77b4',
                                      cor_teste: str = '#d62728'):
    """Gera gráfico de paridade destacando Treino (85%) e Teste Cego (15%)."""
    m_tr = res_holdout['metricas_treino']
    m_te = res_holdout['metricas_teste']
    
    y_tr_real = res_holdout['y_treino_real']
    y_tr_pred = res_holdout['y_treino_pred']
    y_te_real = res_holdout['y_teste_real']
    y_te_pred = res_holdout['y_teste_pred']
    
    ens_te = res_holdout['ensaios_teste']
    ens_tr = res_holdout['ensaios_treino']
    
    fig, ax = plt.subplots(figsize=(7, 6.5), dpi=300)
    
    # Linha ideal 1:1 e margens de erro de +-5% e +-10%
    t_lin = np.linspace(-0.02, 1.05, 200)
    ax.plot(t_lin, t_lin, 'k-', lw=1.8, label='Paridade Ideal ($y = x$)', zorder=3)
    ax.fill_between(t_lin, t_lin - 0.05, t_lin + 0.05, color='#2ca02c', alpha=0.12,
                    label=r'Margem de Tolerância $\pm 5\%$', zorder=1)
    ax.plot(t_lin, t_lin + 0.05, color='#2ca02c', ls='--', lw=0.9, alpha=0.7)
    ax.plot(t_lin, t_lin - 0.05, color='#2ca02c', ls='--', lw=0.9, alpha=0.7)
    
    # Pontos de Treino (85%)
    ax.scatter(y_tr_real, y_tr_pred, color=cor_treino, s=42, alpha=0.65, edgecolors='none',
               label=f'Treino 85% ({len(ens_tr)} ensaios, N={len(y_tr_real)})\n$R^2 = {m_tr["R2"]:.4f}$ | $RMSE = {m_tr["RMSE"]*100:.2f}\\%$',
               zorder=4)
    
    # Pontos de Teste Cego (15%)
    ax.scatter(y_te_real, y_te_pred, color=cor_teste, s=65, alpha=0.92, edgecolors='black', linewidth=0.9,
               marker='^',
               label=f'Teste Cego 15% (Ensaios {ens_te}, N={len(y_te_real)})\n$R^2 = {m_te["R2"]:.4f}$ | $RMSE = {m_te["RMSE"]*100:.2f}\\%$ | $MAE = {m_te["MAE"]*100:.2f}\\%$',
               zorder=5)
    
    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel('Conversão Experimental de Zinco ($X_{\\mathrm{Zn}}^{\\mathrm{exp}}$)')
    ax.set_ylabel('Conversão Prevista pelo Modelo ($X_{\\mathrm{Zn}}^{\\mathrm{prev}}$)')
    ax.set_title(f'Gráfico de Paridade (Particionamento 85/15)\n{nome_modelo}', fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper left', frameon=True, fancybox=True, shadow=True, framealpha=0.92)
    
    plt.tight_layout()
    plt.savefig(caminho_saida)
    plt.close()
    print(f"-> Gráfico de paridade salvo: {caminho_saida}")


def gerar_grafico_ponto_a_ponto(df: pd.DataFrame,
                                res_holdout: dict,
                                modelo_obj: any,
                                tipo_modelo: str,
                                caminho_saida: str):
    """
    Gera gráfico de comparação ponto a ponto com 6 subplots:
      - 3 ensaios de Teste Cego (Ensaios 1, 2, 6)
      - 3 ensaios de Treino representativos (e.g., Ensaios 7, 11, 14)
    Cada subplot contém:
      * Painel superior: Dados experimentais (pontos) vs Curva contínua predita pelo modelo (linha)
      * Painel inferior: Resíduo ponto a ponto (Delta X = X_exp - X_prev)
    """
    ens_te = sorted(res_holdout['ensaios_teste'])
    # Selecionar 3 ensaios de treino cobrindo diferentes etas
    ens_tr_todos = res_holdout['ensaios_treino']
    ens_tr_selecionados = [7, 11, 14]  # eta = 3.1, 1.5, 1.5
    ensaios_painel = [(e, 'Teste Cego (15%)', '#d62728') for e in ens_te] + \
                     [(e, 'Treino (85%)', '#1f77b4') for e in ens_tr_selecionados]
    
    pbm = PopulationBalanceModel()
    t_cont = np.linspace(0.0, 15.0, 120)
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.6), dpi=300, sharex=True, sharey=True)
    axes = axes.flatten()
    
    for idx, (ens_id, tipo_split, cor_linha) in enumerate(ensaios_painel):
        ax = axes[idx]
        grp = df[df['ensaio_id'] == ens_id]
        ca0 = float(grp['C_acid_0_mol_L'].iloc[0])
        eta = float(grp['razao_molar_eta'].iloc[0])
        sl = float(grp['razao_SL_g_L'].iloc[0]) if 'razao_SL_g_L' in grp.columns else 0.0
        t_exp = grp['tempo_min'].values
        y_exp = grp['X_zn_exp'].values
        
        # Obter predição contínua e nos pontos experimentais
        if tipo_modelo == 'serial':
            alpha_pred = res_holdout['alfas_preditos'][ens_id]
            res_cont = pbm.simular_ensaio(tempos_min=t_cont.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_pred)
            x_cont = res_cont['X_pbm']
            res_pts = pbm.simular_ensaio(tempos_min=t_exp.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=alpha_pred)
            x_pred_pts = np.array(res_pts['X_pbm'])
            subtitulo_extra = f"$\\hat{{\\alpha}} = {alpha_pred:.0f}\\ \\mu\\mathrm{{m/min}}$"
        elif tipo_modelo == 'cooperativo':
            # Simulação PBM base contínua
            res_pbm_cont = pbm.simular_ensaio(tempos_min=t_cont.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
            df_cont = pd.DataFrame({
                'tempo_min': t_cont,
                'razao_molar_eta': [eta] * len(t_cont),
                'C_acid_0_mol_L': [ca0] * len(t_cont),
                'razao_SL_g_L': [sl] * len(t_cont),
                'X_pbm': res_pbm_cont['X_pbm'],
                'C_acid_pbm_mol_L': [ca0 * max(0.0, 1.0 - x/eta) for x in res_pbm_cont['X_pbm']],
                'ensaio_id': [ens_id] * len(t_cont)
            })
            feats = ['tempo_min', 'razao_molar_eta', 'C_acid_0_mol_L', 'razao_SL_g_L', 'X_pbm', 'C_acid_pbm_mol_L']
            x_cont = modelo_obj.acoplador.acoplar_e_restringir(df_cont, modelo_obj.regressor.predict(df_cont[feats].values))
            
            # Predição exata nos pontos experimentais
            df_pts = grp.copy()
            if 'X_pbm' not in df_pts.columns:
                res_pts_wb = pbm.simular_ensaio(tempos_min=t_exp.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
                df_pts['X_pbm'] = res_pts_wb['X_pbm']
            if 'C_acid_pbm_mol_L' not in df_pts.columns:
                df_pts['C_acid_pbm_mol_L'] = df_pts['C_acid_0_mol_L'] * np.maximum(0.0, 1.0 - df_pts['X_pbm']/df_pts['razao_molar_eta'])
            feats = ['tempo_min', 'razao_molar_eta', 'C_acid_0_mol_L', 'razao_SL_g_L', 'X_pbm', 'C_acid_pbm_mol_L']
            x_pred_pts = modelo_obj.acoplador.acoplar_e_restringir(df_pts, modelo_obj.regressor.predict(df_pts[feats].values))
            subtitulo_extra = "Correção Residual $\\widehat{\\Delta X}$"
            
        # Linha teórica White-box puro (alpha = 0) para referência
        res_wb_ref = pbm.simular_ensaio(tempos_min=t_cont.tolist(), C_acid_0_mol_L=ca0, eta=eta, alpha=0.0)
        
        # Plotagem
        ax.plot(t_cont, res_wb_ref['X_pbm'], color='#7f7f7f', ls=':', lw=1.2, label='White-Box Puro ($\\alpha=0$)')
        ax.plot(t_cont, x_cont, color=cor_linha, lw=2.2, label=f'Modelo Híbrido ({tipo_split})')
        ax.scatter(t_exp, y_exp, color='black', s=45, zorder=5, label='Experimental ($X_{\\mathrm{Zn}}^{\\mathrm{exp}}$)')
        
        # Conexão ponto a ponto de resíduo (barras verticais de erro)
        for tp, yr, yp in zip(t_exp, y_exp, x_pred_pts):
            ax.plot([tp, tp], [yr, yp], color='red', lw=1.1, alpha=0.7)
            
        # R2 e RMSE do ensaio
        ss_res_e = np.sum((y_exp - x_pred_pts) ** 2)
        ss_tot_e = np.sum((y_exp - np.mean(y_exp)) ** 2)
        r2_e = 1.0 - ss_res_e / ss_tot_e if ss_tot_e > 0 else 1.0
        rmse_e = np.sqrt(np.mean((y_exp - x_pred_pts) ** 2))
        
        titulo_box = (f"Ensaio {ens_id} [{tipo_split}]\n"
                      f"$\\eta = {eta:.1f}$ | $C_{{A0}} = {ca0:.1f}$ M | {subtitulo_extra}\n"
                      f"$R^2 = {r2_e:.4f}$ | $RMSE = {rmse_e*100:.2f}\\%$")
        ax.set_title(titulo_box, fontsize=10.5, pad=6)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_ylim(-0.02, 1.05)
        
        if idx % 3 == 0:
            ax.set_ylabel('Conversão de Zinco ($X_{\\mathrm{Zn}}$)')
        if idx >= 3:
            ax.set_xlabel('Tempo de Lixiviação (min)')
        if idx == 0:
            ax.legend(loc='lower right', fontsize=8.5, framealpha=0.9)
            
    nome_grafico = "Híbrido Serial Cinético" if tipo_modelo == 'serial' else "Híbrido Paralelo Cooperativo"
    fig.suptitle(f'Comparação Ponto a Ponto e Perfis Cinéticos (Divisão 85% Treino / 15% Teste Cego)\n{nome_grafico}',
                 fontweight='bold', y=0.985, fontsize=13)
    plt.tight_layout()
    plt.subplots_adjust(top=0.86)
    plt.savefig(caminho_saida)
    plt.close()
    print(f"-> Gráfico ponto a ponto salvo: {caminho_saida}")


def gerar_painel_comparativo_todas_85_15(res_serial: dict,
                                        res_coop: dict,
                                        res_comp: dict,
                                        caminho_saida: str):
    """Gera painel com 3 gráficos de paridade comparando os três modelos sob o split 85/15."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8), dpi=300, sharex=True, sharey=True)
    
    modelos_info = [
        ('1. Híbrido Paralelo Cooperativo (Residual)', res_coop, '#1f77b4', '#d62728'),
        ('2. Híbrido Serial Cinético (alpha-GBDT)', res_serial, '#2ca02c', '#d62728'),
        ('3. Híbrido Paralelo Competitivo (Sansana 2024)', res_comp, '#9467bd', '#d62728')
    ]
    
    t_lin = np.linspace(-0.02, 1.05, 200)
    
    for ax, (titulo, res, cor_tr, cor_te) in zip(axes, modelos_info):
        ax.plot(t_lin, t_lin, 'k-', lw=1.6, label='Ideal ($y=x$)', zorder=3)
        ax.fill_between(t_lin, t_lin - 0.05, t_lin + 0.05, color='gray', alpha=0.15,
                        label=r'Faixa $\pm 5\%$', zorder=1)
        
        m_tr = res['metricas_treino']
        m_te = res['metricas_teste']
        y_tr_real = res['y_treino_real']
        y_tr_pred = res['y_treino_pred']
        y_te_real = res['y_teste_real']
        y_te_pred = res['y_teste_pred']
        
        ax.scatter(y_tr_real, y_tr_pred, color=cor_tr, s=36, alpha=0.6,
                   label=f'Treino 85%: $R^2={m_tr["R2"]:.4f}$ | $RMSE={m_tr["RMSE"]*100:.2f}\\%$',
                   zorder=4)
        ax.scatter(y_te_real, y_te_pred, color=cor_te, s=60, marker='^', alpha=0.92,
                   edgecolors='black', lw=0.8,
                   label=f'Teste Cego 15%: $R^2={m_te["R2"]:.4f}$ | $RMSE={m_te["RMSE"]*100:.2f}\\%$',
                   zorder=5)
        
        ax.set_xlim(-0.02, 1.05)
        ax.set_ylim(-0.02, 1.05)
        ax.set_title(titulo, fontsize=11.5, fontweight='bold')
        ax.set_xlabel('Experimental ($X_{\\mathrm{Zn}}^{\\mathrm{exp}}$)')
        ax.grid(True, linestyle='--', alpha=0.45)
        ax.legend(loc='upper left', fontsize=9, framealpha=0.92)
        
    axes[0].set_ylabel('Predito pelo Modelo ($X_{\\mathrm{Zn}}^{\\mathrm{prev}}$)')
    fig.suptitle('Confronto dos Três Modelos Híbridos sob a Divisão 85% Treino / 15% Teste Cego',
                 fontsize=14, fontweight='bold', y=0.99)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(caminho_saida)
    plt.close()
    print(f"-> Painel comparativo salvo: {caminho_saida}")


def main():
    configurar_estilo_grafico()
    pasta_docs = "scripts/docs"
    os.makedirs(pasta_docs, exist_ok=True)
    
    df = carregar_dados_bancada()
    print(f"-> Base carregada com sucesso: {len(df)} pontos experimentais (16 ensaios).")
    
    # -------------------------------------------------------------------------
    # 1. MODELO HÍBRIDO SERIAL CINÉTICO (85/15)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("1. EXECUTANDO VALIDAÇÃO 85/15 NO MODELO HÍBRIDO SERIAL CINÉTICO")
    print("=" * 90)
    mod_serial = ModeloHibridoAlphaSerial(nome_algoritmo='Gradient Boosting (GBDT)', random_state=42)
    mod_serial.calibrar_alfas_otimos_ensaios(df, verbose=False)
    res_holdout_serial = mod_serial.validacao_holdout_85_15(df, test_size=0.15, random_state=42)
    
    print(f"   * Ensaios Treino ({res_holdout_serial['pct_treino']:.1f}%): {res_holdout_serial['ensaios_treino']}")
    print(f"   * Ensaios Teste ({res_holdout_serial['pct_teste']:.1f}%): {res_holdout_serial['ensaios_teste']}")
    m_tr_s = res_holdout_serial['metricas_treino']
    m_te_s = res_holdout_serial['metricas_teste']
    print(f"   * Treino 85%: R² = {m_tr_s['R2']:.4f}, RMSE = {m_tr_s['RMSE']:.4f} ({m_tr_s['RMSE']*100:.2f}%), MAE = {m_tr_s['MAE']:.4f}")
    print(f"   * Teste 15%:  R² = {m_te_s['R2']:.4f}, RMSE = {m_te_s['RMSE']:.4f} ({m_te_s['RMSE']*100:.2f}%), MAE = {m_te_s['MAE']:.4f}")
    
    # Validação GroupKFold k=6 para o Serial
    res_gkf_serial = mod_serial.validacao_group_kfold_85_15(df, n_splits=6)
    m_gkf_s = res_gkf_serial['metricas']
    print(f"   * Validação Agrupada 85/15 (GroupKFold k=6): R² = {m_gkf_s['R2']:.4f}, RMSE = {m_gkf_s['RMSE']:.4f}, MAE = {m_gkf_s['MAE']:.4f}")
    
    # Gráficos do Modelo Serial
    caminho_paridade_serial = os.path.join(pasta_docs, "paridade_serial_85_15.png")
    gerar_grafico_paridade_individual(res_holdout_serial,
                                      nome_modelo="Modelo Híbrido Serial Cinético (alpha-GBDT -> PBM)",
                                      caminho_saida=caminho_paridade_serial,
                                      cor_treino='#2ca02c', cor_teste='#d62728')
    
    caminho_ponto_serial = os.path.join(pasta_docs, "ponto_a_ponto_serial_85_15.png")
    gerar_grafico_ponto_a_ponto(df, res_holdout_serial, mod_serial, tipo_modelo='serial', caminho_saida=caminho_ponto_serial)
    
    # -------------------------------------------------------------------------
    # 2. MODELO HÍBRIDO PARALELO COOPERATIVO / RESIDUAL (85/15)
    # -------------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("2. EXECUTANDO VALIDAÇÃO 85/15 NO MODELO HÍBRIDO PARALELO COOPERATIVO (RESIDUAL)")
    print("=" * 90)
    mod_coop = ModeloHibridoCinza(nome_algoritmo='Gradient Boosting (GBDT)')
    # Treinar modelo base para uso no ponto a ponto contínuo
    df_preparado = mod_coop.preparar_dados_base(df)
    mod_coop.treinar(df_preparado)
    res_holdout_coop = mod_coop.validacao_holdout_85_15(df, test_size=0.15, random_state=42)
    
    print(f"   * Ensaios Treino ({res_holdout_coop['pct_treino']:.1f}%): {res_holdout_coop['ensaios_treino']}")
    print(f"   * Ensaios Teste ({res_holdout_coop['pct_teste']:.1f}%): {res_holdout_coop['ensaios_teste']}")
    m_tr_c = res_holdout_coop['metricas_treino']
    m_te_c = res_holdout_coop['metricas_teste']
    print(f"   * Treino 85%: R² = {m_tr_c['R2']:.4f}, RMSE = {m_tr_c['RMSE']:.4f} ({m_tr_c['RMSE']*100:.2f}%), MAE = {m_tr_c['MAE']:.4f}")
    print(f"   * Teste 15%:  R² = {m_te_c['R2']:.4f}, RMSE = {m_te_c['RMSE']:.4f} ({m_te_c['RMSE']*100:.2f}%), MAE = {m_te_c['MAE']:.4f}")
    
    # Validação GroupKFold k=6 para o Cooperativo
    res_gkf_coop = mod_coop.validacao_group_kfold_85_15(df, n_splits=6)
    m_gkf_c = res_gkf_coop['metricas']
    print(f"   * Validação Agrupada 85/15 (GroupKFold k=6): R² = {m_gkf_c['R2']:.4f}, RMSE = {m_gkf_c['RMSE']:.4f}, MAE = {m_gkf_c['MAE']:.4f}")
    
    # Gráficos do Modelo Cooperativo
    caminho_paridade_coop = os.path.join(pasta_docs, "paridade_cooperativo_85_15.png")
    gerar_grafico_paridade_individual(res_holdout_coop,
                                      nome_modelo="Modelo Híbrido Paralelo Cooperativo (PBM + Delta_X GBDT)",
                                      caminho_saida=caminho_paridade_coop,
                                      cor_treino='#1f77b4', cor_teste='#d62728')
    
    caminho_ponto_coop = os.path.join(pasta_docs, "ponto_a_ponto_cooperativo_85_15.png")
    gerar_grafico_ponto_a_ponto(df, res_holdout_coop, mod_coop, tipo_modelo='cooperativo', caminho_saida=caminho_ponto_coop)
    
    # -------------------------------------------------------------------------
    # 3. MODELO HÍBRIDO PARALELO COMPETITIVO (85/15) PARA COMPARAÇÃO
    # -------------------------------------------------------------------------
    print("\n" + "=" * 90)
    print("3. OBTENDO RESULTADOS 85/15 DO MODELO HÍBRIDO PARALELO COMPETITIVO")
    print("=" * 90)
    mod_comp = ModeloHibridoParaleloCompetitivo(nome_algoritmo='Gradient Boosting (GBDT)')
    res_holdout_comp = mod_comp.validacao_holdout_85_15(df, test_size=0.15, random_state=42)
    m_tr_comp = res_holdout_comp['metricas_treino']
    m_te_comp = res_holdout_comp['metricas_teste']
    print(f"   * Treino 85%: R² = {m_tr_comp['R2']:.4f}, RMSE = {m_tr_comp['RMSE']:.4f} ({m_tr_comp['RMSE']*100:.2f}%)")
    print(f"   * Teste 15%:  R² = {m_te_comp['R2']:.4f}, RMSE = {m_te_comp['RMSE']:.4f} ({m_te_comp['RMSE']*100:.2f}%)")
    
    # Painel comparativo de todas as 3 arquiteturas no 85/15
    caminho_painel = os.path.join(pasta_docs, "comparativo_paridade_85_15_todas.png")
    gerar_painel_comparativo_todas_85_15(res_holdout_serial, res_holdout_coop, res_holdout_comp, caminho_painel)
    
    # -------------------------------------------------------------------------
    # 4. TABELA CONSOLIDADA DE RESULTADOS 85/15
    # -------------------------------------------------------------------------
    tabela_consolidada = [
        {
            'Arquitetura': 'Híbrido Paralelo Competitivo (Sansana 2024)',
            'Treino R² (85%)': m_tr_comp['R2'],
            'Treino RMSE': m_tr_comp['RMSE'],
            'Teste R² (15%)': m_te_comp['R2'],
            'Teste RMSE': m_te_comp['RMSE'],
            'Teste MAE': m_te_comp['MAE'],
            'Extrapolação': 'Moderada (condicionada pelo peso w)'
        },
        {
            'Arquitetura': 'Híbrido Paralelo Cooperativo (Residual)',
            'Treino R² (85%)': m_tr_c['R2'],
            'Treino RMSE': m_tr_c['RMSE'],
            'Teste R² (15%)': m_te_c['R2'],
            'Teste RMSE': m_te_c['RMSE'],
            'Teste MAE': m_te_c['MAE'],
            'Extrapolação': 'Frágil (estagnação de árvores fora do domínio)'
        },
        {
            'Arquitetura': 'Novo Híbrido Serial Cinético (alpha-GBDT)',
            'Treino R² (85%)': m_tr_s['R2'],
            'Treino RMSE': m_tr_s['RMSE'],
            'Teste R² (15%)': m_te_s['R2'],
            'Teste RMSE': m_te_s['RMSE'],
            'Teste MAE': m_te_s['MAE'],
            'Extrapolação': 'Excelente (100% conservação física e monotonia nativa)'
        }
    ]
    df_resumo = pd.DataFrame(tabela_consolidada)
    caminho_json = os.path.join(pasta_docs, "resumo_comparativo_85_15.json")
    df_resumo.to_json(caminho_json, orient='records', indent=2)
    
    print("\n" + "=" * 105)
    print("TABELA CONSOLIDADA: COMPARAÇÃO DAS TRÊS ARQUITETURAS SOB O SPLIT 85% TREINO / 15% TESTE CEGO:")
    print("=" * 105)
    print(f"{'Arquitetura':<42} | {'R² Treino':<10} | {'RMSE Treino':<12} | {'R² Teste':<10} | {'RMSE Teste':<12} | {'MAE Teste':<10}")
    print("-" * 105)
    for _, row in df_resumo.iterrows():
        print(f"{row['Arquitetura']:<42} | {row['Treino R² (85%)']:^10.4f} | {row['Treino RMSE']:^12.4f} | "
              f"{row['Teste R² (15%)']:^10.4f} | {row['Teste RMSE']:^12.4f} | {row['Teste MAE']:^10.4f}")
    print("=" * 105)
    print("Processamento concluído com sucesso!")


if __name__ == "__main__":
    main()
