"""
Script de Treinamento e Avaliação Comparativa do Modelo MLP com Data Augmentation

Compara:
1. MLP Tanh (Original - 128 pontos, alpha=0.1)
2. MLP Tanh (Augmented - 640 pontos, alpha=0.01)
3. SVR Otimizado (Referência - 128 pontos)

Avaliações realizadas:
- Ajuste Completo (R², RMSE, MAE avaliados nos 128 pontos experimentais reais)
- Validação Cruzada Agrupada (GroupKFold com 8 dobras / 85-15):
  * O modelo aumentado treina com 560 pontos gerados em 14 ensaios
  * O teste cego avalia estritamente os 16 pontos experimentais reais dos 2 ensaios deixados de fora
"""

import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.base import clone

diretorio_atual = os.path.dirname(os.path.abspath(__file__))
diretorio_lop = os.path.dirname(diretorio_atual)
diretorio_scripts = os.path.join(diretorio_lop, "scripts")
if diretorio_scripts not in sys.path:
    sys.path.insert(0, diretorio_scripts)

from data_loader import carregar_dados_bancada
from population_balance import PopulationBalanceModel
from hybrid_model import (
    ModeloHibridoCinza,
    preparar_features_e_alvo,
    calcular_metricas,
    AcopladorFisico,
    FEATURES_FISICAS
)
from data_augmentation import gerar_base_aumentada


def construir_pipeline_mlp(alpha: float = 0.01, random_state: int = 42) -> Pipeline:
    return Pipeline([
        ('scaler', StandardScaler()),
        ('mlp', MLPRegressor(
            hidden_layer_sizes=(12, 6),
            activation='tanh',
            solver='lbfgs',
            alpha=alpha,
            max_iter=3000,
            random_state=random_state
        ))
    ])


def executar_benchmark():
    print("=" * 80)
    print("BENCHMARK: MLP (ORIGINAL) vs MLP (COM DATA AUGMENTATION) vs SVR (REFERÊNCIA)")
    print("=" * 80)

    df_exp = carregar_dados_bancada()
    pbm = PopulationBalanceModel()
    acoplador = AcopladorFisico()

    # 1. Base real com simulações White-box
    h_temp = ModeloHibridoCinza('SVR (Otimizado)')
    df_base_real = h_temp.preparar_dados_base(df_exp)

    # 2. Base aumentada (640 pontos)
    df_aug = gerar_base_aumentada(salvar_csv=True)

    ensaios = df_base_real['ensaio_id'].unique()
    gkf = GroupKFold(n_splits=8)

    # -------------------------------------------------------------
    # A) VALIDAÇÃO CRUZADA GROUPKFOLD (8 DOBRAS / 85-15)
    # -------------------------------------------------------------
    print("\n-> Executando Validação Cruzada GroupKFold (8 Dobras)...")
    
    # MLP Original (sem augmentation, alpha=0.1)
    mlp_orig = ModeloHibridoCinza('Rede Neural (MLP Tanh + L-BFGS)')
    y_cv_mlp_orig = mlp_orig.validacao_cruzada_group_kfold(df_base_real, n_splits=8)
    m_cv_mlp_orig = calcular_metricas(df_base_real['X_zn_exp'].values, y_cv_mlp_orig)

    # SVR Otimizado (Referência sem augmentation)
    svr_ref = ModeloHibridoCinza('SVR (Otimizado)')
    y_cv_svr = svr_ref.validacao_cruzada_group_kfold(df_base_real, n_splits=8)
    m_cv_svr = calcular_metricas(df_base_real['X_zn_exp'].values, y_cv_svr)

    # MLP com Data Augmentation (alpha=0.01)
    y_cv_mlp_aug = np.zeros(len(df_base_real))
    pipe_aug = construir_pipeline_mlp(alpha=0.01)

    for train_ens_idx, val_ens_idx in gkf.split(ensaios, groups=ensaios):
        train_ens = ensaios[train_ens_idx]
        val_ens = ensaios[val_ens_idx]

        # Treino: dados aumentados (14 ensaios x 40 = 560 pontos)
        df_train_fold = df_aug[df_aug['ensaio_id'].isin(train_ens)]
        X_tr = df_train_fold[FEATURES_FISICAS].values
        y_tr = df_train_fold['delta_X'].values

        # Teste cego: dados REAIS puros (2 ensaios x 8 = 16 pontos)
        val_mask = df_base_real['ensaio_id'].isin(val_ens)
        df_val_fold = df_base_real[val_mask].copy()
        X_val = df_val_fold[FEATURES_FISICAS].values

        pipe_fold = clone(pipe_aug)
        pipe_fold.fit(X_tr, y_tr)

        delta_pred_val = pipe_fold.predict(X_val)
        y_cv_mlp_aug[val_mask] = acoplador.acoplar_e_restringir(df_val_fold, delta_pred_val)

    m_cv_mlp_aug = calcular_metricas(df_base_real['X_zn_exp'].values, y_cv_mlp_aug)

    # -------------------------------------------------------------
    # B) AJUSTE COMPLETO (100% DOS DADOS DE TREINO)
    # -------------------------------------------------------------
    print("\n-> Executando Ajuste Completo (avaliado nos 128 pontos reais)...")
    
    # MLP Original Fit
    mlp_orig.treinar(df_base_real)
    y_fit_mlp_orig = mlp_orig.prever(df_base_real)
    m_fit_mlp_orig = calcular_metricas(df_base_real['X_zn_exp'].values, y_fit_mlp_orig)

    # SVR Fit
    svr_ref.treinar(df_base_real)
    y_fit_svr = svr_ref.prever(df_base_real)
    m_fit_svr = calcular_metricas(df_base_real['X_zn_exp'].values, y_fit_svr)

    # MLP com Augmentation Fit
    pipe_aug_full = construir_pipeline_mlp(alpha=0.01)
    pipe_aug_full.fit(df_aug[FEATURES_FISICAS].values, df_aug['delta_X'].values)
    delta_fit_pred = pipe_aug_full.predict(df_base_real[FEATURES_FISICAS].values)
    y_fit_mlp_aug = acoplador.acoplar_e_restringir(df_base_real, delta_fit_pred)
    m_fit_mlp_aug = calcular_metricas(df_base_real['X_zn_exp'].values, y_fit_mlp_aug)

    # -------------------------------------------------------------
    # C) EXIBIÇÃO DA TABELA COMPARATIVA
    # -------------------------------------------------------------
    tabela = [
        {
            "Modelo": "MLP Tanh (Original - 128 pts)",
            "N_treino": "112 (CV) / 128 (Full)",
            "R2_Fit": m_fit_mlp_orig['R2'],
            "MAE_Fit (%)": f"{m_fit_mlp_orig['MAE']*100:.2f}%",
            "R2_CV": m_cv_mlp_orig['R2'],
            "MAE_CV (%)": f"{m_cv_mlp_orig['MAE']*100:.2f}%",
            "RMSE_CV": m_cv_mlp_orig['RMSE']
        },
        {
            "Modelo": "MLP Tanh + Augmentation (640 pts)",
            "N_treino": "560 (CV) / 640 (Full)",
            "R2_Fit": m_fit_mlp_aug['R2'],
            "MAE_Fit (%)": f"{m_fit_mlp_aug['MAE']*100:.2f}%",
            "R2_CV": m_cv_mlp_aug['R2'],
            "MAE_CV (%)": f"{m_cv_mlp_aug['MAE']*100:.2f}%",
            "RMSE_CV": m_cv_mlp_aug['RMSE']
        },
        {
            "Modelo": "SVR Otimizado (Kernel RBF)",
            "N_treino": "112 (CV) / 128 (Full)",
            "R2_Fit": m_fit_svr['R2'],
            "MAE_Fit (%)": f"{m_fit_svr['MAE']*100:.2f}%",
            "R2_CV": m_cv_svr['R2'],
            "MAE_CV (%)": f"{m_cv_svr['MAE']*100:.2f}%",
            "RMSE_CV": m_cv_svr['RMSE']
        }
    ]

    df_res = pd.DataFrame(tabela)
    print("\n" + "=" * 80)
    print("RESULTADOS CONSOLIDADOS DO EXPERIMENTO DE AUGMENTATION")
    print("=" * 80)
    print(df_res.to_string(index=False))

    return {
        'df_res': df_res,
        'pipe_aug_full': pipe_aug_full,
        'df_base_real': df_base_real,
        'df_aug': df_aug
    }


if __name__ == '__main__':
    executar_benchmark()
