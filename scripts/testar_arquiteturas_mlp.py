import sys, os
diretorio_scripts = os.path.dirname(os.path.abspath(__file__))
if diretorio_scripts not in sys.path:
    sys.path.append(diretorio_scripts)

from data_loader import carregar_dados_bancada
from hybrid_model import (
    ModeloHibridoCinza, preparar_features_e_alvo,
    calcular_metricas, AcopladorFisico
)
from sklearn.model_selection import LeaveOneGroupOut, GroupKFold
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import clone
import numpy as np
import pandas as pd

df_exp = carregar_dados_bancada()
hibrido = ModeloHibridoCinza('Rede Neural (MLP Tanh + L-BFGS)')
df_base = hibrido.preparar_dados_base(df_exp)
X, y, grupos = preparar_features_e_alvo(df_base)
acoplador = AcopladorFisico()

propostas = [
    {'nome': 'Arch 1 (Original TCC: 12-6, a=0.10)',      'layers': (12, 6),    'alpha': 0.10},
    {'nome': 'Arch 2 (Rasa 1-Camada: 16, a=0.05)',        'layers': (16,),      'alpha': 0.05},
    {'nome': 'Arch 3 (Rasa Compacta: 12, a=0.03)',        'layers': (12,),      'alpha': 0.03},
    {'nome': 'Arch 4 (2-Camadas Compacta: 10-5, a=0.05)', 'layers': (10, 5),    'alpha': 0.05},
    {'nome': 'Arch 5 (2-Camadas Piramidal: 14-7, a=0.05)','layers': (14, 7),    'alpha': 0.05},
    {'nome': 'Arch 6 (Original Ajustado: 12-6, a=0.05)',  'layers': (12, 6),    'alpha': 0.05},
]

print("=" * 80)
print("TESTE EM LEAVE-ONE-GROUP-OUT (LOGO-CV: 15 Treino / 1 Teste)")
print("=" * 80)

logo = LeaveOneGroupOut()
for p in propostas:
    y_cv = np.zeros(len(df_base))
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        ('mlp', MLPRegressor(
            hidden_layer_sizes=p['layers'],
            activation='tanh',
            solver='lbfgs',
            alpha=p['alpha'],
            max_iter=3000,
            random_state=42
        ))
    ])
    for train_idx, val_idx in logo.split(X, y, grupos):
        m = clone(pipe)
        m.fit(X[train_idx], y[train_idx])
        pred = m.predict(X[val_idx])
        df_val = df_base.iloc[val_idx].copy()
        y_cv[val_idx] = acoplador.acoplar_e_restringir(df_val, pred)
    m_cv = calcular_metricas(df_base['X_zn_exp'].values, y_cv)
    print(f"{p['nome']:<45} | R2 = {m_cv['R2']:.4f} | RMSE = {m_cv['RMSE']:.4f} | MAE = {m_cv['MAE']*100:.2f}%")

print("\n" + "=" * 80)
print("TESTE EM GROUP-KFOLD (8 Dobras: 14 Treino / 2 Teste)")
print("=" * 80)

gkf = GroupKFold(n_splits=8)
for p in propostas:
    y_cv = np.zeros(len(df_base))
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        ('mlp', MLPRegressor(
            hidden_layer_sizes=p['layers'],
            activation='tanh',
            solver='lbfgs',
            alpha=p['alpha'],
            max_iter=3000,
            random_state=42
        ))
    ])
    for train_idx, val_idx in gkf.split(X, y, grupos):
        m = clone(pipe)
        m.fit(X[train_idx], y[train_idx])
        pred = m.predict(X[val_idx])
        df_val = df_base.iloc[val_idx].copy()
        y_cv[val_idx] = acoplador.acoplar_e_restringir(df_val, pred)
    m_cv = calcular_metricas(df_base['X_zn_exp'].values, y_cv)
    print(f"{p['nome']:<45} | R2 = {m_cv['R2']:.4f} | RMSE = {m_cv['RMSE']:.4f} | MAE = {m_cv['MAE']*100:.2f}%")
