import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'scripts')
import numpy as np, pandas as pd
c = pd.read_csv  # noqa
import runpy
ns = runpy.run_path('entregas/graficos_paralelo/graficos_paralelo.py')   # reaproveita curvas e medições
curvas, b, cn, med, treino = ns['curvas'], ns['b'], ns['cn'], ns['med'], ns['treino']
import validacao as V
def met(y, p):
    y, p = np.asarray(y, float), np.asarray(p, float); e = y - p
    return dict(R2=1-(e**2).sum()/((y-y.mean())**2).sum(), RMSE=np.sqrt((e**2).mean()), MAE=np.abs(e).mean(), Erro_max=np.abs(e).max())
# ponto t = 1e-6 com X_exp = 0
p0 = curvas[np.isclose(curvas.tempo_min, 1e-6)][['ensaio_id','x_hat','X_pbm']].assign(tempo_min=1e-6, X_zn_exp=0.0)
B = pd.concat([b[b.tempo_min > 0][['ensaio_id','tempo_min','X_zn_exp','x_hat']], p0[['ensaio_id','tempo_min','X_zn_exp','x_hat']]])
rows=[]
for nome, ids in [('Treino', treino), ('Teste cego', V.ENSAIOS_TESTE_CEGO), ('Global', range(1,17))]:
    p = B[B.ensaio_id.isin(ids)]; rows.append(dict(Conjunto=nome, Ensaios=len(list(ids)), Pontos=len(p), **met(p.X_zn_exp, p.x_hat)))
T1 = pd.DataFrame(rows)
q0 = cn[np.isclose(cn.tempo_min, 1e-6)][['ensaio_id','x_hat','X_pbm']].assign(tempo_min=1e-6, X_relatorio=0.0)
M = pd.concat([med[['ensaio_id','tempo_min','X_relatorio','x_hat','X_pbm','C_acid_0_mol_L','razao_molar_eta']], q0])
rows=[]
for e, g in M.groupby('ensaio_id'):
    rows.append(dict(Ensaio=e, C_A0=g.C_acid_0_mol_L.dropna().iloc[0], eta=g.razao_molar_eta.dropna().iloc[0], Pontos=len(g),
                     **met(g.X_relatorio, g.x_hat), RMSE_wb=met(g.X_relatorio, g.X_pbm)['RMSE'], R2_wb=met(g.X_relatorio, g.X_pbm)['R2']))
T2 = pd.DataFrame(rows)
g = M[M.ensaio_id != 'E1']; T3 = pd.DataFrame([dict(Conjunto='E2-E4', Pontos=len(g), **met(g.X_relatorio, g.x_hat))])
T1.to_csv('entregas/graficos_paralelo/tabela_5_4_1_com_t1e-6.csv', index=False)
T2.to_csv('entregas/graficos_paralelo/tabela_5_4_2_com_t1e-6.csv', index=False)
pd.set_option('display.float_format', lambda v: f'{v:.4f}')
print(T1.to_string(index=False)); print(); print(T2.to_string(index=False)); print(); print(T3.to_string(index=False))
print(); print(p0[p0.ensaio_id.isin([7,8,14])].to_string(index=False)); print(q0.to_string(index=False))
