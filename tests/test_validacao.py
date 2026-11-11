"""Protocolos de validação: ausência de vazamento, estimativa de nu só no treino, métricas e baselines."""
import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import GroupKFold

import validacao as V
from kinetics import ConfigWB, NU_A_PRIORI

CFG = ConfigWB(especiacao='V2', nu=NU_A_PRIORI)
TODOS = set(range(1, 17))


def test_esquemas_particionam_os_ensaios(df):
    for dobras in (V.esquema_gkf8_sklearn(df), V.esquema_logo(df),
                   V.esquema_deixa_nivel(df, 'razao_molar_eta'), V.esquema_deixa_nivel(df, 'C_acid_0_mol_L')):
        testes = [t for d in dobras for t in d.teste]
        assert sorted(testes) == sorted(TODOS), 'cada ensaio é testado exatamente uma vez'
        for d in dobras:
            assert not set(d.treino) & set(d.teste), 'treino e teste não podem compartilhar ensaios'
            assert set(d.treino) | set(d.teste) == TODOS
    rep = V.esquema_gkf8_repetido(df, n_rep=5, seed0=7)
    for r in range(5):
        testes = [t for d in rep if d.rep == r for t in d.teste]
        assert sorted(testes) == sorted(TODOS)
        assert all(len(d.teste) == 2 and len(d.treino) == 14 for d in rep if d.rep == r)


def test_gkf8_sklearn_reproduz_a_particao_do_tcc(df):
    ids = df['ensaio_id'].to_numpy()
    esperado = sorted(tuple(sorted(set(ids[va]))) for _, va in
                      GroupKFold(8).split(np.zeros((len(df), 1)), df['X_zn_exp'].to_numpy(), ids))
    obtido = sorted(d.teste for d in V.esquema_gkf8_sklearn(df))
    assert obtido == esperado


def test_particao_comum_igual_a_do_relatorio(df):
    """Teste cego 7, 8, 14 e as 4 dobras da Tabela B.3: o teste nunca aparece na CV nem no treino final."""
    teste = set(V.ENSAIOS_TESTE_CEGO)
    assert teste == {7, 8, 14}
    treino = set(V.ensaios_treino_comum(df))
    assert treino == TODOS - teste and len(treino) == 13
    cv = V.esquema_cv_tabela_b3(df)
    assert [d.teste for d in cv] == [(1, 5, 11, 16), (4, 10, 15), (3, 9, 13), (2, 6, 12)]
    assert sorted(t for d in cv for t in d.teste) == sorted(treino), 'cada ensaio de treino é validado uma vez'
    for d in cv:
        assert not (set(d.treino) | set(d.teste)) & teste, 'o teste cego não pode entrar na validação cruzada'
        assert set(d.treino) | set(d.teste) == treino and not set(d.treino) & set(d.teste)
    (final,) = V.esquema_teste_cego(df)
    assert set(final.treino) == treino and set(final.teste) == teste


def test_extrapolacao_temporal_separa_por_tempo(df):
    (d,) = V.esquema_extrapolacao_temporal(df, 5.0)
    assert d.t_corte == 5.0 and set(d.treino) == set(d.teste) == TODOS


def test_nuestimado_so_com_ensaios_de_treino(resp):
    assert 1.03 < V.estimar_nu(resp, list(range(1, 17))) < 1.12
    # alterar o ensaio 8 (η = 1) só importa se ele estiver no treino (mediana de poucos ensaios é sensível a ele)
    adulterada = resp.copy()
    adulterada.loc[adulterada['ensaio_id'] == 8, 'C_acid_final_mol_L'] = 0.0
    assert V.estimar_nu(adulterada, [9, 10, 6]) == V.estimar_nu(resp, [9, 10, 6])
    assert V.estimar_nu(adulterada, [8, 9, 10, 6]) != V.estimar_nu(resp, [8, 9, 10, 6])
    # a mediana é robusta a um ensaio discrepante entre vários (propriedade desejada)
    assert V.estimar_nu(adulterada, list(range(1, 17))) == pytest.approx(V.estimar_nu(resp, list(range(1, 17))), abs=0.01)
    assert V.estimar_nu(resp, [2, 7, 12, 16]) == NU_A_PRIORI      # só η = 3,1: sem ensaios elegíveis -> valor a priori


@pytest.mark.parametrize('codigo', [('wb',), ('tabela',), ('ridge', 'B'), ('hibrido', 'extra_trees', 'B'),
                                    ('hibrido', 'gp', 'B'), ('caixa', 'extra_trees'), ('caixa', 'gp')])
def test_previsao_independe_do_alvo_dos_ensaios_de_teste(df, resp, codigo):
    """Vazamento: trocar os dados MEDIDOS dos ensaios de teste (X(t), pH, X final) não pode alterar a previsão."""
    dobra = V.Dobra(treino=tuple(i for i in range(1, 17) if i not in (6, 13)), teste=(6, 13), rep=0, rotulo='teste')
    a, _ = V._rodar_dobra(df, resp, dobra, [codigo], CFG, 'dobra', None, 'vazamento', 42)
    df2, resp2 = df.copy(), resp.copy()
    m = df2['ensaio_id'].isin(dobra.teste)
    df2.loc[m, 'X_zn_exp'] = np.random.default_rng(0).uniform(0.0, 1.0, int(m.sum()))
    resp2.loc[resp2['ensaio_id'].isin(dobra.teste), ['X_zn_final', 'C_acid_final_mol_L', 'pH_final']] = [0.3, 0.9, 0.1]
    b, _ = V._rodar_dobra(df2, resp2, dobra, [codigo], CFG, 'dobra', None, 'vazamento', 42)
    assert np.array_equal(a['x_hat'].to_numpy(), b['x_hat'].to_numpy())
    assert np.array_equal(a['x_bruto'].to_numpy(), b['x_bruto'].to_numpy())
    assert (a['X_exp'].to_numpy() != b['X_exp'].to_numpy()).any()      # (o alvo mudou mesmo)


def test_extrapolacao_temporal_so_usa_tempos_de_treino(df, resp):
    """Treino com t <= 5 min: trocar os pontos de 15 min (alvo do teste) não pode alterar a previsão."""
    (dobra,) = V.esquema_extrapolacao_temporal(df, 5.0)
    codigo = ('hibrido', 'extra_trees', 'B')
    a, _ = V._rodar_dobra(df, resp, dobra, [codigo], CFG, 'a_priori', None, 'tempo', 42)
    df2 = df.copy()
    df2.loc[df2['tempo_min'] == 15.0, 'X_zn_exp'] = 0.123
    b, _ = V._rodar_dobra(df2, resp, dobra, [codigo], CFG, 'a_priori', None, 'tempo', 42)
    assert (a['tempo_min'] == 15.0).all() and len(a) == 16
    assert np.array_equal(a['x_hat'].to_numpy(), b['x_hat'].to_numpy())


def test_tabela_com_nivel_nao_visto_usa_residuo_medio_por_tempo():
    treino = pd.DataFrame({'razao_molar_eta': [1.0, 1.0, 1.5, 1.5], 'tempo_min': [1.0, 2.0, 1.0, 2.0],
                           'delta_X': [-0.1, -0.2, -0.3, -0.4], 'X_pbm': 0.9})
    t = V.TabelaEtaT().ajustar(treino)
    teste = pd.DataFrame({'razao_molar_eta': [1.0, 0.5], 'tempo_min': [1.0, 1.0], 'X_pbm': [0.9, 0.9]})
    assert np.allclose(t.prever_bruto(teste), [0.9 - 0.1, 0.9 - 0.2])      # nível visto; nível novo -> média em t = 1


def test_piso_de_ruido(df):
    pr = V.piso_ruido(df)
    assert pr['sigma'] == pytest.approx(0.0175, abs=0.0005)
    assert pr['graus_liberdade'] == 64
    assert 0.013 < pr['sigma_diferencas'] < 0.019


def _pred_sintetica():
    rng = np.random.default_rng(0)
    linhas = []
    for ens in range(1, 9):
        for t in (0.0, 1.0, 2.0, 5.0):
            x = 0.7 if t > 0 else 0.0
            for modelo, tipo, erro in (('White-Box', 'white-box', 0.10), ('Modelo bom', 'hibrido', 0.02)):
                linhas.append({'esquema': 'S', 'rep': 0, 'dobra': f'e{ens}', 'modelo': modelo, 'tipo': tipo,
                               'ensaio_id': ens, 'tempo_min': t, 'eta_nominal': 1.0, 'eta_real': 1.07, 'ca0': 0.5,
                               'X_exp': x, 'X_pbm': x, 'teto': 1.0, 'x_bruto': x, 'nu': 1.0,
                               'x_hat': x + erro * rng.choice([-1, 1])})
    return pd.DataFrame(linhas)


def test_metricas_ganho_teto_e_bootstrap():
    pred = _pred_sintetica()
    m = V.metricas_por_repeticao(pred, sigma=0.02).set_index('modelo')
    assert m.loc['White-Box', 'RMSE'] == pytest.approx(0.10)
    assert m.loc['Modelo bom', 'RMSE'] == pytest.approx(0.02)
    assert m.loc['Modelo bom', 'ganho'] == pytest.approx(1 - 0.02 ** 2 / 0.10 ** 2)
    assert m.loc['Modelo bom', 'teto_ganho'] == pytest.approx(1 - 0.02 ** 2 / 0.10 ** 2)   # já no ruído: ganho = teto
    diag = pd.DataFrame({'esquema': 'S', 'rep': 0, 'dobra': [f'e{i}' for i in range(1, 9)] * 2,
                         'modelo': ['White-Box'] * 8 + ['Modelo bom'] * 8, 'n_pts': 3, 'pct_acima_teto': 0.0,
                         'pct_abaixo_zero': 0.0, 'pct_quedas': 0.0, 'x0_max_abs': 0.0})
    r = V.resumo_benchmark(pred, diag, sigma=0.02, n_boot=200).set_index('modelo')
    assert r.loc['Modelo bom', 'RMSE_ic95_inf'] <= r.loc['Modelo bom', 'RMSE'] <= r.loc['Modelo bom', 'RMSE_ic95_sup']
    c = V.comparacoes_pareadas(pred, [('Modelo bom', 'White-Box')], n_boot=200)
    assert c.loc[0, 'delta_RMSE'] == pytest.approx(-0.08) and not bool(c.loc[0, 'empate_ic'])


def test_executar_esquema_cobre_cada_ponto_de_teste_uma_vez(df, resp):
    dobras = V.esquema_deixa_nivel(df, 'razao_molar_eta')
    pred, diag = V.executar_esquema(df, resp, 'nivel', dobras, [('wb',), ('tabela',)], CFG, 'a_priori', n_jobs=1)
    por_modelo = pred.groupby('modelo').size()
    assert (por_modelo == 128).all()
    assert pred.groupby(['modelo', 'ensaio_id', 'tempo_min']).size().max() == 1
    assert set(diag['modelo']) == {'White-Box', 'Tabela (η, t)'}
