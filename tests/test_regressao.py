"""
Regressão contra os números anteriores à revisão (congelados em resultados/baseline_antes.csv) e consistência
das constantes de hybrid_model.py com as regras de seleção do pipeline.
"""
import os

import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import GroupKFold

import validacao as V
from experimentos_batelada import escolher_entradas
from hybrid_model import (AcopladorFisico, CONFIG_WB_PRINCIPAL, ENTRADAS_PRINCIPAL, ModeloHibridoCinza,
                          construir_base)
from kinetics import ConfigWB

BASELINE = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'resultados', 'baseline_antes.csv'))


def _ref(item, coluna):
    return float(BASELINE.loc[BASELINE['item'] == item, coluna].iloc[0])


def _rmse_t_pos(base, pred):
    m = base['tempo_min'].to_numpy() > 0
    return float(np.sqrt(np.mean((base['X_zn_exp'].to_numpy()[m] - pred[m]) ** 2)))


@pytest.fixture(scope='module')
def base_legado(df):
    return construir_base(df, ConfigWB.legado())


def test_whitebox_legado_reproduz_r2_documentado(base_legado):
    y, x = base_legado['X_zn_exp'].to_numpy(), base_legado['X_pbm'].to_numpy()
    m = base_legado['tempo_min'].to_numpy() > 0
    r2 = lambda a, p: 1 - np.sum((a - p) ** 2) / np.sum((a - a.mean()) ** 2)
    assert r2(y, x) == pytest.approx(_ref('White-Box legado (PBM bruto)', 'R2_todos'), abs=5e-4)
    assert r2(y[m], x[m]) == pytest.approx(_ref('White-Box legado (PBM bruto)', 'R2_t>0'), abs=5e-4)


def test_formulacao_da_tese_reproduz_r2_documentado(df):
    for alpha, item in ((0.0, 'Tese 2017 (alpha=0)'), (5500.0, 'Tese 2017 (alpha=5500)')):
        b = construir_base(df, ConfigWB.tese(alpha))
        y, x = b['X_zn_exp'].to_numpy(), b['X_pbm'].to_numpy()
        m = b['tempo_min'].to_numpy() > 0
        assert 1 - np.sum((y - x) ** 2) / np.sum((y - y.mean()) ** 2) == pytest.approx(_ref(item, 'R2_todos'), abs=2e-3)
        assert 1 - np.sum((y[m] - x[m]) ** 2) / np.sum((y[m] - y[m].mean()) ** 2) == pytest.approx(_ref(item, 'R2_t>0'), abs=2e-3)


@pytest.mark.parametrize('nome,item,tol', [('Gaussian Process (Kriging)', 'Gaussian Process (Kriging) [GKF8]', 3e-4),
                                           ('Extra Trees', 'Extra Trees [GKF8]', 4e-4)])
def test_protocolo_legado_reproduz_rmse_do_tcc(base_legado, nome, item, tol):
    """Atributos legados + acoplador legado + GroupKFold(8) do sklearn: números de antes da revisão."""
    modelo = ModeloHibridoCinza(nome, entradas='legado')
    p = modelo.validacao_cruzada_group_kfold(base_legado, 8)
    assert _rmse_t_pos(base_legado, p) == pytest.approx(_ref(item, 'RMSE_t>0'), abs=tol)


def test_tabela_legada_reproduz_rmse_do_tcc(base_legado):
    acop = AcopladorFisico(modo='legado')
    ids = base_legado['ensaio_id'].to_numpy()
    p = np.zeros(len(base_legado))
    for tr, va in GroupKFold(8).split(np.zeros((len(base_legado), 1)), groups=ids):
        t = V.TabelaEtaT().ajustar(base_legado.iloc[tr])
        te = base_legado.iloc[va]
        p[va] = acop.restringir(te, t.prever_bruto(te))
    assert _rmse_t_pos(base_legado, p) == pytest.approx(_ref('Tabela media por (eta, t) [GKF8]', 'RMSE_t>0'), abs=3e-4)


def test_teto_fisico_remove_a_correcao_oculta_do_whitebox(base_legado):
    """O teto com η nominal reduzia o RMSE do White-Box de 0,0443 para 0,0415; com η_real/ν isso deixa de ocorrer."""
    zero = np.zeros(len(base_legado))
    legado = AcopladorFisico(modo='legado').acoplar_e_restringir(base_legado, zero)
    fisico = AcopladorFisico(nu=1.0).acoplar_e_restringir(base_legado, zero)
    assert _rmse_t_pos(base_legado, legado) == pytest.approx(_ref('White-Box legado + acoplador legado', 'RMSE_t>0'), abs=3e-4)
    assert _rmse_t_pos(base_legado, fisico) == pytest.approx(_ref('White-Box legado (PBM bruto)', 'RMSE_t>0'), abs=3e-4)


def test_whitebox_principal_e_o_de_laboratorio():
    """Contrato do White-Box principal: força motriz CAf livre, nu = 1, alpha = 0, ks(T) de Balarini, RRB mássica."""
    assert CONFIG_WB_PRINCIPAL == ConfigWB.laboratorio()
    assert (CONFIG_WB_PRINCIPAL.forca_motriz, CONFIG_WB_PRINCIPAL.nu, CONFIG_WB_PRINCIPAL.alpha) == ('H2SO4', 1.0, 0.0)
    assert CONFIG_WB_PRINCIPAL.ks_arrhenius and CONFIG_WB_PRINCIPAL.rrb == 'massica'


def test_whitebox_principal_nao_le_dados_de_conversao_nem_de_ph(df):
    """Sem ajuste aos dados: a simulação não muda quando X medido e pH medido são alterados."""
    alterado = df.copy()
    alterado['X_zn_exp'] = 0.5
    for col in ('pH_0', 'pH_final'):
        if col in alterado.columns:
            alterado[col] = 3.0
    a = construir_base(df, CONFIG_WB_PRINCIPAL, usar_cache=False)
    b = construir_base(alterado, CONFIG_WB_PRINCIPAL, usar_cache=False)
    assert np.array_equal(a['X_pbm'].to_numpy(), b['X_pbm'].to_numpy())
    assert np.array_equal(a['C_acid_pbm_mol_L'].to_numpy(), b['C_acid_pbm_mol_L'].to_numpy())


def test_acido_livre_do_principal_e_o_balanco_do_tcc(df):
    """CAf = CA0 (1 - X/eta): consumo de ácido só o estequiométrico (nu = 1)."""
    b = construir_base(df, CONFIG_WB_PRINCIPAL)
    eta = b['razao_molar_eta_real'].to_numpy()
    esperado = np.maximum(b['C_acid_0_mol_L'].to_numpy() * (1.0 - b['X_pbm'].to_numpy() / eta), 0.0)
    assert np.allclose(b['C_acid_pbm_mol_L'].to_numpy(), esperado, atol=1e-9)


def test_regra_de_escolha_dos_atributos_prefere_o_mais_simples_em_caso_de_empate():
    resumo = pd.DataFrame({'entradas': ['A', 'B', 'C'] * 2, 'modelo': ['m1'] * 3 + ['m2'] * 3,
                           'RMSE': [0.040, 0.021, 0.0205, 0.041, 0.022, 0.0215]})
    rng = np.random.default_rng(0)
    base = rng.normal(0.021, 0.001, 10)
    por_rep = pd.DataFrame({'A': base + 0.02, 'B': base + 0.0003 + rng.normal(0, 0.001, 10), 'C': base})
    sel = escolher_entradas(resumo, por_rep)
    assert sel['melhor_por_rmse'] == 'C' and sel['escolhido'] == 'B'
    # se B for claramente pior que C em todas as repetições, escolhe C
    por_rep['B'] = base + 0.004
    assert escolher_entradas(resumo, por_rep)['escolhido'] == 'C'
