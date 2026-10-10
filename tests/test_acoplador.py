"""Acoplador físico: limites, X(0) = 0, monotonia isotônica e invariância à ordem das linhas."""
import numpy as np
import pandas as pd
import pytest

from hybrid_model import AcopladorFisico, construir_base
from kinetics import ConfigWB, NU_A_PRIORI

CFG = ConfigWB(especiacao='V2', nu=NU_A_PRIORI)


@pytest.fixture(scope='module')
def base(df):
    return construir_base(df, CFG)


def _bruto_ruidoso(base, semente=0, amp=0.08):
    rng = np.random.default_rng(semente)
    return base['X_pbm'].to_numpy() + rng.normal(0.0, amp, len(base))


def test_saida_respeita_teto_zero_e_monotonia(base):
    acop = AcopladorFisico(nu=CFG.nu)
    x = acop.restringir(base, _bruto_ruidoso(base))
    teto = np.minimum(1.0, base['razao_molar_eta_real'].to_numpy() / CFG.nu)
    assert np.all(x >= 0.0) and np.all(x <= teto + 1e-12)
    assert np.all(x[base['tempo_min'].to_numpy() == 0.0] == 0.0)
    for _, g in base.assign(x=x).groupby('ensaio_id'):
        assert np.all(np.diff(g.sort_values('tempo_min')['x'].to_numpy()) >= -1e-12)


def test_invariante_a_ordem_das_linhas(base):
    acop = AcopladorFisico(nu=CFG.nu)
    bruto = _bruto_ruidoso(base, 1)
    ref = acop.restringir(base, bruto)
    perm = np.random.default_rng(3).permutation(len(base))
    emb = acop.restringir(base.iloc[perm].reset_index(drop=True), bruto[perm])
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    assert np.allclose(emb[inv], ref, atol=1e-12)


def test_isotonica_nao_e_viesada_para_cima_como_o_maximo_acumulado():
    """Ruído simétrico em torno de um patamar: a projeção isotônica não infla a média, o cummax infla."""
    t = np.arange(0.0, 20.0)
    rng = np.random.default_rng(0)
    bruto = 0.8 + rng.normal(0.0, 0.02, len(t))
    bruto[0] = 0.0
    q = pd.DataFrame({'ensaio_id': 1, 'tempo_min': t, 'razao_molar_eta_real': 2.0, 'razao_molar_eta': 2.0,
                      'X_pbm': 0.0})
    iso = AcopladorFisico(modo='fisico').restringir(q, bruto)
    cummax = np.maximum.accumulate(np.clip(bruto, 0, 1))
    assert abs(iso[5:].mean() - 0.8) < abs(cummax[5:].mean() - 0.8)


def test_modo_legado_usa_eta_nominal_e_maximo_acumulado():
    q = pd.DataFrame({'ensaio_id': 1, 'tempo_min': [0.0, 1.0, 2.0, 3.0], 'razao_molar_eta': 0.5,
                      'razao_molar_eta_real': 0.535, 'X_pbm': 0.0})
    x = AcopladorFisico(modo='legado').restringir(q, np.array([0.0, 0.45, 0.40, 0.70]))
    assert np.allclose(x, [0.0, 0.45, 0.45, 0.50])          # teto 0,5 (nominal) e cummax


def test_teto_fisico_usa_eta_real_sobre_nu():
    q = pd.DataFrame({'ensaio_id': [1, 1], 'tempo_min': [0.0, 5.0], 'razao_molar_eta': 0.5,
                      'razao_molar_eta_real': 0.535, 'X_pbm': 0.0})
    assert AcopladorFisico(nu=1.07).teto(q)[0] == pytest.approx(0.5)
    assert AcopladorFisico(nu=1.0).teto(q)[0] == pytest.approx(0.535)


def test_diagnostico_conta_violacoes_da_saida_bruta():
    q = pd.DataFrame({'ensaio_id': 1, 'tempo_min': [0.0, 1.0, 2.0, 3.0, 4.0], 'razao_molar_eta': 1.0,
                      'razao_molar_eta_real': 1.0, 'X_pbm': 0.0})
    bruto = np.array([0.03, 0.5, 0.45, 1.2, -0.1])     # acima do teto em t=3; abaixo de 0 em t=4; quedas 1->2 e 3->4
    d = AcopladorFisico(nu=1.0).diagnosticar(q, bruto)
    assert d['pct_acima_teto'] == pytest.approx(25.0)
    assert d['pct_abaixo_zero'] == pytest.approx(25.0)
    assert d['pct_quedas'] == pytest.approx(50.0)
    assert d['x0_max_abs'] == pytest.approx(0.03)


def test_exige_eta_real(base):
    sem = base.drop(columns=['razao_molar_eta_real'])
    with pytest.raises(KeyError):
        AcopladorFisico(nu=1.0).restringir(sem, np.zeros(len(sem)))
