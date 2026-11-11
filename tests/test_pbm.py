"""Invariantes e equivalências do PBM (White-Box)."""
import os

import numpy as np
import pandas as pd
import pytest

from kinetics import ConfigWB, NU_A_PRIORI
from population_balance import PopulationBalanceModel

DIR_DADOS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dados')
CONFIGS = {
    'laboratorio': ConfigWB.laboratorio(),
    'legado': ConfigWB.legado(),
    'V1_nu': ConfigWB(especiacao='V1', nu=NU_A_PRIORI),
    'V2_nu': ConfigWB(especiacao='V2', nu=NU_A_PRIORI),
    'tese': ConfigWB.tese(0.0),
}
TEMPOS = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0]


def _sim(cfg, ca0, eta, rapido=True, tempos=TEMPOS):
    pbm = PopulationBalanceModel.para_config(cfg)
    f = pbm.simular_ensaio_rapido if rapido else pbm.simular_ensaio
    return f(tempos, ca0, eta, config=cfg)


@pytest.mark.parametrize('nome', list(CONFIGS))
def test_invariantes_fisicas(df, nome):
    """X(0) = 0; X não decresce; X <= min(1, eta/nu); CAf >= 0 e decresce; H+ >= 0."""
    cfg = CONFIGS[nome]
    for _, g in df.groupby('ensaio_id'):
        ca0, eta = float(g['C_acid_0_mol_L'].iloc[0]), float(g['razao_molar_eta_real'].iloc[0])
        r = _sim(cfg, ca0, eta)
        x, caf = r['X_pbm'], r['C_acid_mol_L']
        assert x[0] == 0.0
        assert np.all(np.diff(x) >= -1e-9), 'X deve ser não decrescente'
        assert x.max() <= min(1.0, eta / cfg.nu) + 2e-3, 'X não pode exceder o limite estequiométrico'
        assert np.all(caf >= 0.0) and np.all(np.diff(caf) <= 1e-9)
        if cfg.forca_motriz == 'H+':
            assert np.all(r['C_H_mol_L'] >= 0.0)


@pytest.mark.parametrize('nome', ['laboratorio', 'legado', 'V2_nu', 'tese'])
def test_integrador_rapido_equivale_a_referencia(df, nome):
    """O integrador por tabela difere da referência em menos de 1e-4 em X."""
    cfg = CONFIGS[nome]
    pior = 0.0
    for _, g in df.groupby('ensaio_id'):
        ca0, eta = float(g['C_acid_0_mol_L'].iloc[0]), float(g['razao_molar_eta_real'].iloc[0])
        a = _sim(cfg, ca0, eta, rapido=False)['X_pbm']
        b = _sim(cfg, ca0, eta, rapido=True)['X_pbm']
        pior = max(pior, float(np.abs(a - b).max()))
    assert pior < 1e-4, pior


def test_referencia_reproduz_resultados_legados():
    """Com a configuração legada, a integração de referência reproduz o CSV congelado do PBM antigo."""
    ref = pd.read_csv(os.path.join(DIR_DADOS, 'pbm_legado_referencia.csv'))
    pbm = PopulationBalanceModel()
    for ens, g in ref.groupby('ensaio_id', sort=False):
        r = pbm.simular_ensaio(g['tempo_min'].tolist(), float(g['C_acid_0_mol_L'].iloc[0]),
                               float(g['razao_molar_eta_real'].iloc[0]))
        assert np.allclose(r['X_pbm'], g['X_pbm'].to_numpy(), atol=1e-9)
        assert np.allclose(r['C_acid_mol_L'], g['C_acid_pbm_mol_L'].to_numpy(), atol=1e-9)


def test_nu_limita_a_conversao_com_acido_limitante():
    """Com eta < 1 e nu > 1 o patamar é eta/nu (ácido totalmente consumido), não eta."""
    eta = 0.535
    x_nu1 = _sim(ConfigWB.legado(), 0.5, eta, tempos=[0.0, 15.0, 60.0])['X_pbm'][-1]
    x_nu = _sim(ConfigWB(nu=1.07), 0.5, eta, tempos=[0.0, 15.0, 60.0])['X_pbm'][-1]
    assert abs(x_nu1 - eta) < 0.01
    assert abs(x_nu - eta / 1.07) < 0.01


def test_tempos_fora_de_ordem_sao_devolvidos_na_ordem_pedida():
    cfg = CONFIGS['legado']
    a = _sim(cfg, 0.5, 1.07, tempos=[15.0, 0.0, 2.0, 0.5])
    b = _sim(cfg, 0.5, 1.07, tempos=[0.0, 0.5, 2.0, 15.0])
    assert np.allclose(a['X_pbm'][[1, 3, 2, 0]], b['X_pbm'])


def test_config_incompativel_com_instancia_levanta_erro():
    with pytest.raises(ValueError):
        PopulationBalanceModel(formulacao='tese').simular_ensaio([0.0, 1.0], 0.5, 1.0, config=ConfigWB.legado())
