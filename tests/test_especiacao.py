"""Especiação do H2SO4: balanços, limites e reprodução do [H+] de Balarini."""
import numpy as np
import pytest

from kinetics import ConfigWB, especiar_escalar, k2_hso4, k_znso4, nu_a_priori, ph_de_h


def test_k2_reproduz_h_de_balarini_a_30c():
    """6 g/L de H2SO4 (0,0612 mol/L) a 30 °C dá [H+] = 0,0689 mol/L (Balarini); tolerância de 3%."""
    s = 0.0612
    for esp in ('V0', 'V1'):
        h = especiar_escalar(s, 2 * s, k2_hso4(30.0, esp))['h']
        assert abs(h - 0.0689) / 0.0689 < 0.03, (esp, h)


def test_k2_van_t_hoff():
    assert k2_hso4(25.0, 'V1') == pytest.approx(10 ** -1.96, rel=1e-9)
    assert k2_hso4(40.0, 'V1') < k2_hso4(30.0, 'V1') < k2_hso4(25.0, 'V1')   # dissociação exotérmica
    assert k2_hso4(40.0, 'V0') == 0.0099
    assert k_znso4(40.0) > k_znso4(25.0)                                      # associação endotérmica


@pytest.mark.parametrize('esp', ['V0', 'V1', 'V2'])
def test_ph_inicial_dos_16_ensaios(resp, esp):
    """pH_0 previsto fica a menos de 0,15 do medido (a solução lixiviante não tem Zn)."""
    cfg = ConfigWB(especiacao=esp)
    for _, r in resp.iterrows():
        ca0 = float(r['C_acid_0_mol_L'])
        ph = ph_de_h(especiar_escalar(ca0, 2 * ca0, cfg.k2)['h'])
        assert abs(ph - r['pH_0']) < 0.15, (r['ensaio_id'], ph, r['pH_0'])


@pytest.mark.parametrize('S,H,zn', [(0.5, 0.4, 0.3), (1.0, 0.3, 0.8), (1.5, 2.9, 0.1), (0.1, 0.02, 0.05)])
def test_balancos_e_equilibrios_com_par_ionico(S, H, zn):
    k2, kz = 0.0072, 226.0
    e = especiar_escalar(S, H, k2, zn, kz)
    assert e['h'] + e['hso4'] == pytest.approx(H, rel=1e-9)                    # prótons
    assert e['hso4'] + e['so4'] + e['znso4'] == pytest.approx(S, rel=1e-9)     # sulfato
    assert e['h'] * e['so4'] / e['hso4'] == pytest.approx(k2, rel=1e-6)        # 2ª dissociação
    zn_livre = zn - e['znso4']
    assert e['znso4'] == pytest.approx(kz * zn_livre * e['so4'], rel=1e-6)     # par iônico
    assert min(e.values()) >= 0.0


def test_sem_zinco_v2_igual_v1():
    a = especiar_escalar(1.0, 0.6, 0.0072)
    b = especiar_escalar(1.0, 0.6, 0.0072, 0.0, 226.0)
    assert a['h'] == pytest.approx(b['h'], rel=1e-12)


def test_par_ionico_eleva_h_livre():
    """Sequestrar sulfato como ZnSO4 desloca HSO4- -> H+ + SO4(2-): [H+] não pode diminuir."""
    sem = especiar_escalar(1.0, 0.4, 0.0072)['h']
    com = especiar_escalar(1.0, 0.4, 0.0072, 0.8, 226.0)['h']
    assert com > sem


def test_valores_limite():
    assert especiar_escalar(1.0, 0.0, 0.0072)['h'] == 0.0           # sem ácido livre
    assert especiar_escalar(0.0, 0.5, 0.0072)['h'] == 0.0
    assert np.isnan(ph_de_h(0.0))


def test_nu_a_priori_e_plausivel():
    assert 1.03 <= nu_a_priori() <= 1.09
