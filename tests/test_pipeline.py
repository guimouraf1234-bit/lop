"""Pipeline: execução reduzida, arquivos esperados e determinismo (duas execuções idênticas)."""
import hashlib
import json
import os

import pytest

import experimentos_batelada as E
from hybrid_model import CONFIG_WB_PRINCIPAL, ModeloHibridoCinza, obter_especificacoes


def _hash_pasta(pasta):
    return {f: hashlib.sha256(open(os.path.join(pasta, f), 'rb').read()).hexdigest() for f in sorted(os.listdir(pasta))}


def test_etapas_rapidas_geram_arquivos_e_sao_deterministicas(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    for d in (a, b):
        E.main(['--rapido', '--etapas', 'dados,whitebox', '--saida', str(d), '--n-jobs', '1'])
    esperados = {'piso_ruido.csv', 'escada_whitebox.csv', 'selecao_whitebox.csv', 'termodinamica_por_ensaio.csv',
                 'config.json'}
    assert esperados <= set(os.listdir(a))
    assert _hash_pasta(str(a)) == _hash_pasta(str(b))
    cfg = json.load(open(a / 'config.json', encoding='utf-8'))
    assert cfg['white_box_principal']['forca_motriz'] == 'H2SO4'
    assert cfg['white_box_principal']['nu'] == 1.0


def test_conjuntos_disponiveis_treinam_com_o_principal(df):
    """Regressão: a etapa 'residuo' treinava o conjunto 'C' (precisa de [H+]) com o White-Box de laboratório e quebrava."""
    assert 'C' not in E.CONJUNTOS_DISPONIVEIS
    for ent in E.CONJUNTOS_DISPONIVEIS:
        m = ModeloHibridoCinza(obter_especificacoes()['gp'].nome, entradas=ent, config_wb=CONFIG_WB_PRINCIPAL)
        m.treinar(m.preparar_dados_base(df))


@pytest.mark.slow
def test_benchmark_reduzido_e_deterministico(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    for d in (a, b):
        E.main(['--rapido', '--etapas', 'dados,benchmark,aprendizado', '--saida', str(d)])
    ha, hb = _hash_pasta(str(a)), _hash_pasta(str(b))
    assert ha == hb
    assert {'benchmark_batelada.csv', 'comparacoes_pareadas.csv', 'predicoes_oof.csv.gz', 'curva_aprendizado.csv'} <= set(ha)
