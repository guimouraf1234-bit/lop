"""Configuração comum dos testes: caminho dos módulos, dados e marcadores."""
import os
import sys
import warnings

import pytest

DIR_TESTES = os.path.dirname(os.path.abspath(__file__))
DIR_SCRIPTS = os.path.abspath(os.path.join(DIR_TESTES, '..', 'scripts'))
if DIR_SCRIPTS not in sys.path:
    sys.path.insert(0, DIR_SCRIPTS)
os.environ['PYTHONPATH'] = DIR_SCRIPTS + os.pathsep + os.environ.get('PYTHONPATH', '')
warnings.filterwarnings('ignore')

from data_loader import carregar_dados_bancada, carregar_respostas_finais  # noqa: E402


@pytest.fixture(scope='session')
def df():
    return carregar_dados_bancada()


@pytest.fixture(scope='session')
def resp():
    return carregar_respostas_finais()


def pytest_collection_modifyitems(config, items):
    if os.environ.get('RUN_SLOW') == '1':
        return
    pular = pytest.mark.skip(reason='teste lento: defina RUN_SLOW=1')
    for item in items:
        if 'slow' in item.keywords:
            item.add_marker(pular)
