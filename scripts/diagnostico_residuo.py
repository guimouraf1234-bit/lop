"""
diagnostico_residuo.py
----------------------
Exporta a tabela do resíduo do White-Box, ΔX = X_exp - X_pbm, para uma configuração dada (padrão: o White-Box
principal). É só um EXPORTADOR para inspeção: o pipeline não lê este arquivo e sempre recalcula o PBM.

    python scripts/diagnostico_residuo.py                  # White-Box principal -> resultados/dados_com_residuo.csv
    python scripts/diagnostico_residuo.py --legado         # White-Box do TCC antes da revisão (V0, nu = 1)
"""

import argparse
import os
import sys

DIR_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
if DIR_SCRIPTS not in sys.path:
    sys.path.insert(0, DIR_SCRIPTS)

from data_loader import carregar_dados_bancada  # noqa: E402
from hybrid_model import CONFIG_WB_PRINCIPAL, construir_base  # noqa: E402
from kinetics import ConfigWB  # noqa: E402


def calcular_residuos(config: ConfigWB = CONFIG_WB_PRINCIPAL, caminho_saida: str = None):
    """Quadro com X_pbm, C_acid_pbm_mol_L, H_pbm_mol_L e delta_X (e grava o CSV se `caminho_saida`)."""
    df = construir_base(carregar_dados_bancada(), config)
    if caminho_saida:
        os.makedirs(os.path.dirname(os.path.abspath(caminho_saida)), exist_ok=True)
        df.to_csv(caminho_saida, index=False, float_format='%.8g')
        print(f'Resíduos ({config.rotulo()}) salvos em: {caminho_saida}')
    return df


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--legado', action='store_true', help='usa o White-Box anterior à revisão (V0, nu = 1)')
    ap.add_argument('--saida', default=os.path.join(DIR_SCRIPTS, '..', 'resultados', 'dados_com_residuo.csv'))
    a = ap.parse_args()
    calcular_residuos(ConfigWB.legado() if a.legado else CONFIG_WB_PRINCIPAL, a.saida)
