"""Avalia E1–E4 sem ajuste aos ensaios novos, na base de 0,611 g Zn/g de sólido e com η corrigida pelas massas."""
from pathlib import Path
import json
import warnings
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_config, dump
import gerar_resultados_gp_relatorio as G
from data_loader import carregar_dados_bancada
from hybrid_model import construir_base, construir_base_condicao, CONFIG_WB_PRINCIPAL as CFG
from experimentos_batelada import _metricas_pp
import validacao as V


def ler_celulas(path):
    ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
        return {c.attrib['r']: float(c.find('s:v', ns).text)
                for c in root.findall('.//s:sheetData/s:row/s:c', ns)
                if c.find('s:v', ns) is not None and c.attrib.get('t') not in ('s', 'str')}


def treinar(seed, tr):
    warnings.filterwarnings('ignore')
    return V.criar_estimador(('hibrido', 'gp', 'B'), seed).ajustar(tr)


def main():
    warnings.filterwarnings('ignore')
    root = Path(__file__).resolve().parents[1]
    cells = ler_celulas(root / 'entregas/conversao_ensaios_novos.xlsx')
    # Identificação E2/E3 já trocada no arquivo da branch, conforme nota em A48.
    condicoes = [('E1', .25, .75, 6.7, 11), ('E2', .75, 1.25, 12.1, 21),
                 ('E3', .75, 2.25, 6.7, 31), ('E4', .25, .40, 12.6, 41)]
    teor_zn = .611
    malhas, medidas, condicoes_efetivas = [], [], []
    for eid, ca0, nominal, massa, inicio in condicoes:
        idx = len(malhas) + 1
        assert cells[f'F{idx+1}'] == ca0 and cells[f'G{idx+1}'] == massa
        eta_efetiva = .2 * ca0 * 65.38 / (massa * teor_zn)
        condicoes_efetivas.append({'ensaio': eid, 'eta_planejada_Zn_total': nominal, 'eta_ZnO': eta_efetiva})
        malhas.append(construir_base_condicao(eta_efetiva, ca0,
                                              np.arange(0, 15.001, .25), CFG,
                                              ensaio_id=eid, eta_nominal=eta_efetiva))
        for row in range(inicio, inicio + 6):
            # Zn(g/L) = leitura × massa após diluição / massa da alíquota.
            zn = cells[f'E{row}'] * cells[f'C{row}'] / cells[f'B{row}']
            medidas.append({'ensaio_id': eid, 'tempo_min': cells[f'A{row}'],
                             'C_acid_0_mol_L': ca0, 'razao_molar_eta': eta_efetiva,
                             'eta_planejada_Zn_total': nominal,
                             'massa_solido_g': massa, 'Zn_g_L': zn,
                             'X_relatorio': zn * .2 / (massa * teor_zn)})
    malha = pd.concat(malhas, ignore_index=True)
    base = construir_base(carregar_dados_bancada(), CFG)
    treino = V.ensaios_treino_comum(base)
    print('GP: mesmas dez sementes, apenas os 13 ensaios originais no ajuste.', flush=True)
    with parallel_config(backend='loky', n_jobs=4, inner_max_num_threads=1):
        modelos = Parallel()(delayed(treinar)(seed, base[base.ensaio_id.isin(treino)]) for seed in G.SEEDS)
    dump(modelos, G.OUT / 'modelos_gp_13_ensaios.joblib')
    # Confere reprodução do comitê que gerou os resultados anteriores.
    anterior = pd.read_csv(G.OUT / 'curvas_gp.csv')
    bruto_ant = np.mean([m.prever_bruto(anterior) for m in modelos], axis=0)
    assert np.max(np.abs(bruto_ant - anterior.x_bruto)) < 1e-8
    malha['x_bruto'] = np.mean([m.prever_bruto(malha) for m in modelos], axis=0)
    malha['x_hat'] = G.ACOP.restringir(malha, malha.x_bruto)
    malha.to_csv(G.OUT / 'curvas_ensaios_novos_gp.csv', index=False)
    p = pd.DataFrame(medidas)
    p['x_hat'] = G.pontos_da_malha(p, malha, malha.x_hat)
    p['X_pbm'] = G.pontos_da_malha(p, malha, malha.X_pbm)
    p.to_csv(G.OUT / 'predicoes_ensaios_novos_gp.csv', index=False)
    linhas = []
    for versao, teor in [('relatorio', teor_zn)]:
        d = p.assign(X_zn_exp=p['X_' + versao])
        G.painel(d, malha, ['E1', 'E2', 'E3', 'E4'], (2, 2), 'figura_4_ensaios_novos_gp_' + versao)
        import shutil
        for ext in ('png', 'pdf'):
            shutil.copy2(G.OUT / f'figura_4_ensaios_novos_gp_{versao}.{ext}',
                         G.OUT / f'figura_4_ensaios_novos_gp_0611_corrigido.{ext}')
        for e, g in d.groupby('ensaio_id'):
            linhas.append({'base': versao, 'teor_Zn_g_g': teor, 'ensaio': e,
                           **_metricas_pp(g.X_zn_exp, g.x_hat),
                           'eta_ZnO': float(g.iloc[0].razao_molar_eta),
                           'medicoes_acima_teto': int((g.X_zn_exp > min(1, g.iloc[0].razao_molar_eta)).sum()),
                           'RMSE_whitebox_pp': _metricas_pp(g.X_zn_exp, g.X_pbm)['RMSE_pp'],
                           'medida_15min': float(g.iloc[-1].X_zn_exp),
                           'prevista_15min': float(g.iloc[-1].x_hat)})
    pd.DataFrame(linhas).to_csv(G.OUT / 'metricas_ensaios_novos_gp.csv', index=False)
    (G.OUT / 'config_ensaios_novos.json').write_text(json.dumps({
        'whitebox': CFG.para_dict(), 'treino': list(treino), 'sementes': list(G.SEEDS),
        'ajuste_aos_ensaios_novos': False, 'volume_L': .2,
        'teor_zn_g_g': teor_zn, 'teor_Zn_total_usado_na_pesagem': .650,
        'base_eta_avaliacao': 'V * CA0 * M_Zn / (massa * 0,611)',
        'condicoes_efetivas': condicoes_efetivas,
        'condicoes': condicoes, 'rotulos': 'E2 e E3 conforme arquivo já trocado, nota A48',
        'unidade_leituras_assumida': 'g/L conforme tratamento do relatório',
        'projecao': 'malha 0,25 min; métricas extraídas da mesma curva'},
        ensure_ascii=False, indent=2), encoding='utf-8')
    print(pd.DataFrame(linhas).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
