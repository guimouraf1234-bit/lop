"""
tabelas_batelada.py
-------------------
Converte os CSV de lop/resultados/ em tabelas Markdown (vírgula decimal) para os documentos do TCC, de modo que
nenhum número seja digitado à mão:

    python scripts/tabelas_batelada.py [--resultados DIR] [--saida ARQUIVO.md]

Grava resultados/tabelas.md. Só lê os CSV (e recalcula o White-Box, que é determinístico, para a Tabela do PBM).
"""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

DIR_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
if DIR_SCRIPTS not in sys.path:
    sys.path.insert(0, DIR_SCRIPTS)

from data_loader import carregar_dados_bancada          # noqa: E402
from hybrid_model import CONFIG_WB_PRINCIPAL, construir_base  # noqa: E402

DIR_RESULTADOS = os.path.abspath(os.path.join(DIR_SCRIPTS, '..', 'resultados'))
ROT = {'gp': 'GP', 'extra_trees': 'Extra Trees', 'gbdt': 'GBDT', 'random_forest': 'Random Forest', 'svr': 'SVR',
       'mlp_tanh': 'MLP tanh', 'ridge': 'Ridge'}


def num(v, casas=4):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return '–'
    return f'{v:.{casas}f}'.replace('.', ',')


def md(df: pd.DataFrame, alinhar=None) -> str:
    cols = list(df.columns)
    alinhar = alinhar or (['---'] + ['---:'] * (len(cols) - 1))
    linhas = ['| ' + ' | '.join(str(c) for c in cols) + ' |', '| ' + ' | '.join(a.replace('---', ':---', 1) if a == '---' else a for a in alinhar) + ' |']
    for _, r in df.iterrows():
        linhas.append('| ' + ' | '.join(str(r[c]) for c in cols) + ' |')
    return '\n'.join(linhas)


def nome_curto(m: str) -> str:
    h = re.match(r'Híbrido (\w+) \[(\w)\]', m)
    if h:
        return f'Híbrido {ROT[h.group(1)]} [{h.group(2)}]'
    c = re.match(r'Caixa-preta (\w+)', m)
    if c:
        return f'Caixa-preta {ROT[c.group(1)]}'
    return m


def ler(R, n):
    return pd.read_csv(os.path.join(R, n))


def t_piso_ruido(R):
    p = ler(R, 'piso_ruido.csv')
    p['valor'] = [num(v, 4) if g.startswith('sigma') else str(int(v)) for g, v in zip(p['grandeza'], p['valor'])]
    return '### T01. Piso de ruído experimental\n\n' + md(p.rename(columns={'grandeza': 'Grandeza', 'valor': 'Valor'}))


def t_escada_wb(R):
    e = ler(R, 'escada_whitebox.csv')
    out = pd.DataFrame({'Degrau': e['degrau'], 'RMSE (t>0)': e['RMSE_t>0'].map(num), 'R² (t>0)': e['R2_t>0'].map(num),
                        'R² (todos)': e['R2_todos'].map(num), 'RMSE pH₀': e['RMSE_pH0'].map(lambda v: num(v, 3)),
                        'RMSE pH_f (espec.)': e['RMSE_pHf_espec'].map(lambda v: num(v, 3)),
                        'RMSE pH_f (ponta a ponta)': e['RMSE_pHf_ponta'].map(lambda v: num(v, 3))})
    return '### T02. Escada do White-Box (X em 112 pontos; pH em 16 e 12 ensaios)\n\n' + md(out)


def t_selecao_wb(R):
    s = ler(R, 'selecao_whitebox.csv')
    out = pd.DataFrame({'Especiação': s['especiacao'], 'RMSE pH₀': s['RMSE_pH0'].map(lambda v: num(v, 3)),
                        'RMSE pH_f': s['RMSE_pHf_espec'].map(lambda v: num(v, 3)),
                        'viés pH_f': s['vies_pHf_espec'].map(lambda v: num(v, 3)),
                        'RMSE pH combinado': s['RMSE_pH_combinado_espec'].map(lambda v: num(v, 3)),
                        'RMSE X (t>0)': s['RMSE_t>0'].map(num), 'R² X (t>0)': s['R2_t>0'].map(num)})
    return '### T03. Seleção do White-Box por erro combinado de pH (ν a priori)\n\n' + md(out)


def t_ph_por_ensaio(R):
    t = ler(R, 'termodinamica_por_ensaio.csv')
    p = t[t['especiacao'] == 'V2'][['ensaio_id', 'eta_nominal', 'ca0', 'pH0_med', 'pHf_med', 'pHf_comparavel']].copy()
    for esp in ('V0', 'V1', 'V2'):
        g = t[t['especiacao'] == esp].set_index('ensaio_id')
        p[f'pH₀ prev. {esp}'] = p['ensaio_id'].map(g['pH0_prev']).map(lambda v: num(v, 2))
        p[f'pH_f prev. {esp}'] = p.apply(lambda r: num(g.loc[r['ensaio_id'], 'pHf_espec'], 2) if r['pHf_comparavel'] else '–', axis=1)
    p['pH₀ med.'] = p['pH0_med'].map(lambda v: num(v, 1))
    p['pH_f med.'] = p['pHf_med'].map(lambda v: num(v, 1))
    out = p[['ensaio_id', 'eta_nominal', 'ca0', 'pH₀ med.', 'pH₀ prev. V0', 'pH₀ prev. V2', 'pH_f med.', 'pH_f prev. V0',
             'pH_f prev. V1', 'pH_f prev. V2']].rename(columns={'ensaio_id': 'Ensaio', 'eta_nominal': 'η', 'ca0': 'C_A0 (mol/L)'})
    out['η'] = out['η'].map(lambda v: num(v, 1))
    out['C_A0 (mol/L)'] = out['C_A0 (mol/L)'].map(lambda v: num(v, 1))
    return '### T04. pH medido e previsto (especiação sobre o estado final medido; ν a priori)\n\n' + md(out)


def t_benchmark(R):
    b = ler(R, 'benchmark_batelada.csv')
    blocos = []
    for esq, g in b.groupby('esquema', sort=False):
        g = g.sort_values('RMSE')
        out = pd.DataFrame({
            'Modelo': g['modelo'].map(nome_curto), 'Tipo': g['tipo'],
            'RMSE': g['RMSE'].map(num),
            'IC95%': [f'[{num(a, 4)}; {num(c, 4)}]' for a, c in zip(g['RMSE_ic95_inf'], g['RMSE_ic95_sup'])],
            'dp entre partições': g['RMSE_dp_reps'].map(lambda v: num(v, 4)),
            'R² (t>0)': g['R2'].map(num), 'Ganho vs WB': g['ganho'].map(lambda v: num(v, 3)),
            '% do teto': (100.0 * g['ganho'] / g['teto_ganho']).map(lambda v: num(v, 0)),
            'bruto > teto (%)': g['pct_acima_teto'].map(lambda v: num(v, 1)),
            'bruto quedas (%)': g['pct_quedas'].map(lambda v: num(v, 1))})
        blocos.append(f'**{esq}** (teto do ganho = {num(g["teto_ganho"].iloc[0], 3)})\n\n' + md(out))
    return '### T05. Benchmark por protocolo de validação (RMSE em t>0; ganho = 1 − MSE/MSE_WB)\n\n' + '\n\n'.join(blocos)


def t_por_nivel_eta(R):
    b = ler(R, 'benchmark_batelada.csv')
    g = b[b['esquema'].str.startswith('Deixa-um-η')].sort_values('RMSE')
    out = pd.DataFrame({'Modelo': g['modelo'].map(nome_curto), 'RMSE total': g['RMSE'].map(num),
                        'η = 0,5 ausente': g['RMSE_eta0.5'].map(num), 'η = 1,0 ausente': g['RMSE_eta1'].map(num),
                        'η = 1,5 ausente': g['RMSE_eta1.5'].map(num), 'η = 3,1 ausente': g['RMSE_eta3.1'].map(num)})
    return '### T06. Extrapolação em η: RMSE no nível deixado fora do treino\n\n' + md(out)


def t_pareadas(R):
    c = ler(R, 'comparacoes_pareadas.csv')
    sel = c[c['A'].str.startswith('Híbrido') & c['B'].str.startswith('Caixa-preta')]
    sel = sel[sel.apply(lambda r: r['A'].split()[1] == r['B'].split()[1], axis=1)]
    out = pd.DataFrame({'Protocolo': sel['esquema'], 'Híbrido': sel['A'].map(nome_curto), 'Caixa-preta': sel['B'].map(nome_curto),
                        'Δ RMSE (híbrido − caixa-preta)': sel['delta_RMSE'].map(lambda v: num(v, 4)),
                        'IC95%': [f'[{num(a, 4)}; {num(b, 4)}]' for a, b in zip(sel['ic95_inf'], sel['ic95_sup'])],
                        'Empate (IC contém 0)': sel['empate_ic'].map(lambda v: 'sim' if v else 'não')})
    return '### T07. Híbrido contra a caixa-preta do MESMO modelo (diferença pareada, bootstrap de ensaios)\n\n' + md(out)


def t_pareadas_melhor(R):
    c = ler(R, 'comparacoes_pareadas.csv')
    sel = c[c['A'] == c['melhor_hibrido']]
    sel = sel[sel['B'].isin(['White-Box', 'Tabela (η, t)', 'Ridge (resíduo)']) | sel['B'].str.startswith('Caixa-preta') | sel['B'].str.startswith('Híbrido')]
    out = pd.DataFrame({'Protocolo': sel['esquema'], 'Melhor híbrido': sel['A'].map(nome_curto), 'Contra': sel['B'].map(nome_curto),
                        'Δ RMSE': sel['delta_RMSE'].map(lambda v: num(v, 4)),
                        'IC95%': [f'[{num(a, 4)}; {num(b, 4)}]' for a, b in zip(sel['ic95_inf'], sel['ic95_sup'])],
                        'Empate': sel['empate_ic'].map(lambda v: 'sim' if v else 'não')})
    return '### T08. Melhor híbrido de cada protocolo contra os demais modelos (Δ < 0: o melhor híbrido erra menos)\n\n' + md(out)


def t_aprendizado(R):
    c = ler(R, 'curva_aprendizado.csv')
    p = c.pivot(index='modelo', columns='n_treino', values='RMSE')
    p = p.loc[[m for m in sorted(p.index, key=lambda m: ('White' not in m, 'Tabela' not in m, 'Híbrido' not in m, m))]]
    out = p.map(num).reset_index()
    out['modelo'] = out['modelo'].map(nome_curto)
    out.columns = ['Modelo'] + [f'n = {n}' for n in p.columns]
    return '### T09. Curva de aprendizado: RMSE em 2 ensaios fora do treino, por nº de ensaios de treino (30 sorteios)\n\n' + md(out)


def t_ablacao(R):
    a = ler(R, 'ablacao_entradas.csv')
    a['slug'] = a['modelo'].str.extract(r'Híbrido (\w+) \[')[0].map(ROT)
    p = a.pivot_table(index=['esquema_ablacao', 'slug'], columns='entradas', values='RMSE', sort=False).reset_index()
    out = pd.DataFrame({'Protocolo': p['esquema_ablacao'], 'Modelo': p['slug'], 'A (t, η, C_A0)': p['A'].map(num),
                        'B (A + X_PBM)': p['B'].map(num), 'C (B + C_Af, [H⁺])': p['C'].map(num)})
    return '### T10. Ablação dos atributos do ML (RMSE em t>0)\n\n' + md(out)


def t_nu(R):
    s = ler(R, 'sensibilidade_nu.csv')
    s = s[s['modelo'].isin(['White-Box', 'Híbrido gp [B]', 'Híbrido extra_trees [B]'])]
    p = s.pivot_table(index=['esquema', 'modelo'], columns='nu_modo', values='RMSE', sort=False).reset_index()
    cols = [c for c in p.columns if c not in ('esquema', 'modelo')]
    out = pd.DataFrame({'Protocolo': p['esquema'], 'Modelo': p['modelo'].map(nome_curto)})
    for c in cols:
        out[c] = p[c].map(num)
    return '### T11. Sensibilidade ao fator ν (RMSE em t>0)\n\n' + md(out)


def t_escada_hibrido(R):
    e = ler(R, 'escada_hibrido.csv')
    p = e.pivot_table(index=['esquema', 'modelo'], columns='degrau_whitebox', values='RMSE', sort=False).reset_index()
    cols = [c for c in p.columns if c not in ('esquema', 'modelo')]
    out = pd.DataFrame({'Protocolo': p['esquema'], 'Modelo': p['modelo'].map(nome_curto)})
    for c in cols:
        out[c] = p[c].map(num)
    return '### T12. O mesmo ML sobre degraus do White-Box (RMSE em t>0)\n\n' + md(out)


def t_cobertura(R):
    c = ler(R, 'cobertura_intervalos.csv')
    out = pd.DataFrame({'Protocolo': c['esquema'], 'Modelo': c['modelo'].map(nome_curto), 'Método': c['metodo'],
                        'Cobertura': c['cobertura'].map(lambda v: num(100 * v, 0) + '%'),
                        'Largura média': c['largura_media'].map(num), 'n': c['n']})
    return '### T13. Cobertura empírica dos intervalos de 95% (t>0, pontos fora-da-dobra)\n\n' + md(out)


def t_residuo_patamar(R):
    p = ler(R, 'residuo_patamar_por_ensaio.csv')
    out = pd.DataFrame({'Ensaio': p['ensaio_id'], 'η': p['razao_molar_eta'].map(lambda v: num(v, 1)),
                        'η real': p['razao_molar_eta_real'].map(lambda v: num(v, 3)),
                        'C_A0': p['C_acid_0_mol_L'].map(lambda v: num(v, 1)), 'X_exp (15 min)': p['X_zn_exp'].map(lambda v: num(v, 2)),
                        'X_PBM (15 min)': p['X_pbm'].map(num), 'ΔX': p['delta_X'].map(lambda v: num(v, 3))})
    return '### T14. Resíduo do White-Box principal aos 15 min\n\n' + md(out)


def t_escalas_gp(R):
    g = ler(R, 'gp_escalas.csv')
    out = pd.DataFrame({'Atributos': g['entradas'], 'Entrada': g['atributo'], 'Escala de comprimento (σ das entradas)': g['escala_de_comprimento'].map(lambda v: num(v, 3)),
                        'Ruído branco (norm.)': g['ruido_branco'].map(lambda v: num(v, 3))})
    return '### T15. Escalas de comprimento do GP (ARD) ajustadas com todos os dados\n\n' + md(out)


def t_pbm(df):
    base = construir_base(df, CONFIG_WB_PRINCIPAL)
    tempos = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0]
    blocos = []
    for nome, col in (('X_PBM (conversão do White-Box principal)', 'X_pbm'), ('ΔX = X_exp − X_PBM', 'delta_X')):
        p = base.pivot(index='ensaio_id', columns='tempo_min', values=col)
        meta = base.groupby('ensaio_id')[['razao_molar_eta', 'C_acid_0_mol_L']].first()
        out = pd.DataFrame({'Ensaio': p.index, 'η': meta['razao_molar_eta'].map(lambda v: num(v, 1)).values,
                            'C_A0': meta['C_acid_0_mol_L'].map(lambda v: num(v, 1)).values})
        for t in tempos:
            out[f'{t:g} min'.replace('.', ',')] = p[t].map(lambda v: num(v, 2)).values
        blocos.append(f'**{nome}**\n\n' + md(out))
    return '### T16. Conversão do White-Box principal e resíduo, por ensaio\n\n' + '\n\n'.join(blocos)


def main(argv=None):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--resultados', default=DIR_RESULTADOS)
    ap.add_argument('--saida', default=None)
    a = ap.parse_args(argv)
    R = a.resultados
    secoes = [t_piso_ruido(R), t_escada_wb(R), t_selecao_wb(R), t_ph_por_ensaio(R), t_benchmark(R), t_por_nivel_eta(R),
              t_pareadas(R), t_pareadas_melhor(R), t_aprendizado(R), t_ablacao(R), t_nu(R), t_escada_hibrido(R),
              t_cobertura(R), t_residuo_patamar(R), t_escalas_gp(R), t_pbm(carregar_dados_bancada())]
    cab = ('# Tabelas dos resultados da batelada\n\nGerado por `scripts/tabelas_batelada.py` a partir de `resultados/*.csv`. '
           'Não edite à mão.\n\n')
    saida = a.saida or os.path.join(R, 'tabelas.md')
    with open(saida, 'w', encoding='utf-8') as f:
        f.write(cab + '\n\n'.join(secoes) + '\n')
    print(f'-> {saida}')


if __name__ == '__main__':
    main()
