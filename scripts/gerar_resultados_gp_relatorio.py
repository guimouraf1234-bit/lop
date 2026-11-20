"""Recalcula exclusivamente o GP adotado e gera figuras e texto para o relatório.

Executar na raiz: python scripts/gerar_resultados_gp_relatorio.py
Saídas independentes das análises antigas: entregas/resultados_gp/.
"""
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from joblib import Parallel, delayed, parallel_config
from data_loader import carregar_dados_bancada
from hybrid_model import CONFIG_WB_PRINCIPAL as CFG, construir_base, construir_base_condicao, AcopladorFisico
import validacao as V
from experimentos_batelada import _metricas_pp, _pct_inconsistente

OUT = Path(__file__).resolve().parents[1] / 'entregas' / 'resultados_gp'
SEEDS = tuple(range(42, 52))
TESTE = V.ENSAIOS_TESTE_CEGO
ACOP = AcopladorFisico(nu=CFG.nu)


def ajustar(seed, tr, malha):
    warnings.filterwarnings('ignore')
    m = V.criar_estimador(('hibrido', 'gp', 'B'), seed).ajustar(tr)
    return m.prever_bruto(malha)


def pontos_da_malha(pontos, malha, valores):
    indice = pd.MultiIndex.from_frame(malha[['ensaio_id', 'tempo_min']])
    pos = indice.get_indexer(pd.MultiIndex.from_frame(pontos[['ensaio_id', 'tempo_min']]))
    assert (pos >= 0).all()
    return np.asarray(valores)[pos]


def salvar_fig(fig, nome):
    fig.savefig(OUT / (nome + '.png'), dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT / (nome + '.pdf'), bbox_inches='tight', facecolor='white')
    plt.close(fig)


# Paleta da modelagem paralela no relatório (Figuras 5.4.1 e 5.4.2): white-box azul pontilhado, GP azul-claro
# tracejado, limite estequiométrico cinza pontilhado e pontos experimentais vazados.
COR_WB, COR_GP, COR_LIMITE = '#0072B2', '#56B4E9', '#666666'
TRACO_GP = (0, (6, 3))


def painel(base, curvas, ensaios, shape, nome):
    from matplotlib.ticker import FuncFormatter
    grande = len(ensaios) == 16
    fonte = 11 if grande else 14
    plt.rcParams.update({'font.family': 'Times New Roman', 'mathtext.fontset': 'stix',
                         'font.size': fonte, 'axes.linewidth': 1.2})
    fig, axes = plt.subplots(*shape, figsize=(17, 13) if grande else (12.8, 9.8),
                             sharex=True, sharey=True, squeeze=False)
    fmt = lambda v, n: f'{v:.{n}f}'.replace('.', ',')
    for k, (ax, e) in enumerate(zip(axes.flat, ensaios)):
        lin, col = divmod(k, shape[1])
        g = curvas[curvas.ensaio_id == e]
        todos = base[base.ensaio_id == e]
        p = todos[todos.tempo_min > 0]
        ax.plot(g.tempo_min, g.X_pbm, ls=':', color=COR_WB, lw=2.6)
        ax.plot(g.tempo_min, g.x_hat, ls=TRACO_GP, color=COR_GP, lw=3)
        ax.axhline(min(1, p.razao_molar_eta.iloc[0]), color=COR_LIMITE, ls=':', lw=1.3)
        ax.scatter(todos.tempo_min, todos.X_zn_exp, s=44 if grande else 62, facecolors='white',
                   edgecolors='black', linewidths=1.4, zorder=5)
        r = p.iloc[0]
        ax.text(.03, .97, f'Ensaio {e}{"*" if e in TESTE else ""}',
                transform=ax.transAxes, ha='left', va='top', fontsize=fonte+3, fontweight='bold')
        met = _metricas_pp(p.X_zn_exp, p.x_hat)
        wb = _metricas_pp(p.X_zn_exp, p.X_pbm)
        nota = (f'$C_{{A0}}$ = {fmt(r.C_acid_0_mol_L, 2)} mol/L, $\\eta$ = {fmt(r.razao_molar_eta, 2)}'
                + '\n' + f'RMSE: white-box {fmt(wb["RMSE_pp"] / 100, 3)}; GP {fmt(met["RMSE_pp"] / 100, 3)}')
        ax.text(.97, .05, nota, transform=ax.transAxes, ha='right', va='bottom', fontsize=fonte,
                linespacing=1.5)
        ax.set(xlim=(-.3, 15.5), ylim=(0, 1.12), xticks=[0, 5, 10, 15], yticks=np.arange(0, 1.01, .2))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: fmt(v, 1)))
        ax.tick_params(axis='both', direction='in', length=6, width=1.1, pad=5,
                       labelbottom=lin == shape[0] - 1, labelleft=col == 0)
        if lin == shape[0] - 1:
            ax.set_xlabel('Tempo de reação, $t$ (min)', fontsize=fonte+2)
        if col == 0:
            ax.set_ylabel('Conversão de zinco, $X_{Zn}$ (–)', fontsize=fonte+2)
        ax.grid(color='#bfbfbf', lw=.9, ls=':')
        ax.set_axisbelow(True)
    handles = [Line2D([], [], color=COR_WB, ls=':', lw=2.6, marker='^', markersize=8,
                      label='White-box (novo)'),
               Line2D([], [], color=COR_GP, ls=TRACO_GP, lw=3, marker='P', markersize=9,
                      label='Híbrido paralelo (GP)'),
               Line2D([], [], color=COR_LIMITE, ls=':', lw=1.3,
                      label='Limite estequiométrico, min(1, $\\eta$)'),
               Line2D([], [], color='black', marker='o', markerfacecolor='white',
                      markeredgewidth=1.4, markersize=8, ls='none', label='Experimental')]
    fig.legend(handles=handles, loc='upper center', ncol=4, bbox_to_anchor=(.53, .995),
               fontsize=fonte+1, frameon=False, handlelength=3.0, columnspacing=1.6)
    fig.subplots_adjust(left=.08, right=.985, bottom=.08, top=.92, wspace=.06, hspace=.06)
    salvar_fig(fig, nome)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings('ignore')
    b = construir_base(carregar_dados_bancada(), CFG).sort_values(['ensaio_id', 'tempo_min']).reset_index(drop=True)
    malha = pd.concat([construir_base_condicao(float(g.razao_molar_eta_real.iloc[0]),
                         float(g.C_acid_0_mol_L.iloc[0]), np.arange(0, 15.001, .25), CFG,
                         ensaio_id=int(e), eta_nominal=float(g.razao_molar_eta.iloc[0]))
                      for e, g in b.groupby('ensaio_id')], ignore_index=True)
    treino = V.ensaios_treino_comum(b)
    dobras = V.esquema_cv_tabela_b3(b)
    tarefas = [(seed, d) for seed in SEEDS for d in dobras]
    print('Validação cruzada: GP fixado, 4 dobras e 10 sementes.', flush=True)
    with parallel_config(backend='loky', n_jobs=4, inner_max_num_threads=1):
        bruto_cv = Parallel()(delayed(ajustar)(seed, b[b.ensaio_id.isin(d.treino)],
                                               malha[malha.ensaio_id.isin(d.teste)]) for seed, d in tarefas)
    cv = []
    for (seed, d), bruto in zip(tarefas, bruto_cv):
        m = malha[malha.ensaio_id.isin(d.teste)]
        p = b[(b.ensaio_id.isin(d.teste)) & (b.tempo_min > 0)].copy()
        p['x_hat'] = pontos_da_malha(p, m, ACOP.restringir(m, bruto))
        p['semente'], p['dobra'] = seed, d.rotulo
        cv.append(p)
    cv = pd.concat(cv, ignore_index=True)
    cv.to_csv(OUT / 'predicoes_validacao_cruzada.csv', index=False)
    cv_seed = cv.groupby('semente').apply(lambda p: _metricas_pp(p.X_zn_exp, p.x_hat)['RMSE_pp'], include_groups=False)
    cv_dobra = cv.groupby(['semente', 'dobra']).apply(lambda p: _metricas_pp(p.X_zn_exp, p.x_hat)['RMSE_pp'], include_groups=False)
    print('Treinamento final: 13 ensaios, média das 10 sementes e projeção na malha comum.', flush=True)
    with parallel_config(backend='loky', n_jobs=4, inner_max_num_threads=1):
        raws = np.array(Parallel()(delayed(ajustar)(s, b[b.ensaio_id.isin(treino)], malha) for s in SEEDS))
    curvas = malha.copy()
    curvas['x_bruto'] = raws.mean(axis=0)
    curvas['x_hat'] = ACOP.restringir(malha, curvas.x_bruto)
    b['x_hat'] = pontos_da_malha(b, curvas, curvas.x_hat)
    b['x_bruto'] = pontos_da_malha(b, curvas, curvas.x_bruto)
    b['conjunto'] = np.where(b.ensaio_id.isin(TESTE), 'Teste', 'Treino')
    b.to_csv(OUT / 'predicoes_gp.csv', index=False)
    curvas.to_csv(OUT / 'curvas_gp.csv', index=False)
    linhas = []
    for nome, ids in [('Treino', treino), ('Teste', TESTE), ('Global', range(1, 17))]:
        p = b[b.ensaio_id.isin(ids) & (b.tempo_min > 0)]
        c = curvas[curvas.ensaio_id.isin(ids)]
        k = (p.X_zn_exp - p.x_hat).abs().idxmax()
        linhas.append({'conjunto': nome, 'ensaios': len(ids), **_metricas_pp(p.X_zn_exp, p.x_hat),
                       'erro_max_ensaio': int(p.loc[k, 'ensaio_id']), 'erro_max_t': float(p.loc[k, 'tempo_min']),
                       'inconsistencia_bruta_pct': _pct_inconsistente(c, c.x_bruto, ACOP.teto(c)),
                       'inconsistencia_final_pct': _pct_inconsistente(c, c.x_hat, ACOP.teto(c))})
    metricas = pd.DataFrame(linhas).set_index('conjunto')
    metricas.to_csv(OUT / 'metricas_gp.csv')
    por_ensaio = pd.DataFrame([{'ensaio': int(e), **_metricas_pp(g.X_zn_exp, g.x_hat)}
                               for e, g in b[b.tempo_min > 0].groupby('ensaio_id')])
    por_ensaio.to_csv(OUT / 'metricas_por_ensaio.csv', index=False)
    pd.DataFrame({'RMSE_pp': cv_seed}).to_csv(OUT / 'cv_por_semente.csv')
    cv_dobra.rename('RMSE_pp').to_csv(OUT / 'cv_por_dobra.csv')
    config = {'whitebox': CFG.para_dict(), 'teste': list(TESTE), 'treino': list(treino), 'sementes': list(SEEDS),
              'regressor': 'GP Matérn 3/2 (adotado pelo autor)', 'passo_min': .25,
              'projecao': 'Uma projeção isotônica na malha comum; métricas extraídas dessa curva.',
              'CV_RMSE_pp_media': float(cv_seed.mean()), 'CV_RMSE_pp_dp_sementes': float(cv_seed.std())}
    (OUT / 'config.json').write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding='utf-8')
    plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 10, 'axes.spines.top': True,
                         'axes.spines.right': True, 'figure.facecolor': 'white'})
    painel(b, curvas, [1, 7, 8, 14], (2, 2), 'figura_1_treino_teste_gp')
    painel(b, curvas, list(range(1, 17)), (4, 4), 'figura_2_dezesseis_ensaios_gp')
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family': 'Times New Roman', 'mathtext.fontset': 'stix', 'font.size': 14})
    fig, ax = plt.subplots(figsize=(6, 5.5))
    for nome, cor, marker in [('Treino', COR_WB, '^'), ('Teste', COR_GP, 'o')]:
        p = b[(b.conjunto == nome) & (b.tempo_min > 0)]
        ax.scatter(p.X_zn_exp, p.x_hat, facecolors='white', edgecolors=cor, marker=marker, s=40, label=nome, linewidths=1.3)
    ax.plot([0, 1.05], [0, 1.05], ':', color='#666666', lw=1.5)
    ax.set(xlim=(0, 1.05), ylim=(0, 1.05), xlabel='Conversão medida (–)', ylabel='Conversão prevista (–)')
    ax.tick_params(direction='in', top=True, right=True, length=6)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v:.1f}'.replace('.', ',')))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v:.1f}'.replace('.', ',')))
    ax.set_aspect('equal'); ax.grid(ls=':', color='#bfbfbf', alpha=.8)
    ax.set_axisbelow(True); ax.legend(frameon=False)
    fig.tight_layout(); salvar_fig(fig, 'figura_3_paridade_gp')
    print(metricas.to_string(), flush=True)
    print(json.dumps(config, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
