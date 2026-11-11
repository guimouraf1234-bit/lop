# CLAUDE.md — TCC: modelagem híbrida da lixiviação de ustulado de zinco

Projeto de TCC em grupo (Engenharia Química, UFMG). Este repositório é a abordagem **híbrida paralela cooperativa**:
White-Box (balanço populacional, `population_balance.py` + `kinetics.py`) + regressor de ML que aprende o resíduo
ΔX = X_exp − X_pbm + acoplador físico (teto estequiométrico, X(0) = 0, regressão isotônica). Outros membros do grupo
fizeram o híbrido **serial** (rede prevê a taxa de retração), um black-box e uma PINN — os resultados deles estão no
relatório Word, não neste repositório. Textos, nomes e comentários em português.

## Ambiente
```bash
pip install -r requirements.txt
python -m pytest -q          # ~35 s; RUN_SLOW=1 inclui o teste de determinismo
```
Rodar os comandos a partir de `lop/` (raiz do repositório). Arquivos `.py` usam fim de linha CRLF.

## Pipeline (todo número, tabela e figura sai daqui)
```bash
python scripts/experimentos_batelada.py                      # todas as etapas -> resultados/*.csv (~15 min)
python scripts/experimentos_batelada.py --etapas particao_comum --saida /tmp/res   # uma etapa
python scripts/figuras_batelada.py [--figuras f13,f13r,f14]  # lê os CSV -> resultados/figuras/
python scripts/tabelas_batelada.py                           # lê os CSV -> resultados/tabelas.md
```
**Atenção:** o pipeline sempre regrava `config.json` na pasta de saída. Ao rodar só algumas etapas, use `--saida`
numa pasta temporária e copie apenas os CSV daquela etapa para `resultados/` — senão o `config.json` fica parcial.

## Partição comum com os outros modelos do grupo (etapa `particao_comum`)
Para comparar com o híbrido serial e o black-box do relatório, usa-se a mesma partição (Apêndice B do relatório):
- **teste cego:** Ensaios 7, 8 e 14 (21 pontos com t > 0), consultado uma única vez;
- **treino:** os outros 13 ensaios (91 pontos); **CV** com as 4 dobras fixas da Tabela B.3:
  {1, 5, 11, 16}, {4, 10, 15}, {3, 9, 13}, {2, 6, 12} (`validacao.DOBRAS_TABELA_B3`);
- 10 sementes (42–51); modelo final = média das 10 sementes; escolha do regressor **só pela CV**;
- métricas em p.p. de conversão, só em t > 0: R², RMSE, MAE, erro máximo; inconsistência física na malha de 0,25 min
  antes e depois do acoplador.

Regressores testados como ML do híbrido e como caixa-preta: os 6 do catálogo principal
(`hybrid_model.obter_especificacoes`) + 6 adicionais (`obter_especificacoes_adicionais`: k-NN, HistGBDT, AdaBoost,
Kernel Ridge, Ridge polinomial, GP RBF). Para testar outro regressor: criar a fábrica e o `ModeloSpec` em
`obter_especificacoes_adicionais()`, incluir o slug em `SLUGS_PARTICAO_COMUM` (`experimentos_batelada.py`) e em
`MARC`/`ROT` (`figuras_batelada.py`).

Saídas: `resultados/particao_comum_*.csv`, tabela T17 em `resultados/tabelas.md`, figuras
`resultados/figuras/f13_teste_cego_particao_comum.png`, `f14_resumo_regressores_particao_comum.png` e uma figura
por regressor em `resultados/figuras/particao_comum/` (essas são versionadas; as demais figuras são ignoradas pelo git).

## Regras
- Nunca usar os Ensaios 7, 8 e 14 para escolher modelo, hiperparâmetro ou critério (teste em `tests/test_validacao.py`).
- Não alterar o benchmark principal (`obter_especificacoes`, `codigos_benchmark`) sem rodar o pipeline completo.
- Branch de trabalho: `revisao/validade-hibrido`.
