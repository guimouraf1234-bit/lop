# README -- Pipeline único (`experimentos_batelada.py`, `figuras_batelada.py`, `tabelas_batelada.py`)

Todo número, tabela e figura do trabalho sai destes três comandos (a partir da raiz de `lop/`):

```powershell
& ..\.venv\Scripts\python.exe scripts/experimentos_batelada.py       # ~11 min (12 núcleos); protocolo completo -> resultados/*.csv, config.json
& ..\.venv\Scripts\python.exe scripts/figuras_batelada.py            # lê os CSV -> resultados/figuras/*.png (12 figuras)
& ..\.venv\Scripts\python.exe scripts/tabelas_batelada.py            # lê os CSV -> resultados/tabelas.md (T01 a T16)
& ..\.venv\Scripts\python.exe scripts/experimentos_batelada.py --rapido   # versão reduzida (~3 min), para desenvolvimento
& ..\.venv\Scripts\python.exe -m pytest                                # 56 testes (RUN_SLOW=1 inclui o teste de determinismo do benchmark)
```
Opções: `--etapas whitebox,benchmark` (subconjunto), `--saida DIR`, `--n-jobs N`.

## Etapas e arquivos
| Etapa | Arquivos em `resultados/` | O que responde |
| :--- | :--- | :--- |
| `dados` | `piso_ruido.csv` | σ do erro experimental |
| `whitebox` | `escada_whitebox.csv`, `selecao_whitebox.csv`, `termodinamica_por_ensaio.csv` | qual White-Box é coerente com o pH medido |
| `ablacao` | `ablacao_entradas.csv` | quais atributos o ML precisa (A, B ou C) |
| `benchmark` | `benchmark_batelada.csv`, `comparacoes_pareadas.csv`, `predicoes_oof.csv.gz` | híbrido × baselines × caixa-preta em 6 protocolos |
| `nu` | `sensibilidade_nu.csv` | efeito da proveniência de ν |
| `aprendizado` | `curva_aprendizado.csv` | valor da física com poucos dados |
| `intervalos` | `cobertura_intervalos.csv`, `cobertura_intervalos_pontos.csv.gz` | calibração dos intervalos de 95% |
| `residuo` | `residuo_aprendido.csv`, `gp_escalas.csv`, `residuo_patamar_por_ensaio.csv` | o que o ML aprendeu |
| `dominio` | `dominio_incerteza_eta.csv` | discordância e incerteza entre os níveis de η |
| `escada_hibrido` | `escada_hibrido.csv` | o híbrido sobre degraus do White-Box |
| `curvas` | `curvas_logo.csv` | curvas contínuas fora-da-dobra (figura F01) |

`config.json` registra versões, sementes, o White-Box e os atributos principais, o σ e se as regras de seleção (`whitebox`, `entradas`)
concordam com as constantes `CONFIG_WB_PRINCIPAL` e `ENTRADAS_PRINCIPAL` de `hybrid_model.py` (um teste falha se divergirem).

## Reprodutibilidade
Sementes fixas (modelos 42; partições 0…19; curva de aprendizado 1000…; bootstrap 0), floats com 8 algarismos significativos, ordem estável das
dobras; duas execuções geram CSV idênticos (`tests/test_pipeline.py`). Sem cache em disco: o PBM é recalculado a cada execução.

## Figuras (`figuras_batelada.py`)
Cor segue a **entidade** (híbrido azul, caixa-preta laranja, White-Box tinta, baselines cinza); o modelo é dado pelo marcador (identidade nunca só por cor).
Rampa azul ordinal para A < B < C e $C_{A0}$. Títulos quantitativos são montados dos dados. Cada figura tem seu gêmeo em tabela nos CSV.

| Figura | Conteúdo |
| :--- | :--- |
| f01 | cinética dos 16 ensaios, curvas **fora-da-dobra** (LOGO) com faixa de 95% do GP |
| f02 | RMSE por protocolo e modelo, com IC95% e as linhas do White-Box e de σ |
| f03 | curva de aprendizado |
| f04, f05 | escada do White-Box e pH final previsto × medido |
| f06 | ablação dos atributos |
| f07 | extrapolação por nível de η |
| f08 | cobertura dos intervalos |
| f09, f10 | resíduo do White-Box e aprendido; domínio de aplicabilidade |
| f11 | o ML sobre degraus do White-Box |
| f12 | paridade fora-da-dobra |
