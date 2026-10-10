# README -- Validação (`validacao.py`)

## 1. Princípios
* **Por ensaio inteiro.** Nenhum ponto de um ensaio de teste entra no treino. O PBM é **recalculado em cada dobra**: o único parâmetro do White-Box que depende de dados, ν, é estimado só com ensaios de treino (`estimar_nu`).
* **Métricas só em $t > 0$.** Em $t = 0$ o resíduo é trivial e infla o R².
* **Piso de ruído** $\sigma$ (`piso_ruido`): desvio-padrão combinado nos patamares das curvas (médias de triplicatas): **0,0175**
  (von Neumann, que ignora deriva lenta: 0,0156), 64 graus de liberdade. Como o RMSE do White-Box é $e_{WB}$, o ganho máximo possível sobre ele é $1 - \sigma^2/e_{WB}^2$ (`teto_ganho`).
* **Baselines obrigatórios:** White-Box sozinho, tabela média por (η, t), Ridge e **caixa-preta** (o mesmo regressor prevendo X direto de t, η real, $C_{A0}$, com as mesmas restrições físicas). O híbrido só é "válido" onde supera todos.
* **Incerteza:** IC95% por **bootstrap de ensaios** do RMSE e das diferenças pareadas; `empate` = IC da diferença contém 0.

## 2. Protocolos (`esquema_*`)
| Protocolo | Dobras | Mede |
| :--- | :--- | :--- |
| `esquema_gkf8_sklearn` | 8 dobras de 2 ensaios (a partição do TCC) | interpolação (a mesma dos números antigos) |
| `esquema_gkf8_repetido` | 20 partições aleatórias em 8 dobras | interpolação, com dispersão entre partições |
| `esquema_logo` | 16 dobras (15 treino, 1 teste) | interpolação, uso máximo dos dados |
| `esquema_deixa_nivel(df, 'razao_molar_eta')` | 4 dobras: um nível de η inteiro fora | **extrapolação em η** |
| `esquema_deixa_nivel(df, 'C_acid_0_mol_L')` | 4 dobras | extrapolação em $C_{A0}$ |
| `esquema_extrapolacao_temporal` | treina t ≤ 5 min, prevê 15 min | extrapolação temporal (usa ν a priori) |
| `dobras_curva_aprendizado` | n ∈ {4,…,14} ensaios de treino, 2 de teste, 30 sorteios | valor da física com poucos dados |

## 3. Estimadores (`criar_estimador(codigo)`)
`('wb',)`, `('tabela',)`, `('ridge', entradas)`, `('hibrido', slug, entradas)`, `('caixa', slug)`. As classes `WhiteBoxPuro`, `TabelaEtaT`, `ResidualML`, `CaixaPretaML`
têm `ajustar(base_treino)` e `prever_bruto(base_teste)`; o acoplador é aplicado fora (igual para todos).

## 4. Execução e métricas
```python
import validacao as V
pred, diag = V.executar_esquema(df, resp, 'LOGO', V.esquema_logo(df), [('wb',), ('hibrido', 'gp', 'B'), ('caixa', 'gp')],
                                cfg_base=CONFIG_WB_PRINCIPAL, nu_modo='dobra', n_jobs=-1)
resumo = V.resumo_benchmark(pred, diag, sigma=V.piso_ruido(df)['sigma'])
pares  = V.comparacoes_pareadas(pred, [('Híbrido gp [B]', 'Caixa-preta gp')])
```
`resumo_benchmark` devolve por (protocolo, modelo): RMSE (média e dp entre repetições, IC95%), MAE, viés, R², **ganho** $= 1 - \mathrm{MSE}/\mathrm{MSE}_{WB}$,
teto do ganho, RMSE por nível de η e a **auditoria da saída bruta** (violações antes do acoplador).

## 5. Validação termodinâmica do White-Box
`tabela_termodinamica(resp, config)`: pH₀ previsto (16 ensaios) e pH final de duas formas — por **especiação sobre o estado final medido** (testa só a química)
e **ponta a ponta** (PBM em 15 min). `escada_whitebox` percorre os degraus da tese até V2 + ν; `selecionar_whitebox` aplica a regra: menor erro combinado de pH
(pH₀ e pH_f dos 12 ensaios com η ≥ 1) com ν a priori; o ajuste de X só desempata.

## 6. Intervalos
`cobertura_intervalos`: GP (desvio analítico) e **conformes** (quantil do erro fora-do-ensaio no treino, LOGO interno), avaliados em pontos fora-da-dobra (`resumo_cobertura`).

## 7. Testes
`tests/test_validacao.py` inclui o **teste de vazamento**: trocar `X_zn_exp` e as respostas finais dos ensaios de teste não pode alterar nenhuma previsão
(verificado com um mutante que treina com todos os ensaios: o teste o detecta).
