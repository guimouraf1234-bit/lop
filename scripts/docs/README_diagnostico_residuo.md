# README -- Módulo de Diagnóstico do Resíduo (`diagnostico_residuo.py`)

## 1. Visão Geral e Objetivo
O módulo `diagnostico_residuo.py` realiza a ponte quantitativa entre o **Modelo Mecanicista (\textit{White-Box})** e o **Modelo de Aprendizado de Máquina (\textit{Black-Box})**.

Ele executa a simulação do Balanço Populacional (`population_balance.py`) para todos os 16 ensaios de bancada da tese de Bortot Coelho (2017), compara com os 128 pontos experimentais reais e calcula o **Resíduo Físico ($\Delta X$)**:
$$\Delta X = X_{\text{experimental}} - X_{\text{PBM}}$$

Ao final, gera e salva o arquivo de dados consolidado `data/processed/dados_com_residuo.csv`, que serve como conjunto de dados de entrada para o treinamento do modelo híbrido (*Grey-Box*).

---

## 2. Fundamentação Teórica: O Significado do Resíduo

O modelo de primeiros princípios assume partículas esféricas homogêneas e ácido sulfúrico em solução ideal. Na prática experimental de laboratório, ocorrem três fenômenos secundários complexos:

1. **Desvios Iniciais Positivos ($\Delta X > 0$ nos primeiros minutos):**
   * O concentrado ustulado contém uma fração elevada de partículas ultrafinas com arestas vivas e alta rugosidade superficial.
   * Essas partículas dissolvem-se quase que instantaneamente no contato inicial com o ácido, fazendo com que a conversão real nos primeiros 60 segundos seja ligeiramente superior à de esferas lisas ideais.

2. **Desvios Terminais Negativos ($\Delta X < 0$ nos minutos finais):**
   * Conforme a lixiviação avança, acumulam-se grandes quantidades de íons $\text{Zn}^{2+}$, $\text{Fe}^{2+}/\text{Fe}^{3+}$ e $\text{SO}_4^{2-}$ em solução, reduzindo a atividade termodinâmica do próton $\text{H}^+$ (**efeito de íon comum**).
   * Silicatos solúveis presentes na calcina polimerizam em meio ácido, formando um filme superficial de gel de sílica amorfa (**passivação difusional**).
   * A reação desacelera antes do previsto pela física pura.

O resíduo $\Delta X$ não é ruído aleatório gaussiano: **ele carrega uma assinatura físico-química determinística**, perfeitamente adequada para ser aprendida por algoritmos de regressão supervisionada.

---

## 3. Desempenho do Modelo Mecanicista Puro (White-Box)

A avaliação estatística dos 128 pontos de bancada confirma a alta precisão da física fundamental:

* **Coeficiente de Determinação ($R^2$):** **$0{,}9253$**  
  *(A física pura sozinha explica $92{,}5\%$ da variância de todos os experimentos)*.
* **Raiz do Erro Quadrático Médio ($RMSE$):** **$0{,}0886$** ($8{,}86\%$).
* **Erro Médio Absoluto ($MAE$):** **$0{,}0507$** ($5{,}07\%$).

---

## 4. Estrutura do Código e Dados Gerados

### `calcular_residuos() -> pd.DataFrame`
1. Chama `carregar_dados_bancada()` para obter os dados experimentais.
2. Agrupa por `ensaio_id` e executa `pbm.simular_ensaio` para cada condição operacional $(C_{A0}, \eta, t)$.
3. Adiciona as colunas:
   * `X_pbm`: Conversão teórica mecanicista calculada.
   * `delta_X`: Resíduo ($X_{\text{exp}} - X_{\text{pbm}}$).
4. Exporta para `data/processed/dados_com_residuo.csv`.

---

## 5. Exemplo de Saída dos Dados Gerados

Trecho dos primeiros registros com o resíduo calculado para o Ensaio 1 ($C_{A0} = 0{,}1\text{ mol/L}$, $\eta = 0{,}5$):

| ensaio_id | tempo_min | X_zn_exp | X_pbm | delta_X |
| :---: | :---: | :---: | :---: | :---: |
| 1 | 0.0 | 0.00 | 0.0000 | 0.0000 |
| 1 | 0.5 | 0.42 | 0.3239 | +0.0961 |
| 1 | 1.0 | 0.46 | 0.4216 | +0.0384 |
| 1 | 2.0 | 0.46 | 0.4810 | -0.0210 |
| 1 | 3.0 | 0.49 | 0.4950 | -0.0050 |
| 1 | 4.0 | 0.51 | 0.4987 | +0.0113 |
| 1 | 5.0 | 0.49 | 0.4996 | -0.0096 |
| 1 | 15.0 | 0.50 | 0.5000 | 0.0000 |

---

## 6. Como Executar
Execute no terminal:
```powershell
& .venv\Scripts\python.exe scripts/diagnostico_residuo.py
```
*Saída:*
```text
Resíduos calculados e salvos com sucesso em: ...\data\processed\dados_com_residuo.csv
```
O arquivo `data/processed/dados_com_residuo.csv` estará pronto para alimentar o **Modelo Híbrido Grey-Box**!
