# Experimento de Data Augmentation Temporal (PCHIP + PBM)

Este diretório contém a implementação, os experimentos e os resultados gráficos da ampliação de dados temporal (*Temporal Data Augmentation*) aplicada tanto à Rede Neural **Multilayer Perceptron (MLP Tanh + L-BFGS)** quanto ao **Support Vector Regression (SVR Otimizado)**.

---

## 1. Motivação Científica

Na modelagem original com $N=128$ pontos (16 ensaios $\times$ 8 amostras de tempo):
* A validação cruzada por ensaios inteiros (*GroupKFold* com 8 dobras / proporção 85-15) disponibilizava apenas **112 pontos de treino** por rodada.
* Avaliou-se a hipótese de que o adensamento temporal de 128 para 640 pontos através de interpolação física guiada pelo modelo PBM poderia aprimorar o aprendizado das superfícies contínuas.

---

## 2. Metodologia de Augmentation Físico-Termodinâmica

Para expandir a base de dados de **128 para 640 pontos** sem introduzir ruídos ou violações físicas:
1. **Preservação de Monotonicidade ($dX/dt \ge 0$):**
   * Aplicação prévia de acúmulo monotônico $\max(X)$ para eliminar ruídos analíticos de titulação.
2. **Respeito ao Limite Estequiométrico ($\text{teto} = \min(1{,}0; \eta)$):**
   * Nenhum ponto sintético ultrapassa a conversão máxima permitida pelo reagente limitante.
3. **Interpolação PCHIP (*Piecewise Cubic Hermite Interpolating Polynomial*):**
   * Interpolação cúbica com derivadas contínuas que preserva monotonicidade estrita sem oscilações espúrias (efeito Runge).
4. **Solução Mecanicista PBM em Todos os Instantes Sintéticos:**
   * O Balanço Populacional analítico é integrado na grade densa de 40 tempos por ensaio, gerando os valores exatos de $X_{\text{PBM}}$ e $C_{A,\text{teórico}}$.
5. **Rigor Científico na Validação Cruzada (Blind Test):**
   * O adensamento de pontos foi aplicado **exclusivamente nas partições de treino**.
   * Os 2 ensaios deixados de fora em cada uma das 8 dobras foram avaliados **estritamente em seus pontos experimentais reais brutos** ($N=16$ por dobra), garantindo ausência total de vazamento de dados (*sem data leakage*).

---

## 3. Resultados Quantitativos Consolidados

Avaliados rigorosamente nos **pontos experimentais reais**:

| Modelo | Amostras no Treino | $R^2$ Ajuste | MAE Ajuste (%) | $R^2_{\text{CV}}$ (85/15) | MAE$_{\text{CV}}$ (%) | RMSE$_{\text{CV}}$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **MLP Tanh (Original - 128 pts)** | 112 (CV) / 128 (Full) | 0,9855 | 2,27% | 0,9702 | 2,80% | 0,0559 |
| **MLP Tanh + Augmentation (640 pts)** | 560 (CV) / 640 (Full) | 0,9981 | 1,07% | 0,9798 | 2,93% | 0,0461 |
| **SVR Otimizado (Original - 128 pts)** | 112 (CV) / 128 (Full) | **0,9979** | **0,93%** | **0,9865** | **2,47%** | **0,0376** |
| **SVR Otimizado + Augmentation (640 pts)** | 560 (CV) / 640 (Full) | **0,9983** | **0,90%** | 0,9446 | 4,54% | 0,0763 |

---

## 4. Análise Crítica dos Resultados (Por que o SVR caiu em CV?)

1. **Hiper-sensibilidade dos Vetores de Suporte:**
   * O SVR seleciona como *support vectors* os pontos na margem $\pm \epsilon$. Ao adicionar 640 pontos artificiais, a densidade de vetores de suporte aumentou dramaticamente, memorizando trajetórias específicas dos 14 ensaios de treino.
   * Quando confrontado com ensaios com parâmetros operacionais distintos ($\eta$ e $C_{A0}$ diferentes) no teste cego, o SVR perdeu capacidade de extrapolação generalizada, caindo de $R^2_{\text{CV}} = 0{,}9865$ para $0{,}9446$.
2. **Suavidade Visual:**
   * O SVR com Kernel RBF ($C^\infty$) sofre menos com quebras e cotovelos do que a rede neural MLP, mantendo uma curvatura assintótica lisa.
   * Contudo, o **SVR Original (sem augmentation)** é superior tanto visualmente quanto em termos de métricas estatísticas e integridade física.

---

## 5. Arquivos Gerados neste Diretório

* [cinetica_16_ensaios_svr_augmentation.png](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/cinetica_16_ensaios_svr_augmentation.png): Grade 4x4 do SVR treinado com 640 pontos.
* [comparativo_svr_com_vs_sem_augmentation.png](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/comparativo_svr_com_vs_sem_augmentation.png): Comparação visual direta (SVR Sem vs Com Augmentation).
* [cinetica_16_ensaios_mlp_tanh_augmentation.png](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/cinetica_16_ensaios_mlp_tanh_augmentation.png): Grade 4x4 do MLP com Data Augmentation.
* [comparativo_com_vs_sem_augmentation.png](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/comparativo_com_vs_sem_augmentation.png): Comparação visual direta do MLP.
* [gerar_graficos_svr_augmentation.py](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/gerar_graficos_svr_augmentation.py): Script de renderização para o SVR aumentado.
* [gerar_graficos_mlp_augmentation.py](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/gerar_graficos_mlp_augmentation.py): Script de renderização para o MLP aumentado.
* [data_augmentation.py](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/data_augmentation.py): Módulo gerador da base de 640 pontos PCHIP + PBM.
* [treinar_e_avaliar_mlp.py](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/treinar_e_avaliar_mlp.py): Script de validação cruzada comparativa.
* [dados_aumentados_640pts.csv](file:///c:/Users/GuilhermeMouradeSouz/Documents/lop_git/teste_augmentation/dados_aumentados_640pts.csv): Base com 640 instâncias aumentadas.
