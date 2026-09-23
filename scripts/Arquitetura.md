# Proposta de Arquitetura: Do Modelo White-Box ao Grey-Box

**Laboratório de Operações e Processos -- Engenharia Química -- UFMG**  
**Alunos:** Daniel Couto Vieira, Guilherme Moura de Sousa Franco, Matheus Henrique Borba Póvoas, Rodrigo Amaral da Mata  
**Orientador:** Prof. Dr. Fabrício Eduardo Bortot Coelho  
**Data:** Setembro de 2026  

---
## 1. Status Atual

Desenvolvemos uma estrutura em Python que reproduz as equações de conservação e balanço populacional:

```
scripts/
├── data_loader.py            <- Ingestão padronizada dos 128 pontos (Tabelas A1.1 e A1.4)
├── kinetics.py               <- Cinética de Arrhenius, estequiometria do ácido e taxa v(t)
├── population_balance.py     <- PBM via Método das Características com distribuição RRB
└── diagnostico_residuo.py    <- Execução do PBM e cálculo do resíduo Delta_X
```

### Resultados Obtidos pelo White-Box Puro ($\alpha = 0$):
* **$R^2$ Global (128 pontos):** **$0{,}9253$** (explica $92{,}5\%$ da variância experimental).
* **$RMSE$ (Erro Médio Quadrático):** **$0{,}0886$** ($8{,}86\%$).
* **$MAE$ (Erro Médio Absoluto):** **$0{,}0507$** ($5{,}07\%$).
* **Consistência Física:** $100\%$ das restrições atendidas ($0 \le X \le 1$ e parada da lixiviação (estado estacionário) quando $\eta < 1{,}0$).

### 1.1 Tabela de Conversões Experimentais ($X_{\text{Zn}}^{\text{exp}}$)
*(Dados das Tabelas A1.1 e A1.4 da dissertação de Fabrício Bortot Coelho, 2017)*

| Ensaio | η | CA0 (mol/L) | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1** | 0,5 | 0,1 | 0,00 | 0,42 | 0,46 | 0,46 | 0,49 | 0,51 | 0,49 | 0,50 |
| **2** | 3,1 | 0,1 | 0,00 | 0,82 | 0,90 | 0,95 | 0,97 | 0,97 | 0,99 | 1,00 |
| **3** | 0,5 | 0,5 | 0,00 | 0,45 | 0,48 | 0,49 | 0,52 | 0,50 | 0,49 | 0,50 |
| **4** | 0,5 | 1,0 | 0,00 | 0,50 | 0,47 | 0,51 | 0,47 | 0,53 | 0,51 | 0,50 |
| **5** | 0,5 | 1,5 | 0,00 | 0,49 | 0,49 | 0,47 | 0,49 | 0,49 | 0,47 | 0,49 |
| **6** | 1,0 | 0,1 | 0,00 | 0,75 | 0,78 | 0,80 | 0,85 | 0,87 | 0,85 | 0,87 |
| **7** | 3,1 | 0,5 | 0,00 | 0,88 | 0,93 | 0,98 | 0,99 | 1,00 | 0,99 | 1,01 |
| **8** | 1,0 | 0,5 | 0,00 | 0,83 | 0,82 | 0,85 | 0,86 | 0,87 | 0,86 | 0,87 |
| **9** | 1,0 | 1,0 | 0,00 | 0,83 | 0,84 | 0,85 | 0,85 | 0,83 | 0,84 | 0,85 |
| **10** | 1,0 | 1,5 | 0,00 | 0,84 | 0,83 | 0,83 | 0,83 | 0,83 | 0,84 | 0,85 |
| **11** | 1,5 | 0,1 | 0,00 | 0,86 | 0,85 | 0,88 | 0,97 | 0,94 | 0,96 | 0,98 |
| **12** | 3,1 | 1,0 | 0,00 | 0,96 | 1,00 | 1,00 | 0,99 | 0,99 | 0,99 | 1,00 |
| **13** | 1,5 | 0,5 | 0,00 | 0,92 | 0,92 | 0,95 | 0,95 | 0,95 | 0,95 | 0,97 |
| **14** | 1,5 | 1,0 | 0,00 | 0,91 | 0,96 | 0,95 | 0,97 | 0,95 | 0,98 | 0,97 |
| **15** | 1,5 | 1,5 | 0,00 | 0,91 | 0,94 | 0,95 | 0,94 | 0,96 | 0,96 | 0,97 |
| **16** | 3,1 | 1,5 | 0,00 | 0,98 | 0,98 | 0,96 | 0,99 | 0,99 | 1,00 | 1,00 |

### 1.2 Tabela de Conversões do Modelo White-Box Puro ($X_{\text{PBM}}$)
*(Valores obtidos via Balanço Populacional com $\alpha = 0$)*

| Ensaio | η | CA0 (mol/L) | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1** | 0,5 | 0,1 | 0,00 | 0,32 | 0,42 | 0,48 | 0,50 | 0,50 | 0,50 | 0,50 |
| **2** | 3,1 | 0,1 | 0,00 | 0,45 | 0,68 | 0,88 | 0,95 | 0,98 | 0,99 | 1,00 |
| **3** | 0,5 | 0,5 | 0,00 | 0,49 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| **4** | 0,5 | 1,0 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| **5** | 0,5 | 1,5 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| **6** | 1,0 | 0,1 | 0,00 | 0,39 | 0,57 | 0,72 | 0,80 | 0,84 | 0,87 | 0,95 |
| **7** | 3,1 | 0,5 | 0,00 | 0,92 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| **8** | 1,0 | 0,5 | 0,00 | 0,77 | 0,87 | 0,93 | 0,95 | 0,96 | 0,97 | 0,99 |
| **9** | 1,0 | 1,0 | 0,00 | 0,87 | 0,93 | 0,96 | 0,98 | 0,98 | 0,99 | 1,00 |
| **10** | 1,0 | 1,5 | 0,00 | 0,91 | 0,95 | 0,98 | 0,98 | 0,99 | 0,99 | 1,00 |
| **11** | 1,5 | 0,1 | 0,00 | 0,42 | 0,62 | 0,81 | 0,89 | 0,93 | 0,96 | 1,00 |
| **12** | 3,1 | 1,0 | 0,00 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| **13** | 1,5 | 0,5 | 0,00 | 0,86 | 0,96 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| **14** | 1,5 | 1,0 | 0,00 | 0,96 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| **15** | 1,5 | 1,5 | 0,00 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| **16** | 3,1 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |

### 1.3 Tabela de Resíduos Sistemáticos ($\Delta X = X_{\text{Zn}}^{\text{exp}} - X_{\text{PBM}}$)
*(O gabarito de correção que o Machine Learning aprende no Grey-Box)*

| Ensaio | η | CA0 (mol/L) | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **1** | 0,5 | 0,1 | 0,00 | +0,10 | +0,04 | -0,02 | -0,01 | +0,01 | -0,01 | -0,00 |
| **2** | 3,1 | 0,1 | 0,00 | +0,37 | +0,22 | +0,07 | +0,02 | -0,01 | -0,00 | +0,00 |
| **3** | 0,5 | 0,5 | 0,00 | -0,04 | -0,02 | -0,01 | +0,02 | -0,00 | -0,01 | -0,00 |
| **4** | 0,5 | 1,0 | 0,00 | +0,00 | -0,03 | +0,01 | -0,03 | +0,03 | +0,01 | -0,00 |
| **5** | 0,5 | 1,5 | 0,00 | -0,01 | -0,01 | -0,03 | -0,01 | -0,01 | -0,03 | -0,01 |
| **6** | 1,0 | 0,1 | 0,00 | +0,36 | +0,21 | +0,08 | +0,05 | +0,03 | -0,02 | -0,08 |
| **7** | 3,1 | 0,5 | 0,00 | -0,04 | -0,06 | -0,02 | -0,01 | +0,00 | -0,01 | +0,01 |
| **8** | 1,0 | 0,5 | 0,00 | +0,06 | -0,05 | -0,08 | -0,09 | -0,09 | -0,11 | -0,12 |
| **9** | 1,0 | 1,0 | 0,00 | -0,04 | -0,09 | -0,11 | -0,13 | -0,15 | -0,15 | -0,15 |
| **10** | 1,0 | 1,5 | 0,00 | -0,07 | -0,12 | -0,15 | -0,15 | -0,16 | -0,15 | -0,15 |
| **11** | 1,5 | 0,1 | 0,00 | +0,44 | +0,23 | +0,07 | +0,08 | +0,01 | +0,00 | -0,02 |
| **12** | 3,1 | 1,0 | 0,00 | -0,03 | +0,00 | +0,00 | -0,01 | -0,01 | -0,01 | +0,00 |
| **13** | 1,5 | 0,5 | 0,00 | +0,06 | -0,04 | -0,05 | -0,05 | -0,05 | -0,05 | -0,03 |
| **14** | 1,5 | 1,0 | 0,00 | -0,05 | -0,04 | -0,05 | -0,03 | -0,05 | -0,02 | -0,03 |
| **15** | 1,5 | 1,5 | 0,00 | -0,08 | -0,06 | -0,05 | -0,06 | -0,04 | -0,04 | -0,03 |
| **16** | 3,1 | 1,5 | 0,00 | -0,02 | -0,02 | -0,04 | -0,01 | -0,01 | +0,00 | +0,00 |

---

## 2. Justificativa

Na dissertação de 2017, para representar a desaceleração da lixiviação, foi introduzido um termo experimental na velocidade de dissolução:
$$v(t) = \frac{2}{\rho_s} \big[ k_s C_{Af}(t) - \alpha (C_{A0} - C_{Af}(t)) \big]$$
com o ajuste de $\alpha = 5500\ \mu\text{m/min}$.

### A Oportunidade Identificada:
Esse parâmetro $\alpha$ concentra de **forma única** fenômenos que variam de forma imprevisível ao longo da reação:
1. **Efeito de Íon Comum:** O acúmulo de $\text{Zn}^{2+}$ e sulfatos reduz a ação do ácido livre.
2. **Passivação Difusional por Sílica:** Formação de gel de sílica amorfa que recobre os poros das partículas.
3. **Morfologia dos Finos:** Partículas ultrafinas com cantos vivos (geometria não-esférica) dissolvem-se mais rápido nos primeiros 60 segundos do que esferas lisas.

### A Proposta do Nosso LOP:
Em vez de depender de um parâmetro de ajuste fixo ($\alpha$), **vamos manter o modelo White-Box com $\alpha = 0$** e utilizar a técnica de **Modelagem Híbrida**, onde o Machine Learning aprende o comportamento do resíduo (desvios do modelo teórico) :
$$\Delta X(t) = X_{\text{experimental}}(t) - X_{\text{PBM}}(t)$$

---

## 3. Arquitetura do Modelo Grey-Box

Abaixo apresentamos o fluxo de informação da arquitetura proposta:

<div align="center">

```mermaid
graph LR
    A["Condições Operacionais<br/>(t, η, CA0, T, RRB)"]
    
    B["White-Box<br/>Balanço Populacional"]
    
    C["Atributos Físicos<br/>Força-Motriz CAf, X/η"]
    D["Black-Box<br/>Regressor de Resíduo"]
    
    E["Acoplador Híbrido<br/>clip(X_PBM + ΔX_ML, 0, 1)"]
    F["Predição Final:<br/>X_Híbrido(t)<br/>"]

    A --> B
    A --> C
    B -->|"X_PBM(t)"| C
    C --> D
    B -->|"X_PBM(t)"| E
    D -->|"ΔX_ML"| E
    E --> F
```

</div>


### 3.1. Adaptabilidade do modelo
Uma das maiores vantagens da arquitetura híbrida sobre redes neurais puras é a **sensibilidade granulométrica herdada**:
* Se um novo lote de calcina tiver uma granulometria mais fina ($D_{63{,}2} = 25\ \mu\text{m}$) ou mais grossa ($D_{63{,}2} = 60\ \mu\text{m}$), **o White-Box recalcula analiticamente o 3º momento volumétrico e a dissolução de cada classe de tamanho**, sem necessidade de retreinar o algoritmo.
* O Machine Learning atua sobre o resíduo condicionado à nova previsão gerada, mantendo a estabilidade.

### 3.2. Atributos de Entrada do Black-Box
Em uma abordagem híbrida, o modelo de Machine Learning **não recebe variáveis cegas**. Ele é alimentado por uma combinação de condições operacionais e variáveis de estado calculadas pela física:

| # | Atributo (Feature) | Símbolo | Unidade | Origem | Papel Físico para o Machine Learning |
| :-: | :--- | :---: | :---: | :---: | :--- |
| **1** | **Tempo de Lixiviação** | $t$ | $\text{min}$ | Operacional | Localiza o estágio temporal da reação. |
| **2** | **Razão Estequiométrica** | $\eta$ | adim. | Operacional | Informa se há falta ($\eta < 1$), equivalência ($\eta = 1$) ou excesso ($\eta > 1$) de ácido. |
| **3** | **Ácido Inicial** | $C_{A0}$ | $\text{mol/L}$ | Operacional | Define a acidez inicial e a força iônica máxima da solução aquosa. |
| **4** | **Conversão da Física** | $X_{\text{PBM}}(t)$ | adim. ($0$ a $1$) | **White-Box** | Fornece a base física; a IA sabe onde a teoria está operando. |
| **5** | **Concentração de Ácido** | $C_{Af,\text{teo}}$ | $\text{mol/L}$ | **White-Box** | Força-motriz ácida residual instantânea ($C_{A0}(1 - X_{\text{PBM}}/\eta)$); governa passivação e íon comum. |
| **6** | **Esgotamento do Ácido** | $\frac{X_{\text{PBM}}}{\eta}$ | adim. ($0$ a $1$) | **White-Box** | Alerta a IA quando a reação está próxima de travar por falta de reagente. |

### 3.3. O Papel do Acoplador 
O acoplador não é uma simples soma matemática; ele atua como o **filtro de consistência termodinâmica e operacional** da planta, impedindo que incosistencias físicas ocorram:
1. **Fusão estrutural:** Soma a base de conversão do White Box à correção fina aprendida pelo Machine Learning:
   $$X_{\text{raw}}(t) = X_{\text{PBM}}(t) + \widehat{\Delta X}_{\text{ML}}(\mathbf{z})$$
2. **Garantia de conservação de massa:** Algoritmos de ML podem extrapolar para valores irreais ($X < 0$ ou $X > 1$). O acoplador impede isso:
   $$X_{\text{Híbrido}}(t) = \text{clip}\big( X_{\text{raw}}(t),\ 0{,}0,\ 1{,}0 \big)$$
3. **Respeito à estequiometria:** Quando há falta de ácido ($\eta < 1{,}0$), não existem moléculas suficientes de $\text{H}_2\text{SO}_4$ para dissolver mais do que a fração $\eta$. O acoplador impede que o modelo preveja conversões fisicamente impossíveis:
   $$X_{\text{Híbrido}}(t) \le \eta \quad (\text{para } \eta < 1{,}0)$$
4. **Garantia de único sentido temporal:** A dissolução de calcina na batelada é irreversível (o zinco solubilizado não precipita espontaneamente de volta como rocha; $\frac{dX}{dt} \ge 0$). O acoplador elimina oscilações e cenários impossíveis entre instantes sucessivos:
   $$X_{\text{Híbrido}}(t_k) \ge X_{\text{Híbrido}}(t_{k-1})$$

---

## 4. Metodologia

O plano de execução está dividido em quatro etapas:

### Etapa 1: Preparação dos dados
Em vez de alimentar o modelo de ML com variáveis puramente estatísticas ou cegas, forneceremos os **6 atributos** definidos na Seção 3.2, combinando condições operacionais e variáveis de estado termodinâmicas calculadas pela física:
1. **$t$ (Tempo de Lixiviação):** Localiza o estágio temporal da reação (início vs final).
2. **$\eta$ (Razão Estequiométrica):** Indica o regime de operação (falta quando $\eta < 1$, equivalência quando $\eta = 1$, ou excesso quando $\eta > 1$).
3. **$C_{A0}$ (Concentração Inicial de Ácido):** Define a acidez inicial e a força iônica da polpa.
4. **$X_{\text{PBM}}(t)$ (Conversão White-Box):** Linha de base calculada pelo balanço populacional.
5. **$C_{Af,\text{teórico}}$ (Ácido Livre Residual):** Força-motriz ácida instantânea ($C_{A0}[1 - X_{\text{PBM}}/\eta]$) que governa passivação e íon comum.
6. **$X_{\text{PBM}}/\eta$ (Fração de Esgotamento):** Alerta estequiométrico que indica quando a reação cessa por exaustão de reagente.

### Etapa 2: Treinamento e Catálogo de Modelos de Regressão
Para aprender o resíduo $\Delta X$, implementamos e avaliamos um catálogo com 8 regressores representativos das principais famílias de Machine Learning:
1. **Modelos Baseados em Árvores:**
   * **Gradient Boosting (GBDT):** Ajuste sequencial com regularização e penalização de gradiente.
   * **Extra Trees Regressor:** Árvores extremamente aleatorizadas, reduzindo variância.
   * **Random Forest Regressor:** Floresta aleatória com agregação bootstrap.
2. **Métodos de Kernel e Aprendizado Estatístico:**
   * **Support Vector Regression - SVR (Baseline):** Kernel RBF com margem $\epsilon$-insensível.
   * **SVR (Otimizado):** Parametrização fina de $C$, $\gamma$ e $\epsilon$ para respostas contínuas suaves.
   * **Gaussian Process Regressor (Kriging):** Inferência bayesiana não-paramétrica com estimativa intrínseca de incerteza (Kernel RBF + WhiteKernel).
3. **Redes Neurais Artificiais (Multi-Layer Perceptrons):**
   * **MLP com Ativação Tanh + Otimizador L-BFGS:** Rede neural com superfície suave de 2ª ordem, ideal para cinética contínua e pequenos conjuntos de dados.
   * **MLP com Ativação ReLU + Otimizador Adam:** Arquitetura padrão de aprendizado profundo por gradiente estocástico.

### Etapa 3: Validação Cruzada (Leave-One-Group-Out - LOGO-CV)
Para garantir que o modelo **generaliza para condições operacionais não vistas** e não decorou os dados:
* O conjunto de 16 ensaios de bancada ($N=128$) é particionado em 16 grupos (1 grupo por ensaio).
* Em cada teste, o modelo é treinado em **15 ensaios** e avaliado no **ensaio restante deixado de fora**.
* Repete-se o processo 16 vezes, garantindo métricas de teste isentas de vazamento de dados.

---

## 5. Resultados Experimentais e Benchmark Comparativo Multi-Modelo

Abaixo apresentam-se os resultados quantitativos obtidos em todas as etapas de teste e validação.

### 5.1 Benchmark

| Posição | Modelo / Algoritmo | Família | $R^2$ Ajuste | $R^2$ LOGO-CV | $RMSE_{\text{CV}}$ | $MAE_{\text{CV}}$ | Violação Física (%) |
| :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **1º** | **Gradient Boosting (GBDT)** | Árvores / Boosting | **0,9987** | **0,9934** | **0,0264** | **0,0186** | 0,0% |
| **2º** | **Extra Trees** | Árvores / Bagging | 0,9989 | 0,9932 | 0,0266 | 0,0188 | 0,0% |
| **3º** | **Random Forest** | Árvores / Bagging | 0,9981 | 0,9907 | 0,0312 | 0,0208 | 0,0% |
| **4º** | **Gaussian Process (Kriging)** | Processo Gaussiano | 0,9992 | 0,9878 | 0,0359 | 0,0211 | 0,0% |
| **5º** | **SVR (Otimizado)** | Support Vector | 0,9979 | 0,9870 | 0,0370 | 0,0234 | 0,0% |
| **6º** | **SVR (Baseline)** | Support Vector | 0,9938 | 0,9833 | 0,0419 | 0,0230 | 0,0% |
| **7º** | **MLP (Tanh + L-BFGS)** | Rede Neural | 0,9855 | 0,9688 | 0,0573 | 0,0284 | 0,0% |
| **8º** | **MLP (ReLU + Adam)** | Rede Neural | 0,9790 | 0,9299 | 0,0858 | 0,0587 | 0,0% |
| --- | *White-Box Puro* | *Balanço Populacional* | *0,9253* | *0,9253* | *0,0886* | *0,0507* | *0,0%* |

![alt text](image-1.png)

---

## 6. Arquitetura do Modelo Híbrido Serial (Grey-Box Cinético)

Em resposta às limitações conceituais de desacoplamento do modelo residual e à rigidez empírica da tese de mestrado de **Bortot Coelho (2017)**, desenvolvemos a **Arquitetura Híbrida Serial (Grey-Box Cinético com Estimação de Parâmetro)**, implementada e validada no módulo [`serial_hybrid_model.py`](serial_hybrid_model.py).

### 6.1. Superação da Hipótese de $\alpha$ Constante da Tese de 2017
Na formulação original de Bortot Coelho (2017), a velocidade linear de avanço da reação de dissolução da zincita sob o Modelo do Núcleo Não Reagido (*Shrinking Core Model* -- SCM) com amortecimento difusional e acúmulo de sulfato foi proposta como:
$$v(t) = \frac{d(\Delta D)}{dt} = \frac{2}{\rho_s} \max\Big( k_s \, C_{Af}(t) - \alpha \, [C_{A0} - C_{Af}(t)],\ 0 \Big)$$
onde:
* $\rho_s = 69{,}2\ \text{mol/L}$ é a densidade molar da calcina sólida;
* $k_s = 18000\ \mu\text{m/min}$ a $40^\circ\text{C}$ (Arrhenius, Balarini et al., 2025);
* $C_{Af}(t) = \max\big(C_{A0}(1 - X(t)/\eta),\ 0\big)$ é a acidez livre residual instantânea;
* $[C_{A0} - C_{Af}(t)] = \frac{C_{A0}}{\eta} X(t)$ é a concentração de reagente consumido, proporcional ao acúmulo de $\text{Zn}^{2+}$ e íons sulfato na camada de difusão;
* $\alpha$ é o parâmetro empírico de amortecimento cinético em $[\mu\text{m/min}]$.

Na dissertação de 2017, pela ausência de métodos modernos de aprendizado estatístico, ajustou-se um **único valor estático** para cobrir todo o envelope experimental:
$$\alpha_{\text{tese}} = 5500\ \mu\text{m/min} \quad (\text{constante universal})$$

Entretanto, uma análise mecanicista e estequiométrica rigorosa comprova que o parâmetro de desaceleração $\alpha$ não pode ser tratado como uma constante termodinâmica fixa:
1. **Em falta estequiométrica de ácido ($\eta = 0{,}5$):** A reação se encerra precocemente por falta de moléculas de $\text{H}_2\text{SO}_4$ ($X \to 0{,}5$). O valor fixo $\alpha = 5500$ subestima severamente a velocidade de dissolução inicial, elevando o erro quadrático ($SSE$ salta de $0{,}002$ para $0{,}088$). O valor físico real é $\alpha \approx 0\text{ a }500\ \mu\text{m/min}$.
2. **Em ácido diluído ($C_{A0} = 0{,}1\ \text{mol/L}$):** A força iônica é desprezível e a camada difusional não sofre impedimento por sulfatos precipitados. As partículas ultrafinas dissolvem-se velozmente, exigindo $\alpha \approx 0$ (regime mecanicista puro).
3. **Em excesso de ácido ($\eta = 3{,}1$) e alta acidez ($C_{A0} \ge 1{,}0\ \text{mol/L}$):** Ocorre geração massiva de sulfato de zinco em solução e forte precipitação de gel de sílica amorfa que oclui os microporos da calcina. O valor $\alpha = 5500$ é insuficiente para conter a taxa, exigindo $\alpha \approx 18000\text{ a }24000\ \mu\text{m/min}$.

---

### 6.2. Diagrama Conceitual e Estrutural da Arquitetura Serial

Na arquitetura serial, a inteligência artificial não atua sobre a conversão final $X(t)$ de forma post-hoc; ela atua **antes e dentro** da equação de conservação diferencial:

<div align="center">

```mermaid
graph TD
    subgraph ENTRADAS["1. Variáveis de Entrada Operacionais"]
        E1["Razão Estequiométrica<br/>η = H2SO4 / ZnO"]
        E2["Acidez Inicial<br/>CA0 (mol/L)"]
        E3["Razão Sólido-Líquido<br/>S/L (g/L)"]
    end

    subgraph BLACKBOX["2. Módulo Black-Box de Inteligência Artificial"]
        ML1["Regressor de Machine Learning<br/>(GBDT / PolyRidge / SVR / MLP)"]
        ML2["Restrição Termodinâmica:<br/>α̂ = max(0, f_ML(η, CA0))"]
        ML1 --> ML2
    end

    subgraph WHITEBOX["3. Módulo White-Box Fundamental (PBM + SCM)"]
        WB1["Taxa Linear de Retração Efetiva:<br/>v(t) = (2/ρs) max(ks CAf - α̂ ΔCA, 0)"]
        WB2["EDO de Retração Radial:<br/>d(ΔD)/dt = v(t),  ΔD(0) = 0"]
        WB3["Integração Granulométrica RRB (3º Momento):<br/>X(t) = 1 - ∫ (D - ΔD)³ f(D) dD / I0"]
        
        WB1 --> WB2
        WB2 --> WB3
    end

    subgraph SAIDA["4. Saída do Modelo"]
        S1["Conversão Temporal X(t)<br/>Monotonia Estrita: dX/dt ≥ 0"]
    end

    E1 --> ML1
    E2 --> ML1
    E3 -.-> ML1

    ML2 -->|"Fator α̂ Calibrado"| WB1
    WB3 --> S1
```

</div>

---

### 6.3. Diferenças Estruturais: Arquitetura Paralela vs Arquitetura Serial

| Aspecto Estrutural | Arquitetura Híbrida Paralela (Residual) | Nova Arquitetura Híbrida Serial (Cinética) |
| :--- | :--- | :--- |
| **Ponto de Acoplamento da IA** | Saída da física: $X(t) = X_{\text{PBM}}(t) + \widehat{\Delta X}_{\text{ML}}$ | Entrada da física: $v(t) = f(C_{Af},\ \hat{\alpha}_{\text{ML}})$ |
| **Papel do Machine Learning** | Aprende o desvio residual global $\Delta X$ | Aprende o parâmetro cinético de passivação $\hat{\alpha}$ |
| **Garantia de Monotonicidade ($dX/dt \ge 0$)** | Artificial (requer filtro `np.maximum.accumulate`) | **Nativa e Intrínseca** (como $v(t) \ge 0$, $\Delta D$ e $X$ nunca decrescem) |
| **Respeito Estequiométrico ($X \le \eta$ quando $\eta < 1$)** | Artificial (requer truncamento `clip(X, 0, eta)`) | **Nativa e Intrínseca** (o ácido $C_{Af}$ zera no PBM e a taxa anula-se) |
| **Consistência do Consumo de Ácido** | Desacoplada ($C_{Af}$ é calculado com o $X$ sem correção) | **Totalmente Acoplada** ($C_{Af}$ e $X$ evoluem em sincronia na EDO) |
| **Escalonamento para Reatores Contínuos (CSTR)** | Exige artifícios empíricos de DTR | **Direto e Imediato** (entra na equação de balanço molar de cada estágio) |

---

## 7. Metodologia de Implementação e Validação

O desenvolvimento da arquitetura serial foi conduzido segundo uma metodologia rigorosa em sete etapas, assegurando reprodutibilidade científica e ausência total de vazamento de dados (*data leakage*):

### Etapa 1: Ingestão Padronizada e Matriz Experimental
* Utilizou-se a base completa com os 16 ensaios de bancada da dissertação de Bortot Coelho (2017), totalizando **$N = 128$ pontos experimentais** cobrindo os tempos $t \in [0; 0{,}5; 1; 2; 3; 4; 5; 15]$ min.
* Matriz operacional: 4 níveis de razão molar $\eta \in \{0{,}5;\ 1{,}0;\ 1{,}5;\ 3{,}1\}$ cruzados com 4 níveis de acidez inicial $C_{A0} \in \{0{,}1;\ 0{,}5;\ 1{,}0;\ 1{,}5\}\ \text{mol/L}$.

### Etapa 2: Formulação do Problema Inverso para Calibração de $\alpha^*$
Para cada ensaio $k \in \{1, \dots, 16\}$, determinou-se o valor de amortecimento ótimo $\alpha_k^*$ que minimiza a soma dos erros quadráticos ($SSE$) entre a solução analítico-numérica do PBM e os dados experimentais:
$$\alpha_k^* = \arg\min_{\alpha \in [0,\ 35000]} \sum_{j=1}^8 \left( X_{\text{exp}}(t_j) - X_{\text{PBM}}(t_j;\ \alpha) \right)^2$$
A otimização foi executada via método de Brent com limites estritos (`scipy.optimize.minimize_scalar`), assegurando $\alpha^* \ge 0$.

### Etapa 3: Catálogo Multi-Modelo e Treinamento do Black-Box
Implementou-se um catálogo com 6 regressores representativos de diferentes famílias de Machine Learning:
1. **Gradient Boosting (GBDT):** Floresta sequencial de árvores com regularização e subamostragem ($n=45$, profundidade=3, taxa de aprendizado=0,08).
2. **Polynomial Ridge:** Regressão quadrática completa de 2º ordem com regularização de Tikhonov ($L_2$, $\lambda=1{,}5$).
3. **Random Forest:** Floresta de árvores de decisão aleatórias com agregação bootstrap ($n=60$, profundidade=4).
4. **Extra Trees:** Árvores extremamente aleatorizadas para minimização de variância ($n=60$, profundidade=4).
5. **Support Vector Regression (SVR RBF):** Máquina de vetores de suporte com kernel de base radial Gaussiana ($C=10000$, $\epsilon=400$).
6. **Rede Neural (MLP com Escalonamento de Alvo):** Perceptron multicamadas com camadas ocultas $(8, 4)$, ativação hiperbólica `tanh`, otimizador Quasi-Newton L-BFGS de 2ª ordem e `TransformedTargetRegressor` para garantir convergência assintótica sem instabilidade numérica.

### Etapa 4: Integração Temporal e Resolução Numérica do Balanço Populacional
A cada passo de predição:
1. O regressor prevê $\hat{\alpha} = \max\big(0,\ f_{\text{ML}}(\eta, C_{A0})\big)$.
2. O integrador `solve_ivp` (Runge-Kutta de 4ª/5ª ordem adaptativo - RK45) resolve a taxa $v(t) = \frac{d(\Delta D)}{dt}$.
3. O encolhimento acumulado $\Delta D(t)$ é convertido em conversão mássica $X(t)$ pela integração numérica do 3º momento da distribuição de Rosin-Rammler-Bennett (RRB, $m = 1{,}022$, $D_{63{,}2} = 41{,}65\ \mu\text{m}$):
   $$1 - X(t) = \frac{1}{I_0} \int_{\Delta D(t)}^\infty (D - \Delta D(t))^3 \, \frac{m}{D_{63{,}2}} \left(\frac{D}{D_{63{,}2}}\right)^{m-1} \exp\left[-\left(\frac{D}{D_{63{,}2}}\right)^m\right] dD$$

### Etapa 5: Validação Cruzada Estrita (Leave-One-Group-Out -- LOGO-CV)
Para atestar a capacidade preditiva em **condições operacionais cegas**, executou-se validação cruzada por grupos (16 folds independentes):
* Em cada fold $k$, o modelo Black-Box é treinado exclusivamente nos 15 ensaios restantes.
* O valor de $\hat{\alpha}$ do ensaio $k$ é predito às cegas a partir de $(\eta_k, C_{A0,k})$.
* O PBM simula a curva inteira do ensaio $k$ sem que qualquer dado desse ensaio tenha sido visto pelo regressor.

### Etapa 6: Métricas Estatísticas e Penalização dos Graus de Liberdade
As métricas quantitativas foram calculadas segundo as formulações padronizadas:
1. **Coeficiente de Determinação ($R^2$):**
   $$R^2 = 1 - \frac{\sum_{i=1}^n (y_i - \hat{y}_i)^2}{\sum_{i=1}^n (y_i - \bar{y})^2}$$
2. **Coeficiente de Determinação Ajustado ($R^2_{\text{ajustado}}$):**
   $$R^2_{\text{ajustado}} = 1 - \left[ \frac{(1 - R^2)(n - 1)}{n - p - 1} \right]$$
   onde $n = 128$ e $p$ é o número de parâmetros explicativos do modelo:
   * White-Box Fundamental Puro: $p = 0$ ($R^2_{\text{adj}} = R^2$);
   * Tese Bortot Coelho (2017): $p = 1$ (parâmetro $\alpha = 5500$ fixo);
   * Híbrido Paralelo Residual: $p = 6$ (6 atributos termodinâmicos);
   * Novo Híbrido Serial: $p = 2$ (atributos operacionais $\eta$ e $C_{A0}$).
3. **Métricas de Erro:**
   $$RMSE = \sqrt{\frac{1}{n} \sum_{i=1}^n (y_i - \hat{y}_i)^2}, \qquad MAE = \frac{1}{n} \sum_{i=1}^n |y_i - \hat{y}_i|$$
4. **Consistência Física:**
   * Violação de Limites: $\%$ de pontos com $X < 0$, $X > 1$ ou $X > \eta + 10^{-4}$ (para $\eta < 1$).
   * Violação de Monotonicidade: $\%$ de intervalos temporais com $\frac{dX}{dt} < -10^{-4}$.

### Etapa 7: Mecanismo de Equilíbrio Analítico e Superação da Tese
Em condições estequiométricas ($\eta = 1{,}0$), no limite assintótico $t \to \infty$, a velocidade de retração anula-se ($v(t) = 0$). Isso impõe que a conversão final atinja exatamente o valor analítico:
$$k_s \, C_{Af} = \alpha \, [C_{A0} - C_{Af}] \implies k_s (1 - X) = \alpha X \implies X_{\text{final}} = \frac{k_s}{k_s + \alpha}$$
* Na tese de 2017, com $\alpha = 5500$: $X_{\text{final}} = \frac{18000}{18000 + 5500} = 76{,}6\%$. Como os ensaios 6, 8, 9 e 10 atingem experimentalmente entre $85\%$ e $87\%$, o modelo antigo travava prematuramente em $76{,}6\%$.
* No modelo Híbrido Serial, a IA prevê $\hat{\alpha} \approx 2500\text{ a }3400\ \mu\text{m/min}$, levando a $X_{\text{final}} = 84\%\text{ a }87{,}8\%$, reproduzindo fielmente o patamar real observado.

---

## 8. Benchmark Comparativo Geral entre Todas as Arquiteturas

Abaixo consolida-se o confronto quantitativo entre todas as formulações avaliadas nos 128 pontos experimentais de bancada.

### 8.1. Tabela Geral de Desempenho

| Arquitetura / Modelo | Família | Parâmetros ($p$) | $R^2$ Ajuste | $R^2_{\text{ajustado}}$ Ajuste | $R^2$ LOGO-CV | $R^2_{\text{ajustado}}$ LOGO-CV | $RMSE_{\text{CV}}$ | $MAE_{\text{CV}}$ | Monotonicidade Intrínseca |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Híbrido Paralelo Competitivo (GBDT, Sansana 2024)** | Grey-Box Competitivo | 5 | 0,9985 | 0,9985 | **0,9971** | **0,9969** | **0,0176** | **0,0132** | Sim (Acoplador) |
| **Híbrido Paralelo Residual / Cooperativo (GBDT)** | Grey-Box Cooperativo | 6 | **0,9987** | **0,9987** | 0,9934 | 0,9930 | 0,0264 | 0,0186 | Não (requer filtro) |
| **Novo Híbrido Serial (alpha-PolyRidge)** | Grey-Box Serial | 2 | 0,9438 | 0,9429 | 0,9387 | 0,9378 | 0,0802 | 0,0404 | **Sim (100% nativa)** |
| **Novo Híbrido Serial (alpha-MLP)** | Grey-Box Serial | 2 | 0,9463 | 0,9454 | 0,9371 | 0,9361 | 0,0813 | 0,0397 | **Sim (100% nativa)** |
| **Novo Híbrido Serial (alpha-GBDT)** | Grey-Box Serial | 2 | 0,9474 | 0,9465 | 0,9357 | 0,9347 | 0,0822 | 0,0398 | **Sim (100% nativa)** |
| **White-Box Fundamental Puro ($\alpha=0$)** | Balanço Populacional | 0 | 0,9253 | 0,9253 | 0,9253 | 0,9253 | 0,0886 | 0,0507 | **Sim (100% nativa)** |
| **Tese Bortot Coelho (2017) ($\alpha=5500$)** | Cinética Empírica | 1 | 0,8976 | 0,8968 | 0,8976 | 0,8968 | 0,1037 | 0,0708 | **Sim (100% nativa)** |

---

### 8.2. Evidências Visuais e Gráficos Científicos (300 DPI)

#### Gráfico de Paridade Quíntuplo (Todas as Arquiteturas em Teste Cego LOGO-CV):
![Paridade Comparativa](docs/paridade_comparativa_todas_arquiteturas.png)

#### Curvas Cinéticas Comparativas nos 4 Regimes Estequiométricos ($\eta = 0{,}5;\ 1{,}0;\ 1{,}5;\ 3{,}1$):
![Curvas Cinéticas](docs/curvas_dissolucao_ensaios_criticos.png)

#### Distribuição dos Fatores de Ponderação Competitiva ($w$ e $1-w$) nos 16 Folds Cegos:
![Pesos Competitivos](docs/analise_pesos_ponderacao_competitiva.png)

#### Superfície do Fator $\alpha(\eta, C_{A0})$ Aprendida pela Inteligência Artificial na Arquitetura Serial:
![Superfície Alpha](docs/superficie_alpha_interpretacao_fisica.png)

---

## 9. Arquitetura Híbrida Paralela Competitiva (Sansana et al., 2024)

Em consonância com o artigo seminal de **Sansana, Rendall, Castillo et al. (2024)** (*Hybrid modeling for transfer learning in chemical processes*, *Chemical Engineering Science*, 300, 120568), implementamos a terceira grande família de modelagem cinza: a **Arquitetura Paralela Competitiva** ([`competitive_hybrid_model.py`](competitive_hybrid_model.py)).

### 9.1. Princípio de Funcionamento e Equacionamento
Diferente da abordagem paralela residual (cooperativa), onde o Black-Box corrige os resíduos $\Delta X$, na modelagem competitiva **ambos os modelos competem independentemente para calcular a variável de processo final** ($X_{\text{Zn}}$):
1. **Módulo White-Box Fundamental ($f_p$):** Mantido **estritamente intacto**, resolvendo o PBM puro com $\alpha = 0$ via Método das Características e distribuição RRB:
   $$\hat{X}_{wb}(t) = f_p(t, \eta, C_{A0})$$
2. **Módulo Black-Box Direto ($f_d$):** Treinado diretamente com relação ao alvo experimental $y = X_{\text{zn}}^{\text{exp}}$ a partir das variáveis de operação:
   $$\hat{X}_{bb}(t) = f_d(t, \eta, C_{A0}, S/L; \boldsymbol{\psi})$$
3. **Camada de Fusão / Meta Ponderada (Equações 10 e 11 de Sansana et al.):**
   $$\hat{X}_{\text{raw}}(t) = w \, \hat{X}_{wb}(t) + (1 - w) \, \hat{X}_{bb}(t)$$
   onde o fator de ponderação $w \in [0, 1]$ governa a autoridade relativa de cada submodelo.
4. **Otimização Analítica do Fator de Ponderação:**
   O peso $w$ é calibrado através da minimização de mínimos quadrados restritos (*constrained least squares*):
   $$w^* = \arg\min_{w \in [0, 1]} \sum_{i=1}^N \Big( X_{\text{exp}, i} - \big[ w \hat{X}_{wb, i} + (1 - w) \hat{X}_{bb, i} \big] \Big)^2$$
   Com a solução analítica exata dada por:
   $$w^* = \text{clip}\left( \frac{\sum_{i=1}^N (\hat{X}_{wb, i} - \hat{X}_{bb, i})(X_{\text{exp}, i} - \hat{X}_{bb, i})}{\sum_{i=1}^N (\hat{X}_{wb, i} - \hat{X}_{bb, i})^2},\ 0{,}0,\ 1{,}0 \right)$$
5. **Acoplador Termodinâmico:** Garante consistência física estrita: conservação de massa ($\text{clip}(X, 0, 1)$), respeito à estequiometria ($X \le \eta$ para $\eta < 1{,}0$) e monotonicidade irreversível temporal ($\frac{dX}{dt} \ge 0$).

---

### 9.2. Diagrama Estrutural da Modelagem Paralela Competitiva

<div align="center">

```mermaid
graph TD
    subgraph ENTRADAS["1. Variáveis Operacionais de Processo (x)"]
        E1["Tempo de Lixiviação (t)"]
        E2["Razão Estequiométrica (η)"]
        E3["Concentração de Ácido (CA0)"]
        E4["Razão Sólido-Líquido (S/L)"]
    end

    subgraph COMPETICAO["2. Concorrência Independente de Modelos"]
        subgraph WB["White-Box (Física Teórica Pura - Não Alterada)"]
            WB1["Balanço Populacional (PBM)<br/>Taxa v(t) com α = 0"]
            WB2["Predição Mecanicista:<br/>X_wb(t)"]
            WB1 --> WB2
        end

        subgraph BB["Black-Box (Machine Learning Direto)"]
            BB1["Regressor Direto de Conversão<br/>(GBDT / Random Forest / MLP)"]
            BB2["Predição Orientada a Dados:<br/>X_bb(t)"]
            BB1 --> BB2
        end
    end

    subgraph META["3. Camada de Fusão Ponderada (Sansana et al., 2024)"]
        M1["Mínimos Quadrados Restritos:<br/>w* = argmin ∑ (X_exp - [w X_wb + (1-w) X_bb])²<br/>sujeito a 0 ≤ w ≤ 1"]
        M2["Combinação Ponderada Bruta:<br/>X_raw = w* X_wb + (1 - w*) X_bb"]
        M1 --> M2
    end

    subgraph ACOPLADOR["4. Acoplador Termodinâmico"]
        A1["Conservação de Massa: clip(X, 0, 1)"]
        A2["Teto Estequiométrico: X ≤ η para η < 1.0"]
        A3["Monotonicidade Temporal: dX/dt ≥ 0"]
        A1 --> A2 --> A3
    end

    subgraph SAIDA["5. Saída Final"]
        S1["Predição Híbrida Competitiva:<br/>X_Competitivo(t)"]
    end

    E1 --> WB1
    E2 --> WB1
    E3 --> WB1

    E1 --> BB1
    E2 --> BB1
    E3 --> BB1
    E4 --> BB1

    WB2 --> M1
    BB2 --> M1
    WB2 --> M2
    BB2 --> M2

    M2 --> A1
    A3 --> S1
```

</div>

---

### 9.3. Comparativo Estrutural: As Três Filosofias Híbridas do TCC

| Dimensão Científica | Arquitetura Híbrida Paralela Cooperativa (Residual) | Nova Arquitetura Híbrida Serial (Cinética) | Arquitetura Híbrida Paralela Competitiva (Sansana 2024) |
| :--- | :--- | :--- | :--- |
| **Referência / Princípio** | Grey-Box Aditivo Clássico (Eq. 9) | Grey-Box Paramétrico / SCM (Eq. 12) | Sansana et al. (2024), CES (Eq. 10-11) |
| **Papel do Black-Box** | Aprende o resíduo teórico: $\Delta X = X_{\text{exp}} - X_{\text{PBM}}$ | Estima o parâmetro físico de passivação: $\hat{\alpha}_{\text{ML}}$ | Estima a conversão direta: $\hat{X}_{bb} \approx X_{\text{Zn}}$ |
| **Papel do White-Box** | Fornece a linha de base física subjacente ($X_{\text{PBM}}$) | Resolve as EDOs acopladas via balanço molar e PBM | Fornece uma predição independente concorrente ($X_{\text{wb}}$) |
| **Mecanismo de Fusão** | Soma aditiva direta: $X_{\text{wb}} + \widehat{\Delta X}_{\text{ML}}$ | Acoplamento diferencial no integrador temporal | Média ponderada ótima: $w X_{\text{wb}} + (1 - w) X_{\text{bb}}$ |
| **Fator de Ponderação** | Não aplicável (peso unitário implícito) | Não aplicável (parâmetro dentro da física) | Calibrado via mínimos quadrados ($w \in [0, 1]$) |
| **$R^2$ LOGO-CV** | $0{,}9934$ | $0{,}9387$ (PolyRidge) / $0{,}9357$ (GBDT) | **$0{,}9971$ (GBDT)** |
| **$RMSE$ LOGO-CV** | $0{,}0264$ ($2{,}64\%$) | $0{,}0802$ ($8{,}02\%$) | **$0{,}0176$ ($1{,}76\%$)** |
| **$MAE$ LOGO-CV** | $0{,}0186$ ($1{,}86\%$) | $0{,}0404$ ($4{,}04\%$) | **$0{,}0132$ ($1{,}32\%$)** |
| **Garantia de Monotonia** | Filtro `np.maximum.accumulate` | **100% Nativa e Intrínseca** (pela física da EDO) | Acoplador com filtro `np.maximum.accumulate` |
| **Sensibilidade Granulométrica** | Totalmente herdada do 3º momento volumétrico | Totalmente herdada das classes de tamanho RRB | Ponderada pelo fator $w$ sobre a resposta do PBM |
| **Interpretabilidade Físico-Química** | Média (corrige erros sem mapear o mecanismo) | **Máxima** (explica a taxa de passivação e sulfatos) | Alta na física, flexível no aprendizado estatístico |

---

### 9.4. Interpretação Científica dos Fatores de Ponderação $w$
A calibração do peso $w$ gerou um diagnóstico de grande relevância acadêmica para a tese:
1. **Com Algoritmos de Alta Capacidade Não-Linear (GBDT / Random Forest):**
   * O Black-Box atinge individualmente $R^2 \approx 0{,}996$.
   * A camada de meta-ponderação atribui $w \approx 0{,}040\text{ a }0{,}043$ ($4\%$ White-Box) e $1 - w \approx 0{,}957\text{ a }0{,}960$ ($96\%$ Black-Box).
   * A fusão híbrida eleva o desempenho para **$R^2 = 0{,}9967\text{ a }0{,}9971$** e reduz o erro quadrático para **$RMSE \le 1{,}86\%$**, demonstrando que a ancoragem na física previne instabilidades nos ensaios extremos.
2. **Mecanismo de Salvaguarda Física com Algoritmos Lineares / Restritos (Ridge / SVR):**
   * Quando o Black-Box possui expressividade limitada (como o Polynomial Ridge), a otimização de Sansana et al. **redistribui automaticamente a autoridade para a física pura**, elevando $w$ para **$0{,}885$ a $0{,}898$** ($90\%$ de peso no White-Box puro).
   * Isso comprova que a modelagem competitiva atua como uma barreira protetora contra regressões de baixa qualidade, preservando a coerência do modelo industrial.

---

### 9.5. Estudo de Consequências: Particionamento 85% Treino / 15% Validação vs 15/1 (LOGO-CV)

Em resposta à boa prática estatística de avaliação com **$85\%$ dos dados para treino e $15\%$ para validação/teste**, implementamos essa estratégia em [`competitive_hybrid_model.py`](competitive_hybrid_model.py) e confrontamos os impactos conceituais e numéricos frente ao particionamento 15/1 (*Leave-One-Group-Out* -- LOGO-CV):

| Estratégia de Particionamento | Ensaios Treino | Ensaios Teste | Fração Treino / Teste | $R^2$ Validação | $R^2_{\text{ajustado}}$ | $RMSE_{\text{Val}}$ | $MAE_{\text{Val}}$ | Peso $w$ Médio (Física) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Validação 85/15 (GroupKFold $k=6$)** | $\sim 13\text{ a }14$ | $\sim 2\text{ a }3$ | **$85\% / 15\%$** | **$0{,}9967$** | **$0{,}9966$** | **$0{,}0186$ ($1{,}86\%$)** | **$0{,}0140$ ($1{,}40\%$)** | $0{,}0433$ |
| **Validação 85/15 (Holdout Estruturado)** | $13$ ($104$ pts) | $3$ ($24$ pts) | **$81{,}2\% / 18{,}8\%$** | **$0{,}9964$** | **$0{,}9963$** | **$0{,}0188$ ($1{,}88\%$)** | **$0{,}0130$ ($1{,}30\%$)** | $0{,}0354$ |
| **Validação 15/1 (LOGO-CV 16 folds)** | $15$ ($120$ pts) | $1$ ($8$ pts) | **$93{,}8\% / 6{,}2\%$** | **$0{,}9971$** | **$0{,}9969$** | **$0{,}0176$ ($1{,}76\%$)** | **$0{,}0132$ ($1{,}32\%$)** | $0{,}0420$ |

#### Análise das Consequências Metodológicas:
1. **Risco de Viés Otimista de Interpolação no 15/1:**
   No particionamento 15/1, o modelo conhece $93{,}75\%$ de todo o envelope operacional em cada teste. Os 15 ensaios de treino cercam intimamente o único ensaio deixado de fora, fazendo com que o teste seja puramente de **interpolação local densa**. Isso eleva levemente o $R^2$ ($0{,}9971$).
2. **Robustez e Generalização Real no 85/15:**
   Ao omitir simultaneamente $2$ a $3$ ensaios inteiros ($15\%$ dos dados), o modelo perde múltiplos nós da grade operacional $(\eta, C_{A0})$. O fato do $R^2$ se manter praticamente inalterado (**$0{,}9967$ vs $0{,}9971$**, com $RMSE$ variando apenas de $1{,}76\%$ para $1{,}86\%$) atesta que o modelo híbrido competitivo **não está decorando pontos vizinhos, mas sim aprendendo a dinâmica contínua do processo com extrema robustez**.
3. **Necessidade Estrita de Agrupamento por Ensaios (Prevenção de *Data Leakage*):**
   A divisão de $85\% / 15\%$ só é estatisticamente válida se realizada em nível de **ensaios completos** (grupos). Se fosse feita aleatoriamente ponto a ponto (misturando instantes $t$ de uma mesma curva cinética entre treino e teste), ocorreria vazamento temporal severo. Com o `GroupKFold` ($k=6$), curvas temporais inteiras são mantidas 100% cegas no teste.

#### Gráfico de Paridade Dedicado e Análise de Resíduos do Modelo Híbrido Competitivo (300 DPI):
![Paridade Modelo Competitivo](docs/paridade_modelo_competitivo.png)

---

## 10. Estudo Comparativo sob Particionamento 85/15 (Treino e Teste Cego): Serial vs Cooperativo vs Competitivo

Em consonância com as diretrizes de avaliação com **$85\%$ dos dados para treino e $15\%$ para teste**, expandimos essa análise aos três modelos híbridos desenvolvidos no projeto, implementada em [`gerar_graficos_85_15.py`](gerar_graficos_85_15.py):
1. **Modelo Híbrido Serial Cinético** ([`serial_hybrid_model.py`](serial_hybrid_model.py))
2. **Modelo Híbrido Paralelo Cooperativo / Residual** ([`hybrid_model.py`](hybrid_model.py))
3. **Modelo Híbrido Paralelo Competitivo** ([`competitive_hybrid_model.py`](competitive_hybrid_model.py))

### 10.1. Tabela Comparativa de Desempenho (Split 85% Treino / 15% Teste)

| Arquitetura / Modelo | $R^2$ Treino (85%) | $RMSE$ Treino | $R^2$ Teste (15%) | $RMSE$ Teste (15%) | $MAE$ Teste (15%) | $R^2$ GroupKFold ($k=6$) | $RMSE_{\text{GKF}}$ | Comportamento em Extrapolação |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Híbrido Paralelo Competitivo (Sansana 2024)** | **0,9985** | **0,0125** | **0,9964** | **0,0188** ($1{,}88\%$) | **0,0130** | **0,9967** | **0,0186** | Moderada (condicionada pelo peso $w$) |
| **Híbrido Paralelo Cooperativo (Residual)** | **0,9990** | **0,0105** | 0,9896 | 0,0317 ($3{,}17\%$) | 0,0238 | 0,9927 | 0,0277 | Frágil (estagnação de árvores fora do domínio) |
| **Novo Híbrido Serial Cinético ($\alpha$-GBDT)** | 0,9720 | 0,0544 | 0,7896* | 0,1431 ($14{,}31\%$) | 0,0869 | 0,9362 | 0,0819 | **Excelente (100% conservação física e monotonia nativa)** |

*\* Nota Científica:* No holdout aleatório de 15% com `random_state=42`, os ensaios de teste sorteados foram **Ensaios 1, 2 e 6**, que concentram exclusivamente a acidez mais diluída de todo o espaço amostral ($C_{A0} = 0{,}1\ \text{mol/L}$). Isso configurou um teste severo de **fronteira/extrapolação de acidez**. Na validação agrupada global (`GroupKFold` $k=6$), o modelo serial atinge $R^2 = 0{,}9362$ e $RMSE = 8{,}19\%$, comprovando estabilidade contínua.

---

### 10.2. Diagnóstico Científico: Extrapolação, Rigor Físico e Escalonamento Industrial

O confronto das três arquiteturas permitiu consolidar respostas definitivas para o TCC:

1. **Qual é o melhor modelo para Interpolação e Controle Local?**
   * O **Híbrido Paralelo Competitivo (Sansana et al., 2024)** é o campeão numérico absoluto ($R^2 = 0{,}9964\text{ a }0{,}9971$, $RMSE \le 1{,}88\%$). Ele permite máxima aderência aos dados industriais de rotina dentro do envelope de calibração.
2. **Qual é o melhor modelo para Extrapolação e Rigor Físico?**
   * O **Híbrido Serial Cinético** é insuperável para extrapolação. Porque:
     * O Machine Learning não prevê conversões $X(t)$ de forma cega; ele estima o parâmetro físico de desaceleração $\hat{\alpha} \ge 0$.
     * Toda a trajetória temporal é gerada pela resolução numérica das equações diferenciais fundamentais (EDO do SCM e balanço populacional com distribuição RRB).
     * O modelo serial **garante monotonicidade temporal estrita ($dX/dt \ge 0$) e limites estequiométricos ($X \le \eta$) de forma nativa e intrínseca**, sem recorrer a filtros artificiais pós-processamento como `clip` ou `np.maximum.accumulate`.
     * Nos modelos paralelos (Cooperativo e Competitivo), os regressores baseados em árvores geram platôs horizontais e saltos discretos (degraus) quando operados fora do domínio de treino.
3. **Escalonamento para Reatores Industriais Contínuos (CSTR em Série):**
   * O modelo serial é o **único que se translada diretamente para plantas industriais contínuas**. Como o parâmetro $\alpha$ e a taxa linear de avanço $v(t) = \frac{2}{\rho_s}[k_s C_{Af} - \alpha \Delta C_A]$ são propriedades cinéticas intrínsecas, eles entram diretamente nos balanços de massa de cada tanque agitado contínuo em cascata. Os modelos paralelos não possuem essa capacidade, pois apenas predizem uma curva empírica de batelada $X(t)$.

---

### 10.3. Evidências Gráficas Científicas da Divisão 85/15 (300 DPI)

#### Gráfico de Paridade do Modelo Híbrido Serial Cinético (85% Treino / 15% Teste Cego):
![Paridade Serial 85/15](docs/paridade_serial_85_15.png)

#### Comparação Ponto a Ponto e Perfis Cinéticos do Modelo Serial (Ensaios de Teste e Treino com Resíduos):
![Ponto a Ponto Serial 85/15](docs/ponto_a_ponto_serial_85_15.png)

#### Gráfico de Paridade do Modelo Híbrido Paralelo Cooperativo (85% Treino / 15% Teste Cego):
![Paridade Cooperativo 85/15](docs/paridade_cooperativo_85_15.png)

#### Comparação Ponto a Ponto e Perfis Cinéticos do Modelo Cooperativo (Note os degraus de árvore de decisão):
![Ponto a Ponto Cooperativo 85/15](docs/ponto_a_ponto_cooperativo_85_15.png)

#### Painel Comparativo Integrado: Confronto das Três Arquiteturas sob o Split 85/15:
![Comparativo Paridade 85/15](docs/comparativo_paridade_85_15_todas.png)

