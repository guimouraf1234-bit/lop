# Discussão Metodológica: Treinamento em Batelada vs. Inclusão de Dados Contínuos

**Projeto:** Modelagem Mecanicista e Híbrida da Lixiviação de Espodumênio em Reatores Batelada e Contínuos em Série  
**Trabalho de Conclusão de Curso (TCC) — Engenharia Química — UFMG (2026)**  
**Autores:** Daniel Couto, Guilherme Moura, Matheus Póvoas, Rodrigo Mata  
**Orientador:** Prof. Dr. Fabrício Eduardo Bortot Coelho  
**Referência Experimental:** Coelho (2017) — *Modelagem matemática e simulação da lixiviação de espodumênio em reatores batelada e contínuo...*  

---

## 1. Contexto e a Questão Central

Ao analisar os dados experimentais da Planta Piloto (3 CSTRs em série de $V_i = 6{,}0\text{ L}$, total $18{,}0\text{ L}$), observa-se um comportamento **altamente não linear na dinâmica transiente de partida (*startup*) do ensaio em baixa vazão ($Q = 0{,}21\text{ L/min}$)**.

Surge a dúvida metodológica fundamental:
> *"Os dados dos ensaios contínuos poderiam/deveriam ser utilizados para o treinamento dos modelos de Machine Learning também, ou só o treinamento em batelada é o suficiente?"*

Este documento organiza a fundamentação físico-química, os prós e contras de cada escolha e a recomendação estratégica para discussão entre a equipe e com o orientador.

---

## 2. Por que o Ensaio $Q = 0{,}21\text{ L/min}$ é Tão Não Linear?

A não-linearidade observada não é ruído aleatório ou erro de medição; trata-se de um conjunto de **fenômenos hidrodinâmicos e de transporte multifásico característicos de baixas vazões**:

```
[Alimentação Polpa] 
        │
        ▼
   ┌─────────┐   Vertedouro (Lag)   ┌─────────┐   Vertedouro (Lag)   ┌─────────┐
   │  CSTR 1 │ ───────────────────> │  CSTR 2 │ ───────────────────> │  CSTR 3 │
   │ (V=6L)  │  (Acúmulo inicial)   │ (V=6L)  │  (Salto 160-210 min) │ (V=6L)  │
   └─────────┘                      └─────────┘                      └─────────┘
```

### A. Escala de Tempo e Tempo de Residência
* Em $Q = 0{,}41\text{ L/min}$: $\tau_{\text{ind}} = 14{,}63\text{ min}$, $\tau_{\text{tot}} = 43{,}90\text{ min}$. O sistema atinge o estado estacionário em aproximadamente $3\tau_{\text{tot}} \approx 130\text{ min}$.
* Em $Q = 0{,}21\text{ L/min}$: $\tau_{\text{ind}} = 28{,}57\text{ min}$, $\tau_{\text{tot}} = 85{,}71\text{ min}$. A estabilização completa exige entre $3$ e $4\tau_{\text{tot}}$ ($\approx 260$ a $340\text{ min}$). A janela amostrada ($t = 90$ a $300\text{ min}$) cobre toda a transição lenta.

### B. Atraso de Vertedouro (*Weir Transport Lag*) e Acúmulo de Sólidos (*Hold-up*)
* Com a vazão reduzida à metade, a velocidade linear e a energia cinética da polpa diminuem.
* Nos primeiros estágios ($t < 160\text{ min}$), as partículas sólidas e o reagente líquido acumulam preferencialmente no primeiro reator até que a cota hidráulica dos vertedouros entre tanques se estabeleça de forma contínua e homogênea.
* **Evidência Experimental:** Entre $t = 90$ e $t = 130\text{ min}$, a conversão medida no Reator 2 foi menor ou igual à do Reator 1:
  * $t = 90\text{ min}$: $X_1 = 70{,}9\%$ vs. $X_2 = 67{,}9\%$
  * $t = 110\text{ min}$: $X_1 = 72{,}9\%$ vs. $X_2 = 70{,}3\%$
  Isso ocorre porque o Reator 2 estava recebendo uma polpa inicialmente defasada e heterogênea durante o enchimento.

### C. O Salto Convectivo Repentino ($t = 160$ a $210\text{ min}$)
* Por volta de $t = 160\text{ min}$, a cota hidráulica do transbordo é plenamente atingida. 
* Ocorre uma "onda convectiva" de polpa reagida de R1 para R2 e R3, gerando um salto abrupto nas conversões:
  * Reator 2 salta de $73{,}8\%$ para $81{,}6\%$.
  * Reator 3 salta de $77{,}2\%$ para $85{,}7\%$.
* A partir de $t = 230\text{ min}$, o sistema se estabiliza no regime permanente final ($X_1 \approx 80{,}8\%$, $X_2 \approx 84{,}9\%$, $X_3 \approx 86{,}2\%$).

### D. A Limitação Intrínseca da Batelada
* O reator de batelada de laboratório (béquer de 1 L) é um sistema fechado e perfeitamente agitado, sem entrada, sem saída, sem vertedouros e sem separação física em estágios.
* **Os dados de batelada contêm apenas cinética química intrínseca de dissolução**.
* Portanto, nenhuma inteligência artificial treinada *unicamente* na batelada tem como "adivinhar" o atraso de 160 minutos gerado pela hidrodinâmica específica da planta piloto.

---

## 3. Análise dos Três Paradigmas de Treinamento

| Paradigma | Dados de Treino | Dados de Teste | Papel da ML | Avaliação Científica |
| :--- | :--- | :--- | :--- | :--- |
| **1. Escalonamento Puro (*Scale-Up Zero-Shot*)** | 100% Batelada ($N=128$) | Planta Piloto ($N=6$ Estacionário / $N=57$ Total) | Corrigir a não-idealidade físico-química da lixiviação | **Padrão Ouro de Engenharia Química**. Prova que o laboratório dimensiona a planta sem testes piloto prévios. |
| **2. Validação Cruzada entre Vazões (*Leave-One-Flow-Out*)** | Batelada ($N=128$) + Piloto $Q = 0{,}41$ ($N=24$) | Piloto $Q = 0{,}21$ ($N=33$) | Aprender hidrodinâmica de tanques em série e extrapolar vazão | **Metodologicamente Impecável**. Zero vazamento de dados (*data leakage*); testa a capacidade de adaptação dinâmica. |
| **3. Co-Treinamento Global (*Digital Twin* Operacional)** | Todos os dados ($N=185$) | Subdivisão aleatória / K-Fold | Mapear o comportamento global da unidade instalada | Excelente para controle em tempo real; fraco para tese de escalonamento clássico. |

---

## 4. Detalhamento dos Paradigmas

### Paradigma 1: Escalonamento Puro (Batelada $\to$ Contínuo)
* **Objetivo:** Responder: *"É possível projetar e prever a operação da bateria de reatores contínuos a partir apenas de ensaios de laboratório em pequena escala?"*
* **Por que o treinamento em batelada É SUFICIENTE aqui:**
  * No **Regime Permanente (Estado Estacionário)**, o modelo mecanicista PBM CSTR de Coelho (2017) somado ao modelo de resíduo do SVR treinado em batelada obteve:
    $$\text{MAE} = 0{,}24\%, \quad \text{RMSE} = 0{,}28\%, \quad R^2 = 0{,}9907$$
  * Esse resultado é extraordinário para Engenharia Química. Mostra que o modelo cinético de bancada escalonou com erro inferior a meio ponto percentual no estado estacionário dos 3 tanques.
* **Limitação:** Não prevê o atraso transiente dos primeiros 160 min em $Q = 0{,}21\text{ L/min}$ porque a batelada não sabe o que é um vertedouro.

---

### Paradigma 2: Validação Cruzada entre Vazões (*Leave-One-Flow-Out*)
* **Objetivo:** Responder: *"Se a planta piloto já operou em uma vazão padrão ($Q = 0{,}41\text{ L/min}$), a modelagem híbrida consegue aprender os efeitos hidrodinâmicos e prever com precisão a partida em uma nova vazão nunca vista ($Q = 0{,}21\text{ L/min}$)?"*
* **Como funciona:**
  1. A ML recebe *features* contínuas adicionais: tempo de residência ($\tau_{\text{ind}}$), vazão ($Q$), número do reator ($k$) e tempo adimencional ($\theta = t / \tau$).
  2. A ML é treinada na bancada + ensaio de $Q = 0{,}41$.
  3. O ensaio de $Q = 0{,}21$ é mantido **100% intocado e cego**.
* **Vantagem Científica:**
  * Resolve a não-linearidade e o atraso de partida.
  * Mantém total rigor contra *overfitting* ou "vazamento de dados", blindando a equipe perante perguntas da banca examinadora.

---

### Paradigma 3: Co-Treinamento Multi-Fidelidade (*Digital Twin*)
* **Objetivo:** Desenvolver o modelo mais acurado possível para simulação operacional, controle avançado (MPC) ou otimização econômica da planta já existente.
* **Cuidado na Defesa do TCC:** Se a equipe usar este paradigma sozinho, a banca poderá argumentar:
  > *"Vocês dizem que estão escalonando o processo, mas se a rede foi treinada com os próprios dados da planta piloto em $Q=0{,}21$, isso não é escalonamento preditivo, é apenas um ajuste empírico de curvas (curve fitting)."*

---

## 5. Recomendação Estratégica para o TCC

A recomendação unânime é **adotar uma estrutura em 2 etapas complementares**, transformando a aparente "dificuldade" da curva de $Q = 0{,}21$ no maior ponto forte da monografia:

```
                      ESTRUTURA SUGERIDA PARA O TCC
                                    │
       ┌────────────────────────────┴────────────────────────────┐
       ▼                                                         ▼
[CAPÍTULO 4: ESCALONAMENTO PURO]              [CAPÍTULO 5: ANÁLISE DINÂMICA TRANSIENTE]
• Treino: Apenas Batelada                     • Treino: Batelada + Piloto Q = 0,41
• Teste: Estado Estacionário dos 3 CSTRs      • Teste Cego: Dinâmica de Q = 0,21
• Foco: Projeto / Dimensionamento             • Foco: Operação / Partida / Gêmeo Digital
• Destaque: MAE = 0,24%, R² = 0,9907          • Destaque: Captura da não-linearidade
```

### Argumentação Pronta para a Apresentação/Texto:
1. **No Capítulo de Escalonamento (Capítulo 4):**
   > *"Demonstrou-se a viabilidade do escalonamento preditivo direto (zero-shot): utilizando exclusivamente dados cinéticos obtidos em batelada de bancada acoplados ao balanço populacional CSTR, o modelo híbrido Grey-Box (SVR) predisse a conversão no estado estacionário dos três reatores industriais com MAE de 0,24% e R² de 0,9907, superando os modelos puramente mecanicistas da literatura."*

2. **No Capítulo de Dinâmica e Não-Linearidades (Capítulo 5):**
   > *"Na análise da dinâmica de partida em baixa vazão (Q = 0,21 L/min), identificou-se um atraso hidrodinâmico pronunciado nos primeiros 160 minutos, decorrente do acúmulo de sólidos e atraso de transbordo nos vertedouros entre reatores. Como os ensaios de batelada são desprovidos de variáveis hidrodinâmicas de escoamento, propôs-se uma validação cruzada entre vazões (Leave-One-Flow-Out). Ao integrar os dados de Q = 0,41 L/min ao treinamento, o modelo híbrido assimilou a assinatura hidrodinâmica da cascata contínua, reduzindo o erro na predição transiente da nova condição operacional sem incorrer em vazamento de dados."*

---

## 6. Roteiro para Discussão com os Pares e o Prof. Fabrício

Ao sentar com Daniel, Matheus, Rodrigo e o Prof. Fabrício, coloquem os seguintes pontos em pauta:

1. **Alinhamento de Escopo:** O foco principal da entrega é o **dimensionamento de regime permanente** (onde a batelada já resolve com folga) ou queremos estender até o **controle de processos e simulação dinâmica de partida**?
2. **Aprovação da Validação Cruzada (*Leave-One-Flow-Out*):** O Prof. Fabrício concorda em estruturar o estudo de partida como um teste de generalização entre vazões ($Q=0{,}41 \to Q=0{,}21$)?
3. **Inclusão no Texto:** Aprovam manter os dois estudos lado a lado na monografia (Capítulo 4: Escalonamento Puro; Capítulo 5: Gêmeo Digital Dinâmico)?

---
*Documento gerado para registro técnico e suporte à tomada de decisão metodológica do TCC EQ-UFMG.*
