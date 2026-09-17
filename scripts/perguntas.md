# Perguntas Frequentes & Respostas Técnicas para o TCC

## 1. E se a granulometria do minério for diferente?

### Resposta Rápida:
O modelo absorve a mudança de granulometria **diretamente pela física do Balanço Populacional (PBM)**, sem precisar retreinar o algoritmo de Machine Learning.

### Detalhamento Técnico:
1. **Onde a granulometria entra no modelo:**  
   A granulometria da calcina é definida pelos parâmetros da distribuição de Rosin-Rammler-Bennett (RRB):
   * $D_{63{,}2}$: diâmetro característico de corte em massa ($63{,}2\%$).
   * $m$: coeficiente de uniformidade da dispersão granulométrica.
2. **O que acontece se o lote for mais fino ($D_{63{,}2} < 41{,}65\ \mu\text{m}$):**  
   * A classe `PopulationBalanceModel(m, d632)` recalcula analiticamente o 3º momento inicial $I_0 = D_{63{,}2}^3 \cdot \Gamma(1 + 3/m)$.
   * Como a área superficial específica aumenta proporcionalmente a $1/D$, o encolhimento de diâmetro consome as partículas menores muito mais rápido nos primeiros 2 minutos.
   * A conversão $X_{\text{PBM}}(t)$ sobe com maior inclinação inicial, respeitando a estequiometria do ácido livre.
3. **Por que isso é um diferencial do Modelo Híbrido (Grey-Box)?**  
   * Uma Rede Neural pura (Black-Box) erraria completamente ao receber um minério com tamanho diferente, pois ela não sabe calcular integrais de momento nem conhece conservação de massa.
   * No modelo Híbrido, **a física carrega $90\%$ do efeito da nova granulometria**, e o Machine Learning apenas ajusta o pequeno resíduo $\Delta X$.
