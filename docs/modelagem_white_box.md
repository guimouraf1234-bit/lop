# Documentação Teórica e Equacionamento do Modelo Mecanicista (White-Box) para a Lixiviação de Concentrado Ustulado de Zinco

**Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG**  
**Discentes:** Daniel Couto Vieira, Guilherme Moura de Sousa Franco, Matheus Henrique Borba Póvoas, Rodrigo Amaral da Mata  
**Orientador:** Prof. Dr. Fabrício Eduardo Bortot Coelho  
*Departamento de Engenharia Química -- Universidade Federal de Minas Gerais (UFMG)*  
*Setembro de 2026*

---

## Resumo
Este documento apresenta a fundamentação teórica rigorosa, as hipóteses simplificadoras e o equacionamento matemático do Modelo Mecanicista Baseado em Primeiros Princípios (*White-Box*) para o processo de lixiviação ácida atmosférica de concentrado ustulado de zinco (calcina) em reator batelada. O modelo acopla o Balanço Populacional unidimensional com taxa de retração governada pelo Método das Características à distribuição granulométrica inicial de Rosin-Rammler-Bennett (RRB), considerando a cinética química interfacial heterogênea de dissolução da zincita ($\text{ZnO}$) e a estequiometria do consumo de ácido sulfúrico livre. A formulação descrita estabelece a base física do modelo híbrido cinza (*Grey-Box*), fornecendo a referência mecanicista sobre a qual o aprendizado de máquina atua na predição de desvios e fenômenos secundários não modelados.

---

## 1. Introdução e Hipóteses do Modelo

A lixiviação atmosférica de concentrado ustulado de zinco (calcina) em meio de ácido sulfúrico diluído constitui a etapa primordial da rota hidrometalúrgica de recuperação do zinco metálico. O concentrado ustulado é composto majoritariamente por zincita ($\text{ZnO}$), mineral prontamente solúvel, e ferrita de zinco ($\text{ZnFe}_2\text{O}_4$), fase refratária de dissolução lenta em condições brandas.

Para o desenvolvimento do modelo de primeiros princípios (*White-Box*), estabeleceram-se as seguintes premissas físico-químicas e hidrodinâmicas:
1. **Regime Hidrodinâmico:** O reator batelada opera sob agitação mecânica intensa ($\ge 1000\text{ rpm}$), garantindo que a resistência difusional na camada limite fluida ao redor da partícula seja desprezível ($Sh \to \infty$).
2. **Mecanismo Controlador:** A dissolução do sólido é controlada pela reação química superficial na interface sólido-líquido, descrita pelo Modelo de Núcleo Não Reagido (*Shrinking Core Model* -- SCM) com encolhimento contínuo do raio da partícula.
3. **Morfologia das Partículas:** Assume-se geometria esférica isotrópica com fator de forma volumétrico constante ($\psi = 1$) para todas as frações granulométricas.
4. **Fenômenos Populacionais:** Ausência de quebra mecânica de partículas e ausência de aglomeração/agregação no meio reacional ($B = D = 0$).
5. **Isotermicidade:** O sistema é mantido isotérmico na temperatura controlada do banho termostático ($T = 40^\circ\text{C} = 313{,}15\text{ K}$).
6. **Seletividade Mineralógica:** Em lixiviação neutra/branda ($T \le 40^\circ\text{C}$), apenas a zincita solúvel ($\text{ZnO}$) participa ativamente da reação no período monitorado ($t \le 15\text{ min}$). A fração de zinco associada à ferrita ($13\%$) permanece inerte.

---

## 2. Estequiometria e Balanço de Massa no Reator Batelada

A reação de dissolução da zincita por ácido sulfúrico aquoso obedece à estequiometria mássica e molar unitária (1:1):

$$\text{ZnO}_{(s)} + \text{H}_2\text{SO}_{4(aq)} \longrightarrow \text{ZnSO}_{4(aq)} + \text{H}_2\text{O}_{(l)}$$

Sejam $n_{A0}$ o número inicial de moles de reagente fluido ($\text{H}_2\text{SO}_4$) e $n_{B0}$ o número inicial de moles de reagente sólido solúvel ($\text{ZnO}$) carregados no reator de volume líquido $V$:

$$n_{A0} = C_{A0} \, V$$

$$n_{B0} = \frac{m_{B0} \, T_{\text{ZnO}}}{M_{\text{ZnO}}}$$

onde $C_{A0}$ é a concentração inicial de ácido ($\text{mol}\cdot\text{L}^{-1}$), $m_{B0}$ é a massa total de concentrado ustulado alimentada ($\text{g}$), $T_{\text{ZnO}}$ é o teor em massa de zincita na amostra ($76{,}1\% = 0{,}761$), e $M_{\text{ZnO}} = 81{,}38\text{ g}\cdot\text{mol}^{-1}$ é a massa molar da zincita.

A razão molar estequiométrica $\eta$ entre o ácido sulfúrico e a zincita é definida por:

$$\eta = \frac{n_{A0}}{n_{B0}} = \frac{C_{A0} \, V \, M_{\text{ZnO}}}{m_{B0} \, T_{\text{ZnO}}}$$

Pelo balanço molar de conservação de espécies em sistema fechado, a concentração instantânea de ácido sulfúrico livre em solução, $C_{Af}(t)$, relaciona-se diretamente com a conversão fracionária da zincita, $X(t)$:

$$C_{Af}(t) = \frac{n_{A0} - n_{B0} \, X(t)}{V} = C_{A0} \left( 1 - \frac{n_{B0}}{n_{A0}} X(t) \right)$$

Substituindo a razão estequiométrica $\eta$, obtém-se a lei de consumo de ácido:

$$C_{Af}(t) = C_{A0} \left( 1 - \frac{X(t)}{\eta} \right)$$

### Análise dos Três Regimes Operacionais:
* **Deficiência de Ácido ($\eta < 1{,}0$):** O ácido livre esgota-se ($C_{Af} \to 0$) antes da solubilização completa do zinco. A reação cessa no patamar estequiométrico máximo teórico:
  $$X_{\text{máx}} = \eta \quad (\text{para } \eta < 1{,}0)$$
  *(Para $\eta = 0{,}5$, a conversão máxima estabiliza-se rigorosamente em $X \approx 0{,}50$)*.
* **Proporção Estequiométrica ($\eta = 1{,}0$):** Ácido e sólido em proporções equivalentes. Conforme $X \to 1$, $C_{Af} \to 0$, gerando acentuada desaceleração cinética terminal.
* **Excesso de Ácido ($\eta > 1{,}0$):** Permanece ácido residual no licor final ($C_{Af} > 0$), permitindo que a conversão da zincita atinja $X \to 1{,}0$ ($100\%$).

---

## 3. Cinética Química Heterogênea e Taxa de Retração

Para uma partícula de diâmetro esférico $D$, o número de moles de zincita remanescente é:

$$N_B = \rho_s \, \left( \frac{\pi}{6} D^3 \right)$$

Diferenciando em relação ao tempo:

$$\frac{dN_B}{dt} = \frac{\pi}{2} \rho_s \, D^2 \frac{dD}{dt}$$

A taxa intrínseca de reação heterogênea na interface sólido-líquido de área superficial $A_p = \pi D^2$ é de primeira ordem em relação ao ácido livre:

$$\left( -\frac{1}{A_p} \frac{dN_B}{dt} \right) = k_s \, C_{Af}(t)$$

onde $k_s$ é a constante cinética superficial ($\mu\text{m}\cdot\text{min}^{-1}$) e $\rho_s = 69{,}2\text{ mol}\cdot\text{L}^{-1}$ de sólido. Igualando as taxas, deduz-se a velocidade linear de retração do diâmetro da partícula:

$$v(t) \equiv -\frac{dD}{dt} = \frac{2 \, k_s}{\rho_s} \, C_{Af}(t)$$

### Dependência Térmica (Arrhenius)
$$k_s(T) = k_{s,\text{ref}} \exp \left[ -\frac{E_a}{R} \left( \frac{1}{T} - \frac{1}{T_{\text{ref}}} \right) \right]$$

onde $E_a = 13{,}45\text{ kJ}\cdot\text{mol}^{-1}$ (Balarini et al., 2025), $T_{\text{ref}} = 313{,}15\text{ K}$ ($40^\circ\text{C}$), e $k_{s,\text{ref}} = 1{,}8 \times 10^4\ \mu\text{m}\cdot\text{min}^{-1}$.

### Termo de Amortecimento Empírico (Tese Bortot Coelho, 2017)
$$v(t) = \frac{2}{\rho_s} \max \Big( k_s \, C_{Af}(t) - \alpha \, [C_{A0} - C_{Af}(t)],\ 0 \Big)$$
com $\alpha = 5{,}5 \times 10^3\ \mu\text{m}\cdot\text{min}^{-1}$.

---

## 4. Balanço Populacional e Método das Características

A equação diferencial parcial populacional para a densidade $\psi(D,t)$ em sistema perfeitamente misturado é:

$$\frac{\partial \psi(D,t)}{\partial t} + \frac{\partial \big[ -v(t) \, \psi(D,t) \big]}{\partial D} = 0$$

Pelo Método das Características, como $v(t)$ não depende de $D$:

$$\frac{dD}{dt} = -v(t) \implies D(t) = \max \big( D_0 - \Delta D(t),\ 0 \big)$$

onde o encolhimento acumulado do diâmetro no intervalo $[0, t]$ é:

$$\Delta D(t) = \int_0^t v(t') \, dt'$$

---

## 5. Distribuição Granulométrica de Rosin-Rammler-Bennett (RRB)

A distribuição granulométrica cumulativa $F_0(D)$ e a densidade de probabilidade $f_0(D)$ inicial da calcina de Três Marias são dadas por:

$$F_0(D) = 1 - \exp \left[ -\left( \frac{D}{D_{63{,}2}} \right)^m \right]$$

$$f_0(D) = \frac{m}{D_{63{,}2}} \left( \frac{D}{D_{63{,}2}} \right)^{m-1} \exp \left[ -\left( \frac{D}{D_{63{,}2}} \right)^m \right]$$

com $m = 1{,}022$ e $D_{63{,}2} = 41{,}65\ \mu\text{m}$.

---

## 6. Conversão Fracionária ($X_{\text{PBM}}$) pelo 3º Momento Volumétrico

A conversão volumétrica/mássica da zincita é calculada pela razão entre o 3º momento volumétrico residual e o inicial:

$$X_{\text{PBM}}(t) = 1 - \frac{\mu_3(t)}{\mu_3(0)} = 1 - \frac{\displaystyle \int_{\Delta D(t)}^\infty \big[ D - \Delta D(t) \big]^3 f_0(D) \, dD}{\displaystyle \int_0^\infty D^3 f_0(D) \, dD}$$

### Solução Analítica Exata do 3º Momento Inicial ($I_0$)
Para a distribuição RRB/Weibull:

$$I_0 \equiv \int_0^\infty D^3 f_0(D) \, dD = D_{63{,}2}^3 \, \Gamma \left( 1 + \frac{3}{m} \right)$$

Para $m = 1{,}022$ e $D_{63{,}2} = 41{,}65\ \mu\text{m}$:
$$\Gamma \left( 1 + \frac{3}{1{,}022} \right) \approx 5{,}6645 \implies I_0 \approx 4{,}094 \times 10^5\ \mu\text{m}^3$$

### Sistema Dinâmico Acoplado (Integrado via RK45)
$$\begin{cases}
\displaystyle \frac{d(\Delta D)}{dt} = v \big( C_{Af}(t) \big) \\[1em]
\displaystyle X_{\text{PBM}}(t) = 1 - \frac{1}{I_0} \int_{\Delta D(t)}^{D_{\text{máx}}} \big[ D - \Delta D(t) \big]^3 f_0(D) \, dD \\[1em]
\displaystyle C_{Af}(t) = C_{A0} \left( 1 - \frac{X_{\text{PBM}}(t)}{\eta} \right)
\end{cases}$$

com $\Delta D(0) = 0{,}0\ \mu\text{m}$, $X_{\text{PBM}}(0) = 0{,}0$ e $C_{Af}(0) = C_{A0}$.

---

## 7. Tabela de Parâmetros Físico-Químicos

| Símbolo | Descrição | Valor | Unidade | Fonte |
| :--- | :--- | :---: | :---: | :--- |
| $k_{s,\text{ref}}$ | Constante cinética superficial a $40^\circ\text{C}$ | $1{,}8 \times 10^4$ | $\mu\text{m}\cdot\text{min}^{-1}$ | Balarini et al. (2025) |
| $E_a$ | Energia de ativação aparente | $13{,}45$ | $\text{kJ}\cdot\text{mol}^{-1}$ | Balarini et al. (2025) |
| $\rho_s$ | Densidade molar da fase sólida | $69{,}2$ | $\text{mol}\cdot\text{L}^{-1}$ | Bortot Coelho (2017) |
| $m$ | Índice de uniformidade RRB | $1{,}022$ | adim. | Bortot Coelho (2017) |
| $D_{63{,}2}$ | Diâmetro característico RRB | $41{,}65$ | $\mu\text{m}$ | Bortot Coelho (2017) |
| $T_{\text{ZnO}}$ | Teor de zincita na calcina bruta | $76{,}1$ | $\% \text{ m/m}$ | Bortot Coelho (2017) |
| $f_{\text{ZnO}}$ | Fração de zinco solúvel (zincita) | $0{,}87$ | adim. | Bortot Coelho (2017) |
| $M_{\text{ZnO}}$ | Massa molar da zincita | $81{,}38$ | $\text{g}\cdot\text{mol}^{-1}$ | Lide (2005) |
| $M_{\text{H}_2\text{SO}_4}$ | Massa molar do ácido sulfúrico | $98{,}08$ | $\text{g}\cdot\text{mol}^{-1}$ | Lide (2005) |
| $V$ | Volume reacional de bancada | $0{,}40$ | $\text{L}$ | Bortot Coelho (2017) |
| $T$ | Temperatura de operação | $40{,}0$ | $^\circ\text{C}$ | Bortot Coelho (2017) |
| $\alpha$ | Amortecimento empírico da tese | $5{,}5 \times 10^3$ | $\mu\text{m}\cdot\text{min}^{-1}$ | Bortot Coelho (2017) |

---

## 8. Justificativa para a Modelagem Híbrida Grey-Box

Embora o modelo *White-Box* apresente sólida consistência física ($R^2 = 0{,}9253$), a existência do resíduo sistemático $\Delta X = X_{\text{exp}} - X_{\text{PBM}}$ decorre de fenômenos não modelados:
1. **Efeito de Íon Comum:** O acúmulo de $\text{Zn}^{2+}$ e $\text{SO}_4^{2-}$ reduz o coeficiente de atividade do ácido livre $\gamma_{\text{H}^+}$, diminuindo a velocidade real da reação.
2. **Passivação Difusional por Sílica:** Silicatos solúveis na calcina precipitam um gel de sílica amorfa que oclui microporos das partículas.
3. **Morfologia dos Finos:** Partículas ultrafinas com cantos vivos e rugosidade dissolvem mais rápido no início ($t \le 1\text{ min}$) do que esferas lisas ideais.

A modelagem **Grey-Box** combina a confiabilidade estequiométrica do PBM à flexibilidade do Machine Learning, que prevê esse resíduo e eleva o coeficiente de determinação para $R^2 = 0{,}9979$.

---

## Referências Bibliográficas
1. BALARINI, J. C.; ARAÚJO, E. M. R.; COELHO, F. E. B.; POLLI, L. O.; KONZEN, C.; MIRANDA, T. L. S.; SALUM, A. Leaching kinetics of roasted zinc concentrate in sulphuric acid solutions. *Revista Observatorio de la Economía Latinoamericana*, v. 23, n. 11, 2025. DOI: 10.55905/oelv23n11-188.
2. BORTOT COELHO, F. E. *Desenvolvimento e Validação de um Modelo de Balanço Populacional para a Lixiviação de um Concentrado Ustulado de Zinco em Bancada e em Escala Piloto*. 2017. 248 f. Dissertação (Mestrado em Engenharia Química) -- Universidade Federal de Minas Gerais, Belo Horizonte, 2017.
3. CRUNDWELL, F. K.; BRYSON, A. W. The calculation of the rate of dissolution of minerals in batch and continuous reactors. *Hydrometallurgy*, v. 29, n. 1-3, p. 275-295, 1992.
4. FOGLER, H. S. *Elements of Chemical Reaction Engineering*. 4. ed. Boston: Prentice Hall, 2006.
5. LEVENSPIEL, O. *Chemical Reaction Engineering*. 3. ed. New York: John Wiley & Sons, 1999.
6. LIDE, D. R. *CRC Handbook of Chemistry and Physics*. 86. ed. Boca Raton: CRC Press, 2005.
