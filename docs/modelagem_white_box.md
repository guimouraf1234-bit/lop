# Documentação Teórica e Equacionamento do Modelo Mecanicista (White-Box) para a Lixiviação de Concentrado Ustulado de Zinco

**Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG**  
**Discentes:** Daniel Couto Vieira, Guilherme Moura de Sousa Franco, Matheus Henrique Borba Póvoas, Rodrigo Amaral da Mata  
**Orientador:** Prof. Dr. Fabrício Eduardo Bortot Coelho  
*Departamento de Engenharia Química -- Universidade Federal de Minas Gerais (UFMG)*  
*Setembro de 2026*

---

## Resumo
Este documento apresenta as hipóteses e o equacionamento do Modelo Mecanicista (*White-Box*) para a lixiviação ácida atmosférica de concentrado ustulado de zinco (calcina) em reator batelada. O ponto de partida é o modelo de **Bortot Coelho (2017)** (item 4.3, Tabela 5.9 e Apêndice A2.1 da dissertação): Balanço Populacional unidimensional resolvido pelo Método das Características, distribuição granulométrica de Rosin-Rammler-Bennett (RRB), cinética de reação química interfacial (Eq. 3.42) e consumo estequiométrico de ácido sulfúrico (Eq. 3.51). Sobre ele aplicam-se três correções de consistência física: (i) a RRB, ajustada a frações mássicas, é integrada em base mássica; (ii) a constante $k_s$ é referida à temperatura em que foi estimada (30 °C); (iii) a força motriz é a concentração de $\text{H}^+$, a mesma base em que $k_s$ foi estimada, calculada pela especiação $\text{HSO}_4^-/\text{SO}_4^{2-}$. O parâmetro empírico $\alpha$ da Eq. 4.1 é fixado em zero: o desvio que ele corrigia passa a ser aprendido pelo modelo de Machine Learning do híbrido *Grey-Box*. A formulação original da tese continua disponível no código (`formulacao='tese'`).

---

## 1. Introdução e Hipóteses do Modelo

A lixiviação atmosférica de concentrado ustulado de zinco (calcina) em meio de ácido sulfúrico diluído constitui a etapa primordial da rota hidrometalúrgica de recuperação do zinco metálico. O concentrado ustulado é composto majoritariamente por zincita ($\text{ZnO}$, 76,1% m/m), mineral prontamente solúvel, e por ferrita de zinco ($\text{ZnFe}_2\text{O}_4$), willemita e esfalerita, fases de dissolução lenta em condições brandas.

Premissas do modelo (Bortot Coelho, 2017, item 4.3):
1. **Regime Hidrodinâmico:** agitação intensa ($1000\text{ rpm}$), com resistência difusional na camada limite desprezível.
2. **Mecanismo Controlador:** reação química superficial de primeira ordem em $[\text{H}^+]$ (Balarini, 2009), com encolhimento contínuo das partículas.
3. **Morfologia das Partículas:** partículas esféricas.
4. **Fenômenos Populacionais:** sem quebra nem aglomeração ($B = D = 0$).
5. **Isotermicidade:** $T = 40^\circ\text{C}$.
6. **Seletividade Mineralógica:** apenas a reação da zincita é considerada (Eq. 4.2), tanto na cinética quanto no consumo de ácido. A conversão experimental é normalizada pelo zinco da zincita, descontando o zinco da ferrita (Eq. 4.5).

---

## 2. Estequiometria e Balanço de Massa no Reator Batelada

$$\text{ZnO}_{(s)} + \text{H}_2\text{SO}_{4(aq)} \longrightarrow \text{ZnSO}_{4(aq)} + \text{H}_2\text{O}_{(l)} \qquad (4.2)$$

A razão molar é definida sobre a zincita (Eq. 4.3):

$$\eta = \frac{n_{A0}}{n_{B0}} = \frac{V \, C_{A0}}{m_{B0} \, T_{\text{ZnO}} / M_{\text{ZnO}}} \qquad (4.3)$$

com $T_{\text{ZnO}} = 0{,}761$ e $M_{\text{ZnO}} = 81{,}38\text{ g}\cdot\text{mol}^{-1}$. O modelo usa o $\eta$ calculado com as massas de cada ensaio (coluna `razao_molar_eta_real`): $0{,}535 / 1{,}069 / 1{,}58\text{–}1{,}61 / 3{,}19\text{–}3{,}29$ para os níveis nominais $0{,}5 / 1{,}0 / 1{,}5 / 3{,}1$ do planejamento. É com esses valores que se reproduz a soma quadrática de $\approx 0{,}05$ relatada no ajuste de $\alpha$ da tese (item 5.2.4).

A conversão experimental (Eq. 4.5) também está na base da zincita:

$$X_{\text{Zn}} = \frac{D \, V \, (C_{\text{Zn},f} - C_{\text{Fe},f}/2)}{m_{B0} \, T_{\text{ZnO}}} \qquad (4.5)$$

Balanço de ácido em sistema fechado (Eq. 3.51, Herbst, 1979):

$$C_{Af}(t) = C_{A0} \left( 1 - \frac{X(t)}{\eta} \right) \qquad (3.51)$$

> O consumo de ácido medido é de $1{,}07$ a $1{,}13$ mol por mol de ZnO dissolvido, e não $1{,}00$: outros óxidos do concentrado também consomem ácido. Por coerência com a hipótese 6, esse consumo não entra no modelo; é a causa do desvio de $-0{,}04$ no patamar com $\eta = 0{,}5$ (Seção 9).

---

## 3. Cinética Química e Taxa de Variação do Diâmetro

Para controle por reação química, a taxa de variação do diâmetro (Eq. 3.42 da tese, LeBlanc e Fogler, 1987) é:

$$v \equiv -\frac{dD}{dt} = \frac{2 \, k_s}{\rho_s} \, C_A \qquad (3.42)$$

com $\rho_s = 69{,}2\ \text{mol}\cdot\text{L}^{-1}$.

### 3.1. Base e temperatura de $k_s$
Balarini (2009) estimou $k_s = 3{,}00 \times 10^{-2}\ \text{cm}\cdot\text{s}^{-1} = 1{,}8 \times 10^4\ \mu\text{m}\cdot\text{min}^{-1}$ usando, na Eq. 3.42, a **concentração hidrogeniônica** $C_A = [\text{H}^+]$ (item 3.3.2 da tese), a **$30^\circ\text{C}$** (Balarini et al., 2025, Fig. 8: $T = 30^\circ\text{C}$, $6\ \text{g/L}$ de $\text{H}_2\text{SO}_4$, $[\text{H}^+] = 0{,}0689\ \text{mol/L}$). A tese aplica esse valor a $40^\circ\text{C}$ e sobre a molaridade de $\text{H}_2\text{SO}_4$. Para manter a coerência com a estimação:

$$k_s(T) = k_{s,30^\circ\text{C}} \exp \left[ -\frac{E_a}{R} \left( \frac{1}{T} - \frac{1}{303{,}15} \right) \right], \qquad E_a = 13{,}45\ \text{kJ}\cdot\text{mol}^{-1}$$

o que dá $k_s(40^\circ\text{C}) = 2{,}13 \times 10^4\ \mu\text{m}\cdot\text{min}^{-1}$, e a força motriz é $C_A = [\text{H}^+]$.

### 3.2. Concentração de $\text{H}^+$ por especiação
O ácido sulfúrico tem a primeira dissociação completa e a segunda parcial:

$$\text{HSO}_4^- \rightleftharpoons \text{H}^+ + \text{SO}_4^{2-}, \qquad K_2 = \frac{[\text{H}^+][\text{SO}_4^{2-}]}{[\text{HSO}_4^-]} = 0{,}0099\ \text{mol/L}$$

Cada mol de ZnO dissolvido consome 2 mol de $\text{H}^+$, mas o sulfato permanece em solução como $\text{ZnSO}_4$. Assim, com prótons totais $H_T = 2\,C_{Af}$ e sulfato total $S_T = C_{A0}$:

$$[\text{H}^+] + [\text{HSO}_4^-] = 2\,C_{Af}, \qquad [\text{HSO}_4^-] = \frac{[\text{H}^+]\,\big(C_{A0} - [\text{HSO}_4^-]\big)}{K_2}$$

$$[\text{HSO}_4^-] = \frac{(S_T + H_T + K_2) - \sqrt{(S_T + H_T + K_2)^2 - 4\,S_T H_T}}{2}, \qquad [\text{H}^+] = H_T - [\text{HSO}_4^-]$$

O valor $K_2 = 0{,}0099$ não é ajustado: é o que reproduz o $[\text{H}^+] = 0{,}0689\ \text{mol/L}$ usado por Balarini para $6\ \text{g/L}$ de $\text{H}_2\text{SO}_4$ ($0{,}0612\ \text{mol/L}$). Fisicamente, o sulfato liberado pela dissolução tampona o $\text{H}^+$ e a reação desacelera mesmo com ácido livre disponível. A especiação usa concentrações (sem coeficientes de atividade), uma aproximação nas forças iônicas mais altas ($C_{A0} = 1{,}5\ \text{mol/L}$).

### 3.3. Parâmetro $\alpha$ da tese
Bortot Coelho (2017) adicionou um parâmetro ajustável à Eq. 3.42 para corrigir a conversão final (Eq. 4.1):

$$v = \frac{2}{\rho_s} \Big[ k_s \, C_{Af} - \alpha \, (C_{A0} - C_{Af}) \Big] \qquad (4.1)$$

com $\alpha = 5{,}5 \times 10^3\ \mu\text{m}\cdot\text{min}^{-1}$. **Neste trabalho, $\alpha = 0$**: o termo é um ajuste empírico das conversões finais e, quando o ácido é limitante, impede o seu esgotamento (Seção 8). No modelo híbrido, esse papel é assumido pelo regressor de Machine Learning. O modelo da tese com $\alpha = 5500$ é mantido nos gráficos como referência ("Tese 2017").

---

## 4. Balanço Populacional e Método das Características

Equação do balanço populacional em batelada (Eq. 3.38):

$$\frac{\partial \psi(D,t)}{\partial t} = -\frac{\partial}{\partial D} \big[ -v \, \psi(D,t) \big] \qquad (3.38)$$

Como $v$ não depende de $D$, a solução pelo Método das Características é um deslocamento da distribuição inicial:

$$\psi(D,t) = \psi\big(D + \Delta D(t),\ 0\big), \qquad \Delta D(t) = \int_0^t v(t') \, dt'$$

A tese resolve a EDP pelo Método Numérico das Linhas (`NDSolve` do Mathematica, Apêndice A2.1). A solução pelas características é a mesma solução exata, integrada aqui com RK45 (`scipy.integrate.solve_ivp`).

---

## 5. Distribuição Granulométrica Inicial (RRB)

A RRB foi ajustada às frações **mássicas** passantes do peneiramento e da difração a laser (Tabela 5.4 da tese). Sua derivada é, portanto, a densidade mássica $q_3(D)$:

$$q_3(D) = \frac{m}{D_{63{,}2}} \left( \frac{D}{D_{63{,}2}} \right)^{m-1} \exp \left[ -\left( \frac{D}{D_{63{,}2}} \right)^m \right]$$

com $m = 1{,}022$ e $D_{63{,}2} = 41{,}65\ \mu\text{m}$, no domínio $D_{\text{mín}} = 0{,}1\ \mu\text{m}$ a $D_{\text{máx}} = 297\ \mu\text{m}$ (peneira 50#) (Tabela 5.9). Verificação: $F(38\ \mu\text{m}) = 0{,}598$, próximo das frações mássicas medidas ($53{,}6\%$ passante na peneira de $38\ \mu\text{m}$; $55{,}9\%$ abaixo de $36{,}5\ \mu\text{m}$ na difração a laser).

---

## 6. Conversão Fracionária da Zincita

Como $q_3$ já é ponderada pela massa, uma partícula de diâmetro inicial $D$ retém a fração $(1 - \Delta D/D)^3$ da sua massa:

$$X(t) = 1 - \frac{\int_{\max(\Delta D,\,D_{\text{mín}})}^{D_{\text{máx}}} \left( 1 - \frac{\Delta D}{D} \right)^3 q_3(D) \, dD}{\int_{D_{\text{mín}}}^{D_{\text{máx}}} q_3(D) \, dD}$$

> A tese usa a Eq. 3.49, $X = 1 - \int (D - \Delta D)^3 f\,dD \,/ \int D^3 f\,dD$, que só vale se $f$ for densidade **em número**. Aplicada à RRB mássica, pondera a massa duas vezes e torna a população efetiva muito mais grossa: com $\Delta D = 5\ \mu\text{m}$, a conversão sai $0{,}12$ em vez de $0{,}51$.

### Sistema Dinâmico Acoplado
$$\begin{cases}
\displaystyle \frac{d(\Delta D)}{dt} = \frac{2 \, k_s(T)}{\rho_s} \, [\text{H}^+](t) \\[1em]
\displaystyle X(t) = 1 - \frac{\int_{\max(\Delta D,\,D_{\text{mín}})}^{D_{\text{máx}}} (1 - \Delta D/D)^3 \, q_3(D) \, dD}{\int_{D_{\text{mín}}}^{D_{\text{máx}}} q_3(D) \, dD} \\[1em]
\displaystyle C_{Af}(t) = C_{A0} \left( 1 - \frac{X(t)}{\eta} \right) \\[1em]
\displaystyle [\text{H}^+](t) = f_{\text{especiação}}\big(C_{A0},\ C_{Af}(t);\ K_2\big)
\end{cases}$$

com $\Delta D(0) = 0$, $X(0) = 0$ e $C_{Af}(0) = C_{A0}$.

---

## 7. Tabela de Parâmetros

| Símbolo | Descrição | Valor | Unidade | Fonte |
| :--- | :--- | :---: | :---: | :--- |
| $k_{s,30^\circ\text{C}}$ | Constante de velocidade (base $[\text{H}^+]$, $30^\circ\text{C}$) | $1{,}8 \times 10^4$ | $\mu\text{m}\cdot\text{min}^{-1}$ | Balarini (2009); Balarini et al. (2025) |
| $E_a$ | Energia de ativação | $13{,}45$ | $\text{kJ}\cdot\text{mol}^{-1}$ | Balarini (2009) |
| $K_2$ | Constante da 2ª dissociação do $\text{H}_2\text{SO}_4$ | $0{,}0099$ | $\text{mol}\cdot\text{L}^{-1}$ | Implícita em Balarini et al. (2025) |
| $\alpha$ | Parâmetro da Eq. 4.1 (tese: $5{,}5 \times 10^3$) | $0$ | $\mu\text{m}\cdot\text{min}^{-1}$ | Este trabalho |
| $\rho_s$ | Densidade molar da zincita | $69{,}2$ | $\text{mol}\cdot\text{L}^{-1}$ | Perry e Green (1984) |
| $T_{\text{ZnO}}$ | Teor de zincita no concentrado ustulado | $76{,}1$ | $\% \text{ m/m}$ | Bortot Coelho (2017) |
| $M_{\text{ZnO}}$ | Massa molar da zincita | $81{,}38$ | $\text{g}\cdot\text{mol}^{-1}$ | Bortot Coelho (2017) |
| $m$ | Índice de uniformidade RRB (base mássica) | $1{,}022$ | adim. | Bortot Coelho (2017) |
| $D_{63{,}2}$ | Diâmetro característico RRB (base mássica) | $41{,}65$ | $\mu\text{m}$ | Bortot Coelho (2017) |
| $D_{\text{mín}}$ | Limite inferior da distribuição | $0{,}1$ | $\mu\text{m}$ | Bortot Coelho (2017) |
| $D_{\text{máx}}$ | Limite superior da distribuição | $297$ | $\mu\text{m}$ | Bortot Coelho (2017) |
| $V$ | Volume reacional de bancada | $0{,}40$ | $\text{L}$ | Bortot Coelho (2017) |
| $T$ | Temperatura de operação | $40{,}0$ | $^\circ\text{C}$ | Bortot Coelho (2017) |

---

## 8. Pontos da Formulação da Tese e Correções Adotadas

| # | Ponto da tese | Tratamento neste trabalho |
| :-: | :--- | :--- |
| 1 | **Base da RRB.** A RRB foi ajustada a frações mássicas (Tabela 5.4), mas é usada como distribuição em número (Eq. 3.45) e ponderada por $D^3$ (Eq. 3.49), contando a massa duas vezes. | **Corrigido:** integração em base mássica (Seção 6). |
| 2 | **Temperatura de $k_s$.** Estimado a $30^\circ\text{C}$ (Balarini et al., 2025), aplicado a $40^\circ\text{C}$ sem correção. | **Corrigido:** Arrhenius a partir de $30^\circ\text{C}$ (Seção 3.1). |
| 3 | **Base de $k_s$.** Estimado com $C_A = [\text{H}^+]$ (item 3.3.2 da tese), aplicado sobre $[\text{H}_2\text{SO}_4]$. | **Corrigido:** força motriz $[\text{H}^+]$ por especiação (Seção 3.2). Um fator fixo de conversão ($1{,}126$) não basta, porque o tamponamento pelo sulfato varia ao longo da reação. |
| 4 | **Parâmetro $\alpha$.** Pela Eq. 5.3, com $\eta = 0{,}535$ o modelo para em $X = 0{,}41$ e deixa 23% do ácido sem reagir, enquanto os ensaios 1, 3, 4 e 5 medem $X = 0{,}50$ e ácido residual $\le 5 \times 10^{-6}\ \text{mol/L}$. O texto e a Tabela 5.9 informam $5{,}5 \times 10^3$, o Apêndice A2.1 usa $5\,396$ e o mínimo da Eq. 5.4 está em $\approx 3{,}7 \times 10^3$. | **Removido** ($\alpha = 0$); o desvio é aprendido pelo ML. |
| 5 | **Consumo de ácido por outros minerais.** O ácido consumido é de $1{,}07$ a $1{,}13$ mol por mol de ZnO dissolvido. | **Mantido fora do modelo** (hipótese 6); fica no resíduo. |
| 6 | **Detalhes menores.** $\rho_B$ com unidade $\text{mol}\cdot\text{m}^{-3}$ na Tabela 5.9 (o correto é $\text{mol}\cdot\text{L}^{-1}$); $D_{\text{mín}} = 0{,}01\ \mu\text{m}$ no Apêndice A2.1 contra $0{,}1\ \mu\text{m}$ na Tabela 5.9; $D_{63{,}2} = 42{,}09\ \mu\text{m}$ no ajuste linear da Tabela 5.6 contra $41{,}65\ \mu\text{m}$ no modelo. | Adotados os valores da Tabela 5.9. |

Efeito de cada correção sobre o $R^2$ excluindo $t = 0$ (128 pontos, $\alpha = 0$):

| Formulação | $R^2$ ($t > 0$) | $RMSE$ |
| :--- | :---: | :---: |
| Tese ($\alpha = 0$) | $0{,}7475$ | $0{,}098$ |
| + RRB em base mássica | $0{,}8348$ | $0{,}079$ |
| + $k_s$ referido a $30^\circ\text{C}$ | $0{,}8269$ | $0{,}081$ |
| + força motriz $[\text{H}^+]$ por especiação (**modelo adotado**) | $\mathbf{0{,}9486}$ | $\mathbf{0{,}044}$ |
| Tese com $\alpha = 5500$ (referência) | $0{,}7697$ | $0{,}094$ |

---

## 9. Justificativa para a Modelagem Híbrida Grey-Box

O modelo *White-Box* corrigido apresenta $R^2 = 0{,}9837$ ($RMSE = 0{,}0414$, $MAE = 0{,}0309$; $R^2 = 0{,}9486$ excluindo os pontos triviais em $t = 0$). O resíduo $\Delta X = X_{\text{exp}} - X_{\text{PBM}}$ remanescente tem três assinaturas:
1. **Patamar incompleto com $\eta = 1$:** o modelo chega a $0{,}96\text{–}0{,}97$ e o experimento para em $0{,}85\text{–}0{,}87$ (resíduo de $-0{,}11$ em 15 min), com ácido livre residual medido de $0{,}056 / 0{,}12 / 0{,}23\ \text{mol/L}$ (previsto: $0{,}048 / 0{,}095 / 0{,}142$). A tese (item 5.2.2) explica esse platô pela queda da força motriz: cerca de $90\%$ do ácido é consumido no início da reação, quando a taxa é alta, e a conversão não é total mesmo com ácido e zincita em proporção estequiométrica. A especiação reproduz esse mecanismo e leva o modelo de $1{,}00$ a $0{,}96	ext{–}0{,}97$, mas não até $0{,}85	ext{–}0{,}87$. A parte restante pode vir de efeitos não incluídos (coeficientes de atividade, concentração de ZnSO₄ no licor, que a própria tese cita em 4.3) e dos demais minerais (ferrita, willemita, esfalerita), a que a tese atribui os $98\%$ com $\eta = 1{,}5$; isso não foi testado.
2. **Transientes iniciais ($t = 0{,}5$ min):** resíduo médio de $+0{,}09$ com $\eta = 1$ e de $-0{,}05$ com $\eta = 3{,}1$ (dentro de $\pm 0{,}10$ em todos os ensaios).
3. **Patamar com $\eta = 0{,}5$:** o modelo para em $X = \eta = 0{,}535$ e o experimento em $0{,}50$ ($-0{,}04$), pelo consumo de ácido por outros óxidos (Seção 8, item 5). Com $\eta \ge 1{,}5$, a conversão final é reproduzida ($-0{,}03$ a $0{,}00$).

A modelagem **Grey-Box** combina a estrutura estequiométrica do PBM com a flexibilidade do Machine Learning, que prevê esse resíduo. Em validação cruzada agrupada (GroupKFold, 8 dobras), os melhores regressores são o Gaussian Process ($R^2_{\text{CV}} = 0{,}9966$; $0{,}9893$ excluindo $t = 0$) e o Extra Trees ($0{,}9965$; $0{,}9891$), empatados também em 10 partições aleatórias dos ensaios ($0{,}9883 \pm 0{,}0013$ e $0{,}9888 \pm 0{,}0013$ com $t > 0$).

---

## Referências Bibliográficas
1. BALARINI, J. C. *Estudo Cinético do Ustulado de Zinco da Votorantim Metais: Mecanismo Não Catalítico Associado à Caracterização e à Determinação de Tempos de Residência*. Universidade Federal de Minas Gerais, 2009.
2. BALARINI, J. C.; ARAÚJO, E. M. R.; COELHO, F. E. B.; POLLI, L. O.; KONZEN, C.; MIRANDA, T. L. S.; SALUM, A. Leaching kinetics of roasted zinc concentrate in sulphuric acid solutions. *Revista Observatorio de la Economía Latinoamericana*, v. 23, n. 11, 2025. DOI: 10.55905/oelv23n11-188.
3. BORTOT COELHO, F. E. *Desenvolvimento e Validação de um Modelo de Balanço Populacional para a Lixiviação de um Concentrado Ustulado de Zinco em Bancada e em Escala Piloto*. 2017. 248 f. Dissertação (Mestrado em Engenharia Química) -- Universidade Federal de Minas Gerais, Belo Horizonte, 2017.
4. CRUNDWELL, F. K.; BRYSON, A. W. The calculation of the rate of dissolution of minerals in batch and continuous reactors. *Hydrometallurgy*, v. 29, n. 1-3, p. 275-295, 1992.
5. HERBST, J. A. Rate Processes in Multiparticle Metallurgical Systems. In: SOHN, H. Y.; WADSWORTH, M. E. (Ed.). *Rate Processes of Extractive Metallurgy*. New York: Plenum Press, 1979.
6. LEBLANC, S. E.; FOGLER, H. S. Population balance modeling of the dissolution of polydisperse solids: rate limiting regimes. *AIChE Journal*, v. 33, n. 1, p. 54-63, 1987.
