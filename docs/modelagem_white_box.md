# Documentação Teórica e Equacionamento do Modelo Mecanicista (White-Box) para a Lixiviação de Concentrado Ustulado de Zinco

**Trabalho de Conclusão de Curso (TCC) -- Engenharia Química -- UFMG**  
**Discentes:** Daniel Couto Vieira, Guilherme Moura de Sousa Franco, Matheus Henrique Borba Póvoas, Rodrigo Amaral da Mata  
**Orientador:** Prof. Dr. Fabrício Eduardo Bortot Coelho  
*Departamento de Engenharia Química -- Universidade Federal de Minas Gerais (UFMG)*  
*Outubro de 2026*

---

## Resumo
Este documento apresenta as hipóteses e o equacionamento do Modelo Mecanicista (*White-Box*) para a lixiviação ácida atmosférica de concentrado ustulado de zinco (calcina) em reator batelada. O ponto de partida é o modelo de **Bortot Coelho (2017)** (item 4.3, Tabela 5.9 e Apêndice A2.1 da dissertação): Balanço Populacional unidimensional resolvido pelo Método das Características, distribuição granulométrica de Rosin-Rammler-Bennett (RRB), cinética de reação química interfacial (Eq. 3.42) e consumo estequiométrico de ácido sulfúrico (Eq. 3.51). Sobre ele aplicam-se correções de consistência física: (i) a RRB, ajustada a frações mássicas, é integrada em base mássica; (ii) a constante $k_s$ é referida à temperatura em que foi estimada (30 °C) e corrigida por Arrhenius até 40 °C; (iii) a força motriz é a concentração de ácido livre $C_{Af}$ (Eq. 3.51), com consumo estequiométrico de ácido ($\nu = 1$). O parâmetro empírico $\alpha$ da Eq. 4.1 é fixado em zero: o desvio que ele corrigia passa a ser aprendido pelo modelo de Machine Learning do híbrido *Grey-Box*. Esse é o **White-Box principal** (`ConfigWB.laboratorio()`), formulado como se fôssemos ao laboratório, sem nenhum dado de conversão nem de pH. As variantes com especiação por pH (V0 a V2, com par iônico $\text{ZnSO}_4^0$) e com $\nu$ estimado pela composição (1,06) são **comparativas**: explicam o pH medido, mas não entram no principal, pela mesma razão de não ajustar nada aos dados. A formulação original da tese continua disponível no código (`ConfigWB.tese()`). Os resultados do híbrido e a validação estão em `validade_modelo_hibrido.md`.

---

## 1. Introdução e Hipóteses do Modelo

A lixiviação atmosférica de concentrado ustulado de zinco (calcina) em meio de ácido sulfúrico diluído constitui a etapa primordial da rota hidrometalúrgica de recuperação do zinco metálico. O concentrado ustulado é composto majoritariamente por zincita ($\text{ZnO}$, 76,1% m/m), mineral prontamente solúvel, e por ferrita de zinco ($\text{ZnFe}_2\text{O}_4$), willemita e esfalerita, fases de dissolução lenta em condições brandas.

Premissas do modelo (Bortot Coelho, 2017, item 4.3):
1. **Regime Hidrodinâmico:** agitação intensa ($1000\text{ rpm}$), com resistência difusional na camada limite desprezível.
2. **Mecanismo Controlador:** reação química superficial de primeira ordem (Balarini, 2009), com encolhimento contínuo das partículas. A força motriz do principal é $C_{Af}$; a variante comparativa usa $[\text{H}^+]$, a base em que Balarini estimou $k_s$.
3. **Morfologia das Partículas:** partículas esféricas.
4. **Fenômenos Populacionais:** sem quebra nem aglomeração ($B = D = 0$).
5. **Isotermicidade:** $T = 40^\circ\text{C}$.
6. **Seletividade Mineralógica:** a cinética considera apenas a dissolução da zincita (Eq. 4.2). O consumo de ácido pelos demais constituintes entra apenas pelo fator $\nu$ (Seção 2). A conversão experimental é normalizada pelo zinco da zincita, descontando o zinco da ferrita (Eq. 4.5).

---

## 2. Estequiometria e Balanço de Massa no Reator Batelada

$$\text{ZnO}_{(s)} + \text{H}_2\text{SO}_{4(aq)} \longrightarrow \text{ZnSO}_{4(aq)} + \text{H}_2\text{O}_{(l)} \qquad (4.2)$$

A razão molar é definida sobre a zincita (Eq. 4.3):

$$\eta = \frac{n_{A0}}{n_{B0}} = \frac{V \, C_{A0}}{m_{B0} \, T_{\text{ZnO}} / M_{\text{ZnO}}} \qquad (4.3)$$

com $T_{\text{ZnO}} = 0{,}761$ e $M_{\text{ZnO}} = 81{,}38\text{ g}\cdot\text{mol}^{-1}$. O modelo usa o $\eta$ calculado com as massas de cada ensaio (coluna `razao_molar_eta_real`): $0{,}535 / 1{,}069 / 1{,}58\text{–}1{,}61 / 3{,}19\text{–}3{,}29$ para os níveis nominais $0{,}5 / 1{,}0 / 1{,}5 / 3{,}1$ do planejamento. É com esses valores que se reproduz a soma quadrática de $\approx 0{,}05$ relatada no ajuste de $\alpha$ da tese (item 5.2.4).

A conversão experimental (Eq. 4.5) também está na base da zincita:

$$X_{\text{Zn}} = \frac{D \, V \, (C_{\text{Zn},f} - C_{\text{Fe},f}/2)}{m_{B0} \, T_{\text{ZnO}}} \qquad (4.5)$$

Balanço de ácido em sistema fechado (Eq. 3.51, Herbst, 1979), generalizado pelo fator $\nu$ (mol de $\text{H}_2\text{SO}_4$ consumidos por mol de ZnO dissolvido):

$$C_{Af}(t) = C_{A0} \left( 1 - \frac{\nu \, X(t)}{\eta} \right) \qquad (3.51)$$

**Fator $\nu$.** Nos ensaios com $\eta = 0{,}5$ o ácido é totalmente consumido (pH final 5 a 6) e a conversão final é $0{,}50$, não $\eta = 0{,}535$: logo $\nu = \eta/X_{final} = 1{,}07$ (1,09 no ensaio 5). Outros constituintes da calcina também consomem ácido (óxidos de Ca, Cu, Cd e Mn; ferrita, cujo zinco é descontado da conversão pela Eq. 4.5). O valor **a priori**, calculado só pela composição, é $\nu = 1{,}06 \pm 0{,}02$ (`kinetics.nu_a_priori`). Ele é **comparativo**: o White-Box principal usa $\nu = 1$ (Eq. 3.51 original), para que nenhum dado entre no fator de fechamento. O teto estequiométrico da conversão passa a ser $\min(1;\ \eta/\nu)$.

> **Atenção:** o "ácido livre final" da Tabela A1.1 **não é uma titulação**. A tese (itens 4.4.2 e 4.7.2 e Apêndice A4) o calcula a partir do pH medido, por especiação de ácido sulfúrico **puro** (balanço de cargas sem $\text{Zn}^{2+}$). Reconstruído do pH, ele reproduz o valor tabelado (razão mediana 0,97; o resto é arredondamento do pH). Os consumos de ácido de 1,07 a 1,13 mol/mol que vêm dele não são independentes do pH; o valor robusto é o dos ensaios com $\eta = 0{,}5$.

---

## 3. Cinética Química e Taxa de Variação do Diâmetro

Para controle por reação química, a taxa de variação do diâmetro (Eq. 3.42 da tese, LeBlanc e Fogler, 1987) é:

$$v \equiv -\frac{dD}{dt} = \frac{2 \, k_s}{\rho_s} \, C_A \qquad (3.42)$$

com $\rho_s = 69{,}2\ \text{mol}\cdot\text{L}^{-1}$.

### 3.1. Base e temperatura de $k_s$
Balarini (2009) estimou $k_s = 3{,}00 \times 10^{-2}\ \text{cm}\cdot\text{s}^{-1} = 1{,}8 \times 10^4\ \mu\text{m}\cdot\text{min}^{-1}$ usando, na Eq. 3.42, a **concentração hidrogeniônica** $C_A = [\text{H}^+]$ (item 3.3.2 da tese), a **30 °C** (Balarini et al., 2025, Fig. 8: $T = 30^\circ\text{C}$, $6\ \text{g/L}$ de $\text{H}_2\text{SO}_4$, $[\text{H}^+] = 0{,}0689\ \text{mol/L}$). A tese aplica esse valor a $40^\circ\text{C}$ e sobre a molaridade de $\text{H}_2\text{SO}_4$. Para referir $k_s$ à temperatura em que foi estimado:

$$k_s(T) = k_{s,30^\circ\text{C}} \exp \left[ -\frac{E_a}{R} \left( \frac{1}{T} - \frac{1}{303{,}15} \right) \right], \qquad E_a = 13{,}45\ \text{kJ}\cdot\text{mol}^{-1}$$

o que dá $k_s(40^\circ\text{C}) = 2{,}13 \times 10^4\ \mu\text{m}\cdot\text{min}^{-1}$.

**Força motriz do principal.** O principal aplica esse $k_s(40^\circ\text{C})$ sobre $C_A = C_{Af}$ (Eq. 3.51), como a tese, mas com a correção de Arrhenius. Isso mistura bases: $k_s$ foi estimado com $[\text{H}^+]$ ($0{,}0689$ mol/L para $0{,}0612$ mol/L de $\text{H}_2\text{SO}_4$, razão $1{,}126$). A inconsistência é **assumida**: nenhum dado de conversão ou de pH entra na construção do modelo. O que a troca de base custa na conversão aparece na escada da Seção 8: o RMSE do White-Box vai de $0{,}081$ com $C_{Af}$ para $0{,}044$ com $[\text{H}^+]$ (V0).

### 3.2. Concentração de $\text{H}^+$ por especiação (variantes comparativas)
O principal **não usa** esta seção. As variantes abaixo trocam a força motriz $C_{Af}$ por $[\text{H}^+]$ e existem para comparar o modelo com o pH medido e medir o efeito da especiação sobre a conversão.

O ácido sulfúrico tem a primeira dissociação completa e a segunda parcial, $\text{HSO}_4^- \rightleftharpoons \text{H}^+ + \text{SO}_4^{2-}$, com $K_2 = [\text{H}^+][\text{SO}_4^{2-}]/[\text{HSO}_4^-]$. Cada mol de ZnO dissolvido consome 2 mol de $\text{H}^+$, mas o sulfato permanece em solução. Com prótons totais $H_T = 2\,C_{Af}$ e sulfato total $S_T = C_{A0}$:

$$[\text{H}^+] + [\text{HSO}_4^-] = H_T, \qquad [\text{HSO}_4^-] + [\text{SO}_4^{2-}] + [\text{ZnSO}_4^0] = S_T$$

Três variantes foram avaliadas (`ConfigWB.especiacao`), todas em concentrações (sem coeficientes de atividade):

| Variante | $K_2$ | Par iônico $\text{ZnSO}_4^0$ | Fonte |
| :--- | :--- | :--- | :--- |
| **V0** (versão anterior) | $0{,}0099$ mol/L, constante | não | reproduz o $[\text{H}^+] = 0{,}0689$ mol/L de Balarini (6 g/L, 30 °C) |
| **V1** | van't Hoff: $0{,}00720$ mol/L a 40 °C | não | Apêndice A4 da tese: $\text{p}K_{a2} = 1{,}96$ a 25 °C, $\Delta H^\circ = -5{,}2$ kcal/mol |
| **V2** (melhor pH) | igual a V1 | sim: $K_z = [\text{ZnSO}_4^0]/([\text{Zn}^{2+}][\text{SO}_4^{2-}]) = 226$ L/mol a 40 °C | IUPAC (Powell et al., 2013, Parte 5): $\log K^\circ = 2{,}304$, $\Delta_r H^\circ = +6{,}0$ kJ/mol |

Sem par iônico, $[\text{HSO}_4^-] = \frac{(S_T + H_T + K_2) - \sqrt{(S_T + H_T + K_2)^2 - 4\,S_T H_T}}{2}$ e $[\text{H}^+] = H_T - [\text{HSO}_4^-]$. Com o par iônico (Zn dissolvido $= X\,C_{A0}/\eta$), resolve-se numericamente a equação monótona $K_2\,b = (H_T - b)\,s$ em $s = [\text{SO}_4^{2-}]$.

**Validação independente (pH).** O $K_2$ constante de V0 reproduz o $[\text{H}^+]$ de Balarini, mas **não** o pH final medido: com o estado final medido ($X_{final}$ e $\nu$), V0 prevê o pH 0,82 unidade alto demais (RMSE 0,99) e subestima o $[\text{H}^+]$ em até 56 vezes (ensaio 10: pH 2,35 previsto contra 0,6 medido), porque o sulfato do $\text{ZnSO}_4$ tampona o $\text{H}^+$ como $\text{HSO}_4^-$. O par iônico (V2) reduz o RMSE do pH final para **0,26** (viés +0,20) sem ajustar nenhum parâmetro; o pH inicial é reproduzido por todas as variantes (RMSE 0,052). Por isso V2 é a variante comparativa mais coerente com o pH; o principal não usa especiação (Seção 3).

**Ressalvas.** (i) $K_z$ é de força iônica zero e aplicada **sem correção de atividade**: V2 é uma cota superior do pareamento, e o viés residual é compatível com coeficientes de atividade não modelados (uma variante V3 é o passo seguinte). (ii) Os valores IUPAC foram obtidos de resumo secundário; conferir no relatório antes de citar. (iii) A tese lista $\text{p}K_{a2} = 1{,}96$ e também $K_{a2} = 12{,}383 \times 10^{-3}$ ($\text{p}K_{a2} = 1{,}91$); usa-se 1,96, que reproduz Balarini a 30 °C (−0,4%). Para o efeito de cada variante sobre a conversão, ver a Seção 8.

### 3.3. Parâmetro $\alpha$ da tese
Bortot Coelho (2017) adicionou um parâmetro ajustável à Eq. 3.42 para corrigir a conversão final (Eq. 4.1):

$$v = \frac{2}{\rho_s} \Big[ k_s \, C_{Af} - \alpha \, (C_{A0} - C_{Af}) \Big] \qquad (4.1)$$

com $\alpha = 5{,}5 \times 10^3\ \mu\text{m}\cdot\text{min}^{-1}$. **Neste trabalho, $\alpha = 0$**: o termo é um ajuste empírico das conversões finais e, quando o ácido é limitante, impede o seu esgotamento (Seção 8). No modelo híbrido, esse papel é assumido pelo regressor de Machine Learning. O modelo da tese com $\alpha = 5500$ é mantido como referência ("Tese 2017").

---

## 4. Balanço Populacional e Método das Características

Equação do balanço populacional em batelada (Eq. 3.38):

$$\frac{\partial \psi(D,t)}{\partial t} = -\frac{\partial}{\partial D} \big[ -v \, \psi(D,t) \big] \qquad (3.38)$$

Como $v$ não depende de $D$, a solução pelo Método das Características é um deslocamento da distribuição inicial:

$$\psi(D,t) = \psi\big(D + \Delta D(t),\ 0\big), \qquad \Delta D(t) = \int_0^t v(t') \, dt'$$

A tese resolve a EDP pelo Método Numérico das Linhas (`NDSolve` do Mathematica, Apêndice A2.1). A solução pelas características é a mesma solução exata, integrada aqui por RK45 (integrador de referência) e, para os experimentos, por LSODA sobre uma tabela $X(\Delta D)$ (diferença máxima < $10^{-4}$ em X).

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
\displaystyle \frac{d(\Delta D)}{dt} = \frac{2 \, k_s(T)}{\rho_s} \, C_{Af}(t) \\[1em]
\displaystyle X(t) = 1 - \frac{\int_{\max(\Delta D,\,D_{\text{mín}})}^{D_{\text{máx}}} (1 - \Delta D/D)^3 \, q_3(D) \, dD}{\int_{D_{\text{mín}}}^{D_{\text{máx}}} q_3(D) \, dD} \\[1em]
\displaystyle C_{Af}(t) = C_{A0} \left( 1 - \frac{\nu\,X(t)}{\eta} \right), \qquad \nu = 1
\end{cases}$$

com $\Delta D(0) = 0$, $X(0) = 0$ e $C_{Af}(0) = C_{A0}$. Nas variantes comparativas, a primeira equação passa a usar $[\text{H}^+](t) = f_{\text{especiação}}\big(C_{A0},\ C_{Af}(t),\ X(t)\,C_{A0}/\eta;\ K_2(T),\ K_z(T)\big)$ no lugar de $C_{Af}(t)$, e $\nu$ pode valer $1{,}06$.

---

## 7. Tabela de Parâmetros

| Símbolo | Descrição | Valor | Unidade | Fonte |
| :--- | :--- | :---: | :---: | :--- |
| $k_{s,30^\circ\text{C}}$ | Constante de velocidade (base $[\text{H}^+]$, $30^\circ\text{C}$) | $1{,}8 \times 10^4$ | $\mu\text{m}\cdot\text{min}^{-1}$ | Balarini (2009); Balarini et al. (2025) |
| $E_a$ | Energia de ativação | $13{,}45$ | $\text{kJ}\cdot\text{mol}^{-1}$ | Balarini (2009) |
| $K_2(T)$ | 2ª dissociação do $\text{H}_2\text{SO}_4$ (V1 e V2; comparativo) | $0{,}00720$ a $40^\circ\text{C}$ | $\text{mol}\cdot\text{L}^{-1}$ | Apêndice A4 da tese (van't Hoff) |
| $K_z(T)$ | Associação $\text{Zn}^{2+} + \text{SO}_4^{2-} \to \text{ZnSO}_4^0$ (V2; comparativo) | $226$ a $40^\circ\text{C}$ | $\text{L}\cdot\text{mol}^{-1}$ | IUPAC (Powell et al., 2013); *a conferir* |
| $\nu$ | Consumo de ácido por mol de ZnO | $1$ (principal); $1{,}06$ (comparativo) | adim. | estequiometria do ZnO (principal); composição da calcina (a priori, comparativo); $\eta/X_{final} = 1{,}07$ em $\eta = 0{,}5$ |
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

Nenhum parâmetro desta tabela é ajustado aos 128 pontos de conversão. O principal usa apenas $k_s$, $E_a$, $\rho_s$, a composição, a RRB e a geometria do ensaio; $k_s$, $E_a$ e a RRB vêm de ajustes anteriores (Balarini; granulometria medida), não dos pontos de conversão.

---

## 8. Pontos da Formulação da Tese e Correções Adotadas

| # | Ponto da tese | Tratamento neste trabalho |
| :-: | :--- | :--- |
| 1 | **Base da RRB.** A RRB foi ajustada a frações mássicas (Tabela 5.4), mas é usada como distribuição em número (Eq. 3.45) e ponderada por $D^3$ (Eq. 3.49), contando a massa duas vezes. | **Corrigido:** integração em base mássica (Seção 6). |
| 2 | **Temperatura de $k_s$.** Estimado a $30^\circ\text{C}$ (Balarini et al., 2025), aplicado a $40^\circ\text{C}$ sem correção. | **Corrigido:** Arrhenius a partir de $30^\circ\text{C}$ (Seção 3.1). |
| 3 | **Base de $k_s$.** Estimado com $C_A = [\text{H}^+]$ (item 3.3.2 da tese), aplicado sobre $[\text{H}_2\text{SO}_4]$. | **Mantido no principal, com variante comparativa:** o principal aplica $k_s$ sobre $C_{Af}$ (com Arrhenius, Seção 3.1). A força motriz $[\text{H}^+]$ por especiação (Seção 3.2) é a variante comparativa e ajusta melhor a conversão. Um fator fixo de conversão ($1{,}126$) não basta, porque o tamponamento pelo sulfato varia ao longo da reação. |
| 4 | **Parâmetro $\alpha$.** Pela Eq. 5.3, com $\eta = 0{,}535$ o modelo para em $X = 0{,}41$ e deixa 23% do ácido sem reagir, enquanto os ensaios 1, 3, 4 e 5 medem $X = 0{,}50$ e ácido residual $\le 5 \times 10^{-6}\ \text{mol/L}$. O texto e a Tabela 5.9 informam $5{,}5 \times 10^3$, o Apêndice A2.1 usa $5\,396$ e o mínimo da Eq. 5.4 está em $\approx 3{,}7 \times 10^3$. | **Removido** ($\alpha = 0$); o desvio é aprendido pelo ML. |
| 5 | **Consumo de ácido por outros constituintes.** O ácido consumido é de $1{,}07$ a $1{,}13$ mol por mol de ZnO dissolvido. | **Comparativo:** $\nu = 1$ no principal; $\nu = 1{,}06$ (a priori, pela composição) nas variantes comparativas (Seção 2). |
| 6 | **Detalhes menores.** $\rho_B$ com unidade $\text{mol}\cdot\text{m}^{-3}$ na Tabela 5.9 (o correto é $\text{mol}\cdot\text{L}^{-1}$); $D_{\text{mín}} = 0{,}01\ \mu\text{m}$ no Apêndice A2.1 contra $0{,}1\ \mu\text{m}$ na Tabela 5.9; $D_{63{,}2} = 42{,}09\ \mu\text{m}$ no ajuste linear da Tabela 5.6 contra $41{,}65\ \mu\text{m}$ no modelo. | Adotados os valores da Tabela 5.9. |
| 7 | **Especiação do H₂SO₄** (Apêndice A4) com constantes de 25 °C e sem $\text{Zn}^{2+}$ no balanço de cargas. | **Comparativo:** $K_2(T)$ por van't Hoff e par iônico $\text{ZnSO}_4^0$ (V2), validados pelo pH final medido, mas fora do principal (Seção 3.2). |
| 8 | **"Ácido livre final" da Tabela A1.1** apresentado como medida. | É **derivado do pH** por especiação ideal de ácido puro; o pH é a única medida independente (Seção 2). |

Efeito de cada passo sobre a conversão (métricas excluindo $t = 0$; 112 pontos) e sobre o pH (16 ensaios para pH₀, 12 com $\eta \ge 1$ para o pH final, estado final medido):

| Formulação | $R^2$ ($t > 0$) | RMSE | RMSE pH final |
| :--- | :---: | :---: | :---: |
| Tese ($\alpha = 0$) | $0{,}7475$ | $0{,}098$ | -- |
| Tese com $\alpha = 5500$ (referência) | $0{,}7697$ | $0{,}094$ | -- |
| + RRB em base mássica | $0{,}8347$ | $0{,}079$ | -- |
| + $k_s(T)$ de Arrhenius: **White-Box principal ($C_{Af}$, $\nu = 1$)** | $\mathbf{0{,}8269}$ | $\mathbf{0{,}081}$ | -- |
| + força motriz $[\text{H}^+]$, V0 ($K_2 = 0{,}0099$; versão anterior) | $0{,}9486$ | $0{,}044$ | $0{,}86$ |
| V1: $K_2(T)$ | $0{,}9524$ | $0{,}043$ | $0{,}96$ |
| V2: V1 + par iônico $\text{ZnSO}_4^0$ | $0{,}8525$ | $0{,}075$ | $\mathbf{0{,}13}$ |
| V0 + $\nu$ | $0{,}9593$ | $0{,}039$ | $0{,}99$ |
| V1 + $\nu$ | $0{,}9544$ | $0{,}042$ | $1{,}08$ |
| V2 + $\nu$ (comparativo; antigo principal) | $0{,}8992$ | $0{,}062$ | $0{,}26$ |

O White-Box principal é o degrau em negrito e é o que **menos** usa: nenhuma constante além de $k_s$, $E_a$ e da granulometria, sem $K_2$, par iônico nem $\nu$. Em troca, ajusta a conversão **pior** que todas as variantes com $[\text{H}^+]$ ($0{,}081$ contra $0{,}039$ a $0{,}075$). O $K_2$ constante de V0 reduz o $[\text{H}^+]$ ao longo da reação e, sem ter sido ajustado a esses dados, absorve parte do patamar em $\eta \approx 1$ que a física não explica. V2 não esconde esse patamar e é a única variante coerente com o pH, mas a escolha do principal não usa o pH. O efeito sobre o desempenho do híbrido (a extrapolação em $\eta$ piora com o principal) está em `validade_modelo_hibrido.md` (Seções 6.0 e 6.7).

---

## 9. Justificativa para a Modelagem Híbrida Grey-Box

O modelo *White-Box* principal ($C_{Af}$, $\nu = 1$) apresenta $R^2 = 0{,}9450$ nos 128 pontos e $0{,}8269$ excluindo $t = 0$ (RMSE $0{,}0813$ contra o piso de ruído experimental de $0{,}0175$). O resíduo $\Delta X = X_{\text{exp}} - X_{\text{PBM}}$ remanescente (Tabelas T14 e T16 de `resultados/tabelas.md`) tem as assinaturas:
1. **Patamar incompleto com $\eta = 1$:** aos 15 min o modelo chega a $0{,}997$--$1{,}000$ e o experimento para em $0{,}85$--$0{,}87$ (resíduo de $-0{,}13$ a $-0{,}15$, maior em módulo com $C_{A0} \ge 1$ mol/L). Nos ensaios 8 a 10 o pH final medido ($0{,}6$ a $1{,}2$) indica $[\text{H}^+]$ de $0{,}06$ a $0{,}25$ mol/L, com o qual o PBM dissolve tudo: o patamar **não** parece limitado pela acidez. A tese (item 5.2.2) o atribui à queda da força motriz, o que nem $C_{Af}$ nem a especiação com par iônico (V2) reproduzem. Hipóteses, **nenhuma testada**: zinco ocluído ou refratário (ferrita, silicatos), passivação por $\text{CaSO}_4$/$\text{PbSO}_4$ (1,4% de Ca e 2,3% de Pb na calcina), efeito de íon comum e atividade iônica, e erro da base de cálculo da Eq. 4.5.
2. **Transientes iniciais ($t = 0{,}5$ min):** o modelo é mais rápido que o experimento em 15 dos 16 ensaios, com resíduos de $-0{,}145$ ($\eta = 1{,}07$, $C_{A0} = 1{,}5$ mol/L, ensaio 10) a $+0{,}029$ ($\eta = 1{,}6$, $C_{A0} = 0{,}1$ mol/L, ensaio 11), bem acima do ruído experimental ($\sigma = 0{,}0175$); a causa não foi investigada.
3. **Patamar com $\eta = 0{,}5$:** com $\nu = 1$ o modelo para em $\eta = 0{,}535$ e o experimento em $0{,}49$--$0{,}50$, resíduo de $-0{,}035$ a $-0{,}045$. É o consumo extra de ácido ($\nu \approx 1{,}07$, Seção 2) que o principal deliberadamente não inclui; com $\nu = 1{,}06$ o resíduo era de $-0{,}004$ a $-0{,}014$. Com $\eta \ge 1{,}5$ o resíduo vai de $-0{,}03$ a $+0{,}01$.

A modelagem **Grey-Box** combina a estrutura estequiométrica do PBM com a flexibilidade do Machine Learning, que prevê esse resíduo. A validação por ensaio inteiro, contra o White-Box, uma tabela sem ML, o Ridge e uma caixa-preta (só ML), em seis protocolos, está em `validade_modelo_hibrido.md`. Em resumo: em interpolação o híbrido chega a $0{,}0215$--$0{,}0228$ de RMSE (próximo ao ruído), mas uma caixa-preta já se aproxima (a diferença é significativa só no Extra Trees); na **extrapolação em $\eta$** os híbridos erram $0{,}064$--$0{,}070$ (14 a 21% menos que o White-Box, $0{,}0813$), porém empatam com a tabela sem ML ($0{,}0637$). Dentre as caixas-pretas, Extra Trees ($0{,}0764$) e GBDT ($0{,}0794$) ficam um pouco melhores que o White-Box, e GP ($0{,}0978$) e SVR ($0{,}49$) piores; os híbridos superam todas. O domínio de validade são os quatro níveis de $\eta$ ensaiados.

---

## Referências Bibliográficas
1. BALARINI, J. C. *Estudo Cinético do Ustulado de Zinco da Votorantim Metais: Mecanismo Não Catalítico Associado à Caracterização e à Determinação de Tempos de Residência*. Universidade Federal de Minas Gerais, 2009.
2. BALARINI, J. C.; ARAÚJO, E. M. R.; COELHO, F. E. B.; POLLI, L. O.; KONZEN, C.; MIRANDA, T. L. S.; SALUM, A. Leaching kinetics of roasted zinc concentrate in sulphuric acid solutions. *Revista Observatorio de la Economía Latinoamericana*, v. 23, n. 11, 2025. DOI: 10.55905/oelv23n11-188.
3. BORTOT COELHO, F. E. *Desenvolvimento e Validação de um Modelo de Balanço Populacional para a Lixiviação de um Concentrado Ustulado de Zinco em Bancada e em Escala Piloto*. 2017. 248 f. Dissertação (Mestrado em Engenharia Química) -- Universidade Federal de Minas Gerais, Belo Horizonte, 2017.
4. CRUNDWELL, F. K.; BRYSON, A. W. The calculation of the rate of dissolution of minerals in batch and continuous reactors. *Hydrometallurgy*, v. 29, n. 1-3, p. 275-295, 1992.
5. HERBST, J. A. Rate Processes in Multiparticle Metallurgical Systems. In: SOHN, H. Y.; WADSWORTH, M. E. (Ed.). *Rate Processes of Extractive Metallurgy*. New York: Plenum Press, 1979.
6. LEBLANC, S. E.; FOGLER, H. S. Population balance modeling of the dissolution of polydisperse solids: rate limiting regimes. *AIChE Journal*, v. 33, n. 1, p. 54-63, 1987.
7. POWELL, K. J. et al. Chemical speciation of environmentally significant metals with inorganic ligands. Part 5: The Zn²⁺ + OH⁻, Cl⁻, CO₃²⁻, SO₄²⁻, and PO₄³⁻ systems (IUPAC Technical Report). *Pure and Applied Chemistry*, v. 85, n. 12, p. 2249-2311, 2013. DOI: 10.1351/PAC-REP-13-06-03. *(valores conferidos apenas em resumo secundário)*
8. SANSANA, J. et al. (2024), conforme a revisão bibliográfica do grupo (`data/ENTREGA PARCIAL II - Referência Bibliográfica.docx.md`), sobre arquiteturas de modelagem híbrida.
