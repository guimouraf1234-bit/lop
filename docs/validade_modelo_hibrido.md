# Validade do modelo híbrido da lixiviação em batelada

**Escopo:** só a batelada (16 ensaios, 128 pontos). A planta piloto e o híbrido serial ficam fora (seção 9).
**Atualização (2026-10-11):** o White-Box principal passou a ser o de laboratório (CAf, ν = 1); ver seções 1 e 6.0.
**Fonte de todos os números:** `resultados/*.csv`, geradas por `scripts/experimentos_batelada.py`; as tabelas `T##` estão em
`resultados/tabelas.md` e as figuras `F##` em `resultados/figuras/` (regeneráveis; ver `scripts/docs/README_experimentos_batelada.md`).

---

## 1. Resumo

> **Atualização de 2026-10-11.** O White-Box principal é o **de laboratório** (`ConfigWB.laboratorio()`): força motriz $C_{Af}$ (Eq. 3.51), ν = 1, α = 0, $k_s(T)$ de Balarini com correção de Arrhenius e RRB mássica medida, sem dado de conversão nem de pH na construção. A versão anterior (V2 + ν = 1,06, escolhida pelo pH) passa a ser comparativa. Os números abaixo são do principal atual; a comparação com a versão anterior está na seção 6.0. As seções 4.2, 4.4 e 6.7 analisam a versão anterior e estão marcadas.

1. **Interpolação: os híbridos não mudam.** RMSE de **0,0215 a 0,0228** (GKF8 e LOGO), dentro do mesmo intervalo da versão anterior (0,0208 a 0,0230). O White-Box sozinho sobe de 0,0597 para **0,0813**. Contra a caixa-preta, o ganho é pequeno: o Extra Trees híbrido é significativamente melhor (Δ = −0,0043, IC 95% [−0,0082; −0,0004], GKF8 × 20), mas GP e GBDT empatam com a caixa-preta correspondente.
2. **Extrapolação em η: piora, e o híbrido deixa de superar a tabela.** Com a deixa-um-η-fora, os híbridos vão a **0,064 a 0,070** (Extra Trees 0,0644; GBDT 0,0698; GP 0,0672), contra 0,043 a 0,054 da versão anterior. A tabela (η, t), que é um baseline só de dados, sai de 0,0527 para 0,0637 porque o acoplador passou a usar o teto η/ν com ν = 1. Assim, o Extra Trees híbrido empata com a tabela, enquanto antes ficava 15% abaixo dela. Contra a caixa-preta, o híbrido continua melhor (Extra Trees: Δ = −0,0120, IC 95% [−0,0206; −0,0053]).
3. **Poucos dados (4 ensaios de treino):** híbridos 0,050 (GP) e 0,052 (Extra Trees), contra 0,043 a 0,047 antes e 0,065 a 0,087 da caixa-preta. Com 8 ensaios, o GP híbrido (0,032) fica atrás da caixa-preta GP (0,028).
4. **O que explica a perda.** A força motriz [H⁺] é o que faz o White-Box ajustar a conversão: com [H⁺] e K₂ constante (V0), o RMSE cai de 0,0813 para 0,0443. CAf não tem esse efeito. A variante V1 (K₂ por van't Hoff com as constantes do Apêndice A4 da tese, também sem ajuste aos dados e com ν = 1) dá RMSE de 0,0426 na conversão, mas erro de pH combinado de 0,63, contra 0,17 da V2 + ν.
5. **Intervalos:** cobertura de 90 a 96% em LOGO (nominal 95%) e de **36 a 60%** quando se extrapola em η (antes, 49 a 71%).
6. **Patamar em η ≈ 1:** o principal dissolve tudo em 15 min (X = 1,000), e o experimento para em 0,85 a 0,87. O resíduo fica em −0,13 a −0,15 em η = 1,07. É o mesmo problema da versão anterior, com causa ainda desconhecida (seção 7).
7. **Domínio:** a discordância entre modelos é, em média, 0,016 nos níveis de η ensaiados (η real ≈ 0,54; 1,07; 1,6; 3,2) e 0,030 fora deles. Em η ≈ 0,85, a média é 0,043 (máximo 0,055); o máximo geral é 0,081.

---

## 2. O que mudou em relação à versão anterior

| Ponto | Antes | Agora | Por quê |
| :--- | :--- | :--- | :--- |
| White-Box principal | V2 + ν = 1,06, escolhido pelo pH medido | CAf livre, ν = 1, sem especiação e sem dado de conversão nem de pH (`ConfigWB.laboratorio()`) | decisão do usuário (2026-10-11): o principal não pode usar dados para ajustar parâmetros |
| Métricas | R² global (inclui 16 pontos triviais em t = 0), "violação 0%" depois do clip | RMSE/MAE/viés/R² **só em t > 0**, ganho sobre o White-Box e seu **teto pelo ruído**, violações **na saída bruta** | o R² global é dominado pela variância entre níveis de η, que o White-Box já explica |
| Validação | GroupKFold(8) com η nominal igual no treino e no teste | 6 protocolos: GKF8 repetido (20), LOGO, **deixa-um-η-fora**, deixa-um-$C_{A0}$-fora, **extrapolação temporal**, curva de aprendizado | o GroupKFold só mede interpolação em $C_{A0}$, que quase não importa |
| Baselines | nenhum | White-Box, tabela por (η, t), Ridge, **caixa-preta** | sem eles não se sabe o que a física acrescenta |
| Estatística | R² a 4 casas, ranking | **IC95% por bootstrap de ensaios**, diferenças pareadas, empates declarados, 20 partições | a ordem entre os melhores não é significativa |
| Hiperparâmetros | escolhidos com o mesmo CV que se reporta | **CV aninhada** (SVR, MLP) | evita viés otimista |
| Acoplador | teto com **η nominal**, máximo acumulado | teto **η_real/ν**, X(0) = 0, **regressão isotônica** | o teto antigo era uma correção oculta (equivalia a ν = 1,07) |
| White-Box | K₂ = 0,0099 constante | K₂(T) + par iônico, ν; escolha por pH | o [H⁺] antigo ficava até 56 vezes abaixo do medido |
| Cache | CSV do resíduo lido por `preparar_dados_base` (ignorava o DataFrame recebido) | PBM **sempre recalculado** | risco de dados obsoletos |
| Código | 4 laços de gráfico, 8 pastas de figuras, 3 scripts quebrados | pipeline único, testes, `arquivo/` | os números divergiam entre figuras |

Itens do parecer anterior corrigidos nesta revisão: (a) o "ácido livre final" da Tabela A1.1 **não é uma titulação**, é calculado a partir do pH
(seção 4.3), logo não é medida independente; (b) o consumo de ácido "medido" ν = 1,07–1,13 só é robusto nos ensaios com η = 0,5; (c) o fator
"30 a 110×" de [H⁺] era ponta a ponta; no estado final medido o pior caso é 56×.

---

## 3. Arquitetura e proveniência dos parâmetros

**Paralela cooperativa** (Sansana et al., 2024): $X = \mathcal{R}[X_{PBM} + \widehat{\Delta X}_{ML}(t,\ \eta_{real},\ C_{A0},\ X_{PBM})]$.
O critério de Sansana et al. para essa arquitetura (modelo de primeiros princípios com limitação estrutural mas bom comportamento geral, em
**domínio conhecido e controlado**) é exatamente o que os resultados mostram: ela funciona dentro do planejamento e perde força fora dele.

| Parâmetro do White-Box | Valor | Proveniência | Ajustado aos 128 pontos? |
| :--- | :--- | :--- | :--- |
| RRB (m, D₆₃,₂, limites) | 1,022; 41,65 µm; 0,1 a 297 µm | Tabela 5.9, ajustada a frações **mássicas**; integrada em base mássica | não |
| $k_s$, $E_a$ | 1,8·10⁴ µm/min a 30 °C; 13,45 kJ/mol | Balarini (2009); Arrhenius de 30 para 40 °C | não |
| $K_2(T)$ | 0,00720 mol/L a 40 °C | Apêndice A4 da tese: pKa₂ = 1,96 (25 °C), ΔH° = −5,2 kcal/mol | não |
| $K_z(T)$ (ZnSO₄⁰) | 226 L/mol a 40 °C | IUPAC (Powell et al., 2013): log K° = 2,304, ΔrH° = +6,0 kJ/mol (**valores de resumo secundário: conferir**) | não |
| ν | 1,06 | composição da calcina (óxidos de Ca, Cu, Cd, Mn e ferrita); η = 0,5 dá 1,07 | não (a priori) |
| $T_{ZnO}$, $M_{ZnO}$, $\rho_s$ | 0,761; 81,38 g/mol; 69,2 mol/L | tese | não |

Atributos do ML: conjunto **B** = (t, η real, $C_{A0}$, $X_{PBM}$), escolhido pela ablação (seção 6.5). **O White-Box e os atributos foram escolhidos por regras
declaradas e reproduzíveis** (`validacao.selecionar_whitebox`, `experimentos_batelada.escolher_entradas`); um teste falha se as constantes divergirem.

---

## 4. White-Box

> **Nota (2026-10-11):** esta seção analisa a versão anterior do White-Box (V2 + ν = 1,06, escolhida pelo pH). O principal atual não usa especiação nem ν; a escada da 4.1 continua válida como comparação, e o resultado do principal atual está na seção 6.0.

### 4.1. Escada de formulações (T02, F04)
A correção da RRB para base mássica e a especiação de H⁺ reproduzem os valores já documentados (R² em t > 0: tese 0,7475 → 0,8347 → 0,8269 → **0,9486**).
Os novos degraus são:

| Degrau | RMSE de X (t>0) | RMSE pH_f |
| :--- | ---: | ---: |
| V0 (K₂ = 0,0099, TCC) | 0,0443 | 0,865 |
| V1 (K₂(T)) | 0,0426 | 0,959 |
| V2 (V1 + par iônico) | 0,0750 | **0,128** |
| V0 + ν | **0,0394** | 0,987 |
| V1 + ν | 0,0417 | 1,085 |
| **V2 + ν (principal)** | 0,0620 | 0,256 |

### 4.2. Por que V2 (T03, T04, F05)

> Comparativo. Esta análise justificava V2 como principal pelo pH medido; esse critério foi abandonado em 2026-10-11 (seção 1).
Com K₂ constante, a especiação ideal prevê o pH final **0,8 a 0,9 unidade alcalino demais** (viés +0,82 para V0; até 1,75 unidade no ensaio 10): o sulfato do ZnSO₄
tampona o H⁺ como HSO₄⁻ e o modelo subestima o H⁺ livre em até **56 vezes**. O par iônico libera o H⁺ (sequestra sulfato como ZnSO₄⁰) e reduz o erro
de pH final de 0,99 para 0,26 unidade (viés +0,20), **sem ajustar nada**. O pH inicial é reproduzido em qualquer variante (RMSE 0,052).

**Ressalvas.** (i) A constante IUPAC é de força iônica zero e foi usada sem correção de atividade; V2 é uma **cota superior** do pareamento. O erro que sobra (V2
ainda prevê pH alto demais, viés +0,20) é compatível com coeficientes de atividade não modelados. (ii) A comparação supõe o pH medido em 40 °C (a tese diz que o termopar
estava acoplado ao medidor) e ideal (sem atividade). (iii) O pH ponta a ponta de V2 (erro 0,89) é pior que o por especiação (0,26): o PBM superestima a conversão em η ≈ 1 e consome ácido demais.

### 4.3. O "ácido livre final" não é medido
A tese (itens 4.4.2 e 4.7.2, Apêndice A4) calcula `C_acid_final` a partir do pH medido, por especiação de **ácido sulfúrico puro** (balanço de cargas sem Zn²⁺). Reconstruindo-o
do pH com as constantes do Apêndice A4 obtém-se 0,97 do valor tabelado (mediana; o resto é arredondamento do pH em 1 casa). Consequências: o pH é a **única** medida independente
do estado ácido, e os ν "medidos" de 1,07–1,13 vêm de uma especiação que ignora o Zn²⁺. O valor robusto vem dos ensaios com η = 0,5: o ácido se esgota (pH 5 a 6), X final = 0,50 e
η real = 0,535, logo ν = 1,07 (1,09 no ensaio 5). A estimativa a priori (1,06) concorda.

### 4.4. ν (T11)

> Comparativo. O principal atual usa ν = 1; esta seção mede o efeito de ν na versão anterior.
Estimar ν dentro de cada dobra (mediana dos ensaios de treino com η ≤ 1,5) ou usar o valor a priori dá resultados equivalentes (White-Box 0,0597 contra 0,0620; Extra Trees em extrapolação 0,0450 contra 0,0471).
**Sem ν (ν = 1)** o White-Box piora para 0,0750 e o Extra Trees em extrapolação em η para 0,0626: o ν importa.

### 4.5. O patamar em η ≈ 1 (T14, T16, F09)
Aos 15 min o resíduo do White-Box principal é ≈ 0 em η = 0,5 (−0,004 a −0,014) e η = 3,1, −0,02 a −0,03 em η = 1,5 e **−0,10 a −0,14 em η = 1,0**, crescendo com $C_{A0}$.
Nos ensaios 8 a 10 o pH final medido (0,6 a 1,2) indica [H⁺] de 0,06 a 0,25 mol/L, com o qual o PBM dissolve tudo; portanto o patamar **não** parece ser limitado pela acidez.

---

## 5. Protocolos de validação
Ver `scripts/docs/README_validacao.md`. Pontos de rigor:
* partições **por ensaio inteiro**, com o PBM recalculado e ν estimado só no treino em cada dobra; um **teste de vazamento** (trocar os dados medidos dos ensaios de teste não pode mudar a previsão) passa, e um mutante que treina com tudo é detectado;
* métricas só em t > 0; **σ = 0,0175** (64 graus de liberdade; 0,0156 por von Neumann) estimado nos patamares das curvas; o ganho máximo possível sobre o White-Box é 0,913;
* IC95% por **bootstrap de ensaios**; comparações pareadas; "empate" = IC da diferença contém 0.

---

## 6. Resultados

> **Nota (2026-10-11):** os números das subseções 6.1 a 6.7 são da versão anterior (V2 + ν = 1,06). O resultado do principal atual está na 6.0, com a comparação lado a lado.

### 6.0. Principal atual (CAf, ν = 1) e comparação com a versão anterior (V2 + ν)

RMSE em t > 0, no formato **anterior → atual**. O baseline de tabela (η, t) também muda, porque o acoplador usa o teto η/ν.

| Modelo | GKF8 × 20 | LOGO | Deixa-um-η-fora | Poucos dados (n = 4) |
| :--- | :---: | :---: | :---: | :---: |
| White-Box | 0,0594 → **0,0813** | 0,0597 → **0,0813** | 0,0595 → **0,0813** | 0,0571 → **0,0796** |
| Tabela (η, t) | 0,0306 → 0,0299 | 0,0303 → 0,0296 | 0,0527 → **0,0637** | — |
| Híbrido GP [B] | 0,0224 → 0,0223 | 0,0208 → 0,0216 | 0,0536 → **0,0672** | 0,0467 → **0,0501** |
| Híbrido Extra Trees [B] | 0,0230 → 0,0228 | 0,0224 → 0,0218 | 0,0450 → **0,0644** | 0,0431 → **0,0524** |
| Híbrido GBDT [B] | 0,0209 → 0,0220 | 0,0211 → 0,0215 | 0,0432 → **0,0698** | — |
| Caixa-preta GP | 0,0234 → 0,0239 | 0,0236 → 0,0242 | 0,0950 → 0,0978 | 0,0606 → 0,0650 |
| Caixa-preta Extra Trees | 0,0270 → 0,0271 | 0,0264 → 0,0265 | 0,0727 → 0,0764 | 0,0822 → 0,0873 |

Leitura: a interpolação não muda. A extrapolação em η piora, e o híbrido de árvore deixa de ficar abaixo da tabela. A perda vem da força motriz: a escada da seção 4.1 mostra que [H⁺] (V0, com K₂ constante) leva o White-Box de 0,0813 para 0,0443, e CAf não tem esse efeito. Os dados de origem são os mesmos; só a física do principal mudou.

### 6.1. Interpolação (GKF8 do TCC, GKF8 × 20, LOGO, $C_{A0}$ ausente; T05, F02, F12)
| Modelo | GKF8 (TCC) | GKF8 × 20 (± dp) | LOGO | $C_{A0}$ ausente |
| :--- | ---: | ---: | ---: | ---: |
| White-Box | 0,0597 | 0,0594 ± 0,0004 | 0,0597 | 0,0597 |
| Tabela (η, t) | 0,0303 | 0,0306 ± 0,0009 | 0,0303 | 0,0303 |
| Ridge | 0,0384 | 0,0388 ± 0,0013 | 0,0389 | 0,0355 |
| Híbrido MLP | **0,0185** | **0,0196 ± 0,0017** | **0,0173** | 0,0240 |
| Híbrido GBDT | 0,0211 | 0,0209 ± 0,0009 | 0,0211 | 0,0239 |
| Híbrido GP | 0,0227 | 0,0224 ± 0,0012 | 0,0208 | **0,0225** |
| Híbrido Extra Trees | 0,0226 | 0,0230 ± 0,0007 | 0,0224 | 0,0233 |
| Caixa-preta GP | 0,0236 | 0,0234 ± 0,0012 | 0,0236 | 0,0274 |
| Caixa-preta GBDT | 0,0240 | 0,0238 ± 0,0010 | 0,0236 | 0,0288 |

* Em GKF8 e LOGO todos os híbridos (exceto o SVR, 1,35σ) ficam entre 1,0 e 1,3σ; o melhor (MLP, LOGO 0,0173) está no próprio piso de ruído (ganho/teto = 100%).
* **Hierarquia estatística:** a MLP é significativamente melhor que a maioria dos modelos (por 0,002 a 0,009; T08), mas empata com GP, GBDT e Random Forest no GKF8 × 20; a diferença prática é pequena. Entre GBDT, RF, GP e ET há empates.
* **A física acrescenta pouco aqui:** híbrido contra a caixa-preta do mesmo modelo (GKF8 × 20), Δ RMSE = −0,0010 (GP), −0,0029 (GBDT), −0,0040 (ET); só o ET exclui o zero, por pouco (T07).
* A caixa-preta SVR (kernel isotrópico) é muito pior (0,0855 a 0,1006), por limitação do kernel com um único γ para t, η e $C_{A0}$; os comparadores relevantes são GP e GBDT.
* A **tabela por (η, t) sem ML** (0,0303) é melhor que o White-Box e o Ridge, mas perde para o melhor híbrido por 0,008 a 0,013 (significativo; T08).
* **Auditoria da saída bruta** (antes do acoplador): a caixa-preta GP viola o teto em 16% dos pontos e tem 24% de passos decrescentes; o híbrido GP, 9,6% e 3,7%. A estrutura física melhora a saída bruta, embora o acoplador a corrija.

### 6.2. Extrapolação em η (T05, T06, F07)
| Modelo | RMSE | ganho/teto |
| :--- | ---: | ---: |
| Híbrido Random Forest | **0,0417** | 56% |
| Híbrido GBDT | 0,0432 | 52% |
| Híbrido Extra Trees | 0,0450 | 47% |
| Híbrido MLP | 0,0467 | 42% |
| Ridge | 0,0507 | 30% |
| Tabela (η, t) | 0,0527 | 24% |
| Híbrido GP | 0,0536 | 21% |
| White-Box | 0,0595 | 0 |
| Caixa-preta ET / GBDT / GP | 0,0727 / 0,0758 / 0,0950 | negativo |

* O melhor híbrido supera o White-Box em Δ = −0,0178 (IC [−0,0252; −0,0080]), a tabela em −0,0110 e o Ridge em −0,0091: **todos significativos**. Contra a caixa-preta do mesmo modelo, a diferença é de −0,028 a −0,041 (GP, ET, GBDT), também significativa.
* Por nível ausente (T06), o melhor híbrido (RF) ganha do White-Box em η = 1,0 (0,065 contra 0,100), 1,5 (0,033 contra 0,050) e 3,1 (0,034 contra 0,035, marginal); **perde em η = 0,5** (0,023 contra 0,021). Em η = 1,5 a tabela e o Ridge (0,027) ganham dos híbridos.
* **GP e SVR são os piores híbridos aqui** (0,0536 e 0,0582): longe dos dados o GP reverte à média a priori do resíduo (zero) e devolve o White-Box, enquanto as árvores repetem o resíduo do nível vizinho. A caixa-preta GP falha por outro motivo: com só (t, η, $C_{A0}$) o ARD aprende escala de comprimento de 0,044 desvios-padrão em η (T15), isto é, trata cada nível de η isoladamente.
* O ceiling estequiométrico η_real/ν salva a caixa-preta em η = 0,5 (0,024) e em nenhum outro nível.

### 6.3. Extrapolação temporal e em $C_{A0}$ (T05)
* Treinando com t ≤ 5 min e prevendo 15 min: híbridos RF/GBDT/ET 0,011 a 0,012, **empatados com a caixa-preta ET** (0,0121); o GP híbrido 0,0177 contra 0,200 da caixa-preta GP. O teste é fácil (a curva está no patamar) e a tabela não extrapola (= White-Box).
* $C_{A0}$ ausente: híbridos 0,0225 a 0,0248, caixa-pretas 0,0274 a 0,0288; empatados em 2 dos 3 pares (GP e GBDT; só o Extra Trees exclui o zero, por pouco).

### 6.4. Curva de aprendizado (T09, F03)
| Ensaios de treino | 4 | 8 | 14 |
| :--- | ---: | ---: | ---: |
| Melhor híbrido | 0,0390 (GBDT) | 0,0280 (GBDT) | 0,0214 (GP) |
| Melhor caixa-preta | 0,0861 (GP) | 0,0311 (GBDT) | 0,0232 (GP) |
| Tabela (η, t) | 0,0466 | 0,0365 | 0,0301 |
| White-Box | 0,0581 | 0,0568 | 0,0591 |

A vantagem da física é de ~2× com 4 ensaios e cai para ~10% com 14: é um argumento forte para usar o híbrido **enquanto os dados forem escassos**.

### 6.5. Atributos (T10, F06, T15)
Informar a conversão do White-Box ($X_{PBM}$, conjunto B) reduz o erro frente ao conjunto A, sobretudo em extrapolação: GBDT 0,0324 → 0,0210 (GKF8) e 0,0634 → 0,0432 (η ausente); Extra Trees 0,0638 → 0,0450; SVR 0,0849 → 0,0582.
O conjunto C (mais $C_{Af}$ e [H⁺]) empata com B; **por regra, vale o mais simples (B)**. O GP com C degenera (ruído branco → 0 e escalas de 10⁵).

### 6.6. Intervalos (T13, F08)
Cobertura de 95% nominal: **89 a 91% em LOGO** (GP analítico 89%, conformes 90–91%), **49 a 71% com η ausente** (GP analítico 71%, conforme 49 a 66%). Os intervalos calibrados por interpolação não valem para extrapolar.

### 6.7. O mesmo ML sobre degraus do White-Box (T12, F11): o trade-off de D2

> Resolvido em 2026-10-11: o usuário escolheu o principal de laboratório (CAf, ν = 1). A análise abaixo vale para a escolha anterior, e a escada da F11 agora inclui o principal atual como degrau 4.
| RMSE (Extra Trees / GP) | White-Box do TCC (V0, ν = 1) | **V0 + ν** | **V2 + ν (principal)** |
| :--- | ---: | ---: | ---: |
| GKF8 × 10 | 0,0200 / 0,0235 | 0,0196 / 0,0237 | 0,0229 / 0,0226 |
| η ausente | 0,0341 / 0,0455 | **0,0290 / 0,0392** | 0,0450 / 0,0536 |
| n = 4 ensaios | 0,0351 / 0,0335 | **0,0296 / 0,0315** | 0,0431 / 0,0467 |
| n = 8 ensaios | 0,0225 / 0,0279 | **0,0207 / 0,0261** | 0,0271 / 0,0291 |

Mais física nem sempre melhora o híbrido: **V0 + ν é a melhor base preditiva** (10 a 36% menos erro em extrapolação em η e com poucos dados), mas é o menos coerente com o pH. O K₂ constante
funciona como um parâmetro efetivo que absorve parte do patamar em η ≈ 1 (sem ter sido ajustado a esses dados). **Em interpolação não há diferença relevante.** A escolha entre as duas bases é
do grupo e deve ser **declarada e justificada**; este trabalho segue a regra D2 (V2) e mostra a alternativa.

---

## 7. Limitações e hipóteses (não testadas)
* **Domínio:** 4 níveis de η, 4 de $C_{A0}$, t ≤ 15 min, 40 °C, 1 000 rpm. Fora disso o resultado é extrapolação (F10).
* **Ruído:** os dados são médias de triplicatas; sem as réplicas individuais não há σ por ponto (pedido em `docs/pedido_replicas_prof_fabricio.md`). O piso de ruído é uma estimativa (±10%): por isso o melhor híbrido pode "atingir" 100% do ganho.
* **Hiperparâmetros:** ET, RF e GBDT usam os valores do TCC (não aninhados); SVR e MLP, CV aninhada leve (grade pequena).
* **IUPAC:** constantes de resumo secundário; **atividade iônica** não modelada (um V3 com coeficientes de atividade é a continuação natural).
* **Patamar em η ≈ 1 — hipóteses:** zinco ocluído ou refratário (ferrita, silicatos), passivação por CaSO₄/PbSO₄ (1,4% de Ca e 2,3% de Pb na calcina), efeito de íon comum e atividade, erro na base de cálculo da Eq. 4.5. Nenhuma foi testada; o ML só o descreve empiricamente.
* **Teste cego real:** todos os "testes cegos" são fora-da-dobra sobre os mesmos 16 ensaios. Só ensaios novos (por exemplo η = 0,75 e 1,25) seriam um teste independente.

## 8. Recomendações
1. Apresentar o **GBDT** como híbrido de referência (nunca pior que o 3º em 6 protocolos) e o **GP** quando se precisar de incerteza; declarar a MLP como a melhor em interpolação.
2. Declarar o **domínio** e o fato de que a vantagem do híbrido sobre uma caixa-preta é **grande em extrapolação e com poucos dados, pequena em interpolação**.
3. Decidir a base: o principal atual é o de laboratório (CAf, ν = 1), decidido em 2026-10-11. Se a extrapolação em η for prioridade, V0 + ν (seção 6.7) continua sendo a alternativa comparativa.
4. Enviar o pedido de réplicas; planejar 2 a 3 ensaios em η intermediários.

## 9. Fora do escopo
Planta piloto (`continuous_reactor_model.py`: fator 0,05 arbitrário; volumes 18 L no JSON contra 6 L na dissertação) e o **híbrido serial** (ML estimando parâmetros dentro do PBM), que garantiria X(0) = 0,
monotonia e estequiometria por construção; a branch `feat/modelagem-alpha` (do grupo) já explora uma versão serial para α, e o resultado de extrapolação aqui o motiva.

## 10. Reprodução
```powershell
& ..\.venv\Scripts\python.exe -m pytest                           # 56 testes
& ..\.venv\Scripts\python.exe scripts/experimentos_batelada.py    # ~11 min
& ..\.venv\Scripts\python.exe scripts/figuras_batelada.py
& ..\.venv\Scripts\python.exe scripts/tabelas_batelada.py
```
