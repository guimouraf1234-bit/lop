# Tabelas dos resultados da batelada

Gerado por `scripts/tabelas_batelada.py` a partir de `resultados/*.csv`. Não edite à mão.

### T01. Piso de ruído experimental

| Grandeza | Valor |
| :--- | ---: |
| sigma (patamares) | 0,0175 |
| sigma (diferenças, von Neumann) | 0,0156 |
| graus de liberdade | 64 |
| sigma eta=0.5 | 0,0180 |
| sigma eta=1 | 0,0162 |
| sigma eta=1.5 | 0,0222 |
| sigma eta=3.1 | 0,0097 |

### T02. Escada do White-Box (X em 112 pontos; pH em 16 e 12 ensaios)

| Degrau | RMSE (t>0) | R² (t>0) | R² (todos) | RMSE pH₀ | RMSE pH_f (espec.) | RMSE pH_f (ponta a ponta) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| Tese (α = 0) | 0,0981 | 0,7475 | 0,9198 | – | – | – |
| Tese (α = 5500) | 0,0937 | 0,7697 | 0,9268 | – | – | – |
| + RRB em base mássica | 0,0794 | 0,8347 | 0,9475 | – | – | – |
| + ks(T) de Arrhenius | 0,0813 | 0,8269 | 0,9450 | – | – | – |
| + força motriz [H+] (V0: K2 = 0,0099, TCC) | 0,0443 | 0,9486 | 0,9837 | 0,052 | 0,865 | 1,075 |
| V1: K2(T) por van't Hoff | 0,0426 | 0,9524 | 0,9849 | 0,052 | 0,959 | 1,145 |
| V2: V1 + par iônico ZnSO4 | 0,0750 | 0,8525 | 0,9532 | 0,052 | 0,128 | 0,453 |
| V0 + nu | 0,0394 | 0,9593 | 0,9871 | 0,052 | 0,987 | 1,213 |
| V1 + nu | 0,0417 | 0,9544 | 0,9855 | 0,052 | 1,085 | 1,270 |
| V2 + nu | 0,0620 | 0,8992 | 0,9680 | 0,052 | 0,256 | 0,890 |

### T03. Seleção do White-Box por erro combinado de pH (ν a priori)

| Especiação | RMSE pH₀ | RMSE pH_f | viés pH_f | RMSE pH combinado | RMSE X (t>0) | R² X (t>0) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| V2 | 0,052 | 0,256 | 0,202 | 0,172 | 0,0620 | 0,8992 |
| V0 | 0,052 | 0,987 | 0,817 | 0,647 | 0,0394 | 0,9593 |
| V1 | 0,052 | 1,085 | 0,905 | 0,711 | 0,0417 | 0,9544 |

### T04. pH medido e previsto (especiação sobre o estado final medido; ν a priori)

| Ensaio | η | C_A0 (mol/L) | pH₀ med. | pH₀ prev. V0 | pH₀ prev. V2 | pH_f med. | pH_f prev. V0 | pH_f prev. V1 | pH_f prev. V2 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 1,0 | 0,97 | 0,97 | 6,0 | – | – | – |
| 2 | 3,1 | 0,1 | 1,0 | 0,97 | 0,97 | 1,2 | 1,29 | 1,31 | 1,20 |
| 3 | 0,5 | 0,5 | 0,2 | 0,29 | 0,29 | 5,4 | – | – | – |
| 4 | 0,5 | 1,0 | 0,0 | -0,00 | -0,00 | 5,3 | – | – | – |
| 5 | 0,5 | 1,5 | -0,2 | -0,18 | -0,18 | 5,0 | – | – | – |
| 6 | 1,0 | 0,1 | 1,0 | 0,97 | 0,97 | 2,0 | 2,50 | 2,62 | 2,08 |
| 7 | 3,1 | 0,5 | 0,3 | 0,29 | 0,29 | 0,5 | 0,72 | 0,74 | 0,55 |
| 8 | 1,0 | 0,5 | 0,2 | 0,29 | 0,29 | 1,2 | 2,44 | 2,57 | 1,59 |
| 9 | 1,0 | 1,0 | 0,0 | -0,00 | -0,00 | 0,9 | 2,35 | 2,49 | 1,27 |
| 10 | 1,0 | 1,5 | -0,2 | -0,18 | -0,18 | 0,6 | 2,35 | 2,48 | 1,12 |
| 11 | 1,5 | 0,1 | 1,0 | 0,97 | 0,97 | 1,5 | 1,90 | 1,99 | 1,59 |
| 12 | 3,1 | 1,0 | 0,1 | -0,00 | -0,00 | 0,2 | 0,44 | 0,44 | 0,25 |
| 13 | 1,5 | 0,5 | 0,2 | 0,29 | 0,29 | 0,7 | 1,68 | 1,80 | 0,95 |
| 14 | 1,5 | 1,0 | 0,0 | -0,00 | -0,00 | 0,5 | 1,65 | 1,78 | 0,67 |
| 15 | 1,5 | 1,5 | -0,2 | -0,18 | -0,18 | 0,2 | 1,63 | 1,76 | 0,50 |
| 16 | 3,1 | 1,5 | -0,2 | -0,18 | -0,18 | -0,1 | 0,27 | 0,28 | 0,07 |

### T05. Benchmark por protocolo de validação (RMSE em t>0; ganho = 1 − MSE/MSE_WB)

**GKF8 (partição do TCC)** (teto do ganho = 0,914)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido MLP tanh [B] | hibrido | 0,0185 | [0,0153; 0,0218] | – | 0,9910 | 0,904 | 99 | 7,1 | 5,4 |
| Híbrido GBDT [B] | hibrido | 0,0211 | [0,0175; 0,0251] | – | 0,9883 | 0,875 | 96 | 4,5 | 10,7 |
| Híbrido Random Forest [B] | hibrido | 0,0221 | [0,0175; 0,0273] | – | 0,9872 | 0,863 | 94 | 8,9 | 9,8 |
| Híbrido Extra Trees [B] | hibrido | 0,0226 | [0,0175; 0,0274] | – | 0,9866 | 0,857 | 94 | 5,4 | 3,6 |
| Híbrido GP [B] | hibrido | 0,0227 | [0,0167; 0,0282] | – | 0,9865 | 0,855 | 94 | 8,9 | 2,7 |
| Híbrido SVR [B] | hibrido | 0,0234 | [0,0191; 0,0274] | – | 0,9857 | 0,847 | 93 | 10,7 | 19,6 |
| Caixa-preta GP | caixa-preta | 0,0236 | [0,0197; 0,0276] | – | 0,9854 | 0,844 | 92 | 15,2 | 22,3 |
| Caixa-preta GBDT | caixa-preta | 0,0240 | [0,0193; 0,0286] | – | 0,9849 | 0,838 | 92 | 10,7 | 5,4 |
| Caixa-preta Extra Trees | caixa-preta | 0,0264 | [0,0229; 0,0297] | – | 0,9818 | 0,805 | 88 | 8,9 | 8,9 |
| Tabela (η, t) | baseline | 0,0303 | [0,0223; 0,0386] | – | 0,9759 | 0,742 | 81 | 6,2 | 17,0 |
| Ridge (resíduo) | baseline | 0,0384 | [0,0304; 0,0464] | – | 0,9613 | 0,585 | 64 | 0,0 | 0,0 |
| White-Box | white-box | 0,0597 | [0,0409; 0,0771] | – | 0,9067 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,0855 | [0,0632; 0,1082] | – | 0,8084 | -1,053 | -115 | 33,9 | 33,9 |

**GKF8 repetido (20 partições)** (teto do ganho = 0,913)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido MLP tanh [B] | hibrido | 0,0196 | [0,0168; 0,0225] | 0,0017 | 0,9898 | 0,890 | 97 | 7,3 | 8,3 |
| Híbrido GBDT [B] | hibrido | 0,0209 | [0,0174; 0,0250] | 0,0009 | 0,9885 | 0,875 | 96 | 6,5 | 8,0 |
| Híbrido Random Forest [B] | hibrido | 0,0219 | [0,0174; 0,0269] | 0,0009 | 0,9874 | 0,864 | 95 | 8,3 | 10,3 |
| Híbrido GP [B] | hibrido | 0,0224 | [0,0176; 0,0272] | 0,0012 | 0,9868 | 0,857 | 94 | 9,6 | 3,7 |
| Híbrido Extra Trees [B] | hibrido | 0,0230 | [0,0182; 0,0276] | 0,0007 | 0,9861 | 0,850 | 93 | 4,9 | 3,7 |
| Caixa-preta GP | caixa-preta | 0,0234 | [0,0197; 0,0273] | 0,0012 | 0,9856 | 0,844 | 92 | 16,0 | 23,8 |
| Híbrido SVR [B] | hibrido | 0,0236 | [0,0209; 0,0265] | 0,0034 | 0,9851 | 0,839 | 92 | 15,8 | 16,1 |
| Caixa-preta GBDT | caixa-preta | 0,0238 | [0,0193; 0,0280] | 0,0010 | 0,9851 | 0,839 | 92 | 12,4 | 5,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0270 | [0,0234; 0,0306] | 0,0012 | 0,9809 | 0,793 | 87 | 8,6 | 9,7 |
| Tabela (η, t) | baseline | 0,0306 | [0,0226; 0,0389] | 0,0009 | 0,9754 | 0,734 | 80 | 7,0 | 18,2 |
| Ridge (resíduo) | baseline | 0,0388 | [0,0307; 0,0468] | 0,0013 | 0,9605 | 0,572 | 63 | 0,0 | 0,0 |
| White-Box | white-box | 0,0594 | [0,0410; 0,0766] | 0,0004 | 0,9076 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,1006 | [0,0713; 0,1427] | 0,0312 | 0,7103 | -2,129 | -233 | 31,4 | 36,5 |

**LOGO (16 ensaios)** (teto do ganho = 0,914)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido MLP tanh [B] | hibrido | 0,0173 | [0,0147; 0,0199] | – | 0,9922 | 0,916 | 100 | 5,4 | 8,0 |
| Híbrido GP [B] | hibrido | 0,0208 | [0,0161; 0,0256] | – | 0,9887 | 0,879 | 96 | 9,8 | 1,8 |
| Híbrido GBDT [B] | hibrido | 0,0211 | [0,0173; 0,0252] | – | 0,9883 | 0,875 | 96 | 5,4 | 8,0 |
| Híbrido Random Forest [B] | hibrido | 0,0218 | [0,0173; 0,0269] | – | 0,9876 | 0,867 | 95 | 7,1 | 8,0 |
| Híbrido SVR [B] | hibrido | 0,0222 | [0,0186; 0,0257] | – | 0,9871 | 0,862 | 94 | 13,4 | 15,2 |
| Híbrido Extra Trees [B] | hibrido | 0,0224 | [0,0174; 0,0272] | – | 0,9868 | 0,859 | 94 | 3,6 | 2,7 |
| Caixa-preta GBDT | caixa-preta | 0,0236 | [0,0190; 0,0279] | – | 0,9854 | 0,844 | 92 | 11,6 | 4,5 |
| Caixa-preta GP | caixa-preta | 0,0236 | [0,0197; 0,0277] | – | 0,9854 | 0,844 | 92 | 16,1 | 25,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0264 | [0,0224; 0,0304] | – | 0,9818 | 0,805 | 88 | 7,1 | 8,9 |
| Tabela (η, t) | baseline | 0,0303 | [0,0223; 0,0386] | – | 0,9760 | 0,743 | 81 | 5,4 | 17,0 |
| Ridge (resíduo) | baseline | 0,0389 | [0,0306; 0,0472] | – | 0,9603 | 0,575 | 63 | 0,0 | 0,0 |
| White-Box | white-box | 0,0597 | [0,0410; 0,0772] | – | 0,9065 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,0792 | [0,0539; 0,1043] | – | 0,8357 | -0,757 | -83 | 32,1 | 32,1 |

**Deixa-um-η-fora (extrapolação em η)** (teto do ganho = 0,913)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido Random Forest [B] | hibrido | 0,0417 | [0,0305; 0,0523] | – | 0,9545 | 0,510 | 56 | 8,0 | 3,6 |
| Híbrido GBDT [B] | hibrido | 0,0432 | [0,0312; 0,0544] | – | 0,9510 | 0,472 | 52 | 14,3 | 5,4 |
| Híbrido Extra Trees [B] | hibrido | 0,0450 | [0,0348; 0,0551] | – | 0,9470 | 0,428 | 47 | 16,1 | 10,7 |
| Híbrido MLP tanh [B] | hibrido | 0,0467 | [0,0379; 0,0547] | – | 0,9429 | 0,385 | 42 | 40,2 | 30,4 |
| Ridge (resíduo) | baseline | 0,0507 | [0,0396; 0,0610] | – | 0,9326 | 0,273 | 30 | 0,0 | 0,0 |
| Tabela (η, t) | baseline | 0,0527 | [0,0394; 0,0652] | – | 0,9273 | 0,216 | 24 | 0,0 | 22,3 |
| Híbrido GP [B] | hibrido | 0,0536 | [0,0430; 0,0643] | – | 0,9247 | 0,188 | 21 | 0,0 | 9,8 |
| Híbrido SVR [B] | hibrido | 0,0582 | [0,0433; 0,0716] | – | 0,9114 | 0,044 | 5 | 5,4 | 33,9 |
| White-Box | white-box | 0,0595 | [0,0406; 0,0770] | – | 0,9073 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0727 | [0,0549; 0,0885] | – | 0,8614 | -0,494 | -54 | 25,0 | 8,9 |
| Caixa-preta GBDT | caixa-preta | 0,0758 | [0,0589; 0,0898] | – | 0,8492 | -0,626 | -69 | 25,0 | 1,8 |
| Caixa-preta GP | caixa-preta | 0,0950 | [0,0710; 0,1161] | – | 0,7635 | -1,550 | -170 | 47,3 | 17,0 |
| Caixa-preta SVR | caixa-preta | 0,4914 | [0,2519; 0,6830] | – | -5,3298 | -67,250 | -7366 | 42,9 | 25,0 |

**Deixa-um-C_A0-fora** (teto do ganho = 0,914)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GP [B] | hibrido | 0,0225 | [0,0173; 0,0278] | – | 0,9868 | 0,858 | 94 | 8,9 | 1,8 |
| Híbrido Extra Trees [B] | hibrido | 0,0233 | [0,0177; 0,0286] | – | 0,9858 | 0,848 | 93 | 2,7 | 1,8 |
| Híbrido GBDT [B] | hibrido | 0,0239 | [0,0175; 0,0314] | – | 0,9850 | 0,840 | 92 | 2,7 | 10,7 |
| Híbrido MLP tanh [B] | hibrido | 0,0240 | [0,0183; 0,0303] | – | 0,9849 | 0,838 | 92 | 14,3 | 12,5 |
| Híbrido Random Forest [B] | hibrido | 0,0248 | [0,0173; 0,0330] | – | 0,9839 | 0,828 | 91 | 4,5 | 8,9 |
| Caixa-preta Extra Trees | caixa-preta | 0,0274 | [0,0233; 0,0316] | – | 0,9803 | 0,789 | 86 | 7,1 | 10,7 |
| Caixa-preta GP | caixa-preta | 0,0274 | [0,0208; 0,0339] | – | 0,9803 | 0,789 | 86 | 16,1 | 25,0 |
| Híbrido SVR [B] | hibrido | 0,0287 | [0,0201; 0,0368] | – | 0,9784 | 0,769 | 84 | 9,8 | 10,7 |
| Caixa-preta GBDT | caixa-preta | 0,0288 | [0,0217; 0,0352] | – | 0,9783 | 0,768 | 84 | 10,7 | 3,6 |
| Tabela (η, t) | baseline | 0,0303 | [0,0223; 0,0386] | – | 0,9760 | 0,743 | 81 | 5,4 | 17,0 |
| Ridge (resíduo) | baseline | 0,0355 | [0,0277; 0,0428] | – | 0,9670 | 0,647 | 71 | 0,0 | 0,0 |
| White-Box | white-box | 0,0597 | [0,0410; 0,0772] | – | 0,9065 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,1813 | [0,1141; 0,2369] | – | 0,1384 | -8,213 | -899 | 19,6 | 31,2 |

**Extrapolação temporal (t<=5 -> 15 min)** (teto do ganho = 0,928)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido Random Forest [B] | hibrido | 0,0108 | [0,0079; 0,0132] | – | 0,9971 | 0,973 | 105 | 0,0 | 0,0 |
| Híbrido GBDT [B] | hibrido | 0,0112 | [0,0076; 0,0146] | – | 0,9969 | 0,970 | 105 | 0,0 | 0,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0121 | [0,0095; 0,0147] | – | 0,9964 | 0,966 | 104 | 0,0 | 0,0 |
| Híbrido Extra Trees [B] | hibrido | 0,0122 | [0,0093; 0,0148] | – | 0,9963 | 0,965 | 104 | 0,0 | 0,0 |
| Híbrido GP [B] | hibrido | 0,0177 | [0,0122; 0,0228] | – | 0,9922 | 0,927 | 100 | 0,0 | 0,0 |
| Híbrido MLP tanh [B] | hibrido | 0,0202 | [0,0122; 0,0272] | – | 0,9899 | 0,904 | 97 | 12,5 | 0,0 |
| Caixa-preta GBDT | caixa-preta | 0,0213 | [0,0139; 0,0287] | – | 0,9888 | 0,894 | 96 | 6,2 | 0,0 |
| Ridge (resíduo) | baseline | 0,0559 | [0,0275; 0,0777] | – | 0,9227 | 0,270 | 29 | 50,0 | 0,0 |
| Híbrido SVR [B] | hibrido | 0,0571 | [0,0413; 0,0715] | – | 0,9193 | 0,238 | 26 | 12,5 | 0,0 |
| White-Box | white-box | 0,0654 | [0,0327; 0,0902] | – | 0,8941 | 0,000 | 0 | 0,0 | 0,0 |
| Tabela (η, t) | baseline | 0,0654 | [0,0327; 0,0902] | – | 0,8941 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta GP | caixa-preta | 0,2000 | [0,1543; 0,2391] | – | 0,0100 | -8,349 | -900 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,2892 | [0,2295; 0,3430] | – | -1,0700 | -18,548 | -1999 | 12,5 | 0,0 |

### T06. Extrapolação em η: RMSE no nível deixado fora do treino

| Modelo | RMSE total | η = 0,5 ausente | η = 1,0 ausente | η = 1,5 ausente | η = 3,1 ausente |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Híbrido Random Forest [B] | 0,0417 | 0,0228 | 0,0649 | 0,0326 | 0,0338 |
| Híbrido GBDT [B] | 0,0432 | 0,0236 | 0,0673 | 0,0353 | 0,0338 |
| Híbrido Extra Trees [B] | 0,0450 | 0,0237 | 0,0672 | 0,0428 | 0,0344 |
| Híbrido MLP tanh [B] | 0,0467 | 0,0289 | 0,0630 | 0,0464 | 0,0418 |
| Ridge (resíduo) | 0,0507 | 0,0541 | 0,0713 | 0,0276 | 0,0391 |
| Tabela (η, t) | 0,0527 | 0,0473 | 0,0815 | 0,0265 | 0,0387 |
| Híbrido GP [B] | 0,0536 | 0,0417 | 0,0775 | 0,0518 | 0,0326 |
| Híbrido SVR [B] | 0,0582 | 0,0790 | 0,0664 | 0,0463 | 0,0270 |
| White-Box | 0,0595 | 0,0211 | 0,0997 | 0,0504 | 0,0351 |
| Caixa-preta Extra Trees | 0,0727 | 0,0237 | 0,1129 | 0,0782 | 0,0415 |
| Caixa-preta GBDT | 0,0758 | 0,0237 | 0,0953 | 0,1078 | 0,0419 |
| Caixa-preta GP | 0,0950 | 0,0237 | 0,1034 | 0,0610 | 0,1453 |
| Caixa-preta SVR | 0,4914 | 0,0215 | 0,1207 | 0,0696 | 0,9727 |

### T07. Híbrido contra a caixa-preta do MESMO modelo (diferença pareada, bootstrap de ensaios)

| Protocolo | Híbrido | Caixa-preta | Δ RMSE (híbrido − caixa-preta) | IC95% | Empate (IC contém 0) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | Híbrido GP [B] | Caixa-preta GP | -0,0009 | [-0,0072; 0,0044] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0038 | [-0,0076; 0,0004] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0029 | [-0,0082; 0,0020] | sim |
| GKF8 (partição do TCC) | Híbrido SVR [B] | Caixa-preta SVR | -0,0621 | [-0,0860; -0,0400] | não |
| GKF8 repetido (20 partições) | Híbrido GP [B] | Caixa-preta GP | -0,0010 | [-0,0057; 0,0033] | sim |
| GKF8 repetido (20 partições) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0040 | [-0,0076; -0,0002] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0029 | [-0,0079; 0,0016] | sim |
| GKF8 repetido (20 partições) | Híbrido SVR [B] | Caixa-preta SVR | -0,0813 | [-0,1196; -0,0479] | não |
| LOGO (16 ensaios) | Híbrido GP [B] | Caixa-preta GP | -0,0028 | [-0,0084; 0,0019] | sim |
| LOGO (16 ensaios) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0039 | [-0,0074; -0,0000] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0025 | [-0,0077; 0,0020] | sim |
| LOGO (16 ensaios) | Híbrido SVR [B] | Caixa-preta SVR | -0,0570 | [-0,0828; -0,0319] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido GP [B] | Caixa-preta GP | -0,0414 | [-0,0666; -0,0141] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0277 | [-0,0387; -0,0162] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0326 | [-0,0485; -0,0179] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido SVR [B] | Caixa-preta SVR | -0,4333 | [-0,6297; -0,1862] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GP | -0,0049 | [-0,0118; 0,0014] | sim |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0041 | [-0,0081; -0,0001] | não |
| Deixa-um-C_A0-fora | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0049 | [-0,0115; 0,0012] | sim |
| Deixa-um-C_A0-fora | Híbrido SVR [B] | Caixa-preta SVR | -0,1526 | [-0,2074; -0,0884] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GP [B] | Caixa-preta GP | -0,1823 | [-0,2176; -0,1393] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | 0,0001 | [-0,0021; 0,0025] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0100 | [-0,0160; -0,0039] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido SVR [B] | Caixa-preta SVR | -0,2321 | [-0,2825; -0,1781] | não |

### T08. Melhor híbrido de cada protocolo contra os demais modelos (Δ < 0: o melhor híbrido erra menos)

| Protocolo | Melhor híbrido | Contra | Δ RMSE | IC95% | Empate |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | White-Box | -0,0411 | [-0,0602; -0,0217] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Tabela (η, t) | -0,0118 | [-0,0189; -0,0052] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Ridge (resíduo) | -0,0199 | [-0,0289; -0,0116] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Híbrido GP [B] | -0,0042 | [-0,0078; -0,0002] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Híbrido Extra Trees [B] | -0,0041 | [-0,0073; -0,0007] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Híbrido GBDT [B] | -0,0026 | [-0,0050; -0,0001] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Híbrido Random Forest [B] | -0,0036 | [-0,0063; -0,0005] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Híbrido SVR [B] | -0,0048 | [-0,0086; -0,0014] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Caixa-preta GP | -0,0050 | [-0,0091; -0,0014] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Caixa-preta Extra Trees | -0,0079 | [-0,0120; -0,0039] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Caixa-preta GBDT | -0,0055 | [-0,0101; -0,0011] | não |
| GKF8 (partição do TCC) | Híbrido MLP tanh [B] | Caixa-preta SVR | -0,0670 | [-0,0902; -0,0445] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | White-Box | -0,0397 | [-0,0584; -0,0204] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Tabela (η, t) | -0,0109 | [-0,0183; -0,0040] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Ridge (resíduo) | -0,0191 | [-0,0283; -0,0106] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Híbrido GP [B] | -0,0027 | [-0,0057; 0,0005] | sim |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Híbrido Extra Trees [B] | -0,0033 | [-0,0056; -0,0008] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Híbrido GBDT [B] | -0,0013 | [-0,0039; 0,0017] | sim |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Híbrido Random Forest [B] | -0,0022 | [-0,0057; 0,0016] | sim |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Híbrido SVR [B] | -0,0041 | [-0,0064; -0,0023] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Caixa-preta GP | -0,0037 | [-0,0063; -0,0010] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Caixa-preta Extra Trees | -0,0073 | [-0,0103; -0,0038] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Caixa-preta GBDT | -0,0041 | [-0,0071; -0,0009] | não |
| GKF8 repetido (20 partições) | Híbrido MLP tanh [B] | Caixa-preta SVR | -0,0854 | [-0,1238; -0,0517] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | White-Box | -0,0425 | [-0,0609; -0,0233] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Tabela (η, t) | -0,0130 | [-0,0203; -0,0060] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Ridge (resíduo) | -0,0217 | [-0,0309; -0,0134] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Híbrido GP [B] | -0,0035 | [-0,0065; -0,0005] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Híbrido Extra Trees [B] | -0,0052 | [-0,0091; -0,0014] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Híbrido GBDT [B] | -0,0038 | [-0,0068; -0,0011] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Híbrido Random Forest [B] | -0,0045 | [-0,0084; -0,0011] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Híbrido SVR [B] | -0,0049 | [-0,0081; -0,0021] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Caixa-preta GP | -0,0064 | [-0,0105; -0,0023] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Caixa-preta Extra Trees | -0,0091 | [-0,0137; -0,0045] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Caixa-preta GBDT | -0,0063 | [-0,0108; -0,0020] | não |
| LOGO (16 ensaios) | Híbrido MLP tanh [B] | Caixa-preta SVR | -0,0619 | [-0,0876; -0,0365] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | White-Box | -0,0178 | [-0,0252; -0,0080] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Tabela (η, t) | -0,0110 | [-0,0156; -0,0057] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Ridge (resíduo) | -0,0091 | [-0,0162; -0,0034] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido GP [B] | -0,0120 | [-0,0175; -0,0071] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido Extra Trees [B] | -0,0033 | [-0,0071; -0,0008] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido GBDT [B] | -0,0016 | [-0,0027; -0,0004] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido SVR [B] | -0,0165 | [-0,0329; -0,0021] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido MLP tanh [B] | -0,0050 | [-0,0126; 0,0011] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta GP | -0,0533 | [-0,0752; -0,0296] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta Extra Trees | -0,0311 | [-0,0428; -0,0185] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta GBDT | -0,0342 | [-0,0496; -0,0194] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta SVR | -0,4498 | [-0,6432; -0,2087] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | White-Box | -0,0373 | [-0,0569; -0,0177] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Tabela (η, t) | -0,0078 | [-0,0145; -0,0019] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Ridge (resíduo) | -0,0130 | [-0,0221; -0,0053] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido Extra Trees [B] | -0,0008 | [-0,0039; 0,0018] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido GBDT [B] | -0,0014 | [-0,0051; 0,0023] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido Random Forest [B] | -0,0023 | [-0,0067; 0,0023] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido SVR [B] | -0,0062 | [-0,0135; 0,0011] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido MLP tanh [B] | -0,0016 | [-0,0042; 0,0012] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GP | -0,0049 | [-0,0118; 0,0014] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta Extra Trees | -0,0049 | [-0,0107; 0,0006] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GBDT | -0,0063 | [-0,0127; 0,0001] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta SVR | -0,1588 | [-0,2166; -0,0905] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | White-Box | -0,0547 | [-0,0792; -0,0228] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Tabela (η, t) | -0,0547 | [-0,0792; -0,0228] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Ridge (resíduo) | -0,0451 | [-0,0667; -0,0171] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Híbrido GP [B] | -0,0070 | [-0,0116; -0,0023] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Híbrido Extra Trees [B] | -0,0015 | [-0,0029; -0,0001] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Híbrido GBDT [B] | -0,0005 | [-0,0021; 0,0015] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Híbrido SVR [B] | -0,0464 | [-0,0606; -0,0308] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Híbrido MLP tanh [B] | -0,0095 | [-0,0175; -0,0007] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Caixa-preta GP | -0,1893 | [-0,2281; -0,1439] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Caixa-preta Extra Trees | -0,0014 | [-0,0040; 0,0012] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Caixa-preta GBDT | -0,0105 | [-0,0173; -0,0036] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Random Forest [B] | Caixa-preta SVR | -0,2785 | [-0,3330; -0,2184] | não |

### T09. Curva de aprendizado: RMSE em 2 ensaios fora do treino, por nº de ensaios de treino (30 sorteios)

| Modelo | n = 4 | n = 6 | n = 8 | n = 10 | n = 12 | n = 14 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| White-Box | 0,0581 | 0,0571 | 0,0568 | 0,0576 | 0,0589 | 0,0591 |
| Tabela (η, t) | 0,0466 | 0,0367 | 0,0365 | 0,0331 | 0,0309 | 0,0301 |
| Híbrido Extra Trees [B] | 0,0421 | 0,0323 | 0,0288 | 0,0259 | 0,0225 | 0,0220 |
| Híbrido GBDT [B] | 0,0390 | 0,0284 | 0,0280 | 0,0237 | 0,0213 | 0,0217 |
| Híbrido GP [B] | 0,0462 | 0,0358 | 0,0315 | 0,0278 | 0,0221 | 0,0214 |
| Caixa-preta Extra Trees | 0,0876 | 0,0445 | 0,0357 | 0,0289 | 0,0283 | 0,0267 |
| Caixa-preta GBDT | 0,0990 | 0,0413 | 0,0311 | 0,0264 | 0,0252 | 0,0239 |
| Caixa-preta GP | 0,0861 | 0,0524 | 0,0377 | 0,0270 | 0,0248 | 0,0232 |

### T10. Ablação dos atributos do ML (RMSE em t>0)

| Protocolo | Modelo | A (t, η, C_A0) | B (A + X_PBM) | C (B + C_Af, [H⁺]) |
| :--- | ---: | ---: | ---: | ---: |
| GKF8 repetido (10) | GP | 0,0251 | 0,0226 | 0,0235 |
| GKF8 repetido (10) | Extra Trees | 0,0266 | 0,0229 | 0,0214 |
| GKF8 repetido (10) | GBDT | 0,0324 | 0,0210 | 0,0212 |
| GKF8 repetido (10) | SVR | 0,0450 | 0,0239 | 0,0235 |
| Deixa-um-η-fora | GP | 0,0616 | 0,0536 | 0,0511 |
| Deixa-um-η-fora | Extra Trees | 0,0638 | 0,0450 | 0,0447 |
| Deixa-um-η-fora | GBDT | 0,0634 | 0,0432 | 0,0467 |
| Deixa-um-η-fora | SVR | 0,0849 | 0,0582 | 0,0545 |

### T11. Sensibilidade ao fator ν (RMSE em t>0)

| Protocolo | Modelo | nu estimado na dobra | nu a priori (1.06) | nu = 1 (sem consumo extra de ácido) |
| :--- | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | White-Box | 0,0597 | 0,0620 | 0,0750 |
| GKF8 (partição do TCC) | Híbrido GP [B] | 0,0227 | 0,0213 | 0,0220 |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | 0,0226 | 0,0226 | 0,0218 |
| Deixa-um-η-fora | White-Box | 0,0595 | 0,0620 | 0,0750 |
| Deixa-um-η-fora | Híbrido GP [B] | 0,0536 | 0,0557 | 0,0522 |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | 0,0450 | 0,0471 | 0,0626 |

### T12. O mesmo ML sobre degraus do White-Box (RMSE em t>0)

| Protocolo | Modelo | 1. Tese (α = 0) | 2. White-Box do TCC (V0, nu = 1) | 3. V0 + nu | 4. V2 + nu (principal) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 repetido (10) | White-Box | 0,0981 | 0,0443 | 0,0393 | 0,0593 |
| GKF8 repetido (10) | Híbrido GP [B] | 0,0437 | 0,0235 | 0,0237 | 0,0226 |
| GKF8 repetido (10) | Híbrido Extra Trees [B] | 0,0346 | 0,0200 | 0,0196 | 0,0229 |
| GKF8 repetido (10) | Caixa-preta GP | 0,0237 | 0,0237 | 0,0231 | 0,0231 |
| GKF8 repetido (10) | Caixa-preta Extra Trees | 0,0270 | 0,0270 | 0,0269 | 0,0269 |
| Deixa-um-η-fora | White-Box | 0,0981 | 0,0443 | 0,0394 | 0,0595 |
| Deixa-um-η-fora | Híbrido GP [B] | 0,0837 | 0,0455 | 0,0392 | 0,0536 |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | 0,0661 | 0,0341 | 0,0290 | 0,0450 |
| Deixa-um-η-fora | Caixa-preta GP | 0,0978 | 0,0978 | 0,0950 | 0,0950 |
| Deixa-um-η-fora | Caixa-preta Extra Trees | 0,0764 | 0,0764 | 0,0727 | 0,0727 |
| Poucos dados (n = 4 ensaios de treino) | White-Box | 0,0904 | 0,0450 | 0,0403 | 0,0571 |
| Poucos dados (n = 4 ensaios de treino) | Híbrido GP [B] | 0,0668 | 0,0335 | 0,0315 | 0,0467 |
| Poucos dados (n = 4 ensaios de treino) | Híbrido Extra Trees [B] | 0,0639 | 0,0351 | 0,0296 | 0,0431 |
| Poucos dados (n = 4 ensaios de treino) | Caixa-preta GP | 0,0650 | 0,0650 | 0,0606 | 0,0606 |
| Poucos dados (n = 4 ensaios de treino) | Caixa-preta Extra Trees | 0,0873 | 0,0873 | 0,0822 | 0,0822 |
| Poucos dados (n = 8 ensaios de treino) | White-Box | 0,0904 | 0,0450 | 0,0395 | 0,0553 |
| Poucos dados (n = 8 ensaios de treino) | Híbrido GP [B] | 0,0555 | 0,0279 | 0,0261 | 0,0291 |
| Poucos dados (n = 8 ensaios de treino) | Híbrido Extra Trees [B] | 0,0478 | 0,0225 | 0,0207 | 0,0271 |
| Poucos dados (n = 8 ensaios de treino) | Caixa-preta GP | 0,0276 | 0,0276 | 0,0273 | 0,0273 |
| Poucos dados (n = 8 ensaios de treino) | Caixa-preta Extra Trees | 0,0344 | 0,0344 | 0,0310 | 0,0310 |

### T13. Cobertura empírica dos intervalos de 95% (t>0, pontos fora-da-dobra)

| Protocolo | Modelo | Método | Cobertura | Largura média | n |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | Conforme (LOGO interno) | 62% | 0,0790 | 112 |
| Deixa-um-η-fora | Híbrido GBDT [B] | Conforme (LOGO interno) | 66% | 0,0732 | 112 |
| Deixa-um-η-fora | Híbrido GP [B] | Conforme (LOGO interno) | 49% | 0,0896 | 112 |
| Deixa-um-η-fora | Híbrido GP [B] | GP (analítico) | 71% | 0,1376 | 112 |
| LOGO (16 ensaios) | Híbrido Extra Trees [B] | Conforme (LOGO interno) | 90% | 0,0737 | 112 |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Conforme (LOGO interno) | 90% | 0,0752 | 112 |
| LOGO (16 ensaios) | Híbrido GP [B] | Conforme (LOGO interno) | 91% | 0,0820 | 112 |
| LOGO (16 ensaios) | Híbrido GP [B] | GP (analítico) | 89% | 0,0615 | 112 |

### T14. Resíduo do White-Box principal aos 15 min

| Ensaio | η | η real | C_A0 | X_exp (15 min) | X_PBM (15 min) | ΔX |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,535 | 0,1 | 0,50 | 0,5044 | -0,004 |
| 2 | 3,1 | 3,290 | 0,1 | 1,00 | 1,0000 | 0,000 |
| 3 | 0,5 | 0,535 | 0,5 | 0,50 | 0,5044 | -0,004 |
| 4 | 0,5 | 0,535 | 1,0 | 0,50 | 0,5044 | -0,004 |
| 5 | 0,5 | 0,535 | 1,5 | 0,49 | 0,5044 | -0,014 |
| 6 | 1,0 | 1,069 | 0,1 | 0,87 | 0,9719 | -0,102 |
| 7 | 3,1 | 3,192 | 0,5 | 1,01 | 1,0000 | 0,010 |
| 8 | 1,0 | 1,069 | 0,5 | 0,87 | 0,9883 | -0,118 |
| 9 | 1,0 | 1,069 | 1,0 | 0,85 | 0,9918 | -0,142 |
| 10 | 1,0 | 1,069 | 1,5 | 0,85 | 0,9933 | -0,143 |
| 11 | 1,5 | 1,584 | 0,1 | 0,98 | 1,0000 | -0,020 |
| 12 | 3,1 | 3,216 | 1,0 | 1,00 | 1,0000 | 0,000 |
| 13 | 1,5 | 1,608 | 0,5 | 0,97 | 1,0000 | -0,030 |
| 14 | 1,5 | 1,602 | 1,0 | 0,97 | 1,0000 | -0,030 |
| 15 | 1,5 | 1,604 | 1,5 | 0,97 | 1,0000 | -0,030 |
| 16 | 3,1 | 3,208 | 1,5 | 1,00 | 1,0000 | 0,000 |

### T15. Escalas de comprimento do GP (ARD) ajustadas com todos os dados

| Atributos | Entrada | Escala de comprimento (σ das entradas) | Ruído branco (norm.) |
| :--- | ---: | ---: | ---: |
| A | tempo_min | 0,668 | 0,078 |
| A | razao_molar_eta_real | 0,044 | 0,078 |
| A | C_acid_0_mol_L | 3,114 | 0,078 |
| B | tempo_min | 4,200 | 0,102 |
| B | razao_molar_eta_real | 0,577 | 0,102 |
| B | C_acid_0_mol_L | 6,079 | 0,102 |
| B | X_pbm | 0,364 | 0,102 |
| C | tempo_min | 0,579 | 0,000 |
| C | razao_molar_eta_real | 100000,000 | 0,000 |
| C | C_acid_0_mol_L | 1,617 | 0,000 |
| C | X_pbm | 0,103 | 0,000 |
| C | C_acid_pbm_mol_L | 0,763 | 0,000 |
| C | H_pbm_mol_L | 3288,253 | 0,000 |

### T16. Conversão do White-Box principal e resíduo, por ensaio

**X_PBM (conversão do White-Box principal)**

| Ensaio | η | C_A0 | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 0,00 | 0,48 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 2 | 3,1 | 0,1 | 0,00 | 0,88 | 0,96 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 |
| 3 | 0,5 | 0,5 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 4 | 0,5 | 1,0 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 5 | 0,5 | 1,5 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 6 | 1,0 | 0,1 | 0,00 | 0,72 | 0,81 | 0,87 | 0,90 | 0,92 | 0,93 | 0,97 |
| 7 | 3,1 | 0,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 8 | 1,0 | 0,5 | 0,00 | 0,86 | 0,91 | 0,94 | 0,95 | 0,96 | 0,97 | 0,99 |
| 9 | 1,0 | 1,0 | 0,00 | 0,90 | 0,93 | 0,96 | 0,97 | 0,97 | 0,98 | 0,99 |
| 10 | 1,0 | 1,5 | 0,00 | 0,92 | 0,94 | 0,96 | 0,97 | 0,98 | 0,98 | 0,99 |
| 11 | 1,5 | 0,1 | 0,00 | 0,80 | 0,90 | 0,96 | 0,98 | 0,99 | 1,00 | 1,00 |
| 12 | 3,1 | 1,0 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 13 | 1,5 | 0,5 | 0,00 | 0,97 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 14 | 1,5 | 1,0 | 0,00 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 15 | 1,5 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 16 | 3,1 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |

**ΔX = X_exp − X_PBM**

| Ensaio | η | C_A0 | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 0,00 | -0,06 | -0,04 | -0,04 | -0,01 | 0,01 | -0,01 | -0,00 |
| 2 | 3,1 | 0,1 | 0,00 | -0,06 | -0,06 | -0,04 | -0,03 | -0,03 | -0,01 | 0,00 |
| 3 | 0,5 | 0,5 | 0,00 | -0,05 | -0,02 | -0,01 | 0,02 | -0,00 | -0,01 | -0,00 |
| 4 | 0,5 | 1,0 | 0,00 | -0,00 | -0,03 | 0,01 | -0,03 | 0,03 | 0,01 | -0,00 |
| 5 | 0,5 | 1,5 | 0,00 | -0,01 | -0,01 | -0,03 | -0,01 | -0,01 | -0,03 | -0,01 |
| 6 | 1,0 | 0,1 | 0,00 | 0,03 | -0,03 | -0,07 | -0,05 | -0,05 | -0,08 | -0,10 |
| 7 | 3,1 | 0,5 | 0,00 | -0,12 | -0,07 | -0,02 | -0,01 | 0,00 | -0,01 | 0,01 |
| 8 | 1,0 | 0,5 | 0,00 | -0,03 | -0,09 | -0,09 | -0,09 | -0,09 | -0,11 | -0,12 |
| 9 | 1,0 | 1,0 | 0,00 | -0,07 | -0,09 | -0,11 | -0,12 | -0,14 | -0,14 | -0,14 |
| 10 | 1,0 | 1,5 | 0,00 | -0,08 | -0,11 | -0,13 | -0,14 | -0,15 | -0,14 | -0,14 |
| 11 | 1,5 | 0,1 | 0,00 | 0,06 | -0,05 | -0,08 | -0,01 | -0,05 | -0,04 | -0,02 |
| 12 | 3,1 | 1,0 | 0,00 | -0,04 | 0,00 | 0,00 | -0,01 | -0,01 | -0,01 | 0,00 |
| 13 | 1,5 | 0,5 | 0,00 | -0,05 | -0,07 | -0,05 | -0,05 | -0,05 | -0,05 | -0,03 |
| 14 | 1,5 | 1,0 | 0,00 | -0,08 | -0,04 | -0,05 | -0,03 | -0,05 | -0,02 | -0,03 |
| 15 | 1,5 | 1,5 | 0,00 | -0,09 | -0,06 | -0,05 | -0,06 | -0,04 | -0,04 | -0,03 |
| 16 | 3,1 | 1,5 | 0,00 | -0,02 | -0,02 | -0,04 | -0,01 | -0,01 | 0,00 | 0,00 |
