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
| Laboratório: CAf, nu = 1 (principal) | 0,0813 | 0,8269 | 0,9450 | – | – | – |

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

**GKF8 (partição do TCC)** (teto do ganho = 0,953)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido Extra Trees [B] | hibrido | 0,0215 | [0,0163; 0,0263] | – | 0,9879 | 0,930 | 98 | 1,8 | 0,9 |
| Híbrido GBDT [B] | hibrido | 0,0215 | [0,0173; 0,0262] | – | 0,9879 | 0,930 | 98 | 0,9 | 1,8 |
| Híbrido SVR [B] | hibrido | 0,0218 | [0,0180; 0,0258] | – | 0,9876 | 0,928 | 97 | 7,1 | 12,5 |
| Híbrido GP [B] | hibrido | 0,0221 | [0,0176; 0,0265] | – | 0,9872 | 0,926 | 97 | 4,5 | 4,5 |
| Híbrido Random Forest [B] | hibrido | 0,0226 | [0,0164; 0,0288] | – | 0,9866 | 0,922 | 97 | 0,0 | 6,2 |
| Híbrido MLP tanh [B] | hibrido | 0,0231 | [0,0181; 0,0279] | – | 0,9860 | 0,919 | 96 | 6,2 | 2,7 |
| Caixa-preta GP | caixa-preta | 0,0241 | [0,0205; 0,0279] | – | 0,9848 | 0,912 | 96 | 5,4 | 22,3 |
| Caixa-preta GBDT | caixa-preta | 0,0243 | [0,0197; 0,0287] | – | 0,9845 | 0,911 | 96 | 0,0 | 5,4 |
| Caixa-preta Extra Trees | caixa-preta | 0,0266 | [0,0232; 0,0298] | – | 0,9815 | 0,893 | 94 | 2,7 | 8,9 |
| Tabela (η, t) | baseline | 0,0296 | [0,0219; 0,0381] | – | 0,9771 | 0,868 | 91 | 2,7 | 18,8 |
| Ridge (resíduo) | baseline | 0,0424 | [0,0333; 0,0514] | – | 0,9529 | 0,728 | 76 | 0,0 | 0,0 |
| White-Box | white-box | 0,0813 | [0,0552; 0,1047] | – | 0,8269 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,0873 | [0,0670; 0,1090] | – | 0,8002 | -0,154 | -16 | 29,5 | 33,9 |

**GKF8 repetido (20 partições)** (teto do ganho = 0,953)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GBDT [B] | hibrido | 0,0220 | [0,0179; 0,0265] | 0,0008 | 0,9873 | 0,926 | 97 | 0,5 | 3,3 |
| Híbrido GP [B] | hibrido | 0,0223 | [0,0180; 0,0269] | 0,0010 | 0,9869 | 0,924 | 97 | 2,1 | 4,6 |
| Híbrido Extra Trees [B] | hibrido | 0,0228 | [0,0178; 0,0275] | 0,0013 | 0,9863 | 0,921 | 97 | 2,4 | 2,4 |
| Híbrido Random Forest [B] | hibrido | 0,0229 | [0,0167; 0,0291] | 0,0006 | 0,9862 | 0,920 | 97 | 0,0 | 6,9 |
| Híbrido MLP tanh [B] | hibrido | 0,0233 | [0,0197; 0,0270] | 0,0019 | 0,9857 | 0,917 | 96 | 5,7 | 5,0 |
| Híbrido SVR [B] | hibrido | 0,0235 | [0,0207; 0,0266] | 0,0022 | 0,9854 | 0,916 | 96 | 7,3 | 11,2 |
| Caixa-preta GP | caixa-preta | 0,0239 | [0,0204; 0,0275] | 0,0012 | 0,9850 | 0,913 | 96 | 5,1 | 23,8 |
| Caixa-preta GBDT | caixa-preta | 0,0241 | [0,0197; 0,0281] | 0,0010 | 0,9847 | 0,912 | 96 | 0,0 | 5,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0271 | [0,0236; 0,0306] | 0,0012 | 0,9807 | 0,888 | 93 | 2,5 | 9,7 |
| Tabela (η, t) | baseline | 0,0299 | [0,0223; 0,0384] | 0,0007 | 0,9766 | 0,865 | 91 | 2,5 | 19,4 |
| Ridge (resíduo) | baseline | 0,0430 | [0,0336; 0,0523] | 0,0014 | 0,9515 | 0,720 | 76 | 0,6 | 0,0 |
| White-Box | white-box | 0,0813 | [0,0552; 0,1047] | 0,0000 | 0,8269 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,1018 | [0,0736; 0,1433] | 0,0309 | 0,7044 | -0,708 | -74 | 23,8 | 36,5 |

**LOGO (16 ensaios)** (teto do ganho = 0,953)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GBDT [B] | hibrido | 0,0215 | [0,0170; 0,0260] | – | 0,9879 | 0,930 | 98 | 0,9 | 2,7 |
| Híbrido GP [B] | hibrido | 0,0216 | [0,0169; 0,0261] | – | 0,9878 | 0,930 | 98 | 3,6 | 3,6 |
| Híbrido Extra Trees [B] | hibrido | 0,0218 | [0,0166; 0,0265] | – | 0,9876 | 0,928 | 97 | 2,7 | 0,9 |
| Híbrido Random Forest [B] | hibrido | 0,0228 | [0,0159; 0,0292] | – | 0,9864 | 0,921 | 97 | 0,0 | 5,4 |
| Caixa-preta GBDT | caixa-preta | 0,0239 | [0,0195; 0,0281] | – | 0,9850 | 0,914 | 96 | 0,0 | 4,5 |
| Híbrido SVR [B] | hibrido | 0,0241 | [0,0199; 0,0283] | – | 0,9847 | 0,912 | 96 | 6,2 | 10,7 |
| Caixa-preta GP | caixa-preta | 0,0242 | [0,0205; 0,0281] | – | 0,9847 | 0,912 | 96 | 4,5 | 25,0 |
| Híbrido MLP tanh [B] | hibrido | 0,0248 | [0,0197; 0,0297] | – | 0,9839 | 0,907 | 95 | 5,4 | 4,5 |
| Caixa-preta Extra Trees | caixa-preta | 0,0265 | [0,0226; 0,0305] | – | 0,9816 | 0,893 | 94 | 2,7 | 8,9 |
| Tabela (η, t) | baseline | 0,0296 | [0,0219; 0,0381] | – | 0,9771 | 0,868 | 91 | 2,7 | 18,8 |
| Ridge (resíduo) | baseline | 0,0429 | [0,0335; 0,0523] | – | 0,9517 | 0,721 | 76 | 0,9 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,0811 | [0,0579; 0,1052] | – | 0,8276 | 0,004 | 0 | 24,1 | 32,1 |
| White-Box | white-box | 0,0813 | [0,0552; 0,1047] | – | 0,8269 | 0,000 | 0 | 0,0 | 0,0 |

**Deixa-um-η-fora (extrapolação em η)** (teto do ganho = 0,953)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ridge (resíduo) | baseline | 0,0579 | [0,0450; 0,0708] | – | 0,9120 | 0,492 | 52 | 8,9 | 0,0 |
| Híbrido MLP tanh [B] | hibrido | 0,0615 | [0,0507; 0,0722] | – | 0,9007 | 0,427 | 45 | 23,2 | 17,0 |
| Tabela (η, t) | baseline | 0,0637 | [0,0453; 0,0805] | – | 0,8938 | 0,386 | 41 | 0,0 | 17,0 |
| Híbrido Extra Trees [B] | hibrido | 0,0644 | [0,0493; 0,0781] | – | 0,8913 | 0,372 | 39 | 0,0 | 6,2 |
| Híbrido GP [B] | hibrido | 0,0672 | [0,0540; 0,0798] | – | 0,8817 | 0,317 | 33 | 0,0 | 5,4 |
| Híbrido GBDT [B] | hibrido | 0,0698 | [0,0535; 0,0840] | – | 0,8723 | 0,262 | 28 | 0,0 | 5,4 |
| Híbrido Random Forest [B] | hibrido | 0,0711 | [0,0530; 0,0868] | – | 0,8676 | 0,235 | 25 | 0,0 | 3,6 |
| Híbrido SVR [B] | hibrido | 0,0737 | [0,0607; 0,0853] | – | 0,8574 | 0,177 | 19 | 1,8 | 26,8 |
| Caixa-preta Extra Trees | caixa-preta | 0,0764 | [0,0612; 0,0907] | – | 0,8470 | 0,116 | 12 | 25,0 | 8,9 |
| Caixa-preta GBDT | caixa-preta | 0,0794 | [0,0646; 0,0917] | – | 0,8348 | 0,046 | 5 | 25,0 | 1,8 |
| White-Box | white-box | 0,0813 | [0,0552; 0,1047] | – | 0,8269 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta GP | caixa-preta | 0,0978 | [0,0763; 0,1177] | – | 0,7491 | -0,449 | -47 | 47,3 | 17,0 |
| Caixa-preta SVR | caixa-preta | 0,4919 | [0,2529; 0,6835] | – | -5,3409 | -35,629 | -3737 | 41,1 | 25,0 |

**Deixa-um-C_A0-fora** (teto do ganho = 0,953)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido Extra Trees [B] | hibrido | 0,0210 | [0,0150; 0,0264] | – | 0,9885 | 0,933 | 98 | 2,7 | 2,7 |
| Híbrido GBDT [B] | hibrido | 0,0232 | [0,0177; 0,0287] | – | 0,9859 | 0,918 | 96 | 0,9 | 3,6 |
| Híbrido MLP tanh [B] | hibrido | 0,0237 | [0,0178; 0,0294] | – | 0,9853 | 0,915 | 96 | 4,5 | 4,5 |
| Híbrido Random Forest [B] | hibrido | 0,0237 | [0,0155; 0,0313] | – | 0,9853 | 0,915 | 96 | 0,0 | 5,4 |
| Híbrido GP [B] | hibrido | 0,0239 | [0,0180; 0,0300] | – | 0,9851 | 0,914 | 96 | 5,4 | 4,5 |
| Híbrido SVR [B] | hibrido | 0,0250 | [0,0182; 0,0313] | – | 0,9836 | 0,905 | 95 | 4,5 | 10,7 |
| Caixa-preta Extra Trees | caixa-preta | 0,0276 | [0,0235; 0,0318] | – | 0,9801 | 0,885 | 93 | 2,7 | 10,7 |
| Caixa-preta GP | caixa-preta | 0,0279 | [0,0217; 0,0341] | – | 0,9795 | 0,882 | 92 | 6,2 | 25,0 |
| Caixa-preta GBDT | caixa-preta | 0,0290 | [0,0220; 0,0353] | – | 0,9780 | 0,873 | 92 | 0,0 | 3,6 |
| Tabela (η, t) | baseline | 0,0296 | [0,0219; 0,0381] | – | 0,9771 | 0,868 | 91 | 2,7 | 18,8 |
| Ridge (resíduo) | baseline | 0,0383 | [0,0297; 0,0465] | – | 0,9615 | 0,777 | 82 | 0,0 | 0,0 |
| White-Box | white-box | 0,0813 | [0,0552; 0,1047] | – | 0,8269 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,1815 | [0,1144; 0,2370] | – | 0,1368 | -3,986 | -418 | 16,1 | 31,2 |

**Extrapolação temporal (t<=5 -> 15 min)** (teto do ganho = 0,943)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido Extra Trees [B] | hibrido | 0,0114 | [0,0086; 0,0142] | – | 0,9968 | 0,976 | 103 | 0,0 | 0,0 |
| Híbrido GBDT [B] | hibrido | 0,0114 | [0,0078; 0,0148] | – | 0,9968 | 0,976 | 103 | 0,0 | 0,0 |
| Híbrido Random Forest [B] | hibrido | 0,0116 | [0,0085; 0,0144] | – | 0,9967 | 0,975 | 103 | 0,0 | 0,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0121 | [0,0095; 0,0147] | – | 0,9964 | 0,973 | 103 | 0,0 | 0,0 |
| Híbrido GP [B] | hibrido | 0,0134 | [0,0107; 0,0157] | – | 0,9956 | 0,967 | 103 | 0,0 | 0,0 |
| Híbrido MLP tanh [B] | hibrido | 0,0165 | [0,0090; 0,0233] | – | 0,9932 | 0,950 | 101 | 0,0 | 0,0 |
| Caixa-preta GBDT | caixa-preta | 0,0213 | [0,0139; 0,0287] | – | 0,9888 | 0,916 | 97 | 0,0 | 0,0 |
| Ridge (resíduo) | baseline | 0,0648 | [0,0379; 0,0869] | – | 0,8962 | 0,226 | 24 | 50,0 | 0,0 |
| Híbrido SVR [B] | hibrido | 0,0655 | [0,0488; 0,0822] | – | 0,8939 | 0,210 | 22 | 12,5 | 0,0 |
| White-Box | white-box | 0,0736 | [0,0417; 0,0996] | – | 0,8658 | 0,000 | 0 | 0,0 | 0,0 |
| Tabela (η, t) | baseline | 0,0736 | [0,0417; 0,0996] | – | 0,8658 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta GP | caixa-preta | 0,2000 | [0,1543; 0,2391] | – | 0,0100 | -6,377 | -676 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,2892 | [0,2295; 0,3430] | – | -1,0700 | -14,425 | -1529 | 12,5 | 0,0 |

### T06. Extrapolação em η: RMSE no nível deixado fora do treino

| Modelo | RMSE total | η = 0,5 ausente | η = 1,0 ausente | η = 1,5 ausente | η = 3,1 ausente |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Ridge (resíduo) | 0,0579 | 0,0461 | 0,0872 | 0,0486 | 0,0365 |
| Híbrido MLP tanh [B] | 0,0615 | 0,0537 | 0,0818 | 0,0630 | 0,0399 |
| Tabela (η, t) | 0,0637 | 0,0289 | 0,1040 | 0,0310 | 0,0599 |
| Híbrido Extra Trees [B] | 0,0644 | 0,0287 | 0,0942 | 0,0746 | 0,0363 |
| Híbrido GP [B] | 0,0672 | 0,0248 | 0,0940 | 0,0767 | 0,0521 |
| Híbrido GBDT [B] | 0,0698 | 0,0291 | 0,0916 | 0,0950 | 0,0349 |
| Híbrido Random Forest [B] | 0,0711 | 0,0229 | 0,0932 | 0,0992 | 0,0339 |
| Híbrido SVR [B] | 0,0737 | 0,0753 | 0,0938 | 0,0733 | 0,0436 |
| Caixa-preta Extra Trees | 0,0764 | 0,0525 | 0,1129 | 0,0782 | 0,0415 |
| Caixa-preta GBDT | 0,0794 | 0,0525 | 0,0953 | 0,1078 | 0,0419 |
| White-Box | 0,0813 | 0,0516 | 0,1396 | 0,0544 | 0,0361 |
| Caixa-preta GP | 0,0978 | 0,0525 | 0,1034 | 0,0610 | 0,1453 |
| Caixa-preta SVR | 0,4919 | 0,0464 | 0,1207 | 0,0696 | 0,9727 |

### T07. Híbrido contra a caixa-preta do MESMO modelo (diferença pareada, bootstrap de ensaios)

| Protocolo | Híbrido | Caixa-preta | Δ RMSE (híbrido − caixa-preta) | IC95% | Empate (IC contém 0) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0051 | [-0,0095; -0,0005] | não |
| GKF8 (partição do TCC) | Híbrido GP [B] | Caixa-preta GP | -0,0020 | [-0,0069; 0,0029] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0028 | [-0,0072; 0,0016] | sim |
| GKF8 (partição do TCC) | Híbrido SVR [B] | Caixa-preta SVR | -0,0655 | [-0,0878; -0,0457] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0021 | [-0,0062; 0,0017] | sim |
| GKF8 repetido (20 partições) | Híbrido GP [B] | Caixa-preta GP | -0,0016 | [-0,0060; 0,0027] | sim |
| GKF8 repetido (20 partições) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0043 | [-0,0082; -0,0004] | não |
| GKF8 repetido (20 partições) | Híbrido SVR [B] | Caixa-preta SVR | -0,0826 | [-0,1205; -0,0498] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0024 | [-0,0068; 0,0016] | sim |
| LOGO (16 ensaios) | Híbrido GP [B] | Caixa-preta GP | -0,0026 | [-0,0073; 0,0019] | sim |
| LOGO (16 ensaios) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0047 | [-0,0086; -0,0006] | não |
| LOGO (16 ensaios) | Híbrido SVR [B] | Caixa-preta SVR | -0,0570 | [-0,0817; -0,0341] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido GP [B] | Caixa-preta GP | -0,0307 | [-0,0546; -0,0064] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0120 | [-0,0206; -0,0053] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0096 | [-0,0156; -0,0043] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido SVR [B] | Caixa-preta SVR | -0,4181 | [-0,6134; -0,1720] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0066 | [-0,0133; -0,0009] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GP | -0,0041 | [-0,0108; 0,0020] | sim |
| Deixa-um-C_A0-fora | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0057 | [-0,0112; -0,0007] | não |
| Deixa-um-C_A0-fora | Híbrido SVR [B] | Caixa-preta SVR | -0,1565 | [-0,2132; -0,0897] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0008 | [-0,0026; 0,0007] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GP [B] | Caixa-preta GP | -0,1867 | [-0,2251; -0,1412] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0099 | [-0,0153; -0,0040] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido SVR [B] | Caixa-preta SVR | -0,2237 | [-0,2741; -0,1712] | não |

### T08. Melhor híbrido de cada protocolo contra os demais modelos (Δ < 0: o melhor híbrido erra menos)

| Protocolo | Melhor híbrido | Contra | Δ RMSE | IC95% | Empate |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | White-Box | -0,0598 | [-0,0860; -0,0331] | não |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Tabela (η, t) | -0,0081 | [-0,0145; -0,0027] | não |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Ridge (resíduo) | -0,0209 | [-0,0321; -0,0112] | não |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Híbrido GP [B] | -0,0006 | [-0,0035; 0,0022] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Híbrido GBDT [B] | -0,0001 | [-0,0027; 0,0019] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Híbrido Random Forest [B] | -0,0012 | [-0,0044; 0,0021] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Híbrido SVR [B] | -0,0003 | [-0,0048; 0,0038] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Híbrido MLP tanh [B] | -0,0016 | [-0,0062; 0,0026] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta GP | -0,0026 | [-0,0064; 0,0012] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0051 | [-0,0095; -0,0005] | não |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta GBDT | -0,0028 | [-0,0067; 0,0008] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta SVR | -0,0658 | [-0,0892; -0,0451] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | White-Box | -0,0592 | [-0,0833; -0,0334] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Tabela (η, t) | -0,0078 | [-0,0144; -0,0028] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Ridge (resíduo) | -0,0210 | [-0,0308; -0,0119] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido GP [B] | -0,0003 | [-0,0029; 0,0024] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido Extra Trees [B] | -0,0008 | [-0,0028; 0,0014] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido Random Forest [B] | -0,0009 | [-0,0042; 0,0022] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido SVR [B] | -0,0016 | [-0,0044; 0,0011] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido MLP tanh [B] | -0,0014 | [-0,0049; 0,0022] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GP | -0,0019 | [-0,0055; 0,0014] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta Extra Trees | -0,0051 | [-0,0093; -0,0009] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0021 | [-0,0062; 0,0017] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta SVR | -0,0842 | [-0,1224; -0,0510] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | White-Box | -0,0598 | [-0,0839; -0,0340] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Tabela (η, t) | -0,0081 | [-0,0146; -0,0030] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Ridge (resíduo) | -0,0215 | [-0,0313; -0,0124] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido GP [B] | -0,0001 | [-0,0030; 0,0030] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido Extra Trees [B] | -0,0003 | [-0,0022; 0,0019] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido Random Forest [B] | -0,0013 | [-0,0046; 0,0020] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido SVR [B] | -0,0027 | [-0,0062; 0,0010] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido MLP tanh [B] | -0,0034 | [-0,0083; 0,0012] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GP | -0,0027 | [-0,0066; 0,0010] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta Extra Trees | -0,0051 | [-0,0096; -0,0006] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0024 | [-0,0068; 0,0016] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta SVR | -0,0596 | [-0,0850; -0,0354] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | White-Box | -0,0197 | [-0,0354; -0,0007] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Tabela (η, t) | -0,0021 | [-0,0142; 0,0129] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Ridge (resíduo) | 0,0036 | [-0,0011; 0,0096] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Híbrido GP [B] | -0,0056 | [-0,0131; 0,0035] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Híbrido Extra Trees [B] | -0,0028 | [-0,0101; 0,0069] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Híbrido GBDT [B] | -0,0083 | [-0,0180; 0,0031] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Híbrido Random Forest [B] | -0,0095 | [-0,0203; 0,0033] | sim |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Híbrido SVR [B] | -0,0122 | [-0,0173; -0,0076] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Caixa-preta GP | -0,0363 | [-0,0635; -0,0102] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Caixa-preta Extra Trees | -0,0148 | [-0,0264; -0,0033] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Caixa-preta GBDT | -0,0178 | [-0,0297; -0,0057] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido MLP tanh [B] | Caixa-preta SVR | -0,4303 | [-0,6237; -0,1863] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | White-Box | -0,0603 | [-0,0859; -0,0336] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Tabela (η, t) | -0,0086 | [-0,0137; -0,0043] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Ridge (resíduo) | -0,0174 | [-0,0276; -0,0085] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Híbrido GP [B] | -0,0029 | [-0,0054; -0,0003] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Híbrido GBDT [B] | -0,0023 | [-0,0048; -0,0001] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Híbrido Random Forest [B] | -0,0027 | [-0,0054; 0,0002] | sim |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Híbrido SVR [B] | -0,0040 | [-0,0081; 0,0001] | sim |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Híbrido MLP tanh [B] | -0,0027 | [-0,0066; 0,0005] | sim |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta GP | -0,0070 | [-0,0145; -0,0009] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0066 | [-0,0133; -0,0009] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta GBDT | -0,0080 | [-0,0152; -0,0019] | não |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta SVR | -0,1605 | [-0,2181; -0,0931] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | White-Box | -0,0623 | [-0,0889; -0,0299] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Tabela (η, t) | -0,0623 | [-0,0889; -0,0299] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Ridge (resíduo) | -0,0534 | [-0,0764; -0,0260] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Híbrido GP [B] | -0,0020 | [-0,0046; 0,0002] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Híbrido GBDT [B] | -0,0000 | [-0,0017; 0,0018] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Híbrido Random Forest [B] | -0,0002 | [-0,0019; 0,0014] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Híbrido SVR [B] | -0,0541 | [-0,0689; -0,0388] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Híbrido MLP tanh [B] | -0,0052 | [-0,0130; 0,0031] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta GP | -0,1887 | [-0,2271; -0,1433] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0008 | [-0,0026; 0,0007] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta GBDT | -0,0100 | [-0,0148; -0,0045] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta SVR | -0,2779 | [-0,3313; -0,2182] | não |

### T09. Curva de aprendizado: RMSE em 2 ensaios fora do treino, por nº de ensaios de treino (30 sorteios)

| Modelo | n = 4 | n = 6 | n = 8 | n = 10 | n = 12 | n = 14 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| White-Box | 0,0811 | 0,0811 | 0,0811 | 0,0811 | 0,0811 | 0,0811 |
| Tabela (η, t) | 0,0509 | 0,0374 | 0,0371 | 0,0323 | 0,0304 | 0,0296 |
| Híbrido Extra Trees [B] | 0,0530 | 0,0390 | 0,0319 | 0,0255 | 0,0219 | 0,0216 |
| Híbrido GBDT [B] | 0,0573 | 0,0381 | 0,0315 | 0,0242 | 0,0214 | 0,0214 |
| Híbrido GP [B] | 0,0521 | 0,0404 | 0,0376 | 0,0268 | 0,0225 | 0,0217 |
| Caixa-preta Extra Trees | 0,0897 | 0,0475 | 0,0388 | 0,0295 | 0,0295 | 0,0269 |
| Caixa-preta GBDT | 0,1007 | 0,0433 | 0,0332 | 0,0271 | 0,0258 | 0,0242 |
| Caixa-preta GP | 0,0880 | 0,0537 | 0,0386 | 0,0277 | 0,0256 | 0,0237 |

### T10. Ablação dos atributos do ML (RMSE em t>0) (conjunto C omitido: exige [H⁺], indisponível com a força motriz C_Af)

| Protocolo | Modelo | A (t, η, C_A0) | B (A + X_PBM) |
| :--- | ---: | ---: | ---: |
| GKF8 repetido (10) | GP | 0,0255 | 0,0223 |
| GKF8 repetido (10) | Extra Trees | 0,0265 | 0,0229 |
| GKF8 repetido (10) | GBDT | 0,0312 | 0,0221 |
| GKF8 repetido (10) | SVR | 0,0529 | 0,0242 |
| Deixa-um-η-fora | GP | 0,0739 | 0,0672 |
| Deixa-um-η-fora | Extra Trees | 0,0752 | 0,0644 |
| Deixa-um-η-fora | GBDT | 0,0830 | 0,0698 |
| Deixa-um-η-fora | SVR | 0,1280 | 0,0737 |

### T11. Sensibilidade ao fator ν (RMSE em t>0)

| Protocolo | Modelo | nu estimado na dobra (comparativo) | nu a priori (1.06, comparativo) | nu = 1 (White-Box principal) |
| :--- | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | White-Box | 0,0722 | 0,0739 | 0,0813 |
| GKF8 (partição do TCC) | Híbrido GP [B] | 0,0207 | 0,0209 | 0,0221 |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | 0,0231 | 0,0235 | 0,0215 |
| Deixa-um-η-fora | White-Box | 0,0720 | 0,0739 | 0,0813 |
| Deixa-um-η-fora | Híbrido GP [B] | 0,0648 | 0,0655 | 0,0672 |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | 0,0593 | 0,0624 | 0,0644 |

### T12. O mesmo ML sobre degraus do White-Box (RMSE em t>0)

| Protocolo | Modelo | 1. Tese (α = 0) | 2. White-Box do TCC (V0, nu = 1) | 3. V0 + nu (estimado na dobra, comparativo) | 4. CAf, nu = 1 (principal) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 repetido (10) | White-Box | 0,0981 | 0,0443 | 0,0393 | 0,0813 |
| GKF8 repetido (10) | Híbrido GP [B] | 0,0437 | 0,0235 | 0,0237 | 0,0223 |
| GKF8 repetido (10) | Híbrido Extra Trees [B] | 0,0346 | 0,0200 | 0,0196 | 0,0229 |
| GKF8 repetido (10) | Caixa-preta GP | 0,0237 | 0,0237 | 0,0231 | 0,0237 |
| GKF8 repetido (10) | Caixa-preta Extra Trees | 0,0270 | 0,0270 | 0,0269 | 0,0270 |
| Deixa-um-η-fora | White-Box | 0,0981 | 0,0443 | 0,0394 | 0,0813 |
| Deixa-um-η-fora | Híbrido GP [B] | 0,0837 | 0,0455 | 0,0392 | 0,0672 |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | 0,0661 | 0,0341 | 0,0290 | 0,0644 |
| Deixa-um-η-fora | Caixa-preta GP | 0,0978 | 0,0978 | 0,0950 | 0,0978 |
| Deixa-um-η-fora | Caixa-preta Extra Trees | 0,0764 | 0,0764 | 0,0727 | 0,0764 |
| Poucos dados (n = 4 ensaios de treino) | White-Box | 0,0904 | 0,0450 | 0,0403 | 0,0796 |
| Poucos dados (n = 4 ensaios de treino) | Híbrido GP [B] | 0,0668 | 0,0335 | 0,0315 | 0,0501 |
| Poucos dados (n = 4 ensaios de treino) | Híbrido Extra Trees [B] | 0,0639 | 0,0351 | 0,0296 | 0,0524 |
| Poucos dados (n = 4 ensaios de treino) | Caixa-preta GP | 0,0650 | 0,0650 | 0,0606 | 0,0650 |
| Poucos dados (n = 4 ensaios de treino) | Caixa-preta Extra Trees | 0,0873 | 0,0873 | 0,0822 | 0,0873 |
| Poucos dados (n = 8 ensaios de treino) | White-Box | 0,0904 | 0,0450 | 0,0395 | 0,0796 |
| Poucos dados (n = 8 ensaios de treino) | Híbrido GP [B] | 0,0555 | 0,0279 | 0,0261 | 0,0321 |
| Poucos dados (n = 8 ensaios de treino) | Híbrido Extra Trees [B] | 0,0478 | 0,0225 | 0,0207 | 0,0273 |
| Poucos dados (n = 8 ensaios de treino) | Caixa-preta GP | 0,0276 | 0,0276 | 0,0273 | 0,0276 |
| Poucos dados (n = 8 ensaios de treino) | Caixa-preta Extra Trees | 0,0344 | 0,0344 | 0,0310 | 0,0344 |

### T13. Cobertura empírica dos intervalos de 95% (t>0, pontos fora-da-dobra)

| Protocolo | Modelo | Método | Cobertura | Largura média | n |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | Conforme (LOGO interno) | 43% | 0,0797 | 112 |
| Deixa-um-η-fora | Híbrido GBDT [B] | Conforme (LOGO interno) | 46% | 0,0815 | 112 |
| Deixa-um-η-fora | Híbrido GP [B] | Conforme (LOGO interno) | 36% | 0,0953 | 112 |
| Deixa-um-η-fora | Híbrido GP [B] | GP (analítico) | 60% | 0,1409 | 112 |
| LOGO (16 ensaios) | Híbrido Extra Trees [B] | Conforme (LOGO interno) | 94% | 0,0844 | 112 |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Conforme (LOGO interno) | 96% | 0,0805 | 112 |
| LOGO (16 ensaios) | Híbrido GP [B] | Conforme (LOGO interno) | 96% | 0,0852 | 112 |
| LOGO (16 ensaios) | Híbrido GP [B] | GP (analítico) | 90% | 0,0716 | 112 |

### T14. Resíduo do White-Box principal aos 15 min

| Ensaio | η | η real | C_A0 | X_exp (15 min) | X_PBM (15 min) | ΔX |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,535 | 0,1 | 0,50 | 0,5347 | -0,035 |
| 2 | 3,1 | 3,290 | 0,1 | 1,00 | 1,0000 | 0,000 |
| 3 | 0,5 | 0,535 | 0,5 | 0,50 | 0,5347 | -0,035 |
| 4 | 0,5 | 0,535 | 1,0 | 0,50 | 0,5347 | -0,035 |
| 5 | 0,5 | 0,535 | 1,5 | 0,49 | 0,5347 | -0,045 |
| 6 | 1,0 | 1,069 | 0,1 | 0,87 | 0,9967 | -0,127 |
| 7 | 3,1 | 3,192 | 0,5 | 1,01 | 1,0000 | 0,010 |
| 8 | 1,0 | 1,069 | 0,5 | 0,87 | 1,0000 | -0,130 |
| 9 | 1,0 | 1,069 | 1,0 | 0,85 | 1,0000 | -0,150 |
| 10 | 1,0 | 1,069 | 1,5 | 0,85 | 1,0000 | -0,150 |
| 11 | 1,5 | 1,584 | 0,1 | 0,98 | 1,0000 | -0,020 |
| 12 | 3,1 | 3,216 | 1,0 | 1,00 | 1,0000 | 0,000 |
| 13 | 1,5 | 1,608 | 0,5 | 0,97 | 1,0000 | -0,030 |
| 14 | 1,5 | 1,602 | 1,0 | 0,97 | 1,0000 | -0,030 |
| 15 | 1,5 | 1,604 | 1,5 | 0,97 | 1,0000 | -0,030 |
| 16 | 3,1 | 3,208 | 1,5 | 1,00 | 1,0000 | 0,000 |

### T15. Escalas de comprimento do GP (ARD) ajustadas com todos os dados

| Atributos | Entrada | Escala de comprimento (σ das entradas) | Ruído branco (norm.) |
| :--- | ---: | ---: | ---: |
| A | tempo_min | 0,706 | 0,063 |
| A | razao_molar_eta_real | 0,065 | 0,063 |
| A | C_acid_0_mol_L | 3,655 | 0,063 |
| B | tempo_min | 4,267 | 0,069 |
| B | razao_molar_eta_real | 0,549 | 0,069 |
| B | C_acid_0_mol_L | 6,235 | 0,069 |
| B | X_pbm | 0,288 | 0,069 |

### T16. Conversão do White-Box principal e resíduo, por ensaio

**X_PBM (conversão do White-Box principal)**

| Ensaio | η | C_A0 | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 0,00 | 0,52 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 |
| 2 | 3,1 | 0,1 | 0,00 | 0,89 | 0,97 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 3 | 0,5 | 0,5 | 0,00 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 |
| 4 | 0,5 | 1,0 | 0,00 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 |
| 5 | 0,5 | 1,5 | 0,00 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 | 0,53 |
| 6 | 1,0 | 0,1 | 0,00 | 0,77 | 0,86 | 0,92 | 0,95 | 0,96 | 0,97 | 1,00 |
| 7 | 3,1 | 0,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 8 | 1,0 | 0,5 | 0,00 | 0,94 | 0,97 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 |
| 9 | 1,0 | 1,0 | 0,00 | 0,97 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 10 | 1,0 | 1,5 | 0,00 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 11 | 1,5 | 0,1 | 0,00 | 0,83 | 0,93 | 0,98 | 0,99 | 1,00 | 1,00 | 1,00 |
| 12 | 3,1 | 1,0 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 13 | 1,5 | 0,5 | 0,00 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 14 | 1,5 | 1,0 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 15 | 1,5 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 16 | 3,1 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |

**ΔX = X_exp − X_PBM**

| Ensaio | η | C_A0 | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 0,00 | -0,10 | -0,07 | -0,07 | -0,04 | -0,02 | -0,04 | -0,03 |
| 2 | 3,1 | 0,1 | 0,00 | -0,07 | -0,07 | -0,05 | -0,03 | -0,03 | -0,01 | 0,00 |
| 3 | 0,5 | 0,5 | 0,00 | -0,08 | -0,05 | -0,04 | -0,01 | -0,03 | -0,04 | -0,03 |
| 4 | 0,5 | 1,0 | 0,00 | -0,03 | -0,06 | -0,02 | -0,06 | -0,00 | -0,02 | -0,03 |
| 5 | 0,5 | 1,5 | 0,00 | -0,04 | -0,04 | -0,06 | -0,04 | -0,04 | -0,06 | -0,04 |
| 6 | 1,0 | 0,1 | 0,00 | -0,02 | -0,08 | -0,12 | -0,10 | -0,09 | -0,12 | -0,13 |
| 7 | 3,1 | 0,5 | 0,00 | -0,12 | -0,07 | -0,02 | -0,01 | 0,00 | -0,01 | 0,01 |
| 8 | 1,0 | 0,5 | 0,00 | -0,11 | -0,15 | -0,14 | -0,14 | -0,13 | -0,14 | -0,13 |
| 9 | 1,0 | 1,0 | 0,00 | -0,14 | -0,15 | -0,15 | -0,15 | -0,17 | -0,16 | -0,15 |
| 10 | 1,0 | 1,5 | 0,00 | -0,15 | -0,17 | -0,17 | -0,17 | -0,17 | -0,16 | -0,15 |
| 11 | 1,5 | 0,1 | 0,00 | 0,03 | -0,08 | -0,10 | -0,02 | -0,06 | -0,04 | -0,02 |
| 12 | 3,1 | 1,0 | 0,00 | -0,04 | 0,00 | 0,00 | -0,01 | -0,01 | -0,01 | 0,00 |
| 13 | 1,5 | 0,5 | 0,00 | -0,07 | -0,08 | -0,05 | -0,05 | -0,05 | -0,05 | -0,03 |
| 14 | 1,5 | 1,0 | 0,00 | -0,09 | -0,04 | -0,05 | -0,03 | -0,05 | -0,02 | -0,03 |
| 15 | 1,5 | 1,5 | 0,00 | -0,09 | -0,06 | -0,05 | -0,06 | -0,04 | -0,04 | -0,03 |
| 16 | 3,1 | 1,5 | 0,00 | -0,02 | -0,02 | -0,04 | -0,01 | -0,01 | 0,00 | 0,00 |

### T17. Partição comum do relatório (teste cego 7, 8, 14; CV nas 4 dobras da Tabela B.3; 10 sementes)

Modelo selecionado só pela CV: **Híbrido GBDT [B]**. Erros em p.p. de conversão, t > 0; inconsistência = % dos nós da malha de 0,25 min com conversão decrescente ou acima do teto, antes (bruta) e depois do acoplador.

| Conjunto | Ensaios | Pontos | R² | RMSE (p.p.) | MAE (p.p.) | Erro máximo (p.p.) | Inconsistência bruta (%) | Inconsistência final (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Treino | 13 | 91 | 0,997 | 1,12 | 0,85 | 3,51 | 2,3 | 0,0 |
| Teste cego | 3 | 21 | 0,878 | 2,11 | 1,40 | 7,27 | 2,2 | 0,0 |
| Global | 16 | 112 | 0,995 | 1,36 | 0,95 | 7,27 | 2,3 | 0,0 |

Todos os modelos, ordenados pelo erro da CV (o teste cego não participa da escolha):

| Modelo | Tipo | CV fora da dobra (p.p.) | Teste RMSE (p.p.) | Teste R² | E7 (p.p.) | E8 (p.p.) | E14 (p.p.) | Inconsist. bruta (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GBDT [B] | hibrido | 1,97 ± 0,00 | 2,11 | 0,878 | 2,99 | 1,32 | 1,64 | 2,2 |
| Híbrido GP [B] | hibrido | 1,98 ± 0,00 | 2,72 | 0,796 | 3,99 | 1,80 | 1,77 | 1,1 |
| Híbrido HistGBDT [B] | hibrido | 2,06 ± 0,00 | 2,08 | 0,881 | 2,84 | 1,73 | 1,38 | 5,0 |
| Híbrido Random Forest [B] | hibrido | 2,08 ± 0,03 | 2,36 | 0,847 | 3,31 | 1,89 | 1,49 | 1,7 |
| Caixa-preta MLP tanh | caixa-preta | 2,23 ± 0,03 | 1,90 | 0,901 | 1,30 | 2,83 | 1,04 | 13,9 |
| Híbrido Extra Trees [B] | hibrido | 2,32 ± 0,06 | 2,29 | 0,855 | 3,29 | 1,14 | 1,91 | 0,0 |
| Híbrido AdaBoost [B] | hibrido | 2,34 ± 0,06 | 2,45 | 0,835 | 3,47 | 1,91 | 1,50 | 0,0 |
| Caixa-preta GBDT | caixa-preta | 2,57 ± 0,00 | 2,03 | 0,886 | 2,30 | 2,05 | 1,71 | 0,0 |
| Caixa-preta GP | caixa-preta | 2,62 ± 0,00 | 2,32 | 0,852 | 2,50 | 2,44 | 1,98 | 34,4 |
| Caixa-preta HistGBDT | caixa-preta | 2,65 ± 0,00 | 1,93 | 0,897 | 2,32 | 1,86 | 1,54 | 14,4 |
| Caixa-preta GP (RBF) | caixa-preta | 2,80 ± 0,05 | 2,80 | 0,785 | 3,64 | 2,44 | 2,07 | 27,8 |
| Caixa-preta Random Forest | caixa-preta | 2,81 ± 0,05 | 2,32 | 0,852 | 2,49 | 2,66 | 1,68 | 0,0 |
| Híbrido MLP tanh [B] | hibrido | 2,95 ± 0,23 | 2,27 | 0,859 | 3,21 | 0,92 | 2,06 | 0,0 |
| Caixa-preta Extra Trees | caixa-preta | 3,00 ± 0,11 | 2,74 | 0,793 | 2,10 | 3,24 | 2,78 | 0,0 |
| Tabela (η, t) | baseline | 3,19 ± 0,00 | 2,31 | 0,853 | 3,39 | 1,18 | 1,79 | 11,1 |
| Híbrido Kernel Ridge [B] | hibrido | 3,21 ± 0,00 | 1,79 | 0,912 | 2,38 | 1,30 | 1,52 | 0,0 |
| Caixa-preta AdaBoost | caixa-preta | 3,35 ± 0,04 | 3,36 | 0,689 | 4,74 | 2,54 | 2,24 | 0,0 |
| Híbrido GP (RBF) [B] | hibrido | 3,47 ± 1,30 | 2,99 | 0,755 | 3,89 | 2,56 | 2,25 | 6,7 |
| Híbrido k-NN [B] | hibrido | 3,67 ± 0,00 | 5,59 | 0,140 | 3,59 | 5,67 | 6,98 | 4,4 |
| Ridge (resíduo) | baseline | 4,02 ± 0,00 | 4,17 | 0,521 | 3,82 | 4,69 | 3,96 | 0,0 |
| Híbrido SVR [B] | hibrido | 4,32 ± 0,00 | 1,78 | 0,913 | 2,22 | 1,55 | 1,46 | 2,2 |
| Híbrido Ridge polinomial [B] | hibrido | 4,65 ± 0,00 | 3,95 | 0,571 | 3,49 | 4,21 | 4,11 | 5,6 |
| Caixa-preta Ridge polinomial | caixa-preta | 7,16 ± 0,00 | 7,12 | -0,394 | 7,13 | 6,88 | 7,34 | 45,0 |
| White-Box | white-box | 7,96 ± 0,00 | 8,81 | -1,134 | 5,28 | 13,45 | 4,89 | 0,0 |
| Caixa-preta Kernel Ridge | caixa-preta | 9,44 ± 0,00 | 5,45 | 0,183 | 5,18 | 6,29 | 4,76 | 45,6 |
| Caixa-preta k-NN | caixa-preta | 18,59 ± 0,00 | 10,94 | -2,292 | 5,59 | 13,09 | 12,50 | 6,1 |
| Caixa-preta SVR | caixa-preta | 19,85 ± 0,00 | 8,71 | -1,088 | 7,73 | 9,90 | 8,37 | 42,8 |
