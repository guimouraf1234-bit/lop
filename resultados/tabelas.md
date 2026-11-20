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
| Laboratório: CAf, nu = 1 (principal) | 0,0702 | 0,8710 | 0,9590 | – | – | – |

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

**GKF8 (partição do TCC)** (teto do ganho = 0,937)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Caixa-preta GBDT | caixa-preta | 0,0218 | [0,0182; 0,0252] | – | 0,9875 | 0,903 | 96 | 8,0 | 6,2 |
| Híbrido GBDT [B] | hibrido | 0,0219 | [0,0174; 0,0264] | – | 0,9875 | 0,903 | 96 | 2,7 | 9,8 |
| Caixa-preta GP | caixa-preta | 0,0221 | [0,0186; 0,0259] | – | 0,9871 | 0,900 | 96 | 17,0 | 21,4 |
| Híbrido Random Forest [B] | hibrido | 0,0227 | [0,0178; 0,0281] | – | 0,9865 | 0,895 | 96 | 6,2 | 14,3 |
| Híbrido GP [B] | hibrido | 0,0231 | [0,0175; 0,0291] | – | 0,9860 | 0,891 | 95 | 11,6 | 4,5 |
| Híbrido SVR [B] | hibrido | 0,0237 | [0,0190; 0,0281] | – | 0,9853 | 0,886 | 95 | 13,4 | 17,0 |
| Híbrido MLP tanh [B] | hibrido | 0,0244 | [0,0195; 0,0294] | – | 0,9844 | 0,879 | 94 | 10,7 | 4,5 |
| Híbrido Extra Trees [B] | hibrido | 0,0247 | [0,0176; 0,0313] | – | 0,9840 | 0,876 | 93 | 3,6 | 5,4 |
| Caixa-preta Extra Trees | caixa-preta | 0,0256 | [0,0227; 0,0282] | – | 0,9829 | 0,867 | 93 | 6,2 | 8,0 |
| Tabela (η, t) | baseline | 0,0338 | [0,0231; 0,0452] | – | 0,9701 | 0,768 | 82 | 5,4 | 19,6 |
| Ridge (resíduo) | baseline | 0,0447 | [0,0347; 0,0546] | – | 0,9475 | 0,593 | 63 | 0,0 | 0,0 |
| White-Box | white-box | 0,0702 | [0,0441; 0,0930] | – | 0,8710 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,0853 | [0,0662; 0,1034] | – | 0,8093 | -0,478 | -51 | 36,6 | 37,5 |

**GKF8 repetido (20 partições)** (teto do ganho = 0,937)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GBDT [B] | hibrido | 0,0217 | [0,0173; 0,0263] | 0,0005 | 0,9877 | 0,904 | 96 | 2,4 | 9,3 |
| Caixa-preta GP | caixa-preta | 0,0221 | [0,0185; 0,0259] | 0,0010 | 0,9871 | 0,900 | 96 | 15,6 | 24,2 |
| Caixa-preta GBDT | caixa-preta | 0,0222 | [0,0186; 0,0257] | 0,0008 | 0,9871 | 0,900 | 96 | 11,5 | 6,6 |
| Híbrido Random Forest [B] | hibrido | 0,0227 | [0,0179; 0,0281] | 0,0007 | 0,9864 | 0,895 | 95 | 6,2 | 13,6 |
| Híbrido GP [B] | hibrido | 0,0229 | [0,0178; 0,0283] | 0,0011 | 0,9862 | 0,893 | 95 | 9,7 | 3,3 |
| Híbrido MLP tanh [B] | hibrido | 0,0239 | [0,0193; 0,0285] | 0,0015 | 0,9850 | 0,884 | 94 | 10,0 | 8,1 |
| Híbrido SVR [B] | hibrido | 0,0241 | [0,0209; 0,0274] | 0,0023 | 0,9846 | 0,881 | 94 | 13,7 | 13,9 |
| Híbrido Extra Trees [B] | hibrido | 0,0254 | [0,0185; 0,0319] | 0,0009 | 0,9830 | 0,869 | 93 | 3,7 | 5,8 |
| Caixa-preta Extra Trees | caixa-preta | 0,0262 | [0,0230; 0,0290] | 0,0011 | 0,9820 | 0,861 | 92 | 6,5 | 10,5 |
| Tabela (η, t) | baseline | 0,0342 | [0,0234; 0,0456] | 0,0010 | 0,9693 | 0,762 | 81 | 5,5 | 20,7 |
| Ridge (resíduo) | baseline | 0,0451 | [0,0348; 0,0553] | 0,0017 | 0,9467 | 0,587 | 63 | 0,1 | 0,0 |
| White-Box | white-box | 0,0702 | [0,0441; 0,0930] | 0,0000 | 0,8710 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,1007 | [0,0723; 0,1427] | 0,0343 | 0,7050 | -1,286 | -137 | 31,6 | 36,7 |

**LOGO (16 ensaios)** (teto do ganho = 0,937)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GBDT [B] | hibrido | 0,0215 | [0,0170; 0,0261] | – | 0,9878 | 0,906 | 97 | 0,0 | 9,8 |
| Caixa-preta GP | caixa-preta | 0,0217 | [0,0175; 0,0258] | – | 0,9877 | 0,904 | 96 | 15,2 | 25,9 |
| Caixa-preta GBDT | caixa-preta | 0,0220 | [0,0185; 0,0254] | – | 0,9873 | 0,901 | 96 | 8,0 | 6,2 |
| Híbrido Random Forest [B] | hibrido | 0,0225 | [0,0179; 0,0276] | – | 0,9867 | 0,897 | 96 | 5,4 | 14,3 |
| Híbrido GP [B] | hibrido | 0,0233 | [0,0176; 0,0289] | – | 0,9858 | 0,890 | 95 | 10,7 | 2,7 |
| Híbrido Extra Trees [B] | hibrido | 0,0249 | [0,0176; 0,0316] | – | 0,9837 | 0,874 | 93 | 3,6 | 5,4 |
| Híbrido SVR [B] | hibrido | 0,0249 | [0,0198; 0,0299] | – | 0,9837 | 0,874 | 93 | 7,1 | 11,6 |
| Caixa-preta Extra Trees | caixa-preta | 0,0255 | [0,0223; 0,0285] | – | 0,9830 | 0,868 | 93 | 6,2 | 8,9 |
| Híbrido MLP tanh [B] | hibrido | 0,0257 | [0,0186; 0,0325] | – | 0,9827 | 0,866 | 92 | 8,0 | 6,2 |
| Tabela (η, t) | baseline | 0,0338 | [0,0231; 0,0452] | – | 0,9701 | 0,768 | 82 | 5,4 | 19,6 |
| Ridge (resíduo) | baseline | 0,0449 | [0,0346; 0,0552] | – | 0,9471 | 0,590 | 63 | 0,0 | 0,0 |
| White-Box | white-box | 0,0702 | [0,0441; 0,0930] | – | 0,8710 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,0795 | [0,0572; 0,1002] | – | 0,8345 | -0,282 | -30 | 31,2 | 34,8 |

**Deixa-um-η-fora (extrapolação em η)** (teto do ganho = 0,937)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido Random Forest [B] | hibrido | 0,0476 | [0,0350; 0,0597] | – | 0,9405 | 0,539 | 57 | 2,7 | 14,3 |
| Híbrido GBDT [B] | hibrido | 0,0505 | [0,0378; 0,0629] | – | 0,9331 | 0,481 | 51 | 25,0 | 13,4 |
| Híbrido MLP tanh [B] | hibrido | 0,0535 | [0,0403; 0,0661] | – | 0,9250 | 0,419 | 45 | 15,2 | 29,5 |
| Híbrido Extra Trees [B] | hibrido | 0,0593 | [0,0451; 0,0724] | – | 0,9078 | 0,286 | 30 | 1,8 | 7,1 |
| Ridge (resíduo) | baseline | 0,0610 | [0,0462; 0,0747] | – | 0,9025 | 0,244 | 26 | 0,0 | 0,0 |
| Tabela (η, t) | baseline | 0,0630 | [0,0445; 0,0799] | – | 0,8959 | 0,193 | 21 | 0,0 | 25,0 |
| Híbrido GP [B] | hibrido | 0,0691 | [0,0536; 0,0842] | – | 0,8750 | 0,031 | 3 | 0,0 | 3,6 |
| White-Box | white-box | 0,0702 | [0,0441; 0,0930] | – | 0,8710 | 0,000 | 0 | 0,0 | 0,0 |
| Híbrido SVR [B] | hibrido | 0,0711 | [0,0577; 0,0833] | – | 0,8673 | -0,028 | -3 | 5,4 | 29,5 |
| Caixa-preta Extra Trees | caixa-preta | 0,0787 | [0,0587; 0,0964] | – | 0,8375 | -0,259 | -28 | 25,0 | 12,5 |
| Caixa-preta GP | caixa-preta | 0,1289 | [0,0820; 0,1661] | – | 0,5646 | -2,374 | -253 | 46,4 | 17,9 |
| Caixa-preta GBDT | caixa-preta | 0,1834 | [0,1059; 0,2484] | – | 0,1181 | -5,834 | -622 | 25,0 | 6,2 |
| Caixa-preta SVR | caixa-preta | 0,4916 | [0,2521; 0,6831] | – | -5,3331 | -48,079 | -5129 | 42,9 | 25,9 |

**Deixa-um-C_A0-fora** (teto do ganho = 0,937)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GP [B] | hibrido | 0,0222 | [0,0162; 0,0284] | – | 0,9870 | 0,899 | 96 | 7,1 | 0,9 |
| Híbrido MLP tanh [B] | hibrido | 0,0233 | [0,0169; 0,0298] | – | 0,9858 | 0,890 | 95 | 3,6 | 4,5 |
| Caixa-preta GP | caixa-preta | 0,0254 | [0,0201; 0,0303] | – | 0,9831 | 0,869 | 93 | 17,0 | 24,1 |
| Híbrido Extra Trees [B] | hibrido | 0,0261 | [0,0181; 0,0333] | – | 0,9821 | 0,861 | 92 | 3,6 | 6,2 |
| Caixa-preta Extra Trees | caixa-preta | 0,0265 | [0,0231; 0,0295] | – | 0,9817 | 0,858 | 92 | 6,2 | 11,6 |
| Híbrido SVR [B] | hibrido | 0,0275 | [0,0186; 0,0351] | – | 0,9802 | 0,847 | 90 | 8,0 | 14,3 |
| Híbrido GBDT [B] | hibrido | 0,0280 | [0,0190; 0,0360] | – | 0,9795 | 0,841 | 90 | 4,5 | 9,8 |
| Caixa-preta GBDT | caixa-preta | 0,0288 | [0,0233; 0,0340] | – | 0,9782 | 0,831 | 89 | 12,5 | 3,6 |
| Híbrido Random Forest [B] | hibrido | 0,0305 | [0,0197; 0,0401] | – | 0,9756 | 0,811 | 87 | 4,5 | 10,7 |
| Tabela (η, t) | baseline | 0,0338 | [0,0231; 0,0452] | – | 0,9701 | 0,768 | 82 | 5,4 | 19,6 |
| Ridge (resíduo) | baseline | 0,0403 | [0,0312; 0,0488] | – | 0,9575 | 0,670 | 72 | 0,0 | 0,0 |
| White-Box | white-box | 0,0702 | [0,0441; 0,0930] | – | 0,8710 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,1774 | [0,1139; 0,2298] | – | 0,1755 | -5,390 | -575 | 17,9 | 31,2 |

**Extrapolação temporal (t<=5 -> 15 min)** (teto do ganho = 0,933)

| Modelo | Tipo | RMSE | IC95% | dp entre partições | R² (t>0) | Ganho vs WB | % do teto | bruto > teto (%) | bruto quedas (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GBDT [B] | hibrido | 0,0105 | [0,0078; 0,0132] | – | 0,9973 | 0,976 | 105 | 0,0 | 0,0 |
| Híbrido Extra Trees [B] | hibrido | 0,0106 | [0,0075; 0,0136] | – | 0,9972 | 0,976 | 105 | 0,0 | 0,0 |
| Híbrido Random Forest [B] | hibrido | 0,0113 | [0,0077; 0,0148] | – | 0,9968 | 0,972 | 104 | 0,0 | 0,0 |
| Caixa-preta Extra Trees | caixa-preta | 0,0123 | [0,0094; 0,0153] | – | 0,9962 | 0,967 | 104 | 6,2 | 0,0 |
| Híbrido GP [B] | hibrido | 0,0139 | [0,0105; 0,0172] | – | 0,9952 | 0,958 | 103 | 0,0 | 0,0 |
| Caixa-preta GBDT | caixa-preta | 0,0198 | [0,0145; 0,0248] | – | 0,9903 | 0,915 | 98 | 12,5 | 0,0 |
| Híbrido MLP tanh [B] | hibrido | 0,0206 | [0,0127; 0,0279] | – | 0,9895 | 0,908 | 97 | 37,5 | 0,0 |
| Híbrido SVR [B] | hibrido | 0,0656 | [0,0495; 0,0803] | – | 0,8935 | 0,069 | 7 | 6,2 | 0,0 |
| Ridge (resíduo) | baseline | 0,0663 | [0,0340; 0,0913] | – | 0,8913 | 0,050 | 5 | 68,8 | 0,0 |
| White-Box | white-box | 0,0680 | [0,0337; 0,0940] | – | 0,8856 | 0,000 | 0 | 0,0 | 0,0 |
| Tabela (η, t) | baseline | 0,0680 | [0,0337; 0,0940] | – | 0,8856 | 0,000 | 0 | 0,0 | 0,0 |
| Caixa-preta GP | caixa-preta | 0,1964 | [0,1502; 0,2353] | – | 0,0455 | -7,344 | -787 | 0,0 | 0,0 |
| Caixa-preta SVR | caixa-preta | 0,2895 | [0,2294; 0,3441] | – | -1,0743 | -17,132 | -1836 | 12,5 | 0,0 |

### T06. Extrapolação em η: RMSE no nível deixado fora do treino

| Modelo | RMSE total | η = 0,5 ausente | η = 1,0 ausente | η = 1,5 ausente | η = 3,1 ausente |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Híbrido Random Forest [B] | 0,0476 | 0,0225 | 0,0741 | 0,0452 | 0,0323 |
| Híbrido GBDT [B] | 0,0505 | 0,0258 | 0,0782 | 0,0483 | 0,0330 |
| Híbrido MLP tanh [B] | 0,0535 | 0,0257 | 0,0811 | 0,0572 | 0,0306 |
| Híbrido Extra Trees [B] | 0,0593 | 0,0292 | 0,0833 | 0,0713 | 0,0345 |
| Ridge (resíduo) | 0,0610 | 0,0655 | 0,0895 | 0,0392 | 0,0323 |
| Tabela (η, t) | 0,0630 | 0,0520 | 0,1021 | 0,0278 | 0,0446 |
| Híbrido GP [B] | 0,0691 | 0,0455 | 0,1025 | 0,0703 | 0,0394 |
| White-Box | 0,0702 | 0,0245 | 0,1234 | 0,0526 | 0,0330 |
| Híbrido SVR [B] | 0,0711 | 0,0795 | 0,0856 | 0,0687 | 0,0433 |
| Caixa-preta Extra Trees | 0,0787 | 0,0258 | 0,1243 | 0,0841 | 0,0399 |
| Caixa-preta GP | 0,1289 | 0,0258 | 0,1017 | 0,0583 | 0,2281 |
| Caixa-preta GBDT | 0,1834 | 0,0258 | 0,3471 | 0,1072 | 0,0440 |
| Caixa-preta SVR | 0,4916 | 0,0234 | 0,1229 | 0,0687 | 0,9727 |

### T07. Híbrido contra a caixa-preta do MESMO modelo (diferença pareada, bootstrap de ensaios)

| Protocolo | Híbrido | Caixa-preta | Δ RMSE (híbrido − caixa-preta) | IC95% | Empate (IC contém 0) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta GBDT | 0,0000 | [-0,0039; 0,0032] | sim |
| GKF8 (partição do TCC) | Híbrido GP [B] | Caixa-preta GP | 0,0010 | [-0,0036; 0,0056] | sim |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0008 | [-0,0076; 0,0053] | sim |
| GKF8 (partição do TCC) | Híbrido SVR [B] | Caixa-preta SVR | -0,0616 | [-0,0805; -0,0432] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0005 | [-0,0048; 0,0030] | sim |
| GKF8 repetido (20 partições) | Híbrido GP [B] | Caixa-preta GP | 0,0007 | [-0,0035; 0,0052] | sim |
| GKF8 repetido (20 partições) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0007 | [-0,0068; 0,0050] | sim |
| GKF8 repetido (20 partições) | Híbrido SVR [B] | Caixa-preta SVR | -0,0819 | [-0,1196; -0,0486] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0005 | [-0,0047; 0,0031] | sim |
| LOGO (16 ensaios) | Híbrido GP [B] | Caixa-preta GP | 0,0016 | [-0,0029; 0,0060] | sim |
| LOGO (16 ensaios) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0006 | [-0,0065; 0,0048] | sim |
| LOGO (16 ensaios) | Híbrido SVR [B] | Caixa-preta SVR | -0,0545 | [-0,0762; -0,0323] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido GP [B] | Caixa-preta GP | -0,0598 | [-0,1048; -0,0075] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0194 | [-0,0331; -0,0073] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,1329 | [-0,1866; -0,0653] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido SVR [B] | Caixa-preta SVR | -0,4204 | [-0,6153; -0,1739] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GP | -0,0032 | [-0,0086; 0,0024] | sim |
| Deixa-um-C_A0-fora | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0003 | [-0,0073; 0,0057] | sim |
| Deixa-um-C_A0-fora | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0008 | [-0,0064; 0,0038] | sim |
| Deixa-um-C_A0-fora | Híbrido SVR [B] | Caixa-preta SVR | -0,1499 | [-0,2046; -0,0861] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0093 | [-0,0128; -0,0057] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GP [B] | Caixa-preta GP | -0,1825 | [-0,2190; -0,1385] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido Extra Trees [B] | Caixa-preta Extra Trees | -0,0018 | [-0,0038; 0,0000] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido SVR [B] | Caixa-preta SVR | -0,2239 | [-0,2779; -0,1698] | não |

### T08. Melhor híbrido de cada protocolo contra os demais modelos (Δ < 0: o melhor híbrido erra menos)

| Protocolo | Melhor híbrido | Contra | Δ RMSE | IC95% | Empate |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | White-Box | -0,0483 | [-0,0727; -0,0212] | não |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Tabela (η, t) | -0,0119 | [-0,0211; -0,0038] | não |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Ridge (resíduo) | -0,0229 | [-0,0342; -0,0123] | não |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Híbrido GP [B] | -0,0013 | [-0,0032; 0,0012] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Híbrido Extra Trees [B] | -0,0029 | [-0,0077; 0,0016] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Híbrido Random Forest [B] | -0,0008 | [-0,0025; 0,0009] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Híbrido SVR [B] | -0,0018 | [-0,0055; 0,0023] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Híbrido MLP tanh [B] | -0,0025 | [-0,0058; 0,0001] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta GP | -0,0003 | [-0,0039; 0,0031] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta Extra Trees | -0,0037 | [-0,0084; 0,0009] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta GBDT | 0,0000 | [-0,0039; 0,0032] | sim |
| GKF8 (partição do TCC) | Híbrido GBDT [B] | Caixa-preta SVR | -0,0634 | [-0,0826; -0,0435] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | White-Box | -0,0485 | [-0,0722; -0,0216] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Tabela (η, t) | -0,0125 | [-0,0216; -0,0047] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Ridge (resíduo) | -0,0234 | [-0,0346; -0,0127] | não |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido GP [B] | -0,0012 | [-0,0031; 0,0010] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido Extra Trees [B] | -0,0037 | [-0,0091; 0,0007] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido Random Forest [B] | -0,0010 | [-0,0022; 0,0002] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido SVR [B] | -0,0025 | [-0,0053; 0,0007] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Híbrido MLP tanh [B] | -0,0022 | [-0,0063; 0,0014] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GP | -0,0005 | [-0,0043; 0,0031] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta Extra Trees | -0,0045 | [-0,0092; 0,0001] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0005 | [-0,0048; 0,0030] | sim |
| GKF8 repetido (20 partições) | Híbrido GBDT [B] | Caixa-preta SVR | -0,0844 | [-0,1220; -0,0499] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | White-Box | -0,0486 | [-0,0726; -0,0216] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Tabela (η, t) | -0,0123 | [-0,0213; -0,0044] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Ridge (resíduo) | -0,0234 | [-0,0348; -0,0126] | não |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido GP [B] | -0,0017 | [-0,0037; 0,0009] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido Extra Trees [B] | -0,0034 | [-0,0087; 0,0013] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido Random Forest [B] | -0,0010 | [-0,0024; 0,0003] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido SVR [B] | -0,0034 | [-0,0078; 0,0010] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Híbrido MLP tanh [B] | -0,0042 | [-0,0103; 0,0010] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GP | -0,0002 | [-0,0039; 0,0035] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta Extra Trees | -0,0039 | [-0,0085; 0,0006] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0005 | [-0,0047; 0,0031] | sim |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Caixa-preta SVR | -0,0579 | [-0,0800; -0,0351] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | White-Box | -0,0225 | [-0,0342; -0,0078] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Tabela (η, t) | -0,0154 | [-0,0225; -0,0060] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Ridge (resíduo) | -0,0134 | [-0,0222; -0,0048] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido GP [B] | -0,0214 | [-0,0273; -0,0146] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido Extra Trees [B] | -0,0117 | [-0,0195; -0,0049] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido GBDT [B] | -0,0029 | [-0,0045; -0,0014] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido SVR [B] | -0,0235 | [-0,0376; -0,0115] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Híbrido MLP tanh [B] | -0,0058 | [-0,0114; -0,0007] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta GP | -0,0812 | [-0,1220; -0,0332] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta Extra Trees | -0,0311 | [-0,0427; -0,0184] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta GBDT | -0,1358 | [-0,1904; -0,0678] | não |
| Deixa-um-η-fora (extrapolação em η) | Híbrido Random Forest [B] | Caixa-preta SVR | -0,4439 | [-0,6393; -0,2017] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | White-Box | -0,0479 | [-0,0735; -0,0204] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Tabela (η, t) | -0,0115 | [-0,0209; -0,0036] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Ridge (resíduo) | -0,0180 | [-0,0285; -0,0090] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido Extra Trees [B] | -0,0039 | [-0,0106; 0,0009] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido GBDT [B] | -0,0057 | [-0,0108; 0,0001] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido Random Forest [B] | -0,0083 | [-0,0143; -0,0004] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido SVR [B] | -0,0052 | [-0,0119; 0,0016] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Híbrido MLP tanh [B] | -0,0010 | [-0,0028; 0,0008] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GP | -0,0032 | [-0,0086; 0,0024] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta Extra Trees | -0,0042 | [-0,0102; 0,0018] | sim |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta GBDT | -0,0066 | [-0,0120; -0,0018] | não |
| Deixa-um-C_A0-fora | Híbrido GP [B] | Caixa-preta SVR | -0,1551 | [-0,2102; -0,0919] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | White-Box | -0,0575 | [-0,0843; -0,0230] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Tabela (η, t) | -0,0575 | [-0,0843; -0,0230] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Ridge (resíduo) | -0,0558 | [-0,0816; -0,0231] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Híbrido GP [B] | -0,0034 | [-0,0077; -0,0000] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Híbrido Extra Trees [B] | -0,0001 | [-0,0016; 0,0020] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Híbrido Random Forest [B] | -0,0008 | [-0,0027; 0,0011] | sim |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Híbrido SVR [B] | -0,0551 | [-0,0708; -0,0382] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Híbrido MLP tanh [B] | -0,0101 | [-0,0177; -0,0016] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta GP | -0,1859 | [-0,2248; -0,1402] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta Extra Trees | -0,0018 | [-0,0036; -0,0000] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta GBDT | -0,0093 | [-0,0128; -0,0057] | não |
| Extrapolação temporal (t<=5 -> 15 min) | Híbrido GBDT [B] | Caixa-preta SVR | -0,2790 | [-0,3335; -0,2184] | não |

### T09. Curva de aprendizado: RMSE em 2 ensaios fora do treino, por nº de ensaios de treino (30 sorteios)

| Modelo | n = 4 | n = 6 | n = 8 | n = 10 | n = 12 | n = 14 |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| White-Box | 0,0700 | 0,0700 | 0,0700 | 0,0700 | 0,0700 | 0,0700 |
| Tabela (η, t) | 0,0546 | 0,0425 | 0,0422 | 0,0372 | 0,0344 | 0,0338 |
| Híbrido Extra Trees [B] | 0,0506 | 0,0369 | 0,0326 | 0,0274 | 0,0236 | 0,0236 |
| Híbrido GBDT [B] | 0,0518 | 0,0344 | 0,0311 | 0,0254 | 0,0210 | 0,0211 |
| Híbrido GP [B] | 0,0541 | 0,0414 | 0,0339 | 0,0266 | 0,0225 | 0,0220 |
| Caixa-preta Extra Trees | 0,0874 | 0,0462 | 0,0355 | 0,0282 | 0,0282 | 0,0260 |
| Caixa-preta GBDT | 0,1240 | 0,0573 | 0,0519 | 0,0251 | 0,0238 | 0,0224 |
| Caixa-preta GP | 0,0924 | 0,0533 | 0,0369 | 0,0260 | 0,0241 | 0,0219 |

### T10. Ablação dos atributos do ML (RMSE em t>0) (conjunto C omitido: exige [H⁺], indisponível com a força motriz C_Af)

| Protocolo | Modelo | A (t, η, C_A0) | B (A + X_PBM) |
| :--- | ---: | ---: | ---: |
| GKF8 repetido (10) | GP | 0,0279 | 0,0229 |
| GKF8 repetido (10) | Extra Trees | 0,0290 | 0,0254 |
| GKF8 repetido (10) | GBDT | 0,0346 | 0,0218 |
| GKF8 repetido (10) | SVR | 0,0574 | 0,0245 |
| Deixa-um-η-fora | GP | 0,0871 | 0,0691 |
| Deixa-um-η-fora | Extra Trees | 0,0801 | 0,0593 |
| Deixa-um-η-fora | GBDT | 0,0893 | 0,0505 |
| Deixa-um-η-fora | SVR | 0,1200 | 0,0711 |

### T11. Sensibilidade ao fator ν (RMSE em t>0)

| Protocolo | Modelo | nu estimado na dobra (comparativo) | nu a priori (1.06, comparativo) | nu = 1 (White-Box principal) |
| :--- | ---: | ---: | ---: | ---: |
| GKF8 (partição do TCC) | White-Box | 0,0521 | 0,0542 | 0,0702 |
| GKF8 (partição do TCC) | Híbrido GP [B] | 0,0259 | 0,0250 | 0,0231 |
| GKF8 (partição do TCC) | Híbrido Extra Trees [B] | 0,0269 | 0,0260 | 0,0247 |
| Deixa-um-η-fora | White-Box | 0,0524 | 0,0542 | 0,0702 |
| Deixa-um-η-fora | Híbrido GP [B] | 0,0527 | 0,0509 | 0,0691 |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | 0,0360 | 0,0381 | 0,0593 |

### T12. O mesmo ML sobre degraus do White-Box (RMSE em t>0)

| Protocolo | Modelo | 1. Tese (α = 0) | 2. White-Box do TCC (V0, nu = 1) | 3. V0 + nu (estimado na dobra, comparativo) | 4. CAf, nu = 1 (principal) |
| :--- | ---: | ---: | ---: | ---: | ---: |
| GKF8 repetido (10) | White-Box | 0,0981 | 0,0443 | 0,0393 | 0,0702 |
| GKF8 repetido (10) | Híbrido GP [B] | 0,0434 | 0,0236 | 0,0236 | 0,0229 |
| GKF8 repetido (10) | Híbrido Extra Trees [B] | 0,0351 | 0,0201 | 0,0197 | 0,0254 |
| GKF8 repetido (10) | Caixa-preta GP | 0,0225 | 0,0225 | 0,0219 | 0,0218 |
| GKF8 repetido (10) | Caixa-preta Extra Trees | 0,0263 | 0,0263 | 0,0261 | 0,0261 |
| Deixa-um-η-fora | White-Box | 0,0981 | 0,0443 | 0,0394 | 0,0702 |
| Deixa-um-η-fora | Híbrido GP [B] | 0,0794 | 0,0418 | 0,0392 | 0,0691 |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | 0,0661 | 0,0341 | 0,0285 | 0,0593 |
| Deixa-um-η-fora | Caixa-preta GP | 0,1309 | 0,1309 | 0,1288 | 0,1289 |
| Deixa-um-η-fora | Caixa-preta Extra Trees | 0,0820 | 0,0820 | 0,0786 | 0,0787 |
| Poucos dados (n = 4 ensaios de treino) | White-Box | 0,0904 | 0,0450 | 0,0403 | 0,0668 |
| Poucos dados (n = 4 ensaios de treino) | Híbrido GP [B] | 0,0702 | 0,0325 | 0,0321 | 0,0506 |
| Poucos dados (n = 4 ensaios de treino) | Híbrido Extra Trees [B] | 0,0632 | 0,0346 | 0,0295 | 0,0498 |
| Poucos dados (n = 4 ensaios de treino) | Caixa-preta GP | 0,0703 | 0,0703 | 0,0660 | 0,0655 |
| Poucos dados (n = 4 ensaios de treino) | Caixa-preta Extra Trees | 0,0871 | 0,0871 | 0,0821 | 0,0816 |
| Poucos dados (n = 8 ensaios de treino) | White-Box | 0,0904 | 0,0450 | 0,0395 | 0,0668 |
| Poucos dados (n = 8 ensaios de treino) | Híbrido GP [B] | 0,0545 | 0,0261 | 0,0261 | 0,0293 |
| Poucos dados (n = 8 ensaios de treino) | Híbrido Extra Trees [B] | 0,0476 | 0,0225 | 0,0206 | 0,0305 |
| Poucos dados (n = 8 ensaios de treino) | Caixa-preta GP | 0,0272 | 0,0272 | 0,0268 | 0,0266 |
| Poucos dados (n = 8 ensaios de treino) | Caixa-preta Extra Trees | 0,0335 | 0,0335 | 0,0301 | 0,0302 |

### T13. Cobertura empírica dos intervalos de 95% (t>0, pontos fora-da-dobra)

| Protocolo | Modelo | Método | Cobertura | Largura média | n |
| :--- | ---: | ---: | ---: | ---: | ---: |
| Deixa-um-η-fora | Híbrido Extra Trees [B] | Conforme (LOGO interno) | 57% | 0,0940 | 112 |
| Deixa-um-η-fora | Híbrido GBDT [B] | Conforme (LOGO interno) | 70% | 0,0939 | 112 |
| Deixa-um-η-fora | Híbrido GP [B] | Conforme (LOGO interno) | 45% | 0,0996 | 112 |
| Deixa-um-η-fora | Híbrido GP [B] | GP (analítico) | 57% | 0,1439 | 112 |
| LOGO (16 ensaios) | Híbrido Extra Trees [B] | Conforme (LOGO interno) | 88% | 0,0914 | 112 |
| LOGO (16 ensaios) | Híbrido GBDT [B] | Conforme (LOGO interno) | 90% | 0,0715 | 112 |
| LOGO (16 ensaios) | Híbrido GP [B] | Conforme (LOGO interno) | 89% | 0,0767 | 112 |
| LOGO (16 ensaios) | Híbrido GP [B] | GP (analítico) | 88% | 0,0652 | 112 |

### T14. Resíduo do White-Box principal aos 15 min

| Ensaio | η | η real | C_A0 | X_exp (15 min) | X_PBM (15 min) | ΔX |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,535 | 0,1 | 0,50 | 0,5000 | -0,000 |
| 2 | 3,1 | 3,290 | 0,1 | 1,00 | 1,0000 | 0,000 |
| 3 | 0,5 | 0,535 | 0,5 | 0,50 | 0,5000 | -0,000 |
| 4 | 0,5 | 0,535 | 1,0 | 0,50 | 0,5000 | -0,000 |
| 5 | 0,5 | 0,535 | 1,5 | 0,49 | 0,5000 | -0,010 |
| 6 | 1,0 | 1,069 | 0,1 | 0,87 | 0,9774 | -0,107 |
| 7 | 3,1 | 3,192 | 0,5 | 1,01 | 1,0000 | 0,010 |
| 8 | 1,0 | 1,069 | 0,5 | 0,87 | 0,9946 | -0,125 |
| 9 | 1,0 | 1,069 | 1,0 | 0,85 | 0,9972 | -0,147 |
| 10 | 1,0 | 1,069 | 1,5 | 0,85 | 0,9981 | -0,148 |
| 11 | 1,5 | 1,584 | 0,1 | 0,98 | 1,0000 | -0,020 |
| 12 | 3,1 | 3,216 | 1,0 | 1,00 | 1,0000 | 0,000 |
| 13 | 1,5 | 1,608 | 0,5 | 0,97 | 1,0000 | -0,030 |
| 14 | 1,5 | 1,602 | 1,0 | 0,97 | 1,0000 | -0,030 |
| 15 | 1,5 | 1,604 | 1,5 | 0,97 | 1,0000 | -0,030 |
| 16 | 3,1 | 3,208 | 1,5 | 1,00 | 1,0000 | 0,000 |

### T15. Escalas de comprimento do GP (ARD) ajustadas com todos os dados

| Atributos | Entrada | Escala de comprimento (σ das entradas) | Ruído branco (norm.) |
| :--- | ---: | ---: | ---: |
| A | tempo_min | 1,228 | 0,137 |
| A | razao_molar_eta | 0,050 | 0,137 |
| A | C_acid_0_mol_L | 3,746 | 0,137 |
| B | tempo_min | 4,381 | 0,074 |
| B | razao_molar_eta | 0,510 | 0,074 |
| B | C_acid_0_mol_L | 6,254 | 0,074 |
| B | X_pbm | 0,317 | 0,074 |

### T16. Conversão do White-Box principal e resíduo, por ensaio

**X_PBM (conversão do White-Box principal)**

| Ensaio | η | C_A0 | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 0,00 | 0,49 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 2 | 3,1 | 0,1 | 0,00 | 0,85 | 0,95 | 0,99 | 1,00 | 1,00 | 1,00 | 1,00 |
| 3 | 0,5 | 0,5 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 4 | 0,5 | 1,0 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 5 | 0,5 | 1,5 | 0,00 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 | 0,50 |
| 6 | 1,0 | 0,1 | 0,00 | 0,72 | 0,82 | 0,89 | 0,92 | 0,93 | 0,94 | 0,98 |
| 7 | 3,1 | 0,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 8 | 1,0 | 0,5 | 0,00 | 0,90 | 0,94 | 0,97 | 0,98 | 0,98 | 0,99 | 0,99 |
| 9 | 1,0 | 1,0 | 0,00 | 0,94 | 0,97 | 0,98 | 0,99 | 0,99 | 0,99 | 1,00 |
| 10 | 1,0 | 1,5 | 0,00 | 0,96 | 0,98 | 0,99 | 0,99 | 0,99 | 0,99 | 1,00 |
| 11 | 1,5 | 0,1 | 0,00 | 0,79 | 0,90 | 0,97 | 0,99 | 0,99 | 1,00 | 1,00 |
| 12 | 3,1 | 1,0 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 13 | 1,5 | 0,5 | 0,00 | 0,98 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 14 | 1,5 | 1,0 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 15 | 1,5 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |
| 16 | 3,1 | 1,5 | 0,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 | 1,00 |

**ΔX = X_exp − X_PBM**

| Ensaio | η | C_A0 | 0 min | 0,5 min | 1 min | 2 min | 3 min | 4 min | 5 min | 15 min |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0,5 | 0,1 | 0,00 | -0,07 | -0,04 | -0,04 | -0,01 | 0,01 | -0,01 | -0,00 |
| 2 | 3,1 | 0,1 | 0,00 | -0,03 | -0,05 | -0,04 | -0,03 | -0,03 | -0,01 | 0,00 |
| 3 | 0,5 | 0,5 | 0,00 | -0,05 | -0,02 | -0,01 | 0,02 | -0,00 | -0,01 | -0,00 |
| 4 | 0,5 | 1,0 | 0,00 | -0,00 | -0,03 | 0,01 | -0,03 | 0,03 | 0,01 | -0,00 |
| 5 | 0,5 | 1,5 | 0,00 | -0,01 | -0,01 | -0,03 | -0,01 | -0,01 | -0,03 | -0,01 |
| 6 | 1,0 | 0,1 | 0,00 | 0,03 | -0,04 | -0,09 | -0,07 | -0,06 | -0,09 | -0,11 |
| 7 | 3,1 | 0,5 | 0,00 | -0,12 | -0,07 | -0,02 | -0,01 | 0,00 | -0,01 | 0,01 |
| 8 | 1,0 | 0,5 | 0,00 | -0,07 | -0,12 | -0,12 | -0,12 | -0,11 | -0,13 | -0,12 |
| 9 | 1,0 | 1,0 | 0,00 | -0,11 | -0,13 | -0,13 | -0,14 | -0,16 | -0,15 | -0,15 |
| 10 | 1,0 | 1,5 | 0,00 | -0,12 | -0,15 | -0,16 | -0,16 | -0,16 | -0,15 | -0,15 |
| 11 | 1,5 | 0,1 | 0,00 | 0,07 | -0,05 | -0,09 | -0,02 | -0,05 | -0,04 | -0,02 |
| 12 | 3,1 | 1,0 | 0,00 | -0,04 | 0,00 | 0,00 | -0,01 | -0,01 | -0,01 | 0,00 |
| 13 | 1,5 | 0,5 | 0,00 | -0,06 | -0,08 | -0,05 | -0,05 | -0,05 | -0,05 | -0,03 |
| 14 | 1,5 | 1,0 | 0,00 | -0,09 | -0,04 | -0,05 | -0,03 | -0,05 | -0,02 | -0,03 |
| 15 | 1,5 | 1,5 | 0,00 | -0,09 | -0,06 | -0,05 | -0,06 | -0,04 | -0,04 | -0,03 |
| 16 | 3,1 | 1,5 | 0,00 | -0,02 | -0,02 | -0,04 | -0,01 | -0,01 | 0,00 | 0,00 |

### T17. Partição comum do relatório (teste cego 7, 8, 14; CV nas 4 dobras da Tabela B.3; 10 sementes)

Modelo selecionado só pela CV: **Híbrido GP [B]**. Erros em p.p. de conversão, t > 0; inconsistência = % dos nós da malha de 0,25 min com conversão decrescente ou acima do teto, antes (bruta) e depois do acoplador.

| Conjunto | Ensaios | Pontos | R² | RMSE (p.p.) | MAE (p.p.) | Erro máximo (p.p.) | Inconsistência bruta (%) | Inconsistência final (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Treino | 13 | 91 | 0,997 | 1,16 | 0,85 | 3,64 | 0,9 | 0,0 |
| Teste cego | 3 | 21 | 0,817 | 2,58 | 1,76 | 8,51 | 1,1 | 0,0 |
| Global | 16 | 112 | 0,994 | 1,53 | 1,02 | 8,51 | 0,9 | 0,0 |

Todos os modelos, ordenados pelo erro da CV (o teste cego não participa da escolha):

| Modelo | Tipo | CV fora da dobra (p.p.) | Teste RMSE (p.p.) | Teste R² | E7 (p.p.) | E8 (p.p.) | E14 (p.p.) | Inconsist. bruta (%) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Híbrido GP [B] | hibrido | 1,97 ± 0,00 | 2,58 | 0,817 | 3,77 | 1,62 | 1,75 | 1,1 |
| Híbrido HistGBDT [B] | hibrido | 2,06 ± 0,00 | 1,98 | 0,892 | 2,89 | 1,25 | 1,37 | 6,1 |
| Híbrido GBDT [B] | hibrido | 2,10 ± 0,00 | 2,22 | 0,865 | 3,16 | 1,80 | 1,24 | 2,8 |
| Híbrido Random Forest [B] | hibrido | 2,16 ± 0,03 | 2,39 | 0,843 | 3,43 | 1,78 | 1,48 | 3,9 |
| Caixa-preta MLP tanh | caixa-preta | 2,17 ± 0,03 | 1,92 | 0,899 | 1,30 | 2,87 | 1,05 | 14,4 |
| Híbrido AdaBoost [B] | hibrido | 2,32 ± 0,03 | 2,89 | 0,770 | 4,12 | 2,20 | 1,79 | 1,1 |
| Caixa-preta GBDT | caixa-preta | 2,45 ± 0,00 | 2,12 | 0,876 | 2,51 | 1,90 | 1,89 | 0,0 |
| Caixa-preta HistGBDT | caixa-preta | 2,52 ± 0,00 | 2,06 | 0,883 | 2,61 | 1,73 | 1,71 | 16,7 |
| Híbrido Extra Trees [B] | hibrido | 2,59 ± 0,08 | 2,36 | 0,847 | 3,50 | 1,20 | 1,74 | 0,0 |
| Caixa-preta Random Forest | caixa-preta | 2,60 ± 0,04 | 2,39 | 0,843 | 2,44 | 2,87 | 1,73 | 0,0 |
| Caixa-preta GP | caixa-preta | 2,67 ± 0,00 | 1,96 | 0,894 | 1,31 | 2,47 | 1,92 | 33,9 |
| Caixa-preta GP (RBF) | caixa-preta | 2,79 ± 0,09 | 1,99 | 0,891 | 1,30 | 2,44 | 2,06 | 25,0 |
| Caixa-preta Extra Trees | caixa-preta | 3,06 ± 0,09 | 2,67 | 0,804 | 1,37 | 3,22 | 3,02 | 0,0 |
| Híbrido MLP tanh [B] | hibrido | 3,08 ± 0,25 | 2,23 | 0,863 | 3,30 | 0,85 | 1,83 | 3,9 |
| Híbrido k-NN [B] | hibrido | 3,27 ± 0,00 | 5,79 | 0,076 | 3,61 | 5,54 | 7,55 | 6,1 |
| Caixa-preta AdaBoost | caixa-preta | 3,34 ± 0,04 | 3,47 | 0,668 | 4,82 | 2,68 | 2,39 | 0,0 |
| Tabela (η, t) | baseline | 3,74 ± 0,00 | 2,48 | 0,831 | 3,71 | 1,12 | 1,84 | 11,1 |
| Híbrido GP (RBF) [B] | hibrido | 3,98 ± 1,58 | 2,95 | 0,761 | 3,71 | 3,07 | 1,69 | 10,6 |
| Híbrido SVR [B] | hibrido | 4,20 ± 0,00 | 1,81 | 0,910 | 2,44 | 1,28 | 1,48 | 0,0 |
| Ridge (resíduo) | baseline | 4,30 ± 0,00 | 4,00 | 0,559 | 3,82 | 4,97 | 2,96 | 0,0 |
| Híbrido Kernel Ridge [B] | hibrido | 4,48 ± 0,00 | 1,81 | 0,909 | 2,53 | 1,07 | 1,54 | 0,0 |
| Híbrido Ridge polinomial [B] | hibrido | 4,59 ± 0,00 | 3,73 | 0,616 | 3,72 | 3,82 | 3,66 | 3,3 |
| White-Box | white-box | 6,82 ± 0,00 | 7,80 | -0,674 | 5,21 | 11,49 | 4,84 | 0,0 |
| Caixa-preta Ridge polinomial | caixa-preta | 7,78 ± 0,00 | 6,92 | -0,316 | 6,40 | 6,87 | 7,44 | 45,0 |
| Caixa-preta Kernel Ridge | caixa-preta | 9,15 ± 0,00 | 5,34 | 0,215 | 4,71 | 6,35 | 4,81 | 46,1 |
| Caixa-preta k-NN | caixa-preta | 18,51 ± 0,00 | 11,13 | -2,408 | 5,60 | 13,52 | 12,55 | 6,7 |
| Caixa-preta SVR | caixa-preta | 19,66 ± 0,00 | 8,87 | -1,165 | 7,94 | 10,11 | 8,42 | 42,8 |
