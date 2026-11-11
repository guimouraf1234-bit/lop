# Trecho movido da metodologia (4.5.6) para a Seção 5 (Resultados)

> Texto removido de METODOLOGIA_HIBRIDA_PARALELA.md, sem alterações. Ajuste a numeração das equações e a intro ao colar.

O LOGO mede a capacidade de interpolar entre condições vizinhas. O deixa-um-$\eta$-fora mede o que o híbrido deve entregar na prática, prever uma razão molar não ensaiada, e é onde a arquitetura cooperativa é mais exigida. A curva de aprendizado verifica a expectativa da revisão de que a informação física pesa mais quando os dados são escassos. Os demais protocolos (GKF8 repetido, deixa-um-$C_{A0}$-fora e extrapolação temporal) estão no Apêndice B.

A métrica principal é o RMSE em $t>0$, e não o $R^2$, porque o $R^2$ global é dominado pela variância entre níveis de $\eta$, que o White-Box já explica (Equações 3.25 a 3.27; Chicco; Warrens; Jurman, 2021). O piso de ruído experimental foi estimado nos patamares das curvas, $\sigma = 0{,}0175$, e o ganho sobre o White-Box e seu teto são

$$G = 1 - \frac{MSE_{modelo}}{MSE_{WB}}, \qquad G_{max} = 1 - \frac{\sigma^2}{MSE_{WB}} \tag{4.X.X}$$

de modo que um erro abaixo de $\sigma$ não seja lido como melhoria real. A incerteza das métricas foi estimada por *bootstrap* (Efron; Tibshirani, 1993) com o ensaio como unidade (2000 reamostragens), e as comparações entre modelos usaram diferenças pareadas: dois modelos são empatados quando o IC 95% da diferença contém zero.

### Intervalos de predição

O intervalo de predição vem do desvio-padrão analítico do GP (Equação 4.X.X). A cobertura empírica é comparada com o nível nominal de 95%, em LOGO e em extrapolação em $\eta$.
