# Modelagem híbrida paralela

Esta seção apresenta os resultados do modelo híbrido paralelo cooperativo, no qual um Processo Gaussiano (GP), com núcleo Matérn 3/2, estima o resíduo entre a conversão experimental e a previsão do modelo fenomenológico. O white-box emprega a razão molar nominal, a constante cinética fixa de 18.000 µm/min e α = 0, em conformidade com a referência comum do grupo. A previsão final corresponde à média das saídas brutas de dez modelos, com sementes de 42 a 51, seguida da aplicação do acoplador físico. Os erros são expressos em pontos percentuais (p.p.) de conversão e calculados apenas nos instantes com medição experimental e t > 0.

Nos 13 ensaios de treinamento, o modelo apresentou RMSE de 1,18 p.p., erro absoluto médio de 0,86 p.p. e R² de 0,997 (Tabela 1). Essas métricas caracterizam o ajuste aos ensaios utilizados na construção do modelo e não constituem, isoladamente, uma medida de sua capacidade de previsão. O RMSE global de 1,57 p.p. também reúne pontos de treino e de teste, devendo ser interpretado como uma descrição conjunta do desempenho.

Tabela 1 – Desempenho do modelo híbrido paralelo com Processo Gaussiano na previsão da conversão de zinco, considerando somente os instantes experimentais com t > 0.

| Conjunto | Ensaios | Pontos | R² (–) | RMSE (p.p.) | MAE (p.p.) | Erro máximo (p.p.) | Inconsistência física final (%) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Treino | 13 | 91 | 0,997 | 1,18 | 0,86 | 3,64 | 0,0 |
| Teste | 3 | 21 | 0,805 | 2,66 | 1,79 | 8,89 | 0,0 |
| Global | 16 | 112 | 0,994 | 1,57 | 1,04 | 8,89 | 0,0 |

Fonte: Elaborado pelos autores.

Nos Ensaios 7, 8 e 14, reservados para teste e ausentes do treinamento, o RMSE foi de 2,66 p.p., com R² de 0,805, MAE de 1,79 p.p. e erro máximo de 8,89 p.p. O maior desvio ocorreu no Ensaio 7, em 0,5 min. Os RMSE individuais foram de 3,91 p.p. no Ensaio 7, 1,69 p.p. no Ensaio 8 e 1,76 p.p. no Ensaio 14. Frente ao RMSE de 7,80 p.p. do white-box nos mesmos 21 pontos, a correção residual reduziu o RMSE em aproximadamente 66%.

A Figura 1 compara as previsões com as medições no Ensaio 1, de treinamento, e nos três ensaios de teste. No Ensaio 8, de razão molar nominal unitária, a correção residual reproduziu o patamar abaixo da conversão completa: em 15 min, a conversão prevista foi de 0,852, frente a 0,870 medida. O white-box, por sua vez, aproxima-se da conversão completa nessa condição. No Ensaio 14, a previsão em 15 min foi de 0,972, próxima da medição de 0,970. No Ensaio 7, o modelo manteve uma superestimação nos primeiros instantes, que concentrou o maior erro de teste.

Figura 1 – Conversão de zinco medida e prevista pelo modelo híbrido paralelo com Processo Gaussiano e pelo white-box comum nos Ensaios 1, 7, 8 e 14; o asterisco identifica os ensaios de teste.

Fonte: Elaborado pelos autores.

Como a avaliação de teste reúne apenas três ensaios, o desempenho foi examinado também por validação cruzada agrupada por ensaio, com as quatro dobras comuns aos modelos serial e black-box e apenas os 13 ensaios de treinamento. O RMSE fora da dobra, calculado nos 91 pontos experimentais e depois promediado entre as dez sementes, foi de 1,97 p.p. As médias por dobra variaram de 1,17 a 2,49 p.p. O maior erro no teste, em relação à validação cruzada, reforça a necessidade de distinguir o desempenho entre condições operacionais.

Antes da aplicação do acoplador, 0,56% dos nós das curvas de teste apresentaram inconsistência física segundo o critério implementado. Após a imposição da condição inicial, do limite estequiométrico e da monotonicidade, essa fração foi nula. As restrições foram aplicadas em uma malha comum de 0,25 min, da qual foram extraídas as previsões nos instantes de medição. Assim, a admissibilidade física das curvas finais decorre da ação do acoplador, enquanto a redução do erro em relação ao white-box expressa o resultado da correção residual em conjunto com esse acoplamento.

As curvas dos 16 ensaios são apresentadas na Figura 2, e a comparação entre conversões medidas e previstas, na Figura 3. O maior afastamento da identidade entre os pontos de teste está associado à etapa inicial do Ensaio 7. O bom ajuste global deve ser considerado juntamente com essa limitação localizada e com o número reduzido de condições independentes de teste.

Figura 2 – Conversão de zinco medida e prevista pelo modelo híbrido paralelo com Processo Gaussiano e pelo white-box comum nos 16 ensaios originais; o asterisco identifica os ensaios de teste.

Figura 3 – Comparação entre a conversão medida e a prevista pelo modelo híbrido paralelo com Processo Gaussiano nos 112 pontos experimentais com t > 0, distinguindo os conjuntos de treino e teste.

Fonte das Figuras 2 e 3: Elaborado pelos autores.

## Ensaios novos de laboratório

O modelo foi confrontado com os quatro ensaios novos E1 a E4, com volume de solução de 200 mL e amostragem em 0,5; 1; 2; 3; 5 e 15 min. Foram utilizadas as mesmas configurações e sementes dos modelos treinados nos 13 ensaios originais, sem ajuste aos ensaios novos. As conversões foram calculadas com teor de 0,611 g de Zn associado à zincita por grama de sólido. As massas dos ensaios, por sua vez, haviam sido definidas com teor de zinco total de 0,650 g/g. Para manter a mesma base na avaliação, a razão molar foi recalculada pelas massas efetivamente alimentadas, usando η = V CA0 M_Zn / (m × 0,611). As razões planejadas foram conservadas apenas como identificação das condições de preparo.

A identificação de E2 e E3 segue a correção já registrada na planilha. A conversão das leituras conserva a hipótese de unidade em g/L e não inclui as correções por ferro e traçador empregadas na Equação 4.5 da dissertação original. Sua interpretação depende da confirmação dessas hipóteses e da representatividade da composição adotada para o material dos novos ensaios.

Tabela 2 – Condições dos ensaios novos na base de zinco associado à zincita (0,611 g/g), conversões em 15 min e erros nos seis instantes de amostragem.

| Ensaio | CA0 (mol/L) | η planejada na base Zn total | η efetiva na base zincita | Limite de X | Medida em 15 min | Prevista em 15 min | RMSE GP (p.p.) | Erro máximo GP (p.p.) | RMSE white-box (p.p.) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| E1 | 0,25 | 0,75 | 0,799 | 0,799 | 0,809 | 0,770 | 5,97 | 9,74 | 4,08 |
| E2 | 0,75 | 1,25 | 1,327 | 1,000 | 0,936 | 0,910 | 2,95 | 4,11 | 11,86 |
| E3 | 0,75 | 2,25 | 2,396 | 1,000 | 0,971 | 0,986 | 5,36 | 9,43 | 7,19 |
| E4 | 0,25 | 0,40 | 0,425 | 0,425 | 0,447 | 0,370 | 6,81 | 7,79 | 1,43 |

Fonte: Elaborado pelos autores.

A Figura 4 apresenta as curvas correspondentes, com as razões molares efetivas. Em E2 e E3, o RMSE do híbrido foi inferior ao do white-box, passando de 11,86 para 2,95 p.p. e de 7,19 para 5,36 p.p., respectivamente. Em E3, o maior erro ocorreu no primeiro instante de amostragem, enquanto a previsão em 15 min diferiu da medição em aproximadamente 1,53 p.p.

A correção da razão molar elevou os limites de E1 e E4 para 0,799 e 0,425. Em ambos, a medição em 0,5 min ficou abaixo do limite calculado, enquanto os cinco pontos seguintes o ultrapassaram. Portanto, a mistura entre a base de zinco total e a de zincita explicava parte da aparente incompatibilidade, mas não toda. Em E1, a conversão medida também diminuiu após 1 min. Essas diferenças exigem a conferência das condições efetivas, da composição do sólido e do tratamento das leituras, antes de atribuir integralmente os desvios à modelagem.

O GP apresentou erro superior ao do white-box em E1 e E4. Em E4, previu um patamar de aproximadamente 0,370, abaixo do limite de 0,425. Essa condição permanece abaixo do menor nível nominal de treinamento (0,50) e evidencia uma limitação da correção residual em extrapolação. O acoplador assegura a admissibilidade física da curva, mas não garante sua acurácia. As diferenças de desempenho não devem ser interpretadas apenas como efeito do volume de 200 mL, pois a disponibilidade de ácido e o tratamento experimental também diferem.

Figura 4 – Conversão de zinco medida e prevista pelo modelo híbrido paralelo com Processo Gaussiano e pelo white-box nos ensaios novos E1 a E4, com conversões e razões molares na base de 0,611 g de Zn associado à zincita por grama de sólido. As razões molares são calculadas pelas massas efetivamente alimentadas; E2 e E3 seguem a identificação corrigida da planilha.

Fonte: Elaborado pelos autores.

# Sugestão de conclusão

A modelagem híbrida paralela com Processo Gaussiano permitiu corrigir os desvios do modelo fenomenológico nas condições avaliadas, reduzindo o RMSE dos três ensaios de teste de 7,80 para 2,66 p.p. A correção residual reproduziu o patamar de conversão inferior à unidade próximo da condição estequiométrica, enquanto o acoplador físico assegurou o atendimento à condição inicial, à monotonicidade e ao limite de disponibilidade de ácido. Persistiu, entretanto, uma superestimação da conversão no início do Ensaio 7. Nos ensaios novos, o híbrido melhorou a previsão de E2 e E3 em relação ao white-box, mas apresentou pior desempenho em E1 e E4, sendo E4 uma condição de extrapolação em razão molar. A interpretação de E1 e E4 também depende da conferência das conversões experimentais e da composição do sólido. A análise utiliza conversões e razões molares na mesma base de 0,611 g/g; E1 e E4 ainda apresentam cinco de seis medições acima dos limites calculados. Assim, o bom desempenho na partição original não assegura acurácia em todas as condições novas.
