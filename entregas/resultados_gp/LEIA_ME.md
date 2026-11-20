# Resultados do GP com o white-box comum

Esta pasta foi recalculada após a alteração do white-box. Os arquivos antigos em resultados/ não foram substituídos.

- Script reproduzível: scripts/gerar_resultados_gp_relatorio.py.
- GP adotado pelo autor; não selecionado pelo menor erro de teste.
- Núcleo Matérn 3/2, ARD, ruído branco, normalização do alvo, cinco reinícios da verossimilhança marginal. Pré-processamento ajustado somente no conjunto de treino de cada dobra.
- Modelo final treinado nos 13 ensaios; dez sementes 42–51; média das saídas brutas antes da projeção.
- White-box: η nominal, ks = 18000 µm/min sem Arrhenius, α = 0, ν = 1, força motriz CAf.
- Atributos do GP: tempo, η nominal, CA0 e X_PBM (conjunto B), na mesma base de η do white-box e dos demais modelos.
- Teste: ensaios 7, 8 e 14. CV: quatro dobras da partição comum. O histórico anterior já incluiu avaliações desses ensaios; não se deve alegar uma primeira consulta inédita nesta execução.
- Projeção isotônica realizada na mesma malha de 0,25 min na CV e na avaliação final. Métricas extraídas dessa curva nos instantes medidos, t > 0. Isso corrige a diferença entre projetar separadamente os pontos e a curva na rotina antiga.
- Inconsistência: porcentagem dos nós t > 0 com queda, conversão negativa ou acima do teto; tolerância de 0,001 em conversão (0,1 p.p.). Valores finais são também verificados contra os limites, sem essa tolerância.
- As figuras não incluem faixas de incerteza ou barras de erro experimental. A comparação usa as médias medidas disponíveis.
- Este pacote inclui os 16 ensaios originais e a avaliação dos quatro ensaios novos E1–E4. O script scripts/gerar_ensaios_novos_gp.py reproduz essa extensão, usando o mesmo ajuste aos 13 ensaios e verificando a reprodução das previsões originais.

Para o corpo do relatório, usar a Figura 1 e a Tabela 1. A Figura 2 pode ir ao apêndice, e a Figura 3 é complementar. PNG em 300 dpi e PDF vetorial são fornecidos.

Estilo atualizado conforme a referência visual fornecida: Times New Roman, white-box azul pontilhado, GP púrpura, pontos experimentais vazados, teto cinza pontilhado, valores dos eixos visíveis em todos os painéis.

A base adotada para E1–E4 é 0,611 g Zn/g de sólido. O teor de 0,650 serviu ao preparo/pesagem, e não à conversão usada na comparação. A conversão continua na base de 0,611 g Zn/g. A razão molar usada na entrada do GP, no white-box e no teto do acoplador é a nominal do planejamento (0,75; 1,25; 2,25; 0,40), na mesma base dos 16 ensaios originais e igual à usada pelo white-box comum do grupo. A coluna razao_molar_eta contém esse valor, eta_ZnO guarda a razão na base da zincita, V CA0 × 65,38 / (massa × 0,611), e eta_planejada_Zn_total conserva o valor planejado. Os ensaios originais permanecem na partição e configuração comum já avaliadas. A identificação E2/E3 já foi corrigida na planilha e não foi trocada novamente. As unidades das leituras e a contribuição de outras fases requerem confirmação; não foram aplicadas as correções por ferro e traçador da dissertação original.

A versão `figura_4_ensaios_novos_gp_0611_corrigido` é a atual. A versão antiga `_planilha` permanece apenas como arquivo histórico local e foi excluída do pacote para entrega.
