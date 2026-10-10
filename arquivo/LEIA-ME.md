# Arquivo (material superado, mantido só para consulta)

Nada aqui é usado pelo pipeline atual (`scripts/experimentos_batelada.py` e `scripts/figuras_batelada.py`).
Foi movido, não apagado, na revisão de validade do modelo híbrido (batelada). Motivos:

| Item | Por que saiu |
| :--- | :--- |
| `scripts_antigos/plotar_cinetica.py`, `plotar_cinetica_teste_cego.py`, `gerar_graficos_logo.py`, `gerar_graficos_todos_modelos.py` | Quatro laços de gráfico quase idênticos, cada um com sua própria forma de validação cruzada e de calcular métricas, o que gerava números diferentes entre figuras. Substituídos por `figuras_batelada.py`, que só lê os CSV do pipeline. Importam funções que não existem mais (`plotar_resultados`, o acoplador estático). |
| `scripts_antigos/testar_arquiteturas_mlp.py`, `plotar_comparativo_arquiteturas_mlp.py`, `plotar_mlp_original_vs_calibrada.py` | Escolha de arquitetura da MLP com o mesmo CV usado para reportar (otimista) e com o nome de modelo `'Rede Neural (MLP Tanh + L-BFGS)'`, que já não existe no catálogo (ValueError). A MLP agora é ajustada por CV aninhada dentro do pipeline. |
| `teste_augmentation/` | O aumento de dados (PCHIP + PBM) simulava o PBM com o η **nominal** (o resíduo-alvo do treino ficava em −0,003 em η = 0,5, contra −0,037 no teste), usava uma versão antiga do White-Box e interpola o próprio dado: não traz informação sobre condições novas. Os números do README (por exemplo SVR com CV 0,9865 contra 0,9446) não são reproduzíveis com o código atual. |
| `figuras_antigas/` | Pastas de gráficos duplicadas (com e sem acento, cópias byte a byte) e as figuras de `scripts/docs/`: restam 59 PNGs distintas. Duas das oito pastas originais (`graficos2` e `gráficos2` da raiz de `lop/`) não existiam mais quando a arrumação chegou a elas; o conteúdo era idêntico ao de `scripts/graficos2` e `scripts/gráficos2`, que estão aqui, e nenhuma figura única se perdeu. As figuras atuais são geradas em `resultados/figuras/`. |
| `previsoes_pontos_teste.csv` | Previsões de GP, ET e ML puro em condições não ensaiadas (η = 0,75), geradas por um script não versionado e com o White-Box antigo. Substituídas pela análise de domínio (`resultados/dominio_incerteza_eta.csv`, figura f10). |
| `dados_com_residuo_legado.csv` | Cache do resíduo do White-Box **anterior** à revisão (K₂ = 0,0099, ν = 1). O pipeline não lê mais cache: o PBM é recalculado a cada execução. Uma cópia reduzida serve de referência de regressão em `tests/dados/`. |

Para reproduzir os números "antes" da revisão use `resultados/baseline_antes.csv` e `tests/test_regressao.py`.
