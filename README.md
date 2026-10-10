# LOP -- Modelagem híbrida da lixiviação de calcina de zinco (batelada)

TCC de Engenharia Química (UFMG, 2026). White-Box (balanço populacional) + Machine Learning residual, validado por ensaio inteiro,
com baselines (tabela, Ridge, caixa-preta), protocolos de extrapolação e domínio de aplicabilidade declarado.

**Leia primeiro:** [`docs/validade_modelo_hibrido.md`](docs/validade_modelo_hibrido.md) (o que foi feito, os resultados e as limitações).

## Estrutura
```
scripts/
  kinetics.py              ConfigWB, cinética, K2(T), par iônico ZnSO4, nu                (README_kinetics.md)
  population_balance.py    PBM (integrador de referência e rápido)                         (README_population_balance.md)
  data_loader.py           128 pontos, respostas finais, parâmetros                         (README_data_loader.md)
  hybrid_model.py          atributos A/B/C, acoplador físico, catálogo, ModeloHibridoCinza (README_hybrid_model.md)
  validacao.py             protocolos, baselines, métricas, bootstrap, intervalos           (README_validacao.md)
  experimentos_batelada.py PIPELINE ÚNICO -> resultados/*.csv                              (README_experimentos_batelada.md)
  figuras_batelada.py      lê os CSV -> resultados/figuras/*.png
  tabelas_batelada.py      lê os CSV -> resultados/tabelas.md
  diagnostico_residuo.py   exporta ΔX = X_exp − X_PBM                                       (README_diagnostico_residuo.md)
  continuous_reactor_model.py   planta piloto (FORA do escopo da revisão; roda em modo legado)
tests/                     56 testes (invariantes físicos, vazamento, regressão contra os números antigos, pipeline)
resultados/                saídas do pipeline (CSV versionados; figuras regeneráveis)
docs/                      relatório de validade, white-box, arquitetura
arquivo/                   material superado (ver arquivo/LEIA-ME.md)
```

## Reproduzir tudo
```powershell
& ..\.venv\Scripts\python.exe -m pytest                              # 56 testes, ~20 s
& ..\.venv\Scripts\python.exe scripts/experimentos_batelada.py       # protocolo completo, ~11 min
& ..\.venv\Scripts\python.exe scripts/figuras_batelada.py
& ..\.venv\Scripts\python.exe scripts/tabelas_batelada.py
```
Dependências: `requirements.txt` (numpy, scipy, pandas, scikit-learn, joblib, matplotlib, pypdf, pytest).

## Decisões de projeto que mudam o texto do TCC
1. **White-Box principal = V2** (K₂(T) por van't Hoff + par iônico ZnSO₄⁰ + ν a priori), escolhido pelo **pH medido** (única medida independente do estado ácido), não pelo ajuste de X.
2. **Atributos do ML = B** (t, η real, $C_{A0}$, $X_{PBM}$).
3. **Validade:** o híbrido é avaliado contra White-Box, tabela e caixa-preta em seis protocolos; o domínio são os **quatro níveis de η ensaiados** (entre eles não há informação).
4. **Fora do escopo:** planta piloto (volume 18 L × 6 L no JSON, fator 0,05 arbitrário no código); híbrido serial (ver branch `feat/modelagem-alpha`).
