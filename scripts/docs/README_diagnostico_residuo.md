# README -- Exportador do resíduo (`diagnostico_residuo.py`)

## 1. Papel
Exporta a tabela do **resíduo do White-Box**, $\Delta X = X_{exp} - X_{PBM}$, para a configuração escolhida (padrão: o White-Box
principal). É só uma ferramenta de inspeção: o pipeline **não lê** o arquivo e sempre recalcula o PBM (não há cache).

```powershell
& .venv\Scripts\python.exe scripts/diagnostico_residuo.py              # White-Box principal -> resultados/dados_com_residuo.csv
& .venv\Scripts\python.exe scripts/diagnostico_residuo.py --legado     # White-Box anterior à revisão (V0, nu = 1)
```
```python
from diagnostico_residuo import calcular_residuos
df = calcular_residuos()          # DataFrame com X_pbm, C_acid_pbm_mol_L, H_pbm_mol_L, delta_X
```

## 2. O que o resíduo contém (fatos medidos; as causas são hipóteses)
Resíduo aos 15 min do White-Box principal (Tabela T14 de `resultados/tabelas.md`):

* **η = 0,5:** ≈ 0 (−0,004 a −0,014): o ν corrige o patamar; a versão antiga previa 0,535 contra 0,50 medido (resíduo −0,035).
* **η = 1,0:** **−0,10 a −0,14**. O experimento satura em 0,85–0,87 com ácido e zincita em proporção estequiométrica; o PBM prevê 0,97–0,99.
* **η = 1,5:** −0,02 a −0,03; **η = 3,1:** ≈ 0.
* O desvio concentra-se em η ≈ 1 e cresce com $C_{A0}$ (−0,10 a 0,1 mol/L; −0,14 a 1,5 mol/L).

Nesse patamar a cinética parece **não** ser o limitante: nos ensaios 8 a 10 o pH final medido (0,6 a 1,2) indica [H⁺] de 0,06 a
0,25 mol/L e, com o $k_s$ de Balarini ($v = 617\,[\text{H}^+]$ µm/min por mol/L), qualquer partícula da distribuição se dissolveria em
poucos minutos (é o que o PBM prevê: X = 0,99). Hipóteses a testar (nenhuma foi testada neste trabalho): fração de zinco ocluída ou refratária
(ferrita, silicatos), passivação por CaSO₄/PbSO₄ (a calcina tem 1,4% de Ca e 2,3% de Pb), efeito de íon comum e atividade iônica,
erro de base de cálculo da Eq. 4.5. O ML aprende esse desvio **como função empírica de (t, η, $C_{A0}$)**, só onde há ensaios.

## 3. Estrutura
`calcular_residuos(config=CONFIG_WB_PRINCIPAL, caminho_saida=None)` usa `hybrid_model.construir_base`: simula o PBM rápido em cada
ensaio com `config` e devolve o quadro com `X_pbm`, `C_acid_pbm_mol_L`, `H_pbm_mol_L` e `delta_X`.
