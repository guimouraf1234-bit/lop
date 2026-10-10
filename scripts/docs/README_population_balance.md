# README -- Balanço Populacional (`population_balance.py`)

## 1. Papel
`population_balance.py` implementa o **Balanço Populacional unidimensional (PBM)** da lixiviação de calcina de zinco em
batelada (Bortot Coelho, 2017: item 4.3, Tabela 5.9 e Apêndice A2.1). Integra o encolhimento das partículas e entrega a
**conversão da zincita** $X_{PBM}$, o **ácido livre** $C_{Af}$, o **[H⁺]** e o pH do licor. É a base física (White-Box) do híbrido.

A configuração física é um `kinetics.ConfigWB` (ver `README_kinetics.md`). Sem `config`, `simular_ensaio` reproduz exatamente
o comportamento anterior à revisão (`formulacao='corrigida'` ou `'tese'`).

| | `ConfigWB.legado()` (TCC antes da revisão) | **`CONFIG_WB_PRINCIPAL`** (V2, ν = 1,06) | `ConfigWB.tese(α)` |
| :--- | :--- | :--- | :--- |
| RRB | base mássica | base mássica | numérica, peso $D^3$ (Eq. 3.49) |
| $k_s$ | 18000 µm/min a 30 °C + Arrhenius | idem | 18000 µm/min fixo |
| Força motriz | [H⁺], $K_2 = 0{,}0099$ | [H⁺], $K_2(T)$ + par iônico ZnSO₄⁰ | $C_{Af}$ |
| $C_{Af}$ | $C_{A0}(1 - X/\eta)$ | $C_{A0}(1 - \nu X/\eta)$ | $C_{A0}(1 - X/\eta)$ |
| RMSE / R² de X (t > 0) | 0,0443 / 0,9486 | 0,0620 / 0,8992 | 0,0981 / 0,7475 (α = 0) |

O White-Box principal ajusta **pior** a conversão que o legado: ele não "esconde" com o $K_2$ constante o patamar que
a física não explica, e em troca é o único coerente com o pH medido (ver `docs/validade_modelo_hibrido.md`).

## 2. Equações (resumo)
* Balanço (Eq. 3.38) com $B = D = 0$ e $v$ independente de $D$: solução por características, $\psi(D,t) = \psi(D + \Delta D(t), 0)$, $\Delta D = \int_0^t v\,dt'$.
* RRB (Tabela 5.9, ajustada a frações **mássicas**): $q_3(D) = \tfrac{m}{D_{63{,}2}}(D/D_{63{,}2})^{m-1}\exp[-(D/D_{63{,}2})^m]$, $m = 1{,}022$, $D_{63{,}2} = 41{,}65\ \mu$m, $D \in [0{,}1;\ 297]\ \mu$m.
* Conversão (base mássica): $X = 1 - \dfrac{\int_{\max(\Delta D, D_{min})}^{D_{max}}(1-\Delta D/D)^3\,q_3\,dD}{\int_{D_{min}}^{D_{max}} q_3\,dD}$.
* $\eta = V C_{A0} / (m_{B0} T_{ZnO}/M_{ZnO})$ (Eq. 4.3; $T_{ZnO} = 0{,}761$, $M_{ZnO} = 81{,}38$ g/mol), calculado com as massas de cada ensaio (`razao_molar_eta_real`).
* EDO acoplada: $\dfrac{d\Delta D}{dt} = \dfrac{2k_s(T)}{\rho_s}\,C(t)$, com $C$ = [H⁺] (ou $C_{Af}$ na tese) obtida do estado químico para a conversão atual.

## 3. Dois integradores equivalentes
| Método | Como | Custo (16 ensaios) |
| :--- | :--- | :--- |
| `simular_ensaio(...)` | **referência**: recalcula a quadratura a cada passo da EDO (RK45) | ~1–3 s |
| `simular_ensaio_rapido(...)` | usa a tabela $X(\Delta D)$ (3501 pontos, calculada uma vez por instância) e LSODA | ~0,5 s |

A diferença máxima em X é < 1e-4 (medida: ≤ 5e-6; `tests/test_pbm.py`). O pipeline usa o rápido; a referência reproduz o CSV
legado a 1e-9. Tempos fora de ordem são aceitos (a saída segue a ordem pedida).

## 4. API
```python
PopulationBalanceModel(m=1.022, d632=41.65, d_min=0.1, d_max=297.0, formulacao='corrigida')   # 'corrigida' | 'tese'
PopulationBalanceModel.para_config(cfg)          # instância compartilhada com a base de RRB de cfg
pbm.simular_ensaio(tempos_min, C_acid_0_mol_L, eta, T_celsius=40.0, alpha=0.0, config=None)
pbm.simular_ensaio_rapido(...)                   # mesma interface
calcular_eta(volume_L, C_acid_0_mol_L, massa_solido_g)
```
Saída (dict): `tempo_min`, `delta_D_um`, `X_pbm`, `C_acid_mol_L` (com ν), `C_H_mol_L` e `pH_pbm` (NaN na formulação da tese).
Uma `config` com `rrb` incompatível com a instância levanta `ValueError`.

## 5. Exemplo (ensaio 6: $C_{A0} = 0{,}10$ mol/L, 4,0 g, η = 1,069)
```python
from kinetics import ConfigWB, NU_A_PRIORI
from population_balance import PopulationBalanceModel, calcular_eta

eta = calcular_eta(volume_L=0.4, C_acid_0_mol_L=0.10, massa_solido_g=4.0)            # 1,0694
tempos = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0]
for nome, cfg in (('legado', ConfigWB.legado()), ('principal', ConfigWB(especiacao='V2', nu=NU_A_PRIORI))):
    r = PopulationBalanceModel.para_config(cfg).simular_ensaio_rapido(tempos, 0.10, eta, config=cfg)
    print(nome, r['X_pbm'].round(3))
```
```text
legado    [0.    0.666 0.75  0.823 0.86  0.883 0.9   0.961]
principal [0.    0.719 0.805 0.872 0.902 0.92  0.932 0.972]
experimento: 0, 0.75, 0.78, 0.80, 0.85, 0.87, 0.85, 0.87
```
Aos 15 min o experimento mede 0,87; o White-Box principal prevê 0,972 (resíduo −0,10). Nos quatro ensaios com η = 1 os resíduos aos 15 min vão de −0,10 a −0,14 (Tabela T14 de `resultados/tabelas.md`).

## 6. Testes e reprodução
```powershell
& .venv\Scripts\python.exe -m pytest tests/test_pbm.py -q                  # invariantes e equivalência dos integradores
& .venv\Scripts\python.exe scripts/diagnostico_residuo.py                  # exporta o resíduo ΔX em resultados/dados_com_residuo.csv
```
