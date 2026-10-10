# README -- Módulo de Balanço Populacional (`population_balance.py`)

## 1. Visão Geral e Objetivo
O módulo `population_balance.py` implementa o **Balanço Populacional Unidimensional (PBM)** da lixiviação de concentrado ustulado de zinco em reator batelada, a partir do modelo de Bortot Coelho (2017) (item 4.3, Tabela 5.9 e Apêndice A2.1 da dissertação).

Ele é o componente central da modelagem mecanicista (\textit{White-Box}): integra o encolhimento das partículas ao longo do tempo e calcula a **conversão da zincita ($X_{\text{PBM}}$)**, o **ácido livre ($C_{Af}$)** e o **[H⁺]** do licor.

Duas formulações estão disponíveis (argumento `formulacao`):

| | `'corrigida'` (padrão) | `'tese'` (Apêndice A2.1) |
| :--- | :--- | :--- |
| RRB | densidade mássica $q_3$ | densidade numérica com peso $D^3$ (Eq. 3.49) |
| $k_s$ | 18000 µm/min a 30 °C + Arrhenius (21344 a 40 °C) | 18000 µm/min direto |
| Força motriz | $[\text{H}^+]$ por especiação | $C_{Af}$ (molaridade de H₂SO₄) |
| $R^2$ ($t > 0$), $\alpha = 0$ | 0,9486 | 0,7475 |

Em ambas, $\alpha = 0$ por padrão; `alpha=5500` reproduz a Eq. 4.1 da tese.

---

## 2. Fundamentação Teórica e Equacionamento

### 2.1. Equação do Balanço Populacional (Eq. 3.38)
Sem quebra nem aglomeração ($B = D = 0$):
$$\frac{\partial \psi(D,t)}{\partial t} = -\frac{\partial}{\partial D} \big[ -v \, \psi(D,t) \big]$$

### 2.2. Solução pelo Método das Características
Como $v$ não depende de $D$, a distribuição apenas se desloca:
$$\psi(D,t) = \psi\big(D + \Delta D(t),\ 0\big), \qquad \Delta D(t) = \int_0^t v(t') \, dt'$$

### 2.3. Distribuição Granulométrica Inicial (RRB, Tabela 5.9)
$$q_3(D) = \frac{m}{D_{63{,}2}} \left( \frac{D}{D_{63{,}2}} \right)^{m-1} \exp \left[ -\left( \frac{D}{D_{63{,}2}} \right)^m \right]$$
com $m = 1{,}022$, $D_{63{,}2} = 41{,}65\ \mu\text{m}$, $D_{\text{mín}} = 0{,}1\ \mu\text{m}$ e $D_{\text{máx}} = 297\ \mu\text{m}$. A RRB foi ajustada a frações mássicas (Tabela 5.4 da tese).

### 2.4. Conversão da Zincita
Formulação corrigida (cada classe retém a fração $(1 - \Delta D/D)^3$ da sua massa):
$$X(t) = 1 - \frac{\int_{\max(\Delta D,\,D_{\text{mín}})}^{D_{\text{máx}}} (1 - \Delta D/D)^3 \, q_3(D) \, dD}{\int_{D_{\text{mín}}}^{D_{\text{máx}}} q_3(D) \, dD}$$
Formulação da tese (Eq. 3.49): $X = 1 - \int (D - \Delta D)^3 f\,dD \,/ \int D^3 f\,dD$.

### 2.5. Razão Molar (Eq. 4.3)
$$\eta = \frac{V \, C_{A0}}{m_{B0} \, T_{\text{ZnO}} / M_{\text{ZnO}}}, \qquad T_{\text{ZnO}} = 0{,}761,\ M_{\text{ZnO}} = 81{,}38\ \text{g/mol}$$
calculada por `calcular_eta` com as massas de cada ensaio (coluna `razao_molar_eta_real` gerada pelo `data_loader`).

### 2.6. Sistema Dinâmico Acoplado (formulação corrigida)
$$\begin{cases}
\displaystyle \frac{d(\Delta D)}{dt} = \frac{2 \, k_s(T)}{\rho_s} \, [\text{H}^+] \\[0.8em]
\displaystyle C_{Af}(t) = C_{A0} \left( 1 - \frac{X(t)}{\eta} \right) \\[0.8em]
\displaystyle [\text{H}^+] = f_{\text{especiação}}(C_{A0},\ C_{Af};\ K_2 = 0{,}0099)
\end{cases}$$
resolvido por Runge-Kutta de passo adaptativo (RK45, `scipy.integrate.solve_ivp`).

---

## 3. Estrutura do Módulo

### `calcular_eta(volume_L, C_acid_0_mol_L, massa_solido_g) -> float`
Razão molar H₂SO₄ / ZnO pela Eq. 4.3.

### `PopulationBalanceModel(m=1.022, d632=41.65, d_min=0.1, d_max=297.0, formulacao='corrigida')`
Inicializa a RRB e a integral de normalização (`I0`).

### `rrb_pdf(D: float) -> float`
Avalia a densidade RRB em $\mu\text{m}^{-1}$.

### `calcular_conversao_zincita(delta_D: float) -> float`
Conversão da zincita para o encolhimento acumulado $\Delta D$, por quadratura de Gauss-Kronrod (`scipy.integrate.quad`).

### `simular_ensaio(tempos_min, C_acid_0_mol_L, eta, T_celsius=40.0, alpha=0.0) -> dict`
Integra o sistema acoplado. Retorna:
* `tempo_min`: tempos avaliados.
* `delta_D_um`: encolhimento acumulado $\Delta D(t)$ em $\mu\text{m}$.
* `X_pbm`: conversão da zincita (0 a 1), comparável a $X_{\text{Zn}}^{\text{exp}}$ (Eq. 4.5).
* `C_acid_mol_L`: ácido livre $C_{Af}(t)$ em mol/L.
* `C_H_mol_L`: $[\text{H}^+](t)$ em mol/L (NaN na formulação da tese).

---

## 4. Exemplo de Uso

```python
from population_balance import PopulationBalanceModel, calcular_eta

pbm = PopulationBalanceModel()                       # formulação corrigida
eta = calcular_eta(volume_L=0.4, C_acid_0_mol_L=0.10, massa_solido_g=4.0)  # 1,069 (ensaio 6)

tempos = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0]
res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=0.10, eta=eta)

for t, x, c, h in zip(res['tempo_min'], res['X_pbm'], res['C_acid_mol_L'], res['C_H_mol_L']):
    print(f"t = {t:4.1f} min | X = {x*100:5.1f}% | [H2SO4] = {c:6.4f} | [H+] = {h:6.4f} mol/L")
```

Saída: a conversão sobe para $66{,}6\%$ em 0,5 min, $90{,}0\%$ em 5 min e $96{,}1\%$ em 15 min, com $0{,}0101$ mol/L de ácido livre e $[\text{H}^+] = 0{,}0022$ mol/L ao final. O experimento mede $75\%$, $85\%$ e $87\%$.

---

## 5. Como Testar
```powershell
& .venv\Scripts\python.exe scripts/diagnostico_residuo.py
```
Recalcula `X_pbm`, `C_acid_pbm_mol_L` e o resíduo `delta_X` dos 128 pontos e grava `data/processed/dados_com_residuo.csv`.
