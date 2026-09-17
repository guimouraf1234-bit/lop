# README -- Módulo de Balanço Populacional (`population_balance.py`)

## 1. Visão Geral e Objetivo
O módulo `population_balance.py` implementa a solução analítico-numérica do **Balanço Populacional Unidimensional (PBM)** para a lixiviação de concentrado ustulado de zinco em reator batelada.

Ele é o componente central da modelagem mecanicista (\textit{White-Box}), responsável por integrar o encolhimento de todas as frações de tamanho da distribuição granulométrica ao longo do tempo e calcular a **conversão global de zinco ($X_{\text{PBM}}$)** e o **perfil de consumo de ácido ($C_{Af}$)**.

---

## 2. Fundamentação Teórica e Equacionamento

### 2.1. A Equação do Balanço Populacional (PBM)
Para um reator batelada com esfericidade constante e sem fenômenos de quebra ou aglomeração de partículas ($B = D = 0$), a evolução da densidade de distribuição de tamanhos $\psi(D,t)$ é descrita por:
$$\frac{\partial \psi(D,t)}{\partial t} + \frac{\partial \big[ -v(t) \, \psi(D,t) \big]}{\partial D} = 0$$

### 2.2. Solução pelo Método das Características
Como a velocidade de retração $v(t)$ é puramente química e independe do tamanho $D$, todas as partículas encolhem à mesma taxa em um dado instante:
$$\frac{dD}{dt} = -v(t) \implies D(t) = \max \big( D_0 - \Delta D(t),\ 0 \big)$$

onde $\Delta D(t)$ representa o encolhimento linear acumulado do diâmetro:
$$\Delta D(t) = \int_0^t v(t') \, dt' \quad [\mu\text{m}]$$

### 2.3. Distribuição Granulométrica Inicial (Rosin-Rammler-Bennett -- RRB)
A calcina de zinco de Três Marias segue a distribuição RRB:
$$F_0(D) = 1 - \exp \left[ -\left( \frac{D}{D_{63{,}2}} \right)^m \right]$$
$$f_0(D) = \frac{m}{D_{63{,}2}} \left( \frac{D}{D_{63{,}2}} \right)^{m-1} \exp \left[ -\left( \frac{D}{D_{63{,}2}} \right)^m \right]$$
com $m = 1{,}022$ e $D_{63{,}2} = 41{,}65\ \mu\text{m}$ (Bortot Coelho, 2017).

### 2.4. Conversão pelo 3º Momento Volumétrico
A conversão fracionária da zincita, $X_{\text{PBM}}(t)$, é a fração de volume de sólido solubilizado:
$$X_{\text{PBM}}(t) = 1 - \frac{\mu_3(t)}{\mu_3(0)} = 1 - \frac{\displaystyle \int_{\Delta D(t)}^\infty \big( D - \Delta D(t) \big)^3 f_0(D) \, dD}{\displaystyle \int_0^\infty D^3 f_0(D) \, dD}$$

* **Solução Analítica Exata de $I_0$ (3º momento inicial):**
  $$I_0 \equiv \int_0^\infty D^3 f_0(D) \, dD = D_{63{,}2}^3 \cdot \Gamma \left( 1 + \frac{3}{m} \right)$$
  Para $m = 1{,}022$, $\Gamma(1 + 3/1{,}022) \approx 5{,}6645 \implies I_0 \approx 4{,}094 \times 10^5\ \mu\text{m}^3$.

### 2.5. Sistema Dinâmico Acoplado
O acoplamento mútuo entre a retração da partícula e o esgotamento do ácido resulta no sistema de EDOs:
$$\begin{cases}
\displaystyle \frac{d(\Delta D)}{dt} = v \big( C_{Af}(t) \big) \\[0.8em]
\displaystyle X_{\text{PBM}}(t) = 1 - \frac{1}{I_0} \int_{\Delta D(t)}^{D_{\text{máx}}} \big( D - \Delta D(t) \big)^3 f_0(D) \, dD \\[0.8em]
\displaystyle C_{Af}(t) = C_{A0} \left( 1 - \frac{X_{\text{PBM}}(t)}{\eta} \right)
\end{cases}$$
resolvido numericamente pelo algoritmo de Runge-Kutta de passo adaptativo (RK45 - `scipy.integrate.solve_ivp`).

---

## 3. Estrutura da Classe `PopulationBalanceModel`

### `__init__(m=1.022, d632=41.65, f_zno=0.87)`
Inicializa os parâmetros granulométricos e calcula analiticamente o 3º momento inicial $I_0$ via função Gamma.

### `rrb_pdf(D: float) -> float`
Avalia a densidade de probabilidade de tamanho de partícula $f_0(D)$ em $\mu\text{m}^{-1}$.

### `calcular_conversao_zincita(delta_D: float) -> float`
Calcula a conversão fracionária da zincita a partir do encolhimento linear acumulado $\Delta D$ via quadratura de Gauss-Kronrod (`scipy.integrate.quad`).

### `simular_ensaio(tempos_min, C_acid_0_mol_L, eta, T_celsius=40.0, alpha=0.0) -> dict`
Integra o sistema acoplado para as condições operacionais informadas. Retorna um dicionário com:
* `tempo_min`: Vetor de tempos avaliados.
* `delta_D_um`: Retração acumulada $\Delta D(t)$ em $\mu\text{m}$.
* `X_pbm`: Conversão da zincita ($0$ a $1$).
* `C_acid_mol_L`: Concentração de ácido livre residual $C_{Af}(t)$ em $\text{mol/L}$.

---

## 4. Exemplo de Uso

```python
from population_balance import PopulationBalanceModel

pbm = PopulationBalanceModel()

# Simular Ensaio 6 (CA0 = 0.1 mol/L, eta = 1.0)
tempos = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 15.0]
res = pbm.simular_ensaio(tempos_min=tempos, C_acid_0_mol_L=0.10, eta=1.0)

for t, x, c in zip(res['tempo_min'], res['X_pbm'], res['C_acid_mol_L']):
    print(f"t = {t:4.1f} min | X = {x*100:5.1f}% | [H2SO4] = {c:6.4f} mol/L")
```

---

## 5. Como Testar
Execute no terminal:
```powershell
& .venv\Scripts\python.exe scripts/population_balance.py
```
*(Ou através de um script de teste)*.
