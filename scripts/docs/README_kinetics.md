# README -- Módulo de Cinética Química Heterogênea (`kinetics.py`)

## 1. Visão Geral e Papel Físico
O módulo `kinetics.py` implementa as leis físico-químicas fundamentais de dissolução da zincita ($\text{ZnO}$) por ácido sulfúrico ($\text{H}_2\text{SO}_4$).

Ele atua como o **motor cinético elementar** da modelagem mecanicista (\textit{White-Box}), convertendo o estado termodinâmico do fluido (temperatura e concentração de ácido livre) na **taxa linear de encolhimento do diâmetro da partícula** ($v(t) = dD/dt$).

---

## 2. Fundamentação Teórica e Equacionamento

### 2.1. Estequiometria e Balanço de Massa do Ácido
A reação de dissolução da zincita:
$$\text{ZnO}_{(s)} + \text{H}_2\text{SO}_{4(aq)} \longrightarrow \text{ZnSO}_{4(aq)} + \text{H}_2\text{O}_{(l)}$$
obedece à estequiometria molar unitária ($1:1$).

Pelo balanço molar no reator em batelada, a concentração de ácido sulfúrico livre instantâneo no licor, $C_{Af}(t)$, em função da conversão do sólido $X_{\text{Zn}}(t)$ e da razão molar estequiométrica $\eta$, é:
$$C_{Af}(t) = C_{A0} \left( 1 - \frac{X_{\text{Zn}}(t)}{\eta} \right)$$

* **Garantia Física:** A concentração nunca pode ser negativa: $C_{Af}(t) \ge 0$.
* **Patamar de Esgotamento ($\eta < 1{,}0$):** O ácido esgota-se ($C_{Af} \to 0$) quando a conversão atinge $X = \eta$.

### 2.2. Velocidade Linear de Retração do Diâmetro ($v(t)$)
De acordo com o Modelo do Núcleo Não Reagido (*Shrinking Core Model* -- SCM) com controle interfacial por reação química superficial (Balarini et al., 2025; Bortot Coelho, 2017), a velocidade com que o diâmetro da partícula decresce é:
$$v(t) \equiv -\frac{dD}{dt} = \frac{2 \, k_s}{\rho_s} \, C_{Af}(t) \quad [\mu\text{m/min}]$$

onde:
* $k_s = 18.000\ \mu\text{m/min}$ ($3{,}0 \times 10^{-4}\text{ m/s}$ ou $0{,}03\text{ cm/s}$ a $40^\circ\text{C}$).
* $\rho_s = 69{,}2\text{ mol/L}$ de sólido (densidade molar da zincita).

### 2.3. Formulação com Amortecimento Empírico (Tese Bortot Coelho, 2017)
Para representar a desaceleração provocada pelo acúmulo de sulfatos e passivação, a tese de 2017 introduziu o parâmetro ajustável $\alpha$:
$$v(t) = \frac{2}{\rho_s} \max \Big( k_s \, C_{Af}(t) - \alpha \, [C_{A0} - C_{Af}(t)],\ 0 \Big)$$
* Para $\alpha = 0{,}0$: Modelo mecanicista puro de primeiros princípios.
* Para $\alpha = 5500{,}0\ \mu\text{m/min}$: Modelo com amortecimento empírico calibrado na tese de 2017.

### 2.4. Dependência Térmica (Equação de Arrhenius)
$$k_s(T) = k_{s,\text{ref}} \exp \left[ -\frac{E_a}{R} \left( \frac{1}{T} - \frac{1}{T_{\text{ref}}} \right) \right]$$
* $E_a = 13{,}45\text{ kJ/mol}$ (Energia de ativação aparente; Balarini et al., 2025).
* $T_{\text{ref}} = 313{,}15\text{ K}$ ($40^\circ\text{C}$).
* $R = 8{,}31446\text{ J/(mol K)}$.

---

## 3. Funções Disponíveis

### `calcular_ks(temperatura_celsius: float) -> float`
Calcula o valor corrigido de $k_s$ em $\mu\text{m/min}$ para qualquer temperatura operacional.

### `concentracao_acido_livre(C_acid_0: float, X_zn: float, eta: float) -> float`
Calcula a concentração instantânea de $\text{H}_2\text{SO}_4$ livre em $\text{mol/L}$, garantindo que o valor permaneça $\ge 0$.

### `taxa_retricao_linear_v(C_acid_mol_L: float, ks_um_min: float, alpha: float, C_acid_0_mol_L: float) -> float`
Calcula a taxa linear $v(t) = dD/dt$ em $\mu\text{m/min}$.

---

## 4. Exemplo de Uso

```python
from kinetics import concentracao_acido_livre, taxa_retricao_linear_v, calcular_ks

# 1. Ajuste térmico para 50 °C
ks_50c = calcular_ks(temperatura_celsius=50.0)
print(f"ks a 50 °C: {ks_50c:.1f} um/min")

# 2. Avaliação de ácido livre para CA0 = 0.5 M, eta = 1.0 e conversão de 60%
caf = concentracao_acido_livre(C_acid_0=0.5, X_zn=0.60, eta=1.0)
print(f"Ácido livre residual: {caf:.4f} mol/L")

# 3. Taxa linear de dissolução da partícula
v = taxa_retricao_linear_v(C_acid_mol_L=caf, alpha=0.0)
print(f"Velocidade de retração do diâmetro: {v:.2f} um/min")
```

---

## 5. Como Testar
Execute no terminal:
```powershell
& .venv\Scripts\python.exe scripts/kinetics.py
```
*Saída esperada:*
```text
=== TESTE CINÉTICO BEM SUCEDIDO ===
Para CA0 = 0.5 mol/L, eta = 1.0 e X = 50.0%:
-> Ácido livre residual CAf: 0.2500 mol/L
-> Taxa de encolhimento (Física Pura): 130.06 um/min
-> Taxa de encolhimento (Tese alpha=5500): 90.32 um/min
```
