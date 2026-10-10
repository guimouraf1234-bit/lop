# README -- Módulo de Cinética Química Heterogênea (`kinetics.py`)

## 1. Visão Geral e Papel Físico
O módulo `kinetics.py` implementa as leis cinéticas de dissolução da zincita ($\text{ZnO}$) por ácido sulfúrico ($\text{H}_2\text{SO}_4$), a partir do modelo de Bortot Coelho (2017) (item 4.3 e Tabela 5.9 da dissertação), com as correções de consistência descritas em `docs/modelagem_white_box.md`: $k_s$ referido a 30 °C e força motriz $[\text{H}^+]$ por especiação.

Ele atua como o **motor cinético elementar** da modelagem mecanicista (\textit{White-Box}), convertendo o estado do licor na **taxa de variação do diâmetro da partícula** ($v = -dD/dt$).

---

## 2. Fundamentação Teórica e Equacionamento

### 2.1. Estequiometria e Balanço de Massa do Ácido
$$\text{ZnO}_{(s)} + \text{H}_2\text{SO}_{4(aq)} \longrightarrow \text{ZnSO}_{4(aq)} + \text{H}_2\text{O}_{(l)} \qquad (4.2)$$

Pelo balanço molar no reator batelada (Eq. 3.51), com $X_{\text{Zn}}$ a conversão da zincita (Eq. 4.5) e $\eta$ a razão molar H₂SO₄ / ZnO calculada com as massas do ensaio (Eq. 4.3):
$$C_{Af}(t) = C_{A0} \left( 1 - \frac{X_{\text{Zn}}(t)}{\eta} \right), \qquad C_{Af} \ge 0$$

### 2.2. Concentração de H⁺ por Especiação
$\text{HSO}_4^- \rightleftharpoons \text{H}^+ + \text{SO}_4^{2-}$, com $K_2 = 0{,}0099\ \text{mol/L}$ (`K2_HSO4_MOL_L`), valor que reproduz o $[\text{H}^+] = 0{,}0689\ \text{mol/L}$ usado por Balarini para 6 g/L de H₂SO₄. Com prótons totais $H_T = 2\,C_{Af}$ e sulfato total $S_T = C_{A0}$ (o ZnSO₄ formado mantém o sulfato em solução):
$$[\text{HSO}_4^-] = \frac{(S_T + H_T + K_2) - \sqrt{(S_T + H_T + K_2)^2 - 4\,S_T H_T}}{2}, \qquad [\text{H}^+] = H_T - [\text{HSO}_4^-]$$

### 2.3. Taxa de Variação do Diâmetro
$$v(t) \equiv -\frac{dD}{dt} = \frac{2}{\rho_s} \Big[ k_s \, C - \alpha \, \big(C_{A0} - C_{Af}(t)\big) \Big] \quad [\mu\text{m/min}]$$

com $C = [\text{H}^+]$ (argumento `C_motriz_mol_L`) ou, se ele não for informado, $C = C_{Af}$ (formulação da tese).

| Parâmetro | Valor | Fonte |
| :--- | :---: | :--- |
| $k_s$ a 30 °C, base [H⁺] (`KS_REF_UM_MIN`) | $1{,}8 \times 10^4\ \mu\text{m/min}$ | Balarini (2009) |
| $E_a$ (`EA_KJ_MOL`) | $13{,}45\ \text{kJ/mol}$ | Balarini (2009) |
| $\rho_s$ (`RHO_S_MOL_L`) | $69{,}2\ \text{mol/L}$ | Tabela 5.9 |
| $\alpha$ da tese (`ALPHA_TESE_UM_MIN`) | $5{,}5 \times 10^3\ \mu\text{m/min}$ | Tabela 5.9 (só referência) |

* $\alpha = 0$ (padrão): Eq. 3.42, White-Box do modelo híbrido.
* $\alpha = 5500$: Eq. 4.1 da tese, usada como referência nos gráficos.

### 2.4. Dependência Térmica
$$k_s(T) = k_{s,30^\circ\text{C}} \exp \left[ -\frac{E_a}{R} \left( \frac{1}{T} - \frac{1}{303{,}15} \right) \right]$$
A 40 °C: $k_s = 21\,344\ \mu\text{m/min}$.

---

## 3. Funções Disponíveis

### `calcular_ks(temperatura_celsius: float) -> float`
Retorna $k_s$ em $\mu\text{m/min}$ na base [H⁺] (18000 a 30 °C; 21344 a 40 °C).

### `concentracao_acido_livre(C_acid_0: float, X_zn: float, eta: float) -> float`
Concentração de $\text{H}_2\text{SO}_4$ livre em mol/L (Eq. 3.51), garantindo valor $\ge 0$.

### `concentracao_h_livre(C_acid_0: float, C_acid_livre: float, K2: float = 0.0099) -> float`
Concentração de $\text{H}^+$ em mol/L pela especiação $\text{HSO}_4^-/\text{SO}_4^{2-}$.

### `taxa_retricao_linear_v(C_acid_mol_L, ks_um_min, alpha=0.0, C_acid_0_mol_L=0.0, C_motriz_mol_L=None) -> float`
Calcula $v$ em $\mu\text{m/min}$. O termo de $\alpha$ só atua quando `C_acid_0_mol_L` é informado.

---

## 4. Exemplo de Uso

```python
from kinetics import calcular_ks, concentracao_acido_livre, concentracao_h_livre, taxa_retricao_linear_v

ks = calcular_ks(40.0)                                     # 21344 um/min
caf = concentracao_acido_livre(C_acid_0=0.5, X_zn=0.50, eta=1.0)
h = concentracao_h_livre(C_acid_0=0.5, C_acid_livre=caf)
v = taxa_retricao_linear_v(C_acid_mol_L=caf, ks_um_min=ks, C_acid_0_mol_L=0.5, C_motriz_mol_L=h)
print(f"CAf = {caf:.3f} mol/L | [H+] = {h:.4f} mol/L | v = {v:.2f} um/min")
```

---

## 5. Como Testar
```powershell
& .venv\Scripts\python.exe -c "import sys; sys.path.insert(0,'scripts'); from kinetics import *; ks=calcular_ks(40.0); c=concentracao_acido_livre(0.5,0.5,1.0); h=concentracao_h_livre(0.5,c); print(ks, c, h, taxa_retricao_linear_v(c, ks, C_acid_0_mol_L=0.5, C_motriz_mol_L=h))"
```
*Saída esperada (CA0 = 0,5 mol/L, η = 1,0, X = 50%):*
```text
ks(40 °C) ≈ 21344,1 um/min | CAf = 0,25 mol/L | [H+] ≈ 0,0656 mol/L | v ≈ 40,46 um/min
```
Com ZnSO₄ = 0,25 mol/L no licor, só 13% dos prótons do ácido livre restante estão como H⁺: é o tamponamento pelo sulfato.
