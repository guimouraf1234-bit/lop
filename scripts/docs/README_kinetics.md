# README -- Cinética, estequiometria e especiação (`kinetics.py`)

## 1. Papel
`kinetics.py` reúne as leis do White-Box que não dependem da granulometria: a velocidade de dissolução, o balanço de
ácido e a **especiação do H₂SO₄**. Todo parâmetro tem proveniência declarada no código e nenhum é ajustado aos 128
pontos de conversão. A configuração física completa é um objeto congelado, `ConfigWB`, gravado junto de todo resultado.

## 2. `ConfigWB` (configuração do White-Box)

| Campo | Valores | Significado |
| :--- | :--- | :--- |
| `rrb` | `'massica'` (padrão), `'numerica'` | RRB integrada em base mássica (correta) ou como a Eq. 3.49 da tese |
| `ks_arrhenius` | `True` / `False` | $k_s$ referido a 30 °C e corrigido por Arrhenius, ou fixo em 1,8·10⁴ µm/min |
| `forca_motriz` | `'H+'` (padrão), `'H2SO4'` | força motriz = [H⁺] da especiação, ou a molaridade do ácido (tese) |
| `especiacao` | `'V0'`, `'V1'`, `'V2'` | ver seção 3 (só vale com `forca_motriz = 'H+'`) |
| `nu` | float ≥ 0 | mol de H₂SO₄ consumidos por mol de ZnO dissolvido (seção 4) |
| `T_celsius`, `alpha` | 40, 0 | temperatura e parâmetro $\alpha$ da Eq. 4.1 (0 = Eq. 3.42) |

Fábricas: `ConfigWB.legado()` (White-Box do TCC antes da revisão: V0, ν = 1), `ConfigWB.tese(alpha)` (Apêndice A2.1),
`cfg.com_nu(nu)`. O White-Box **principal** do trabalho é `hybrid_model.CONFIG_WB_PRINCIPAL = ConfigWB(especiacao='V2', nu=1.06)`.

## 3. Especiação do H₂SO₄ (concentrações, sem coeficientes de atividade)
Com prótons totais $H_T = 2\,C_{Af}$ e sulfato total $S_T = C_{A0}$ (o sulfato do ZnSO₄ formado fica em solução):

$$[\text{H}^+] + [\text{HSO}_4^-] = H_T,\qquad [\text{HSO}_4^-]+[\text{SO}_4^{2-}]+[\text{ZnSO}_4^0] = S_T,\qquad
K_2 = \frac{[\text{H}^+][\text{SO}_4^{2-}]}{[\text{HSO}_4^-]},\qquad [\text{ZnSO}_4^0] = K_z\,[\text{Zn}^{2+}][\text{SO}_4^{2-}]$$

| Variante | $K_2$ | Par iônico ZnSO₄⁰ | Fonte |
| :--- | :--- | :--- | :--- |
| **V0** (TCC) | 0,0099 mol/L, constante | não | reproduz o [H⁺] = 0,0689 mol/L de Balarini (6 g/L, 30 °C) |
| **V1** | van't Hoff: 0,01096 (25 °C), 0,00949 (30 °C), **0,00720 (40 °C)** | não | Apêndice A4 da tese: pKa₂ = 1,96 e ΔH° = −5,2 kcal/mol |
| **V2** | igual a V1 | sim: $K_z$ = 201 (25 °C), 226 (40 °C) L/mol | IUPAC (Powell et al., *Pure Appl. Chem.* 85, 2249, 2013, Parte 5): log₁₀K° = 2,304 ± 0,042 e ΔrH° = +6,0 ± 0,5 kJ/mol |

Sem par iônico, a solução é a raiz da quadrática (`concentracao_h_livre`). Com ele, resolve-se a equação monótona
$f(s) = K_2 b - (H_T - b)s$ em $s = [\text{SO}_4^{2-}]$ por `brentq` (`especiar_escalar`).

**Ressalvas honestas**
* Os valores IUPAC foram obtidos de um resumo secundário (o texto primário não pôde ser aberto): **conferir no relatório antes de citar**.
* $K_z$ é de força iônica zero e é aplicada **sem correção de atividade**; em licor concentrado (I > 1 mol/L) a constante aparente é bem menor. V2 é, portanto, uma **cota superior** do efeito do par iônico, e V0/V1 a cota inferior (sem par).
* A tese lista pKa₂ = 1,96 em um ponto e Ka₂ = 12,383·10⁻³ (pKa₂ = 1,91) em outro; usa-se 1,96, que reproduz Balarini a 30 °C (diferença de −0,4%).

**Validação independente (pH medido).** O pH é a única medida independente do estado ácido (o "ácido livre final" da Tabela A1.1
é calculado a partir do pH, por especiação ideal de ácido puro, ver `docs/validade_modelo_hibrido.md`). V2 é a única
variante coerente com o pH final: ver `resultados/selecao_whitebox.csv`.

## 4. Fator de fechamento estequiométrico ν
$C_{Af}(t) = C_{A0}\,(1 - \nu X/\eta)$ (Eq. 3.51 com ν = 1). Nos ensaios com η = 0,5 o ácido é totalmente consumido (pH final 5–6)
com $X_{final} = 0{,}50 < \eta_{real} = 0{,}535$, logo $\nu = \eta_{real}/X_{final} \approx 1{,}07$: há outro consumo de ácido
(óxidos de Ca, Cu, Cd e Mn; ferrita) ou erro de fechamento do balanço. `nu_a_priori()` = **1,06**, estimado só pela composição da
calcina (ver docstring); `validacao.estimar_nu` o estima dentro de cada dobra de CV. O teto do acoplador passa a ser $\min(1, \eta_{real}/\nu)$.

## 5. Funções
| Função | Descrição |
| :--- | :--- |
| `calcular_ks(T)` | $k_s$ em µm/min (18000 a 30 °C; 21344 a 40 °C) |
| `k2_hso4(T, especiacao)`, `k_znso4(T)` | constantes de equilíbrio por van't Hoff |
| `concentracao_acido_livre(C_A0, X, eta, nu=1.0)` | $C_{Af}$ (≥ 0) |
| `concentracao_h_livre(C_A0, C_Af, K2)` | [H⁺] sem par iônico (V0/V1) |
| `especiar_escalar(S, H, K2, zn, kz)`, `especiar(C_A0, C_Af, cfg, zn)` | especiação completa (h, HSO₄⁻, SO₄²⁻, ZnSO₄⁰) |
| `ph_de_h(h)` | −log₁₀[H⁺] (NaN se h = 0) |
| `estado_quimico(C_A0, X, eta, cfg)` | $(C_{Af},\ \text{força motriz},\ [\text{H}^+])$ usado pelo PBM |
| `taxa_retricao_linear_v(...)` | $v = -dD/dt = \tfrac{2}{\rho_s}[k_s C - \alpha(C_{A0}-C_{Af})]$ |

## 6. Exemplo (CA0 = 1,0 mol/L, X = 0,85, η = 1,069, ν = 1,06; o ensaio 9 mediu pH_f = 0,9)
```python
from kinetics import ConfigWB, concentracao_acido_livre, especiar, ph_de_h
for esp in ('V0', 'V1', 'V2'):
    cfg = ConfigWB(especiacao=esp, nu=1.06)
    caf = concentracao_acido_livre(1.0, 0.85, 1.069, cfg.nu)
    e = especiar(1.0, caf, cfg, zn_dissolvido=0.85 * 1.0 / 1.069)
    print(esp, f"CAf = {caf:.4f}  [H+] = {e['h']:.4f}  pH = {ph_de_h(e['h']):.2f}")
```
```text
V0 CAf = 0.1572  [H+] = 0.0044  pH = 2.35
V1 CAf = 0.1572  [H+] = 0.0033  pH = 2.49
V2 CAf = 0.1572  [H+] = 0.0541  pH = 1.27
```
Sem o par iônico o [H⁺] fica ~12× abaixo do valor de V2 (e ~30× abaixo do medido, 10⁻⁰·⁹ = 0,126); com ele, o sulfato (0,705 mol/L) fica quase todo
como ZnSO₄⁰ e o H⁺ deixa de ser tamponado como HSO₄⁻.

## 7. Testes
`tests/test_especiacao.py` (balanços, limites, [H⁺] de Balarini, pH₀ dos 16 ensaios) e `tests/test_pbm.py`.
```powershell
& .venv\Scripts\python.exe -m pytest tests/test_especiacao.py -q
```
