# README -- Modelo híbrido (`hybrid_model.py`)

## 1. Arquitetura: paralela cooperativa (Sansana et al., 2024)
$$X(t) = \mathcal{R}\big[\,X_{PBM}(t) + \widehat{\Delta X}_{ML}(\mathbf{z})\,\big],\qquad \Delta X = X_{exp} - X_{PBM}$$
* **Base física** (`kinetics.ConfigWB` + `population_balance`): White-Box principal = V2, ν = 1,06.
* **ML residual:** aprende ΔX a partir dos atributos $\mathbf{z}$ (abaixo).
* **Acoplador** $\mathcal{R}$ (`AcopladorFisico`): impõe as restrições físicas à saída **bruta** do ML.

Quando usar: o White-Box descreve o comportamento geral mas tem uma limitação estrutural (aqui, o patamar em η ≈ 1) e o domínio é conhecido e controlado.
O desempenho em condição **nova** depende do domínio (ver "Domínio de aplicabilidade").

## 2. Atributos do ML (`CONJUNTOS_ENTRADA`)
| Conjunto | Colunas | Observação |
| :--- | :--- | :--- |
| `'A'` | `tempo_min`, `razao_molar_eta_real`, `C_acid_0_mol_L` | só as três variáveis independentes do planejamento |
| **`'B'`** (principal) | A + `X_pbm` | o ML "vê" onde a teoria opera |
| `'C'` | B + `C_acid_pbm_mol_L`, `H_pbm_mol_L` | estado químico do licor segundo o PBM |
| `'legado'` | `tempo_min`, η nominal, `C_acid_0_mol_L`, `razao_SL_g_L`, `X_pbm`, `C_acid_pbm_mol_L` | só para reproduzir os números antigos |

O conjunto B foi escolhido pela **ablação** (`resultados/ablacao_entradas.csv`) com regra declarada: menor RMSE médio; havendo empate
(diferença pareada menor que o desvio-padrão da diferença entre repetições), vale o conjunto mais simples. $S/L$ foi retirado: é função
exata de $C_{A0}/\eta$. Usa-se η **real** (Eq. 4.3), não o nível nominal.

## 3. Acoplador físico (`AcopladorFisico(nu, modo='fisico')`)
1. $0 \le X \le X_{teto}$, com **$X_{teto} = \min(1, \eta_{real}/\nu)$** (mesma estequiometria do PBM). O teto antigo usava o η nominal (0,5) contra o η real do PBM (0,535): equivalia a ν = 1,07 **sem declarar**, e reduzia o RMSE do White-Box de 0,0443 para 0,0415.
2. $X(0) = 0$ (condição inicial conhecida).
3. Monotonia $dX/dt \ge 0$ por ensaio, pela **projeção L2 (regressão isotônica)** em ordem de tempo. O máximo acumulado usado antes é viesado para cima e só funcionava com as linhas agrupadas e ordenadas.

`diagnosticar(df, x_bruto)` audita a **saída bruta**: % de pontos acima do teto, abaixo de zero, passos com queda e $|X(0)|$ máximo. A
"violação 0%" da versão anterior era medida **depois** do clip, portanto tautológica.

## 4. Catálogo (`obter_especificacoes()`, chave = slug)
`gp` (C·Matérn 3/2 ARD + ruído, log(1+t)), `extra_trees`, `gbdt`, `random_forest`, `svr` (C, γ, ε por CV aninhada, alvo padronizado), `mlp_tanh` (arquitetura e α por CV aninhada).
Nome, slug, título e fábrica ficam num só lugar (`ModeloSpec`). `AjusteAninhado` escolhe hiperparâmetros por CV **interna** agrupada por ensaio,
usando só os ensaios de treino; ET, RF e GBDT usam os hiperparâmetros herdados do TCC (podem estar levemente otimistas: ver `sensibilidade` no relatório).
`obter_catalogo_modelos()` é o catálogo **legado** (9 regressores, 6 atributos), mantido para `entradas='legado'` e para o módulo do contínuo.

## 5. API
```python
from hybrid_model import ModeloHibridoCinza, construir_base, construir_base_condicao
from data_loader import carregar_dados_bancada

df = carregar_dados_bancada()
m = ModeloHibridoCinza('Gaussian Process (Kriging)')     # entradas='B', White-Box principal
base = m.preparar_dados_base(df)                         # recalcula o PBM; sem cache em disco
m.treinar(base)
x_hat = m.prever(base)                                   # com restrições físicas; avisa fora do domínio
x_hat, sd = m.prever_com_incerteza(base)                 # só GP: desvio-padrão do resíduo
dX = m.prever_residuo(base)                              # só o ML
cond = construir_base_condicao(eta_real=1.25, C_acid_0_mol_L=0.75, tempos_min=[0.5, 2, 15], config=m.config_wb)
m.prever(cond)                                           # condição nova (aqui η entre níveis: UserWarning de extrapolação)
```
`construir_base(df, config)` devolve **cópia** com `X_pbm`, `C_acid_pbm_mol_L`, `H_pbm_mol_L`, `delta_X`. Seu cache em memória guarda só a saída do PBM
(chave = entradas do PBM); `delta_X` sempre sai dos dados recebidos.

## 6. Domínio de aplicabilidade (`fora_do_dominio`)
$t \le 15$ min, η real $\in [0{,}53;\ 3{,}30]$, $C_{A0} \in [0{,}1;\ 1{,}5]$ mol/L, T = 40 °C. Fora dele `prever` emite `UserWarning`. **Entre** os 4
níveis ensaiados de η o ML também não tem informação (as validações de extrapolação mostram o efeito); o domínio verdadeiro é, portanto, os 4 níveis de η.

## 7. Compatibilidade
`continuous_reactor_model.py` (fora do escopo da revisão) usa `ModeloHibridoCinza(..., entradas='legado')` e roda como antes.

## 8. Testes
`tests/test_acoplador.py`, `tests/test_regressao.py` (reproduz os números anteriores com `entradas='legado'` + acoplador legado).
