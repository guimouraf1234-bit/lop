# README -- Carregamento de dados (`data_loader.py`)

## 1. Papel
Camada de dados do projeto: lê os JSON de `data/extracted/` e devolve DataFrames/dicionários padronizados para o PBM, o
híbrido e o pipeline de validação.

## 2. Dados

| Arquivo | Conteúdo | Origem |
| :--- | :--- | :--- |
| `ensaios_bancada.json` | 16 ensaios × 8 tempos = **128 pontos** $(t, X_{Zn})$ e a resposta final de cada ensaio | Tabelas A1.1 e A1.4 de Bortot Coelho (2017); médias de triplicatas |
| `parametros_processo.json` | constantes cinéticas, propriedades da calcina, RRB | tese e Balarini et al. (2025) |
| `ensaios_planta_piloto.json` | cascata de 3 CSTR (fora do escopo desta revisão) | Tabelas A1.6 e A1.7 |

Delineamento: $\eta \in \{0{,}5;\ 1{,}0;\ 1{,}5;\ 3{,}1\}$ × $C_{A0} \in \{0{,}1;\ 0{,}5;\ 1{,}0;\ 1{,}5\}$ mol/L, T = 40 °C, 1000 rpm, 400 mL, tempos
0; 0,5; 1; 2; 3; 4; 5; 15 min.

> **Atenção (revisão):** `C_acid_final_mol_L` **não é uma titulação**. A tese (item 4.4.2 e Apêndice A4) o calcula a partir
> do **pH medido**, por especiação ideal de ácido sulfúrico puro (balanço de cargas sem Zn²⁺). O **pH é a única medida independente**
> do estado ácido; por isso a validação termodinâmica do White-Box usa pH, não $C_{Af}$. As duas colunas só fornecem a mesma informação.
>
> `ensaios_planta_piloto.json` informa 18 L por reator (volume_total 54 L), mas a dissertação (item 4.7.2) diz **6 L** por reator
> (τ = 14,6 e 28,6 min por reator); o código do contínuo usa 6 L. Corrigir o JSON antes de usá-lo.

## 3. Funções
| Função | Retorno |
| :--- | :--- |
| `obter_caminhos()` | caminhos absolutos dos dados (a partir da posição do script) |
| `carregar_dados_bancada()` | DataFrame de 128 linhas: operação, `X_zn_exp` e `razao_molar_eta_real` (Eq. 4.3, pelas massas; o `razao_molar_eta` é o nível nominal 0,5/1,0/1,5/3,1) |
| `carregar_respostas_finais()` | uma linha por ensaio: `pH_0`, `pH_final`, `C_acid_final_mol_L` (derivado do pH), `Zn_licor_final_g_L`, `X_zn_final`, `X_fe_final`, η nominal e real |
| `carregar_parametros()` | dicionário de `parametros_processo.json` |
| `carregar_dados_piloto()` | dicionário do piloto |
| `carregar_dados_com_residuo()` | **obsoleto**: cache do resíduo do White-Box antigo, em `arquivo/`. O pipeline recalcula o PBM a cada execução |

## 4. Exemplo
```python
from data_loader import carregar_dados_bancada, carregar_respostas_finais

df = carregar_dados_bancada()
print(len(df), df['ensaio_id'].nunique())                     # 128 16
print(df[df['razao_molar_eta'] == 1.0][['ensaio_id', 'tempo_min', 'X_zn_exp']].head(3))
resp = carregar_respostas_finais()
print(resp[['ensaio_id', 'pH_0', 'pH_final', 'X_zn_final']].head(3))
```

## 5. Teste
```powershell
& .venv\Scripts\python.exe -m pytest tests -q          # usa estes carregadores em todos os testes
```
