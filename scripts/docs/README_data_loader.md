# README -- Módulo de Carregamento de Dados (`data_loader.py`)

## 1. Visão Geral e Objetivo
O módulo `data_loader.py` é responsável pela ingestão, estruturação e disponibilização em memória de todos os dados experimentais e propriedades físico-químicas utilizados no projeto de modelagem da lixiviação de concentrado ustulado de zinco.

Ele atua como a **camada de dados (Data Layer)** do projeto, garantindo que os scripts mecanicistas e de aprendizado de máquina acessem as informações experimentais de forma padronizada via DataFrames do Pandas e dicionários tipados.

---

## 2. Estrutura dos Dados Carregados

O módulo consome arquivos JSON padronizados localizados no diretório `data/extracted/`:

1. **`ensaios_bancada.json` (16 ensaios de bancada em triplicata):**
   * Contém 128 registros planos $(t, X_{\text{Zn}})$.
   * Variáveis operacionais monitoradas:
     * `tempo_min`: Tempos discretos de amostragem ($0{,}0;\ 0{,}5;\ 1{,}0;\ 2{,}0;\ 3{,}0;\ 4{,}0;\ 5{,}0;\ 15{,}0\text{ min}$).
     * `X_zn_exp`: Conversão fracionária de zinco medida por Espectrometria de Absorção Atômica (EAA).
     * `razao_molar_eta`: Razão estequiométrica molar entre $\text{H}_2\text{SO}_4$ e $\text{ZnO}$ ($\eta \in \{0{,}5;\ 1{,}0;\ 1{,}5;\ 3{,}1\}$).
     * `C_acid_0_mol_L`: Concentração inicial de ácido ($C_{A0} \in \{0{,}1;\ 0{,}5;\ 1{,}0;\ 1{,}5\}\text{ mol/L}$).
     * `razao_SL_g_L`: Razão sólido/líquido em $\text{g/L}$.
     * `C_acid_final_mol_L` e `pH_final`: Concentração final de ácido e pH medidos.
   * *Fonte:* Tabelas A1.1 e A1.4 da dissertação de Bortot Coelho (2017).

2. **`parametros_processo.json`:**
   * Constantes cinéticas ($k_s$, $E_a$), propriedades do sólido ($\rho_s = 69{,}2\text{ mol/L}$, teor de $\text{ZnO} = 76{,}1\%$, fração de zinco solúvel $f_{\text{ZnO}} = 0{,}87$), e parâmetros granulométricos RRB ($m = 1{,}022$, $D_{63{,}2} = 41{,}65\ \mu\text{m}$).

3. **`ensaios_planta_piloto.json`:**
   * Dados contínuos da cascata de 3 reatores CSTR em série ($18\text{ L}$ cada) sob vazões de $0{,}41\text{ L/min}$ e $0{,}21\text{ L/min}$ (Tabelas A1.6 e A1.7).

Além disso, gerencia o acesso aos dados pré-calculados/processados em `data/processed/`:
4. **`dados_com_residuo.csv`:**
   * Dataset de 128 pontos contendo as variáveis de processo, a conversão física do PBM ($X_{\text{PBM}}$), a concentração ácida remanescente ($C_{A,\text{PBM}}$) e o resíduo experimental ($\Delta X$).

---

## 3. Funções Disponíveis

### `obter_caminhos() -> dict`
Descobre dinamicamente a localização absoluta dos diretórios e arquivos em `data/extracted/` e `data/processed/` a partir da posição do script, evitando erros de caminho relativo em diferentes sistemas operacionais.

### `carregar_dados_bancada() -> pd.DataFrame`
Carrega os 16 ensaios de bancada e retorna um `pd.DataFrame` contendo 128 linhas e 16 colunas com todos os atributos operacionais e respostas analíticas.

### `carregar_parametros() -> dict`
Retorna um dicionário hierárquico com todas as propriedades termodinâmicas, granulométricas e cinéticas da calcina de zinco e da solução lixiviante.

### `carregar_dados_com_residuo() -> Optional[pd.DataFrame]`
Carrega o conjunto consolidado de dados processados com resíduos (`data/processed/dados_com_residuo.csv`) para uso direto no treinamento e validação do modelo híbrido (*Grey-Box*).

---

## 4. Exemplo de Uso

```python
from data_loader import carregar_dados_bancada, carregar_parametros

# 1. Carrega os ensaios de bancada
df = carregar_dados_bancada()
print(f"Total de observações carregadas: {len(df)}")

# Filtrar ensaios com razão estequiométrica unitária (eta = 1.0)
df_eta1 = df[df["razao_molar_eta"] == 1.0]
print(df_eta1[["ensaio_id", "tempo_min", "X_zn_exp"]].head())

# 2. Carrega as constantes do processo
params = carregar_parametros()
ks = params["parametros_cineticos_reacao_quimica"]["constante_cinetica_ks"]["valor_um_min"]
print(f"Constante cinética de referência: {ks} um/min")
```

---

## 5. Como Testar
No terminal PowerShell com o ambiente virtual ativado:
```powershell
& .venv\Scripts\python.exe scripts/data_loader.py
```
*Saída esperada:* Confirmação de 128 observações e exibição das 5 primeiras linhas do DataFrame.
