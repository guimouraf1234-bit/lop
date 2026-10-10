import json
import os
import pandas as pd

def obter_caminhos():
    diretorio_atual = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.abspath(os.path.join(diretorio_atual, "..", "data", "extracted"))
    pasta_processados = os.path.abspath(os.path.join(diretorio_atual, "..", "data", "processed"))
    
    return {
        "bancada": os.path.join(pasta_dados, "ensaios_bancada.json"),
        "parametros": os.path.join(pasta_dados, "parametros_processo.json"),
        "piloto": os.path.join(pasta_dados, "ensaios_planta_piloto.json"),
        "pasta_processados": pasta_processados,
        "dados_com_residuo": os.path.join(pasta_processados, "dados_com_residuo.csv")
    }

def carregar_dados_bancada():
    caminhos = obter_caminhos()
    with open(caminhos["bancada"], "r", encoding="utf-8") as arquivo:
        dados_brutos = json.load(arquivo)
    df = pd.DataFrame(dados_brutos["registros_planos"])
    # Razão molar recalculada pelas massas (Eq. 4.3 da tese), usada no PBM;
    # 'razao_molar_eta' guarda o valor nominal do planejamento (0,5 / 1,0 / 1,5 / 3,1)
    from population_balance import calcular_eta
    df["razao_molar_eta_real"] = [
        calcular_eta(v, c, m) for v, c, m in
        zip(df["volume_solucao_L"], df["C_acid_0_mol_L"], df["massa_solido_g"])
    ]
    return df

def carregar_parametros():
    caminhos = obter_caminhos()
    with open(caminhos["parametros"], "r", encoding="utf-8") as arquivo:
        parametros = json.load(arquivo)
    return parametros

def carregar_dados_com_residuo():
    caminhos = obter_caminhos()
    caminho_csv = caminhos["dados_com_residuo"]
    if os.path.exists(caminho_csv):
        return pd.read_csv(caminho_csv)
    caminho_legado = os.path.abspath(os.path.join(os.path.dirname(__file__), "dados_com_residuo.csv"))
    if os.path.exists(caminho_legado):
        return pd.read_csv(caminho_legado)
    return None

def carregar_dados_piloto():
    caminhos = obter_caminhos()
    with open(caminhos["piloto"], "r", encoding="utf-8") as arquivo:
        dados_piloto = json.load(arquivo)
    return dados_piloto
