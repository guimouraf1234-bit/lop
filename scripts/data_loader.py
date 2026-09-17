import json
import os
import pandas as pd

def obter_caminhos():
    diretorio_atual = os.path.dirname(os.path.abspath(__file__))
    pasta_dados = os.path.abspath(os.path.join(diretorio_atual, "..", "data", "extracted"))
    
    return {
        "bancada": os.path.join(pasta_dados, "ensaios_bancada.json"),
        "parametros": os.path.join(pasta_dados, "parametros_processo.json"),
        "piloto": os.path.join(pasta_dados, "ensaios_planta_piloto.json")
    }

def carregar_dados_bancada():
    caminhos = obter_caminhos()
    with open(caminhos["bancada"], "r", encoding="utf-8") as arquivo:
        dados_brutos = json.load(arquivo)
    df = pd.DataFrame(dados_brutos["registros_planos"])
    return df

def carregar_parametros():
    caminhos = obter_caminhos()
    with open(caminhos["parametros"], "r", encoding="utf-8") as arquivo:
        parametros = json.load(arquivo)
    return parametros
