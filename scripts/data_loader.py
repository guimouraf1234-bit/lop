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
    """
    128 observações (16 ensaios x 8 tempos) das Tabelas A1.1 e A1.4 da tese. Colunas relevantes:
    razao_molar_eta (nível nominal do planejamento), razao_molar_eta_real (Eq. 4.3, pelas massas),
    C_acid_0_mol_L, razao_SL_g_L, tempo_min, X_zn_exp.
    Atenção: 'C_acid_final_mol_L' NÃO é uma titulação: a tese o calcula a partir do pH medido, por especiação
    ideal de ácido sulfúrico puro (item 4.4.2 e Apêndice A4). O pH é a única medida independente.
    """
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
    """OBSOLETO: cache do resíduo do White-Box legado. O pipeline atual recalcula o PBM a cada uso."""
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


def carregar_respostas_finais():
    """
    Resposta final de cada ensaio (uma linha por ensaio), da Tabela A1.1 da tese:
    pH_0 (solução lixiviante), pH_final, X_zn_final, X_fe_final, Zn_licor_final_g_L e C_acid_final_mol_L.
    C_acid_final_mol_L é derivado do pH_final (ver `carregar_dados_bancada`), não é medida independente.
    Inclui eta nominal, eta real (Eq. 4.3), massa, volume e CA0.
    """
    caminhos = obter_caminhos()
    with open(caminhos["bancada"], "r", encoding="utf-8") as arquivo:
        dados_brutos = json.load(arquivo)
    from population_balance import calcular_eta
    linhas = []
    for chave, ens in dados_brutos["ensaios"].items():
        cond = ens["condicoes_operacionais"]
        resp = ens["respostas_finais"]
        linhas.append({
            "ensaio_id": int(chave),
            "razao_molar_eta": cond["razao_molar_eta"],
            "C_acid_0_mol_L": cond["C_acid_0_mol_L"],
            "massa_solido_g": cond["massa_solido_g"],
            "volume_solucao_L": cond["volume_solucao_L"],
            "razao_molar_eta_real": calcular_eta(cond["volume_solucao_L"], cond["C_acid_0_mol_L"], cond["massa_solido_g"]),
            "pH_0": resp["pH_0"],
            "pH_final": resp["pH_final"],
            "C_acid_final_mol_L": resp["C_acid_final_mol_L"],
            "Zn_licor_final_g_L": resp["Zn_licor_final_g_L"],
            "X_zn_final": resp["X_zn_final"],
            "X_fe_final": resp["X_fe_final"],
        })
    return pd.DataFrame(linhas).sort_values("ensaio_id").reset_index(drop=True)
