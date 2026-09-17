import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

def gerar_diagrama_arquitetura():
    # Cria figura com proporção ampla e alta resolução
    fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8.5)
    ax.axis('off')
    
    # Fundo branco puro
    fig.patch.set_facecolor('#ffffff')
    ax.set_facecolor('#ffffff')

    # Cores da paleta profissional
    cor_entrada = "#f1f3f4"
    borda_entrada = "#5f6368"
    
    cor_white = "#e8f0fe"
    borda_white = "#1a73e8"
    
    cor_feat = "#fef7e0"
    borda_feat = "#f9ab00"
    
    cor_black = "#fce8e6"
    borda_black = "#d93025"
    
    cor_hybrid = "#e6f4ea"
    borda_hybrid = "#137333"

    # Função auxiliar para desenhar caixas arredondadas
    def desenhar_bloco(x, y, w, h, titulo, subtitulo, cor_fundo, cor_borda, fontsize_tit=10, fontsize_sub=8.5):
        box = patches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.15,rounding_size=0.2",
            facecolor=cor_fundo,
            edgecolor=cor_borda,
            linewidth=1.8,
            zorder=2
        )
        ax.add_patch(box)
        
        # Título
        ax.text(x + w/2, y + h - 0.35, titulo, ha='center', va='center',
                fontsize=fontsize_tit, fontweight='bold', color='#202124', zorder=3)
        
        # Subtítulo / Equações
        ax.text(x + w/2, y + (h - 0.4)/2, subtitulo, ha='center', va='center',
                fontsize=fontsize_sub, color='#3c4043', multialignment='center', zorder=3)

    # 1. BLOCO DE ENTRADAS (Esquerda)
    desenhar_bloco(
        0.5, 2.5, 2.4, 4.0,
        "VARIÁVEIS DE PROCESSO",
        "• Tempo: t [min]\n• Razão Molar: η\n• Ácido Inicial: CA0 [mol/L]\n• Temperatura: T [°C]\n• Granulometria: RRB\n  (m = 1.022, D63.2 = 41.65 µm)",
        cor_entrada, borda_entrada, 10, 8.5
    )

    # 2. RAMO SUPERIOR: WHITE-BOX (Física Pura - PBM)
    desenhar_bloco(
        4.0, 5.2, 4.6, 2.6,
        "1. RAMO WHITE-BOX (Mecanicista - PBM)",
        "Taxa Linear: v(t) = (2 ks / ρs) · CAf(t)\nDistribuição RRB: f0(D)\nEDO do Balanço Populacional (Método das Características):\n" +
        r"$\mathbf{X_{PBM}(t) = 1 - \frac{\int (D - \Delta D)^3 f_0(D) dD}{I_0}}$",
        cor_white, borda_white, 10, 8.5
    )

    # 3. ENGENHARIA DE ATRIBUTOS FÍSICOS (Domain-Knowledge Features)
    desenhar_bloco(
        4.0, 0.8, 4.6, 2.8,
        "2. ENGENHARIA DE ATRIBUTOS FÍSICOS",
        "Atributos termodinâmicos alimentados ao ML:\n" +
        r"• Força-Motriz Residual: $C_{Af,teo} = C_{A0}(1 - X_{PBM}/\eta)$" + "\n" +
        r"• Esgotamento Estequiométrico: $X_{PBM} / \eta$" + "\n" +
        r"• Retração Acumulada do Diâmetro: $\Delta D(t)$",
        cor_feat, borda_feat, 10, 8.5
    )

    # 4. RAMO INFERIOR: BLACK-BOX (Machine Learning)
    desenhar_bloco(
        9.5, 0.8, 4.0, 2.8,
        "3. RAMO BLACK-BOX (Machine Learning)",
        "Modelagem do Resíduo Físico:\n" +
        "Random Forest / MLP / SVR\n\n" +
        r"Alvo: $\mathbf{\Delta X = X_{exp} - X_{PBM}}$" + "\n" +
        r"Predição: $\mathbf{\widehat{\Delta X}_{ML}(z)}$",
        cor_black, borda_black, 10, 8.5
    )

    # 5. ACOPLADOR HÍBRIDO GREY-BOX (Direita)
    desenhar_bloco(
        9.5, 5.0, 4.0, 2.8,
        "4. ACOPLAMENTO HÍBRIDO GREY-BOX",
        "Fusão Física + Aprendizado de Máquina:\n" +
        r"$\mathbf{X_{raw} = X_{PBM}(t) + \widehat{\Delta X}_{ML}(z)}$" + "\n\n" +
        "Restrição de Conservação de Massa:\n" +
        r"$\mathbf{X_{Híbrido} = clip(X_{raw},\ 0.0,\ 1.0)}$" + "\n\n" +
        r"$\mathbf{R^2 = 0{,}9979}$ | $\mathbf{RMSE = 0{,}0147}$",
        cor_hybrid, borda_hybrid, 10, 8.5
    )

    # SETAS E CONEXÕES
    # Entrada -> White-Box
    ax.annotate("", xy=(4.0, 6.5), xytext=(2.9, 5.2),
                arrowprops=dict(arrowstyle="-|>", color="#1a73e8", lw=2, mutation_scale=15))
    ax.text(3.4, 6.1, "Física", color="#1a73e8", fontsize=8.5, fontweight="bold")

    # Entrada -> Feature Engineering
    ax.annotate("", xy=(4.0, 2.2), xytext=(2.9, 3.8),
                arrowprops=dict(arrowstyle="-|>", color="#f9ab00", lw=2, mutation_scale=15))
    ax.text(3.2, 2.8, "Condições", color="#b06000", fontsize=8.5, fontweight="bold")

    # White-Box -> Feature Engineering (XPBM desce)
    ax.annotate("", xy=(6.3, 3.6), xytext=(6.3, 5.2),
                arrowprops=dict(arrowstyle="-|>", color="#1a73e8", lw=2, mutation_scale=15))
    ax.text(6.4, 4.4, r"$X_{PBM}(t)$", color="#1a73e8", fontsize=9, fontweight="bold")

    # White-Box -> Acoplador Híbrido (XPBM vai direto para a soma)
    ax.annotate("", xy=(9.5, 6.5), xytext=(8.6, 6.5),
                arrowprops=dict(arrowstyle="-|>", color="#1a73e8", lw=2.5, mutation_scale=18))
    ax.text(9.05, 6.75, r"$X_{PBM}$", ha="center", color="#1a73e8", fontsize=9, fontweight="bold")

    # Feature Engineering -> Black-Box
    ax.annotate("", xy=(9.5, 2.2), xytext=(8.6, 2.2),
                arrowprops=dict(arrowstyle="-|>", color="#d93025", lw=2, mutation_scale=15))
    ax.text(9.05, 2.45, "Features", ha="center", color="#d93025", fontsize=8.5, fontweight="bold")

    # Black-Box -> Acoplador Híbrido (Delta_X sobe para a soma)
    ax.annotate("", xy=(11.5, 5.0), xytext=(11.5, 3.6),
                arrowprops=dict(arrowstyle="-|>", color="#d93025", lw=2.5, mutation_scale=18))
    ax.text(11.6, 4.3, r"$\widehat{\Delta X}_{ML}$", color="#d93025", fontsize=10, fontweight="bold")

    # Título Principal do Diagrama
    plt.suptitle("ARQUITETURA DA MODELAGEM HÍBRIDA GREY-BOX (PBM + RESÍDUO PARALELO)",
                 fontsize=13, fontweight='bold', color='#1a0dab', y=0.96)
    plt.title("Framework de Engenharia Química para Lixiviação de Concentrado Ustulado de Zinco -- TCC UFMG",
              fontsize=9.5, color='#5f6368', pad=10)

    # Salva figura em alta resolução
    os.makedirs("scripts", exist_ok=True)
    caminho_salvar = "scripts/arquitetura_greybox.png"
    plt.savefig(caminho_salvar, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Diagrama salvo com sucesso em: {caminho_salvar}")

if __name__ == "__main__":
    gerar_diagrama_arquitetura()
