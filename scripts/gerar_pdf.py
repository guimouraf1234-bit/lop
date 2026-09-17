"""
gerar_pdf.py
Converte scripts/Arquitetura.md para um PDF com padrão acadêmico sóbrio (estilo ABNT / artigo científico)
utilizando Microsoft Edge Headless com KaTeX, Mermaid.js e Marked.js.
Elimina a estética de "IA / template web" (sem azuis saturados, sem bordas coloridas).
"""

import os
import subprocess
import sys
from pypdf import PdfReader

def converter_markdown_para_pdf():
    diretorio_scripts = os.path.dirname(os.path.abspath(__file__))
    caminho_md = os.path.join(diretorio_scripts, "Arquitetura.md")
    caminho_html = os.path.join(diretorio_scripts, "_temp_arquitetura.html")
    caminho_pdf = os.path.join(diretorio_scripts, "Arquitetura.pdf")

    with open(caminho_md, "r", encoding="utf-8") as f:
        conteudo_md = f.read()

    # Prepara o template HTML de estilo acadêmico (estilo artigo/monografia da UFMG)
    html_template = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <title>Proposta de Arquitetura: White-Box ao Grey-Box</title>
    
    <!-- Marked.js para Markdown -->
    <script src="https://cdn.jsdelivr.net/npm/marked@9.1.6/marked.min.js"></script>
    
    <!-- Mermaid.js para diagramas -->
    <script src="https://cdn.jsdelivr.net/npm/mermaid@10.6.1/dist/mermaid.min.js"></script>
    
    <!-- KaTeX para renderização matemática fiel ao LaTeX -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.css">
    <script src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/katex.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/katex@0.16.9/dist/contrib/auto-render.min.js"></script>

    <style>
        @page {
            size: A4;
            margin: 15mm 18mm 15mm 18mm;
        }
        
        body {
            font-family: 'Times New Roman', 'Times', 'Cambria', 'Georgia', serif;
            font-size: 10pt;
            line-height: 1.42;
            color: #111111;
            background: #ffffff;
            margin: 0;
            padding: 0;
        }

        /* Títulos sóbrios no padrão de relatório acadêmico/ABNT */
        h1 {
            font-size: 15pt;
            font-weight: bold;
            color: #000000;
            text-align: center;
            margin-top: 0;
            margin-bottom: 6px;
            line-height: 1.25;
            border-bottom: none;
        }

        /* Subtítulo institucional / autores */
        h1 + p {
            text-align: center;
            font-size: 9.5pt;
            color: #222222;
            margin-bottom: 12px;
            line-height: 1.35;
        }

        h2 {
            font-size: 12pt;
            font-weight: bold;
            color: #000000;
            margin-top: 16px;
            margin-bottom: 6px;
            text-transform: uppercase;
            letter-spacing: 0.3px;
            border-bottom: 0.5pt solid #000000;
            padding-bottom: 2px;
            page-break-after: avoid;
        }

        h3 {
            font-size: 10.5pt;
            font-weight: bold;
            color: #000000;
            margin-top: 10px;
            margin-bottom: 4px;
            page-break-after: avoid;
        }

        p {
            margin-top: 0;
            margin-bottom: 6px;
            text-align: justify;
            text-indent: 12mm;
        }

        /* Parágrafos após títulos ou tabelas não precisam de recuo de primeira linha */
        h1 + p, h2 + p, h3 + p, .tabela-bloco p {
            text-indent: 0;
        }

        strong {
            font-weight: bold;
            color: #000000;
        }

        hr {
            border: 0;
            border-top: 0.5pt solid #888888;
            margin: 12px 0;
        }

        /* =========================================================
           TABELAS NO PADRÃO ABNT / BOOKTABS (Artigos Científicos)
           Regras: 
           - Linha superior e inferior espessas (1.5pt)
           - Linha de cabeçalho (1pt)
           - SEM bordas verticais nas laterais (clássico ABNT/Elsevier)
           - Fundo branco/neutro, sem cores azuis
           ========================================================= */
        .tabela-bloco {
            page-break-inside: avoid !important;
            break-inside: avoid !important;
            margin: 10px 0;
        }

        .tabela-bloco h3 {
            margin-top: 6px !important;
            margin-bottom: 1px !important;
            font-size: 9.5pt !important;
            font-weight: bold;
            text-align: left;
            text-indent: 0;
        }

        .tabela-bloco p {
            margin-bottom: 4px !important;
            font-size: 8pt !important;
            font-style: italic;
            color: #444444;
            text-align: left;
            text-indent: 0;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin: 4px 0 6px 0;
            font-size: 7.6pt;
            font-family: 'Times New Roman', serif;
            border-top: 1.2pt solid #000000;
            border-bottom: 1.2pt solid #000000;
            page-break-inside: avoid;
        }

        th {
            border-top: 1.2pt solid #000000;
            border-bottom: 0.8pt solid #000000;
            border-left: none;
            border-right: none;
            padding: 3px 2px;
            font-weight: bold;
            color: #000000;
            background-color: #ffffff;
            text-align: center;
        }

        td {
            border-top: none;
            border-bottom: 0.4pt solid #e0e0e0;
            border-left: none;
            border-right: none;
            padding: 2.2px 2px;
            text-align: center;
            color: #000000;
        }

        tr:last-child td {
            border-bottom: 1.2pt solid #000000;
        }

        /* =========================================================
           BLOCOS DE CÓDIGO E CAIXAS TÉCNICAS
           ========================================================= */
        pre, code {
            font-family: 'Consolas', 'Courier New', monospace;
            font-size: 8pt;
        }

        pre {
            background-color: #fbfbfb;
            border: 0.5pt solid #cccccc;
            border-radius: 0;
            padding: 6px 8px;
            overflow-x: auto;
            margin: 6px 0;
            line-height: 1.3;
        }

        code {
            background-color: #f5f5f5;
            padding: 1px 3px;
            border: 0.5pt solid #dddddd;
            border-radius: 0;
        }

        pre code {
            background-color: transparent;
            padding: 0;
            border: none;
        }

        /* =========================================================
           DIAGRAMA MERMAID SÓBRIO / MONOCROMÁTICO
           ========================================================= */
        .mermaid {
            display: flex;
            justify-content: center;
            margin: 12px 0;
            page-break-inside: avoid;
        }

        .mermaid svg {
            max-width: 90% !important;
            height: auto !important;
        }

        /* Equações KaTeX */
        .katex {
            font-size: 1.0em;
        }
        
        .katex-display {
            margin: 6px 0;
            text-align: center;
        }

        /* Listas */
        ul, ol {
            margin-top: 3px;
            margin-bottom: 6px;
            padding-left: 22px;
        }

        li {
            margin-bottom: 2px;
            text-align: justify;
        }
    </style>
</head>
<body>
    <div id="content"></div>

    # pyrefly: ignore [parse-error]
    <script id="raw-markdown" type="text/plain">""" + conteudo_md + """</script>

    <script>
        async function processarERenderizar() {
            const raw = document.getElementById('raw-markdown').textContent;
            
            // 1. Proteger blocos de math para que o Marked não altere caracteres especiais
            const mathTokens = [];
            let protegido = raw.replace(/\\$\\$([\\s\\S]*?)\\$\\$/g, function(match) {
                mathTokens.push(match);
                return '%%%MATHBLOCK' + (mathTokens.length - 1) + '%%%';
            });
            protegido = protegido.replace(/\\$([^\\$\\n]+?)\\$/g, function(match) {
                mathTokens.push(match);
                return '%%%MATHBLOCK' + (mathTokens.length - 1) + '%%%';
            });

            // 2. Parser do Markdown
            let html = marked.parse(protegido);

            // 3. Restaurar equações intactas
            html = html.replace(/%%%MATHBLOCK(\\d+)%%%/g, function(match, id) {
                return mathTokens[parseInt(id)];
            });

