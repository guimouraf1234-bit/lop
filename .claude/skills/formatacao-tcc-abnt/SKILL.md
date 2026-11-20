---
name: formatacao-tcc-abnt
description: Formata, revisa e corrige TCC, monografia, dissertação, tese ou relatório acadêmico segundo as normas ABNT (NBR 14724 trabalhos acadêmicos, NBR 6023 referências, NBR 10520 citações, NBR 6024 numeração progressiva, NBR 6027 sumário, NBR 6028 resumo, NBR 6034 índice, NBR 12225 lombada e normas tabulares do IBGE). Use sempre que o usuário pedir para formatar, normalizar, padronizar ou revisar um trabalho "nas normas ABNT", montar capa, folha de rosto, sumário, resumo/abstract, lista de figuras, legendas de figura e tabela, citações (autor-data, et al., apud) ou referências bibliográficas, ou verificar margens, fonte, espaçamento e paginação de um .docx, LaTeX/abnTeX2 ou Markdown — mesmo que ele não cite a ABNT explicitamente, mas fale em "TCC", "monografia", "relatório final", "normalização" ou "referência no formato certo".
---

# Formatação de TCC segundo a ABNT

Esta skill reúne as regras da ABNT para trabalhos acadêmicos e um verificador automático de `.docx`.
O objetivo é que o trabalho saia **consistente e defensável numa banca ou na biblioteca**: cada correção
proposta deve poder ser justificada pela norma (cite a norma e a subseção ao explicar).

## Normas e edições consideradas

| Assunto | Norma | Arquivo de referência |
|---|---|---|
| Estrutura e apresentação de trabalhos acadêmicos | NBR 14724:2024 (substituiu a 2011) | `references/estrutura.md`, `references/apresentacao-grafica.md` |
| Citações | NBR 10520:2023 (substituiu a 2002) | `references/citacoes.md` |
| Referências | NBR 6023:2018 (errata 2020) | `references/referencias.md` |
| Numeração progressiva, alíneas | NBR 6024:2012 | `references/apresentacao-grafica.md` |
| Sumário | NBR 6027:2012 | `references/apresentacao-grafica.md` |
| Resumo e palavras-chave | NBR 6028:2021 | `references/estrutura.md` |
| Índice | NBR 6034:2004 | `references/estrutura.md` |
| Lombada | NBR 12225:2004 | `references/estrutura.md` |
| Tabelas | IBGE, Normas de apresentação tabular, 3. ed., 1993 | `references/apresentacao-grafica.md` |

**Sobre as edições.** A NBR 14724 foi revisada em dez./2024; as regras centrais (A4, margens, fonte 12,
espaçamento 1,5, paginação) não mudaram em relação à 2011. A NBR 10520:2023 mudou bastante as citações (ver
`references/citacoes.md`). Há menções a uma NBR 6023 de 2025 que não puderam ser confirmadas: use a 2018 e
avise o usuário se ele disser que a instituição exige outra edição. **A instituição manda**: quando o guia da
biblioteca da universidade (no caso da UFMG, as *Diretrizes para normalização bibliográfica da produção
acadêmica da UFMG* e o modelo do colegiado do curso) divergir da ABNT, siga o guia e registre a divergência.
A ABNT não fixa família de fonte nem recuo de parágrafo — esses vêm do modelo institucional.

## Como trabalhar

1. **Descubra o formato e o escopo.** Word (`.docx`), LaTeX (abnTeX2) ou Markdown? O usuário quer o trabalho
   todo revisado ou só uma parte (referências, citações, legendas, capa)? Leia apenas os arquivos de
   referência necessários — eles são longos e organizados por assunto.
2. **Para `.docx`, rode o verificador primeiro** — ele encontra rapidamente o que é mecânico (margens, fonte,
   espaçamento, legendas, palavras-chave, citações em caixa alta, numeração de seções):
   ```bash
   python .claude/skills/formatacao-tcc-abnt/scripts/verificar_docx.py caminho/trabalho.docx
   python .claude/skills/formatacao-tcc-abnt/scripts/verificar_docx.py trabalho.docx --json saida.json
   ```
   O relatório separa **ERRO** (contraria regra explícita da norma), **AVISO** (provável problema, confira no
   contexto) e **INFO**. O script é heurístico: confirme cada achado lendo o trecho antes de corrigir, e
   lembre que ele não enxerga tudo (capa, folha de rosto, conteúdo do sumário, ordem dos elementos, notas de
   rodapé, numeração das páginas pré-textuais) — esses itens você verifica com o checklist.
3. **Revise o que o script não cobre** com `references/checklist.md` (lista completa, item a item, com a
   norma de cada regra).
4. **Corrija ou proponha correções.** Em `.docx`, prefira ajustar **estilos** (Normal, Título 1–5, Legenda,
   Citação longa, Referência) em vez de formatação direta parágrafo a parágrafo — assim o sumário automático e
   as listas de figuras/tabelas continuam funcionando. Para editar `.docx` use a skill `docx` se disponível.
   Para referências e citações, reescreva cada item no formato exato da norma (exemplos em
   `references/referencias.md` e `references/citacoes.md`).
5. **Relate** em português: o que foi corrigido, o que ficou pendente e o que depende de decisão do autor ou
   da instituição (ex.: família de fonte, impressão só no anverso ou anverso e verso, sistema autor-data ou
   numérico). Ao citar uma regra, indique a norma e a subseção, p. ex. "NBR 14724:2024, 5.3".

## Regras que mais aparecem em revisão (resumo rápido)

- **Papel e margens:** A4; anverso: superior e esquerda 3 cm, inferior e direita 2 cm (no verso, espelhado:
  superior e direita 3 cm, inferior e esquerda 2 cm).
- **Fonte:** tamanho 12 em todo o trabalho, inclusive capa; tamanho menor e uniforme (em geral 10) para
  citações longas, notas de rodapé, paginação, ficha catalográfica, legendas e fontes de ilustrações e tabelas.
  Texto em preto; cores só nas ilustrações.
- **Espaçamento:** 1,5 no texto; simples nas citações longas, notas, referências, legendas, natureza do
  trabalho e ficha catalográfica. Referências separadas entre si por uma linha em branco simples.
- **Paginação:** pré-textuais contados a partir da folha de rosto, mas sem número impresso; número em
  algarismo arábico a partir da primeira folha textual (Introdução), canto superior direito, 2 cm da borda
  superior e 2 cm da borda direita (no verso, canto superior esquerdo). Apêndices e anexos continuam a
  numeração.
- **Seções:** indicativo numérico arábico, alinhado à esquerda, separado do título por **um espaço** — sem
  ponto, hífen ou travessão (`1 INTRODUÇÃO`, `2.1 Modelo cinético`); até seção quinária; destaque gráfico
  decrescente e idêntico no sumário; seção primária começa em página nova (anverso, se impressão frente e
  verso). Títulos sem número centralizados: ERRATA, AGRADECIMENTOS, LISTA DE ..., RESUMO, ABSTRACT, SUMÁRIO,
  REFERÊNCIAS, GLOSSÁRIO, APÊNDICE, ANEXO, ÍNDICE. Sem título: folha de aprovação, dedicatória, epígrafe.
- **Ilustrações:** identificação **acima** — `Figura 3 – Título` (palavra designativa, número arábico,
  travessão, título); **abaixo**, `Fonte: ...` (obrigatória mesmo se for do autor: "Fonte: elaborado pelo
  autor (2026)."), depois legenda/notas. Citar no texto antes e inserir o mais perto possível.
- **Tabelas (IBGE):** título acima (`Tabela 2 – ...`), bordas horizontais em cima, embaixo e sob o cabeçalho,
  **sem bordas laterais**; fonte e notas abaixo. Quadro (conteúdo textual) é ilustração e é fechado.
- **Equações:** destacadas, centralizadas, numeradas à direita entre parênteses: `(1)`.
- **Siglas:** na primeira ocorrência, nome por extenso seguido da sigla entre parênteses.
- **Citações (NBR 10520:2023):** autor em maiúsculas e minúsculas **dentro e fora** dos parênteses —
  `(Silva, 2020, p. 15)` e `Silva (2020, p. 15)`; até 3 autores separados por `;` nos parênteses e por `e` no
  texto; 4 ou mais pode usar *et al.*; direta curta entre aspas; direta com mais de 3 linhas com recuo de 4 cm,
  fonte menor, espaço simples, sem aspas; página obrigatória na citação direta.
- **Referências (NBR 6023:2018):** alinhadas à esquerda (não justificar), espaço simples, uma linha em branco
  entre elas; `SOBRENOME, Prenome.`; destaque (negrito ou itálico, um só em todo o trabalho) no título da
  obra ou do periódico; `Disponível em: URL. Acesso em: 20 nov. 2026.`; DOI como `DOI: https://doi.org/...`.
- **Resumo (NBR 6028:2021):** parágrafo único, 150 a 500 palavras, voz ativa, 3ª pessoa;
  `Palavras-chave: lixiviação; zinco; modelagem híbrida.` (iniciais minúsculas, `;` entre elas, ponto final).

## LaTeX e Markdown

- **abnTeX2** (`\documentclass[12pt,a4paper,oneside,english,brazil]{abntex2}` + `abntex2cite` com
  `[alf]`) já implementa margens, paginação, sumário e títulos. Atenção: o `abntex2cite` segue a NBR
  10520:2002 e gera a chamada entre parênteses em CAIXA ALTA; para a NBR 10520:2023 use uma versão
  atualizada do `biblatex-abnt` ou ajuste o estilo, e **confira a saída compilada** antes de afirmar que está
  conforme. Use `\caption` acima e `\fonte{}` abaixo da figura/tabela.
- **Markdown/Pandoc**: use um modelo `.docx` de referência (`--reference-doc`) com os estilos ABNT e um CSL
  ABNT (`associacao-brasileira-de-normas-tecnicas.csl`) — confira se o CSL já segue a 10520:2023.

## Arquivos

- `references/estrutura.md` — elementos pré-textuais, textuais e pós-textuais: obrigatoriedade, ordem e
  conteúdo de cada um (capa, lombada, folha de rosto, ficha, errata, aprovação, resumo, listas, sumário,
  apêndices, anexos, glossário, índice).
- `references/apresentacao-grafica.md` — formatação: papel, margens, fonte, espaçamento, notas, paginação,
  numeração progressiva e alíneas, sumário, ilustrações, tabelas (IBGE), equações, siglas.
- `references/citacoes.md` — NBR 10520:2023 completa, com exemplos e as mudanças em relação à 2002.
- `references/referencias.md` — NBR 6023:2018: regras gerais e modelos por tipo de documento.
- `references/checklist.md` — checklist final de conferência, item a item.
- `scripts/verificar_docx.py` — verificador automático de `.docx` (requer `python-docx`).
