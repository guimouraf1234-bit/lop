const fs = require('fs');
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, AlignmentType,
        BorderStyle, ShadingType, PageBreak, VerticalAlign } = require('docx');

const FONT = 'Times New Roman';
const cm = c => Math.round(c * 567);
// texto com trechos [entre colchetes] destacados em amarelo (a preencher)
function runs(texto, opts = {}) {
  return texto.split(/(\[[^\]]*\])/).filter(s => s).map(s => new TextRun({
    text: s, font: FONT, size: opts.size || 24, bold: opts.bold, italics: opts.italics,
    highlight: s.startsWith('[') ? 'yellow' : undefined }));
}
const par = (texto, o = {}) => new Paragraph({
  children: runs(texto, o), alignment: o.align || AlignmentType.JUSTIFIED,
  spacing: { line: o.line || 360, before: o.before || 0, after: o.after ?? 0 },
  indent: o.indent ? { firstLine: cm(1.25) } : undefined, keepNext: o.keepNext });

const borda = { style: BorderStyle.SINGLE, size: 6, color: '000000' };
const nada = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
function tabela(cabec, linhas, larguras) {
  const total = larguras.reduce((a, b) => a + b, 0);
  const cel = (txt, i, header, ultima) => new TableCell({
    width: { size: larguras[i], type: WidthType.DXA }, verticalAlign: VerticalAlign.CENTER,
    margins: { top: 60, bottom: 60, left: 80, right: 80 },
    borders: { top: header ? borda : nada, bottom: (header || ultima) ? borda : nada, left: nada, right: nada },
    children: [new Paragraph({ alignment: header ? AlignmentType.CENTER : AlignmentType.LEFT,
      spacing: { line: 240 }, children: runs(txt, { size: 20, bold: header }) })] });
  return new Table({
    width: { size: total, type: WidthType.DXA }, columnWidths: larguras, alignment: AlignmentType.CENTER,
    rows: [new TableRow({ tableHeader: true, children: cabec.map((t, i) => cel(t, i, true)) }),
           ...linhas.map((l, k) => new TableRow({ cantSplit: true,
             children: l.map((t, i) => cel(t, i, false, k === linhas.length - 1)) }))] });
}
const titulo = t => new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 100, after: 60, line: 240 },
  keepNext: true, children: runs(t, { size: 20, bold: true }) });
const fonte = () => new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 40, after: 120 },
  children: runs('Fonte: Elaborado pelos autores (2026).', { size: 20, italics: true }) });
const secao = t => new Paragraph({ spacing: { before: 200, after: 80, line: 360 }, keepNext: true,
  children: runs(t, { bold: true }) });

const D1 = tabela(['Ferramenta', 'Desenvolvedor', 'Período de uso', 'Membros'], [
  ['Claude / Claude Code', 'Anthropic', 'set.–nov. 2026', 'Guilherme Moura de Souza Franco, [outros]'],
  ['[outras ferramentas]', '[ ]', '[ ]', '[ ]'],
], [2600, 2000, 1900, 2571]);

const D2 = tabela(['Etapa', 'Ferramenta', 'Motivação', 'O que a IA fez', 'O que os autores fizeram'], [
  ['Modelo white-box', 'Claude Code', 'Implementar o modelo de Bortot Coelho (2017)',
   'Gerou o código do balanço populacional e da cinética, e os testes',
   'Escolha das equações e dos parâmetros; comparação com a dissertação; verificação'],
  ['Híbrido paralelo', 'Claude Code', 'Implementar e comparar regressores para o resíduo',
   'Gerou o código do híbrido, da validação cruzada, do acoplador e das figuras',
   'Definição da arquitetura, das entradas, da partição e das restrições físicas; escolha do GP; análise dos resultados'],
  ['Ensaios novos', 'Claude Code', 'Calcular as conversões', 'Gerou a planilha de conversão',
   'Planejamento e execução dos ensaios; medições; correção dos rótulos E2/E3'],
  ['Figuras', 'Claude Code', 'Gerar os gráficos de resultados', 'Gerou os scripts e as figuras',
   'Definição do conteúdo e do layout; conferência dos valores'],
  ['Análise e verificação', 'Claude', 'Apoio à interpretação',
   'Explicou métodos, conferiu tabelas e comparou métodos numéricos', 'Decisões de critério e interpretação final'],
  ['Híbrido serial', '[ ]', '[ ]', '[ ]', '[ ]'],
  ['Black-box (MLP)', '[ ]', '[ ]', '[ ]', '[ ]'],
  ['Texto', '[ ]', '[ ]', '[ ]', 'Redação autoral e revisão crítica'],
], [1450, 1250, 1800, 2200, 2371]);

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 24 } } } },
  sections: [{
    properties: { page: { size: { width: cm(21), height: cm(29.7) },
      margin: { top: cm(3), left: cm(3), bottom: cm(2), right: cm(2) } } },
    children: [
      new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 240 },
        children: runs('APÊNDICE D – DECLARAÇÃO DE USO DE INTELIGÊNCIA ARTIFICIAL', { bold: true }) }),
      par('Este apêndice descreve o uso de ferramentas de Inteligência Artificial generativa no trabalho, conforme o Guia de Inteligência Artificial da UFMG (2026). O conteúdo final foi produzido com curadoria crítica dos autores, que assumem inteira responsabilidade pelos resultados apresentados.', { indent: true }),
      secao('D.1 Ferramentas utilizadas'),
      titulo('Tabela D.1 – Ferramentas de Inteligência Artificial utilizadas no trabalho'), D1, fonte(),
      secao('D.2 Uso por etapa'),
      titulo('Tabela D.2 – Uso de Inteligência Artificial em cada etapa do trabalho'), D2, fonte(),
      new Paragraph({ children: [new PageBreak()] }),
      new Paragraph({ spacing: { after: 240 }, children: runs('Parágrafo para a Metodologia (colar no início da Seção 4 e apagar esta página do apêndice)', { bold: true, italics: true }) }),
      par('Neste trabalho foram utilizadas ferramentas de Inteligência Artificial generativa como apoio à implementação computacional, à geração de figuras e à revisão do texto, conforme o Guia de Inteligência Artificial da UFMG (2026). As formulações físicas, a partição dos dados, os critérios de validação, os ensaios de laboratório e a interpretação dos resultados foram definidos pelos autores, que revisaram o código e os resultados produzidos. As ferramentas, as finalidades e a descrição do uso em cada etapa estão no Apêndice D. O conteúdo final é de inteira responsabilidade dos autores.', { indent: true }),
    ] }] });
Packer.toBuffer(doc).then(b => fs.writeFileSync('Apendice_D_uso_de_IA.docx', b));
