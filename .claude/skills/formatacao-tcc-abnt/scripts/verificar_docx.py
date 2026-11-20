#!/usr/bin/env python3
"""Verificador heurístico de formatação ABNT para trabalhos acadêmicos em .docx.

Verifica (NBR 14724:2024, 10520:2023, 6023:2018, 6024:2012, 6028:2021, IBGE):
página A4 e margens; fonte 12 no texto e menor nas exceções; entrelinha 1,5 / simples;
numeração e alinhamento de títulos; legendas de ilustrações e tabelas (posição, travessão,
"Fonte:", menção no texto); resumo e palavras-chave; chamadas de citação; lista de
referências; siglas; alíneas; paginação no cabeçalho; cor do texto.

Uso:
    python verificar_docx.py trabalho.docx [--json saida.json] [--max-exemplos 5]

Saída: relatório em Markdown no stdout com níveis ERRO / AVISO / INFO. É heurístico:
confirme cada achado no documento antes de corrigir.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field

try:
    import docx
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml.ns import qn
    from docx.shared import Length
    from docx.table import Table
    from docx.text.paragraph import Paragraph
except ImportError:  # pragma: no cover
    sys.exit("python-docx não instalado: pip install python-docx")

CM = 360000  # EMU por cm
TOL_CM = 0.1

TIPOS_ILUSTRACAO = (
    "Figura|Gráfico|Grafico|Quadro|Tabela|Fotografia|Esquema|Fluxograma|Mapa|Desenho|"
    "Organograma|Planta|Retrato|Imagem|Ilustração|Diagrama"
)
RE_LEGENDA = re.compile(rf"^\s*({TIPOS_ILUSTRACAO})\s+((?:[A-Z]\.)?\d+(?:[.\-]\d+)*)\s*(.?)\s*(.*)$")
RE_FONTE = re.compile(r"^\s*Fontes?\s*:", re.I)
RE_TITULO_NUM = re.compile(r"^(\d+(?:\.\d+)*)(\.?)(\s*)([-–—):]?)(\s*)(\S.*)$")

TITULOS_SEM_NUMERO = [
    "ERRATA", "AGRADECIMENTOS", "AGRADECIMENTO", "LISTA DE", "RESUMO", "ABSTRACT", "RESUMEN",
    "RÉSUMÉ", "SUMÁRIO", "REFERÊNCIAS", "GLOSSÁRIO", "APÊNDICE", "APÊNDICES", "ANEXO",
    "ANEXOS", "ÍNDICE",
]
FIM_REFERENCIAS = ("GLOSSÁRIO", "APÊNDICE", "ANEXO", "ÍNDICE")
ROMANOS = {"II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XII", "XIII", "XIV", "XV",
           "XVI", "XVII", "XVIII", "XIX", "XX", "XXI", "CD", "MD"}


@dataclass
class Achado:
    nivel: str          # ERRO, AVISO, INFO
    categoria: str
    mensagem: str
    norma: str = ""
    exemplos: list = field(default_factory=list)


# ---------------------------------------------------------------- utilidades de estilo
def cadeia(estilo):
    while estilo is not None:
        yield estilo
        estilo = estilo.base_style


def _doc_defaults(doc):
    dd = doc.styles.element.find(qn("w:docDefaults"))
    tam, linha = None, None
    if dd is not None:
        sz = dd.find(qn("w:rPrDefault") + "/" + qn("w:rPr") + "/" + qn("w:sz"))
        if sz is not None:
            tam = int(sz.get(qn("w:val"))) / 2
        sp = dd.find(qn("w:pPrDefault") + "/" + qn("w:pPr") + "/" + qn("w:spacing"))
        if sp is not None and sp.get(qn("w:line")):
            regra = sp.get(qn("w:lineRule"), "auto")
            v = int(sp.get(qn("w:line")))
            linha = v / 240 if regra == "auto" else ("pt", v / 20)
    return tam or 10.0, linha or 1.0


class Resolvedor:
    def __init__(self, doc):
        self.doc = doc
        self.tam_padrao, self.linha_padrao = _doc_defaults(doc)

    def tamanho_run(self, par, run):
        if run.font.size:
            return run.font.size.pt
        if run.style is not None:
            for st in cadeia(run.style):
                if st.font.size:
                    return st.font.size.pt
        for st in cadeia(par.style):
            if st.font.size:
                return st.font.size.pt
        return self.tam_padrao

    def tamanho(self, par):
        cont = Counter()
        for r in par.runs:
            n = len(r.text.strip())
            if n:
                cont[self.tamanho_run(par, r)] += n
        if not cont:
            return None
        return cont.most_common(1)[0][0]

    def tamanhos(self, par):
        return {self.tamanho_run(par, r) for r in par.runs if r.text.strip()}

    def fonte_run(self, par, run):
        if run.font.name:
            return run.font.name
        for st in cadeia(par.style):
            if st.font.name:
                return st.font.name
        return None

    def _pf(self, par, attr):
        v = getattr(par.paragraph_format, attr)
        if v is not None:
            return v
        for st in cadeia(par.style):
            v = getattr(st.paragraph_format, attr)
            if v is not None:
                return v
        return None

    def entrelinha(self, par):
        """Retorna ('mult', x) ou ('pt', y)."""
        v = self._pf(par, "line_spacing")
        if v is None:
            d = self.linha_padrao
            return d if isinstance(d, tuple) else ("mult", d)
        if isinstance(v, Length):
            return ("pt", v.pt)
        return ("mult", float(v))

    def entrelinha_mult(self, par):
        tipo, v = self.entrelinha(par)
        if tipo == "mult":
            return v
        tam = self.tamanho(par) or 12
        return v / (tam * 1.15)  # aproximação: 1 linha ≈ 1,15 × tamanho

    def alinhamento(self, par):
        return self._pf(par, "alignment")

    def recuo_esq_cm(self, par):
        v = self._pf(par, "left_indent")
        return (v or 0) / CM

    def espaco_depois_pt(self, par):
        v = self._pf(par, "space_after")
        return v.pt if v is not None else 0.0


# ---------------------------------------------------------------- utilidades de texto
def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def sem_acento(s):
    return "".join(c for c in unicodedata.normalize("NFD", s) if unicodedata.category(c) != "Mn")


def trecho(s, n=70):
    s = norm(s)
    return s if len(s) <= n else s[: n - 1] + "…"


def eh_titulo_estilo(par):
    nome = (par.style.name or "").lower() if par.style is not None else ""
    if nome.startswith(("heading", "título", "titulo", "title")):
        return True
    pprs = [par._element.pPr] + [st.element.pPr for st in cadeia(par.style)]
    for ppr in pprs:
        if ppr is not None:
            ol = ppr.find(qn("w:outlineLvl"))
            if ol is not None:
                return int(ol.get(qn("w:val"))) < 9
    return False


def titulo_sem_numero(texto):
    t = norm(texto).upper()
    return any(t == k or t.startswith(k + " ") for k in TITULOS_SEM_NUMERO)


def tem_imagem(par):
    el = par._element
    return bool(el.findall(".//" + qn("w:drawing")) or el.findall(".//" + qn("w:pict")))


def blocos(doc):
    """Sequência de ('p', Paragraph) e ('t', Table) na ordem do corpo."""
    for el in doc.element.body.iterchildren():
        if el.tag == qn("w:p"):
            yield "p", Paragraph(el, doc)
        elif el.tag == qn("w:tbl"):
            yield "t", Table(el, doc)


def numeracao_automatica(doc, par):
    """Retorna (lvlText, sufixo) da numeração automática do parágrafo, se houver."""
    numpr = None
    for fonte in [par._element.pPr] + [st.element.pPr for st in cadeia(par.style)]:
        if fonte is not None and fonte.numPr is not None:
            numpr = fonte.numPr
            break
    if numpr is None or numpr.numId is None:
        return None
    num_id = numpr.numId.val
    ilvl = numpr.ilvl.val if numpr.ilvl is not None else 0
    try:
        numbering = doc.part.numbering_part.element
    except (KeyError, NotImplementedError):
        return None
    num = numbering.find(f'{qn("w:num")}[@{qn("w:numId")}="{num_id}"]')
    if num is None:
        return None
    abs_id = num.find(qn("w:abstractNumId")).get(qn("w:val"))
    absn = numbering.find(f'{qn("w:abstractNum")}[@{qn("w:abstractNumId")}="{abs_id}"]')
    if absn is None:
        return None
    lvl = absn.find(f'{qn("w:lvl")}[@{qn("w:ilvl")}="{ilvl}"]')
    if lvl is None:
        return None
    txt = lvl.find(qn("w:lvlText"))
    suff = lvl.find(qn("w:suff"))
    fmt = lvl.find(qn("w:numFmt"))
    return (
        txt.get(qn("w:val")) if txt is not None else "",
        suff.get(qn("w:val")) if suff is not None else "tab",
        fmt.get(qn("w:val")) if fmt is not None else "decimal",
    )


# ---------------------------------------------------------------- verificador
class Verificador:
    def __init__(self, caminho, max_ex=5):
        self.doc = docx.Document(caminho)
        self.r = Resolvedor(self.doc)
        self.max_ex = max_ex
        self.achados: list[Achado] = []
        self.pars = self.doc.paragraphs
        self.idx = {p._element: i + 1 for i, p in enumerate(self.pars)}
        self._regioes()

    # --- registro
    def add(self, nivel, cat, msg, norma="", exemplos=None):
        self.achados.append(Achado(nivel, cat, msg, norma, (exemplos or [])[: self.max_ex]))

    def ref(self, par):
        return f"§{self.idx.get(par._element, '?')}: \"{trecho(par.text)}\""

    # --- regiões do documento
    def _regioes(self):
        self.titulos = [p for p in self.pars if p.text.strip() and
                        (eh_titulo_estilo(p) or titulo_sem_numero(p.text))]
        self.i_intro = None
        self.i_ref = None
        self.i_fim_ref = None
        for p in self.titulos:
            t = sem_acento(norm(p.text).upper())
            t_sem_num = re.sub(r"^\d+(\.\d+)*\.?\s*", "", t)
            i = self.idx[p._element]
            if self.i_intro is None and t_sem_num.startswith("INTRODUCAO"):
                self.i_intro = i
            if t_sem_num.startswith(("REFERENCIAS", "BIBLIOGRAFIA")):
                self.i_ref = i
            elif self.i_ref and i > self.i_ref and self.i_fim_ref is None and \
                    t_sem_num.startswith(tuple(sem_acento(x) for x in FIM_REFERENCIAS)):
                self.i_fim_ref = i
        self.ini_texto = self.i_intro or 1
        self.fim_texto = self.i_ref or len(self.pars) + 1

    def regiao(self, par):
        i = self.idx[par._element]
        if self.i_ref and i > self.i_ref and (self.i_fim_ref is None or i < self.i_fim_ref):
            return "referencias"
        if i < self.ini_texto:
            return "pre"
        if i >= self.fim_texto:
            return "pos"
        return "texto"

    def tipo(self, par):
        t = par.text
        if not t.strip():
            return "vazio"
        if eh_titulo_estilo(par) or titulo_sem_numero(t) or (t.isupper() and len(t) < 200):
            return "titulo"
        if RE_LEGENDA.match(t) and len(t) < 300:
            return "legenda"
        if RE_FONTE.match(t) or re.match(r"^\s*Notas?\s*:", t):
            return "fonte"
        if self.regiao(par) == "referencias":
            return "referencia"
        if self.r.recuo_esq_cm(par) >= 3.5:
            return "citacao_longa"
        return "texto"

    # --- verificações
    def pagina(self):
        espelho = self.doc.settings.element.find(qn("w:mirrorMargins")) is not None
        for k, s in enumerate(self.doc.sections, 1):
            w, h = (s.page_width or 0) / CM, (s.page_height or 0) / CM
            if abs(w - 21.0) > TOL_CM or abs(h - 29.7) > TOL_CM:
                self.add("ERRO", "Página", f"Seção {k}: papel {w:.2f} × {h:.2f} cm; deve ser A4 "
                         "(21,0 × 29,7 cm).", "NBR 14724, 5.1")
            m = {"superior": s.top_margin, "inferior": s.bottom_margin,
                 "esquerda": s.left_margin, "direita": s.right_margin}
            m = {k2: (v or 0) / CM for k2, v in m.items()}
            alvo = {"superior": 3.0, "inferior": 2.0, "esquerda": 3.0, "direita": 2.0}
            errados = [f"{k2} {m[k2]:.2f} cm (esperado {alvo[k2]:.0f})" for k2 in alvo
                       if abs(m[k2] - alvo[k2]) > TOL_CM]
            if errados:
                self.add("ERRO", "Margens", f"Seção {k}: " + "; ".join(errados) +
                         (" (margens espelhadas: esquerda = interna, direita = externa)" if espelho else ""),
                         "NBR 14724, 5.1")
            if s.gutter and s.gutter / CM > TOL_CM:
                self.add("AVISO", "Margens", f"Seção {k}: medianiz de {s.gutter / CM:.2f} cm somada à "
                         "margem; a ABNT prevê só 3 cm/2 cm.", "NBR 14724, 5.1")
        if espelho:
            self.add("INFO", "Margens", "Margens espelhadas ativadas (impressão frente e verso).")

    def fontes_e_espacamento(self):
        tam_errado, linha_errada, menores, cores = [], [], defaultdict(list), []
        familias = Counter()
        tam_titulo = []
        for p in self.pars:
            tp = self.tipo(p)
            if tp == "vazio":
                continue
            reg = self.regiao(p)
            tam = self.r.tamanho(p)
            for run in p.runs:
                if run.text.strip():
                    f = self.r.fonte_run(p, run)
                    if f:
                        familias[f] += len(run.text)
                    c = run.font.color
                    if c is not None and c.type is not None and c.rgb is not None and \
                            str(c.rgb) not in ("000000", "auto") and reg != "pre":
                        cores.append(f"{self.ref(p)} cor #{c.rgb}")
            if tp == "texto" and reg in ("texto", "pos", "pre") and tam and abs(tam - 12) > 0.1:
                tam_errado.append(f"{self.ref(p)} → {tam:g} pt")
            if tp == "titulo" and tam and abs(tam - 12) > 0.1:
                tam_titulo.append(f"{self.ref(p)} → {tam:g} pt")
            if tp in ("legenda", "fonte", "citacao_longa") and tam:
                menores[tam].append(self.ref(p))
            if tp == "texto" and reg == "texto":
                m = self.r.entrelinha_mult(p)
                if abs(m - 1.5) > 0.06:
                    linha_errada.append(f"{self.ref(p)} → {m:.2f}")
        if tam_errado:
            self.add("ERRO", "Fonte", f"{len(tam_errado)} parágrafo(s) de texto com tamanho ≠ 12 pt.",
                     "NBR 14724, 5.1", tam_errado)
        if tam_titulo:
            self.add("AVISO", "Fonte", f"{len(tam_titulo)} título(s) com tamanho ≠ 12 pt (a norma pede 12 "
                     "em todo o trabalho; destaque por caixa/negrito).", "NBR 14724, 5.1", tam_titulo)
        if 12.0 in menores:
            self.add("ERRO", "Fonte", f"{len(menores[12.0])} legenda(s), fonte(s) de ilustração ou citação(ões) "
                     "longa(s) em 12 pt; devem estar em tamanho menor.", "NBR 14724, 5.1", menores[12.0])
        tams_menores = sorted(t for t in menores if t < 12)
        if len(tams_menores) > 1:
            self.add("AVISO", "Fonte", "Tamanho menor não uniforme entre legendas/fontes/citações longas: "
                     + ", ".join(f"{t:g} pt ({len(menores[t])}×)" for t in tams_menores) + ".",
                     "NBR 14724, 5.1")
        if linha_errada:
            self.add("ERRO", "Espaçamento", f"{len(linha_errada)} parágrafo(s) do texto com entrelinha ≠ 1,5.",
                     "NBR 14724, 5.2", linha_errada)
        for tp_nome, rotulo in (("citacao_longa", "citação longa"), ("legenda", "legenda"),
                                ("fonte", "fonte/nota de ilustração")):
            ruins = [f"{self.ref(p)} → {self.r.entrelinha_mult(p):.2f}" for p in self.pars
                     if self.tipo(p) == tp_nome and self.r.entrelinha_mult(p) > 1.15]
            if ruins:
                self.add("ERRO", "Espaçamento", f"{len(ruins)} {rotulo}(s) sem espaço simples.",
                         "NBR 14724, 5.2", ruins)
        if familias:
            tot = sum(familias.values())
            principais = [(f, n) for f, n in familias.most_common() if n / tot > 0.02]
            nivel = "AVISO" if len(principais) > 2 else "INFO"
            self.add(nivel, "Fonte", "Famílias de fonte (por caracteres): " +
                     ", ".join(f"{f} {100 * n / tot:.0f}%" for f, n in familias.most_common(5)) +
                     ". A ABNT não fixa a família, mas use uma só.")
        if cores:
            self.add("AVISO", "Cor", f"{len(cores)} trecho(s) de texto em cor diferente de preto "
                     "(cores só em ilustrações).", "NBR 14724, 5.1", cores)

    def titulos_secoes(self):
        ponto, separador, profundo, nao_centralizado, tabs = [], [], [], [], []
        prev_nivel = 0
        for p in self.titulos:
            t = norm(p.text)
            auto = numeracao_automatica(self.doc, p)
            if titulo_sem_numero(t):
                if self.r.alinhamento(p) != WD_ALIGN_PARAGRAPH.CENTER:
                    nao_centralizado.append(self.ref(p))
                if sem_acento(t.upper()) != "REFERENCIAS" and sem_acento(t.upper()).startswith(
                        ("REFERENCIAS", "BIBLIOGRAFIA")):
                    self.add("AVISO", "Títulos", f"{self.ref(p)}: o título normativo é \"REFERÊNCIAS\".",
                             "NBR 14724, 4.3.1")
                if auto and auto[2] != "none":
                    self.add("ERRO", "Títulos", f"{self.ref(p)}: título sem indicativo numérico está "
                             f"numerado automaticamente ({auto[0]}).", "NBR 6024:2012")
                continue
            m = RE_TITULO_NUM.match(t)
            if m:
                num, pto, esp, sinal = m.group(1), m.group(2), m.group(3), m.group(4)
                nivel = num.count(".") + 1
                if pto:
                    ponto.append(self.ref(p))
                if sinal or len(esp) != 1:
                    separador.append(self.ref(p))
                if nivel > 5:
                    profundo.append(self.ref(p))
                if prev_nivel and nivel > prev_nivel + 1:
                    self.add("AVISO", "Títulos", f"{self.ref(p)}: salto de nível (de {prev_nivel} para "
                             f"{nivel}).", "NBR 6024")
                prev_nivel = nivel
            elif auto:
                lvl_txt, suff, fmt = auto
                if fmt != "decimal":
                    self.add("AVISO", "Títulos", f"{self.ref(p)}: numeração automática do tipo '{fmt}'; "
                             "use algarismos arábicos.", "NBR 6024, 4.1")
                if lvl_txt.endswith("."):
                    ponto.append(f"{self.ref(p)} (numeração automática \"{lvl_txt}\")")
                if suff == "tab":
                    tabs.append(self.ref(p))
                prev_nivel = lvl_txt.count("%")
        if ponto:
            self.add("ERRO", "Títulos", "Indicativo de seção seguido de ponto (\"1.\" / \"2.1.\").",
                     "NBR 6024, 4.1", ponto)
        if separador:
            self.add("ERRO", "Títulos", "Indicativo separado do título por sinal (hífen, travessão, "
                     "parêntese, dois-pontos) ou por mais de um espaço; use um único espaço.",
                     "NBR 6024, 4.1", separador)
        if profundo:
            self.add("ERRO", "Títulos", "Seção além da quinária; use alíneas.", "NBR 6024, 4.1", profundo)
        if nao_centralizado:
            self.add("ERRO", "Títulos", "Títulos sem indicativo numérico devem ser centralizados.",
                     "NBR 14724, 5.2.3", nao_centralizado)
        if tabs:
            self.add("AVISO", "Títulos", f"{len(tabs)} título(s) com numeração automática seguida de "
                     "tabulação; a norma prevê um espaço (Definir novo nível → seguir número com: Espaço).",
                     "NBR 6024, 4.1", tabs)

    def legendas(self):
        seq = list(blocos(self.doc))
        por_tipo = defaultdict(list)
        sem_fonte, abaixo, separador, ponto_final, nao_citada = [], [], [], [], []
        texto_corpo = " ".join(p.text for p in self.pars if self.tipo(p) == "texto")
        for k, (kind, obj) in enumerate(seq):
            if kind != "p":
                continue
            m = RE_LEGENDA.match(obj.text)
            if not m or len(obj.text) >= 300 or self.tipo(obj) != "legenda":
                continue
            tipo, num, sep, titulo = m.groups()
            por_tipo[tipo].append((num, obj))
            if sep != "–" or not titulo:
                separador.append(self.ref(obj))
            if titulo.rstrip().endswith("."):
                ponto_final.append(self.ref(obj))
            # posição: o que vem antes/depois
            ant = seq[k - 1] if k > 0 else None
            prox = seq[k + 1] if k + 1 < len(seq) else None
            ant_ilustr = ant and ((ant[0] == "t") or (ant[0] == "p" and tem_imagem(ant[1])))
            prox_ilustr = prox and ((prox[0] == "t") or (prox[0] == "p" and tem_imagem(prox[1])))
            if ant_ilustr and not prox_ilustr:
                abaixo.append(self.ref(obj))
            # fonte depois da ilustração, antes da próxima legenda
            achou = False
            for kind2, obj2 in seq[k + 1: k + 12]:
                if kind2 == "p":
                    if RE_FONTE.match(obj2.text):
                        achou = True
                        break
                    if RE_LEGENDA.match(obj2.text) and self.tipo(obj2) == "legenda":
                        break
            if not achou:
                sem_fonte.append(self.ref(obj))
            if not re.search(rf"\b{tipo}\s+{re.escape(num)}\b", texto_corpo):
                nao_citada.append(self.ref(obj))
        if separador:
            self.add("ERRO", "Ilustrações", "Identificação fora do padrão \"Figura N – Título\" (use "
                     "travessão –, não hífen, ponto ou dois-pontos).", "NBR 14724, 5.8", separador)
        if ponto_final:
            self.add("AVISO", "Ilustrações", "Título de ilustração/tabela terminando em ponto.",
                     "NBR 14724, 5.8 (convenção)", ponto_final)
        if abaixo:
            self.add("ERRO", "Ilustrações", "Identificação abaixo da ilustração/tabela; deve ficar na "
                     "parte superior.", "NBR 14724, 5.8 / IBGE", abaixo)
        if sem_fonte:
            self.add("ERRO", "Ilustrações", "Ilustração/tabela sem \"Fonte:\" abaixo (obrigatória, mesmo "
                     "que seja do autor: \"Fonte: elaborado pelo autor (ano).\").", "NBR 14724, 5.8",
                     sem_fonte)
        if nao_citada:
            self.add("AVISO", "Ilustrações", "Ilustração/tabela não mencionada no texto (\"Figura N\").",
                     "NBR 14724, 5.8", nao_citada)
        for tipo, itens in por_tipo.items():
            nums = [n for n, _ in itens]
            simples = [int(n) for n in nums if n.isdigit()]
            if simples and simples != list(range(simples[0], simples[0] + len(simples))):
                self.add("AVISO", "Ilustrações", f"Numeração de {tipo} fora de sequência: "
                         + ", ".join(nums[:20]) + ".", "NBR 14724, 5.8")
            dup = [n for n, c in Counter(nums).items() if c > 1]
            if dup:
                self.add("ERRO", "Ilustrações", f"{tipo} com número repetido: {', '.join(dup)}.")

    def resumo(self):
        for rotulo, kw_re, nome_kw in (
            ("RESUMO", r"^\s*Palavras?[- ]chave", "Palavras-chave:"),
            ("ABSTRACT", r"^\s*Key[- ]?words", "Keywords:"),
        ):
            ini = next((self.idx[p._element] for p in self.titulos if norm(p.text).upper() == rotulo), None)
            kw = None
            for p in self.pars:
                if re.match(kw_re, p.text, re.I) and (ini is None or self.idx[p._element] > ini):
                    kw = p
                    break
            if ini:
                fim = self.idx[kw._element] if kw else None
                corpo = [p for p in self.pars[ini: (fim - 1) if fim else ini + 6] if p.text.strip()
                         and not eh_titulo_estilo(p)]
                if kw is None:
                    corpo = corpo[:1]
                if len(corpo) > 1:
                    self.add("ERRO", "Resumo", f"{rotulo}: {len(corpo)} parágrafos; deve ser parágrafo único.",
                             "NBR 6028:2021", [self.ref(p) for p in corpo])
                n = sum(len(p.text.split()) for p in corpo)
                if corpo and not 150 <= n <= 500:
                    self.add("AVISO", "Resumo", f"{rotulo} com {n} palavras; trabalhos acadêmicos: 150 a 500.",
                             "NBR 6028:2021")
                elif corpo:
                    self.add("INFO", "Resumo", f"{rotulo} com {n} palavras.")
            if kw is None:
                if ini:
                    self.add("ERRO", "Resumo", f"{rotulo} sem a linha \"{nome_kw}\" logo abaixo.", "NBR 6028:2021")
                continue
            t = norm(kw.text)
            prefixo, _, termos = t.partition(":")
            if not _:
                self.add("ERRO", "Resumo", f"{self.ref(kw)}: falta dois-pontos após \"{nome_kw[:-1]}\".",
                         "NBR 6028:2021")
                continue
            if rotulo == "RESUMO" and prefixo.strip() != "Palavras-chave":
                self.add("AVISO", "Resumo", f"{self.ref(kw)}: grafe \"Palavras-chave:\".", "NBR 6028:2021")
            termos = termos.strip()
            if "," in termos and ";" not in termos:
                self.add("ERRO", "Resumo", f"{self.ref(kw)}: termos separados por vírgula; use ponto e vírgula.",
                         "NBR 6028:2021")
            elif re.search(r"\.\s+\S", termos.rstrip(".")) and ";" not in termos:
                self.add("ERRO", "Resumo", f"{self.ref(kw)}: termos separados por ponto; use ponto e vírgula.",
                         "NBR 6028:2021")
            if not termos.endswith("."):
                self.add("ERRO", "Resumo", f"{self.ref(kw)}: palavras-chave devem terminar com ponto final.",
                         "NBR 6028:2021")
            if rotulo == "RESUMO":
                maiusc = [x.strip() for x in termos.rstrip(".").split(";")
                          if x.strip() and x.strip()[0].isupper() and not x.strip().isupper()]
                if maiusc:
                    self.add("AVISO", "Resumo", f"{self.ref(kw)}: termos com inicial maiúscula ("
                             + "; ".join(maiusc[:6]) + "). Use minúsculas, exceto nomes próprios e "
                             "científicos.", "NBR 6028:2021")

    def citacoes(self):
        caixa_alta, e_dentro, etal, pagina, latinas, p_sem_espaco, sem_pagina = [], [], [], [], [], [], []
        siglas_def = set(re.findall(r"\(([A-Z][A-Z0-9]{1,9})\)", " ".join(p.text for p in self.pars)))
        for p in self.pars:
            if self.tipo(p) not in ("texto", "citacao_longa", "fonte", "legenda"):
                continue
            t = p.text
            for m in re.finditer(r"\(([^()]*?\b(?:1[5-9]|20)\d{2}[a-z]?\b[^()]*)\)", t):
                for parte in m.group(1).split(";"):
                    autor = parte.split(",")[0].strip()
                    autor_limpo = re.sub(r"\bet al\.?", "", autor).strip()
                    palavras = [w for w in re.findall(r"[A-Za-zÀ-ÿ'\-]+", autor_limpo)
                                if w not in ("e", "and", "de", "da", "do", "dos", "das")]
                    if not palavras or not re.search(r"\d{4}", parte):
                        continue
                    if all(w.isupper() for w in palavras) and not \
                            (len(palavras) == 1 and (palavras[0] in siglas_def or len(palavras[0]) <= 4)):
                        caixa_alta.append(f"{self.ref(p)} → ({trecho(m.group(1), 50)})")
                    if re.search(r"\s(e|and|&)\s", autor_limpo) and "," in parte:
                        e_dentro.append(f"{self.ref(p)} → ({trecho(m.group(1), 50)})")
            if re.search(r"\bet\.?\s?al(?!\.)\b|\bet\.\s?al\.", t):
                etal.append(self.ref(p))
            if re.search(r"\b(pg|pág|pag|pp)\.\s*\d", t, re.I):
                pagina.append(self.ref(p))
            if re.search(r"\bp\.\d", t):
                p_sem_espaco.append(self.ref(p))
            if re.search(r"\b(op\. cit|ibid|ibidem|idem|loc\. cit)\b", t, re.I):
                latinas.append(self.ref(p))
            for q in re.finditer(r"[\"“]([^\"”]{25,})[\"”]", t):
                janela = t[max(0, q.start() - 120): q.end() + 60]
                if re.search(r"\b(1[5-9]|20)\d{2}[a-z]?\b", janela) and not \
                        re.search(r"\bp\.\s?\d|não paginad|f\.\s?\d", janela):
                    sem_pagina.append(f"{self.ref(p)} → \"{trecho(q.group(1), 40)}\"")
        if caixa_alta:
            self.add("AVISO", "Citações", "Autor entre parênteses em CAIXA ALTA (forma da NBR 10520:2002). "
                     "A NBR 10520:2023 usa maiúsculas e minúsculas: (Silva, 2020). Siglas ficam em "
                     "maiúsculas. Confira o guia da instituição.", "NBR 10520:2023", caixa_alta)
        if e_dentro:
            self.add("ERRO", "Citações", "Autores ligados por \"e\" dentro dos parênteses; use ponto e "
                     "vírgula: (Silva; Souza, 2020).", "NBR 10520:2023", e_dentro)
        if etal:
            self.add("ERRO", "Citações", "\"et al.\" grafado errado (sem ponto final ou com \"et.\").",
                     "NBR 10520:2023", etal)
        if pagina:
            self.add("ERRO", "Citações", "Página indicada como pg./pág./pp.; use \"p. 15\" ou \"p. 15-20\".",
                     "NBR 10520:2023", pagina)
        if p_sem_espaco:
            self.add("AVISO", "Citações", "\"p.15\" sem espaço; use \"p. 15\".", "NBR 10520:2023", p_sem_espaco)
        if latinas:
            self.add("AVISO", "Citações", "Expressões latinas (op. cit., ibid., idem, loc. cit.) no texto; "
                     "só são permitidas em notas de rodapé (no texto, apenas apud).", "NBR 10520:2023",
                     latinas)
        if sem_pagina:
            self.add("AVISO", "Citações", "Citação direta (entre aspas) sem indicação de página próxima.",
                     "NBR 10520:2023", sem_pagina)

    def referencias(self):
        if not self.i_ref:
            self.add("AVISO", "Referências", "Título \"REFERÊNCIAS\" não encontrado (o verificador procura um "
                     "parágrafo de título com esse texto).", "NBR 14724, 4.3.1")
            return
        refs = [p for p in self.pars if self.regiao(p) == "referencias" and self.tipo(p) == "referencia"]
        if not refs:
            return
        justif = [self.ref(p) for p in refs if self.r.alinhamento(p) == WD_ALIGN_PARAGRAPH.JUSTIFY]
        linha = [self.ref(p) for p in refs if self.r.entrelinha_mult(p) > 1.15]
        recuo = [self.ref(p) for p in refs if (self.r._pf(p, "first_line_indent") or 0) != 0]
        if justif:
            self.add("ERRO", "Referências", f"{len(justif)} referência(s) justificada(s); alinhe à esquerda.",
                     "NBR 6023, 6.3 / NBR 14724, 5.2", justif)
        if linha:
            self.add("ERRO", "Referências", f"{len(linha)} referência(s) sem espaço simples.",
                     "NBR 14724, 5.2", linha)
        if recuo:
            self.add("AVISO", "Referências", f"{len(recuo)} referência(s) com recuo (de primeira linha ou "
                     "pendente); a norma alinha à margem esquerda.", "NBR 6023, 6.3", recuo)
        # separação por linha em branco
        sem_sep = []
        for p in refs:
            i = self.idx[p._element]
            prox = self.pars[i] if i < len(self.pars) else None
            if prox is not None and prox.text.strip() and self.r.espaco_depois_pt(p) < 8:
                sem_sep.append(self.ref(p))
        if sem_sep:
            self.add("AVISO", "Referências", f"{len(sem_sep)} referência(s) sem linha em branco (ou espaço "
                     "depois ≈ 12 pt) antes da seguinte.", "NBR 14724, 5.2", sem_sep)
        # ordem alfabética
        chaves = []
        for p in refs:
            t = norm(p.text)
            if t.startswith("_"):
                continue
            chaves.append((sem_acento(t).upper(), p))
        fora = [f"{self.ref(b)} vem depois de {self.ref(a)}" for (ka, a), (kb, b) in zip(chaves, chaves[1:])
                if re.sub(r"[^A-Z ]", "", kb)[:20] < re.sub(r"[^A-Z ]", "", ka)[:20]]
        numerada = sum(bool(re.match(r"^\[?\d+[\].]?\s", norm(p.text))) for p in refs) > len(refs) / 2
        if fora and not numerada:
            self.add("AVISO", "Referências", f"{len(fora)} par(es) fora da ordem alfabética.",
                     "NBR 6023:2018 (ordenação)", fora)
        padroes = [
            (r"<\s*https?://", "ERRO", "URL entre < >; a NBR 6023:2018 não usa os sinais."),
            (r"Acesso em:\s*\d{1,2}/\d{1,2}/\d{2,4}", "ERRO", "Data de acesso numérica; use \"Acesso em: 20 nov. 2026.\""),
            (r"\[s\.\s?l\.?\]|\[S\.l\.\]|\[s\.n\.\]|\[S\.\s?n\.\]", "ERRO", "Grafe [S. l.] e [s. n.] (com espaço; S maiúsculo em S. l.)."),
            (r"\((?:Org|Orgs|Coord|Ed|Eds|orgs|coords)\.\)", "ERRO", "Use (org.), (coord.), (ed.) — minúsculas, singular."),
            (r"\bpp\.\s*\d", "ERRO", "Use \"p. 10-20\", não \"pp.\"."),
            (r"\bet\.?\s?al(?!\.)\b|\bet\.\s?al\.", "ERRO", "\"et al.\" grafado errado."),
            (r"Acessado em|Acesso:|Disponível:", "ERRO", "Use exatamente \"Disponível em:\" e \"Acesso em:\"."),
        ]
        for rx, nivel, msg in padroes:
            ex = [self.ref(p) for p in refs if re.search(rx, p.text)]
            if ex:
                self.add(nivel, "Referências", f"{msg} ({len(ex)} ocorrência(s))", "NBR 6023:2018", ex)
        sem_acesso = [self.ref(p) for p in refs if "Disponível em" in p.text and "Acesso em" not in p.text]
        if sem_acesso:
            self.add("ERRO", "Referências", "\"Disponível em:\" sem \"Acesso em:\".", "NBR 6023:2018", sem_acesso)
        url_solta = [self.ref(p) for p in refs if re.search(r"https?://", p.text)
                     and "Disponível em" not in p.text and "DOI" not in p.text]
        if url_solta:
            self.add("AVISO", "Referências", "URL sem \"Disponível em:\" (ou \"DOI:\").", "NBR 6023:2018",
                     url_solta)
        self.add("INFO", "Referências", f"{len(refs)} referência(s) encontradas.")

    def siglas(self):
        vistos = {}
        definidas = set()
        for p in self.pars:
            if self.tipo(p) not in ("texto", "citacao_longa") or self.regiao(p) != "texto":
                continue
            t = p.text
            if t.isupper() or re.match(r"^[A-ZÀ-Ý'\-]{2,}(?: [A-ZÀ-Ý'\-]+)*,\s", t):
                continue  # títulos e parágrafos com cara de referência
            # remove chamadas de citação (autor-data), cujos sobrenomes podem estar em caixa alta
            t = re.sub(r"\([^()]*\b(?:1[5-9]|20)\d{2}[a-z]?\b[^()]*\)", lambda m: " " * len(m.group(0)), t)
            for m in re.finditer(r"\b([A-Z]{2,8})s?\b", t):
                s = m.group(1)
                if s in ROMANOS or s in vistos or s in definidas:
                    continue
                depois = t[m.end(): m.end() + 1]
                abre = t.rfind("(", 0, m.start())
                dentro = abre >= 0 and ")" not in t[abre: m.start()]
                if dentro and depois == ")":  # "Nome (SIGLA)" ou "(Expansão, SIGLA)"
                    definidas.add(s)
                else:
                    vistos[s] = p
        nao_def = [f"{s} — 1ª ocorrência {self.ref(p)}" for s, p in vistos.items() if s not in definidas]
        if nao_def:
            self.add("AVISO", "Siglas", f"{len(nao_def)} sigla(s) usada(s) antes de serem definidas como "
                     "\"Nome por extenso (SIGLA)\" (ou nunca definidas no texto).", "NBR 14724, 5.6", nao_def)

    def alineas(self):
        grupos, atual = [], []
        for p in self.pars:
            if re.match(r"^\s*[a-z]{1,2}\)\s", p.text):
                atual.append(p)
            elif re.match(r"^\s*[–—-]\s", p.text) and atual:
                continue  # subalínea
            else:
                if atual:
                    grupos.append(atual)
                atual = []
        if atual:
            grupos.append(atual)
        ruins = []
        for g in grupos:
            for k, p in enumerate(g):
                fim = p.text.rstrip()
                ultima = k == len(g) - 1
                if ultima and not fim.endswith("."):
                    ruins.append(f"{self.ref(p)} (última alínea deve terminar em ponto)")
                elif not ultima and not fim.endswith((";", ":")):
                    ruins.append(f"{self.ref(p)} (alínea deve terminar em ponto e vírgula)")
                if re.match(r"^\s*[a-z]{1,2}\)\s+[A-ZÀ-Ý](?:[a-zà-ÿ]|\s)", p.text):
                    ruins.append(f"{self.ref(p)} (texto da alínea deve começar com minúscula)")
            i0 = self.idx[g[0]._element]
            ant = next((q for q in reversed(self.pars[: i0 - 1]) if q.text.strip()), None)
            if ant is not None and not ant.text.rstrip().endswith(":"):
                ruins.append(f"{self.ref(ant)} (texto antes das alíneas deve terminar em dois-pontos)")
        if ruins:
            self.add("AVISO", "Alíneas", f"{len(ruins)} problema(s) de pontuação em alíneas.", "NBR 6024, 4.2",
                     ruins)

    def paginacao(self):
        def tem_page(el):
            if el is None:
                return False
            for fs in el.iter(qn("w:fldSimple")):
                if "PAGE" in (fs.get(qn("w:instr")) or "").upper().split():
                    return True
            for it in el.iter(qn("w:instrText")):
                if (it.text or "").strip().upper().startswith("PAGE"):
                    return True
            return False

        info = []
        for k, s in enumerate(self.doc.sections, 1):
            cab = None if s.header.is_linked_to_previous and k > 1 else s.header
            rod = None if s.footer.is_linked_to_previous and k > 1 else s.footer
            if cab is not None and tem_page(cab._element):
                pars = [p for p in cab.paragraphs if tem_page(p._element)]
                if pars and self.r.alinhamento(pars[0]) != WD_ALIGN_PARAGRAPH.RIGHT:
                    self.add("AVISO", "Paginação", f"Seção {k}: número de página no cabeçalho não está "
                             "alinhado à direita (se usar tabulação, confira a posição).", "NBR 14724, 5.3")
                info.append(f"seção {k}: número no cabeçalho")
            if rod is not None and tem_page(rod._element):
                self.add("ERRO", "Paginação", f"Seção {k}: número de página no rodapé; deve ficar no canto "
                         "superior direito.", "NBR 14724, 5.3")
            hd = (s.header_distance or 0) / CM
            if hd and abs(hd - 2.0) > 0.15 and cab is not None and tem_page(cab._element):
                self.add("AVISO", "Paginação", f"Seção {k}: cabeçalho a {hd:.2f} cm da borda; o número "
                         "deve ficar a 2 cm da borda superior.", "NBR 14724, 5.3")
            pg = s._sectPr.find(qn("w:pgNumType"))
            if pg is not None and pg.get(qn("w:start")) == "1" and k > 1:
                self.add("AVISO", "Paginação", f"Seção {k}: numeração reiniciada em 1. A contagem começa "
                         "na folha de rosto e não deve ser reiniciada na Introdução.", "NBR 14724, 5.3")
        if not info and not any(a.categoria == "Paginação" for a in self.achados):
            self.add("AVISO", "Paginação", "Nenhum campo de número de página encontrado nos cabeçalhos.",
                     "NBR 14724, 5.3")
        if self.i_intro and self.doc.sections:
            primeira = self.doc.sections[0]
            if tem_page(primeira.header._element) and len(self.doc.sections) == 1:
                self.add("AVISO", "Paginação", "Documento com uma só seção e número de página: os elementos "
                         "pré-textuais ficarão numerados. Insira quebra de seção antes da Introdução.",
                         "NBR 14724, 5.3")

    def executar(self):
        self.pagina()
        self.fontes_e_espacamento()
        self.titulos_secoes()
        self.legendas()
        self.resumo()
        self.citacoes()
        self.referencias()
        self.siglas()
        self.alineas()
        self.paginacao()
        regiao = (f"Introdução no §{self.i_intro}" if self.i_intro else
                  "Introdução não encontrada — todo o documento tratado como texto")
        regiao += f"; Referências no §{self.i_ref}" if self.i_ref else "; Referências não encontradas"
        self.add("INFO", "Documento", f"{len(self.pars)} parágrafos, {len(self.doc.tables)} tabelas, "
                 f"{len(self.doc.sections)} seção(ões). {regiao}.")
        return self.achados


def relatorio_md(achados, caminho):
    ordem = {"ERRO": 0, "AVISO": 1, "INFO": 2}
    cont = Counter(a.nivel for a in achados)
    linhas = [f"# Verificação ABNT — {caminho}", "",
              f"**{cont['ERRO']} erro(s), {cont['AVISO']} aviso(s), {cont['INFO']} informação(ões).** "
              "Verificação heurística: confirme cada item no documento. §n = n-ésimo parágrafo do corpo.", ""]
    for a in sorted(achados, key=lambda a: (ordem[a.nivel], a.categoria)):
        norma = f" _({a.norma})_" if a.norma else ""
        linhas.append(f"- **{a.nivel}** [{a.categoria}] {a.mensagem}{norma}")
        for e in a.exemplos:
            linhas.append(f"    - {e}")
    linhas += ["", "Não verificado automaticamente: capa, folha de rosto, folha de aprovação, ordem dos "
               "elementos, conteúdo do sumário, notas de rodapé, bordas das tabelas, equações, destaque "
               "tipográfico das referências. Use references/checklist.md."]
    return "\n".join(linhas)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docx")
    ap.add_argument("--json", help="grava os achados também em JSON")
    ap.add_argument("--max-exemplos", type=int, default=5)
    a = ap.parse_args(argv)
    v = Verificador(a.docx, a.max_exemplos)
    achados = v.executar()
    print(relatorio_md(achados, a.docx))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump([asdict(x) for x in achados], f, ensure_ascii=False, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
