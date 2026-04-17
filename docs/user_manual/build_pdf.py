"""Build the First Edition PDF from the Markdown source.

Small-surface Markdown-to-PDF renderer using reportlab. Handles the
subset of Markdown used in the manual:

- Headings (# .. ####)
- Paragraphs (blank-line separated)
- Unordered lists (- )
- Ordered lists (1. )
- Tables (GitHub-flavoured | col | col |)
- Code blocks (``` fences)
- Inline: **bold**, *italic*, `code`, __bold__, _italic_

Invocation:
    python docs/user_manual/build_pdf.py
"""

from __future__ import annotations

import re
from html import escape
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    ListFlowable,
    ListItem,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


# ── Unicode font registration ─────────────────────────────────────

# reportlab's built-in Type-1 fonts (Helvetica / Courier) cover only
# WinAnsi / Latin-1. That renders every glyph outside that band —
# α, ², ³, ⁺, ⁻, ⌈, ⌉, ±, µ, °, ÷, × — as a black tofu square.
# Fix: register DejaVu Sans + DejaVu Sans Mono (both shipped with
# matplotlib and covering ~6 000 Unicode glyphs including full Greek,
# mathematical operators, super/subscripts, arrows) as the document
# font family at module load time.

def _register_unicode_fonts() -> tuple[str, str]:
    """Register DejaVu TTFs with reportlab; return (body_family, mono_family).

    Falls back to Helvetica / Courier if the DejaVu TTFs cannot be
    located (the PDF will then show tofu for non-Latin glyphs, but
    the document will at least build).
    """
    try:
        from matplotlib import get_data_path as _mpl_data
    except ImportError:
        return ("Helvetica", "Courier")

    fonts_dir = Path(_mpl_data()) / "fonts" / "ttf"
    variants = {
        ("DejaVuSans", "DejaVuSans.ttf"),
        ("DejaVuSans-Bold", "DejaVuSans-Bold.ttf"),
        ("DejaVuSans-Oblique", "DejaVuSans-Oblique.ttf"),
        ("DejaVuSans-BoldOblique", "DejaVuSans-BoldOblique.ttf"),
        ("DejaVuSansMono", "DejaVuSansMono.ttf"),
        ("DejaVuSansMono-Bold", "DejaVuSansMono-Bold.ttf"),
        ("DejaVuSansMono-Oblique", "DejaVuSansMono-Oblique.ttf"),
        ("DejaVuSansMono-BoldOblique", "DejaVuSansMono-BoldOblique.ttf"),
    }
    registered = []
    for name, fname in variants:
        ttf = fonts_dir / fname
        if ttf.exists():
            try:
                pdfmetrics.registerFont(TTFont(name, str(ttf)))
                registered.append(name)
            except Exception:
                pass

    if "DejaVuSans" not in registered:
        return ("Helvetica", "Courier")

    # Wire the four variants into a family so <b>, <i>, <b><i> work
    # inside reportlab Paragraphs.
    if all(v in registered for v in (
        "DejaVuSans", "DejaVuSans-Bold",
        "DejaVuSans-Oblique", "DejaVuSans-BoldOblique",
    )):
        pdfmetrics.registerFontFamily(
            "DejaVuSans",
            normal="DejaVuSans",
            bold="DejaVuSans-Bold",
            italic="DejaVuSans-Oblique",
            boldItalic="DejaVuSans-BoldOblique",
        )
    if all(v in registered for v in (
        "DejaVuSansMono", "DejaVuSansMono-Bold",
        "DejaVuSansMono-Oblique", "DejaVuSansMono-BoldOblique",
    )):
        pdfmetrics.registerFontFamily(
            "DejaVuSansMono",
            normal="DejaVuSansMono",
            bold="DejaVuSansMono-Bold",
            italic="DejaVuSansMono-Oblique",
            boldItalic="DejaVuSansMono-BoldOblique",
        )

    body = "DejaVuSans"
    mono = "DejaVuSansMono" if "DejaVuSansMono" in registered else "Courier"
    return (body, mono)


_BODY_FONT, _MONO_FONT = _register_unicode_fonts()
_BODY_BOLD = _BODY_FONT + "-Bold" if _BODY_FONT != "Helvetica" else "Helvetica-Bold"
_BODY_ITALIC = _BODY_FONT + "-Oblique" if _BODY_FONT != "Helvetica" else "Helvetica-Oblique"


# ── Styles ────────────────────────────────────────────────────────

def _build_styles() -> dict:
    base = getSampleStyleSheet()
    styles = {}
    styles["body"] = ParagraphStyle(
        "body", parent=base["BodyText"],
        fontName=_BODY_FONT, fontSize=10.5, leading=14,
        spaceAfter=6,
    )
    styles["h1"] = ParagraphStyle(
        "h1", parent=base["Heading1"],
        fontName=_BODY_BOLD, fontSize=22, leading=26,
        spaceBefore=18, spaceAfter=12, textColor=colors.HexColor("#1f2937"),
    )
    styles["h2"] = ParagraphStyle(
        "h2", parent=base["Heading2"],
        fontName=_BODY_BOLD, fontSize=16, leading=20,
        spaceBefore=14, spaceAfter=8, textColor=colors.HexColor("#1f2937"),
    )
    styles["h3"] = ParagraphStyle(
        "h3", parent=base["Heading3"],
        fontName=_BODY_BOLD, fontSize=13, leading=17,
        spaceBefore=10, spaceAfter=6, textColor=colors.HexColor("#374151"),
    )
    styles["h4"] = ParagraphStyle(
        "h4", parent=base["Heading4"],
        fontName=_BODY_BOLD, fontSize=11, leading=14,
        spaceBefore=8, spaceAfter=4, textColor=colors.HexColor("#4b5563"),
    )
    styles["code"] = ParagraphStyle(
        "code", parent=base["Code"],
        fontName=_MONO_FONT, fontSize=8.5, leading=11,
        backColor=colors.HexColor("#f3f4f6"),
        borderColor=colors.HexColor("#d1d5db"),
        borderWidth=0.5, borderPadding=6,
        spaceBefore=4, spaceAfter=6, leftIndent=4, rightIndent=4,
    )
    styles["bullet"] = ParagraphStyle(
        "bullet", parent=styles["body"], leftIndent=6, bulletIndent=0,
    )
    styles["caption"] = ParagraphStyle(
        "caption", parent=base["BodyText"],
        fontName=_BODY_ITALIC, fontSize=9, leading=11,
        textColor=colors.HexColor("#6b7280"), spaceAfter=6,
    )
    return styles


# ── Inline markdown → reportlab mini-HTML ─────────────────────────

_INLINE_CODE_RE = re.compile(r"`([^`]+)`")
_BOLD_STAR_RE = re.compile(r"\*\*([^*]+)\*\*")
_BOLD_UND_RE = re.compile(r"__([^_]+)__")
_ITAL_STAR_RE = re.compile(r"(?<!\*)\*([^*\n]+)\*(?!\*)")
_ITAL_UND_RE = re.compile(r"(?<!_)_([^_\n]+)_(?!_)")


def _inline_md_to_html(text: str) -> str:
    """Convert inline markdown to reportlab's limited HTML subset.

    Code spans are extracted first into placeholders so that their
    contents (which frequently contain underscores / asterisks) are
    not mis-parsed by the bold / italic regexes. After the other
    transforms the placeholders are substituted back.
    """
    placeholders: list[str] = []
    def _stash_code(m):
        idx = len(placeholders)
        placeholders.append(
            f'<font face="{_MONO_FONT}" color="#b91c1c">'
            f'{escape(m.group(1), quote=False)}</font>'
        )
        return f"\x00CODE{idx}\x00"

    # 1. Stash code spans (they must survive bold/italic rewrites)
    tmp = _INLINE_CODE_RE.sub(_stash_code, text)
    # 2. Escape the remaining text
    tmp = escape(tmp, quote=False)
    # 3. Bold, then italic (order matters for overlapping stars)
    tmp = _BOLD_STAR_RE.sub(r"<b>\1</b>", tmp)
    tmp = _BOLD_UND_RE.sub(r"<b>\1</b>", tmp)
    tmp = _ITAL_STAR_RE.sub(r"<i>\1</i>", tmp)
    tmp = _ITAL_UND_RE.sub(r"<i>\1</i>", tmp)
    # 4. Put code spans back
    for idx, html in enumerate(placeholders):
        tmp = tmp.replace(f"\x00CODE{idx}\x00", html)
    return tmp


# ── Block parser ──────────────────────────────────────────────────

def _parse_table_block(lines: list[str]) -> Table:
    """Parse a GitHub-flavoured Markdown table into a reportlab Table."""
    # Filter the separator row (----|----)
    rows = []
    for ln in lines:
        stripped = ln.strip()
        if not stripped.startswith("|"):
            break
        cells = [c.strip() for c in stripped.strip("|").split("|")]
        # skip alignment row
        if all(re.fullmatch(r":?-{3,}:?", c.strip()) for c in cells if c):
            continue
        rows.append(cells)

    if not rows:
        return None

    body_style = ParagraphStyle(
        "tblBody", fontName=_BODY_FONT, fontSize=9, leading=11,
    )
    header_style = ParagraphStyle(
        "tblHead", fontName=_BODY_BOLD, fontSize=9, leading=11,
        textColor=colors.white,
    )

    cell_rows = []
    for i, row in enumerate(rows):
        style = header_style if i == 0 else body_style
        cell_rows.append([Paragraph(_inline_md_to_html(c), style) for c in row])

    n_cols = max(len(r) for r in cell_rows)
    for r in cell_rows:
        while len(r) < n_cols:
            r.append(Paragraph("", body_style))

    tbl = Table(cell_rows, hAlign="LEFT", repeatRows=1)
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f2937")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#9ca3af")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.white, colors.HexColor("#f9fafb")]),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    return tbl


def _md_to_flowables(md: str, styles: dict) -> list:
    """Walk Markdown line-by-line and emit reportlab flowables."""
    flowables: list = []
    lines = md.split("\n")
    i = 0
    n = len(lines)

    def _para_buffer_flush(buf: list[str]):
        if buf:
            text = " ".join(buf).strip()
            if text:
                flowables.append(Paragraph(_inline_md_to_html(text), styles["body"]))
            buf.clear()

    para_buf: list[str] = []

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # Fenced code block
        if stripped.startswith("```"):
            _para_buffer_flush(para_buf)
            i += 1
            code_lines = []
            while i < n and not lines[i].startswith("```"):
                code_lines.append(lines[i])
                i += 1
            i += 1  # consume closing fence
            code_text = "\n".join(code_lines)
            flowables.append(Preformatted(code_text, styles["code"]))
            continue

        # Horizontal rule
        if re.fullmatch(r"\s*-{3,}\s*", line):
            _para_buffer_flush(para_buf)
            flowables.append(Spacer(1, 4))
            flowables.append(
                Table([[""]], colWidths=[170 * mm], rowHeights=[0.5],
                      style=TableStyle([
                          ("LINEBELOW", (0, 0), (-1, -1), 0.5,
                           colors.HexColor("#9ca3af")),
                      ])),
            )
            flowables.append(Spacer(1, 4))
            i += 1
            continue

        # Headings
        m = re.match(r"^(#{1,4})\s+(.*)$", line)
        if m:
            _para_buffer_flush(para_buf)
            level = len(m.group(1))
            text = m.group(2).strip()
            flowables.append(Paragraph(_inline_md_to_html(text),
                                       styles[f"h{level}"]))
            i += 1
            continue

        # Table
        if stripped.startswith("|"):
            _para_buffer_flush(para_buf)
            tbl_start = i
            while i < n and lines[i].strip().startswith("|"):
                i += 1
            tbl = _parse_table_block(lines[tbl_start:i])
            if tbl is not None:
                flowables.append(tbl)
                flowables.append(Spacer(1, 6))
            continue

        # Unordered list
        if re.match(r"^\s*-\s+", line):
            _para_buffer_flush(para_buf)
            items = []
            while i < n:
                m2 = re.match(r"^(\s*)-\s+(.*)$", lines[i])
                if not m2:
                    break
                items.append(ListItem(
                    Paragraph(_inline_md_to_html(m2.group(2)), styles["bullet"]),
                    leftIndent=14,
                ))
                i += 1
            flowables.append(ListFlowable(
                items, bulletType="bullet", leftIndent=12, bulletFontSize=9,
            ))
            continue

        # Ordered list
        if re.match(r"^\s*\d+\.\s+", line):
            _para_buffer_flush(para_buf)
            items = []
            while i < n:
                m2 = re.match(r"^(\s*)\d+\.\s+(.*)$", lines[i])
                if not m2:
                    break
                items.append(ListItem(
                    Paragraph(_inline_md_to_html(m2.group(2)), styles["bullet"]),
                    leftIndent=14,
                ))
                i += 1
            flowables.append(ListFlowable(
                items, bulletType="1", leftIndent=12, bulletFontSize=9,
            ))
            continue

        # Blank line ends the current paragraph
        if not stripped:
            _para_buffer_flush(para_buf)
            flowables.append(Spacer(1, 4))
            i += 1
            continue

        # Default: accumulate into current paragraph
        para_buf.append(stripped)
        i += 1

    _para_buffer_flush(para_buf)
    return flowables


# ── Build ────────────────────────────────────────────────────────

def build(md_path: Path, pdf_path: Path) -> None:
    md = md_path.read_text(encoding="utf-8")
    styles = _build_styles()
    flowables = _md_to_flowables(md, styles)

    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        leftMargin=20 * mm, rightMargin=20 * mm,
        topMargin=20 * mm, bottomMargin=20 * mm,
        title="Polysaccharide-Based Microsphere Emulsification Simulator — First Edition",
        author="EmulSim project",
    )

    def _footer(canvas, doc):
        canvas.saveState()
        canvas.setFont(_BODY_FONT, 8)
        canvas.setFillColor(colors.HexColor("#6b7280"))
        canvas.drawString(
            20 * mm, 12 * mm,
            "EmulSim — Polysaccharide Microsphere Simulator — First Edition",
        )
        canvas.drawRightString(
            doc.pagesize[0] - 20 * mm, 12 * mm,
            f"Page {doc.page}",
        )
        canvas.restoreState()

    doc.build(flowables, onFirstPage=_footer, onLaterPages=_footer)


if __name__ == "__main__":
    here = Path(__file__).parent
    md = here / "polysaccharide_microsphere_simulator_first_edition.md"
    pdf = here / "polysaccharide_microsphere_simulator_first_edition.pdf"
    build(md, pdf)
    size_kb = pdf.stat().st_size / 1024.0
    print(f"Built: {pdf}  ({size_kb:.1f} KB)")
