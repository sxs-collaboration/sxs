"""MkDocs hook that renders the metadata-field documentation as tables.

Any page containing the placeholder

    <!-- metadata-fields-tables -->

has it replaced by one section per entry in `metadata_field_groups`, each with
a table describing the fields in that group, as defined in
`sxs/metadata/metadata_docs.py`.

On such pages, Pandoc-style citations — `[@Key]` or `[@Key1; @Key2]` — are
replaced by numbered links to a "Bibliography" section appended to the end of
the page.  Each key must be present in `docs/references.bib`; otherwise the
build fails.  Each reference links to its DOI and/or arXiv preprint.

`metadata_docs.py` only depends on the standard library, so we load it
directly from its file rather than importing the full `sxs` package (which
would pull in heavy dependencies at docs-build time).

"""

import html
import importlib.util
import re
from math import inf
from pathlib import Path

import bibtexparser
import markdown
from bibtexparser import middlewares
from mkdocs.exceptions import PluginError

PLACEHOLDER = "<!-- metadata-fields-tables -->"

_root = Path(__file__).resolve().parent
_metadata_docs_path = _root / "sxs" / "metadata" / "metadata_docs.py"
_bib_path = _root / "docs" / "references.bib"

# Matches `[@Key]` and `[@Key1; @Key2]`; INSPIRE keys look like `Author:2019abc`.
_citation_re = re.compile(r"\[(@[^\s;\[\]]+(?:\s*;\s*@[^\s;\[\]]+)*)\]")


def _load_metadata_docs():
    spec = importlib.util.spec_from_file_location(
        "_sxs_metadata_docs", _metadata_docs_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _md_to_html(text):
    """Render a field description (markdown) as HTML for use inside a table cell."""
    return markdown.markdown(text)


def _format_type(t):
    name = t.__name__ if isinstance(t, type) and not hasattr(t, "__args__") else repr(t)
    return f"<code>{html.escape(name)}</code>"


def _format_versions(field):
    """Range of metadata versions containing the field, e.g., "0–1" or "2+".

    A field deprecated in version `d` was last present in version `d - 1`.
    """
    introduced = field["introduced"]
    deprecated = field["deprecated"]
    if deprecated == inf:
        return f"{introduced}+"
    last = deprecated - 1
    if last < introduced:
        raise PluginError(
            f"Metadata field {field['name']!r} is deprecated in version {deprecated}, "
            f"but was introduced in version {introduced}"
        )
    if last == introduced:
        return f"{introduced}"
    return f"{introduced}–{last}"


def _format_presence(field, FP):
    present_for = field["present_for"]
    if present_for == FP.ALL:
        return "All"
    names = [member.name for member in FP if member in present_for and member != FP.ALL]
    return ", ".join(html.escape(n) for n in names)


def _breakable(name):
    """Allow long field names to wrap after underscores (and nowhere else)."""
    return html.escape(name).replace("_", "_<wbr>")


def _field_table(fields, FP):
    """Two-column table: field name and its properties | description."""
    rows = []
    for field in fields:
        name = html.escape(field["name"])
        anchor = f"metadata-{name}"
        meta = [
            f"<span>Type: {_format_type(field['intended_type'])}</span>",
            f"<span>Versions: {_format_versions(field)}</span>",
            f"<span>Present for: {_format_presence(field, FP)}</span>",
        ]
        if field["computed"]:
            meta.append("<span><em>Computed by <code>sxs</code></em></span>")
        field_cell = (
            f'<a class="field-name" href="#{anchor}"><code>{_breakable(field["name"])}</code></a>'
            f'<div class="field-meta">{"".join(meta)}</div>'
        )
        rows.append(
            f'<tr id="{anchor}">'
            f'<td class="field">{field_cell}</td>'
            f'<td class="description">{_md_to_html(field["description"])}</td>'
            "</tr>"
        )
    return (
        '<table class="metadata-fields">\n'
        "<thead><tr><th>Field</th><th>Description</th></tr></thead>\n"
        "<tbody>\n" + "\n".join(rows) + "\n</tbody>\n</table>\n"
    )


def _render_tables(md):
    FP = md.FieldPresence
    unknown = {f["group"] for f in md.metadata_fields} - set(md.metadata_field_groups)
    if unknown:
        raise PluginError(
            f"Metadata fields use groups not listed in metadata_field_groups: {sorted(unknown)}"
        )
    sections = []
    for group_key, group_title in md.metadata_field_groups.items():
        fields = [f for f in md.metadata_fields if f["group"] == group_key]
        if not fields:
            continue
        current = [f for f in fields if f["deprecated"] == inf]
        deprecated = [f for f in fields if f["deprecated"] != inf]
        section = f"## {group_title}\n\n"
        if current:
            section += _field_table(current, FP)
        if deprecated:
            # Collapsed by default; the content is still in the page, so the
            # site search indexes it.
            section += (
                '\n<details class="deprecated-fields">\n'
                f"<summary>Deprecated fields ({len(deprecated)})</summary>\n"
                f"{_field_table(deprecated, FP)}"
                "</details>\n"
            )
        sections.append(section)
    return "\n\n".join(sections)


### Citations ###

def _load_bibliography():
    library = bibtexparser.parse_file(
        str(_bib_path),
        append_middleware=[
            middlewares.LatexDecodingMiddleware(),
            middlewares.SeparateCoAuthors(),
            middlewares.SplitNameParts(),
        ],
    )
    if library.failed_blocks:
        problems = "\n".join(
            f"  - {type(block).__name__}: {block.raw[:80]!r}"
            for block in library.failed_blocks
        )
        raise PluginError(f"Failed to parse some entries in {_bib_path}:\n{problems}")
    return library.entries_dict


def _cited_keys(text):
    return [
        key.strip().lstrip("@")
        for match in _citation_re.finditer(text)
        for key in match.group(1).split(";")
    ]


def _format_name(name):
    parts = name.first + name.von + name.last
    full = " ".join(parts)
    if name.jr:
        full += ", " + " ".join(name.jr)
    return full


def _format_authors(authors):
    others = bool(authors) and _format_name(authors[-1]) == "others"
    if others:
        authors = authors[:-1]
    names = [_format_name(a) for a in authors]
    if others or len(names) > 3:
        return f"{names[0]} <em>et al.</em>"
    if len(names) <= 2:
        return " and ".join(names)
    return ", ".join(names[:-1]) + ", and " + names[-1]


def _format_reference(entry):
    f = entry.fields_dict
    get = lambda k: f[k].value if k in f else ""
    esc = html.escape

    parts = []
    if "author" in f:
        parts.append(_format_authors(f["author"].value) + ",")
    if "title" in f:
        parts.append(f"“{esc(get('title'))},”")
    if "journal" in f:
        journal = f"<em>{esc(get('journal'))}</em>"
        if get("volume"):
            journal += f" <strong>{esc(get('volume'))}</strong>"
        if get("pages"):
            journal += ", " + esc(get("pages").replace("--", "–"))
        parts.append(journal)
    if get("year"):
        parts.append(f"({esc(get('year'))})")
    text = " ".join(parts)

    links = []
    if get("doi"):
        doi = get("doi")
        links.append(
            f'<a href="https://doi.org/{esc(doi, quote=True)}">doi:{esc(doi)}</a>'
        )
    if get("eprint"):
        eprint = get("eprint")
        links.append(
            f'<a href="https://arxiv.org/abs/{esc(eprint, quote=True)}">arXiv:{esc(eprint)}</a>'
        )
    if links:
        text += ". " + " · ".join(links)
    return text


def _ref_anchor(key):
    return "ref-" + html.escape(key, quote=True)


def _apply_citations(text, md):
    """Number citations in `text`, link them, and append a Bibliography section."""
    keys = list(dict.fromkeys(_cited_keys(text)))  # unique, in order of appearance
    if not keys:
        return text

    bibliography = _load_bibliography()
    missing = [k for k in keys if k not in bibliography]
    if missing:
        where = {
            k: [f["name"] for f in md.metadata_fields if k in _cited_keys(f["description"])]
            for k in missing
        }
        details = "\n".join(
            f"  - {k}" + (f" (cited in field(s): {', '.join(v)})" if v else "")
            for k, v in where.items()
        )
        raise PluginError(
            f"Citation key(s) not found in {_bib_path.relative_to(_root)}:\n{details}"
        )

    number = {key: i + 1 for i, key in enumerate(keys)}

    def replace(match):
        cited = [k.strip().lstrip("@") for k in match.group(1).split(";")]
        links = ", ".join(
            f'<a href="#{_ref_anchor(k)}" title="{html.escape(k, quote=True)}">{number[k]}</a>'
            for k in cited
        )
        return f'<span class="citation">[{links}]</span>'

    text = _citation_re.sub(replace, text)

    items = "\n".join(
        f'<li id="{_ref_anchor(k)}">{_format_reference(bibliography[k])}</li>'
        for k in keys
    )
    return f'{text}\n\n## Bibliography\n\n<ol class="references">\n{items}\n</ol>\n'


### MkDocs event handlers ###

def on_page_markdown(markdown_text, page, config, files):
    if PLACEHOLDER not in markdown_text:
        return markdown_text
    md = _load_metadata_docs()
    text = markdown_text.replace(PLACEHOLDER, _render_tables(md))
    return _apply_citations(text, md)


def on_serve(server, config, builder):
    # Rebuild on `mkdocs serve` when the field definitions change
    # (references.bib is already watched, since it lives in docs/).
    server.watch(str(_metadata_docs_path))
    return server
