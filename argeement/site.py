"""Summary web page of the latest release (GitHub Pages).

    python -m argeement site --release releases/current/latest --out docs \
        [--figures manuscript/figures] [--config site.yaml]

A single index.html without external dependencies is produced (tables, figures, download links);
the output directory can be published directly as the GitHub Pages source. All numbers come from
the release files.
"""
import html
import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from .build import RELATION_RANK, STRONG_MAX_RANK
from .terms import licence_html

DBS = ["amrfinder", "card", "resfinder", "disinfinder"]
NAME = {"amrfinder": "AMRFinderPlus", "card": "CARD", "resfinder": "ResFinder",
        "disinfinder": "DisinFinder", "pointfinder": "PointFinder", "rgi": "RGI"}
REL_EN = {"identical_protein": "identical protein", "near_identical": "near-identical",
          "length_variant": "length variant", "same_family": "same family", "homolog": "homolog"}
# Figures of the page, in the order of the paper, with the paper's captions (HTML). {name} is
# filled from the macros of the release (macros.json, written by manuscript/make_assets.py --lang en
# and shipped with the package); "asset": static figure in argeement/assets/, "release": PNG made
# for every release by manuscript/figures.R.
FIGURES = [
    ("fig_workflow", "asset",
     "The processing workflow of ARGeement. ResFinder, DisinFinder and PointFinder are databases of "
     "the Center for Genomic Epidemiology (CGE). PointFinder and the mutation part of the databases "
     "enter only the mutation inventory and take no part in the sequence-based mapping."),
    ("fig_hierarchy", "asset",
     "Part of the NCBI gene hierarchy. In allele-rich families the node of an entry is an allele "
     "(e.g. <i>blaTEM-1</i>), and the family is then the parent node (<i>blaTEM</i>); for single "
     "genes (e.g. <i>sul1</i>) the node of the entry is itself the family. Grey nodes indicate "
     "further branches of the hierarchy not detailed in the figure."),
    ("fig1_tartalom", "release",
     "Size of the databases. Number of entries in the source, of entries included in the mapping, of "
     "protein-level entries and of unique proteins. ResFinder and DisinFinder (Center for Genomic "
     "Epidemiology, CGE) are shown together; as some proteins are encoded by both databases, the "
     "number of unique CGE proteins is slightly smaller than the sum for the two databases. The "
     "source count of AMRFinderPlus includes all rows of the Reference Gene Catalog, together with "
     "the rows of mutation, virulence and stress response genes."),
    ("fig2_osztalyok", "release",
     "Number of protein-level entries by harmonised drug class. An entry appears in every class "
     "assigned to it. The figure shows all {nClassAll} harmonised classes, in alphabetical order."),
    ("fig3_atfedes", "release",
     "Overlap between the databases. Each cell shows the percentage of the entries of the database "
     "in the row that have a close counterpart (identical protein, near-identical or length variant) "
     "in the database in the column. The value in parentheses counts only counterparts encoding the "
     "identical protein. For example, the AMRFinderPlus column of the ResFinder row shows the share "
     "of ResFinder entries that have a close or identical-protein counterpart in AMRFinderPlus. The "
     "figure is not symmetrical, because the shares always refer to the entries of the database in "
     "the row, and several entries of one database may correspond to the same entry in the other."),
    ("fig4_nevek", "release",
     "Name concordance. a) Names of identical proteins. Names of cross-database entry pairs "
     "encoding the same protein, by category. b) Proteins of identical names. Relation of the "
     "proteins of cross-database entry pairs with the same name."),
    ("fig6_csaladok", "release",
     "Amino acid identity of all protein pairs of the {famN} NCBI families of the chromosomal "
     "\u03b2-lactamases of <i>Klebsiella</i> species ({famSizesList} proteins), within families "
     "({famPairsWithin} pairs) and between families ({famPairsBetween} pairs), as relative "
     "frequencies. The dashed line is the {thrFam}% identity threshold of the same family; "
     "classification also requires a coverage of at least {thrFamCovPct} on both proteins."),
    ("fig5_mutaciok", "release",
     "Mutation-based determinants by genus in the three sources. The AMRFinderPlus column counts "
     "the listed mutations, the CARD column the point mutations listed in the variant models, and "
     "the PointFinder column the mutation positions. As PointFinder lists several resistant amino "
     "acids at one position as a single item, the numbers are only roughly comparable between the "
     "columns, but all three show the genus coverage exactly."),
    ("fig8_egyezes", "release",
     "Agreement of the gene calls of the three annotation tools on {concGenomes} genomes, by tool "
     "pair, at the loci called by both tools. The bars show the relation of the references of the "
     "two tools in the release, and the diamond the name-based agreement before mapping. The number "
     "in parentheses after each tool pair is the number of loci called by both tools. No locus "
     "occurred at which a reference could not be linked to the release, so this category is not "
     "shown."),
]
FIG_EN = [stem for stem, kind, _ in FIGURES if kind == "release"]   # PNGs copied into the package
ASSETS = Path(__file__).parent / "assets"


class _Macros(dict):
    def __missing__(self, key):
        return "?"


# Figures are shown at their print size times this factor, so that the lettering is the same size
# in every figure (the R figures are PNGs at 200 dpi, the static ones SVGs sized in points).
FIG_SCALE = 1.3


def _display_width(path):
    """Display width (CSS px) of a figure from its physical size; None if it cannot be read."""
    import re
    import struct
    try:
        if path.suffix == ".png":
            b = path.read_bytes()
            width = struct.unpack(">I", b[16:20])[0]
            dpi = 200.0
            i = b.find(b"pHYs")
            if i > 0:
                ppx, _, unit = struct.unpack(">IIB", b[i + 4:i + 13])
                if unit == 1 and ppx:
                    dpi = ppx * 0.0254
            return round(width / dpi * 96 * FIG_SCALE)
        if path.suffix == ".svg":
            m = re.search(r'<svg[^>]*\swidth="([0-9.]+)(pt)?"', path.read_text(errors="ignore")[:2000])
            if m:   # pdftocairo writes the size in points, without a unit
                return round(float(m.group(1)) * 96 / 72 * FIG_SCALE)
    except (OSError, struct.error, ValueError):
        pass
    return None


def _caption(template, macros):
    """The caption with the values of the release; a missing value prints as ?."""
    return template.format_map(_Macros(macros or {}))


def _load_macros(*dirs):
    for d in dirs:
        if d and (Path(d) / "macros.json").exists():
            return json.loads((Path(d) / "macros.json").read_text(encoding="utf-8"))
    return {}



def _fmt(n):
    return f"{int(n):,}".replace(",", "\u2009")


def _pct(a, b, d=1):
    return f"{100 * a / b:.{d}f}%" if b else "–"


def _table(header, rows, cls=""):
    h = "".join(f"<th>{html.escape(str(x))}</th>" for x in header)
    b = "".join("<tr>" + "".join(f"<td>{html.escape(str(x))}</td>" for x in r) + "</tr>"
                for r in rows)
    return f'<table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table>'


def _read(rel, name):
    p = Path(rel) / name
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False) if p.exists() else None


def collect(rel):
    """Summary numbers of a release (as a dict: common source of the page and the JSON output)."""
    rel = Path(rel)
    man = json.loads((rel / "manifest.json").read_text())
    e = _read(rel, "entries.tsv")
    m = _read(rel, "mappings.tsv")
    usable = e[e["usable"] == "True"]
    S = dict(release=rel.resolve().name, created=man.get("created", ""),
             profile=man.get("profile", ""), argeement_version=man.get("argeement_version", ""),
             mmseqs_version=man.get("mmseqs_version", ""),
             database_versions=man.get("database_versions", {}), params=man.get("params", {}))
    src = e[e["parent"] == ""] if "parent" in e.columns else e
    S["databases"] = [dict(db=d, source_entries=int((src["db"] == d).sum()),
                           protein_level_entries=int((usable["db"] == d).sum()),
                           unique_proteins=int(usable.loc[usable["db"] == d, "prot_id"].nunique()))
                      for d in DBS if (e["db"] == d).any()]
    S["totals"] = dict(source_entries=int(len(src)), protein_level_entries=int(len(usable)),
                       unique_proteins=int(usable["prot_id"].nunique()),
                       mappings=int(len(m)) if m is not None else 0)
    # directed close coverage
    if m is not None:
        strong = m[m["relation"].map(RELATION_RANK) <= STRONG_MAX_RANK]
        has = set(zip(strong["db_a"], strong["entry_a"], strong["db_b"]))
        ov = []
        for a in DBS:
            ea = usable[usable["db"] == a]
            if ea.empty:
                continue
            row = dict(db=a, n=int(len(ea)))
            for b in DBS:
                if b != a and (usable["db"] == b).any():
                    row[b] = sum((a, x, b) in has for x in ea["entry_id"]) / len(ea)
            ov.append(row)
        S["overlap"] = ov
        S["relations"] = {k: int(v) for k, v in m["relation"].value_counts().items()}
    nc = man.get("name_concordance") or {}
    S["name_concordance"] = nc
    S["ncbi_validation"] = man.get("ncbi_validation") or {}
    v = rel / "validation_synthetic" / "summary.json"
    if v.exists():
        S["validation"] = json.loads(v.read_text())
    c = rel / "concordance" / "summary.json"
    if c.exists():
        S["concordance"] = json.loads(c.read_text())
    d = rel / "diff" / "diff.tsv"
    if d.exists():
        dd = pd.read_csv(d, sep="\t", dtype=str, keep_default_na=False)
        S["diff"] = {k: int(v) for k, v in dd["change"].value_counts().items()}
    return S


CSS = """
:root{--bg:#fff;--fg:#1d2329;--muted:#5b6670;--line:#dfe3e7;--accent:#0b6e99;--card:#f6f8fa}
@media (prefers-color-scheme:dark){:root{--bg:#121619;--fg:#e6e9ec;--muted:#9aa5ae;
 --line:#2b3238;--accent:#5ab0d6;--card:#1a2025}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
 font:15px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:980px;margin:0 auto;padding:2rem 1.2rem 4rem}
h1{font-size:1.8rem;margin:.2rem 0}h2{font-size:1.25rem;margin-top:2.4rem;
 border-bottom:1px solid var(--line);padding-bottom:.3rem}
.lead{color:var(--muted)}a{color:var(--accent)}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:.8rem;margin:1.2rem 0}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:.8rem 1rem}
.card b{display:block;font-size:1.4rem}.card span{color:var(--muted);font-size:.85rem}
.tw{overflow-x:auto}table{border-collapse:collapse;width:100%;margin:.8rem 0;font-size:.92rem}
th,td{border-bottom:1px solid var(--line);padding:.35rem .6rem;text-align:right;white-space:nowrap}
th:first-child,td:first-child{text-align:left}thead th{color:var(--muted);font-weight:600}
figure{margin:1.6rem 0}figure img{display:block;margin:0 auto .5rem;max-width:100%;height:auto;background:#fff;border-radius:6px;
 border:1px solid var(--line)}figcaption{color:var(--muted);font-size:.88rem}
code{background:var(--card);padding:.1rem .3rem;border-radius:4px}
footer{margin-top:3rem;color:var(--muted);font-size:.85rem}
.cite{padding-left:1rem;border-left:3px solid var(--muted)}
"""


def render(S, cfg, figs):
    t = cfg.get("title", "ARGeement")
    links = []
    for k, lab in (("github", "Source code (GitHub)"), ("zenodo", "Archived releases (Zenodo)"),
                   ("paper", "Paper")):
        if cfg.get(k):
            links.append(f'<a href="{html.escape(cfg[k])}">{lab}</a>')
    tot = S["totals"]
    P = [f"<h1>{html.escape(t)}</h1>",
         f'<p class="lead">{html.escape(cfg.get("subtitle", "Sequence-based, versioned mapping of antimicrobial resistance gene databases"))}</p>',
         f"<p>Latest release <b>{html.escape(S['release'])}</b> (profile <code>{html.escape(S['profile'])}</code>, "
         f"created {html.escape(S['created'][:19].replace('T', ' '))} UTC). " + " · ".join(links) + "</p>",
         '<div class="cards">'
         f"<div class='card'><b>{_fmt(tot['source_entries'])}</b><span>source entries</span></div>"
         f"<div class='card'><b>{_fmt(tot['protein_level_entries'])}</b><span>protein-level entries</span></div>"
         f"<div class='card'><b>{_fmt(tot['unique_proteins'])}</b><span>unique proteins</span></div>"
         f"<div class='card'><b>{_fmt(tot['mappings'])}</b><span>entry–entry relations</span></div>"
         "</div>"]
    if cfg.get("zenodo_doi"):
        P.append(f"<p>Download this release: <a href=\"https://doi.org/{html.escape(cfg['zenodo_doi'])}\">"
                 f"doi:{html.escape(cfg['zenodo_doi'])}</a></p>")
    # versions
    P.append("<h2>Database versions</h2>")
    rows = [[NAME.get(k, k), v] for k, v in sorted(S["database_versions"].items())]
    rows.append(["MMseqs2", S["mmseqs_version"]])
    P.append('<div class="tw">' + _table(["Database / tool", "Version"], rows) + "</div>")
    # tartalom
    P.append("<h2>Content</h2>")
    P.append('<div class="tw">' + _table(
        ["Database", "Source entries", "Protein-level entries", "Unique proteins"],
        [[NAME[d["db"]], _fmt(d["source_entries"]), _fmt(d["protein_level_entries"]),
          _fmt(d["unique_proteins"])] for d in S["databases"]]) + "</div>")
    P.append("<p class='lead'>Protein-level entries are entries translated as a whole plus the "
             "retained open reading frames of entries split into ORFs (operons, reverse-strand "
             "genes).</p>")
    # overlap
    if S.get("overlap"):
        P.append("<h2>Overlap between databases</h2>")
        cols = [d["db"] for d in S["databases"]]
        rows = [[f"{NAME[o['db']]} (n = {_fmt(o['n'])})"] +
                [("–" if c == o["db"] else f"{100 * o.get(c, 0):.1f}%") for c in cols] for o in S["overlap"]]
        P.append('<div class="tw">' + _table(["From \\ to"] + [NAME[c] for c in cols], rows) + "</div>")
        P.append("<p class='lead'>Share of protein-level entries of the row database that have a close "
                 "counterpart (identical protein, near-identical or length variant) in the column "
                 "database.</p>")
    if S.get("relations"):
        P.append('<div class="tw">' + _table(
            ["Relation type", "Directed relations"],
            [[REL_EN.get(k, k), _fmt(v)] for k, v in sorted(S["relations"].items(),
                                                            key=lambda x: RELATION_RANK.get(x[0], 9))])
            + "</div>")
    # konkordancia
    c = S.get("concordance")
    if c:
        P.append("<h2>Real-genome use case</h2>")
        P.append(f"<p>{_fmt(c['genomes'])} complete RefSeq genomes annotated with AMRFinderPlus, RGI and "
                 "ResFinder/DisinFinder (BLAST alignments, database versions of this release).</p>")
        rows = []
        for r in c["pairs"]:
            a, b = r["pair"].split("-")
            mt = r["loci_matched"]
            rows.append([f"{NAME[a]} – {NAME[b]}", _fmt(mt), f"{100 * r['locus_concordance']:.1f}%",
                         f"{100 * r['name_concordance_before']:.1f}%",
                         f"{100 * r['name_concordance_after']:.1f}%"])
        P.append('<div class="tw">' + _table(
            ["Tool pair", "Shared loci", "Locus agreement", "Same name", "Same or close protein"],
            rows) + "</div>")
    # validation
    v = S.get("validation")
    if v:
        o = v["overall"]
        P.append("<h2>Synthetic validation</h2>")
        P.append(f"<p>{_fmt(v['n_cases'])} synthetic variants of {_fmt(v['n_proteins'])} curated proteins "
                 f"(seed {v.get('seed', '')}): original protein among best hits in "
                 f"{100 * o['truth_found']:.1f}% of cases, relation type correct in "
                 f"{100 * o['relation_ok']:.1f}%, number of differences correct in "
                 f"{100 * o['ndiff_ok']:.1f}%.</p>")
    # changes
    if S.get("diff"):
        P.append("<h2>Changes since the previous release</h2>")
        P.append('<div class="tw">' + _table(["Change", "Count"],
                                             [[k.replace("_", " "), _fmt(n)] for k, n in sorted(S["diff"].items())])
                 + "</div>")
    # figures, numbered in the order of the paper; captions are trusted HTML (FIGURES)
    if figs:
        P.append("<h2>Figures</h2>")
        for i, (f, cap) in enumerate(figs, 1):
            alt = html.escape(__import__("re").sub(r"<[^>]+>", "", cap.split(". ")[0]))
            w = _display_width(f)
            wa = f' width="{w}"' if w else ""
            P.append(f'<figure id="fig{i}"><img src="img/{f.name}"{wa} alt="{alt}" loading="lazy">'
                     f"<figcaption><b>Figure {i}.</b> {cap}</figcaption></figure>")
    P.append("<h2>Use</h2><p>Install and run the pipeline as described in the "
             + (f'<a href="{html.escape(cfg["github"])}#readme">README</a>' if cfg.get("github") else "README")
             + ". Translate a gene list with <code>argeement lookup</code>, search new sequences with "
               "<code>argeement query</code>.</p>")
    P.append("<h2>Licence</h2>" + licence_html())
    c = cfg.get("citation") or {}
    if c.get("title") and c.get("doi"):
        doi = html.escape(c["doi"])
        ref = (f"{html.escape(c.get('authors', ''))}. {html.escape(c['title'])}. "
               f"<i>{html.escape(c.get('journal', ''))}</i> {html.escape(c.get('details', ''))}; "
               f'doi: <a href="https://doi.org/{doi}">https://doi.org/{doi}</a>')
        P.append("<h2>Citation</h2><p>If you use ARGeement, please cite:</p>"
                 f'<p class="cite">{ref}</p>')
        if cfg.get("zenodo_doi"):
            z = html.escape(cfg["zenodo_doi"])
            P.append("<p>Please also cite the release you used; this release: "
                     f'<a href="https://doi.org/{z}">doi:{z}</a>.</p>')
    P.append(f"<footer>Generated from <code>{html.escape(S['release'])}</code> by ARGeement "
             f"{html.escape(str(S['argeement_version']))}. Machine-readable summary: "
             "<a href=\"summary.json\">summary.json</a>.</footer>")
    return ("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<title>{html.escape(t)} – latest release</title><style>{CSS}</style></head>"
            f"<body><main>{''.join(P)}</main></body></html>")


def build_site(release=None, out="docs", figures=None, config=None, package_dir=None):
    """The web page from a release directory or from a downloaded package directory
    (package_dir: site/summary.json, site/img/, zenodo.json)."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path(config).read_text()) if config and Path(config).exists() else {}
    if package_dir:
        pk = Path(package_dir)
        S = json.loads((pk / "site" / "summary.json").read_text())
        release = pk
        figures = figures or pk / "site" / "img"
    zj = Path(release) / "zenodo.json"   # output of argeement zenodo: the DOI of this release
    if zj.exists():
        z = json.loads(zj.read_text())
        cfg["zenodo_doi"] = z.get("doi") or cfg.get("zenodo_doi")
        if z.get("conceptdoi"):
            cfg["zenodo"] = f"https://doi.org/{z['conceptdoi']}"
    if not package_dir:
        S = collect(release)
    figs = []
    have = bool(figures and Path(figures).exists())
    macros = _load_macros(figures, Path(figures).parent if figures else None) if have else {}
    (out / "img").mkdir(exist_ok=True)
    for stem, kind, cap in FIGURES:
        src = ASSETS if kind == "asset" else (Path(figures) if have else None)
        if src is None:
            continue
        for ext in (".svg", ".png"):
            f = src / f"{stem}{ext}"
            if f.exists():
                shutil.copy(f, out / "img" / f.name)
                figs.append((out / "img" / f.name, _caption(cap, macros)))
                break
    (out / "index.html").write_text(render(S, cfg, figs), encoding="utf-8")
    (out / "summary.json").write_text(json.dumps(S, indent=2, default=str), encoding="utf-8")
    (out / ".nojekyll").write_text("")
    return out / "index.html"
