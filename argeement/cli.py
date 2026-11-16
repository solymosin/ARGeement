"""Command-line interface: argeement <command> (see argeement --help)."""
import argparse
import json
import csv
import sys
from pathlib import Path

import pandas as pd
import yaml

from .build import RELATION_RANK, build
from .diff import diff
from .export import export
from .package import package
from .query import query
from .validate import run_synthetic
from .report import report
import shutil

from .parsers import FAMRULE_FILES, HIERARCHY_FILES, PARSERS, first_existing, inventory_pointfinder

STD_COLS = ["db", "entry_id", "name", "family", "accession", "aro", "protein",
            "nucleotide", "protein_source", "issues", "note",
            "ncbi_node", "ncbi_family", "lineage", "classes", "mechanism", "source_group",
            "parent", "scope"]


def cmd_parse(a):
    kwargs = {}
    if a.db == "amrfinder" and a.config:
        cfg = yaml.safe_load(Path(a.config).read_text()).get("amrfinder", {})
        for k in ("types", "exclude_subtypes", "scopes", "include_subtypes"):
            if k in cfg:
                kwargs[k] = tuple(cfg[k])
    records, version, inventory = PARSERS[a.db](a.raw, **kwargs)
    if a.expect_version and a.expect_version not in str(version):
        sys.exit(f"ERROR: {a.db} version is '{version}', but the profile expects "
                 f"'{a.expect_version}'. Check the source (local/url/ref).")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STD_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in records:
            w.writerow({**r, "issues": ";".join(r["issues"])})
    out.with_suffix(".version").write_text(str(version) + "\n")
    out.with_suffix(".inventory").write_text(json.dumps(inventory, indent=2, ensure_ascii=False))
    if a.db == "amrfinder":   # the hierarchy and family rules are also copied into the release
        for names, suffix in ((HIERARCHY_FILES, ".ncbi_hierarchy"), (FAMRULE_FILES, ".ncbi_famrules")):
            src = first_existing(a.raw, names)
            if src:
                shutil.copy(src, out.with_suffix(suffix))
    print(f"{a.db}: {len(records)} entries, version: {version}", file=sys.stderr)


def cmd_build(a):
    params = {}
    if a.config:
        params = yaml.safe_load(Path(a.config).read_text()).get("params", {}) or {}
    build(a.std, a.out, params=params, mmseqs_bin=a.mmseqs, threads=a.threads,
          keep_tmp=a.keep_tmp, profile=a.profile)
    # 1. inventory-only databases (e.g. PointFinder): must be in place before the report
    for item in a.extra_inventory:
        name, _, path = item.partition("=")
        shutil.copy(path, Path(a.out) / f"inventory_{name}.json")
        # the inventory version also goes into the machine-readable manifest
        try:
            ver = json.loads(Path(path).read_text()).get("version")
            mf = Path(a.out) / "manifest.json"
            if ver and mf.exists():
                man = json.loads(mf.read_text())
                man.setdefault("database_versions", {})[name] = ver
                mf.write_text(json.dumps(man, indent=2, ensure_ascii=False))
        except (ValueError, OSError) as exc:
            print(f"[warning] the version of {name} was not added to the manifest: {exc}", file=sys.stderr)
    # 2. tartalmi riportok
    try:
        report(a.out, Path(a.out) / "DATABASES.md")
    except Exception as exc:   # a failed report is reported but does not discard the release
        print(f"[warning] DATABASES.md was not created: {exc}", file=sys.stderr)
    # 3. public export (includes the report)
    res = export(a.out)
    print(f"Export: {res['out']} ({res['rows']['relations']} relations instead of the "
          f"{res['rows']['mappings_full']} rows of mappings.tsv; "
          f"{res['size_export'] / 1e6:.1f} MB instead of {res['size_full'] / 1e6:.1f} MB)", file=sys.stderr)


def cmd_inventory(a):
    inv = inventory_pointfinder(a.raw)
    if a.expect_version and a.expect_version not in inv["version"]:
        sys.exit(f"ERROR: pointfinder version is '{inv['version']}', but the profile expects "
                 f"'{a.expect_version}'. Check the source (local/url/ref).")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(inv, indent=2, ensure_ascii=False))
    t = inv["totals"]
    print(f"pointfinder: {t['species']} species, {t['mutation_positions']} mutation positions, "
          f"version: {inv['version']}", file=sys.stderr)


def cmd_export(a):
    res = export(a.release, a.out)
    print(f"Export: {res['out']}", file=sys.stderr)


def cmd_package(a):
    res = package(a.release, a.out, figures=a.figures, no_diff=a.no_diff)
    print(f"Package: {res['archive']} ({res['size'] / 1e6:.1f} MB)\n"
          f"Release notes: {res['notes']}\n"
          f"SHA-256: {res['sha256']}\n"
          f"Publish the package directory with publish_local.sh (Zenodo + GitHub Pages).", file=sys.stderr)


def cmd_query(a):
    res, man = query(a.release, a.fasta, mmseqs_bin=a.mmseqs, threads=a.threads, top=a.top)
    vers = ", ".join(f"{k} {v}" for k, v in sorted(man.get("database_versions", {}).items()))
    prof = f"[{man['profile']}] " if man.get("profile") else ""
    print(f"# {prof}release: {Path(a.release).resolve().name}; versions: {vers}", file=sys.stderr)
    res.to_csv(a.out if a.out else sys.stdout, sep="\t", index=False)


def cmd_validate(a):
    summ, meta = run_synthetic(a.release, a.out, n=a.n, seed=a.seed, mmseqs_bin=a.mmseqs,
                               threads=a.threads)
    print(summ.to_string(index=False), file=sys.stderr)
    print(f"Total: {meta['n_cases']} cases from {meta['n_proteins']} proteins; "
          f"truth found {meta['overall']['truth_found']:.1%}, "
          f"relation type correct {meta['overall']['relation_ok']:.1%}, "
          f"number of differences correct {meta['overall']['ndiff_ok']:.1%}", file=sys.stderr)


def cmd_concordance(a):
    from .concordance import collect, run
    calls = collect(a.release, a.input)
    if calls.empty:
        sys.exit("no readable tool output found")
    S = run(a.release, calls, a.out or Path(a.release) / "concordance", a.min_overlap)
    cols = ["pair", "loci_union", "loci_matched", "locus_concordance", "name_concordance_before",
            "name_concordance_after"]
    print(S[cols].to_string(index=False), file=sys.stderr)


def cmd_zenodo(a):
    from .zenodo import publish
    if not (a.release or a.package):
        sys.exit("--release or --package is required")
    publish(a.release, a.config, sandbox=a.sandbox, dry_run=a.dry_run, package_dir=a.package,
            new_record=a.new)


def cmd_site(a):
    from .site import build_site
    if not (a.release or a.package):
        sys.exit("--release or --package is required")
    p = build_site(a.release, a.out, a.figures, a.config, package_dir=a.package)
    print(f"Web page: {p}", file=sys.stderr)


def cmd_report(a):
    out = a.out or str(Path(a.release) / "DATABASES.md")
    report(a.release, out)
    print(f"Riport: {out}", file=sys.stderr)


def cmd_diff(a):
    counts = diff(a.old, a.new, a.out, all_relations=a.all_relations)
    for k, v in counts.items():
        print(f"{k}\t{v}")


def family_level(q, from_db, src_entries, entries, sub, hier, n_examples=5):
    """Family-level resolution: if the query is a node of the NCBI hierarchy (e.g. AMRFinderPlus
    reports 'blaTEM' for a new allele), the mappings of all entries below the node are
    summarised, one row per target cluster."""
    if hier is None or from_db != "amrfinder" or "lineage" not in src_entries.columns:
        return None
    nodes = set(hier.loc[(hier["node_id"] == q) | (hier.get("symbol", "") == q), "node_id"])
    if not nodes:
        return None
    members = src_entries[src_entries["lineage"].apply(
        lambda l: bool(nodes & set(l.split(">"))) if l else False)]
    label = f"{'/'.join(sorted(nodes))} ({len(members)} tag)"
    if members.empty:
        return [dict(query=q, entry_a=label, status="family_node_without_entries")]
    m = sub[sub["entry_a"].isin(set(members["entry_id"]))]
    if m.empty:
        return [dict(query=q, entry_a=label, status="family_level_no_mapping")]
    tgt = entries.set_index(["db", "entry_id"])
    m = m.assign(fc=[tgt.at[(d, e), "family_cluster"] for d, e in zip(m["db_b"], m["entry_b"])],
                 fam_b=[tgt.at[(d, e), "family"] for d, e in zip(m["db_b"], m["entry_b"])])
    rows = []
    for fc, g in sorted(m.groupby("fc"), key=lambda kv: kv[1]["_r"].min()):
        names = sorted(set(g["name_b"]))
        fams = sorted(set(f for f in g["fam_b"] if f))
        shown = ", ".join(names[:n_examples]) + (f" (+{len(names) - n_examples})" if len(names) > n_examples else "")
        best = g.sort_values(["_r", "n_diff"]).iloc[0]
        rows.append(dict(query=q, entry_a=label, status="family_level", entry_b=fc,
                         name_b=(f"[{'; '.join(fams[:2])}] " if fams else "") + shown,
                         relation=best["relation"], n_diff=best["n_diff"]))
    return rows


def cmd_lookup(a):
    rel = Path(a.release)
    man_file = rel / "manifest.json"
    if man_file.exists():
        man = json.loads(man_file.read_text())
        vers = ", ".join(f"{k} {v}" for k, v in sorted(man.get("database_versions", {}).items()))
        prof = f"[{man['profile']}] " if man.get("profile") else ""
        print(f"# {prof}release: {rel.resolve().name}; versions: {vers}", file=sys.stderr)
    entries = pd.read_csv(rel / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    maps = pd.read_csv(rel / "mappings.tsv", sep="\t", keep_default_na=False,
                       dtype={"entry_a": str, "entry_b": str})
    if a.names:
        queries = [q.strip() for q in a.names.split(",") if q.strip()]
    else:
        src = pd.read_csv(a.input, sep="\t", dtype=str, keep_default_na=False)
        queries = src[a.column].unique().tolist()

    src_entries = entries[entries["db"] == a.from_db]
    hier = None
    if (rel / "ncbi_hierarchy.tsv").exists():
        hier = pd.read_csv(rel / "ncbi_hierarchy.tsv", sep="\t", dtype=str, keep_default_na=False)
        hier.columns = [c.lstrip("#") for c in hier.columns]
    max_rank = RELATION_RANK[a.min_relation]
    sub = maps[(maps["db_a"] == a.from_db) & (maps["db_b"] == a.to_db)].copy()
    sub["_r"] = sub["relation"].map(RELATION_RANK)
    sub = sub[sub["_r"] <= max_rank]
    if a.best_only:
        sub = sub[sub["best_for_a"].astype(str) == "True"]

    out = []
    for q in queries:
        # exact identifier or name match within the source database
        hit = src_entries[(src_entries["entry_id"] == q) | (src_entries["name"] == q)]
        if hit.empty:
            hit = src_entries[src_entries["name"].str.lower() == q.lower()]
        if hit.empty:
            fam_rows = family_level(q, a.from_db, src_entries, entries, sub, hier)
            out.extend(fam_rows or [dict(query=q, status="not_found_in_source")])
            continue
        for e in hit["entry_id"]:
            m = sub[sub["entry_a"] == e].sort_values(["_r", "n_diff"])
            if m.empty:
                out.append(dict(query=q, entry_a=e, status="no_mapping"))
            for r in m.itertuples(index=False):
                out.append(dict(query=q, entry_a=e, status="mapped", entry_b=r.entry_b,
                                name_b=r.name_b, relation=r.relation, n_diff=r.n_diff, pident=r.pident,
                                cov_a=r.cov_a, cov_b=r.cov_b, rbh=r.rbh,
                                n_equally_good_b=r.n_equally_good_b,
                                tie_scope=r.tie_scope))
    cols = ["query", "entry_a", "status", "entry_b", "name_b", "relation", "n_diff", "pident",
            "cov_a", "cov_b", "rbh", "n_equally_good_b", "tie_scope"]
    res = pd.DataFrame(out, columns=cols)
    for c in ("n_equally_good_b", "n_diff"):
        res[c] = pd.to_numeric(res[c]).astype("Int64")
    res.to_csv(sys.stdout, sep="\t", index=False)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="argeement", description=__doc__)
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("parse", help="read one database into the standard format")
    p.add_argument("db", choices=sorted(PARSERS))
    p.add_argument("--raw", required=True, help="directory of the downloaded database")
    p.add_argument("--out", required=True, help="kimeneti TSV")
    p.add_argument("--config")
    p.add_argument("--expect-version", default="",
                   help="stop if the parsed version does not contain this text")
    p.set_defaults(func=cmd_parse)

    p = sp.add_parser("build", help="build the mapping (a release)")
    p.add_argument("--std", nargs="+", required=True, help="TSV files of the parse step")
    p.add_argument("--out", required=True, help="release directory")
    p.add_argument("--config")
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--keep-tmp", action="store_true")
    p.add_argument("--profile", default="", help="profile name (recorded in the manifest)")
    p.add_argument("--extra-inventory", action="append", default=[], metavar="NAME=FILE",
                   help="inventory-only database (e.g. pointfinder=...json) for the report")
    p.set_defaults(func=cmd_build)

    p = sp.add_parser("diff", help="compare two releases")
    p.add_argument("old")
    p.add_argument("new")
    p.add_argument("--out", required=True)
    p.add_argument("--all-relations", action="store_true",
                   help="also compare weak (same_family, homolog) relations")
    p.set_defaults(func=cmd_diff)

    p = sp.add_parser("inventory", help="inventory only, not part of the mapping: PointFinder")
    p.add_argument("db", choices=["pointfinder"])
    p.add_argument("--raw", required=True)
    p.add_argument("--out", required=True, help="output JSON")
    p.add_argument("--expect-version", default="")
    p.set_defaults(func=cmd_inventory)

    p = sp.add_parser("export", help="non-redundant tables of a release for publication")
    p.add_argument("--release", required=True)
    p.add_argument("--out", help="output directory (default: <release>/export)")
    p.set_defaults(func=cmd_export)

    p = sp.add_parser("query", help="closest entries per database for external sequences")
    p.add_argument("--release", required=True)
    p.add_argument("--fasta", required=True, help="protein or nucleotide sequences (FASTA)")
    p.add_argument("--out", help="output TSV (default: stdout)")
    p.add_argument("--top", type=int, default=1,
                   help="1: best only (all if tied); N>1: the N best per database")
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_query)

    p = sp.add_parser("validate", help="synthetic validation with variants of known modifications")
    p.add_argument("--release", required=True)
    p.add_argument("--out", help="output directory (default: <release>/validation_synthetic)")
    p.add_argument("--n", type=int, default=200, help="number of curated proteins sampled")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--mmseqs", default="mmseqs")
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_validate)

    p = sp.add_parser("zenodo", help="publish a release on Zenodo (new version of the record)")
    p.add_argument("--release", help="release directory (packaged here first)")
    p.add_argument("--package", help="finished package directory (e.g. downloaded from the server)")
    p.add_argument("--config", default="zenodo.yaml", help="metadata and id of the first record")
    p.add_argument("--sandbox", action="store_true", help="publish to sandbox.zenodo.org (testing)")
    p.add_argument("--dry-run", action="store_true", help="only print what would be uploaded")
    p.add_argument("--new", action="store_true",
                   help="create the first version of a new record (only when zenodo.yaml has no id)")
    p.set_defaults(func=cmd_zenodo)

    p = sp.add_parser("site", help="summary web page of a release (GitHub Pages)")
    p.add_argument("--release", help="release directory")
    p.add_argument("--package", help="downloaded package directory (site/summary.json, site/img/)")
    p.add_argument("--out", default="docs")
    p.add_argument("--figures", help="directory of PNG figures")
    p.add_argument("--config", default="site.yaml", help="title, GitHub and Zenodo links")
    p.set_defaults(func=cmd_site)

    p = sp.add_parser("concordance",
                      help="agreement of annotations of real genomes before and after mapping")
    p.add_argument("--release", required=True)
    p.add_argument("--input", required=True,
                   help="<genom>/{amrfinder.tsv, rgi.txt, resfinder.json|ResFinder_results_tab.txt}")
    p.add_argument("--out", help="output directory (default: <release>/concordance)")
    p.add_argument("--min-overlap", type=float, default=0.8,
                   help="reciprocal overlap threshold for pairing loci")
    p.set_defaults(func=cmd_concordance)

    p = sp.add_parser("package", help="publication package of a release (archive, notes, web page data)")
    p.add_argument("--release", required=True)
    p.add_argument("--out", help="output directory (default: <release>/package)")
    p.add_argument("--figures", help="directory of PNG figures for the web page")
    p.add_argument("--no-diff", action="store_true",
                   help="omit the changes against the previous build (e.g. first published version)")
    p.set_defaults(func=cmd_package)

    p = sp.add_parser("report", help="database content report (Markdown) of a release")
    p.add_argument("--release", required=True)
    p.add_argument("--out", help="output .md (default: in the release directory)")
    p.set_defaults(func=cmd_report)

    p = sp.add_parser("lookup", help="translate gene names from one database to another")
    p.add_argument("--release", required=True)
    p.add_argument("--from", dest="from_db", required=True, choices=sorted(PARSERS))
    p.add_argument("--to", dest="to_db", required=True, choices=sorted(PARSERS))
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--names", help="comma-separated names or identifiers")
    g.add_argument("--input", help="TSV, e.g. the output of an annotation tool")
    p.add_argument("--column", default="gene", help="name column of the --input file")
    p.add_argument("--min-relation", default="length_variant", choices=list(RELATION_RANK))
    p.add_argument("--best-only", action="store_true")
    p.set_defaults(func=cmd_lookup)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
