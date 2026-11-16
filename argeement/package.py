"""Publication-ready package of a release (for Zenodo and the GitHub Pages summary).

Output (<release>/package/):
  argeement-<profile>-<run>.tar.gz   – the export directory, compressed reproducibly
  RELEASE_NOTES.md                   – description of the release: versions, counts, changes
  <archive>.sha256                   – checksum of the archive
  manifest.json, package.json        – for the Zenodo metadata (argeement zenodo --package)
  site/summary.json, site/img/, site/macros.json – for the web page (argeement site --package)

The package directory alone is sufficient for publication: it can be downloaded from the build
server and published to Zenodo and GitHub Pages from another machine.
"""
import gzip
import hashlib
import io
import json
import tarfile
from collections import Counter
from pathlib import Path

import pandas as pd

from .build import RELATION_RANK
from .export import LINKS, export

CHANGE_ORDER = ["entry_added", "entry_removed", "entry_sequence_changed", "entry_renamed",
                "mapping_added", "mapping_removed", "relation_changed", "family_changed"]


def _deterministic_targz(src_dir, archive, arcroot):
    """Byte-identical archive from identical content: sorted file order, fixed time and owner,
    gzip header without a timestamp, so that the checksum is reproducible as well."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for f in sorted(p for p in Path(src_dir).rglob("*") if p.is_file()):
            info = tarfile.TarInfo(f"{arcroot}/{f.relative_to(src_dir).as_posix()}")
            data = f.read_bytes()
            info.size, info.mtime, info.mode = len(data), 0, 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    with open(archive, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(buf.getvalue())


def _previous_release(rel):
    """Creation time and profile of the release the diff was computed against (diff_summary.md)."""
    f = rel / "diff" / "diff_summary.md"
    if f.exists():
        for line in f.read_text().splitlines():
            if line.startswith("Previous release:"):
                return line.split(":", 1)[1].strip()
    return None


def package(release, out=None, figures=None, no_diff=False):
    """no_diff: leave the changes against the previous build out of the notes and the web page,
    e.g. for the first published version, whose previous build was never published."""
    rel = Path(release)
    man = json.loads((rel / "manifest.json").read_text())
    exp = rel / "export"
    if not all((exp / f).exists() for f in ("SHA256SUMS", "schema.sql", "LICENSE.md")):
        export(rel)          # missing, or written by an earlier version
    out = Path(out) if out else rel / "package"
    out.mkdir(parents=True, exist_ok=True)

    run = rel.resolve().name
    profile = man.get("profile") or "release"
    tag = f"{profile}-{run}"
    name = f"argeement-{tag}"
    archive = out / f"{name}.tar.gz"
    _deterministic_targz(exp, archive, name)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (out / f"{archive.name}.sha256").write_text(f"{digest}  {archive.name}\n")

    # ------------------------------------------------------------ counts
    rd = lambda t: pd.read_csv(exp / f"{t}.tsv", sep="\t", dtype=str, keep_default_na=False)
    db_name = dict(zip(*rd("databases")[["db_id", "db"]].T.values))
    rel_name = dict(zip(*rd("relation_types")[["relation_id", "relation"]].T.values))
    e = rd("entries")
    e["db"] = e["db_id"].map(db_name)
    dbs = [d for d in LINKS if d in set(e["db"])]
    # entries: as in the source (ORF sub-entries counted separately), as in DATABASES.md and the paper
    per_db = e.groupby("db").agg(entries=("parent_eid", lambda s: int((s == "").sum())),
                                 orfs=("parent_eid", lambda s: int((s != "").sum())),
                                 usable=("pid", lambda s: int((s != "").sum())),
                                 proteins=("pid", lambda s: s[s != ""].nunique()))
    rel_counts = Counter(rd("relations")["relation_id"].map(rel_name))
    n_rel = sum(rel_counts.values())

    L = [f"# ARGeement mapping `{tag}`", "",
         f"Sequence-based mapping between antimicrobial resistance reference databases "
         f"(profile `{profile}`, created {man['created']}).", "",
         "| Database | Version | Entries | ORF sub-entries | Usable | Distinct proteins |",
         "|---|---|---|---|---|---|"]
    for d in dbs:
        row = per_db.loc[d]
        L.append(f"| {LINKS[d][0]} | `{man['database_versions'].get(d, '?')}` | {row.entries} | "
                 f"{row.orfs} | {row.usable} | {row.proteins} |")
    L += ["", "*Entries*: included entries of the source database. *ORF sub-entries*: genes obtained "
              "by splitting entries that are not a single coding sequence (e.g. operons). *Usable*: "
              "entries (including ORF sub-entries) with a protein sequence used in the mapping."]
    L += ["", f"**{e.loc[e['pid'] != '', 'pid'].nunique()}** distinct proteins, "
              f"**{n_rel}** relations between them: "
          + ", ".join(f"{k} {v}" for k, v in sorted(rel_counts.items(),
                                                     key=lambda kv: RELATION_RANK.get(kv[0], 9)))
          + ".", ""]

    diff_file = rel / "diff" / "diff.tsv"
    if diff_file.exists() and not no_diff:
        dc = Counter(pd.read_csv(diff_file, sep="\t", dtype=str, keep_default_na=False)["change"])
        prev = _previous_release(rel)
        L += ["## Changes since the previous release of this profile", ""]
        if prev:
            L += [f"Previous release: {prev}.", ""]
        L += [
              "| Change | Count |", "|---|---|"]
        L += [f"| {c} | {dc.get(c, 0)} |" for c in CHANGE_ORDER]
        L += [""]
    else:
        L += ["## Changes", "", "First published release of this profile.", ""]

    L += ["## Files", "",
          f"`{archive.name}` contains the tables of a relational layout linked by integer keys "
          "(`databases.tsv`, `relation_types.tsv`, `families.tsv`, `proteins.tsv`, `entries.tsv`, "
          "`relations.tsv`) with their definitions (`schema.sql`), `manifest.json`, "
          "`DATABASES.md` (content of the databases), `README.md` (schema, linking to the source "
          "databases, example code), `LICENSE.md` and `SHA256SUMS`.", "",
          "Licence of the data: CC-BY 4.0. The source databases have their own terms (commercial "
          "use of CARD itself requires a licence from McMaster University); see `LICENSE.md`.", "",
          f"SHA-256 of the archive: `{digest}`", "",
          "## Source databases", "",
          "Names and identifiers come from the source databases listed above; please cite them "
          "and observe their terms of use.", ""]
    (out / "RELEASE_NOTES.md").write_text("\n".join(L), encoding="utf-8")

    # companion files for publication (usable without the release directory)
    import shutil
    from .site import FIG_EN, collect
    (out / "manifest.json").write_text((rel / "manifest.json").read_text())
    (out / "package.json").write_text(json.dumps(dict(
        tag=tag, run=run, profile=profile, archive=archive.name, sha256=digest), indent=2))
    sd = out / "site"
    (sd / "img").mkdir(parents=True, exist_ok=True)
    summary = collect(rel)
    if no_diff:
        summary.pop("diff", None)
    (sd / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    if figures and Path(figures).exists():
        for stem in FIG_EN:
            f = Path(figures) / f"{stem}.png"
            if f.exists():
                shutil.copy(f, sd / "img" / f.name)
        # values for the figure captions of the web page (manuscript/make_assets.py --lang en)
        if (Path(figures) / "macros.json").exists():
            shutil.copy(Path(figures) / "macros.json", sd / "macros.json")
    return dict(tag=tag, archive=str(archive), size=archive.stat().st_size, sha256=digest,
                notes=str(out / "RELEASE_NOTES.md"))
