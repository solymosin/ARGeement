"""Real-genome use case: agreement of the gene calls of three annotation tools on the same genomes,
before mapping (by name) and after mapping (by the relations of an ARGeement release).

Input per genome: AMRFinderPlus (-n, TSV), RGI main (TXT) and ResFinder (JSON or
ResFinder_results_tab.txt) output. Calls are paired within a contig by reciprocal overlap
(locus-level agreement), and the references of paired calls are looked up in the release.
"""
import json
import re
from pathlib import Path

import pandas as pd

from .build import RELATION_RANK, STRONG_MAX_RANK

TOOLS = ("amrfinder", "rgi", "resfinder")
TOOL_DB = {"amrfinder": "amrfinder", "rgi": "card", "resfinder": "resfinder"}
COLS = ["genome", "tool", "contig", "start", "end", "name", "ref_db", "ref_id"]


def _col(df, *names):
    for n in names:
        if n in df.columns:
            return df[n]
    raise KeyError(f"egyik oszlop sincs meg: {names}")


# ------------------------------------------------------------------ reading tool outputs
def read_amrfinder(path, genome):
    """AMRFinderPlus TSV (3.x and 4.x column names). Acquired genes only (AMR/AMR and
    STRESS/BIOCIDE), without point mutations (POINT*) and HMM-only calls."""
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if d.empty:
        return pd.DataFrame(columns=COLS)
    typ, sub = _col(d, "Type", "Element type"), _col(d, "Subtype", "Element subtype")
    meth = _col(d, "Method")
    keep = (((typ == "AMR") & (sub == "AMR")) | ((typ == "STRESS") & (sub == "BIOCIDE"))) \
        & ~meth.str.startswith("POINT")
    d = d[keep]
    return pd.DataFrame(dict(
        genome=genome, tool="amrfinder", contig=_col(d, "Contig id", "Contig ID"),
        start=_col(d, "Start").astype(int), end=_col(d, "Stop", "End").astype(int),
        name=_col(d, "Element symbol", "Gene symbol"), ref_db="amrfinder",
        ref_id=_col(d, "Closest reference accession", "Accession of closest sequence")))


def read_rgi(path, genome, orf_map=None):
    """RGI main TXT, protein homolog models only (variant and overexpression models are not part
    of the release either). In protein mode (run_tools.sh) the ORF locations come from rgi_orfs.tsv
    (name, contig, start, end, strand); in contig mode the contig is the 'Contig' column without
    its _N suffix, and Start/Stop give the location."""
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if d.empty:
        return pd.DataFrame(columns=COLS)
    d = d[d["Model_type"] == "protein homolog model"]
    aro = d["ARO"].str.replace(r"^(ARO:)?", "ARO:", regex=True)
    if orf_map is not None and Path(orf_map).exists():
        om = pd.read_csv(orf_map, sep="\t", header=None, dtype=str,
                         names=["orf", "contig", "start", "end", "strand"]).set_index("orf")
        orf = d["ORF_ID"].str.split().str[0]
        contig, start, end = orf.map(om["contig"]), orf.map(om["start"]), orf.map(om["end"])
    else:
        contig = d["Contig"].str.replace(r"_\d+$", "", regex=True)
        start, end = d["Start"], d["Stop"]
    return pd.DataFrame(dict(
        genome=genome, tool="rgi", contig=contig, start=start.astype(int), end=end.astype(int),
        name=d["Best_Hit_ARO"], ref_db="card", ref_id=aro))


def read_resfinder(path, genome, entries=None):
    """ResFinder 4: the std_format JSON (seq_regions, 'ref_id' = the FASTA header, e.g.
    blaTEM-1B_1_AY458016), or ResFinder_results_tab.txt, which lacks the variant number, so the
    entry is looked up in the release by name and accession."""
    path = Path(path)
    if path.suffix == ".json":
        js = json.loads(path.read_text())
        rows = []
        for r in js.get("seq_regions", {}).values():
            dbs = " ".join(r.get("ref_database", []) if isinstance(r.get("ref_database"), list)
                           else [str(r.get("ref_database", ""))])
            rows.append(dict(genome=genome, tool="resfinder", contig=str(r["query_id"]).split()[0],
                             start=int(r["query_start_pos"]), end=int(r["query_end_pos"]),
                             name=r.get("name", ""),
                             ref_db="disinfinder" if "disinf" in dbs.lower() else "resfinder",
                             ref_id=r.get("ref_id", "")))
        return pd.DataFrame(rows, columns=COLS)
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if d.empty:
        return pd.DataFrame(columns=COLS)
    pos = d["Position in contig"].str.extract(r"(\d+)\.\.(\d+)").astype(int)
    ref = []
    for g, acc in zip(d["Resistance gene"], d["Accession no."]):
        hit = []
        if entries is not None:
            hit = entries.loc[(entries["db"].isin(["resfinder", "disinfinder"]))
                              & (entries["name"] == g)
                              & entries["entry_id"].str.endswith("_" + acc), ["db", "entry_id"]] \
                .values.tolist()
        ref.append(hit[0] if len(hit) == 1 else ["resfinder", ""])
    return pd.DataFrame(dict(
        genome=genome, tool="resfinder", contig=d["Contig"].str.split().str[0],
        start=pos.min(axis=1), end=pos.max(axis=1), name=d["Resistance gene"],
        ref_db=[r[0] for r in ref], ref_id=[r[1] for r in ref]))


# ------------------------------------------------------------------ merging the calls of one tool
SEP = "|"


def collapse_calls(calls, min_overlap=0.8):
    """Merges the calls of one tool at the same locus into a single call. ResFinder may report
    several equally good references at the same place (e.g. blaTEM-141, -206 and -214 for the
    same gene); these are one locus, the list of references is kept in 'ref_id'
    (separated by |), and the best one is used when the relation is evaluated."""
    out = []
    for (g, t, ctg), d in calls.groupby(["genome", "tool", "contig"], sort=False):
        d = d.sort_values(["start", "end"])
        groups = []
        for r in d.itertuples(index=False):
            for gr in groups:
                ov = min(r.end, gr["end"]) - max(r.start, gr["start"]) + 1
                if ov > 0 and ov >= min_overlap * (r.end - r.start + 1) \
                        and ov >= min_overlap * (gr["end"] - gr["start"] + 1):
                    gr["start"], gr["end"] = min(gr["start"], r.start), max(gr["end"], r.end)
                    gr["refs"].append((r.ref_db, r.ref_id, r.name))
                    break
            else:
                groups.append(dict(start=r.start, end=r.end, refs=[(r.ref_db, r.ref_id, r.name)]))
        for gr in groups:
            refs = list(dict.fromkeys(gr["refs"]))
            out.append(dict(genome=g, tool=t, contig=ctg, start=gr["start"], end=gr["end"],
                            name=SEP.join(dict.fromkeys(x[2] for x in refs)),
                            ref_db=SEP.join(x[0] for x in refs), ref_id=SEP.join(x[1] for x in refs),
                            n_refs=len(refs)))
    return pd.DataFrame(out, columns=COLS + ["n_refs"])


def _refs(r):
    return list(zip(str(r.ref_db).split(SEP), str(r.ref_id).split(SEP)))


# ------------------------------------------------------------------ pairing loci
def _overlap(a, b):
    ov = min(a.end, b.end) - max(a.start, b.start) + 1
    return max(0, ov)


def match_loci(ca, cb, min_overlap=0.8):
    """Pairs the calls of two tools within a genome and contig by reciprocal overlap (the overlap
    is at least min_overlap of both call lengths). Greedy, starting with the largest overlaps."""
    pairs = []
    for (g, ctg), xa in ca.groupby(["genome", "contig"]):
        xb = cb[(cb["genome"] == g) & (cb["contig"] == ctg)]
        cand = []
        for ia, ra in xa.iterrows():
            for ib, rb in xb.iterrows():
                ov = _overlap(ra, rb)
                if ov and ov >= min_overlap * (ra.end - ra.start + 1) \
                        and ov >= min_overlap * (rb.end - rb.start + 1):
                    cand.append((ov, ia, ib))
        used_a, used_b = set(), set()
        for ov, ia, ib in sorted(cand, reverse=True):
            if ia in used_a or ib in used_b:
                continue
            used_a.add(ia)
            used_b.add(ib)
            pairs.append((ia, ib))
    return pairs


# ------------------------------------------------------------------ relations in the release
class Release:
    def __init__(self, rel):
        rel = Path(rel)
        self.entries = pd.read_csv(rel / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
        m = pd.read_csv(rel / "mappings.tsv", sep="\t", keep_default_na=False,
                        dtype={"entry_a": str, "entry_b": str})
        m["_r"] = m["relation"].map(RELATION_RANK)
        self.rel = {}
        for r in m.sort_values(["_r", "n_diff"]).itertuples(index=False):
            self.rel.setdefault((r.db_a, r.entry_a, r.db_b, r.entry_b), (r.relation, int(r.n_diff)))
        usable = self.entries["usable"].astype(str) == "True"
        self.known = set(zip(self.entries["db"], self.entries["entry_id"]))
        self.usable = set(zip(self.entries.loc[usable, "db"], self.entries.loc[usable, "entry_id"]))
        # ORFs of split entries (a ResFinder call may point to a whole operon entry)
        orf = self.entries[(self.entries["parent"] != "") & usable]
        self.children = orf.groupby(["db", "parent"])["entry_id"].apply(list).to_dict()
        # whether there is any close counterpart in another database
        s = m[m["_r"] <= STRONG_MAX_RANK]
        self.strong_to = set(zip(s["db_a"], s["entry_a"], s["db_b"]))

    def _ids(self, db, eid):
        return self.children.get((db, eid), [eid])

    def relation(self, da, ea, db, eb):
        best = None
        for x in self._ids(da, ea):
            for y in self._ids(db, eb):
                r = self.rel.get((da, x, db, y)) or self.rel.get((db, y, da, x))
                if r and (best is None or RELATION_RANK[r[0]] < RELATION_RANK[best[0]]):
                    best = r
        return best

    def has_strong(self, da, ea, db):
        return any((da, x, db) in self.strong_to for x in self._ids(da, ea))

    def status(self, db, eid):
        if not eid:
            return "no_reference"
        ids = self._ids(db, eid)
        if not any((db, x) in self.known for x in ids):
            return "not_in_release"
        if not any((db, x) in self.usable for x in ids):
            return "not_usable"
        return "ok"


def _naive_equal(a, b):
    """Same name, ignoring letter case; with several references, any of them."""
    na = {x.strip().lower() for x in str(a).split(SEP)}
    nb = {x.strip().lower() for x in str(b).split(SEP)}
    return bool(na & nb)


def classify_pair(rel, ra, rb):
    """Classification of a paired locus after mapping (with several references, the best pair)."""
    pa = [x for x in _refs(ra) if rel.status(*x) == "ok"]
    pb = [x for x in _refs(rb) if rel.status(*x) == "ok"]
    if not pa or not pb:
        return "reference_missing", None
    r = None
    for da, ea in pa:
        for db, eb in pb:
            x = rel.relation(da, ea, db, eb)
            if x and (r is None or RELATION_RANK[x[0]] < RELATION_RANK[r[0]]):
                r = x
    if r is None:
        return "unrelated", None
    if r[0] == "identical_protein":
        return "identical", r
    if RELATION_RANK[r[0]] <= STRONG_MAX_RANK:
        return "strong", r
    if r[0] == "same_family":
        return "same_family", r
    return "unrelated", r


def classify_single(rel, r, other_db):
    """A locus called by one tool only: whether its reference has a close counterpart in the other
    database (then the missing call is due to the sensitivity or threshold of the other tool) or
    not (different database content)."""
    ok = [x for x in _refs(r) if rel.status(*x) == "ok"]
    if not ok:
        return "reference_missing"
    return "counterpart_exists" if any(rel.has_strong(d, e, other_db) for d, e in ok) \
        else "no_counterpart"


def run(release, calls, out, min_overlap=0.8):
    rel = Release(release)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    calls.reset_index(drop=True).to_csv(out / "calls_raw.tsv", sep="\t", index=False)
    calls = collapse_calls(calls, min_overlap).reset_index(drop=True)
    calls.to_csv(out / "calls.tsv", sep="\t", index=False)
    loci, summ = [], []
    for i, ta in enumerate(TOOLS):
        for tb in TOOLS[i + 1:]:
            ca, cb = calls[calls["tool"] == ta], calls[calls["tool"] == tb]
            pairs = match_loci(ca, cb, min_overlap)
            ma, mb = {p[0] for p in pairs}, {p[1] for p in pairs}
            for ia, ib in pairs:
                ra, rb = calls.loc[ia], calls.loc[ib]
                cat, r = classify_pair(rel, ra, rb)
                loci.append(dict(pair=f"{ta}-{tb}", genome=ra.genome, contig=ra.contig,
                                 start_a=ra.start, end_a=ra.end, name_a=ra["name"], ref_a=ra.ref_id,
                                 name_b=rb["name"], ref_b=rb.ref_id,
                                 n_refs_a=ra.n_refs, n_refs_b=rb.n_refs,
                                 naive_equal=_naive_equal(ra["name"], rb["name"]), category=cat,
                                 relation="" if r is None else r[0],
                                 n_diff="" if r is None else r[1]))
            for side, cc, matched, other in ((ta, ca, ma, TOOL_DB[tb]), (tb, cb, mb, TOOL_DB[ta])):
                for ix, r in cc[~cc.index.isin(matched)].iterrows():
                    loci.append(dict(pair=f"{ta}-{tb}", genome=r.genome, contig=r.contig,
                                     start_a=r.start, end_a=r.end,
                                     name_a=r["name"] if side == ta else "",
                                     ref_a=r.ref_id if side == ta else "",
                                     name_b=r["name"] if side == tb else "",
                                     ref_b=r.ref_id if side == tb else "",
                                     naive_equal=False,
                                     category=f"only_{side}:" + classify_single(rel, r, other),
                                     relation="", n_diff=""))
    L = pd.DataFrame(loci)
    L.to_csv(out / "loci.tsv", sep="\t", index=False)
    for pair, g in L.groupby("pair", sort=False):
        m = g[~g["category"].str.startswith("only_")]
        n_union = len(g)
        multi = ((m["n_refs_a"].astype(float) > 1) | (m["n_refs_b"].astype(float) > 1)).sum()
        row = dict(pair=pair, genomes=g["genome"].nunique(), loci_union=n_union, loci_matched=len(m),
                   loci_matched_multi_ref=int(multi),
                   locus_concordance=len(m) / n_union if n_union else 0,
                   name_concordance_before=m["naive_equal"].sum() / len(m) if len(m) else 0)
        for c in ("identical", "strong", "same_family", "unrelated", "reference_missing"):
            row[f"after_{c}"] = int((m["category"] == c).sum())
        row["name_concordance_after"] = (row["after_identical"] + row["after_strong"]) / len(m) \
            if len(m) else 0
        for c in sorted(set(g["category"]) - set(m["category"])):
            row[c] = int((g["category"] == c).sum())
        summ.append(row)
    S = pd.DataFrame(summ)
    S.to_csv(out / "summary.tsv", sep="\t", index=False)
    pairs = S.astype(object).where(S.notna(), None).to_dict(orient="records")   # NaN -> null
    (out / "summary.json").write_text(json.dumps(dict(
        release=Path(release).resolve().name, min_overlap=min_overlap,   # not an absolute path
        genomes=int(calls["genome"].nunique()), calls={t: int((calls["tool"] == t).sum()) for t in TOOLS},
        pairs=pairs), indent=2, default=int, allow_nan=False))
    import hashlib
    files = ["calls_raw.tsv", "calls.tsv", "loci.tsv", "summary.tsv", "summary.json"]
    (out / "SHA256SUMS").write_text("".join(
        f"{hashlib.sha256((out / f).read_bytes()).hexdigest()}  {f}\n" for f in files if (out / f).exists()))
    return S


def collect(release, indir):
    """Reads the outputs from a directory tree: <indir>/<genome>/{amrfinder.tsv, rgi.txt,
    resfinder.json or ResFinder_results_tab.txt}."""
    ents = pd.read_csv(Path(release) / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    parts = []
    for gdir in sorted(p for p in Path(indir).iterdir() if p.is_dir()):
        g = gdir.name
        if (f := gdir / "amrfinder.tsv").exists():
            parts.append(read_amrfinder(f, g))
        if (f := gdir / "rgi.txt").exists():
            parts.append(read_rgi(f, g, gdir / "rgi_orfs.tsv"))
        rf = next((f for f in (gdir / "resfinder.json", gdir / "ResFinder_results_tab.txt")
                   if f.exists()), None)
        if rf:
            parts.append(read_resfinder(rf, g, ents))
    parts = [p for p in parts if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS)
