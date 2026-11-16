"""Building the mapping: hash layer, MMseqs2 layer, relations, clusters."""
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml

from . import __version__
from .sequtil import BLOCKING_ISSUES, _best_tiling, n_diff_from_alignment, orf_candidates, seq_hash

CLASS_MAP_FILE = Path(__file__).with_name("class_map.yaml")
EXTRA_COLS = ("ncbi_node", "ncbi_family", "lineage", "classes", "mechanism", "source_group",
              "parent", "scope")


def load_class_map(path=CLASS_MAP_FILE):
    c = yaml.safe_load(Path(path).read_text())
    term = {t.lower(): cls for cls, ts in c["terms"].items() for t in ts}
    for t in c.get("ignore", []):
        term[t.lower()] = None
    return term, {a.lower() for a in c.get("ambiguous", [])}


def harmonize(classes, source_group, term, ambiguous):
    """Native class terms -> common classes. Returns (common classes, unmapped terms).
    If a term is ambiguous or nothing can be mapped, the source file name (ResFinder)
    decides, if it is a known term."""
    out, unmapped, fallback = set(), set(), False
    for c in filter(None, (x.strip() for x in classes.split(";"))):
        k = c.lower()
        if k in term:
            if term[k]:            # None = ignored term
                out.add(term[k])
        elif k in ambiguous:
            fallback = True
        else:
            unmapped.add(c)
    if (fallback or not out) and source_group and term.get(source_group.lower()):
        out.add(term[source_group.lower()])
    return sorted(out), sorted(unmapped)


RELATION_RANK = {"identical_protein": 0, "near_identical": 1, "length_variant": 2,
                 "same_family": 3, "homolog": 4}
STRONG_MAX_RANK = 2  # up to here a relation is "close": all are kept, and the diff tracks them

DEFAULT_PARAMS = {
    # relation type thresholds (identity %, and the SMALLER of the two coverages)
    "near_identical_pident": 98.0, "near_identical_cov": 0.95,
    "same_family_pident": 90.0, "same_family_cov": 0.80,
    # search threshold: no relation is recorded below it
    "search_min_seq_id": 0.50, "search_min_cov": 0.50,
    "search_sensitivity": 7.5,
    # upper limit of hits per query: for large families (e.g. blaPDC) 1000 truncated
    # close hits, so the limit is effectively disabled
    "search_max_seqs": 20000,
    # number of best weak relations kept per entry and target database
    "top_k_weak": 3,
}


# ------------------------------------------------------------ loading
def load_standardized(std_files):
    frames, versions = [], {}
    for f in std_files:
        df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        frames.append(df)
        vfile = Path(f).with_suffix(".version")
        db = df["db"].iloc[0] if len(df) else Path(f).stem
        versions[db] = vfile.read_text().strip() if vfile.exists() else "unknown"
    entries = pd.concat(frames, ignore_index=True)
    for col in EXTRA_COLS:   # older standardised files
        if col not in entries.columns:
            entries[col] = ""
    entries[list(EXTRA_COLS)] = entries[list(EXTRA_COLS)].fillna("")
    entries["key"] = entries["db"] + "::" + entries["entry_id"]
    if entries["key"].duplicated().any():
        dups = entries.loc[entries["key"].duplicated(), "key"].tolist()[:5]
        raise ValueError(f"Duplicate entry identifiers: {dups}")
    return entries, versions


# ------------------------------------------------------------ MMseqs2
def run_mmseqs(fasta, workdir, params, mmseqs_bin, threads):
    out = Path(workdir) / "hits.m8"
    tmp = Path(workdir) / "mmseqs_tmp"
    cols = "query,target,pident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits"
    cmd = [mmseqs_bin, "easy-search", str(fasta), str(fasta), str(out), str(tmp),
           "--min-seq-id", str(params["search_min_seq_id"]),
           "-c", str(params["search_min_cov"]), "--cov-mode", "0",
           "--seq-id-mode", "0",   # identity relative to the alignment length (the default, made explicit)
           "-s", str(params["search_sensitivity"]),
           "--max-seqs", str(params["search_max_seqs"]),
           "--threads", str(threads), "--format-output", cols, "-v", "1"]
    subprocess.run(cmd, check=True)
    hits = pd.read_csv(out, sep="\t", names=cols.split(","))
    hits = hits[hits["query"] != hits["target"]].copy()
    # number of identical positions: MMseqs2's 'nident' field is 0 with default settings,
    # so it is computed back from pident (3 decimals; denominator: alignment length)
    hits["nident"] = (hits["alnlen"] * hits["pident"] / 100).round().astype(int)
    hits["pident"] = hits["pident"].round(2)
    hits["qcov"] = ((hits["qend"] - hits["qstart"] + 1) / hits["qlen"]).round(3)
    hits["tcov"] = ((hits["tend"] - hits["tstart"] + 1) / hits["tlen"]).round(3)
    # number of differences: mismatch or gap columns in the alignment + the unaligned ends
    hits["n_diff"] = n_diff_from_alignment(hits["alnlen"], hits["nident"], hits["qstart"],
                                           hits["qend"], hits["qlen"], hits["tstart"],
                                           hits["tend"], hits["tlen"])
    # best HSP per pair; the full sort key makes the choice independent of the output order of
    # MMseqs2, which varies with the threads and the composition of the database
    hits = (hits.sort_values(["query", "target", "bits", "n_diff", "alnlen", "qstart", "tstart"],
                             ascending=[True, True, False, True, False, True, True], kind="mergesort")
                .drop_duplicates(["query", "target"]))
    return hits


def symmetrize_hits(hits):
    """Assigns every protein pair once.

    The identity and coverage reported by MMseqs2 depend on the search direction, so a pair near a
    threshold could be a close relation from one direction and a weak one from the other, or be
    found in one direction only. Here a single alignment is chosen for every unordered pair (the
    one with fewer differences, on a tie the higher bit score, as in the export) and used in both
    directions with the query- and target-specific fields swapped. Classification, coverage
    summary, family clusters and the export are thus based on the same value per pair.

    If the two directions have the same number of differences and the same (integer) bit score,
    the longer alignment is kept, then the one whose query comes first. Such ties are common: a
    substitution at the last position whose score is zero is included in the local alignment in one
    direction and left out in the other, and the identity of the shorter alignment is higher
    (e.g. 98.1% instead of 97.9% with 8 differences), which can move the pair across the
    near-identical threshold. Without a full sort key the choice depended on the order of the
    MMseqs2 output, so the same pair could be classified differently in two builds.
    """
    if hits.empty:
        return hits
    h = hits.copy()
    a = h[["query", "target"]].min(axis=1)
    b = h[["query", "target"]].max(axis=1)
    h["_a"], h["_b"] = a, b
    h = h.sort_values(["_a", "_b", "n_diff", "bits", "alnlen", "query"],
                      ascending=[True, True, True, False, False, True], kind="mergesort") \
         .drop_duplicates(["_a", "_b"])
    h = h.drop(columns=["_a", "_b"])
    swap = {"query": "target", "target": "query", "qstart": "tstart", "tstart": "qstart",
            "qend": "tend", "tend": "qend", "qlen": "tlen", "tlen": "qlen",
            "qcov": "tcov", "tcov": "qcov"}
    rev = h.rename(columns=swap)[h.columns]
    out = pd.concat([h, rev], ignore_index=True)
    return out.sort_values(["query", "target"]).reset_index(drop=True)


# ------------------------------------------------------------ ORF splitting with curated priority
SUPPORT_BONUS = 100_000   # larger than any possible coverage: a curated candidate always wins
MIN_ORF_AA = 60


def _blocked(issues):
    return bool(set(filter(None, issues.split(";"))) & BLOCKING_ISSUES)


def rescue_split_orfs(entries, p, mmseqs_bin, threads, workdir, curated=("amrfinder", "card")):
    """Repeats the ORF selection of split entries so that candidates related to a curated
    (AMRFinderPlus, CARD) protein take priority.

    At parsing, weighted interval scheduling only considers the covered length, so a longer,
    overlapping foreign ORF could displace the real gene. Here all ORF candidates of the entry are
    grouped by stop codon, the longest member of each group is searched against the curated
    proteins, and candidates of related groups receive a large bonus. The remaining candidates
    fill the remaining space; their support is decided later by orf_support. If an entry has no
    candidate related to a curated protein, or the selection does not change, the ORFs from
    parsing are kept. Returns (entries, keys of the changed entries).
    """
    split = entries[entries["issues"].str.contains("split_into_orfs", regex=False)
                    & (entries["nucleotide"] != "")]
    cur = entries[entries["db"].isin(curated) & (entries["protein"] != "")
                  & ~entries["issues"].map(_blocked)].drop_duplicates("protein")
    if split.empty or cur.empty:
        return entries, []
    cands, reps = {}, {}
    for idx, r in split.iterrows():
        cs = orf_candidates(r["nucleotide"], MIN_ORF_AA)
        for c in cs:
            c["group"] = f"g{idx}_{c['stop_key']}"
            if len(c["protein"]) > len(reps.get(c["group"], "")):
                reps[c["group"]] = c["protein"]
        cands[idx] = cs
    if not reps:
        return entries, []
    work = Path(workdir)
    qfa, tfa, out = work / "orf_reps.faa", work / "curated.faa", work / "orf_reps.m8"
    qfa.write_text("".join(f">{g}\n{aa}\n" for g, aa in reps.items() if aa))
    tfa.write_text("".join(f">c{i}\n{aa}\n" for i, aa in enumerate(cur["protein"])))
    subprocess.run([mmseqs_bin, "easy-search", str(qfa), str(tfa), str(out), str(work / "orf_tmp"),
                    "--min-seq-id", str(p["search_min_seq_id"]), "-c", str(p["search_min_cov"]),
                    "--cov-mode", "0", "--seq-id-mode", "0", "-s", str(p["search_sensitivity"]),
                    "--max-seqs", "100", "--threads", str(threads),
                    "--format-output", "query,target", "-v", "1"], check=True)
    supported = set(pd.read_csv(out, sep="\t", names=["query", "target"])["query"]) \
        if out.stat().st_size else set()

    coord_rx = r"strand ([+-]), nt (\d+)-(\d+)"
    drop, new_rows, changed = [], [], []
    for idx, cs in cands.items():
        if not any(c["group"] in supported for c in cs):
            continue
        for c in cs:
            c["bonus"] = SUPPORT_BONUS if c["group"] in supported else 0
        tiling = _best_tiling(cs)
        par = entries.loc[idx]
        old = entries[(entries["db"] == par["db"]) & (entries["parent"] == par["entry_id"])]
        old_c = sorted(tuple(x) for x in old["note"].str.extract(coord_rx).values.tolist())
        new_c = sorted((o["strand"], str(o["start"] + 1), str(o["end"])) for o in tiling)
        if not len(old) or old_c == new_c:
            continue
        tmpl = old.iloc[0].to_dict()
        for k, o in enumerate(tiling, 1):
            row = dict(tmpl)
            row.update(entry_id=f"{par['entry_id']}#orf{k}", protein=o["protein"], issues="",
                       note=f"ORF {k}/{len(tiling)}, strand {o['strand']}, nt {o['start'] + 1}-{o['end']}"
                            f" (curated-priority)")
            row["key"] = f"{row['db']}::{row['entry_id']}"
            new_rows.append(row)
        drop += list(old.index)
        entries.loc[idx, "note"] = f"split into {len(tiling)} ORF(s): {par['entry_id']}#orf1.."
        changed.append(par["key"])
    if changed:
        entries = pd.concat([entries.drop(index=drop), pd.DataFrame(new_rows)], ignore_index=True)
    return entries, changed


def mmseqs_version(mmseqs_bin):
    try:
        return subprocess.run([mmseqs_bin, "version"], capture_output=True,
                              text=True).stdout.strip()
    except OSError:
        return "unknown"


# ------------------------------------------------------------ kapcsolatok
def classify(pident, qcov, tcov, p):
    cov_min, cov_max = min(qcov, tcov), max(qcov, tcov)
    if pident >= p["near_identical_pident"]:
        if cov_min >= p["near_identical_cov"]:
            return "near_identical"
        # one sequence is almost entirely contained in the other: different reference length
        # (truncated or extended ORF, different start codon)
        if cov_max >= p["near_identical_cov"] and cov_min >= p["same_family_cov"]:
            return "length_variant"
    if pident >= p["same_family_pident"] and cov_min >= p["same_family_cov"]:
        return "same_family"
    return "homolog"


def build_candidates(entries, hits, p):
    """At the level of protein ids (P-...): source P -> target database -> list of candidates."""
    dbs_of_pid = entries[entries["usable"]].groupby("prot_id")["db"].apply(set).to_dict()
    hit_lookup = defaultdict(list)
    for h in hits.itertuples(index=False):
        hit_lookup[h.query].append(h)

    cands = {}  # (pid, target_db) -> list of dict
    for pid, dbs in dbs_of_pid.items():
        per_db = defaultdict(list)
        for tdb in dbs:  # identical sequence in another database
            per_db[tdb].append(dict(target=pid, relation="identical_protein",
                                    pident=100.0, cov_a=1.0, cov_b=1.0, n_diff=0,
                                    bits=float("nan")))
        for h in hit_lookup.get(pid, []):
            rel = classify(h.pident, h.qcov, h.tcov, p)
            for tdb in dbs_of_pid.get(h.target, ()):
                per_db[tdb].append(dict(target=h.target, relation=rel, pident=h.pident,
                                        cov_a=h.qcov, cov_b=h.tcov, n_diff=int(h.n_diff),
                                        bits=h.bits))
        for tdb, lst in per_db.items():
            lst.sort(key=lambda c: (RELATION_RANK[c["relation"]], c["n_diff"],
                                    -(0 if pd.isna(c["bits"]) else c["bits"])))
            strong = [c for c in lst if RELATION_RANK[c["relation"]] <= STRONG_MAX_RANK]
            weak_all = [c for c in lst if RELATION_RANK[c["relation"]] > STRONG_MAX_RANK]
            k = p["top_k_weak"]
            if len(weak_all) > k:
                # candidates tied at the k-th place (same type and number of differences) are all
                # kept, so the truncation does not depend on the order
                cut = (RELATION_RANK[weak_all[k - 1]["relation"]], weak_all[k - 1]["n_diff"])
                weak = [c for c in weak_all
                        if (RELATION_RANK[c["relation"]], c["n_diff"]) <= cut]
            else:
                weak = weak_all
            kept = strong + weak
            best = kept[0]
            for c in kept:
                # tie: same relation type AND exactly the same number of differences
                c["is_best"] = (c["relation"] == best["relation"] and
                                c["n_diff"] == best["n_diff"])
            n_best = sum(c["is_best"] for c in kept)
            for c in kept:
                c["n_equally_good"] = n_best
            cands[(pid, tdb)] = kept
    return cands


def build_mappings(entries, cands, fam, labels=None):
    by_pid_db = entries[entries["usable"]].groupby(["prot_id", "db"])
    members = {k: g for k, g in by_pid_db}
    dbs_of_pid = defaultdict(set)
    for (p_, d) in members:
        dbs_of_pid[p_].add(d)
    best_targets = {k: {c["target"] for c in v if c["is_best"]} for k, v in cands.items()}
    # scope of a tie: none = unambiguous; within_family = the best targets share one
    # family cluster (ambiguous at allele level only); cross_family = different
    # clusters. By construction the latter can only occur for 'homolog' relations.
    def scope(targets):
        if len(targets) <= 1:
            return "none"
        labs = [(labels or {}).get(t, (None, ""))[0] for t in targets]
        keys = labs if all(labs) else [fam.get(t, t) for t in targets]
        return "within_family" if len(set(keys)) == 1 else "cross_family"
    tie_scope = {k: scope(v) for k, v in best_targets.items()}

    rows = []
    for (pid, tdb), lst in cands.items():
        for sdb in dbs_of_pid[pid]:
            if sdb == tdb:
                continue
            src = members[(pid, sdb)]
            for c in lst:
                tgt = members[(c["target"], tdb)]
                rbh = c["is_best"] and pid in best_targets.get((c["target"], sdb), set())
                for a in src.itertuples(index=False):
                    for b in tgt.itertuples(index=False):
                        nt_ident = ("na" if c["relation"] != "identical_protein" or
                                    not a.nt_id or not b.nt_id
                                    else ("yes" if a.nt_id == b.nt_id else "no"))
                        rows.append(dict(
                            db_a=sdb, entry_a=a.entry_id, name_a=a.name,
                            db_b=tdb, entry_b=b.entry_id, name_b=b.name,
                            relation=c["relation"], n_diff=c["n_diff"], pident=c["pident"],
                            cov_a=c["cov_a"], cov_b=c["cov_b"], bits=c["bits"],
                            nt_identical=nt_ident, best_for_a=c["is_best"], rbh=rbh,
                            n_equally_good_b=c["n_equally_good"],
                            tie_scope=tie_scope[(pid, tdb)],
                            prot_a=pid, prot_b=c["target"]))
    df = pd.DataFrame(rows)
    if len(df):
        df["_r"] = df["relation"].map(RELATION_RANK)
        df = (df.sort_values(["db_a", "entry_a", "db_b", "_r", "n_diff"])
                .drop(columns="_r").reset_index(drop=True))
    return df


# ------------------------------------------------------------ klaszterek
def orf_support(entries, hits, curated=("amrfinder", "card"), sole_cover=0.8, operon_cover=0.9):
    """Filtering of the sub-entries obtained by ORF splitting.

    An ORF is kept if its protein is identical to a curated (AMRFinderPlus or CARD) protein or
    related to one above the search threshold. If no ORF of an entry qualifies, the splitting
    pattern decides (see below). Excluded ORFs are flagged 'orf_unsupported'.
    """
    if "parent" not in entries.columns:
        return entries, hits, 0
    is_orf = entries["usable"] & (entries["parent"].fillna("") != "")
    if not is_orf.any():
        return entries, hits, 0
    ref = set(entries.loc[entries["usable"] & entries["db"].isin(curated), "prot_id"])
    linked = set(hits.loc[hits["target"].isin(ref), "query"]) | set(hits.loc[hits["query"].isin(ref), "target"])
    supported = is_orf & (entries["prot_id"].isin(ref) | entries["prot_id"].isin(linked))

    # if no ORF of an entry is supported, the splitting pattern decides:
    #  - a single ORF covering at least sole_cover of the entry: the gene itself
    #  - several ORFs on the same strand covering at least operon_cover together: a unique operon
    #    (e.g. sitABCD, whose genes are absent from the curated databases)
    # if at least one ORF is supported, the unsupported ones (neighbouring, antisense or
    # fragmentary ORFs) are excluded
    nt_len = dict(zip(entries["entry_id"], entries["nucleotide"].fillna("").str.len()))
    orf_rows = entries[is_orf].copy()
    coord = orf_rows["note"].str.extract(r"strand ([+-]), nt (\d+)-(\d+)")
    orf_rows["strand"] = coord[0]
    orf_rows["span"] = coord[2].astype(float) - coord[1].astype(float) + 1
    orf_rows["sup"] = supported[is_orf]
    keep_unsup = set()
    for par, g in orf_rows.groupby("parent"):
        if g["sup"].any():
            continue
        cover = g["span"].sum() / max(nt_len.get(par, 0), 1)
        if (len(g) == 1 and cover >= sole_cover) or \
                (len(g) > 1 and g["strand"].nunique() == 1 and cover >= operon_cover):
            keep_unsup |= set(g.index)
    keep_mask = entries.index.isin(keep_unsup)
    drop = is_orf & ~supported & ~keep_mask
    if drop.any():
        entries.loc[drop, "usable"] = False
        entries.loc[drop, "issues"] = [";".join(filter(None, [i, "orf_unsupported"]))
                                       for i in entries.loc[drop, "issues"]]
        keep = set(entries.loc[entries["usable"], "prot_id"])
        hits = hits[hits["query"].isin(keep) & hits["target"].isin(keep)]
    return entries, hits, int(drop.sum())


def _family_edges(hits, p):
    ok = hits[(hits["pident"] >= p["same_family_pident"]) &
              (hits[["qcov", "tcov"]].min(axis=1) >= p["same_family_cov"])]
    return ok


def family_clusters(entries, hits, p):
    """Representative-based family clusters (without chaining).

    The longest protein not yet assigned becomes the representative; only proteins similar to
    THE REPRESENTATIVE ITSELF above the same_family threshold join its cluster. Unlike single
    linkage, two related families therefore do not merge through an intermediate sequence
    (e.g. SHV – LEN – OKP). The order is deterministic (by length, then by identifier).
    """
    u = entries[entries["usable"]].drop_duplicates("prot_id")
    length = dict(zip(u["prot_id"], u["prot_len"]))
    adj = defaultdict(set)
    for q, t in zip(*(_family_edges(hits, p)[c] for c in ("query", "target"))):
        adj[q].add(t)
        adj[t].add(q)
    rep_of = {}
    for pid in sorted(length, key=lambda x: (-length[x], x)):
        if pid in rep_of:
            continue
        rep_of[pid] = pid
        for nb in adj.get(pid, ()):
            rep_of.setdefault(nb, pid)
    members = defaultdict(list)
    for pid, rep in rep_of.items():
        members[rep].append(pid)
    out = {}
    for ms in members.values():
        # the identifier is computed from the membership: unchanged membership, unchanged ID
        fid = "F-" + hashlib.sha256(",".join(sorted(ms)).encode()).hexdigest()[:12]
        for m in ms:
            out[m] = fid
    return out


def family_labels(entries, hits, p, ambiguous_pp=1.0):
    """NCBI family label per protein.

    'ncbi':     the protein is encoded by an AMRFinderPlus entry: its curated family
    'inferred': the family of the closest (fewest differences) AMRFinderPlus protein, if it is
                related at least at same_family level
    Returns {prot_id: (label, source)}
    """
    a = entries[entries["usable"] & (entries["db"] == "amrfinder") & (entries["ncbi_family"] != "")]
    direct = a.groupby("prot_id")["ncbi_family"].agg(lambda s: ";".join(sorted(set(s)))).to_dict()
    out = {pid: (lab, "ncbi") for pid, lab in direct.items()}
    ok = _family_edges(hits, p)
    ok = ok[ok["target"].isin(direct) & ~ok["query"].isin(direct)]
    if len(ok):
        ok = ok.assign(fam=ok["target"].map(direct))
        # best candidate per family, then the two best different families per protein
        per_fam = ok.sort_values(["query", "n_diff", "bits", "target"],
                                 ascending=[True, True, False, True], kind="mergesort") \
                    .drop_duplicates(["query", "fam"])
        for q, g in per_fam.groupby("query", sort=False):
            first = g.iloc[0]
            source = "inferred"
            if len(g) > 1 and abs(float(first["pident"]) - float(g.iloc[1]["pident"])) < ambiguous_pp:
                # the closest and the second closest family hardly differ
                source = "inferred_ambiguous"
            out[q] = (first["fam"], source)
    return out


# ------------------------------------------------------------ NCBI validation
def hierarchy_validation(entries, famrules_file=None):
    """Comparison of the sequence-based family clusters with the curated NCBI families.

    split:  an NCBI family is split over several clusters (our threshold is stricter,
            or the family is genuinely heterogeneous)
    merged: a cluster merges several NCBI families (our threshold is looser than the curated one)
    """
    a = entries[(entries["db"] == "amrfinder") & entries["usable"] &
                (entries["ncbi_family"] != "") & (entries["family_cluster"] != "")]
    if a.empty:
        return None, None
    fam = a.groupby("ncbi_family").agg(n_entries=("key", "size"),
                                       n_clusters=("family_cluster", "nunique"))
    largest = a.groupby(["ncbi_family", "family_cluster"]).size().groupby(level=0).max()
    fam["largest_cluster_share"] = (largest / fam["n_entries"]).round(3)
    fam["status"] = ["split" if n > 1 else "ok" for n in fam["n_clusters"]]
    if famrules_file and Path(famrules_file).exists():
        fr = pd.read_csv(famrules_file, sep="\t", dtype=str, keep_default_na=False)
        fr.columns = [c.lstrip("#") for c in fr.columns]
        if {"node_id", "blastrule_complete_ident"} <= set(fr.columns):
            rule = dict(zip(fr["node_id"], pd.to_numeric(fr["blastrule_complete_ident"], errors="coerce")))
            fam["ncbi_complete_ident"] = [rule.get(f) if rule.get(f) else "" for f in fam.index]
    cl = a.groupby("family_cluster").agg(
        n_entries=("key", "size"), n_ncbi_families=("ncbi_family", "nunique"),
        ncbi_families=("ncbi_family", lambda x: ",".join(sorted(set(x))[:20])))
    cl["status"] = ["merged" if n > 1 else "ok" for n in cl["n_ncbi_families"]]
    return (fam.sort_values(["status", "n_clusters"], ascending=[False, False]),
            cl.sort_values(["status", "n_ncbi_families"], ascending=[True, False]))


# ------------------------------------------------------------ summary
def coverage_summary(entries, mappings):
    rows = []
    dbs = sorted(entries["db"].unique())
    for a in dbs:
        n_total = int((entries["db"] == a).sum())
        n_usable = int(((entries["db"] == a) & entries["usable"]).sum())
        for b in dbs:
            if a == b:
                continue
            sub = mappings[(mappings["db_a"] == a) & (mappings["db_b"] == b)] if len(mappings) else mappings
            best = sub.groupby("entry_a")["relation"].agg(
                lambda s: min(s, key=RELATION_RANK.get)) if len(sub) else pd.Series(dtype=str)
            counts = best.value_counts().to_dict()
            row = dict(db_a=a, db_b=b, entries_total=n_total, entries_usable=n_usable)
            for rel in RELATION_RANK:
                row[rel] = int(counts.get(rel, 0))
            row["no_match"] = n_usable - int(len(best))
            if len(sub):
                sc = sub[sub["best_for_a"]].groupby("entry_a")["tie_scope"].first().value_counts()
            else:
                sc = pd.Series(dtype=int)
            row["ambiguous_within_family"] = int(sc.get("within_family", 0))
            row["ambiguous_cross_family"] = int(sc.get("cross_family", 0))
            rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------ main function
def build(std_files, outdir, params=None, mmseqs_bin="mmseqs", threads=4, keep_tmp=False,
          profile=""):
    p = {**DEFAULT_PARAMS, **(params or {})}
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    entries, versions = load_standardized(std_files)
    entries["issues"] = entries["issues"].fillna("")
    # ORF splitting with curated priority (revision of the selection made at parsing)
    rwork = Path(tempfile.mkdtemp(prefix="argeement_orf_"))
    try:
        entries, orf_rescued = rescue_split_orfs(entries, p, mmseqs_bin, threads, rwork)
    finally:
        shutil.rmtree(rwork, ignore_errors=True)
    if orf_rescued:
        print(f"ORFs of {len(orf_rescued)} split entries reselected with curated priority: "
              f"{', '.join(orf_rescued)}", file=sys.stderr)
    entries["usable"] = entries.apply(
        lambda r: bool(r["protein"]) and not (set(filter(None, r["issues"].split(";"))) & BLOCKING_ISSUES),
        axis=1)
    entries["prot_id"] = [seq_hash(s, "P") if u else "" for s, u in zip(entries["protein"], entries["usable"])]
    entries["nt_id"] = [seq_hash(s, "N") if s else "" for s in entries["nucleotide"]]
    entries["prot_len"] = entries["protein"].str.len()
    term, ambiguous = load_class_map()
    harm = [harmonize(c, g, term, ambiguous) for c, g in zip(entries["classes"], entries["source_group"])]
    entries["class_h"] = [";".join(h) for h, _ in harm]
    entries["unmapped_classes"] = [";".join(u) for _, u in harm]

    uniq = entries[entries["usable"]].drop_duplicates("prot_id")
    print(f"{len(entries)} entries, {int(entries['usable'].sum())} usable, "
          f"{len(uniq)} unique protein sequences", file=sys.stderr)

    work = Path(tempfile.mkdtemp(prefix="argeement_"))
    try:
        fasta = outdir / "proteins.faa"
        with open(fasta, "w") as fh:
            for r in uniq.itertuples(index=False):
                fh.write(f">{r.prot_id}\n{r.protein}\n")
        hits = symmetrize_hits(run_mmseqs(fasta, work, p, mmseqs_bin, threads))
    finally:
        if not keep_tmp:
            shutil.rmtree(work, ignore_errors=True)

    # support of ORF sub-entries: only ORFs related to a curated AMR protein are kept
    entries, hits, n_unsup = orf_support(entries, hits)
    if n_unsup:
        print(f"{n_unsup} ORF sub-entries not related to any curated AMR protein, excluded", file=sys.stderr)

    cands = build_candidates(entries, hits, p)
    fam = family_clusters(entries, hits, p)
    labels = family_labels(entries, hits, p)
    mappings = build_mappings(entries, cands, fam, labels)
    entries["family_label"] = entries["prot_id"].map(lambda x: labels.get(x, ("", ""))[0])
    entries["family_label_source"] = entries["prot_id"].map(lambda x: labels.get(x, ("", ""))[1])
    entries["family_cluster"] = entries["prot_id"].map(fam).fillna("")

    cols = ["key", "db", "entry_id", "name", "family", "accession", "aro", "prot_id",
            "nt_id", "family_cluster", "prot_len", "protein_source", "usable", "issues", "note",
            "ncbi_node", "ncbi_family", "lineage", "classes", "class_h", "unmapped_classes",
            "mechanism", "source_group", "parent", "family_label", "family_label_source", "scope"]
    entries[cols].to_csv(outdir / "entries.tsv", sep="\t", index=False)
    mappings.to_csv(outdir / "mappings.tsv", sep="\t", index=False)
    hits.to_csv(outdir / "protein_hits.tsv", sep="\t", index=False)
    summary = coverage_summary(entries, mappings)
    summary.to_csv(outdir / "summary.tsv", sep="\t", index=False)
    issues = entries[entries["issues"] != ""][["key", "name", "issues", "usable"]]
    issues.to_csv(outdir / "sequence_issues.tsv", sep="\t", index=False)

    # database inventories into the release (for the report)
    for f in std_files:
        inv = Path(f).with_suffix(".inventory")
        if inv.exists():
            shutil.copy(inv, outdir / f"inventory_{Path(f).stem}.json")

    # NCBI hierarchy into the release (for family-level lookup) and validation
    famrules = None
    for f in std_files:
        for suffix, target in ((".ncbi_hierarchy", "ncbi_hierarchy.tsv"),
                               (".ncbi_famrules", "ncbi_famrules.tsv")):
            src = Path(f).with_suffix(suffix)
            if src.exists():
                shutil.copy(src, outdir / target)
                if target == "ncbi_famrules.tsv":
                    famrules = outdir / target
    # name concordance: names of identical proteins across databases
    from .names import name_concordance
    nc = name_concordance(entries, outdir / "ncbi_hierarchy.tsv")
    nc.to_csv(outdir / "name_concordance.tsv", sep="\t", index=False)
    # reverse direction: same name, different sequence
    from .names import same_name_concordance
    sn = same_name_concordance(entries, hits, p, outdir / "ncbi_hierarchy.tsv")
    sn.to_csv(outdir / "same_name_concordance.tsv", sep="\t", index=False)

    vfam, vcl = hierarchy_validation(entries, famrules)
    validation = {}
    if vfam is not None:
        vfam.to_csv(outdir / "validation_ncbi_families.tsv", sep="\t")
        vcl.to_csv(outdir / "validation_clusters.tsv", sep="\t")
        validation = dict(ncbi_families=int(len(vfam)),
                          ncbi_families_split=int((vfam["status"] == "split").sum()),
                          clusters=int(len(vcl)),
                          clusters_merged=int((vcl["status"] == "merged").sum()))
        print(f"NCBI validation: {validation}", file=sys.stderr)

    manifest = dict(
        argeement_version=__version__, profile=profile, mmseqs_version=mmseqs_version(mmseqs_bin),
        created=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        database_versions=versions, params=p,
        counts=dict(entries=int(len(entries)), usable=int(entries["usable"].sum()),
                    unique_proteins=int(len(uniq)), mappings=int(len(mappings)),
                    family_clusters=int(len(set(fam.values()))),
                    # unambiguous names (the keys above are kept for compatibility with earlier releases)
                    entry_rows_incl_orfs=int(len(entries)),
                    source_entries=int((entries["parent"] == "").sum()),
                    protein_level_entries=int(entries["usable"].sum()),
                    unique_proteins_usable=int(entries.loc[entries["usable"], "prot_id"].nunique()),
                    unique_proteins_all=int(len(uniq))),
        ncbi_validation=validation,
        orf_curated_priority=orf_rescued,
        name_concordance={k: int(v) for k, v in nc["category"].value_counts().items()},
        same_name_concordance={k: int(v) for k, v in sn["relation"].value_counts().items()})
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(summary.to_string(index=False), file=sys.stderr)
    return manifest
