"""Synthetic validation: retrieval of variants made with known modifications.

A random sample of the curated proteins of a release (AMRFinderPlus, CARD) is drawn, and known
variants are made from each (amino-acid substitutions, N-terminal extension, C-terminal truncation,
nucleotide CDS on the forward and reverse strand, without start codon, two-gene operon, gene next
to a foreign ORF), which are then placed with the query command. Metrics per variant type:

  truth_found   the original protein is among the best hits (also in a tie)
  truth_unique  the original protein is the only best hit
  relation_ok   the relation type equals the type expected from the modification
  ndiff_ok      the number of differences equals the size of the modification
"""
import json
import random
import tempfile
from pathlib import Path

import pandas as pd

from .build import DEFAULT_PARAMS, classify
from .query import query

AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"], "C": ["TGT", "TGC"], "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"], "F": ["TTT", "TTC"], "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"], "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "M": ["ATG"], "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "V": ["GTT", "GTC", "GTA", "GTG"], "W": ["TGG"],
    "Y": ["TAT", "TAC"]}
STOPS = ["TAA", "TAG", "TGA"]


def _cds(p, r):
    return "ATG" + "".join(r.choice(CODONS.get(a, ["NNN"])) for a in p[1:]) + r.choice(STOPS)


def _rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _subst(p, k, r):
    """k random amino-acid substitutions; also returns their 1-based positions (the random number
    sequence is unchanged, so earlier runs produce the same cases)."""
    s = list(p)
    poss = r.sample(range(1, len(p)), k)
    for pos in poss:
        s[pos] = r.choice([a for a in AA if a != s[pos]])
    return "".join(s), sorted(x + 1 for x in poss)


def make_cases(proteins, r, params):
    """proteins: [(prot_id, db, seq)] -> list of cases (dict) and FASTA records."""
    cases, fasta = [], []

    def add(kind, pid, db, seq, exp_rel, exp_nd, qid=None, sub=None, subst_pos="", edge_dist=""):
        qid = qid or f"{kind}__{pid}"
        cases.append(dict(query=qid, kind=kind, truth=pid, db=db, exp_relation=exp_rel,
                          exp_ndiff=exp_nd, subquery=sub or qid, subst_pos=subst_pos,
                          edge_dist=edge_dist))
        return qid

    for i, (pid, db, p) in enumerate(proteins):
        L = len(p)
        fasta.append((add("exact_protein", pid, db, p, "identical_protein", 0), p))
        k = r.randint(1, 5)
        sp, poss = _subst(p, k, r)
        # positions of the substitutions and their smallest distance from a protein end (aa)
        fasta.append((add(f"subst_{k}aa", pid, db, None, classify(100 * (L - k) / L, 1, 1, params), k,
                          subst_pos=",".join(map(str, poss)),
                          edge_dist=min(min(x - 1, L - x) for x in poss)), sp))
        e = r.randint(3, 10)
        ext = "M" + "".join(r.choice(AA) for _ in range(e - 1)) + p          # new start e aa upstream
        fasta.append((add("n_extension", pid, db, None, classify(100, L / (L + e), 1, params), e), ext))
        t = r.randint(3, 10)
        fasta.append((add("c_truncation", pid, db, None, classify(100, 1, (L - t) / L, params), t),
                      p[:L - t]))
        cds = _cds(p, r)
        fasta.append((add("cds_forward", pid, db, None, "identical_protein", 0), cds))
        fasta.append((add("cds_reverse", pid, db, None, "identical_protein", 0), _rc(cds)))
        fasta.append((add("cds_no_start", pid, db, None, classify(100, 1, (L - 1) / L, params), 1),
                      cds[3:]))
        # two-gene operon: this protein and the next one of the sample
        p2id, p2db, p2 = proteins[(i + 1) % len(proteins)]
        op = f"operon__{pid}"
        add("operon_gene1", pid, db, None, "identical_protein", 0, qid=op, sub=f"{op}#orf1")
        add("operon_gene2", p2id, p2db, None, "identical_protein", 0, qid=op, sub=f"{op}#orf2")
        fasta.append((op, cds + "AGGAGT" + _cds(p2, r)))
        # gene next to a foreign ORF not similar to any resistance gene
        junk = "M" + "".join(r.choice(AA) for _ in range(89))
        jq = f"gene_junk__{pid}"
        add("gene_with_junk_orf", pid, db, None, "identical_protein", 0, qid=jq, sub=f"{jq}#orf1")
        fasta.append((jq, cds + "AGGAGT" + _cds(junk, r)))
    return cases, fasta


def run_synthetic(release, out=None, n=200, seed=1, mmseqs_bin="mmseqs", threads=4,
                  min_len=100, max_len=1000):
    rel = Path(release)
    man = json.loads((rel / "manifest.json").read_text())
    params = {**DEFAULT_PARAMS, **man.get("params", {})}
    e = pd.read_csv(rel / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    u = e[(e["usable"] == "True") & e["db"].isin(["amrfinder", "card"])]
    seqs = {}
    for line_id, line in _read_faa(rel / "proteins.faa"):
        seqs[line_id] = line
    pool = (u.drop_duplicates("prot_id")
             .assign(seq=lambda x: x["prot_id"].map(seqs))
             .dropna(subset=["seq"]))
    pool = pool[pool["seq"].str.len().between(min_len, max_len) & pool["seq"].str.startswith("M")]
    r = random.Random(seed)
    sample = pool.sample(n=min(n, len(pool)), random_state=seed)
    proteins = list(zip(sample["prot_id"], sample["db"], sample["seq"]))
    cases, fasta = make_cases(proteins, r, params)

    out = Path(out) if out else rel / "validation_synthetic"
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "q.fa"
        fa.write_text("".join(f">{q}\n{s}\n" for q, s in fasta))
        res, _ = query(rel, fa, mmseqs_bin=mmseqs_bin, threads=threads, top=1)

    pid_of = dict(zip(e["db"] + "::" + e["entry_id"], e["prot_id"]))
    res = res[res["status"] == "match"].copy()
    res["prot_b"] = (res["db"] + "::" + res["entry_id"]).map(pid_of)
    # pre-grouping by exact query id and by base id without ORF suffix
    res["base"] = res["query"].str.replace(r"#orf\d+$", "", regex=True)
    by_query = {k: g for k, g in res.groupby(["query", "db"])}
    by_base = {k: g for k, g in res.groupby(["base", "db"])}
    empty = res.iloc[0:0]
    rows = []
    for c in cases:
        sq = c["subquery"]
        if "#orf" in sq:          # a specific ORF is the target (operon)
            hit = by_query.get((sq, c["db"]), empty)
        else:                     # translated as a whole, or as an ORF (e.g. reverse strand)
            hit = by_base.get((sq, c["db"]), empty)
        best_pids = set(hit["prot_b"])
        found = c["truth"] in best_pids
        th = hit[hit["prot_b"] == c["truth"]]
        rel_ = th["relation"].iloc[0] if found else (hit["relation"].iloc[0] if len(hit) else "")
        nd = int(th["n_diff"].iloc[0]) if found and pd.notna(th["n_diff"].iloc[0]) else None
        rows.append(dict(c, truth_found=found, truth_unique=found and len(best_pids) == 1,
                         relation=rel_, n_diff=nd,
                         relation_ok=found and rel_ == c["exp_relation"],
                         ndiff_ok=found and nd == c["exp_ndiff"]))
    df = pd.DataFrame(rows)
    df.to_csv(out / "cases.tsv", sep="\t", index=False)
    kind_order = ["exact_protein", "subst_1aa", "subst_2aa", "subst_3aa", "subst_4aa", "subst_5aa",
                  "n_extension", "c_truncation", "cds_forward", "cds_reverse", "cds_no_start",
                  "operon_gene1", "operon_gene2", "gene_with_junk_orf"]
    summ = (df.groupby("kind")
              .agg(n=("query", "size"), truth_found=("truth_found", "mean"),
                   truth_unique=("truth_unique", "mean"), relation_ok=("relation_ok", "mean"),
                   ndiff_ok=("ndiff_ok", "mean"))
              .reindex([k for k in kind_order if k in set(df["kind"])])
              .reset_index())   # raw proportions: rounding only for display
    summ.to_csv(out / "summary.tsv", sep="\t", index=False)
    meta = dict(release=rel.resolve().name, n_proteins=len(proteins), n_cases=len(df), seed=seed,
                overall=dict(truth_found=float(df["truth_found"].mean()),
                             truth_unique=float(df["truth_unique"].mean()),
                             relation_ok=float(df["relation_ok"].mean()),
                             ndiff_ok=float(df["ndiff_ok"].mean())))
    (out / "summary.json").write_text(json.dumps(meta, indent=2))
    return summ, meta


def _read_faa(path):
    name, buf = None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        yield name, "".join(buf)
