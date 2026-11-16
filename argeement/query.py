"""Placement of external sequences (e.g. a detected gene identified only at family level) in the
databases of a release: the closest entry or entries per database, with the number of differences."""
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from .build import DEFAULT_PARAMS, RELATION_RANK, classify
from .sequtil import choose_translation, clean_protein, n_diff_from_alignment, seq_hash

NT = set("ACGTUN")


def _proteins_from_input(fasta, min_orf_aa=60):
    """Returns [(query id, protein, note)]. Nucleotide sequences are translated and, if needed,
    split into ORFs (as for ResFinder entries)."""
    out = []
    for r in SeqIO.parse(fasta, "fasta"):
        seq = str(r.seq).upper().replace("*", "")
        if set(seq) <= NT:
            # same decision as at parsing: strict whole-sequence translation, otherwise ORFs
            aa, issues, orfs = choose_translation(seq, min_orf_aa)
            if orfs:
                for k, o in enumerate(orfs, 1):
                    out.append((f"{r.id}#orf{k}", o["protein"],
                                f"ORF {k}/{len(orfs)}, strand {o['strand']}, nt {o['start'] + 1}-{o['end']}"))
            elif {"internal_stop", "len_not_multiple_of_3", "empty"} & set(issues):
                out.append((r.id, "", "no ORF found"))
            else:
                out.append((r.id, aa, "translated"))
        else:
            out.append((r.id, clean_protein(seq), "protein"))
    return out


def query(release, fasta, mmseqs_bin="mmseqs", threads=4, top=1):
    rel = Path(release)
    man = json.loads((rel / "manifest.json").read_text())
    p = {**DEFAULT_PARAMS, **man.get("params", {})}
    e = pd.read_csv(rel / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    e = e[e["usable"] == "True"]
    if "family_label" not in e.columns:
        e["family_label"] = ""
    dbs = sorted(set(e["db"]))

    qs = _proteins_from_input(fasta)
    work = Path(tempfile.mkdtemp(prefix="argeement_q_"))
    try:
        qfa = work / "q.faa"
        with open(qfa, "w") as fh:
            for i, (qid, aa, _) in enumerate(qs):
                if aa:
                    fh.write(f">q{i}\n{aa}\n")
        cols = "query,target,pident,alnlen,qlen,tlen,qstart,qend,tstart,tend,evalue,bits"
        hits = pd.DataFrame(columns=cols.split(","))
        if qfa.stat().st_size:
            subprocess.run([mmseqs_bin, "easy-search", str(qfa), str(rel / "proteins.faa"),
                            str(work / "h.m8"), str(work / "tmp"),
                            "--min-seq-id", str(p["search_min_seq_id"]), "-c", str(p["search_min_cov"]),
                            "--cov-mode", "0", "-s", str(p["search_sensitivity"]),
                            "--max-seqs", str(p["search_max_seqs"]), "--threads", str(threads),
                            "--format-output", cols, "-v", "1"], check=True)
            hits = pd.read_csv(work / "h.m8", sep="\t", names=cols.split(","))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    if len(hits):
        hits["qcov"] = (hits["qend"] - hits["qstart"] + 1) / hits["qlen"]
        hits["tcov"] = (hits["tend"] - hits["tstart"] + 1) / hits["tlen"]
        nident = (hits["alnlen"] * hits["pident"] / 100).round()
        hits["n_diff"] = n_diff_from_alignment(hits["alnlen"], nident, hits["qstart"], hits["qend"],
                                               hits["qlen"], hits["tstart"], hits["tend"],
                                               hits["tlen"]).astype(int)

    # pre-indexing: run time linear in the number of queries
    from collections import defaultdict
    entries_of = defaultdict(list)            # (prot_id, db) -> entries
    pid_dbs = defaultdict(set)                # prot_id -> databases
    for r_ in e.itertuples(index=False):
        entries_of[(r_.prot_id, r_.db)].append((r_.entry_id, r_.name, r_.family_label))
        pid_dbs[r_.prot_id].add(r_.db)
    hits_of = defaultdict(list)
    if len(hits):
        for h in hits.itertuples(index=False):
            hits_of[h.query].append(h)

    rows = []
    for i, (qid, aa, note) in enumerate(qs):
        if not aa:
            rows.append(dict(query=qid, note=note, status="not_translatable"))
            continue
        cands = []
        pid = seq_hash(aa, "P")
        if pid in pid_dbs:
            cands.append(dict(target=pid, relation="identical_protein", n_diff=0, pident=100.0,
                              cov_query=1.0, cov_target=1.0))
        for h in hits_of.get(f"q{i}", ()):
            if h.target == pid:
                continue
            cands.append(dict(target=h.target, relation=classify(h.pident, h.qcov, h.tcov, p),
                              n_diff=int(h.n_diff), pident=round(float(h.pident), 2),
                              cov_query=round(h.qcov, 3), cov_target=round(h.tcov, 3)))
        # candidates per database at protein level; entries only for those displayed
        per_db = defaultdict(list)
        for c in cands:
            for d in pid_dbs.get(c["target"], ()):
                per_db[d].append(c)
        found = False
        for d in dbs:
            cd = per_db.get(d)
            if not cd:
                continue
            found = True
            cd.sort(key=lambda c: (RELATION_RANK[c["relation"]], c["n_diff"]))
            best = (RELATION_RANK[cd[0]["relation"]], cd[0]["n_diff"])
            tied = [c for c in cd if (RELATION_RANK[c["relation"]], c["n_diff"]) == best]
            # tie: same relation type and number of differences, DIFFERENT proteins
            n_eq = len({c["target"] for c in tied})
            shown = tied if top == 1 else cd[:top]
            for c in shown:
                for entry_id, name, fam in entries_of[(c["target"], d)]:
                    rows.append(dict(query=qid, note=note, status="match", db=d, entry_id=entry_id,
                                     name=name, relation=c["relation"], n_diff=c["n_diff"],
                                     pident=c["pident"], cov_query=c["cov_query"],
                                     cov_target=c["cov_target"], family_label=fam,
                                     n_equally_good=n_eq))
        if not found:
            rows.append(dict(query=qid, note=note, status="no_match"))
    cols = ["query", "note", "status", "db", "entry_id", "name", "relation", "n_diff", "pident",
            "cov_query", "cov_target", "family_label", "n_equally_good"]
    out = pd.DataFrame(rows, columns=cols)
    for c in ("n_diff", "n_equally_good"):
        out[c] = pd.to_numeric(out[c]).astype("Int64")
    return out, man
