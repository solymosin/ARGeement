"""Name concordance: comparison of the names of entries that encode the same protein.

Every cross-database pair of entries with the same protein (once, db_a < db_b) is assigned to the
first matching category:

  exact                  – identical name
  case                   – differs only in letter case
  formatting             – identical without the 'bla' prefix and punctuation (TEM-1 / blaTEM-1)
  species_prefix         – identical without the CARD species prefix (Ecol_emrE / emrE)
  variant_suffix         – identical without the ResFinder/DisinFinder synonym-variant letter
                           (blaTEM-1A / blaTEM-1)
  ncbi_synonym           – listed as a synonym in the NCBI hierarchy
  allele_number_differs  – same gene name, different allele number (genuine conflict)
  different_name         – everything else (genuine conflict)
"""
import re
from pathlib import Path

import pandas as pd

CATEGORIES = ["exact", "case", "formatting", "species_prefix", "variant_suffix",
              "ncbi_synonym", "allele_number_differs", "different_name"]
CONFLICTS = {"allele_number_differs", "different_name"}

CARD_SPECIES_PREFIX = re.compile(r"^[A-Z][a-z]{3}_")       # Ecol_, Saur_, Abau_ ...
VARIANT_SUFFIX = re.compile(r"(?<=\d)[A-Za-z]$")          # blaTEM-1B -> blaTEM-1
TRAILING_NUMBER = re.compile(r"^(.*?)(\d+)$")


def norm(name):
    s = name.lower()
    s = s[3:] if s.startswith("bla") else s
    return re.sub(r"[^a-z0-9]", "", s)


def _strip(name, db):
    if db == "card":
        name = CARD_SPECIES_PREFIX.sub("", name)
    if db in ("resfinder", "disinfinder"):
        name = VARIANT_SUFFIX.sub("", name)
    return name


def classify(a, db_a, b, db_b, syn_a=frozenset(), syn_b=frozenset()):
    """syn_a / syn_b: normalised NCBI synonyms of side 'a' / 'b' (only AMRFinderPlus entries
    have them). They are always compared with the name on the OTHER side."""
    if a == b:
        return "exact"
    if a.lower() == b.lower():
        return "case"
    na, nb = norm(a), norm(b)
    if na == nb:
        return "formatting"
    sa = CARD_SPECIES_PREFIX.sub("", a) if db_a == "card" else a
    sb = CARD_SPECIES_PREFIX.sub("", b) if db_b == "card" else b
    if norm(sa) == norm(sb):
        return "species_prefix"
    va, vb = norm(_strip(a, db_a)), norm(_strip(b, db_b))
    if va == vb:
        return "variant_suffix"
    if ({nb, vb} & syn_a) or ({na, va} & syn_b):
        return "ncbi_synonym"
    ma, mb = TRAILING_NUMBER.match(va), TRAILING_NUMBER.match(vb)
    if ma and mb and ma.group(1) == mb.group(1) and ma.group(1) and ma.group(2) != mb.group(2):
        return "allele_number_differs"
    return "different_name"


def name_concordance(entries, hierarchy_file=None):
    """entries: entries.tsv of the release (DataFrame). Returns a DataFrame with one row per pair."""
    syn_of_node = {}
    if hierarchy_file and Path(hierarchy_file).exists():
        h = pd.read_csv(hierarchy_file, sep="\t", dtype=str, keep_default_na=False)
        h.columns = [c.lstrip("#") for c in h.columns]
        if "synonyms" in h.columns:
            for node, sym, syn in zip(h["node_id"], h.get("symbol", h["node_id"]), h["synonyms"]):
                # only the listed synonyms; the node's own name is not a synonym
                s = {norm(x) for x in syn.split(",") if x.strip()}
                syn_of_node[node] = frozenset(x for x in s if x)

    u = entries[(entries["usable"].astype(str) == "True") & (entries["prot_id"] != "")]
    if "parent" in u.columns:   # ORF sub-entries carry the operon's name: excluded
        u = u[u["parent"].fillna("") == ""]
    node = dict(zip(u["db"] + "::" + u["entry_id"], u.get("ncbi_node", pd.Series("", index=u.index))))
    rows = []
    for pid, g in u.groupby("prot_id"):
        if g["db"].nunique() < 2:
            continue
        recs = list(g[["db", "entry_id", "name"]].itertuples(index=False))
        for i, (da, ea, na) in enumerate(recs):
            for db_, eb, nb in recs[i + 1:]:
                if da == db_:
                    continue
                (da2, ea2, na2), (db2, eb2, nb2) = sorted([(da, ea, na), (db_, eb, nb)])
                def syn_of(d, e_):
                    if d != "amrfinder":
                        return frozenset()
                    return syn_of_node.get(node.get(f"{d}::{e_}", ""), frozenset())
                rows.append(dict(prot_id=pid, db_a=da2, entry_a=ea2, name_a=na2,
                                 db_b=db2, entry_b=eb2, name_b=nb2,
                                 category=classify(na2, da2, nb2, db2,
                                                   syn_of(da2, ea2), syn_of(db2, eb2))))
    cols = ["prot_id", "db_a", "entry_a", "name_a", "db_b", "entry_b", "name_b", "category"]
    return pd.DataFrame(rows, columns=cols).sort_values(["db_a", "db_b", "category", "name_a"])


# ------------------------------------------------------------ reverse direction
SEQ_RELATIONS = ["identical_protein", "near_identical", "length_variant", "same_family",
                 "homolog", "unrelated"]


def name_key(name, db):
    """Comparison key: without formatting, species-prefix and variant-letter differences
    (the same rules as the first five categories of the first direction)."""
    return norm(_strip(name, db))


def same_name_concordance(entries, hits, params, hierarchy_file=None):
    """Sequence relation of pairs of entries with the same name in different databases.

    relation:   identical_protein ... homolog from the protein hits (the same thresholds as in
                the mapping), 'unrelated' if there is no hit above the search threshold
    name_level: 'allele' if the node of the AMRFinderPlus side is an allele in the NCBI hierarchy;
                'gene' if not; 'unknown' if there is no AMRFinderPlus side or node
    """
    from .build import classify   # avoids a circular import

    allele_flag = {}
    if hierarchy_file and Path(hierarchy_file).exists():
        h = pd.read_csv(hierarchy_file, sep="\t", dtype=str, keep_default_na=False)
        h.columns = [c.lstrip("#") for c in h.columns]
        if "allele" in h.columns:
            allele_flag = dict(zip(h["node_id"], h["allele"]))

    u = entries[(entries["usable"].astype(str) == "True") & (entries["prot_id"] != "")].copy()
    if "parent" in u.columns:
        u = u[u["parent"].fillna("") == ""]
    if "ncbi_node" not in u.columns:
        u["ncbi_node"] = ""
    u["key"] = [name_key(n, d) for n, d in zip(u["name"], u["db"])]

    best = {}
    for q, t, pi, qc, tc, nd in zip(hits["query"], hits["target"], hits["pident"], hits["qcov"],
                                     hits["tcov"], hits["n_diff"]):
        k = (q, t)
        if k not in best or nd < best[k][3]:
            best[k] = (pi, qc, tc, int(nd))

    rows = []
    for key, g in u.groupby("key"):
        if not key or g["db"].nunique() < 2:
            continue
        recs = list(g[["db", "entry_id", "name", "prot_id", "ncbi_node"]].itertuples(index=False))
        for i, ra in enumerate(recs):
            for rb in recs[i + 1:]:
                if ra.db == rb.db:
                    continue
                a, b = sorted([ra, rb])
                if a.prot_id == b.prot_id:
                    rel, nd, pi = "identical_protein", 0, 100.0
                elif (a.prot_id, b.prot_id) in best or (b.prot_id, a.prot_id) in best:
                    pi, qc, tc, nd = best.get((a.prot_id, b.prot_id)) or best[(b.prot_id, a.prot_id)]
                    rel = classify(pi, qc, tc, params)
                else:
                    rel, nd, pi = "unrelated", None, None
                node = a.ncbi_node if a.db == "amrfinder" else (b.ncbi_node if b.db == "amrfinder" else "")
                flag = allele_flag.get(node)
                level = "allele" if flag == "1" else ("gene" if flag == "0" else "unknown")
                rows.append(dict(key=key, db_a=a.db, entry_a=a.entry_id, name_a=a.name,
                                 db_b=b.db, entry_b=b.entry_id, name_b=b.name, name_level=level,
                                 relation=rel, n_diff=nd, pident=pi))
    cols = ["key", "db_a", "entry_a", "name_a", "db_b", "entry_b", "name_b", "name_level",
            "relation", "n_diff", "pident"]
    out = pd.DataFrame(rows, columns=cols)
    out["n_diff"] = pd.to_numeric(out["n_diff"]).astype("Int64")
    return out.sort_values(["db_a", "db_b", "key", "entry_a", "entry_b"])
