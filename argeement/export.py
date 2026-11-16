"""Public export: normalised, non-redundant tables of a release, laid out like the export of a
relational database. The tables refer to each other by integer keys only.

databases.tsv       – db_id -> database, version
relation_types.tsv  – relation_id -> relation type and its thresholds
families.tsv        – family_id -> family cluster (F-...)
proteins.tsv        – pid -> protein (P-...), length, family
entries.tsv         – eid -> database entry, with db_id, pid and parent_eid
relations.tsv       – once per protein pair: pid_a, pid_b, relation_id, differences, identity, coverage
schema.sql          – table definitions with primary and foreign keys (SQLite)
LICENSE.md          – licence of the data (CC-BY 4.0) and terms of the source databases
manifest.json       – versions, parameters
README.md           – schema and usage
DATABASES.md        – content report (if available)
SHA256SUMS          – checksums

The integer keys are assigned per release (in the order of the stable identifiers) and are not
stable across releases; the stable identifiers are `prot_id`, `family_cluster` and
(`db`, `entry_id`).

The full, redundant mappings.tsv stays in the release for checking; every row of the export can be
reconstructed from it (and vice versa), which the tests verify.
"""
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

from .terms import licence_md

LINKS = {
    "amrfinder": ("AMRFinderPlus", "`entry_id` = protein accession, as in the "
                  "`refseq_protein_accession` / `genbank_protein_accession` columns of "
                  "`ReferenceGeneCatalog` and in the headers of `AMRProt.fa`."),
    "card": ("CARD", "`entry_id` = ARO accession (`ARO_accession` in `card.json`, `ARO:` field "
             "of the FASTA headers); `accession` = protein accession of "
             "`protein_fasta_protein_homolog_model.fasta`. If one ARO has several sequences, the "
             "second and later ones get a `#2`, `#3` ... suffix in `entry_id`."),
    "resfinder": ("ResFinder", "`entry_id` = FASTA identifier (`<gene>_<variant>_<accession>`), "
                  "as in the `*.fsa` files and the `Gene_accession no.` column of `phenotypes.txt`."),
    "disinfinder": ("DisinFinder", "`entry_id` = FASTA identifier, as in `disinfectants.fsa` and "
                    "`phenotypes.txt` (same format as ResFinder)."),
}

SCHEMA = """-- ARGeement export: table definitions (SQLite). Every table is a TSV file with a header row;
-- empty fields are NULL. Load with the Python snippet in README.md, or with the sqlite3 shell
-- (foreign keys are switched on only after the empty fields have become NULL):
--   .read schema.sql
--   .mode tabs
--   .import --skip 1 databases.tsv databases
--   .import --skip 1 relation_types.tsv relation_types
--   .import --skip 1 families.tsv families
--   .import --skip 1 proteins.tsv proteins
--   .import --skip 1 entries.tsv entries
--   .import --skip 1 relations.tsv relations
--   UPDATE entries SET pid = NULLIF(pid, ''), parent_eid = NULLIF(parent_eid, ''),
--          name = NULLIF(name, ''), accession = NULLIF(accession, ''),
--          nt_id = NULLIF(nt_id, ''), issues = NULLIF(issues, '');
--   UPDATE proteins SET family_id = NULLIF(family_id, ''), family_label = NULLIF(family_label, ''),
--          family_label_source = NULLIF(family_label_source, '');
--   UPDATE databases SET version = NULLIF(version, '');
--   PRAGMA foreign_keys = ON;
--   PRAGMA foreign_key_check;          -- no output: no violation

CREATE TABLE databases (
    db_id       INTEGER PRIMARY KEY,
    db          TEXT NOT NULL UNIQUE,
    database    TEXT NOT NULL,
    version     TEXT
);

CREATE TABLE relation_types (
    relation_id INTEGER PRIMARY KEY,
    relation    TEXT NOT NULL UNIQUE,
    min_pident  REAL,
    coverage    TEXT,
    note        TEXT
);

CREATE TABLE families (
    family_id      INTEGER PRIMARY KEY,
    family_cluster TEXT NOT NULL UNIQUE,
    n_proteins     INTEGER NOT NULL
);

CREATE TABLE proteins (
    pid                 INTEGER PRIMARY KEY,
    prot_id             TEXT NOT NULL UNIQUE,
    length              INTEGER NOT NULL,
    family_id           INTEGER REFERENCES families (family_id),
    family_label        TEXT,
    family_label_source TEXT
);

CREATE TABLE entries (
    eid         INTEGER PRIMARY KEY,
    db_id       INTEGER NOT NULL REFERENCES databases (db_id),
    entry_id    TEXT NOT NULL,
    name        TEXT,
    accession   TEXT,
    pid         INTEGER REFERENCES proteins (pid),
    nt_id       TEXT,
    parent_eid  INTEGER REFERENCES entries (eid),
    issues      TEXT,
    UNIQUE (db_id, entry_id)
);

CREATE TABLE relations (
    pid_a       INTEGER NOT NULL REFERENCES proteins (pid),
    pid_b       INTEGER NOT NULL REFERENCES proteins (pid),
    relation_id INTEGER NOT NULL REFERENCES relation_types (relation_id),
    n_diff      INTEGER NOT NULL,
    pident      REAL,
    cov_a       REAL,
    cov_b       REAL,
    PRIMARY KEY (pid_a, pid_b),
    CHECK (pid_a < pid_b)
);
"""

# table name -> file, in an order that satisfies the foreign keys
TABLES = ["databases", "relation_types", "families", "proteins", "entries", "relations"]

README = """# ARGeement mapping – release `{release}`

Sequence-based mapping between antimicrobial resistance reference databases.
{profile_line}

| Database | Version |
|---|---|
{versions}

## Files

The tables are laid out like the export of a relational database: each has an integer primary key,
and they refer to each other by these keys only (the full table definitions are in `schema.sql`).

| File | One row per | Columns (primary key first; → foreign key) |
|---|---|---|
| `databases.tsv` | source database | `db_id`, `db`, `database`, `version` |
| `relation_types.tsv` | relation type | `relation_id`, `relation`, `min_pident`, `coverage`, `note` |
| `families.tsv` | family cluster | `family_id`, `family_cluster`, `n_proteins` |
| `proteins.tsv` | distinct protein sequence | `pid`, `prot_id`, `length`, `family_id` → families, `family_label`, `family_label_source` |
| `entries.tsv` | database entry | `eid`, `db_id` → databases, `entry_id`, `name`, `accession`, `pid` → proteins, `nt_id`, `parent_eid` → entries, `issues` |
| `relations.tsv` | pair of distinct proteins (each pair once, `pid_a` < `pid_b`) | `pid_a` → proteins, `pid_b` → proteins, `relation_id` → relation_types, `n_diff`, `pident`, `cov_a`, `cov_b` |
| `schema.sql` | – | table definitions (SQLite) with primary and foreign keys |
| `LICENSE.md` | – | licence of the data (CC-BY 4.0) and terms of the source databases |
| `manifest.json` | – | database versions, thresholds, software versions |
| `DATABASES.md` | – | content of the databases (types, drug classes, overlap) |
| `SHA256SUMS` | – | checksums of the files above |

Empty fields are missing values (NULL): an entry without `pid` could not be used at protein level
(the reason is in `issues`), an entry without `parent_eid` is not an ORF sub-entry.

**The integer keys are valid within this release only.** They are assigned in the order of the
stable identifiers (`pid` in the order of `prot_id`, `family_id` in the order of `family_cluster`,
`eid` in the order of `db`, `entry_id`), so they change whenever proteins or entries are added or
removed. To link releases with each other, use `prot_id`, `family_cluster` and (`db`, `entry_id`).

### Counts in `manifest.json`

| Key | Meaning |
|---|---|
| `counts.source_entries` | entries of the source databases (without ORF sub-entries) |
| `counts.entry_rows_incl_orfs` | all rows of the full entry table, including the ORF sub-entries of split entries |
| `counts.protein_level_entries` | usable protein-level entries: entries translated as a whole plus retained ORF sub-entries |
| `counts.unique_proteins_usable` | distinct proteins of usable entries |
| `counts.unique_proteins_all` | distinct proteins including ORFs excluded for lack of support |
| `counts.entries`, `counts.usable`, `counts.unique_proteins` | legacy keys (same as `entry_rows_incl_orfs`, `protein_level_entries`, `unique_proteins_all`) |

All relations are assigned once per protein pair: of the two search directions the alignment with
fewer differences is used (on a tie, the higher bit score), so the relation of a pair does not
depend on which protein was the query.

### Identifiers

`prot_id` is `P-` followed by the first 16 hexadecimal digits of the SHA-256 hash of the
protein sequence (upper case, without a terminal `*`). Entries that encode the same protein
share the same `prot_id` (and `pid`), so identical proteins need no row in `relations.tsv`.
Nucleotide-only databases (ResFinder, DisinFinder) were translated with the bacterial
genetic code (table 11), with an alternative start codon translated as `M`. `nt_id` is
formed the same way from the nucleotide sequence (`N-`), where the source provides one; it
distinguishes synonymous variants that share a protein.

A ResFinder or DisinFinder entry is translated as a whole only if it has no internal stop codon,
its length is a multiple of three, and it has both a start and a stop codon (or lacks one of them
but no reverse-strand ORF covers at least 80% of it, i.e. a truncated reference). All other
entries (multi-gene operons such as `VanHAX`, genes stored on the reverse strand or in another
reading frame) are split into open reading frames (all six frames, ATG/GTG/TTG start, ≥ 60 amino
acids, chosen to cover the sequence best without overlaps, with ATG preferred unless an
alternative start extends the ORF by at least four codons). An ORF is kept if its protein
is identical or related (above the search threshold) to an AMRFinderPlus or CARD protein. If none
of the ORFs of an entry is related to a curated protein, they are kept when the entry is a single
gene (one ORF covering at least 80% of it) or an operon unique to the database (all ORFs on the
same strand, together covering at least 90% of it). Each kept ORF is a separate
entry `<entry_id>#orf<k>` whose `parent_eid` points to the original entry; the original entry and
the excluded ORFs (`orf_unsupported` in `issues`) have no `pid`.

`name` is the gene name as given by the source database (AMRFinderPlus: allele, or gene
family where no allele is defined; CARD: model name; ResFinder, DisinFinder: gene name part of
the FASTA identifier). Further annotation (drug classes, mechanisms, families) can be joined
from the source databases.

Linking to the source databases (the versions listed above):

{links}

### Relations

`relation_types.tsv` lists the relation types with their thresholds; `relation_id` also orders
them from the closest (0, `identical_protein`) to the weakest. `identical_protein` has no row in
`relations.tsv`: identical proteins share a `pid`. Every protein pair was aligned in both directions
(A as query against B, and B against A); where the two alignments differ slightly, the values of the
direction with fewer differences are given (on a tie, the one with the higher alignment score).
`n_diff` is the number of differing positions: mismatches and gap columns of the alignment plus
unaligned ends. `pident` is the percent identity over the alignment, `cov_a` and `cov_b` the
fraction of the proteins `pid_a` and `pid_b` covered by it. Thresholds:

| relation | identity | coverage |
|---|---|---|
| `near_identical` | ≥ {np}% | both ≥ {nc} |
| `length_variant` | ≥ {np}% | one ≥ {nc}, the other ≥ {fc} |
| `same_family` | ≥ {fp}% | both ≥ {fc} |
| `homolog` | ≥ {sp}% | both ≥ {sc} |

All `near_identical` and `length_variant` relations are included. Of the weaker relations
(`same_family`, `homolog`) only the best {topk} per entry and target database are kept, so
these are not exhaustive.

`family_cluster` groups proteins around a representative (the longest protein): a protein
joins a cluster only if it is related to the representative itself at `same_family` level or
closer, so related families are not chained together through intermediate sequences.
`family_label` is the curated NCBI family (from the AMRFinderPlus gene hierarchy):
`family_label_source` is `ncbi` where the protein is an AMRFinderPlus entry, and `inferred`
where it was taken from the closest AMRFinderPlus protein (at least `same_family` level);
`inferred_ambiguous` marks an inferred label where the closest protein of another NCBI family is
almost as close (identity within 1 percentage point), so the label is uncertain. Both are empty
where no AMRFinderPlus protein is close enough. For reporting families, `family_label` is the
recommended choice.

## Example: best matches of database A in database B

```python
import pandas as pd
rd = lambda f: pd.read_csv(f"{{f}}.tsv", sep="\\t", keep_default_na=False, na_values=[""])
db, e, r = rd("databases"), rd("entries"), rd("relations")
e = e.merge(db[["db_id", "db"]], on="db_id")

# relations in both directions (query q -> target t), plus identical proteins (relation_id 0)
both = pd.concat([r.rename(columns={{"pid_a": "q", "pid_b": "t"}}),
                  r.rename(columns={{"pid_a": "t", "pid_b": "q", "cov_a": "cov_b", "cov_b": "cov_a"}})])
p = e["pid"].dropna().unique()
both = pd.concat([both, pd.DataFrame({{"q": p, "t": p, "relation_id": 0, "n_diff": 0}})])

a, b = "resfinder", "card"
m = (e[e.db == a].merge(both, left_on="pid", right_on="q")
                 .merge(e[e.db == b], left_on="t", right_on="pid", suffixes=("_a", "_b")))
best = m[m.groupby("eid_a")["relation_id"].transform("min") == m["relation_id"]]
best = best[best.groupby("eid_a")["n_diff"].transform("min") == best["n_diff"]]
best = best.merge(rd("relation_types")[["relation_id", "relation"]], on="relation_id")
# several rows per entry_id_a = equally good matches (a tie)
print(best[["entry_id_a", "name_a", "entry_id_b", "name_b", "relation", "n_diff"]].to_string(index=False))
```

## Example: loading into SQLite

```python
import sqlite3, pandas as pd
con = sqlite3.connect("argeement.sqlite")
con.executescript(open("schema.sql").read())
con.execute("PRAGMA foreign_keys = ON")
for t in ["databases", "relation_types", "families", "proteins", "entries", "relations"]:
    pd.read_csv(f"{{t}}.tsv", sep="\\t", keep_default_na=False, na_values=[""]) \\
      .to_sql(t, con, if_exists="append", index=False)
con.commit()
```

## Licence

The data are licensed under CC-BY 4.0. The source databases have their own terms; in particular,
commercial use of CARD itself requires a licence from McMaster University. Of CARD, these files
contain only ARO accessions and terms (CC-BY 4.0), no sequences or other CARD material. Details:
`LICENSE.md`.

Generated by ARGeement {argeement}.
"""


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _key(values):
    """Integer keys 1..n in the order of the (stable) identifiers."""
    return {v: i for i, v in enumerate(sorted(values), start=1)}


def relation_types(p):
    """The relation types with their thresholds; relation_id = rank (0 = closest)."""
    from .build import RELATION_RANK
    cov = {"identical_protein": "identical sequence",
           "near_identical": f"both >= {p['near_identical_cov']}",
           "length_variant": f"one >= {p['near_identical_cov']}, the other >= {p['same_family_cov']}",
           "same_family": f"both >= {p['same_family_cov']}",
           "homolog": f"both >= {p['search_min_cov']}"}
    pid = {"identical_protein": 100, "near_identical": p["near_identical_pident"],
           "length_variant": p["near_identical_pident"], "same_family": p["same_family_pident"],
           "homolog": round(p["search_min_seq_id"] * 100)}
    note = {"identical_protein": "no row in relations.tsv: identical proteins share a pid",
            "near_identical": "probably an allele-level difference",
            "length_variant": "same gene, different reference length",
            "same_family": f"same gene family; only the best {p['top_k_weak']} per entry and "
                           "target database",
            "homolog": f"distant homolog; only the best {p['top_k_weak']} per entry and "
                       "target database"}
    return pd.DataFrame([dict(relation_id=k, relation=r, min_pident=pid[r], coverage=cov[r],
                              note=note[r])
                         for r, k in sorted(RELATION_RANK.items(), key=lambda kv: kv[1])])


def export(release, out=None):
    rel = Path(release)
    out = Path(out) if out else rel / "export"
    out.mkdir(parents=True, exist_ok=True)
    man = json.loads((rel / "manifest.json").read_text())
    p = man["params"]
    e = pd.read_csv(rel / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    m = pd.read_csv(rel / "mappings.tsv", sep="\t", dtype=str, keep_default_na=False)
    if "parent" not in e.columns:
        e["parent"] = ""
    for c in ("family_label", "family_label_source"):
        if c not in e.columns:
            e[c] = ""
    usable = e["usable"] == "True"
    tables = {}

    # ------------------------------------------------------------ databases
    dbs = [d for d in LINKS if d in set(e["db"])] + sorted(set(e["db"]) - set(LINKS))
    db_id = {d: i for i, d in enumerate(dbs, start=1)}
    tables["databases"] = pd.DataFrame(dict(
        db_id=list(db_id.values()), db=dbs,
        database=[LINKS[d][0] if d in LINKS else d for d in dbs],
        version=[man["database_versions"].get(d, "") for d in dbs]))

    # ------------------------------------------------------------ relation types
    rt = relation_types(p)
    tables["relation_types"] = rt
    rel_id = dict(zip(rt["relation"], rt["relation_id"]))

    # ------------------------------------------------------------ proteins, families
    u = e[usable]
    prot = (u.groupby("prot_id")
             .agg(length=("prot_len", "first"), family_cluster=("family_cluster", "first"),
                  family_label=("family_label", "first"),
                  family_label_source=("family_label_source", "first"))
             .reset_index())
    pid = _key(prot["prot_id"])
    fam_id = _key(set(prot["family_cluster"]) - {""})
    prot.insert(0, "pid", prot["prot_id"].map(pid))
    prot.insert(3, "family_id", prot["family_cluster"].map(fam_id).astype("Int64"))
    fams = prot[prot["family_cluster"] != ""].groupby("family_cluster").size()
    tables["families"] = pd.DataFrame(dict(family_id=[fam_id[f] for f in fams.index],
                                           family_cluster=fams.index, n_proteins=fams.values)) \
        .sort_values("family_id")
    tables["proteins"] = prot.drop(columns="family_cluster").sort_values("pid")

    # ------------------------------------------------------------ entries
    ent = e.copy()
    ent["_o"] = ent["db"].map(db_id)
    ent = ent.sort_values(["_o", "entry_id"]).reset_index(drop=True)
    ent.insert(0, "eid", range(1, len(ent) + 1))
    eid = dict(zip(zip(ent["db"], ent["entry_id"]), ent["eid"]))
    ent["db_id"] = ent["_o"]
    ent["pid"] = ent["prot_id"].where(ent["usable"] == "True", "").map(pid).astype("Int64")
    ent["parent_eid"] = pd.Series([eid.get((d, x)) if x else None
                                   for d, x in zip(ent["db"], ent["parent"])],
                                  dtype="Int64")
    missing = ent[(ent["parent"] != "") & ent["parent_eid"].isna()]
    if len(missing):
        raise ValueError(f"parent entry not found: {missing['entry_id'].iloc[0]}")
    tables["entries"] = ent[["eid", "db_id", "entry_id", "name", "accession", "pid", "nt_id",
                             "parent_eid", "issues"]]

    # ------------------------------------------------------------ relations
    r = m[m["relation"] != "identical_protein"][
        ["prot_a", "prot_b", "relation", "n_diff", "pident", "cov_a", "cov_b", "bits"]].copy()
    swap = r["prot_a"] > r["prot_b"]
    r.loc[swap, ["prot_a", "prot_b", "cov_a", "cov_b"]] = \
        r.loc[swap, ["prot_b", "prot_a", "cov_b", "cov_a"]].values
    # every protein pair was aligned in both directions (A->B and B->A); the values of the direction
    # with fewer differences are kept (on a tie, the one with the higher bit score)
    r["_n"] = pd.to_numeric(r["n_diff"], errors="coerce")
    r["_b"] = pd.to_numeric(r["bits"], errors="coerce")
    r = (r.sort_values(["prot_a", "prot_b", "_n", "_b"], ascending=[True, True, True, False])
          .drop_duplicates(["prot_a", "prot_b"]))
    # pid follows the order of prot_id, so prot_a < prot_b implies pid_a < pid_b
    tables["relations"] = pd.DataFrame(dict(
        pid_a=r["prot_a"].map(pid).values, pid_b=r["prot_b"].map(pid).values,
        relation_id=r["relation"].map(rel_id).values, n_diff=r["n_diff"].values,
        pident=r["pident"].values, cov_a=r["cov_a"].values, cov_b=r["cov_b"].values))
    if tables["relations"][["pid_a", "pid_b", "relation_id"]].isna().any().any():
        raise ValueError("relations.tsv: a protein or relation type has no key")

    for t in TABLES:
        tables[t].to_csv(out / f"{t}.tsv", sep="\t", index=False)
    (out / "schema.sql").write_text(SCHEMA, encoding="utf-8")
    (out / "LICENSE.md").write_text(licence_md(rel.resolve().name), encoding="utf-8")

    # ------------------------------------------------------------ manifest, documentation
    shutil.copy(rel / "manifest.json", out / "manifest.json")
    if (rel / "DATABASES.md").exists():
        shutil.copy(rel / "DATABASES.md", out / "DATABASES.md")
    known = [d for d in dbs if d in LINKS]
    (out / "README.md").write_text(README.format(
        release=rel.resolve().name,
        profile_line=f"\nProfile: `{man['profile']}`, created {man['created']}.\n"
        if man.get("profile") else f"\nCreated {man['created']}.\n",
        versions="\n".join(f"| {LINKS[d][0]} | `{man['database_versions'].get(d, '?')}` |" for d in known),
        links="\n".join(f"- **{LINKS[d][0]}:** {LINKS[d][1]}" for d in known),
        np=p["near_identical_pident"], nc=p["near_identical_cov"],
        fp=p["same_family_pident"], fc=p["same_family_cov"],
        sp=round(p["search_min_seq_id"] * 100), sc=p["search_min_cov"], topk=p["top_k_weak"],
        argeement=man.get("argeement_version", "?")), encoding="utf-8")

    files = sorted(f for f in out.iterdir() if f.is_file() and f.name != "SHA256SUMS")
    (out / "SHA256SUMS").write_text("".join(f"{_sha256(f)}  {f.name}\n" for f in files))

    sizes = {f.name: f.stat().st_size for f in files}
    full = (rel / "mappings.tsv").stat().st_size + (rel / "entries.tsv").stat().st_size
    return dict(out=str(out), rows=dict(entries=len(tables["entries"]),
                                        proteins=len(tables["proteins"]),
                                        relations=len(tables["relations"]), mappings_full=len(m)),
                size_export=sum(sizes.values()), size_full=full)
