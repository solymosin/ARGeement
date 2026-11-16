"""Database-specific parsers.

Each returns the same kind of records (dict) with the following fields:
  db, entry_id, name, family, accession, aro, protein, nucleotide,
  protein_source ('native' | 'translated'), issues (list), note

The assumptions about header formats are documented in the docstrings of the functions.
These are the first thing to check with a new database version: if a format changes,
the parser raises an error or a warning instead of silently skipping entries.
"""
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

from .sequtil import choose_translation, clean_nucleotide, clean_protein, find_orfs, translate_cds

# RefSeq (WP_000027057.1) or GenBank (AAA25547.1) protein accession
PROT_ACC = re.compile(r"^(?:[A-Z]{2}_\d+|[A-Z]{3}\d{5,7})\.\d+$")


def warn(msg):
    print(f"[warning] {msg}", file=sys.stderr)


def _record(db, entry_id, name, **kw):
    rec = dict(db=db, entry_id=entry_id, name=name, family="", accession="",
               aro="", protein="", nucleotide="", protein_source="native",
               issues=[], note="", ncbi_node="", ncbi_family="", lineage="",
               classes="", mechanism="", source_group="", parent="", scope="")
    rec.update(kw)
    return rec


def _mut_sources_text(sources):
    prot = [x for x in sources if x.startswith("AMRProt-mutation")]
    dna = [x for x in sources if x.startswith("AMR_DNA-")]
    parts = prot + ([f"{len(dna)} AMR_DNA-* files (rRNA)"] if dna else [])
    return "genes and mutations from " + (" + ".join(parts) if parts else "no mutation file found")


def genus_of(taxon):
    """'Escherichia coli', 'Acinetobacter_baumannii', 'escherichia_coli' -> genus."""
    t = re.split(r"[\s_]+", str(taxon).strip())[0] if str(taxon).strip() else ""
    return t[:1].upper() + t[1:].lower() if t else "(unspecified)"


def _add_mut(table, genus, gene, n_mut):
    g = table.setdefault(genus, {"genes": set(), "mutations": 0})
    g["genes"].add(gene)
    g["mutations"] += n_mut


def _finish_mut(table):
    return {k: {"genes": len(v["genes"]), "mutations": v["mutations"]} for k, v in sorted(table.items())}


# ---------------------------------------------------------------- CARD
def parse_card(raw_dir):
    """CARD 'broadstreet' csomag (https://card.mcmaster.ca/latest/data).

    Files used:
      protein_fasta_protein_homolog_model.fasta
        header: gb|<protein acc>|ARO:<number>|<name> [<species>]
      nucleotide_fasta_protein_homolog_model.fasta   (optional)
        header: gb|<nucl. acc>|<strand>|<coordinates>|ARO:<number>|<name> [<species>]
      card.json  (version: '_version'; gene family: ARO_category)

    Only protein homolog models are used; variant/overexpression models are
    mutation-based, and sequence identity is not informative for them.
    """
    raw = Path(raw_dir)
    version, families = "unknown", {}
    drug_classes, mechanisms = defaultdict(set), defaultdict(set)
    model_types = Counter()
    card_mut = {}
    cj = raw / "card.json"
    if cj.exists():
        with open(cj) as fh:
            card = json.load(fh)
        version = card.get("_version", "unknown")
        for key, model in card.items():
            if key.startswith("_") or not isinstance(model, dict):
                continue
            aro = model.get("ARO_accession")
            mtype = model.get("model_type") or "unknown"
            model_types[mtype] += 1
            if "variant" in mtype:
                seqs = ((model.get("model_sequences") or {}).get("sequence") or {}).values()
                taxa = {(sq.get("NCBI_taxonomy") or {}).get("NCBI_taxonomy_name", "") for sq in seqs}
                snps = ((model.get("model_param") or {}).get("snp") or {}).get("param_value") or {}
                for tx in (taxa or {""}):
                    _add_mut(card_mut, genus_of(tx), model.get("model_name", aro), len(snps))
            for cat in (model.get("ARO_category") or {}).values():
                kind, val = cat.get("category_aro_class_name"), cat.get("category_aro_name", "")
                if kind == "AMR Gene Family":
                    families[f"ARO:{aro}"] = val
                elif kind == "Drug Class":
                    drug_classes[f"ARO:{aro}"].add(val)
                elif kind == "Resistance Mechanism":
                    mechanisms[f"ARO:{aro}"].add(val)
    else:
        warn("CARD: card.json is missing, continuing without version and gene families")

    nt_by_aro = defaultdict(list)
    nt_file = raw / "nucleotide_fasta_protein_homolog_model.fasta"
    if nt_file.exists():
        for r in SeqIO.parse(nt_file, "fasta"):
            parts = r.description.split("|", 5)
            if len(parts) == 6:
                nt_by_aro[parts[4]].append(str(r.seq))

    records, seen = [], defaultdict(int)
    for r in SeqIO.parse(raw / "protein_fasta_protein_homolog_model.fasta", "fasta"):
        parts = r.description.split("|", 3)
        if len(parts) != 4 or not parts[2].startswith("ARO:"):
            raise ValueError(f"CARD: unexpected header format: {r.description[:80]}")
        _, acc, aro, rest = parts
        name = re.sub(r"\s*\[[^\]]*\]\s*$", "", rest).strip()
        seen[aro] += 1
        entry_id = aro if seen[aro] == 1 else f"{aro}#{seen[aro]}"
        nts = nt_by_aro.get(aro, [])
        records.append(_record(
            "card", entry_id, name, accession=acc, aro=aro,
            family=families.get(aro, ""), protein=clean_protein(str(r.seq)),
            nucleotide=clean_nucleotide(nts[0]) if len(nts) == 1 else "",
            classes=";".join(sorted(drug_classes.get(aro, ()))),
            mechanism=";".join(sorted(mechanisms.get(aro, ())))))
    inventory = dict(source_categories={f"model_type: {k}": v for k, v in model_types.most_common()},
                     included_keys=["model_type: protein homolog model"],
                     mutations_by_genus=_finish_mut(card_mut),
                     mutation_unit="variant models (genes) and listed SNPs",
                     included=len(records))
    return records, version, inventory


# ---------------------------------------------------------- AMRFinderPlus
HIERARCHY_FILES = ("ReferenceGeneHierarchy.txt", "ReferenceGeneHierarchy.tsv")
FAMRULE_FILES = ("fam.tsv", "fam.tab", "fam")
CATALOG_FILES = ("ReferenceGeneCatalog.txt", "ReferenceGeneCatalog.tsv")


def first_existing(raw, names):
    return next((Path(raw) / n for n in names if (Path(raw) / n).exists()), None)


def load_hierarchy(path):
    """ReferenceGeneHierarchy: node_id -> dict(parent, symbol, allele, name)."""
    nodes = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rd.fieldnames = [f.lstrip("#").strip() for f in rd.fieldnames]
        if missing := {"node_id", "parent_node_id"} - set(rd.fieldnames):
            raise ValueError(f"Hierarchy: missing column(s): {missing}")
        for row in rd:
            nodes[row["node_id"]] = dict(parent=row.get("parent_node_id", ""),
                                         symbol=row.get("symbol", ""),
                                         allele=row.get("allele", "0"),
                                         name=row.get("name", ""))
    return nodes


def lineage_of(node, nodes, stop=("ALL", "")):
    """Ancestors of a node towards the root, starting with the node itself (cycle-safe)."""
    path, seen = [], set()
    while node and node not in stop and node in nodes and node not in seen:
        seen.add(node)
        path.append(node)
        node = nodes[node]["parent"]
    return path


def ncbi_family_of(node, nodes):
    """For an allele leaf the parent is the family (blaTEM-1 -> blaTEM), otherwise the node itself."""
    if node in nodes and nodes[node]["allele"] == "1":
        return nodes[node]["parent"]
    return node


def parse_amrfinder(raw_dir, types=("AMR",),
                    exclude_subtypes=("POINT", "POINT_DISRUPT", "AMR-SUSCEPTIBLE"),
                    scopes=("core", "plus"), include_subtypes=("BIOCIDE",)):
    """AMRFinderPlus database (NCBI FTP .../AMRFinderPlus/database/latest/).

    Files used:
      ReferenceGeneCatalog.txt/.tsv – source of names and types (TSV with header)
      ReferenceGeneHierarchy.txt/.tsv – NCBI gene hierarchy (optional): gives the node,
        NCBI family and list of ancestors of the entries
      AMRProt.fa (or AMRProt)  – protein sequences
      AMR_CDS.fa (or AMR_CDS)  – nucleotide CDS (optional)
      version.txt

    The exact field order of the FASTA headers is NOT relied upon: the protein accession is
    searched for in every header, and name and type are taken from the catalogue.
    Sequences not in the catalogue are skipped and reported.
    """
    raw = Path(raw_dir)
    version = (" ".join((raw / "version.txt").read_text().split())
               if (raw / "version.txt").exists() else "unknown")

    cat_file = first_existing(raw, CATALOG_FILES)
    if cat_file is None:
        raise FileNotFoundError(f"AMRFinderPlus: no catalogue file in {raw}")
    hfile = first_existing(raw, HIERARCHY_FILES)
    nodes = load_hierarchy(hfile) if hfile else {}
    if not nodes:
        warn("AMRFinderPlus: no hierarchy file, continuing without NCBI families")
    with open(cat_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        reader.fieldnames = [f.lstrip("#").strip() for f in reader.fieldnames]
        required = {"type", "subtype", "refseq_protein_accession"}
        if missing := required - set(reader.fieldnames):
            raise ValueError(f"AMRFinderPlus: missing catalogue column(s): {missing}")
        by_acc = defaultdict(list)   # an accession can have several catalogue rows
        cat_counts = Counter()
        for row in reader:
            cat_counts[f"{row.get('type', '')}/{row.get('subtype', '')}"] += 1
            for col in ("refseq_protein_accession", "genbank_protein_accession"):
                acc = (row.get(col) or "").strip()
                if acc:
                    by_acc[acc].append(row)

    def category_ok(typ, subtype):
        # the configured types without the excluded subtypes, plus separately included subtypes
        # (e.g. STRESS/BIOCIDE for the comparison with DisinFinder)
        return (typ in types and subtype not in exclude_subtypes) or subtype in include_subtypes

    def row_ok(row):
        return (category_ok(row["type"], row["subtype"])
                and (not row.get("scope") or row["scope"] in scopes))

    def first_file(*names):
        return next((raw / n for n in names if (raw / n).exists()), None)

    def find_acc(description):
        for tok in re.split(r"[|\s]", description):
            if PROT_ACC.match(tok) and tok in by_acc:
                return tok
        return None

    nt_by_acc = {}
    cds = first_file("AMR_CDS.fa", "AMR_CDS")
    if cds:
        for r in SeqIO.parse(cds, "fasta"):
            acc = find_acc(r.description)
            if acc:
                nt_by_acc[acc] = str(r.seq)

    prot = first_file("AMRProt.fa", "AMRProt")
    if prot is None:
        raise FileNotFoundError("AMRFinderPlus: AMRProt.fa not found")
    records, not_in_catalog, filtered, seen = [], 0, 0, set()
    for r in SeqIO.parse(prot, "fasta"):
        # mutation references (e.g. WP_...:90:cirA_S90YfsTer15ins2|...|mutation|...):
        # mostly mutant proteins that would only add noise to the sequence comparison
        if "|mutation|" in r.description:
            filtered += 1
            continue
        acc = find_acc(r.description)
        if acc is None:
            not_in_catalog += 1
            continue
        # an entry is included if ANY of its catalogue rows passes the filters
        # (e.g. an acquired gene may have a mutation row besides its own row)
        ok_rows = [row for row in by_acc[acc] if row_ok(row)]
        if not ok_rows:
            filtered += 1
            continue
        row = ok_rows[0]
        node = (row.get("hierarchy_node") or "").strip()
        if acc in seen:
            continue
        seen.add(acc)
        name = (row.get("allele") or "").strip() or (row.get("gene_family") or "").strip() or acc
        records.append(_record(
            "amrfinder", acc, name, accession=acc,
            family=(row.get("gene_family") or "").strip(),
            protein=clean_protein(str(r.seq)),
            nucleotide=clean_nucleotide(nt_by_acc.get(acc, "")),
            note=(row.get("product_name") or "").strip(),
            ncbi_node=node,
            ncbi_family=ncbi_family_of(node, nodes) if node in nodes else "",
            lineage=">".join(lineage_of(node, nodes)),
            classes=";".join(c for c in (row.get("class") or "").split("/") if c),
            # scope: core if any catalogue row of the accession is core
            scope=("core" if any((x.get("scope") or "").strip() == "core" for x in by_acc[acc])
                   else (row.get("scope") or "").strip())))
    if not_in_catalog:
        warn(f"AMRFinderPlus: {not_in_catalog} proteins not found in the catalogue, skipped")
    print(f"AMRFinderPlus: {len(records)} entries, {filtered} filtered out by type/scope/mutation",
          file=sys.stderr)
    amr_mut, mut_sources = {}, []
    mf = first_existing(raw, ("AMRProt-mutation.tsv", "AMRProt-mutation.tab"))
    if mf:
        mut_sources.append(mf.name)
        with open(mf, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            rd.fieldnames = [f.lstrip("#").strip() for f in rd.fieldnames]
            for row in rd:
                sym = row.get("standard_mutation_symbol", "")
                _add_mut(amr_mut, genus_of(row.get("taxgroup", "")), sym.split("_")[0], 1)
    for df in sorted(raw.glob("AMR_DNA-*.tsv")) + sorted(raw.glob("AMR_DNA-*.tab")):
        mut_sources.append(df.name)
        genus = genus_of(df.stem.split("-", 1)[1])
        with open(df, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            rd.fieldnames = [f.lstrip("#").strip() for f in rd.fieldnames]
            for row in rd:
                sym = row.get("standard_mutation_symbol") or row.get("mutation_symbol") or "?"
                _add_mut(amr_mut, genus, sym.split("_")[0], 1)

    scope_counts = Counter()
    for rec in records:
        for row in by_acc[rec["entry_id"]]:
            if row_ok(row):
                scope_counts[row.get("scope") or "?"] += 1
                break
    inventory = dict(source_categories={f"catalog type/subtype: {k}": v
                                        for k, v in cat_counts.most_common()},
                     included_keys=[f"catalog type/subtype: {k}" for k in cat_counts
                                    if category_ok(*k.split("/", 1))],
                     included=len(records), included_scope=dict(scope_counts),
                     mutations_by_genus=_finish_mut(amr_mut),
                     mutation_unit=_mut_sources_text(mut_sources))
    return records, version, inventory


# --------------------------------------------------------------- ResFinder
RESFINDER_ID = re.compile(r"^(?P<name>.+?)_(?P<variant>\d+)_(?P<acc>.+)$")


# ORFs shorter than this (amino acids) are not considered genes when splitting operons
MIN_ORF_AA = 60


def parse_resfinder(raw_dir, db="resfinder"):
    """ResFinder database (git: genomicepidemiology/resfinder_db).

    All *.fsa files are read; identifiers have the format
    <gene name>_<variant number>_<accession>, e.g. blaTEM-1B_1_AY458016.
    Only nucleotide sequences are available, so they are translated (table 11).
    The version comes from VERSION.argeement written by the fetch step (git commit).
    """
    raw = Path(raw_dir)
    vfile = raw / "VERSION.argeement"
    version = " ".join(vfile.read_text().split()) if vfile.exists() else "unknown"

    seqs, src_file = {}, {}
    for fsa in sorted(raw.glob("*.fsa")):
        for r in SeqIO.parse(fsa, "fasta"):
            s = clean_nucleotide(str(r.seq))
            if r.id in seqs and seqs[r.id] != s:
                warn(f"{db}: {r.id} occurs in several files with different sequences")
            seqs.setdefault(r.id, s)
            if fsa.stem != "all":
                src_file.setdefault(r.id, fsa.stem)

    # classes: from phenotypes.txt if available, otherwise from the source file name
    pheno_class = {}
    pfile = raw / "phenotypes.txt"
    if pfile.exists():
        with open(pfile, newline="", encoding="utf-8", errors="replace") as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = [h.strip().lower() for h in next(rd)]
            ci = next((i for i, h in enumerate(header) if "accession" in h), None)
            ki = next((i for i, h in enumerate(header) if h == "class"), None)
            if ci is None or ki is None:
                warn(f"{db}: columns of phenotypes.txt not recognised, using the file names")
            else:
                for row in rd:
                    if len(row) > max(ci, ki):
                        pheno_class[row[ci].strip()] = [c.strip() for c in re.split(r"[,;]", row[ki]) if c.strip()]
    if not seqs:
        raise FileNotFoundError(f"{db}: no *.fsa file in {raw}")

    records = []
    for rid, nt in seqs.items():
        m = RESFINDER_ID.match(rid)
        name, acc = (m["name"], m["acc"]) if m else (rid, "")
        if not m:
            warn(f"{db}: non-standard identifier: {rid}")
        aa, issues = translate_cds(nt)
        cls = pheno_class.get(rid) or ([src_file[rid]] if rid in src_file else [])
        # if it cannot be translated as a single coding sequence (operon, reverse strand, frame shift),
        # the ORFs found in the six frames become sub-entries: <identifier>#orf<k>
        aa, issues, orfs = choose_translation(nt, MIN_ORF_AA)
        if orfs:
            issues = issues + ["split_into_orfs"]   # the split parent does not count as a protein
        for k, o in enumerate(orfs, 1):
            records.append(_record(
                db, f"{rid}#orf{k}", name, accession=acc, protein=o["protein"],
                protein_source="orf", classes=";".join(cls), source_group=src_file.get(rid, ""),
                parent=rid,
                note=f"ORF {k}/{len(orfs)}, strand {o['strand']}, nt {o['start'] + 1}-{o['end']}"))
        records.append(_record(db, rid, name, accession=acc, protein=aa, classes=";".join(cls),
                               source_group=src_file.get(rid, ""),
                               note=f"split into {len(orfs)} ORF(s): {rid}#orf1.." if orfs else "",
                               nucleotide=nt, protein_source="translated",
                               issues=issues))
    file_counts = Counter(src_file.get(r["entry_id"], "?") for r in records if not r["parent"])
    inventory = dict(source_categories={f"file: {k}.fsa": v for k, v in sorted(file_counts.items())},
                     included_keys=[f"file: {k}.fsa" for k in file_counts],
                     included=sum(not r["parent"] for r in records),
                     class_source="phenotypes.txt" if pheno_class else "file names",
                     operon_like=sum("internal_stop" in r["issues"] for r in records),
                     orf_split_parents=len({r["parent"] for r in records if r["parent"]}),
                     orf_entries=sum(bool(r["parent"]) for r in records))
    return records, version, inventory


def parse_disinfinder(raw_dir):
    """DisinFinder (CGE): disinfectant and biocide resistance genes; the format is
    the same as ResFinder (*.fsa + phenotypes.txt)."""
    return parse_resfinder(raw_dir, db="disinfinder")


PARSERS = {"card": parse_card, "amrfinder": parse_amrfinder, "resfinder": parse_resfinder,
           "disinfinder": parse_disinfinder}


# ------------------------------------------------------------ PointFinder (inventory only)
def inventory_pointfinder(raw_dir, top_n=5):
    """PointFinder inventory for the report. NOT part of the mapping: point mutations cannot
    be compared by sequence identity.

    Layout: one directory per species with resistens-overview.txt (header: #Gene_ID ...,
    Res_codon with comma-separated alternatives, Required_mut for combination-dependent
    mutations), genes.txt and RNA_genes.txt (lists of reference genes).
    """
    raw = Path(raw_dir)
    vfile = raw / "VERSION.argeement"
    if vfile.exists():
        version = " ".join(vfile.read_text().split())
    else:   # direct run without the Snakemake fetch step
        v = (raw / "VERSION").read_text().split()[0] if (raw / "VERSION").exists() else "nover"
        commit = "local"
        if (raw / ".git").exists():
            import subprocess
            r = subprocess.run(["git", "-C", str(raw), "rev-parse", "--short", "HEAD"],
                               capture_output=True, text=True)
            commit = r.stdout.strip() or "local"
        version = f"{v} (git {commit})"

    def gene_list(path):
        if not path.exists():
            return []
        return [ln.strip() for ln in path.read_text(errors="replace").splitlines()
                if ln.strip() and not ln.startswith("#")]

    species = []
    for d in sorted(p for p in raw.iterdir() if p.is_dir() and not p.name.startswith(".")):
        ov = d / "resistens-overview.txt"
        if not ov.exists():
            continue
        header, rows = None, []
        for ln in ov.read_text(errors="replace").splitlines():
            if not ln.strip():
                continue
            if ln.startswith("#"):
                if header is None and "gene_id" in ln.lower():
                    header = [h.strip().lstrip("#").lower() for h in ln.split("\t")]
                continue
            if header:
                f = ln.split("\t")
                rows.append(dict(zip(header, f + [""] * (len(header) - len(f)))))
        if header is None:
            warn(f"PointFinder: header of {d.name}/resistens-overview.txt not recognised")
            continue
        antibiotics = Counter(a.strip() for r in rows for a in r.get("resistance", "").split(",")
                              if a.strip())
        species.append(dict(
            species=d.name,
            reference_genes=len(gene_list(d / "genes.txt")),
            rna_genes=len(gene_list(d / "RNA_genes.txt")),
            genes_with_mutations=len({r.get("gene_id", "") for r in rows}),
            mutation_positions=len(rows),
            resistance_alleles=sum(len([x for x in r.get("res_codon", "").split(",") if x.strip()])
                                   for r in rows),
            combination_dependent=sum(bool(r.get("required_mut", "").strip()) for r in rows),
            top_antibiotics=[f"{a} ({n})" for a, n in antibiotics.most_common(top_n)]))
    if not species:
        raise FileNotFoundError(f"PointFinder: no species directory with resistens-overview.txt in {raw}")
    pf_mut = {}
    for x in species:
        g = pf_mut.setdefault(genus_of(x["species"]), {"genes": 0, "mutations": 0})
        g["genes"] += x["genes_with_mutations"]
        g["mutations"] += x["mutation_positions"]
    keys = ("reference_genes", "rna_genes", "mutation_positions", "resistance_alleles",
            "combination_dependent")
    return dict(version=version, species=species, mutations_by_genus=dict(sorted(pf_mut.items())),
                totals={k: sum(s[k] for s in species) for k in keys} | {"species": len(species)})
