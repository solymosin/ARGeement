"""Database content report (Markdown) from a release, for public documentation.

Every number comes from the release files, so the report can be regenerated at each update and is
always consistent with the mapping.
"""
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from .build import RELATION_RANK, STRONG_MAX_RANK, load_class_map
from .names import CATEGORIES, CONFLICTS, SEQ_RELATIONS

DBS = ["amrfinder", "card", "resfinder", "disinfinder"]
NAMES = {"amrfinder": "AMRFinderPlus", "card": "CARD", "resfinder": "ResFinder",
         "disinfinder": "DisinFinder"}

T = {
    "en": dict(
        title="Content of the AMR reference databases",
        intro=("This document describes the antimicrobial resistance (AMR) reference databases "
               "used in this ARGeement release: what types of determinants they contain, how many, "
               "how they are distributed across drug classes, and how they overlap. It is generated "
               "automatically from the release files, so all numbers are consistent with the "
               "sequence-based mapping."),
        release="Release", profile="Profile", created="Created", versions="Database versions",
        database="Database", version="Version",
        s_overview="1. Overview",
        ov_cols=["Database", "Version", "Entries in source", "Included in mapping",
                 "ORF sub-entries", "Usable (protein)", "Unique protein sequences"],
        ov_note=("*Entries in source* counts the categories listed in section 2 (catalog rows, CARD "
                 "models or sequence entries); *Included* is the subset used for the mapping; "
                 "*Usable* excludes entries that could not be used at protein level (see section 6). *ORF "
                 "sub-entries* are genes obtained by splitting entries that are not a single "
                 "coding sequence (e.g. operons); they are counted among the usable entries."),
        s_types="2. Types of determinants",
        ty_cols=["Category in source", "Count", "Included"],
        ty_amr=("AMRFinderPlus catalog rows by `type/subtype`. Point mutations (`POINT`), gene "
                "disruptions (`POINT_DISRUPT`) and susceptibility markers (`AMR-SUSCEPTIBLE`) are "
                "mutation- or context-dependent and cannot be compared by sequence identity. Of the "
                "stress-response genes only biocide resistance genes (`STRESS/BIOCIDE`) are "
                "included, for comparison with DisinFinder and the CARD class *disinfecting agents "
                "and antiseptics*; metal, acid and heat stress and virulence genes are outside the "
                "scope of this mapping."),
        ty_dis=("DisinFinder sequence entries by source file. All entries are genes conferring "
                "resistance or reduced susceptibility to disinfectants and biocides."),
        ty_card=("CARD models by model type. Only protein homolog models (presence of a gene) are "
                 "included; variant, overexpression and knockout models depend on mutations or "
                 "expression and are excluded."),
        ty_res=("ResFinder sequence entries by source file. Chromosomal point mutations are kept in "
                "the separate PointFinder database, which is not part of this mapping."),
        scope="Included AMRFinderPlus entries by scope",
        pf_title="PointFinder (point mutations, not included)",
        pf_intro=("The PointFinder database (version `{v}`) is listed for completeness: it holds "
                  "species-specific chromosomal point mutations, which are not part of the "
                  "sequence-based mapping. *Combination-dependent* mutations confer resistance only "
                  "together with other mutations."),
        pf_cols=["Species", "Reference genes", "RNA genes", "Genes with mutations",
                 "Mutation positions", "Resistance alleles", "Combination-dependent",
                 "Most frequent antibiotics"],
        pf_total="Total",
        mg_title="Mutation-based determinants by genus (not included in the mapping)",
        mg_intro=("Number of genes and mutations per genus in the mutation-based part of each "
                  "source, for comparison of their taxonomic coverage. The units differ by "
                  "source: AMRFinderPlus counts genes and listed mutations ({amr_unit}); CARD counts "
                  "variant models and the SNPs listed in them; PointFinder counts genes with "
                  "mutations and mutation positions. A dash means no data for that genus."),
        mg_cols=["Genus", "AMRFinderPlus genes", "AMRFinderPlus mutations", "CARD models",
                 "CARD SNPs", "PointFinder genes", "PointFinder mutations"],
        mg_genera="genera",
        mech="Resistance mechanisms of the included CARD entries (an entry may have several)",
        s_classes="3. Drug classes",
        cl_intro=("The databases use different drug class vocabularies (e.g. AMRFinderPlus "
                  "`BETA-LACTAM`, CARD `penam`, `cephalosporin`, `carbapenem`). Native terms were "
                  "mapped to common classes with the table in `argeement/class_map.yaml`. An entry "
                  "is counted in every class it is assigned to, so column totals can exceed the "
                  "number of entries. *Only in X* counts entries of database X with no strong match "
                  "(identical, near-identical or length variant) in any of the other databases."),
        cl_cols_only="Only in", no_class="(no class assigned)",
        scope_note=("DisinFinder covers only disinfectant and biocide resistance, and ResFinder "
                    "only antimicrobial (antibiotic) resistance; the two complement each other, so "
                    "classes outside a database's scope are not reported as gaps."),
        native="Native class terms per database",
        unmapped="Terms not covered by the harmonization table",
        unmapped_none="All native class terms are covered by the harmonization table.",
        ignored="Terms ignored by the harmonization table",
        ignored_note=("These terms are deliberately not counted as classes (see the comments in "
                      "`class_map.yaml`); the table shows which gene families carry them, so that "
                      "a change after a database update is visible."),
        families_col="Most frequent gene families",
        unmapped_note=("These entries are listed above under their other (mapped) classes, or under "
                       "*(no class assigned)*. Adding the terms to `class_map.yaml` will include them."),
        s_overlap="4. Overlap between the databases",
        s_names="5. Name concordance of identical proteins",
        nm_intro=("Entries of different databases that encode exactly the same protein are "
                  "compared by name. Each pair (counted once) falls into the first matching "
                  "category: *exact*; *case* (differs only in letter case); *formatting* (identical "
                  "without the `bla` prefix and punctuation, e.g. `TEM-1` / `blaTEM-1`); "
                  "*species prefix* (CARD species prefix, e.g. `Ecol_emrE` / `emrE`); *variant "
                  "suffix* (ResFinder/DisinFinder letter for synonymous variants, e.g. `blaTEM-1A` / "
                  "`blaTEM-1`); *NCBI synonym* (listed as a synonym in the NCBI gene hierarchy); "
                  "*allele number differs* (same gene name, different allele number) and *different "
                  "name*. The last two are genuine naming conflicts: the same protein carries "
                  "incompatible names, so name-based comparison of results fails even where the "
                  "sequences are identical."),
        nm_cols=["Pair", "Pairs", "Exact", "Case", "Formatting", "Species prefix",
                 "Variant suffix", "NCBI synonym", "Allele no. differs", "Different name",
                 "Name not identical¹"],
        nm_foot=("¹ share of the pairs whose names are not identical even ignoring letter case, "
                 "i.e. where a naive name-based comparison fails."),
        nm_ex="Examples of naming conflicts (the full list is in `name_concordance.tsv`)",
        nm_cat={"allele_number_differs": "allele number differs", "different_name": "different name"},
        nm_sub1="Same protein, different name",
        nm_sub2="Same name, different protein",
        sn_intro=("The reverse direction: entries of different databases whose names are equivalent "
                  "(identical after removing formatting, CARD species prefixes and ResFinder "
                  "variant letters) are compared by sequence, using all protein hits of the release "
                  "(*unrelated*: no similarity above the search threshold). Many gene names denote a "
                  "gene with several variants (e.g. `sul1`), so different sequences are expected "
                  "there; but the name-based comparison of results silently treats them as equal. "
                  "Where the name is an *allele* name in the NCBI gene hierarchy, a different "
                  "protein is a genuine conflict."),
        sn_cols=["Pair", "Pairs", "Identical", "Near-identical", "Length variant", "Same family",
                 "Homolog", "Unrelated", "Not identical"],
        sn_allele=("**Allele names:** {n} pairs where the AMRFinderPlus name is an allele; {k} of "
                   "them ({p}) do not encode the identical protein."),
        sn_ex_allele="Allele-name conflicts (most divergent first; full list in `same_name_concordance.tsv`)",
        sn_ex_far="Same name, distant or unrelated proteins",
        olp_intro=("Share of the usable entries of the row database that have a match in the column "
                   "database. *Identical* = identical protein; *strong* = identical, near-identical "
                   "(≥{np}% identity, ≥{nc} coverage) or length variant."),
        olp_cols=["From \\ to", "Identical", "Strong", "Same family", "No strong match"],
        uniq="Entries without a strong match in any of the other databases",
        s_sw="6. Strengths, limitations and differences",
        sw_auto="Observations from the data",
        sw_unique=("**{db}** has {n} {ents} ({p}) without a strong match in any of the other "
                   "databases; most of them in: {top}."),
        sw_unique_none="Every usable **{db}** entry has a strong match in at least one other database.",
        sw_gap=("**{db}** has no entries in the class(es) {cls}, where the other databases have at "
                "least {k} entries."),
        sw_depth=("**{db}**: {n} usable {ents} in {f} gene {fams} "
                  "(median {m} per family; largest: {big})."),
        sw_general="General characteristics",
        general={
            "amrfinder": ("NCBI's curated Reference Gene Database. Protein-based; contains acquired and "
                          "intrinsic resistance genes, taxon-specific point mutations and, in the "
                          "'plus' scope, stress-response and virulence genes. A curated gene "
                          "hierarchy with HMMs allows family-level calls for novel variants, so "
                          "AMRFinderPlus can report a family symbol (e.g. `blaTEM`) instead of an allele."),
            "card": ("Built on the Antibiotic Resistance Ontology (ARO). Uses several model types "
                     "(protein homolog, protein variant, rRNA variant, overexpression, knockout); "
                     "each determinant is annotated with drug classes, resistance mechanism and AMR "
                     "gene family, which makes CARD the richest source of functional annotation."),
            "disinfinder": ("Companion database of ResFinder from the Center for Genomic "
                            "Epidemiology for genes conferring resistance to disinfectants and "
                            "biocides (e.g. quaternary ammonium compounds, aldehydes), stored as "
                            "nucleotide sequences."),
            "resfinder": ("Database of acquired resistance genes from the Center for Genomic "
                          "Epidemiology, stored as nucleotide sequences and designed for assemblies "
                          "and reads. Several alleles may encode the same protein (synonymous "
                          "variants), which is visible here as entries sharing a protein sequence."),
        },
        sw_mapping="Limitations of this mapping",
        lim=[
            "Mutation-based determinants are excluded from all databases; the mapping covers "
            "gene presence only.",
            "Sequence identity does not guarantee functional identity: two near-identical "
            "beta-lactamases can have different substrate spectra. Phenotype interpretation remains "
            "with the source databases.",
            "In allele-rich families, *near-identical* means the closest, but a different, allele.",
            "Sequence-based family clusters are built around representative sequences, so "
            "related families are not chained together; still, {merged} of {clusters} clusters "
            "join several NCBI families, where the curated families are narrower than the "
            "identity threshold. For family-level reporting, the NCBI family (`family_label`, "
            "inferred for other databases from the closest AMRFinderPlus protein) is the "
            "recommended label.",
            "{split_list} could not be translated as a single coding sequence (mostly multi-gene "
            "operons such as `VanHAX`, or sequences in another reading frame); they were split "
            "into open reading frames (`#orf<k>` sub-entries), which take part in the mapping. "
            "In operon-rich classes ({operon_classes}) the entry counts therefore refer to genes, "
            "not to the database's own entries.",
            "{unusable_list} remain unusable at protein level; they are listed in "
            "`sequence_issues.tsv`.",
        ],
        s_methods="7. Methods",
        methods=("Protein sequences were compared all-against-all with MMseqs2 ({mm}). Thresholds: "
                 "near-identical ≥{np}% identity and ≥{nc} coverage of both sequences; same family "
                 "≥{fp}% identity and ≥{fc} coverage. ResFinder nucleotide sequences were translated "
                 "with the bacterial genetic code (table 11). Full parameters are in `manifest.json` "
                 "and `config.yaml` of the release."),
        generated="Generated by ARGeement {v}.",
    ),
}


# declared scope: gaps are only reported for drug classes within the scope of a database
BIOCIDE_LIKE = {"Biocide/antiseptic", "Heat stress", "Fluoride"}
IN_SCOPE = {
    "disinfinder": lambda c: c in BIOCIDE_LIKE,
    "resfinder": lambda c: c not in BIOCIDE_LIKE,
}


def plural(n, one, many, lang="en"):
    return one if n == 1 else many


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def pct(a, b, lang="en"):
    if not b:
        return "–"
    v = f"{100 * a / b:.1f}%"
    return v


def split(s):
    return [x for x in str(s).split(";") if x]


def report(release, out=None, lang="en", gap_min=5, top_n=3):
    t = T[lang]
    P = lambda a, b: pct(a, b, lang)  # noqa: E731
    rel = Path(release)
    man = json.loads((rel / "manifest.json").read_text())
    e = pd.read_csv(rel / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    m = pd.read_csv(rel / "mappings.tsv", sep="\t", dtype=str, keep_default_na=False)
    inv = {db: json.loads((rel / f"inventory_{db}.json").read_text())
           for db in DBS if (rel / f"inventory_{db}.json").exists()}
    for col in ("class_h", "classes", "mechanism", "unmapped_classes", "ncbi_family", "family",
                "parent"):
        if col not in e.columns:
            e[col] = ""
    e["usable_b"] = e["usable"] == "True"
    u = e[e["usable_b"]]
    p = man["params"]
    dbs = [d for d in DBS if d in set(e["db"])]
    L = []

    # ---------------------------------------------------------------- header
    L += [f"# {t['title']}", "", t["intro"], ""]
    L += [f"- **{t['release']}:** `{rel.resolve().name}`"]
    if man.get("profile"):
        L += [f"- **{t['profile']}:** `{man['profile']}`"]
    L += [f"- **{t['created']}:** {man['created']}", ""]
    L += [md_table([t["database"], t["version"]],
                   [[NAMES[d], f"`{man['database_versions'].get(d, '?')}`"] for d in dbs]), ""]

    # ---------------------------------------------------------------- 1. overview
    L += [f"## {t['s_overview']}", ""]
    rows = []
    for d in dbs:
        sub = e[e["db"] == d]
        src = sum(inv.get(d, {}).get("source_categories", {}).values()) or len(sub)
        n_orf = int((sub["parent"] != "").sum())
        rows.append([NAMES[d], f"`{man['database_versions'].get(d, '?')}`", src, len(sub) - n_orf,
                     n_orf, int(sub["usable_b"].sum()), sub.loc[sub["usable_b"], "prot_id"].nunique()])
    L += [md_table(t["ov_cols"], rows), "", t["ov_note"], ""]

    # ---------------------------------------------------------------- 2. types
    L += [f"## {t['s_types']}", ""]
    notes = {"amrfinder": t["ty_amr"], "card": t["ty_card"], "resfinder": t["ty_res"],
             "disinfinder": t["ty_dis"]}
    for d in dbs:
        L += [f"### {NAMES[d]}", "", notes[d], ""]
        cats = inv.get(d, {}).get("source_categories", {})
        if cats:
            rows = []
            inc_keys = inv.get(d, {}).get("included_keys")
            for k, v in cats.items():
                if inc_keys is not None:
                    inc = "✓" if k in inc_keys else ""
                elif d == "amrfinder":
                    ts = k.split(": ", 1)[-1]
                    typ, _, sub_ = ts.partition("/")
                    inc = "✓" if typ == "AMR" and sub_ not in ("POINT", "POINT_DISRUPT", "AMR-SUSCEPTIBLE") else ""
                elif d == "card":
                    inc = "✓" if k.endswith("protein homolog model") else ""
                else:
                    inc = "✓"
                rows.append([f"`{k.split(': ', 1)[-1]}`", v, inc])
            L += [md_table(t["ty_cols"], rows), ""]
        if d == "amrfinder" and inv.get(d, {}).get("included_scope"):
            L += [f"*{t['scope']}:* " + ", ".join(f"{k}: {v}" for k, v in
                                                  sorted(inv[d]["included_scope"].items())), ""]
        if d == "resfinder" and (rel / "inventory_pointfinder.json").exists():
            pf = json.loads((rel / "inventory_pointfinder.json").read_text())
            rows = [[f"*{x['species'].replace('_', ' ')}*", x["reference_genes"], x["rna_genes"],
                     x["genes_with_mutations"], x["mutation_positions"], x["resistance_alleles"],
                     x["combination_dependent"], ", ".join(x["top_antibiotics"][:3])]
                    for x in pf["species"]]
            tt = pf["totals"]
            rows.append([f"**{t['pf_total']}** ({tt['species']})", tt["reference_genes"], tt["rna_genes"],
                         "", tt["mutation_positions"], tt["resistance_alleles"],
                         tt["combination_dependent"], ""])
            L += [f"#### {t['pf_title']}", "", t["pf_intro"].format(v=pf["version"]), "",
                  md_table(t["pf_cols"], rows), ""]
        if d == "card":
            mc = Counter(x for s in e.loc[e["db"] == "card", "mechanism"] for x in split(s))
            if mc:
                L += [f"*{t['mech']}:*", "", md_table([t["ty_cols"][0], t["ty_cols"][1]],
                                                       [[k, v] for k, v in mc.most_common()]), ""]

    # ---------------------------------------------------------------- mutations per genus
    mut_src = {"amrfinder": inv.get("amrfinder", {}).get("mutations_by_genus"),
               "card": inv.get("card", {}).get("mutations_by_genus")}
    if (rel / "inventory_pointfinder.json").exists():
        mut_src["pointfinder"] = json.loads((rel / "inventory_pointfinder.json").read_text()) \
            .get("mutations_by_genus")
    mut_src = {k: v for k, v in mut_src.items() if v}
    if mut_src:
        genera = sorted({g for v in mut_src.values() for g in v},
                        key=lambda g: (-sum(v.get(g, {}).get("mutations", 0) for v in mut_src.values()), g))
        cell = lambda src, g, k: (mut_src[src][g][k] if src in mut_src and g in mut_src[src] else "–")
        rows = [[f"*{g}*"] + [cell(src, g, k) for src in ("amrfinder", "card", "pointfinder")
                              for k in ("genes", "mutations")] for g in genera]
        tot = lambda src, k: sum(v[k] for v in mut_src[src].values()) if src in mut_src else "–"
        rows.append([f"**{t['pf_total']}** ({len(genera)} {t['mg_genera']})"]
                    + [tot(src, k) for src in ("amrfinder", "card", "pointfinder") for k in ("genes", "mutations")])
        amr_unit = inv.get("amrfinder", {}).get("mutation_unit", "–")
        L += [f"### {t['mg_title']}", "", t["mg_intro"].format(amr_unit=amr_unit), "",
              md_table(t["mg_cols"], rows), ""]

    # ---------------------------------------------------------------- close counterparts per entry
    strong = m[m["relation"].map(RELATION_RANK) <= STRONG_MAX_RANK]
    strong_to = defaultdict(set)
    for da, ea, db_ in zip(strong["db_a"], strong["entry_a"], strong["db_b"]):
        strong_to[(da, ea)].add(db_)
    u = u.assign(only=[len(strong_to.get((d, x), set())) == 0 for d, x in zip(u["db"], u["entry_id"])])

    # ---------------------------------------------------------------- 3. drug classes
    L += [f"## {t['s_classes']}", "", t["cl_intro"], ""]
    cls_count = {d: Counter() for d in dbs}
    cls_only = {d: Counter() for d in dbs}
    for d, ch, only in zip(u["db"], u["class_h"], u["only"]):
        cl = split(ch) or [t["no_class"]]
        for c in cl:
            cls_count[d][c] += 1
            if only:
                cls_only[d][c] += 1
    all_cls = sorted({c for d in dbs for c in cls_count[d]},
                     key=lambda c: (c == t["no_class"], -sum(cls_count[d][c] for d in dbs), c))
    header = [""] + [NAMES[d] for d in dbs] + [f"{t['cl_cols_only']} {NAMES[d]}" for d in dbs]
    rows = [[c] + [cls_count[d][c] for d in dbs] + [cls_only[d][c] for d in dbs] for c in all_cls]
    L += [md_table(header, rows), ""]

    L += [f"<details><summary>{t['native']}</summary>", ""]
    for d in dbs:
        nc = Counter(x for s in u.loc[u["db"] == d, "classes"] for x in split(s))
        L += [f"**{NAMES[d]}**", "", md_table([t["ty_cols"][0], t["ty_cols"][1]],
                                              [[f"`{k}`", v] for k, v in nc.most_common()]), ""]
    L += ["</details>", ""]

    L += [f"### {t['unmapped']}", ""]
    um = Counter((d, x) for d, s in zip(u["db"], u["unmapped_classes"]) for x in split(s))
    if um:
        L += [md_table([t["database"], t["ty_cols"][0], t["ty_cols"][1]],
                       [[NAMES[d], f"`{k}`", v] for (d, k), v in um.most_common()]),
              "", t["unmapped_note"], ""]
    else:
        L += [t["unmapped_none"], ""]

    term, _ = load_class_map()
    ign = defaultdict(Counter)
    for d, cs, fam in zip(u["db"], u["classes"], u["family"]):
        for c in split(cs):
            if c.strip().lower() in term and term[c.strip().lower()] is None:
                ign[(d, c.strip())][fam or "–"] += 1
    if ign:
        rows = [[NAMES[d], f"`{c}`", sum(fc.values()),
                 ", ".join(f"{f} ({n})" for f, n in fc.most_common(3))]
                for (d, c), fc in sorted(ign.items(), key=lambda kv: -sum(kv[1].values()))]
        L += [f"### {t['ignored']}", "",
              md_table([t["database"], t["ty_cols"][0], t["ty_cols"][1], t["families_col"]], rows),
              "", t["ignored_note"], ""]

    # ---------------------------------------------------------------- 4. overlap
    L += [f"## {t['s_overlap']}", "",
          t["olp_intro"].format(np=p["near_identical_pident"], nc=p["near_identical_cov"]), ""]
    rank = m["relation"].map(RELATION_RANK)
    best = m.assign(_r=rank).groupby(["db_a", "entry_a", "db_b"])["_r"].min().reset_index()
    rows = []
    for a in dbs:
        n = int((u["db"] == a).sum())
        for b in dbs:
            if a == b:
                continue
            br = best[(best["db_a"] == a) & (best["db_b"] == b)]["_r"]
            ident = int((br == 0).sum())
            strong_n = int((br <= STRONG_MAX_RANK).sum())
            fam = int((br == RELATION_RANK["same_family"]).sum())
            rows.append([f"{NAMES[a]} → {NAMES[b]}", P(ident, n), P(strong_n, n), P(fam, n),
                         P(n - strong_n, n)])
    L += [md_table(t["olp_cols"], rows), ""]
    L += [f"*{t['uniq']}:* " + ", ".join(
        f"{NAMES[d]} {int(u.loc[u['db'] == d, 'only'].sum())} "
        f"({P(int(u.loc[u['db'] == d, 'only'].sum()), int((u['db'] == d).sum()))})" for d in dbs), ""]

    # ---------------------------------------------------------------- 5. name concordance
    ncf = rel / "name_concordance.tsv"
    if ncf.exists():
        nc = pd.read_csv(ncf, sep="\t", dtype=str, keep_default_na=False)
        L += [f"## {t['s_names']}", "", f"### {t['nm_sub1']}", "", t["nm_intro"], ""]
        rows = []
        for (a, b), g in nc.groupby(["db_a", "db_b"], sort=False):
            cc = g["category"].value_counts()
            n = len(g)
            naive_fail = n - int(cc.get("exact", 0)) - int(cc.get("case", 0))
            rows.append([f"{NAMES.get(a, a)} – {NAMES.get(b, b)}", n]
                        + [int(cc.get(c, 0)) for c in CATEGORIES] + [P(naive_fail, n)])
        rows.sort(key=lambda r: -r[1])
        L += [md_table(t["nm_cols"], rows), "", t["nm_foot"], ""]
        conf = nc[nc["category"].isin(CONFLICTS)]
        if len(conf):
            ex = (conf.drop_duplicates(["db_a", "name_a", "db_b", "name_b"])
                      .groupby("category", group_keys=False).head(12))
            L += [f"*{t['nm_ex']}:*", "",
                  md_table([t["nm_cols"][0], "A", "B", ""],
                           [[f"{NAMES.get(r.db_a, r.db_a)} – {NAMES.get(r.db_b, r.db_b)}",
                             f"`{r.name_a}`", f"`{r.name_b}`", t["nm_cat"][r.category]]
                            for r in ex.itertuples(index=False)]), ""]

    snf = rel / "same_name_concordance.tsv"
    if snf.exists():
        sn = pd.read_csv(snf, sep="\t", dtype=str, keep_default_na=False)
        L += [f"### {t['nm_sub2']}", "", t["sn_intro"], ""]
        rows = []
        for (a, b), g in sn.groupby(["db_a", "db_b"], sort=False):
            cc = g["relation"].value_counts()
            rows.append([f"{NAMES.get(a, a)} – {NAMES.get(b, b)}", len(g)]
                        + [int(cc.get(c, 0)) for c in SEQ_RELATIONS]
                        + [P(len(g) - int(cc.get("identical_protein", 0)), len(g))])
        rows.sort(key=lambda r: -r[1])
        L += [md_table(t["sn_cols"], rows), ""]
        al = sn[sn["name_level"] == "allele"]
        if len(al):
            bad = al[al["relation"] != "identical_protein"]
            L += [t["sn_allele"].format(n=len(al), k=len(bad), p=P(len(bad), len(al))), ""]
            if len(bad):
                order = {r: i for i, r in enumerate(reversed(SEQ_RELATIONS))}
                ex = bad.assign(_o=bad["relation"].map(order),
                                _n=pd.to_numeric(bad["n_diff"], errors="coerce").fillna(1e9)) \
                        .sort_values(["_o", "_n"], ascending=[True, False]) \
                        .drop_duplicates(["db_a", "name_a", "db_b", "name_b"]).head(15)
                L += [f"*{t['sn_ex_allele']}:*", "",
                      md_table([t["nm_cols"][0], "A", "B", "relation", "n_diff"],
                               [[f"{NAMES.get(r.db_a, r.db_a)} – {NAMES.get(r.db_b, r.db_b)}",
                                 f"`{r.name_a}`", f"`{r.name_b}`", r.relation, r.n_diff or "–"]
                                for r in ex.itertuples(index=False)]), ""]
        far = sn[(sn["name_level"] != "allele") & sn["relation"].isin(["homolog", "unrelated"])]
        if len(far):
            ex = far.drop_duplicates(["db_a", "name_a", "db_b", "name_b"]).head(12)
            L += [f"*{t['sn_ex_far']}:*", "",
                  md_table([t["nm_cols"][0], "A", "B", "relation", "n_diff"],
                           [[f"{NAMES.get(r.db_a, r.db_a)} – {NAMES.get(r.db_b, r.db_b)}",
                             f"`{r.name_a}`", f"`{r.name_b}`", r.relation, r.n_diff or "–"]
                            for r in ex.itertuples(index=False)]), ""]

    # ---------------------------------------------------------------- 6. strengths and limitations
    L += [f"## {t['s_sw']}", "", f"### {t['sw_auto']}", ""]
    for d in dbs:
        n_only = int(u.loc[u["db"] == d, "only"].sum())
        n_all = int((u["db"] == d).sum())
        top = [c for c, _ in cls_only[d].most_common() if c != t["no_class"]][:top_n]
        top_s = ", ".join(f"{c} ({cls_only[d][c]})" for c in top) or "–"
        if n_only == 0:
            L += ["- " + t["sw_unique_none"].format(db=NAMES[d])]
        else:
            L += ["- " + t["sw_unique"].format(db=NAMES[d], n=n_only, p=P(n_only, n_all), top=top_s,
                                                ents=plural(n_only, "entry", "entries", lang))]
    for d in dbs:
        in_scope = IN_SCOPE.get(d, lambda c: True)
        others = [o for o in dbs if o != d and IN_SCOPE.get(o, lambda c: True)]
        gaps = [c for c in all_cls if c != t["no_class"] and in_scope(c) and cls_count[d][c] == 0
                and others and all(cls_count[o][c] >= gap_min for o in others
                                   if IN_SCOPE.get(o, lambda c: True)(c))]
        if gaps:
            L += ["- " + t["sw_gap"].format(db=NAMES[d], cls=", ".join(gaps), k=gap_min)]
    for d in dbs:
        sub = u[u["db"] == d]
        famcol = "ncbi_family" if d == "amrfinder" and (sub["ncbi_family"] != "").any() else "family"
        f = sub[sub[famcol] != ""][famcol].value_counts()
        if len(f):
            L += ["- " + t["sw_depth"].format(db=NAMES[d], n=len(sub), f=len(f), m=int(f.median()),
                                               ents=plural(len(sub), "entry", "entries", lang),
                                               fams=plural(len(f), "family", "families", lang),
                                               big=", ".join(f"`{k}` ({v})" for k, v in f.head(3).items()))]
    if "disinfinder" in dbs and "resfinder" in dbs:
        L += ["", t["scope_note"]]
    L += ["", f"### {t['sw_general']}", ""]
    L += [f"- **{NAMES[d]}:** {t['general'][d]}" for d in dbs]
    val = man.get("ncbi_validation") or {}
    cge = [d for d in ("resfinder", "disinfinder") if d in dbs]
    split_parents = set(e.loc[e["parent"] != "", "parent"])

    def listing(groups):
        return " and ".join(f"{len(x)} {NAMES[d]} {plural(len(x), 'entry', 'entries', lang)}"
                            for d, x in groups.items() if len(x))

    splitg = {d: e[(e["db"] == d) & e["entry_id"].isin(split_parents)] for d in cge}
    unus = {d: e[(e["db"] == d) & ~e["usable_b"] & ~e["entry_id"].isin(split_parents)] for d in cge}
    split_list, unusable_list = listing(splitg), listing(unus)
    unusable = sum(len(x) for x in unus.values())
    oc = Counter(c for x in splitg.values() for s_ in x["class_h"] for c in split(s_))
    operon_classes = ", ".join(f"{c} ({n})" for c, n in oc.most_common(4)) or "–"
    L += ["", f"### {t['sw_mapping']}", ""]
    for i, s in enumerate(t["lim"]):
        if "{merged}" in s and not val:
            continue
        if "{unusable_list}" in s and not unusable_list:
            continue
        if "{split_list}" in s and not split_list:
            continue
        L += ["- " + s.format(merged=val.get("clusters_merged", "?"), clusters=val.get("clusters", "?"),
                              unusable_list=unusable_list, split_list=split_list,
                              operon_classes=operon_classes)]
    L += [""]

    # ---------------------------------------------------------------- 7. methods
    L += [f"## {t['s_methods']}", "",
          t["methods"].format(mm=man.get("mmseqs_version", "?"), np=p["near_identical_pident"],
                              nc=p["near_identical_cov"], fp=p["same_family_pident"],
                              fc=p["same_family_cov"]), "",
          "---", f"*{t['generated'].format(v=man.get('argeement_version', '?'))}*", ""]

    text = "\n".join(L)
    if out:
        Path(out).write_text(text, encoding="utf-8")
    return text
