# ARGeement – download, parse, build
# Run: ./update.sh [profile file]   or directly:
#   snakemake --cores 4 [--configfile profiles/x.yaml] --config run_id=test
import os
import yaml

configfile: "config.yaml"

PROFILE = config.get("profile", "current")
RUN = config["run_id"]
WORK = f"work/{PROFILE}/{RUN}"
RAW = f"{WORK}/raw"
STD = f"{WORK}/std"
REL = f"releases/{PROFILE}/{RUN}"
DB = config["databases"]
DBS = ["card", "amrfinder", "resfinder"]
if (DB.get("disinfinder") or {}).get("enabled"):
    DBS.append("disinfinder")

# the effective (merged) configuration: passed to every step and stored in the release
CFG = f"{WORK}/config.effective.yaml"
os.makedirs(WORK, exist_ok=True)
with open(CFG, "w") as fh:
    yaml.safe_dump(dict(config), fh, allow_unicode=True, sort_keys=False)


def opt(db, key):
    return (DB.get(db) or {}).get(key) or ""


rule all:
    input: f"{REL}/manifest.json"


rule fetch_card:
    output: directory(f"{RAW}/card")
    log: f"{WORK}/logs/fetch_card.log"
    params:
        local=opt("card", "local"),
        url=opt("card", "url")
    shell:
        """
        mkdir -p {output}
        if [ -n "{params.local}" ]; then
            if [ -d "{params.local}" ]; then cp -r "{params.local}"/. {output}/
            else tar -xjf "{params.local}" -C {output}; fi
        else
            curl -fsSL -o {output}/card.tar.bz2 "{params.url}" 2> {log}
            tar -xjf {output}/card.tar.bz2 -C {output}
            rm {output}/card.tar.bz2
        fi
        # if the archive unpacks into a subdirectory, the required files are moved up
        for f in card.json protein_fasta_protein_homolog_model.fasta \
                 nucleotide_fasta_protein_homolog_model.fasta; do
            [ -e {output}/$f ] || find {output} -name $f -exec mv {{}} {output}/ \\;
        done
        test -s {output}/protein_fasta_protein_homolog_model.fasta
        """


rule fetch_amrfinder:
    output: directory(f"{RAW}/amrfinder")
    log: f"{WORK}/logs/fetch_amrfinder.log"
    params:
        local=opt("amrfinder", "local"),
        url=opt("amrfinder", "url"),
        files=" ".join(["version.txt",
                        "ReferenceGeneCatalog.txt", "ReferenceGeneCatalog.tsv",
                        "ReferenceGeneHierarchy.txt", "ReferenceGeneHierarchy.tsv",
                        "fam.tsv", "fam.tab", "AMRProt.fa", "AMRProt", "AMR_CDS.fa", "AMR_CDS",
                        "AMRProt-mutation.tsv", "AMRProt-mutation.tab", "taxgroup.tsv", "taxgroup.tab"])
    shell:
        """
        mkdir -p {output}
        for f in {params.files}; do
            if [ -n "{params.local}" ]; then
                if [ -e "{params.local}/$f" ]; then cp "{params.local}/$f" {output}/; fi
                for d in "{params.local}"/AMR_DNA-*.tsv "{params.local}"/AMR_DNA-*.tab; do
                    if [ -e "$d" ]; then cp "$d" {output}/; fi
                done
            else
                curl -fsSL -o {output}/$f "{params.url}/$f" 2>> {log} || rm -f {output}/$f
            fi
        done
        # rRNA mutations: AMR_DNA-<taxgroup>.tsv for the taxonomic groups that have a nucleotide
        # reference in the taxgroup table (column 3 > 0)
        if [ -z "{params.local}" ]; then
            for tg in $({{ cat {output}/taxgroup.tsv {output}/taxgroup.tab 2>/dev/null || true; }} | awk -F'\t' '$1 !~ /^#/ && $3 > 0 {{print $1}}'); do
                curl -fsSL -o {output}/AMR_DNA-$tg.tsv "{params.url}/AMR_DNA-$tg.tsv" 2>> {log} || rm -f {output}/AMR_DNA-$tg.tsv
            done
        fi
        test -s {output}/ReferenceGeneCatalog.txt || test -s {output}/ReferenceGeneCatalog.tsv
        test -s {output}/AMRProt.fa || test -s {output}/AMRProt
        """


# git-based databases (ResFinder, DisinFinder, PointFinder): remote or local repository, optional ref
rule fetch_git_db:
    output: directory(f"{RAW}/{{gdb}}")
    wildcard_constraints: gdb="resfinder|pointfinder|disinfinder"
    log: f"{WORK}/logs/fetch_{{gdb}}.log"
    params:
        local=lambda wc: opt(wc.gdb, "local"),
        url=lambda wc: opt(wc.gdb, "url"),
        ref=lambda wc: opt(wc.gdb, "ref")
    shell:
        """
        src="{params.local}"; [ -n "$src" ] || src="{params.url}"
        if [ -d "$src" ] && [ ! -d "$src/.git" ]; then
            # plain local directory (not a git repository): the ref cannot be enforced
            if [ -n "{params.ref}" ]; then
                echo "ERROR: the local source of {wildcards.gdb} ($src) is not a git repository, so the ref ({params.ref}) cannot be checked. Copy it together with its .git directory, or leave the local field empty (the ref is then cloned from the remote repository)." | tee {log} >&2
                exit 1
            fi
            mkdir -p {output}; cp -r "$src"/. {output}/
            commit="local"
        else
            if [ -n "{params.ref}" ]; then
                git clone -q "$src" {output} 2> {log}
                git -C {output} checkout -q "{params.ref}" 2>> {log}
            else
                git clone -q --depth 1 "$src" {output} 2> {log}
            fi
            commit=$(git -C {output} rev-parse --short HEAD)
        fi
        if [ ! -f {output}/VERSION.argeement ]; then
            v=$( [ -f {output}/VERSION ] && head -1 {output}/VERSION | tr -d '\\r' || echo "nover" )
            echo "$v (git $commit)" > {output}/VERSION.argeement
        fi
        """


# PointFinder: inventory for the report only, not part of the mapping
rule inventory_pointfinder:
    input: f"{RAW}/pointfinder"
    output: f"{STD}/pointfinder.inventory.json"
    params:
        expect=lambda wc: opt("pointfinder", "expect_version")
    shell:
        "python -m argeement inventory pointfinder --raw {input} --out {output} "
        "--expect-version '{params.expect}'"


rule parse:
    input: f"{RAW}/{{db}}"
    output:
        tsv=f"{STD}/{{db}}.tsv",
        ver=f"{STD}/{{db}}.version"
    params:
        expect=lambda wc: opt(wc.db, "expect_version")
    shell:
        "python -m argeement parse {wildcards.db} --raw {input} --out {output.tsv} "
        "--config {CFG} --expect-version '{params.expect}'"


PF_ENABLED = bool((DB.get("pointfinder") or {}).get("enabled"))


rule build:
    input:
        std=expand(f"{STD}/{{db}}.tsv", db=DBS),
        extra=[f"{STD}/pointfinder.inventory.json"] if PF_ENABLED else []
    output: f"{REL}/manifest.json"
    threads: config["threads"]
    params:
        extra_inv="--extra-inventory pointfinder=" + f"{STD}/pointfinder.inventory.json" if PF_ENABLED else ""
    shell:
        "python -m argeement build --std {input.std} --out {REL} --config {CFG} "
        "--mmseqs {config[mmseqs]} --threads {threads} --profile {PROFILE} {params.extra_inv} "
        "&& cp {CFG} {REL}/config.yaml"
