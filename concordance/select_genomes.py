"""Genome selection for the real-genome use case: n random complete RefSeq genomes per species,
with a fixed random seed. The selection depends on the NCBI genome set available at the time;
the resulting list (genomes.tsv) is the reference for reproduction.

Usage:
    python concordance/select_genomes.py --n 25 --seed 1 --out concordance/genomes.tsv
Requires the NCBI datasets command-line tool (ncbi-datasets-cli).
"""
import argparse
import json
import random
import subprocess

SPECIES = ["Escherichia coli", "Klebsiella pneumoniae", "Salmonella enterica",
           "Staphylococcus aureus"]


def accessions(species):
    cmd = ["datasets", "summary", "genome", "taxon", species, "--assembly-source", "RefSeq",
           "--assembly-level", "complete", "--as-json-lines"]
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    acc = []
    for line in out.splitlines():
        r = json.loads(line)
        if r.get("accession", "").startswith("GCF_"):
            acc.append(r["accession"])
    return sorted(set(acc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=25, help="number of genomes per species")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--species", nargs="+", default=SPECIES)
    ap.add_argument("--out", default="concordance/genomes.tsv")
    a = ap.parse_args()
    rows = []
    for sp in a.species:
        acc = accessions(sp)
        pick = random.Random(f"{a.seed}:{sp}").sample(acc, min(a.n, len(acc)))
        rows += [(x, sp, len(acc)) for x in sorted(pick)]
        print(f"{sp}: {len(pick)} of {len(acc)} complete genomes")
    with open(a.out, "w") as fh:
        fh.write("accession\tspecies\tpool_size\n")
        fh.writelines(f"{x}\t{s}\t{n}\n" for x, s, n in rows)


if __name__ == "__main__":
    main()
