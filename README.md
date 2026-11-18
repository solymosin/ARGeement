# ARGeement

ARGeement (ARG + *agreement*) maps the entries of AMRFinderPlus, CARD, ResFinder and DisinFinder to
each other by sequence. It does not curate and does not interpret phenotypes: it tells how the
entries of the databases relate to each other **as sequences**. Functional interpretation remains
the task of the source databases.

The latest release is summarised at <https://solymosin.github.io/ARGeement/>; all releases are
archived on Zenodo.

## Installation

```bash
conda env create -f environment.yml
conda activate argeement
tests/run_test.sh          # runs the whole chain on synthetic data
tests/run_pipeline_test.sh # the full Snakemake chain with profiles, from local sources
```

## Usage

```bash
./update.sh                                     # default profile: the latest databases
./update.sh profiles/amrfinder-2026-08-07.yaml  # pinned profile
snakemake --cores 8 --config run_id=test        # build only, without diff

# translating gene names
python -m argeement lookup --release releases/current/latest --from resfinder --to card \
    --names blaTEM-1B,sul1
python -m argeement lookup --release releases/current/latest --from amrfinder --to card \
    --input amrfinder_output.tsv --column "Element symbol"
```

`update.sh` first downloads and parses the databases (cheap) and compares their versions with the
previous release; the expensive build only starts if at least one version has changed (`FORCE=1`
overrides this, e.g. after changing a threshold).

## Pipeline

1. **Download** (`fetch_*`): the CARD broadstreet package, the AMRFinderPlus `latest` directory, the
   ResFinder, DisinFinder and PointFinder git repositories.
2. **Parsing** (`argeement parse`): one standard table. ResFinder/DisinFinder nucleotide sequences
   are translated with genetic code 11, alternative start codons normalised to M. Sequences with
   internal stop codons are excluded from the comparison and listed in `sequence_issues.tsv`. From
   CARD only protein homolog models are used; from AMRFinderPlus only `AMR/AMR` and `STRESS/BIOCIDE`
   entries of core and plus scope, without `POINT`, `POINT_DISRUPT` and `AMR-SUSCEPTIBLE` subtypes and
   without mutation references. An accession may have several catalogue rows; the entry is included
   if any row passes the filters.
3. **Hash layer**: every protein sequence gets a deterministic, version-independent identifier
   (`P-` + first 16 hex digits of SHA-256). Entries with the same `P-` id encode the same protein.
4. **Similarity layer**: MMseqs2 all-against-all search of the unique proteins. Every protein pair
   is assigned once: of the two search directions the alignment with fewer differences (on a tie,
   the higher bit score) is used in both directions.
5. **Relations and clusters**: see below.

## Relation types

| type | condition (default) | meaning |
|---|---|---|
| `identical_protein` | same hash | the same protein |
| `near_identical` | ≥98% identity, both coverages ≥0.95 | probably an allele-level difference |
| `length_variant` | ≥98% identity, one coverage ≥0.95, the other ≥0.80 | same gene, different reference length |
| `same_family` | ≥90% identity, both coverages ≥0.80 | same gene family |
| `homolog` | above the search threshold | more distant relation |

In allele-rich families (beta-lactamases, qnr, mcr, …) `near_identical` means *the closest but a
different allele*: a match at family level, not at allele level. For allele-level agreement use
`identical_protein` or `n_diff == 0`.

The first three are close relations: all are kept, and the diff tracks them. Of the weaker ones,
the best `top_k_weak` per entry and target database are kept. Thresholds are set in `config.yaml`.

## Outputs (`releases/<profile>/<run>/`)

`mappings.tsv` is the main table with directed pairs (both A→B and B→A). Main columns: `relation`,
`n_diff` (number of differences: mismatch or gap positions in the alignment plus the unaligned
ends), `pident`, `cov_a`, `cov_b`, `bits`, `nt_identical` (for identical proteins: whether the
nucleotide sequences are identical too), `best_for_a`, `rbh` (reciprocal best hit), `tie_scope`
(`none`: unambiguous; `within_family`: the tied targets belong to the same family, ambiguous only at
allele level; `cross_family`: different families) and `n_equally_good_b` (number of different target
proteins that are equally good).

`entries.tsv` lists all entries with their protein id (`prot_id`) and family cluster
(`family_cluster`). `summary.tsv` gives the best relation per entry for each pair of databases.
`manifest.json` records database versions, parameters, counts and the MMseqs2 version; `diff/`
holds the changes since the previous release of the profile.

Family clusters (`F-…`) are representative-based: the longest protein is the representative and
only proteins similar to it above the `same_family` threshold join its cluster, so related families
do not chain together. `family_label` is the curated NCBI family: directly for AMRFinderPlus entries
(`ncbi`), and from the closest AMRFinderPlus protein for the other databases (`inferred`, or
`inferred_ambiguous` when the closest protein of another NCBI family is almost as close).

A `homolog` relation in practice means *no counterpart*; reporting the name of a distant homologue
can be misleading.

## NCBI hierarchy

If the AMRFinderPlus database contains `ReferenceGeneHierarchy`, AMRFinderPlus entries get
`ncbi_node`, `ncbi_family` (for an allele leaf its parent, e.g. `blaTEM-1` → `blaTEM`, otherwise the
node itself, e.g. `sul1`) and `lineage`. This enables **family-level lookup** (AMRFinderPlus reports
a family symbol such as `blaTEM` for new alleles; the lookup summarises all entries below the node)
and a **validation against the curated NCBI families** (`validation_ncbi_families.tsv`,
`validation_clusters.tsv`: families split over several clusters and clusters merging families).

## Operons and other reading frames

A ResFinder or DisinFinder entry is translated as a whole only if it has no internal stop codon, a
length divisible by three, and start and stop codons (a truncated reference lacking only one of them
is translated as a whole unless a reverse-strand ORF covers at least 80% of it). All other entries
(e.g. the `VanHAX` operons) are split into ORFs in all six frames (ATG/GTG/TTG start, ≥60 aa; the
selection seeks the best non-overlapping coverage with a cost per ORF and prefers ATG). ORF candidates
related to a curated (AMRFinderPlus or CARD) protein take priority, so that a longer overlapping
foreign ORF cannot displace the real gene. An ORF is kept if it is identical or related to a curated
protein; otherwise ORFs are kept only if the entry is a single gene or a unique operon. Each ORF is
a sub-entry (`<id>#orf<k>`, `parent` column); the split parent does not count as a protein.

## Name concordance

Two directions are analysed: names of identical proteins (`name_concordance.tsv`, eight categories
from identical to different name) and the sequence relation of entries with the same name
(`same_name_concordance.tsv`, with `name_level` telling whether the AMRFinderPlus name is an allele or
a gene name).

## Placing external sequences (`query`)

```bash
python -m argeement query --release releases/current/latest --fasta detected.faa
```

Closest entries per database for protein or nucleotide sequences, with relation type, number of
differences and NCBI family label; ties are all reported (`n_equally_good`).

## Synthetic validation

```bash
python -m argeement validate --release releases/current/latest --n 500 --seed 1
```

Variants with known modifications (1–5 substitutions, N-terminal extension, C-terminal truncation,
CDS on either strand, without start codon, two-gene operon, gene next to a foreign ORF) are made from
a random sample of curated proteins and placed with `query`; per variant type the retrieval of the
original protein, the relation type and the number of differences are measured.

## Real-genome use case (`concordance`)

Agreement of the gene calls of AMRFinderPlus, RGI and ResFinder/DisinFinder on the same genomes,
before mapping (by name) and after mapping (by the relations of the release). The tools run with the
database versions of the release, all with BLAST alignments:

```bash
conda env create -f concordance/env.yaml          # once
conda activate argeement-conc
RUN=work/current/<run> THREADS=36 bash concordance/run_tools.sh
conda activate argeement
python -m argeement concordance --release releases/current/latest \
    --input concordance/calls --out releases/current/latest/concordance
```

The genome list (`concordance/genomes.tsv`) is made with a fixed random seed; the list itself is the
reference for reproduction. Loci are paired within a contig by ≥80% reciprocal overlap.

## Profiles and pinned versions

Several independent release series can coexist. The default profile (`config.yaml`, `current`)
follows the latest databases; a pinned profile reproduces the versions used by an earlier analysis.
A profile file (`profiles/*.yaml`) only contains the differences. Per database: `local` (local
source instead of download), `url`, `ref` (git tag or commit) and `expect_version` (the run stops if
the parsed version does not contain this text). Releases of two profiles can be compared directly:

```bash
python -m argeement diff releases/amrfinder-2026-08-07/latest releases/current/latest --out version_diff
```

## Diff report

The curation interface: `entry_sequence_changed` (different sequence under the same identifier),
`relation_changed` (type of a close relation changed) and `family_changed` (cluster members changed)
are the main candidates for manual review.

## Content report and public export

Every release contains a generated content report (`DATABASES.md`) and a non-redundant public export
(`export/`). The export is laid out like the dump of a relational database: `databases.tsv`,
`relation_types.tsv`, `families.tsv`, `proteins.tsv`, `entries.tsv` and `relations.tsv` (once per
protein pair) refer to each other by integer keys only, and `schema.sql` defines the tables with
their primary and foreign keys, so the export can be loaded into SQLite as it is. The integer keys
are valid within one release; across releases, entries are linked by `db` and `entry_id`, proteins
by `prot_id` and families by `family_cluster`. The export also holds `manifest.json`, `README.md`
(schema and example code) and `SHA256SUMS`. Every row of `mappings.tsv` can be reconstructed from
the export.

## Publication: Zenodo and GitHub Pages

| What | Where | How |
|---|---|---|
| Code | this repository | `git push` |
| Database releases | Zenodo, versioned dataset record | `argeement zenodo` |
| Summary of the latest release | GitHub Pages (`gh-pages` branch) | `argeement site` |

The build runs on a server; publication runs on another machine, the one that holds the git
repository and the credentials. The server has no git repository and no GitHub or Zenodo
credentials: it takes the code from GitHub, and only the package directory (`<release>/package/`)
moves from the server to the publishing machine.

**Build server.** `update_code.sh` takes the code from GitHub (the `main` branch, or a tag,
branch or commit given as argument) and replaces the files of the repository, but never overwrites
`config.yaml` or an existing profile: their new versions are written next to them as `.new`, and
the differences of the settings are shown. Data, releases and local files are left alone.
`server_release.sh` (e.g. weekly from cron) first runs the test suites (`tests/run_test.sh` and
`tests/run_pipeline_test.sh`, log in `logs/`) and stops without building if it fails; otherwise it runs the update, and for a new release
the validation, figures (optional) and `argeement package`.

```bash
./update_code.sh                  # or: ./update_code.sh v1.0.0
tests/run_test.sh
THREADS=36 ./server_release.sh
# 0 3 * * 1  cd /path/to/ARGeement && ARGEEMENT_ENV=/path/to/envs/argeement THREADS=36 ./server_release.sh >> logs/release.log 2>&1
```

**Publishing machine**, from the git repository (code changes are committed and pushed here):

```bash
rsync -a server:/path/to/ARGeement/releases/current/latest/package/ ~/argeement_release/
SANDBOX=1 ./publish_local.sh ~/argeement_release      # test: sandbox.zenodo.org (page not pushed)
./publish_local.sh ~/argeement_release                # production
```

`publish_local.sh` checks the SHA-256 of the archive, publishes a new version of the Zenodo record,
writes the DOI to `zenodo.json` in the package, and pushes the summary page with the new DOI to the
`gh-pages` branch. A release is never published twice. Python with `pyyaml` and `pandas` is enough
on this machine.

One-time setup:

1. *Settings → Pages → Deploy from a branch → `gh-pages` / root* (after the first push of the branch).
2. A Zenodo personal access token with `deposit:write` and `deposit:actions` scopes in
   `~/.config/argeement/zenodo_token` (`chmod 600`); for testing the same on sandbox.zenodo.org in
   `zenodo_sandbox_token`.
3. The first Zenodo version: `NEW=1 ./publish_local.sh <package>` (sandbox: with `SANDBOX=1` too)
   creates the record from the package; the id it prints goes to `deposition_id` (sandbox:
   `sandbox_deposition_id`) in `zenodo.yaml`.
4. Authors, description and links in `zenodo.yaml` and `site.yaml`. The licence notice of the data
   (archive, Zenodo description, web page) comes from `argeement/terms.py`.

## Known limitations

- Mutation-based determinants (gyrA, parC, 23S rRNA, …) are not mapped; they are only inventoried
  in the report.
- Sequence identity does not guarantee functional identity.
- The comparison is at protein level; nucleotide identity is only reported for identical proteins
  (`nt_identical`).

## Licence

The code is under the MIT licence (`LICENSE`). The published data (export archive, Zenodo record)
are under CC-BY 4.0. They contain identifiers, names, public accessions, protein lengths, sequence
hashes and the computed relations, but no sequences of the source databases; of CARD they contain
only ARO accessions and terms, which are under CC-BY 4.0.

The source databases have their own terms, which apply to anyone using them:

| Database | Terms |
|---|---|
| AMRFinderPlus (NCBI) | United States Government Work, public domain |
| CARD (McMaster University) | free for non-commercial research at academic, government and non-profit institutions; commercial use requires a licence from McMaster University; the ARO is under CC-BY 4.0 |
| ResFinder, DisinFinder, PointFinder (CGE, DTU) | Apache License 2.0 |

Running the pipeline downloads CARD; commercial users need a CARD licence for that.

## Citation

If you use ARGeement, please cite:

Vesna Vujić, Georgina Helyes, Adrienn Gréta Tóth, István Csabai, Norbert Solymosi.
ARGeement: sequence-based, versioned mapping of antimicrobial resistance gene databases.
*bioRxiv* 2026.11.17.767227; doi: [https://doi.org/10.64898/2026.11.17.767227](https://doi.org/10.64898/2026.11.17.767227)

Please also cite the release you used: releases are archived on Zenodo, where the concept DOI
[10.5281/zenodo.24278768](https://doi.org/10.5281/zenodo.24278768) always resolves to the latest
release and every release has its own version DOI.


