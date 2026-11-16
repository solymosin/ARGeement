#!/usr/bin/env bash
# Real-genome use case: download genomes and run the three annotation tools with the SAME
# database versions as the release, then compute the agreement.
#
#   conda activate argeement-conc
#   RUN=work/current/2026-11-12_123152 THREADS=36 bash concordance/run_tools.sh
#   conda activate argeement
#   python -m argeement concordance --release releases/current/latest \
#       --input concordance/calls --out releases/current/latest/concordance
#
# RUN: working directory of the run that produced the release; the raw databases are taken from
# there, so the tools use exactly the versions the release was built from.
set -euo pipefail
RUN=${RUN:?give the working directory of the run, e.g. RUN=work/current/2026-11-12_123152}
THREADS=${THREADS:-8}
JOBS=${JOBS:-4}                      # genomes processed in parallel
C=concordance
RAW=$RUN/raw
mkdir -p $C/genomes $C/calls $C/db

# 1. genomes
[[ -s $C/genomes.tsv ]] || python $C/select_genomes.py --out $C/genomes.tsv
if [[ ! -f $C/genomes/.done ]]; then
    tail -n +2 $C/genomes.tsv | cut -f1 > $C/genomes/acc.txt
    datasets download genome accession --inputfile $C/genomes/acc.txt --include genome \
        --filename $C/genomes/genomes.zip
    unzip -o -q $C/genomes/genomes.zip -d $C/genomes
    touch $C/genomes/.done
fi

# 2. databases with the versions of the release
# AMRFinderPlus: the build only downloads the files needed for the mapping (without the HMM
# library), so the complete database is downloaded from NCBI with EXACTLY the version of the
# release (the version directory sits under a directory named after the software version)
if [[ ! -f $C/db/amrfinder/.done ]]; then
    WANT=$(head -1 $RAW/amrfinder/version.txt | tr -d '[:space:]')
    base=https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database
    url=""
    for s in $(curl -fsSL "$base/" | grep -oE 'href="[0-9]+\.[0-9]+/"' | grep -oE '[0-9]+\.[0-9]+' | sort -Vr); do
        if curl -fsSI "$base/$s/$WANT/version.txt" > /dev/null; then url=$base/$s/$WANT; break; fi
    done
    [[ -n $url ]] || { echo "AMRFinderPlus database $WANT not found under $base" >&2; exit 1; }
    echo "AMRFinderPlus database: $url"
    rm -rf $C/db/amrfinder && mkdir -p $C/db/amrfinder
    for f in $(curl -fsSL "$url/" | grep -oE 'href="[^"/?][^"/]*"' | sed 's/^href="//; s/"$//'); do
        curl -fsSL -o "$C/db/amrfinder/$f" "$url/$f"
    done
    got=$(head -1 $C/db/amrfinder/version.txt | tr -d '[:space:]')
    [[ $got == "$WANT" ]] || { echo "Version mismatch: $got (expected: $WANT)" >&2; exit 1; }
    amrfinder_index $C/db/amrfinder
    touch $C/db/amrfinder/.done
fi
if [[ ! -f $C/db/card/.done ]]; then             # RGI: local CARD load with the card.json of the release
    mkdir -p $C/db/card && (cd $C/db/card && rgi load --card_json "$OLDPWD/$RAW/card/card.json" --local)
    touch $C/db/card/.done
fi
# ResFinder/DisinFinder: assembled genomes (-ifa) are searched with BLAST; KMA is not used
# (KMA is meant for raw reads and can miss complete genes in assembled sequences),
# so no KMA index is built
BLASTN=$(command -v blastn) || { echo "blastn not found" >&2; exit 1; }
# ResFinder command: the bioconda package installs run_resfinder.py (recently resfinder);
# "python -m resfinder" only works with the python of the environment (not that of base)
if command -v run_resfinder.py > /dev/null; then RESFINDER_CMD=run_resfinder.py
elif command -v resfinder > /dev/null; then RESFINDER_CMD=resfinder
else RESFINDER_CMD="$(dirname "$(command -v amrfinder)")/python -m resfinder"; fi
$RESFINDER_CMD -h > /dev/null 2>&1 || { echo "ResFinder cannot be started: $RESFINDER_CMD" >&2; exit 1; }
echo "ResFinder: $RESFINDER_CMD"
for g in resfinder disinfinder; do
    if [[ ! -f $C/db/$g/.done ]]; then
        rm -rf $C/db/$g && cp -r $RAW/$g $C/db/$g
        touch $C/db/$g/.done
    fi
done
cat $C/db/amrfinder/version.txt > $C/db/versions.txt 2>/dev/null || true
{ amrfinder --version; rgi main --version; $RESFINDER_CMD --version 2>&1 | head -1; blastn -version | head -1; } \
    >> $C/db/versions.txt 2>&1 || true

# 3. tools per genome
one() {
    set +e                      # an error of one tool does not stop the other genomes
    acc=$1
    fna=$(ls $C/genomes/ncbi_dataset/data/$acc/*.fna | head -1)
    o=$C/calls/$acc
    mkdir -p $o
    [[ -s $o/amrfinder.tsv ]] || amrfinder -n "$fna" -d $C/db/amrfinder --plus \
        --threads $THREADS -o $o/amrfinder.tsv
    if [[ ! -s $o/rgi.txt ]]; then
        # RGI: BLAST alignment (instead of the default DIAMOND), as for the other two tools
        # In contig mode RGI recovers the contig name from the ORF id, which fails for RefSeq names
        # ("NZ_…") and for simplified names too ("Requested rname … does not exist"). ORFs are
        # therefore predicted here with Prodigal, as RGI does in contig mode (single mode), and RGI
        # is run in protein mode. ORF locations are kept in rgi_orfs.tsv
        # (name, contig, start, end, strand).
        prodigal -q -i "$fna" -a "$o/rgi_prodigal.faa" -o /dev/null
        awk -v map="$o/rgi_orfs.tsv" '
            /^>/ {n++; split(substr($0, 2), f, " # "); ctg = f[1]; sub(/_[0-9]+$/, "", ctg);
                  nm = sprintf("orf%06d", n); print nm "\t" ctg "\t" f[2] "\t" f[3] "\t" f[4] > map;
                  print ">" nm; next}
            {gsub(/\*$/, ""); print}' "$o/rgi_prodigal.faa" > "$o/rgi_input.faa"
        (cd $C/db/card && rgi main -i "$OLDPWD/$o/rgi_input.faa" -o "$OLDPWD/$o/rgi" -t protein --local \
            -a BLAST --clean -n $THREADS) || echo "!! RGI error ($acc)" >&2
        [[ -s $o/rgi.txt ]] || echo "!! RGI: no rgi.txt ($acc)" >&2
        rm -f "$o/rgi_input.faa" "$o/rgi_prodigal.faa"
    fi
    if [[ ! -s $o/resfinder.json ]]; then
        $RESFINDER_CMD -ifa "$fna" -o $o/resfinder -b "$BLASTN" -acq -db_res $C/db/resfinder \
            -d -db_disinf $C/db/disinfinder -l 0.6 -t 0.8 -j $o/resfinder/resfinder.json \
            > $o/resfinder.log 2>&1
        js=$(ls $o/resfinder/*.json 2>/dev/null | head -1 || true)
        if [[ -n $js ]]; then
            cp "$js" $o/resfinder.json
        else
            echo "!! ResFinder: no JSON output ($acc), see $o/resfinder.log" >&2
        fi
    fi
}
export -f one
export C THREADS BLASTN RESFINDER_CMD
# the first genome runs alone: RGI (and possibly other tools) builds its BLAST database on the
# first run, and parallel jobs would write it at the same time (corrupted database,
# "Database memory map file error")
first=$(tail -n +2 $C/genomes.tsv | head -1 | cut -f1)
one "$first"
tail -n +2 $C/genomes.tsv | cut -f1 | xargs -P "$JOBS" -I{} bash -c 'one {}'
echo "Done: $C/calls (amrfinder.tsv, rgi.txt, resfinder.json per genome)"
