"""Sequence utilities: translation, cleaning, hashing."""
import hashlib

from Bio.Seq import Seq

# valid start codons of the bacterial genetic code (table 11)
VALID_STARTS = {"ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"}

# these issues exclude a sequence from the protein-level comparison
BLOCKING_ISSUES = {"internal_stop", "empty", "split_into_orfs", "orf_unsupported"}


def clean_protein(seq: str) -> str:
    s = "".join(seq.split()).upper()
    return s[:-1] if s.endswith("*") else s


def clean_nucleotide(seq: str) -> str:
    return "".join(seq.split()).upper().replace("U", "T")


def translate_cds(nt: str):
    """Translates a CDS with genetic code 11. Returns (protein, list of issues).

    An alternative start codon is normalised to M, so that the translated sequence is
    comparable with the sequences stored as proteins by the databases.
    """
    nt = clean_nucleotide(nt)
    if not nt:
        return "", ["empty"]
    issues = []
    nt_use = nt
    if len(nt) % 3:
        issues.append("len_not_multiple_of_3")
        nt_use = nt[: len(nt) - len(nt) % 3]
    if set(nt_use) - set("ACGT"):
        issues.append("ambiguous_bases")
    aa = str(Seq(nt_use).translate(table=11))
    if nt_use[:3] in VALID_STARTS:
        aa = "M" + aa[1:]
    else:
        issues.append("no_start_codon")
    if aa.endswith("*"):
        aa = aa[:-1]
    else:
        issues.append("no_stop_codon")
    if "*" in aa:
        issues.append("internal_stop")
    return aa, issues


def seq_hash(seq: str, prefix: str) -> str:
    """Deterministic, version-independent identifier derived from the sequence."""
    return f"{prefix}-{hashlib.sha256(seq.encode()).hexdigest()[:16]}"


ORF_STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}


def _revcomp(nt):
    return nt.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def find_orfs(nt, min_aa=60):
    """Open reading frames in all six frames (start: ATG/GTG/TTG, up to the stop codon).

    Returns a list of non-overlapping ORFs in sequence order, each a
    dict(start, end, strand, frame, protein), where start/end are 0-based, half-open
    coordinates on the original (+) strand.
    """
    return _best_tiling(orf_candidates(nt, min_aa))


def orf_candidates(nt, min_aa=60):
    """All ORF candidates in all six frames (with every possible start codon), without
    selection. 'stop_key' identifies the stop codon of a candidate; candidates sharing a stop
    codon are versions of the same gene with different start codons."""
    nt = clean_nucleotide(nt)
    n = len(nt)
    cands = []
    for strand, seq in (("+", nt), ("-", _revcomp(nt))):
        for frame in range(3):
            starts = []          # all start codons seen since the previous stop
            for i in range(frame, n - 2, 3):
                cod = seq[i:i + 3]
                if cod in STOPS:
                    # every possible start is a candidate: the selection decides which one
                    # does not run into the neighbouring gene
                    cands += [(strand, frame, st, i + 3) for st in starts
                              if (i - st) // 3 >= min_aa]
                    starts = []
                elif cod in ORF_STARTS:
                    starts.append(i)
            # ORF running to the end of the sequence without a stop (truncated reference)
            end = n - (n - frame) % 3
            cands += [(strand, frame, st, end) for st in starts if (end - st) // 3 >= min_aa]
    orfs = []
    for strand, frame, s, e in cands:
        seq = nt if strand == "+" else _revcomp(nt)
        aa, _ = translate_cds(seq[s:e])
        # coordinates on the (+) strand
        ps, pe = (s, e) if strand == "+" else (n - e, n - s)
        orfs.append(dict(start=ps, end=pe, strand=strand, frame=frame, protein=aa,
                         atg=seq[s:s + 3] == "ATG",
                         stop_key=f"{strand}{pe if strand == '+' else ps}"))
    return orfs


def _best_tiling(orfs, max_overlap=8, orf_cost=150):
    """Weighted interval scheduling: the set of non-overlapping ORFs that together cover the
    most nucleotides. The small overlaps usual in operons (e.g. ATGA) are allowed.
    More robust than 'longest wins': a fortuitously long ORF in another frame cannot displace
    several real genes behind it."""
    if not orfs:
        return []
    o = sorted(orfs, key=lambda x: x["end"])
    ends = [x["end"] for x in o]
    import bisect
    best = [0] * (len(o) + 1)
    take = [False] * len(o)
    prev = []
    for i, x in enumerate(o):
        j = bisect.bisect_right(ends, x["start"] + max_overlap, 0, i)   # compatible predecessors
        prev.append(j)
        # every ORF carries a fixed cost (splitting must not pay off: one complete gene beats two
        # shorter pieces covering a few more nucleotides), and the rarer start codons
        # (GTG, TTG) a 9.5-nucleotide penalty: exact ties are impossible, and if a non-ATG start
        # is at most 3 codons upstream, the ATG wins
        # 'bonus' gives priority to candidates related to a curated protein (build.rescue_split_orfs)
        w = x["end"] - x["start"] - orf_cost - (0 if x.get("atg", True) else 9.5) + x.get("bonus", 0)
        if best[j] + w > best[i]:
            best[i + 1], take[i] = best[j] + w, True
        else:
            best[i + 1] = best[i]
    chosen, i = [], len(o)
    while i > 0:
        if take[i - 1]:
            chosen.append(o[i - 1])
            i = prev[i - 1]
        else:
            i -= 1
    return sorted(chosen, key=lambda x: x["start"])


def choose_translation(nt, min_orf_aa=60, rev_cover=0.8):
    """Whether an entry can be translated as a whole or must be split into ORFs.

    Returns (protein, issues, orfs). If orfs is not empty, the entry is split.
    - Clean (no internal stop, length divisible by 3) with start and stop codon: as a whole.
    - Clean but lacking the start or stop codon (truncated reference): as a whole, unless an ORF
      on the reverse strand covers at least rev_cover of the sequence (gene stored reversed).
    - Everything else: ORF splitting on both strands.
    """
    aa, issues = translate_cds(nt)
    clean = not ({"internal_stop", "len_not_multiple_of_3", "empty"} & set(issues))
    if clean and not ({"no_start_codon", "no_stop_codon"} & set(issues)):
        return aa, issues, []
    orfs = find_orfs(nt, min_orf_aa)
    if clean:
        n = len(clean_nucleotide(nt))
        rev = [o for o in orfs if o["strand"] == "-" and (o["end"] - o["start"]) >= rev_cover * n]
        if not rev:
            return aa, issues, []
    return aa, issues, orfs


def n_diff_from_alignment(alnlen, nident, qstart, qend, qlen, tstart, tend, tlen):
    """Number of differences from an alignment: mismatched and gap positions within the alignment
    plus the unaligned ends. At each end the longer overhang of the two proteins counts, so if the
    alignment ends at the same place on both proteins (e.g. a substitution at the last position),
    that is one difference, not two. Works element-wise with numpy/pandas arrays as well."""
    import numpy as np
    left = np.maximum(qstart - 1, tstart - 1)
    right = np.maximum(qlen - qend, tlen - tend)
    return (alnlen - nident) + left + right
