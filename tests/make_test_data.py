"""Synthetic test data in the formats of the databases (no real sequences).

Usage: python tests/make_test_data.py <output directory> [--v2]
--v2 writes an "updated" version for testing the diff.
"""
import json
import random
import sys
from pathlib import Path

rng = random.Random(42)
AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"], "C": ["TGT", "TGC"], "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"], "F": ["TTT", "TTC"], "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"], "K": ["AAA", "AAG"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "M": ["ATG"], "N": ["AAT", "AAC"],
    "P": ["CCT", "CCC", "CCA", "CCG"], "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"], "V": ["GTT", "GTC", "GTA", "GTG"], "W": ["TGG"],
    "Y": ["TAT", "TAC"]}


def prot(n):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def mutate(p, k, seed):
    r = random.Random(seed)
    s = list(p)
    for pos in r.sample(range(1, len(p)), k):
        s[pos] = r.choice([a for a in AA if a != s[pos]])
    return "".join(s)


def back(p, seed=0, start="ATG"):
    r = random.Random(seed)
    return start + "".join(r.choice(CODONS[a]) for a in p[1:]) + "TAA"


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def wrap(s):
    return "\n".join(s[i:i + 60] for i in range(0, len(s), 60))


v2 = "--v2" in sys.argv
out = Path(sys.argv[1])

TEM1 = prot(286)
TEM116 = mutate(TEM1, 1, 1)
TEM_X = mutate(TEM1, 14, 2)            # ~95%: same family
SUL1 = prot(279)
SUL1_V = mutate(SUL1, 1, 3)
OXA = prot(265)
OXA_EXT = OXA + "".join(rng.choice(AA) for _ in range(25))   # hosszabb referencia
QNR = prot(218)
QNR_R = mutate(QNR, 2, 6)              # ResFinder
QNR_A = mutate(QNR_R, 1, 4)            # two CARD entries, both
QNR_B = mutate(QNR_R, 1, 5)            # exactly 1 difference away: a genuine tie
QNR_C = mutate(QNR_R, 3, 10)           # 3 differences: must NOT be tied
GTGP = prot(180)                       # GTG starttal a ResFinderben
ONLY_CARD = prot(300)
GYRA = prot(875)
QAC = prot(110)
VIM = prot(266)
VANH, VANA, VANX = prot(322), prot(343), prot(202)   # operon a ResFinderben
REV = prot(280)                                       # gene stored on the reverse strand in ResFinder


def _clean_reverse(n, seed0):
    """A protein and coding sequence whose reverse complement has NO stop codon in forward
    frame 0: the pattern of reverse-strand genes that could be mistranslated as a whole."""
    stops = {"TAA", "TAG", "TGA"}
    for sd in range(seed0, seed0 + 100000):
        r_ = random.Random(sd)
        p_ = "M" + "".join(r_.choice(AA) for _ in range(n - 1))
        nt_ = "ATG" + "".join(r_.choice(CODONS[a_]) for a_ in p_[1:]) + "TAA"
        rc = nt_.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if not any(rc[i:i + 3] in stops for i in range(0, len(rc) - 2, 3)):
            return p_, nt_
    raise RuntimeError("no suitable sequence found")


REVC, REVC_NT = _clean_reverse(70, 500)   # reverse strand, no stop in the forward frame
UNIQ = prot(200)                          # ResFinder only, reverse strand, no curated counterpart
EXTRA = prot(90)                          # neighbouring non-resistance ORF in DisinFinder
# chaining: CH_A ~ CH_B and CH_B ~ CH_C (~93%), but CH_A ~ CH_C only ~86%
CH_A = prot(400)
CH_B = mutate(CH_A, 28, 41)
_cb = list(CH_B)
for _pos in random.Random(42).sample([i for i in range(1, 400) if CH_A[i] == CH_B[i]], 28):
    _cb[_pos] = random.Random(_pos).choice([x for x in AA if x != _cb[_pos]])
CH_C = "".join(_cb)
TEM_NEW = mutate(TEM1, 2, 7)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from argeement.sequtil import _best_tiling, choose_translation, find_orfs, orf_candidates  # noqa: E402


def _displaced(n, seed0):
    """A real gene (n amino acids) completely overlapped by a longer foreign ORF on the reverse
    strand, so that selection by coverage alone displaces it. Synonymous codons of the gene are
    chosen so that one frame of the reverse strand has no stop codon."""
    rstop = {"TTA", "CTA", "TCA"}          # read forwards: stop codons on the reverse strand
    for sd in range(seed0, seed0 + 20000):
        r_ = random.Random(sd)
        p = "M" + "".join(r_.choice(AA) for _ in range(n - 1))

        def free(k):                       # k nucleotides without a reverse stop in the triplets at offset 1
            while True:
                x = "".join(r_.choice("ACGT") for _ in range(k))
                if not any(x[i:i + 3] in rstop for i in range(1, k - 2, 3)):
                    return x
        aas = p[1:] + "*"
        opts = [CODONS[a_] if a_ != "*" else ["TAA", "TAG", "TGA"] for a_ in aas]
        # from the end: which codons can continue without a reverse stop
        feas = [None] * len(aas)
        feas[-1] = set(opts[-1])
        for k_ in range(len(aas) - 2, -1, -1):
            feas[k_] = {c_ for c_ in opts[k_] if any((c_[1:] + d_[0]) not in rstop for d_ in feas[k_ + 1])}
        g = "ATG"
        for k_ in range(len(aas)):
            ok = [c_ for c_ in sorted(feas[k_]) if (g[-2:] + c_[0]) not in rstop]
            if not ok:
                break
            g += r_.choice(ok)
        else:
            # left edge: random + reverse stop (TTA at offset 1) + stop in the gene frame
            left = "".join(r_.choice("ACGT") for _ in range(150)) + "T" + "TTA" + "AG" + free(57) + "TAA"
            right = free(120) + "C" + "CAT" + "TT" + "ATT" + "AC" + \
                "".join(r_.choice("ACGT") for _ in range(150))
            s_ = left + g + right
            if (g[-2:] + right[0]) in rstop:
                continue
            cs = orf_candidates(s_, 60)
            grp = {c_["stop_key"] for c_ in cs if c_["protein"] == p}
            if not grp or any(o_["protein"] == p for o_ in find_orfs(s_, 60)) \
                    or not choose_translation(s_, 60)[2]:
                continue
            for c_ in cs:
                c_["bonus"] = 100000 if c_["stop_key"] in grp else 0
            if any(o_["protein"] == p for o_ in _best_tiling(cs)):
                return p, s_
    raise RuntimeError("no suitable sequence found")


QACZ, QACZ_NT = _displaced(90, 700)   # below 100 aa: not sampled by the validation

# ---------------- CARD
c = out / "card"
c.mkdir(parents=True, exist_ok=True)
card = [("ARO:3000873", "AAP20891.1", "TEM-1", TEM1, "TEM beta-lactamase"),
        ("ARO:3000979", "AAB39956.1", "TEM-116", TEM116, "TEM beta-lactamase"),
        ("ARO:3000410", "CAA35633.1", "sul1", SUL1, "sulfonamide resistant sul"),
        ("ARO:3001396", "AAA93527.1", "OXA-1", OXA, "OXA beta-lactamase"),
        ("ARO:3002711", "ABC11111.1", "QnrA1", QNR_A, "quinolone resistance protein (qnr)"),
        ("ARO:3002712", "ABC22222.1", "QnrA2", QNR_B, "quinolone resistance protein (qnr)"),
        ("ARO:3002713", "ABC33333.1", "QnrA3", QNR_C, "quinolone resistance protein (qnr)"),
        ("ARO:3004000", "ABD33333.1", "Ecol_catX", GTGP, "chloramphenicol acetyltransferase (CAT)"),
        ("ARO:3005000", "ABE44444.1", "onlyInCard", ONLY_CARD, "some family"),
        ("ARO:3005001", "ABF55555.1", "qacA", QAC, "MFS antibiotic efflux pump"),
        ("ARO:3005002", "ABF55556.1", "qacZ", QACZ, "SMR antibiotic efflux pump"),
        ("ARO:3002274", "CAB66666.1", "VIM-4", VIM, "VIM beta-lactamase"),
        ("ARO:3002942", "AAA65956.1", "vanH", VANH, "glycopeptide resistance gene cluster"),
        ("ARO:3003999", "ABG77777.1", "REV-1", REV, "REV beta-lactamase"),
        ("ARO:3003998", "ABG88888.1", "aacREVC", REVC, "AAC(3)"),
        ("ARO:3002275", "CAB88888.1", "VIM-2", mutate(VIM, 5, 51), "VIM beta-lactamase"),
        ("ARO:3003900", "ABH99999.1", "MexX", prot(380), "RND efflux pump")]
if v2:
    card[1] = ("ARO:3000979", "AAB39956.1", "TEM-116 beta-lactamase", TEM116, "TEM beta-lactamase")
cj = {"_version": "4.0.2" if v2 else "4.0.1"}
CARD_CLS = {"TEM-1": ["penam", "cephalosporin", "azole antibiotic", "imidazole antibiotic"], "sul1": ["sulfonamide antibiotic"],
            "OXA-1": ["penam"], "Ecol_catX": ["phenicol antibiotic"], "VIM-4": ["carbapenem"], "vanH": ["glycopeptide antibiotic"], "aacREVC": ["aminoglycoside antibiotic"],
            "REV-1": ["cephalosporin"],
            "onlyInCard": ["glycopeptide antibiotic"], "qacA": ["disinfecting agents and antiseptics"]}
CARD_MECH = {"sul1": "antibiotic target replacement", "onlyInCard": "antibiotic target alteration"}
with open(c / "protein_fasta_protein_homolog_model.fasta", "w") as fp, \
        open(c / "nucleotide_fasta_protein_homolog_model.fasta", "w") as fn:
    for i, (aro, acc, name, p, fam) in enumerate(card):
        fp.write(f">gb|{acc}|{aro}|{name} [Escherichia coli]\n{wrap(p)}\n")
        fn.write(f">gb|X{i:05d}.1|+|1-{3 * len(p) + 3}|{aro}|{name} [Escherichia coli]\n"
                 f"{wrap(back(p, 100 + i))}\n")
        cats = {"1": {"category_aro_class_name": "AMR Gene Family", "category_aro_name": fam},
                "2": {"category_aro_class_name": "Resistance Mechanism",
                      "category_aro_name": CARD_MECH.get(name, "antibiotic inactivation")}}
        key = next((k for k in CARD_CLS if name.startswith(k)), None)
        for j, dc in enumerate(CARD_CLS.get(key, ["fluoroquinolone antibiotic"])):
            cats[f"d{j}"] = {"category_aro_class_name": "Drug Class", "category_aro_name": dc}
        cj[str(i)] = {"ARO_accession": aro.split(":")[1], "model_name": name,
                      "model_type": "protein homolog model", "ARO_category": cats}
    # a variant model not in the homolog FASTA: appears in the inventory only
    cj["900"] = {"ARO_accession": "3003926", "model_name": "Escherichia coli gyrA",
                 "model_type": "protein variant model", "ARO_category": {},
                 "model_sequences": {"sequence": {"1": {"NCBI_taxonomy": {"NCBI_taxonomy_name": "Escherichia coli"}}}},
                 "model_param": {"snp": {"param_value": {"1": "S83L", "2": "D87N"}}}}
(c / "card.json").write_text(json.dumps(cj))

# ---------------- AMRFinderPlus
a = out / "amrfinder"
a.mkdir(parents=True, exist_ok=True)
amr = [("WP_000027057.1", "blaTEM-1", "blaTEM", "AMR", "AMR", TEM1),
       ("WP_015058868.1", "blaTEM-116", "blaTEM", "AMR", "AMR", TEM116),
       ("WP_063840001.1", "blaTEM-999", "blaTEM", "AMR", "AMR", TEM_X),
       ("WP_000259031.1", "sul1", "sul1", "AMR", "AMR", SUL1),
       ("WP_001334766.1", "blaOXA-1", "blaOXA-1", "AMR", "AMR", OXA),
       ("WP_000072067.1", "gyrA", "gyrA", "AMR", "POINT", GYRA),
       ("WP_063840002.1", "blaTEM-far", "blaTEM", "AMR", "AMR", mutate(TEM1, 50, 12)),  # NCBI family split
       ("WP_004193055.1", "amrR", "amrR", "AMR", "AMR-SUSCEPTIBLE", prot(210)),
       ("WP_000000001.1", "qacEdelta1", "qacE", "STRESS", "BIOCIDE", QAC),      # bevonva
       ("WP_000000002.1", "merA", "merA", "STRESS", "METAL", prot(560)),       # excluded
       ("WP_000000003.1", "blaVIM-2", "blaVIM", "AMR", "AMR", VIM),
       ("WP_000000004.1", "vanA", "vanA", "AMR", "AMR", VANA),
       ("WP_000000009.1", "vanX-A", "vanX", "AMR", "AMR", VANX),
       ("WP_000000005.1", "blaCHA-1", "blaCHA", "AMR", "AMR", CH_A),
       ("WP_000000006.1", "blaCHB-1", "blaCHB", "AMR", "AMR", CH_B),
       ("WP_000000007.1", "blaCHC-1", "blaCHC", "AMR", "AMR", CH_C),
       ("WP_000000008.1", "mexX", "mexX", "AMR", "AMR", prot(380))]
if v2:
    amr.append(("WP_099999999.1", "blaTEM-300", "blaTEM", "AMR", "AMR", TEM_NEW))
    amr[2] = ("WP_063840001.1", "blaTEM-999", "blaTEM", "AMR", "AMR", mutate(TEM_X, 1, 9))
cat_cols = ["allele", "gene_family", "whitelisted_taxa", "product_name", "scope", "type",
            "subtype", "class", "subclass", "refseq_protein_accession",
            "refseq_nucleotide_accession", "genbank_protein_accession", "hierarchy_node"]
NODE = {"gyrA": ""}   # POINT rows have no node
# .tsv extension in v2: both file name variants are tested
cat_name = "ReferenceGeneCatalog.tsv" if v2 else "ReferenceGeneCatalog.txt"
with open(a / cat_name, "w") as fc, open(a / "AMRProt.fa", "w") as fp, \
        open(a / "AMR_CDS.fa", "w") as fn:
    fc.write("\t".join(cat_cols) + "\n")
    for i, (acc, allele, fam, t, st, p) in enumerate(amr):
        klass = {"blaTEM": "BETA-LACTAM", "blaOXA-1": "BETA-LACTAM", "sul1": "SULFONAMIDE",
                 "gyrA": "QUINOLONE", "qacE": "QUATERNARY AMMONIUM", "blaVIM": "BETA-LACTAM", "vanA": "GLYCOPEPTIDE", "vanX": "GLYCOPEPTIDE",
                 "merA": "MERCURY"}.get(fam, "AMINOGLYCOSIDE")
        fc.write("\t".join([allele, fam, "", f"{allele} product", "core", t, st, klass, klass,
                            acc, f"NG_0{i:05d}.1", "", NODE.get(allele, allele)]) + "\n")
        fp.write(f">0|{acc}|1|1|{allele}|{fam}|hydrolase|2|X|X|{allele} product\n{wrap(p)}\n")
        fn.write(f">{i}|NG_0{i:05d}.1|1|1|{allele}|{fam}|X|{acc}|X\n{wrap(back(p, 100 + i))}\n")
    # fexA-like case: the same accession has a POINT row and then an AMR row
    fexa = prot(475)
    fc.write("\t".join(["fexA_A37V", "fexA", "", "fexA mutation", "core", "AMR", "POINT",
                        "X", "X", "WP_015585966.1", "", "", ""]) + "\n")
    fc.write("\t".join(["fexA", "fexA", "", "fexA efflux", "core", "AMR", "AMR",
                        "PHENICOL/OXAZOLIDINONE", "X", "WP_015585966.1", "", "", "fexA"]) + "\n")
    # fexB: another NCBI family, but close to fexA in sequence -> same cluster
    fc.write("\t".join(["fexB", "fexB", "", "fexB efflux", "core", "AMR", "AMR",
                        "X", "X", "WP_015585967.1", "", "", "fexB"]) + "\n")
    fp.write(f">0|WP_015585967.1|1|1|fexB|fexB|efflux|2|X|X|fexB\n{wrap(mutate(fexa, 15, 13))}\n")
    fp.write(f">0|WP_015585966.1|1|1|fexA|fexA|efflux|2|X|X|fexA\n{wrap(fexa)}\n")
    fp.write(f">WP_015585966.1|1|1|fexA|fexA|mutation|2|||fexA_mut\n{wrap(mutate(fexa, 1, 11))}\n")
    # POINT_DISRUPT with a normal header: must be filtered via the catalogue
    fc.write("\t".join(["mgrB", "mgrB", "", "MgrB", "core", "AMR", "POINT_DISRUPT",
                        "X", "X", "WP_002911375.1", "", "", ""]) + "\n")
    fp.write(f">0|WP_002911375.1|1|1|mgrB|mgrB|x|2|X|X|mgrB\n{wrap(prot(47))}\n")
    # mutation reference in the real AMRProt format: filter, do not warn
    fp.write(f">WP_000489247.1:90:cirA_S90YfsTer15ins2|1|1|cirA|cirA|mutation|2|||"
             f"catecholate_siderophore_receptor_CirA\n{wrap(prot(120))}\n")
    # a protein not in the catalogue
    fp.write(f">0|WP_123456789.1|1|1|mystery|mystery|x\n{wrap(prot(100))}\n")
(a / "AMRProt-mutation.tsv").write_text(
    "#taxgroup\taccession_version\tmutation_position\tstandard_mutation_symbol\treported_mutation_symbol\tclass\tsubclass\tmutated_protein_name\n"
    "Escherichia\tWP_000072067.1\t83\tgyrA_S83L\tgyrA_S83L\tQUINOLONE\tQUINOLONE\tx\n"
    "Escherichia\tWP_000072067.1\t87\tgyrA_D87N\tgyrA_D87N\tQUINOLONE\tQUINOLONE\tx\n"
    "Salmonella\tWP_000072068.1\t57\tparC_S57I\tparC_S57I\tQUINOLONE\tQUINOLONE\tx\n")
(a / "AMR_DNA-Escherichia.tsv").write_text(
    "#accesion_version@gene_name@gene_symbol:start-stop\tmutation_position\tstandard_mutation_symbol\t"
    "reported_mutation_symbol\tclass\tsubclass\tmutated_gene_name\n"
    "NC_000913.3@16S_ribosomal_RNA@16S_rrsB:4166659-4168201\t523\t16S_A523C\t16S_A523C\t"
    "AMINOGLYCOSIDE\tSTREPTOMYCIN\tEscherichia_streptomycin_resistant_16S_rrsB\n")
(a / "version.txt").write_text("2026-08-15.1\n" if not v2 else "2026-10-20.1\n")
# NCBI hierarchy (with the columns of the real file) and family rules
hier = [("ALL", "", "", "0"), ("AMR", "ALL", "", "0"), ("bla-A", "AMR", "", "0"),
        ("blaTEM", "bla-A", "blaTEM", "0"), ("bla-D", "AMR", "", "0"),
        ("blaOXA", "bla-D", "blaOXA", "0"), ("blaOXA-1", "blaOXA", "blaOXA-1", "1"),
        ("sul", "AMR", "sul", "0"), ("sul1", "sul", "sul1", "0"),
        ("fexA", "AMR", "fexA", "0"), ("fexB", "AMR", "fexB", "0"), ("amrR", "AMR", "amrR", "0"),
        ("qacEdelta1", "AMR", "qacEdelta1", "0"), ("blaVIM", "bla-B", "blaVIM", "0"),
        ("bla-B", "AMR", "", "0"), ("blaVIM-2", "blaVIM", "blaVIM-2", "1")]
SYN = {"qacEdelta1": "qacA,qacEdelta"}
hier += [(t, "blaTEM", t, "1") for t in ("blaTEM-1", "blaTEM-116", "blaTEM-999", "blaTEM-300", "blaTEM-far")]
hcols = ["node_id", "parent_node_id", "symbol", "num_prots", "hmm_acc", "prot_acc", "scope", "type",
         "subtype", "class", "subclass", "name", "synonyms", "display_parent", "allele", "db_version"]
with open(a / "ReferenceGeneHierarchy.txt", "w") as fh:
    fh.write("\t".join(hcols) + "\n")
    for node, parent, sym, allele in hier:
        fh.write("\t".join([node, parent, sym, "1", "", "", "core", "AMR", "AMR", "", "", node,
                             SYN.get(node, ""), "0", allele, "x"]) + "\n")
with open(a / "fam.tsv", "w") as fh:
    fh.write("#node_id\tparent_node_id\tgene_symbol\tblastrule_complete_ident\n")
    fh.write("blaTEM\tbla-A\tblaTEM\t0.00\nsul1\tsul\tsul1\t94.00\nfexA\tAMR\tfexA\t90.00\n")

# ---------------- ResFinder
r = out / "resfinder"
r.mkdir(parents=True, exist_ok=True)
with open(r / "beta-lactam.fsa", "w") as f:
    f.write(f">blaTEM-1B_1_AY458016\n{wrap(back(TEM1, 999))}\n")       # synonymous difference
    f.write(f">blaOXA-1_1_HQ170510\n{wrap(back(OXA_EXT, 5))}\n")
with open(r / "sulphonamide.fsa", "w") as f:
    f.write(f">sul1_2_U12338\n{wrap(back(SUL1, 102))}\n")              # = CARD nt
    f.write(f">sul1_3_EU855787\n{wrap(back(SUL1_V, 7))}\n")
with open(r / "quinolone.fsa", "w") as f:
    f.write(f">qnrA1_1_AY070235\n{wrap(back(QNR_R, 8))}\n")
    bad = back(prot(150), 9)
    f.write(f">qnrBroken_1_XX000001\n{wrap(bad[:60] + 'TAA' + bad[63:])}\n")  # internal stop
with open(r / "glycopeptide.fsa", "w") as f:
    op = back(VANH, 31) + "AGGAGT" + back(VANA, 32) + "TT" + back(VANX, 33)
    f.write(f">VanHAX_1_TEST01\n{wrap(op)}\n")
with open(r / "beta-lactam.fsa", "a") as f:
    f.write(f">blaREV-1_1_X01\n{wrap(revcomp(back(REV, 34)))}\n")
with open(r / "aminoglycoside.fsa", "w") as f:
    f.write(f">aacREVC_1_X02\n{wrap(revcomp(REVC_NT))}\n")         # clean forward frame!
    f.write(f">uniqREV_1_X04\n{wrap(revcomp(back(UNIQ, 35)))}\n")
with open(r / "phenicol.fsa", "w") as f:
    f.write(f">catX_1_AB000001\n{wrap(back(GTGP, 106, start='GTG'))}\n")
if v2:   # phenotypes.txt in v2: with an ambiguous and an unknown class term
    (r / "phenotypes.txt").write_text(
        "Gene_accession no.\tClass\tPhenotype\tPMID\tMechanism of resistance\n"
        "blaTEM-1B_1_AY458016\tBeta-lactam\tAmpicillin\t-\tEnzymatic inactivation\n"
        "sul1_2_U12338\tFolate pathway antagonist\tSulfamethoxazole\t-\tReplacement\n"
        "qnrA1_1_AY070235\tQuinolone, Iminophenazine\tCiprofloxacin\t-\tTarget protection\n")
(r / "VERSION.argeement").write_text("2.4.0 (git abc1234)\n" if not v2 else "2.5.0 (git def5678)\n")
# ---------------- DisinFinder
dd = out / "disinfinder"
dd.mkdir(parents=True, exist_ok=True)
with open(dd / "disinfectants.fsa", "w") as f:
    f.write(f">qacA_1_AB566410\n{wrap(back(QAC, 21))}\n")
    f.write(f">formA_1_X73835\n{wrap(back(prot(150), 22))}\n")
    # real gene + a neighbouring non-resistance ORF
    f.write(f">qacX_1_X03\n{wrap(back(QAC, 23) + 'AGGAGT' + back(EXTRA, 24))}\n")
    # unique operon without curated counterpart (like sitABCD): all its genes must be kept
    f.write(f">sitXY_1_X05\n{wrap(back(prot(300), 25) + back(prot(280), 26))}\n")
    # the real gene is displaced by a longer overlapping foreign ORF
    f.write(f">qacZ_1_X06\n{wrap(QACZ_NT)}\n")
(dd / "phenotypes.txt").write_text(
    "Gene_accession no.\tClass\tPhenotype\tPMID\tMechanism of resistance\tNotes\tRequired_gene\n"
    "qacA_1_AB566410\tQuaternary Ammonium Compounds\tBenzylkonium Chloride\t1\tEfflux pump\n"
    "formA_1_X73835\tAldehydes\tFormaldehyde\t1\tEnzymatic degradation\n")
(dd / "VERSION.argeement").write_text("2.0.1 (git 09fc831)\n")

# ---------------- PointFinder (inventory only)
pfd = out / "pointfinder"
for sp, rows in {"campylobacter": [("gyrA", "86", "ACA", "T", "I,A,K,V", "Nalidixic acid, Ciprofloxacin", ""),
                                   ("gyrA", "85", "GAT", "D", "Y", "Nalidixic acid, Ciprofloxacin", "gyrA_T86I"),
                                   ("23S", "2075", "A", "A", "G", "Erythromycin", "")],
                 "salmonella": [("parC", "57", "TCG", "S", "I", "Ciprofloxacin", "")]}.items():
    (pfd / sp).mkdir(parents=True, exist_ok=True)
    with open(pfd / sp / "resistens-overview.txt", "w") as fh:
        fh.write("#Gene_ID\tGene_name\tCodon_pos\tRef_nuc\tRef_codon\tRes_codon\tResistance\tPMID\t"
                 "Mechanism\tNotes\tRequired_mut\n")
        for g, pos, rn, rc, res, ab, req in rows:
            fh.write(f"{g}\t{g}\t{pos}\t{rn}\t{rc}\t{res}\t{ab}\t1\tTarget modification\t\t{req}\n")
    genes = sorted({r[0] for r in rows if r[0] != "23S"})
    (pfd / sp / "genes.txt").write_text("\n".join(genes) + "\n")
    (pfd / sp / "RNA_genes.txt").write_text("23S\n" if any(r[0] == "23S" for r in rows) else "")
(pfd / "docs").mkdir(exist_ok=True)
(pfd / "VERSION").write_text("4.1.1\n")
print(f"tesztadatok: {out}")
