"""Synthetic AMRFinderPlus, RGI and ResFinder output for testing the concordance command.
The references are entries of the test release (make_test_data.py --v2)."""
import json
import sys
from pathlib import Path

out = Path(sys.argv[1]) / "G1"
out.mkdir(parents=True, exist_ok=True)
amr_cols = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name",
            "Scope", "Type", "Subtype", "Class", "Subclass", "Method", "Target length",
            "Reference sequence length", "% Coverage of reference", "% Identity to reference",
            "Alignment length", "Closest reference accession", "Closest reference name",
            "HMM accession", "HMM description"]


def amr(s, e, sym, acc, typ="AMR", sub="AMR", meth="EXACTX"):
    r = dict.fromkeys(amr_cols, "")
    r.update({"Protein id": "NA", "Contig id": "ctg1", "Start": s, "Stop": e, "Strand": "+",
              "Element symbol": sym, "Type": typ, "Subtype": sub, "Method": meth,
              "Closest reference accession": acc})
    return "\t".join(r[c] for c in amr_cols)


(out / "amrfinder.tsv").write_text("\t".join(amr_cols) + "\n" + "\n".join([
    amr("1000", "1861", "blaTEM-1", "WP_000027057.1"),        # all three: identical protein
    amr("3000", "3800", "sul1", "WP_000259031.1"),            # ResFinder: near-identical, RGI: none
    amr("5000", "5800", "blaTEM-999", "WP_063840001.1"),      # RGI TEM-1: same family
    amr("12000", "12800", "blaOXA-1", "WP_001334766.1"),      # AMRFinderPlus only, counterpart exists
    amr("20000", "20300", "gyrA_S83L", "WP_1", sub="POINT", meth="POINTX")]) + "\n")  # kimarad

rgi_cols = ["ORF_ID", "Contig", "Start", "Stop", "Orientation", "Best_Hit_ARO", "ARO", "Model_type"]


def rgi(ctg, s, e, name, aro, mt="protein homolog model"):
    return "\t".join([ctg, ctg, s, e, "+", name, aro, mt])


(out / "rgi.txt").write_text("\t".join(rgi_cols) + "\n" + "\n".join([
    rgi("ctg1_3", "1000", "1861", "TEM-1", "3000873"),
    rgi("ctg1_7", "5010", "5800", "TEM-1", "3000873"),
    rgi("ctg1_9", "10000", "10900", "onlyInCard", "3005000"),            # only in CARD
    rgi("ctg1_11", "20000", "20300", "gyrA", "3000000", "protein variant model")]) + "\n")

(out / "resfinder.json").write_text(json.dumps({"seq_regions": {
    "a": {"name": "blaTEM-1B", "ref_id": "blaTEM-1B_1_AY458016", "query_id": "ctg1",
          "query_start_pos": 1000, "query_end_pos": 1861, "ref_database": ["ResFinder-2.7.0"]},
    # a second, equally good reference at the same locus (how ResFinder reports a tie);
    # must count as one locus
    "a2": {"name": "blaOXA-1", "ref_id": "blaOXA-1_1_HQ170510", "query_id": "ctg1",
           "query_start_pos": 1001, "query_end_pos": 1861, "ref_database": ["ResFinder-2.7.0"]},
    "b": {"name": "sul1", "ref_id": "sul1_3_EU855787", "query_id": "ctg1",
          "query_start_pos": 3000, "query_end_pos": 3800, "ref_database": ["ResFinder-2.7.0"]}}}))
