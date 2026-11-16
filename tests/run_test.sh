#!/usr/bin/env bash
# Runs the whole chain on synthetic data (no downloads).
set -euo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d)
python3 tests/make_test_data.py "$T/v1"
python3 tests/make_test_data.py "$T/v2" --v2
for v in v1 v2; do
    for db in card amrfinder resfinder disinfinder; do
        python3 -m argeement parse $db --raw "$T/$v/$db" --out "$T/$v/std/$db.tsv"
    done
    python3 -m argeement inventory pointfinder --raw "$T/$v/pointfinder" --out "$T/$v/pf.json"
    python3 -m argeement build --std "$T/$v"/std/*.tsv --out "$T/$v/release" --threads 2 \
        --extra-inventory pointfinder="$T/$v/pf.json"
done
python3 -m argeement diff "$T/v1/release" "$T/v2/release" --out "$T/diff"
python3 -m argeement lookup --release "$T/v2/release" --from resfinder --to card \
    --names blaTEM-1B,blaOXA-1,qnrA1,qnrBroken
python3 - "$T" << 'PY'
import os, sys, pandas as pd
T = sys.argv[1]
m = pd.read_csv(f"{T}/v1/release/mappings.tsv", sep="\t")
def rel(a, b):
    r = m[(m.entry_a == a) & (m.entry_b == b)]
    return r.relation.iloc[0] if len(r) else None
assert rel("blaTEM-1B_1_AY458016", "ARO:3000873") == "identical_protein"
assert rel("catX_1_AB000001", "ARO:3004000") == "identical_protein"      # GTG start
assert rel("blaOXA-1_1_HQ170510", "ARO:3001396") == "length_variant"
assert rel("WP_063840001.1", "ARO:3000873") == "same_family"
q = m[(m.entry_a == "qnrA1_1_AY070235") & (m.db_b == "card")]
best = q[q.best_for_a]
assert set(best.name_b) == {"QnrA1", "QnrA2"} and (best.n_diff == 1).all()  # genuine tie
assert not q[q.name_b == "QnrA3"].best_for_a.iloc[0]                    # 3 differences: not best
assert (best.tie_scope == "within_family").all()                        # qnr: within family
assert (m[m.n_equally_good_b == 1].tie_scope == "none").all()
# chaining: A~B and B~C above the family threshold, A~C below -> A and C cannot share a cluster
ch = pd.read_csv(f"{T}/v1/release/entries.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("entry_id")
assert ch.at["WP_000000005.1", "family_cluster"] != ch.at["WP_000000007.1", "family_cluster"]
assert len({ch.at[k, "family_cluster"] for k in ("WP_000000005.1", "WP_000000006.1", "WP_000000007.1")}) == 2
# family label: direct (AMRFinderPlus) and inferred (CARD from the TEM family)
assert ch.at["WP_000027057.1", "family_label"] == "blaTEM" and ch.at["WP_000027057.1", "family_label_source"] == "ncbi"
assert ch.at["ARO:3000873", "family_label"] == "blaTEM" and ch.at["ARO:3000873", "family_label_source"] == "ncbi"   # identical protein
assert ch.at["VanHAX_1_TEST01#orf3", "family_label"] == ""                 # vanX: no AMRFinderPlus relative
assert rel("blaTEM-1B_1_AY458016", "ARO:3000979") == "near_identical"
t = m[(m.entry_a == "blaTEM-1B_1_AY458016") & (m.entry_b == "ARO:3000979")]
assert t.n_diff.iloc[0] == 1                                            # TEM-116: 1 difference
e = pd.read_csv(f"{T}/v1/release/entries.tsv", sep="\t")
assert "WP_000072067.1" not in set(e.entry_id)                          # POINT filtered out
assert not e.entry_id.str.contains("WP_000489247").any()                # mutation header filtered out
assert "WP_002911375.1" not in set(e.entry_id)                          # POINT_DISRUPT filtered out
assert "WP_004193055.1" not in set(e.entry_id)                          # AMR-SUSCEPTIBLE filtered out
t1 = e[e.entry_id == "WP_000027057.1"].iloc[0]
assert t1.ncbi_node == "blaTEM-1" and t1.ncbi_family == "blaTEM"
assert t1.lineage == "blaTEM-1>blaTEM>bla-A>AMR"
assert e[e.entry_id == "WP_000259031.1"].ncbi_family.iloc[0] == "sul1"  # not an allele: its own family
vf = pd.read_csv(f"{T}/v1/release/validation_ncbi_families.tsv", sep="\t").set_index("ncbi_family")
assert vf.at["blaTEM", "status"] == "split" and vf.at["sul1", "status"] == "ok"
assert vf.at["fexA", "ncbi_complete_ident"] == 90.0
vc = pd.read_csv(f"{T}/v1/release/validation_clusters.tsv", sep="\t")
assert (vc[vc.status == "merged"].ncbi_families == "fexA,fexB").any()
import subprocess, io
lk = pd.read_csv(io.StringIO(subprocess.run(
    ["python3", "-m", "argeement", "lookup", "--release", f"{T}/v1/release", "--from", "amrfinder",
     "--to", "card", "--names", "blaTEM,blaTEM-1,nincsilyen"],
    capture_output=True, text=True, check=True).stdout), sep="\t")
famrow = lk[(lk["query"] == "blaTEM") & (lk.status == "family_level")]
assert len(famrow) >= 1 and "TEM beta-lactamase" in famrow.name_b.iloc[0]
assert (lk[lk["query"] == "blaTEM-1"].status == "mapped").all()
assert (lk[lk["query"] == "nincsilyen"].status == "not_found_in_source").all()
rep = open(f"{T}/v2/release/DATABASES.md").read()
assert "| Beta-lactam |" in rep and "`Iminophenazine`" in rep                # harmonised + unmapped term
assert "| `protein variant model` | 1 |  |" in rep                          # excluded model type in the inventory
assert "1 entries" not in rep and "other two" not in rep and "all three" not in rep
from argeement.build import harmonize, load_class_map
tm, amb = load_class_map()
assert harmonize("penicillin beta-lactam", "", tm, amb) == (["Beta-lactam"], [])
assert harmonize("azole antibiotic;cephalosporin;imidazole antibiotic", "", tm, amb) == (["Beta-lactam"], [])
assert "Terms ignored by the harmonization table" in rep and "| CARD | `azole antibiotic` | 2 | TEM beta-lactamase (2) |" in rep
assert "`azole antibiotic`" not in rep.split("### Terms not covered")[1].split("###")[0]
assert harmonize("NONE;BACITRACIN", "", tm, amb) == (["Polymyxin/peptide"], [])   # NONE: not unmapped
assert harmonize("Folate pathway antagonist", "trimethoprim", tm, amb) == (["Trimethoprim"], [])
assert "PointFinder (point mutations, not included)" in rep
assert "WP_000000001.1" in set(e.entry_id) and "WP_000000002.1" not in set(e.entry_id)  # BIOCIDE be, METAL ki
assert rel("qacA_1_AB566410", "WP_000000001.1") == "identical_protein"   # DisinFinder <-> AMRFinderPlus
assert rel("qacA_1_AB566410", "ARO:3005001") == "identical_protein"      # DisinFinder <-> CARD
assert "| `STRESS/BIOCIDE` | 1 | ✓ |" in rep and "| `STRESS/METAL` | 1 |  |" in rep
bio = [l for l in rep.splitlines() if l.startswith("| Biocide/antiseptic |")][0].split("|")
assert [c.strip() for c in bio[2:6]] == ["1", "1", "0", "6"], bio      # AMRF, CARD, ResF, DisinF (qacA, formA, qacX#orf1, sitXY#orf1-2, qacZ#orf)
# displaced gene: ORF selection with curated priority recovers the real gene
import json
qz = m[(m.entry_a.str.startswith("qacZ_1_X06#orf")) & (m.entry_b == "ARO:3005002")]
assert (qz.relation == "identical_protein").any(), qz
man = json.load(open(f"{T}/v1/release/manifest.json"))
assert man["orf_curated_priority"] == ["disinfinder::qacZ_1_X06"], man["orf_curated_priority"]
# (make_test_data.py ensures that selection by coverage alone would displace the gene)
# the AMRFinderPlus scope (core/plus) reaches the release
e1 = pd.read_csv(f"{T}/v1/release/entries.tsv", sep="\t", dtype=str, keep_default_na=False)
assert set(e1.loc[e1.db == "amrfinder", "scope"]) <= {"core", "plus"} and \
    (e1.loc[e1.db == "amrfinder", "scope"] != "").all(), set(e1.loc[e1.db == "amrfinder", "scope"])
assert "| *campylobacter* | 1 | 1 | 2 | 3 | 6 | 1 |" in rep                # 2 genes+23S, 1+4+1 alleles, 1 combination-dependent
ee = pd.read_csv(f"{T}/v2/release/entries.tsv", sep="\t", dtype=str, keep_default_na=False)
assert ee[ee.entry_id == "sul1_2_U12338"].class_h.iloc[0] == "Sulfonamide"  # ambiguous -> file name
# export: relational layout with integer keys, non-redundant, and every row of mappings.tsv can be
# reconstructed from it
X = f"{T}/v2/release/export"
rx = lambda t: pd.read_csv(f"{X}/{t}.tsv", sep="\t", dtype=str, keep_default_na=False)
xdb, xrt, xf, xp, xe, xr = (rx(t) for t in ("databases", "relation_types", "families", "proteins",
                                            "entries", "relations"))
mm = pd.read_csv(f"{T}/v2/release/mappings.tsv", sep="\t", dtype=str, keep_default_na=False)
# primary keys: 1..n, unique; the stable identifiers are unique too
for t, k in ((xdb, "db_id"), (xf, "family_id"), (xp, "pid"), (xe, "eid")):
    assert list(t[k].astype(int)) == list(range(1, len(t) + 1)), k
assert xp.prot_id.is_unique and xf.family_cluster.is_unique and not xe.duplicated(["db_id", "entry_id"]).any()
assert list(xrt.relation) == ["identical_protein", "near_identical", "length_variant", "same_family",
                              "homolog"] and list(xrt.relation_id) == ["0", "1", "2", "3", "4"]
# keys follow the order of the stable identifiers
assert list(xp.prot_id) == sorted(xp.prot_id) and list(xf.family_cluster) == sorted(xf.family_cluster)
# relations: integer keys only, each pair once, pid_a < pid_b, no identical_protein rows
assert list(xr.columns) == ["pid_a", "pid_b", "relation_id", "n_diff", "pident", "cov_a", "cov_b"]
assert not xr.duplicated(["pid_a", "pid_b"]).any() and (xr.pid_a.astype(int) < xr.pid_b.astype(int)).all()
assert set(xr.relation_id) <= {"1", "2", "3", "4"}
# foreign keys and constraints, checked by SQLite after loading with the README snippet
import sqlite3
code = open(f"{X}/README.md").read().split("## Example: loading into SQLite")[1]
code = code.split("```python")[1].split("```")[0].replace('"argeement.sqlite"', '":memory:"')
cwd = os.getcwd(); os.chdir(X)
g = {}; exec(code, g); con = g["con"]
os.chdir(cwd)
con.execute("PRAGMA foreign_keys = ON")
assert con.execute("PRAGMA foreign_key_check").fetchall() == []
assert con.execute("SELECT count(*) FROM relations").fetchone()[0] == len(xr)
assert con.execute("SELECT count(*) FROM entries WHERE pid IS NULL").fetchone()[0] == (xe.pid == "").sum()
assert con.execute("SELECT count(*) FROM entries WHERE parent_eid IS NOT NULL").fetchone()[0] > 0
# ORF sub-entries point to their original entry
orf = con.execute("SELECT c.entry_id, p.entry_id FROM entries c JOIN entries p ON c.parent_eid = p.eid").fetchall()
assert orf and all(c.split("#orf")[0] == p_ for c, p_ in orf)
# family size matches the proteins
assert con.execute("SELECT count(*) FROM families f WHERE n_proteins != "
                   "(SELECT count(*) FROM proteins p WHERE p.family_id = f.family_id)").fetchone()[0] == 0
# the release can be reconstructed: decode the keys back to prot_id and relation names
P = dict(zip(xp.pid, xp.prot_id)); R = dict(zip(xrt.relation_id, xrt.relation))
D = dict(zip(xdb.db_id, xdb.db))
pid = dict(zip(xe.db_id.map(D) + "::" + xe.entry_id, xe.pid.map(P).fillna("")))
ent_rel = pd.read_csv(f"{T}/v2/release/entries.tsv", sep="\t", dtype=str, keep_default_na=False)
assert all(pid[f"{d}::{x}"] == (p_ if u == "True" else "")
           for d, x, p_, u in zip(ent_rel.db, ent_rel.entry_id, ent_rel.prot_id, ent_rel.usable))
rel_ix = {tuple(sorted((P[a_], P[b_]))): (R[r_], int(n_))
          for a_, b_, r_, n_ in zip(xr.pid_a, xr.pid_b, xr.relation_id, xr.n_diff)}
min_nd = {}
for row in mm.itertuples(index=False):
    pa, pb = pid[f"{row.db_a}::{row.entry_a}"], pid[f"{row.db_b}::{row.entry_b}"]
    if row.relation == "identical_protein":
        assert pa == pb, row
    else:
        k = tuple(sorted((pa, pb)))
        assert k in rel_ix and rel_ix[k][0] == row.relation, row
        min_nd[k] = min(min_nd.get(k, 10**9), int(row.n_diff))
# the export keeps the values of the direction with fewer differences
assert all(rel_ix[k][1] == v for k, v in min_nd.items())
used = {tuple(sorted((a_, b_))) for a_, b_ in zip(mm.prot_a, mm.prot_b) if a_ != b_}
assert used == set(rel_ix), "relations.tsv contains a pair that is not in mappings.tsv"
assert len(xr) < len(mm)
import hashlib
for ln in open(f"{X}/SHA256SUMS"):
    h, f = ln.split()
    assert hashlib.sha256(open(f"{X}/{f}", "rb").read()).hexdigest() == h, f
assert {"schema.sql", "LICENSE.md", "databases.tsv", "relation_types.tsv", "families.tsv"} <= \
    {ln.split()[1] for ln in open(f"{X}/SHA256SUMS")}
lic = open(f"{X}/LICENSE.md").read()
assert "CC-BY-4.0" in lic and "McMaster" in lic and "Apache License 2.0" in lic and "public domain" in lic
from argeement.zenodo import metadata
md = metadata(dict(title="t", creators=[]), json.load(open(f"{X}/manifest.json")), "r")
assert md["license"] == "cc-by-4.0" and "McMaster" in md["description"] and "CC-BY 4.0" in md["description"]
# the example code of the README runs in the export directory
code = open(f"{X}/README.md").read().split("```python")[1].split("```")[0]
cwd = os.getcwd(); os.chdir(X)
import contextlib, io as _io
with contextlib.redirect_stdout(_io.StringIO()) as buf:
    exec(code, {})
os.chdir(cwd)
assert "qnrA1_1_AY070235" in buf.getvalue() and "QnrA1" in buf.getvalue()
assert xe.set_index("entry_id").at["WP_000027057.1", "name"] == "blaTEM-1"
# package: reproducible archive and release notes
import subprocess, tarfile, io
pk = lambda out: subprocess.run(["python3", "-m", "argeement", "package", "--release",
                                 f"{T}/v2/release", "--out", out], check=True, capture_output=True, text=True)
r1 = pk(f"{T}/pkg1"); pk(f"{T}/pkg2")
arch = [f for f in os.listdir(f"{T}/pkg1") if f.endswith(".tar.gz")][0]
assert open(f"{T}/pkg1/{arch}", "rb").read() == open(f"{T}/pkg2/{arch}", "rb").read()   # byte-identical
assert "publish_local.sh" in r1.stderr
with tarfile.open(f"{T}/pkg1/{arch}") as tf:
    names = tf.getnames()
root = arch[:-len(".tar.gz")]
assert {f"{root}/{x}" for x in ("databases.tsv", "relation_types.tsv", "families.tsv", "entries.tsv",
                                 "relations.tsv", "proteins.tsv", "schema.sql", "LICENSE.md", "README.md",
                                 "SHA256SUMS", "manifest.json")} <= set(names)
notes = open(f"{T}/pkg1/RELEASE_NOTES.md").read()
assert "| DisinFinder | `2.0.1 (git 09fc831)` |" in notes and "First published release" in notes
# zenodo: first version with --new, then a new version, against a local mock of the deposit API
import http.server, threading, shutil as _sh2
calls = []
class Z(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def _send(self, obj, code=200):
        b = json.dumps(obj).encode(); self.send_response(code)
        self.send_header("Content-Type", "application/json"); self.send_header("Content-Length", str(len(b)))
        self.end_headers(); self.wfile.write(b)
    def _rec(self):
        n = int(self.headers.get("Content-Length") or 0); body = self.rfile.read(n) if n else b""
        calls.append((self.command, self.path, self.headers.get("User-Agent"), body[:200]))
    def do_POST(self):
        self._rec(); B = f"http://127.0.0.1:{srv.server_port}"
        if self.path == "/deposit/depositions":
            self._send({"id": 11, "files": [], "links": {"bucket": f"{B}/files/11"}}, 201)
        elif self.path.endswith("/actions/newversion"):
            self._send({"links": {"latest_draft": f"{B}/deposit/depositions/12"}}, 201)
        elif self.path.endswith("/actions/publish"):
            i = int(self.path.split("/")[3])
            self._send({"id": i, "doi": f"10.5072/zenodo.{i}", "conceptdoi": "10.5072/zenodo.10",
                        "links": {"html": f"{B}/records/{i}"}}, 202)
    def do_GET(self):
        self._rec(); B = f"http://127.0.0.1:{srv.server_port}"
        self._send({"id": 12, "files": [{"links": {"self": f"{B}/old/file"}}], "links": {"bucket": f"{B}/files/12"}})
    def do_PUT(self): self._rec(); self._send({})
    def do_DELETE(self): self._rec(); self.send_response(204); self.end_headers()
srv = http.server.HTTPServer(("127.0.0.1", 0), Z)
threading.Thread(target=srv.serve_forever, daemon=True).start()
os.environ["ZENODO_API"] = f"http://127.0.0.1:{srv.server_port}"; os.environ["ZENODO_TOKEN"] = "x"
from argeement.zenodo import publish
zc = f"{T}/zenodo.yaml"
open(zc, "w").write("deposition_id: null\ntitle: ARGeement\ncreators: [{name: 'A, B'}]\n")
try:
    publish(config=zc, package_dir=f"{T}/pkg1"); raise AssertionError("no id, no --new: must stop")
except SystemExit as zex:
    assert "--new" in str(zex)
r1 = publish(config=zc, package_dir=f"{T}/pkg1", new_record=True)
assert r1["deposition_id"] == 11 and r1["doi"] == "10.5072/zenodo.11"
assert ("POST", "/deposit/depositions") in [c[:2] for c in calls] and not any("newversion" in c[1] for c in calls)
assert all(c[2] and "Python-urllib" not in c[2] for c in calls)      # explicit User-Agent
st = json.load(open(f"{T}/zenodo_state.json")); assert st["deposition_id"] == 11
try:
    publish(config=zc, package_dir=f"{T}/pkg2", new_record=True); raise AssertionError("record exists: --new must stop")
except SystemExit as zex:
    assert "without --new" in str(zex)
_sh2.copy(f"{T}/pkg1/package.json", f"{T}/pkg2/package.json")
pj = json.load(open(f"{T}/pkg2/package.json")); pj["run"] = "second"; json.dump(pj, open(f"{T}/pkg2/package.json", "w"))
calls.clear()
r2 = publish(config=zc, package_dir=f"{T}/pkg2")
assert r2["doi"] == "10.5072/zenodo.12" and any(c[1] == "/deposit/depositions/11/actions/newversion" for c in calls)
assert ("DELETE", "/old/file") in [c[:2] for c in calls]
assert publish(config=zc, package_dir=f"{T}/pkg2")["doi"] == "10.5072/zenodo.12"   # not published twice
srv.shutdown()
# name concordance
from argeement.names import classify
assert classify("blaTEM-1", "amrfinder", "TEM-1", "card") == "formatting"
assert classify("Ecol_emrE", "card", "emrE", "amrfinder") == "species_prefix"
assert classify("blaTEM-1", "amrfinder", "blaTEM-1A", "resfinder") == "variant_suffix"
assert classify("blaVIM-2", "amrfinder", "VIM-4", "card") == "allele_number_differs"
assert classify("mdt", "amrfinder", "Rv1877", "card") == "different_name"
nc = pd.read_csv(f"{T}/v2/release/name_concordance.tsv", sep="\t", dtype=str, keep_default_na=False)
cat = {(r.name_a, r.name_b): r.category for r in nc.itertuples(index=False)}
assert cat[("blaTEM-1", "TEM-1")] == "formatting"
assert cat[("Ecol_catX", "catX")] == "species_prefix"            # card < resfinder
assert cat[("TEM-1", "blaTEM-1B")] == "variant_suffix"
assert cat[("qacEdelta1", "qacA")] == "ncbi_synonym"             # amrfinder < card, NCBI synonym
assert cat[("blaVIM-2", "VIM-4")] == "allele_number_differs"
# the node's own name is not a synonym: this cannot be ncbi_synonym
assert cat[("blaTEM-116", "TEM-116 beta-lactamase")] == "different_name"
assert (nc.category == "ncbi_synonym").sum() == 2                    # qacEdelta1-qacA (CARD and DisinFinder)
assert not ((nc.db_a == nc.db_b) | (nc.db_a > nc.db_b)).any()    # every pair once, ordered
rep2 = open(f"{T}/v2/release/DATABASES.md").read()
assert "## 5. Name concordance" in rep2 and "## 6. Strengths" in rep2 and "## 7. Methods" in rep2
assert "| `blaVIM-2` | `VIM-4` | allele number differs |" in rep2
assert "(see section 6)" in rep2
# ORF-splitting cases
from argeement.sequtil import _best_tiling, choose_translation
# start codon tie: of a GTG 3 codons upstream and the ATG, the ATG wins
ch = _best_tiling([dict(start=0, end=609, atg=False), dict(start=9, end=609, atg=True)])
assert len(ch) == 1 and ch[0]["start"] == 9
# reverse-strand gene with a clean forward frame: must NOT be translated as a whole in the wrong frame
assert rel("aacREVC_1_X02#orf1", "ARO:3003998") == "identical_protein"
pr = e[e.entry_id == "aacREVC_1_X02"].iloc[0]
assert not pr.usable and "split_into_orfs" in pr.issues
# a split parent is never usable (no double counting)
split_par = set(e.loc[e.parent != "", "parent"])
assert not e[e.entry_id.isin(split_par)].usable.any()
# a gene without curated counterpart, recovered as a single ORF, is kept
assert e[e.entry_id == "uniqREV_1_X04#orf1"].usable.iloc[0]
# a non-resistance ORF next to the real gene is excluded
qx = e[e.parent == "qacX_1_X03"].set_index("entry_id")
assert qx.at["qacX_1_X03#orf1", "usable"] and not qx.at["qacX_1_X03#orf2", "usable"]
assert "orf_unsupported" in qx.at["qacX_1_X03#orf2", "issues"]
sx = e[e.parent == "sitXY_1_X05"]
assert len(sx) == 2 and sx.usable.all()                                   # egyedi operon megmarad
# operons and other reading frames: ORF sub-entries
orf = e[e.entry_id.str.startswith("VanHAX_1_TEST01#orf")]
assert len(orf) == 3 and set(orf.parent) == {"VanHAX_1_TEST01"} and (orf.usable == True).all(), orf
assert rel("VanHAX_1_TEST01#orf1", "ARO:3002942") == "identical_protein"          # vanH
assert rel("VanHAX_1_TEST01#orf2", "WP_000000004.1") == "identical_protein"       # vanA
assert rel("blaREV-1_1_X01#orf1", "ARO:3003999") == "identical_protein"           # reverse strand
par = e[e.entry_id == "VanHAX_1_TEST01"].iloc[0]
assert not par.usable and "split into 3 ORF" in par.note
assert not nc.entry_a.str.contains("#orf").any() and not nc.entry_b.str.contains("#orf").any()
assert "| ORF sub-entries |" in rep2.split("\n## 2")[0].replace(" | ORF sub-entries | ", "| ORF sub-entries |")
assert "were split into open reading frames" in rep2
assert "parent_eid" in xe.columns
# query: external sequences
std = pd.read_csv(f"{T}/v1/std/amrfinder.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("entry_id")
vim = std.at["WP_000000003.1", "protein"]
vim_var = vim[:50] + ("A" if vim[50] != "A" else "C") + vim[51:120] + ("W" if vim[120] != "W" else "Y") + vim[121:]
rs = pd.read_csv(f"{T}/v1/std/resfinder.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("entry_id")
with open(f"{T}/q.fa", "w") as fh:
    fh.write(f">tem1\n{std.at['WP_000027057.1', 'protein']}\n>vim_new\n{vim_var}\n")
    fh.write(f">operon_nt\n{rs.at['VanHAX_1_TEST01', 'nucleotide']}\n>random\nMKWQWHHRRPPYYCCWWMMNNQQ\n")
qr = subprocess.run(["python3", "-m", "argeement", "query", "--release", f"{T}/v1/release", "--fasta",
                     f"{T}/q.fa", "--threads", "2"], check=True, capture_output=True, text=True)
q = pd.read_csv(io.StringIO(qr.stdout), sep="\t", dtype=str, keep_default_na=False)
assert "release" in qr.stderr
t1 = q[(q["query"] == "tem1") & (q.db == "card")]
assert (t1.relation == "identical_protein").all() and "TEM-1" in set(t1.name)
vn = q[(q["query"] == "vim_new")].set_index("db")
assert vn.at["amrfinder", "name"] == "blaVIM-2" and vn.at["amrfinder", "n_diff"] == "2"
assert vn.at["card", "name"] == "VIM-4" and vn.at["card", "relation"] == "near_identical"
assert set(q[q["query"].str.startswith("operon_nt#orf")]["query"]) == {f"operon_nt#orf{k}" for k in (1, 2, 3)}
assert (q[q["query"] == "operon_nt#orf2"].relation == "identical_protein").all()          # vanA
assert q[q["query"] == "random"].status.tolist() == ["no_match"]
# mutations per genus, from three sources
assert "Mutation-based determinants by genus" in rep2
assert "| *Escherichia* | 2 | 3 | 1 | 2 | – | – |" in rep2, [l for l in rep2.splitlines() if "Escherichia" in l]  # gyrA + 16S
assert "AMRProt-mutation.tsv + 1 AMR_DNA-* files (rRNA)" in rep2
assert "| *Campylobacter* | – | – | – | – | 2 | 3 |" in rep2
# reverse direction: same name, different protein
sn = pd.read_csv(f"{T}/v2/release/same_name_concordance.tsv", sep="\t", dtype=str, keep_default_na=False)
r_ = {(x.name_a, x.name_b): x for x in sn.itertuples(index=False)}
v2 = r_[("blaVIM-2", "VIM-2")]
assert v2.name_level == "allele" and v2.relation == "near_identical" and v2.n_diff == "5"
assert r_[("mexX", "MexX")].relation == "unrelated"
s1 = sn[(sn.db_a == "amrfinder") & (sn.name_a == "sul1") & (sn.relation == "near_identical")]
assert len(s1) == 1 and s1.name_level.iloc[0] == "gene"                       # gene-level name, expected difference
assert r_[("blaTEM-1", "TEM-1")].relation == "identical_protein"
rep3 = open(f"{T}/v2/release/DATABASES.md").read()
assert "### Same protein, different name" in rep3 and "### Same name, different protein" in rep3
assert "| AMRFinderPlus – CARD | `blaVIM-2` | `VIM-2` | near_identical | 5 |" in rep3
assert "| AMRFinderPlus – CARD | `mexX` | `MexX` | unrelated | – |" in rep3
# number of differences: a substitution at an end is one difference, not two
from argeement.sequtil import n_diff_from_alignment
# 300 aa, substitution at the last position: the alignment ends at position 299 on both
assert n_diff_from_alignment(299, 299, 1, 299, 300, 1, 299, 300) == 1
# N-terminal extension of 5 aa on the query
assert n_diff_from_alignment(300, 300, 6, 305, 305, 1, 300, 300) == 5
# synthetic validation
vr = subprocess.run(["python3", "-m", "argeement", "validate", "--release", f"{T}/v2/release",
                     "--n", "15", "--threads", "2"], check=True, capture_output=True, text=True)
vs = pd.read_csv(f"{T}/v2/release/validation_synthetic/summary.tsv", sep="\t").set_index("kind")
assert (vs["truth_found"] >= 0.95).all(), vs
assert vs.at["cds_reverse", "truth_found"] == 1.0 and vs.at["gene_with_junk_orf", "truth_found"] == 1.0
fx = e[e.entry_id == "WP_015585966.1"]
assert len(fx) == 1 and fx.name.iloc[0] == "fexA"                       # acquired gene is kept
d = pd.read_csv(f"{T}/diff/diff.tsv", sep="\t")
assert {"entry_added", "entry_sequence_changed", "entry_renamed"} <= set(d.change)
# real-genome use case: agreement of tool outputs before and after mapping
subprocess.run(["python3", "tests/make_concordance_data.py", f"{T}/conc_in"], check=True)
subprocess.run(["python3", "-m", "argeement", "concordance", "--release", f"{T}/v2/release",
                "--input", f"{T}/conc_in", "--out", f"{T}/conc"], check=True)
L = pd.read_csv(f"{T}/conc/loci.tsv", sep="\t", keep_default_na=False)
cat = {(r.pair, r.ref_a, r.ref_b): r.category for r in L.itertuples()}
assert cat[("amrfinder-rgi", "WP_000027057.1", "ARO:3000873")] == "identical"
assert cat[("amrfinder-rgi", "WP_063840001.1", "ARO:3000873")] == "same_family"
assert cat[("amrfinder-resfinder", "WP_000259031.1", "sul1_3_EU855787")] == "strong"
assert cat[("amrfinder-rgi", "WP_001334766.1", "")] == "only_amrfinder:counterpart_exists"
assert cat[("amrfinder-rgi", "", "ARO:3005000")] == "only_rgi:no_counterpart"
assert not (L.ref_a == "WP_1").any(), "the point-mutation call is excluded"
tem = L[(L.pair == "amrfinder-resfinder") & (L.ref_a == "WP_000027057.1")]
assert len(tem) == 1 and tem.iloc[0].category == "identical" and int(float(tem.iloc[0].n_refs_b)) == 2, tem
assert not L[L.pair == "amrfinder-resfinder"].category.str.startswith("only_resfinder").any(), \
    "calls at the same locus are merged"
S = pd.read_csv(f"{T}/conc/summary.tsv", sep="\t").set_index("pair")
assert S.loc["amrfinder-resfinder", "name_concordance_before"] == 0.5
assert S.loc["amrfinder-resfinder", "name_concordance_after"] == 1.0

mf = json.load(open(f"{T}/v2/release/manifest.json"))
assert "pointfinder" in mf["database_versions"], mf["database_versions"]
assert mf["counts"]["protein_level_entries"] == mf["counts"]["usable"]

# every protein pair assigned once: hits carry the same values in both directions
ph = pd.read_csv(f"{T}/v2/release/protein_hits.tsv", sep="\t", dtype={"query": str, "target": str})
fw = ph.set_index(["query", "target"])[["pident", "n_diff", "alnlen"]]
bw = ph.rename(columns={"query": "target", "target": "query"}).set_index(["query", "target"])[["pident", "n_diff", "alnlen"]]
assert len(fw) == len(bw) and fw.sort_index().equals(bw.sort_index()), "hits are not symmetric"
qc = ph.set_index(["query", "target"])["qcov"].sort_index()
tc = ph.rename(columns={"query": "target", "target": "query"}).set_index(["query", "target"])["tcov"].sort_index()
assert (qc.values == tc.values).all()

# summary web page (GitHub Pages)
subprocess.run(["python3", "-m", "argeement", "site", "--release", f"{T}/v2/release",
                "--out", f"{T}/site", "--config", "site.yaml"], check=True)
idx = open(f"{T}/site/index.html", encoding="utf-8").read()
assert "<h2>Licence</h2>" in idx and "McMaster" in idx and "creativecommons.org/licenses/by/4.0" in idx
assert "Database versions" in idx and "Overlap between databases" in idx and "Synthetic validation" in idx
sj = json.load(open(f"{T}/site/summary.json"))
assert sj["totals"]["source_entries"] > 0 and sj["database_versions"], sj["totals"]

# update_code.sh: replaces the code, keeps config.yaml, existing profiles and data
import zipfile, shutil as _sh
D = f"{T}/server"
os.makedirs(f"{D}/work/x"); os.makedirs(f"{D}/profiles"); os.makedirs(f"{D}/argeement")
open(f"{D}/work/x/data.txt", "w").write("data")
open(f"{D}/config.yaml", "w").write("threads: 99\n")
open(f"{D}/profiles/amrfinder-2026-08-07.yaml", "w").write("local: /mine\n")
open(f"{D}/profiles/dbs_mine.yaml", "w").write("mine\n")
open(f"{D}/argeement/export.py", "w").write("old\n")
_sh.copy("update_code.sh", f"{D}/update_code.sh")
with zipfile.ZipFile(f"{T}/code.zip", "w") as z:
    for root, dirs, files in os.walk("."):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git", "work", "releases")]
        for f in files:
            z.write(os.path.join(root, f), os.path.join("ARGeement-main", root, f))
r = subprocess.run(["bash", f"{D}/update_code.sh", f"{T}/code.zip"], capture_output=True, text=True)
assert r.returncode == 0, r.stderr
assert open(f"{D}/config.yaml").read() == "threads: 99\n" and os.path.exists(f"{D}/config.yaml.new")
assert open(f"{D}/profiles/amrfinder-2026-08-07.yaml").read() == "local: /mine\n"
assert open(f"{D}/profiles/dbs_mine.yaml").read() == "mine\n"
assert open(f"{D}/work/x/data.txt").read() == "data"
assert open(f"{D}/argeement/export.py").read() == open("argeement/export.py").read()
assert os.path.exists(f"{D}/argeement/terms.py") and "the settings differ" in r.stdout

print("\nAll tests passed.")
PY
rm -rf "$T"
