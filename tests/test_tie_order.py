"""symmetrize_hits must not depend on the order of the MMseqs2 output (tie of two directions)."""
import random
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from argeement.build import symmetrize_hits  # noqa: E402

# A-B: one direction leaves out a zero-score substitution at the last position (379 aa, 98.15%),
# the other includes it (380 aa, 97.89%); same number of differences (8) and bit score
ROWS = [dict(query="A", target="B", pident=98.15, alnlen=379, qlen=380, tlen=380, qstart=1, qend=379,
             tstart=1, tend=379, evalue=0.0, bits=700, qcov=0.997, tcov=0.997, n_diff=8),
        dict(query="B", target="A", pident=97.89, alnlen=380, qlen=380, tlen=380, qstart=1, qend=380,
             tstart=1, tend=380, evalue=0.0, bits=700, qcov=1.0, tcov=1.0, n_diff=8),
        dict(query="C", target="D", pident=99.0, alnlen=300, qlen=300, tlen=300, qstart=1, qend=300,
             tstart=1, tend=300, evalue=0.0, bits=600, qcov=1.0, tcov=1.0, n_diff=3)]

ref = None
for seed in range(100):
    rows = ROWS[:]
    random.Random(seed).shuffle(rows)
    out = symmetrize_hits(pd.DataFrame(rows))
    key = out[["query", "target", "pident", "alnlen", "qcov", "tcov"]].to_csv(index=False)
    ref = ref or key
    assert key == ref, f"order-dependent result (seed {seed})"
ab = out[(out["query"] == "A") & (out["target"] == "B")].iloc[0]
assert ab["alnlen"] == 380 and ab["pident"] == 97.89, "the longer alignment must be kept"
assert len(out) == 4
print("tie order: ok")
