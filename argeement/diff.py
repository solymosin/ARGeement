"""Comparison of two releases: the report shows what is worth reviewing by hand."""
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

STRONG = {"identical_protein", "near_identical", "length_variant"}


def _load(d):
    d = Path(d)
    e = pd.read_csv(d / "entries.tsv", sep="\t", dtype=str, keep_default_na=False)
    m = pd.read_csv(d / "mappings.tsv", sep="\t", dtype=str, keep_default_na=False)
    man = json.loads((d / "manifest.json").read_text())
    return e, m, man


def _family_mates(entries):
    groups = defaultdict(set)
    for k, f in zip(entries["key"], entries["family_cluster"]):
        if f:
            groups[f].add(k)
    return {k: groups[f] - {k} for k, f in zip(entries["key"], entries["family_cluster"]) if f}


def diff(old_dir, new_dir, out_dir, all_relations=False):
    e0, m0, man0 = _load(old_dir)
    e1, m1, man1 = _load(new_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []

    def add(change, key, detail=""):
        rows.append(dict(change=change, key=key, detail=detail))

    # 1. entries
    old = e0.set_index("key")
    new = e1.set_index("key")
    for k in new.index.difference(old.index):
        add("entry_added", k, new.at[k, "name"])
    for k in old.index.difference(new.index):
        add("entry_removed", k, old.at[k, "name"])
    common = old.index.intersection(new.index)
    for k in common:
        if old.at[k, "prot_id"] != new.at[k, "prot_id"]:
            add("entry_sequence_changed", k, f"{old.at[k, 'prot_id']} -> {new.at[k, 'prot_id']}")
        if old.at[k, "name"] != new.at[k, "name"]:
            add("entry_renamed", k, f"{old.at[k, 'name']} -> {new.at[k, 'name']}")

    # 2. kapcsolatok
    def keyed(m):
        if not all_relations:
            m = m[m["relation"].isin(STRONG)]
        return {(r.db_a + "::" + r.entry_a, r.db_b + "::" + r.entry_b): r.relation
                for r in m.itertuples(index=False)}

    k0, k1 = keyed(m0), keyed(m1)
    for pair in k1.keys() - k0.keys():
        if pair[0] in common and pair[1] in common:   # mappings of new entries are not counted twice
            add("mapping_added", pair[0], f"{k1[pair]} -> {pair[1]}")
    for pair in k0.keys() - k1.keys():
        if pair[0] in common and pair[1] in common:
            add("mapping_removed", pair[0], f"{k0[pair]} -> {pair[1]}")
    for pair in k0.keys() & k1.keys():
        if k0[pair] != k1[pair]:
            add("relation_changed", pair[0], f"{pair[1]}: {k0[pair]} -> {k1[pair]}")

    # 3. family cluster membership (only for entries present in both releases)
    f0, f1 = _family_mates(e0), _family_mates(e1)
    cset = set(common)
    for k in common:
        a, b = f0.get(k, set()) & cset, f1.get(k, set()) & cset
        if a != b:
            gained, lost = sorted(b - a), sorted(a - b)
            add("family_changed", k,
                f"+{len(gained)} {gained[:3]} / -{len(lost)} {lost[:3]}")

    df = pd.DataFrame(rows, columns=["change", "key", "detail"])
    df.to_csv(out_dir / "diff.tsv", sep="\t", index=False)

    counts = df["change"].value_counts().to_dict()
    lines = ["# ARGeement diff", "",
             f"Previous release: {man0['created']}" + (f" (profile: {man0['profile']})" if man0.get("profile") else "") + "  ",
             f"New release: {man1['created']}" + (f" (profile: {man1['profile']})" if man1.get("profile") else ""), "",
             "## Database versions", "", "| database | previous | new |", "|---|---|---|"]
    for db in sorted(set(man0["database_versions"]) | set(man1["database_versions"])):
        v0 = " ".join(str(man0["database_versions"].get(db, "-")).split())
        v1 = " ".join(str(man1["database_versions"].get(db, "-")).split())
        lines.append(f"| {db} | {v0} | {v1} |")
    lines += ["", "## Changes", "", "| type | count |", "|---|---|"]
    order = ["entry_added", "entry_removed", "entry_sequence_changed", "entry_renamed",
             "mapping_added", "mapping_removed", "relation_changed", "family_changed"]
    for c in order:
        lines.append(f"| {c} | {counts.get(c, 0)} |")
    lines += ["", "Rows of type `entry_sequence_changed`, `relation_changed` and `family_changed` are",
              "the main candidates for manual review; details: `diff.tsv`."]
    (out_dir / "diff_summary.md").write_text("\n".join(lines) + "\n")
    return counts
