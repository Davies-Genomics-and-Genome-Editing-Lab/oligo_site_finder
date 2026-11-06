"""One TSS per gene from a GENCODE GTF.

Per gene name: the MANE Select transcript, else the longest APPRIS principal_1 transcript, else the longest
transcript ("longest" = genomic span). Only transcript features of the chosen gene/transcript types on the chosen
chromosomes are read. Genes are grouped by NAME, so a symbol with several gene IDs (e.g. readthroughs) gets one TSS.
Pseudoautosomal genes present on both chrX and chrY keep the chrX copy.
"""
from __future__ import annotations

import re
from collections import defaultdict

import pandas as pd

_ATTR = re.compile(r'^(\S+)\s+"([^"]*)"')


def _attributes(s: str) -> dict:
    out = {}
    for field in s.strip().rstrip(";").split(";"):
        m = _ATTR.match(field.strip())
        if m:
            k, v = m.groups()
            if k in out:                        # repeated keys (tag) become lists
                out[k] = (out[k] if isinstance(out[k], list) else [out[k]]) + [v]
            else:
                out[k] = v
    return out


def _tags(attrs: dict) -> set:
    t = attrs.get("tag", [])
    return {t} if isinstance(t, str) else set(t)


def extract_tss(gtf: str, chromosomes: list[str], gene_type: str | None = "protein_coding",
                transcript_type: str | None = "protein_coding") -> pd.DataFrame:
    """Return one row per gene: gene, chr, tss (0-based), strand, source, gene_id, transcript_id."""
    chroms = set(chromosomes)
    by_name = defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript" or f[0] not in chroms:
                continue
            a = _attributes(f[8])
            if gene_type and a.get("gene_type") != gene_type:
                continue
            if transcript_type and a.get("transcript_type") != transcript_type:
                continue
            start, end = int(f[3]), int(f[4])                       # GTF: 1-based, inclusive
            by_name[a.get("gene_name", "")].append(dict(
                chr=f[0], start=start, end=end, strand=f[6], tags=_tags(a), span=end - start + 1,
                gene_id=a.get("gene_id", ""), transcript_id=a.get("transcript_id", "")))

    rows = []
    for gene in sorted(by_name):
        tx = by_name[gene]
        if {t["chr"] for t in tx} >= {"chrX", "chrY"}:                # PAR duplicate: keep chrX
            tx = [t for t in tx if t["chr"] != "chrY"]
        mane = [t for t in tx if "MANE_Select" in t["tags"]]
        appris = [t for t in tx if "appris_principal_1" in t["tags"]]
        if mane:
            chosen, source = mane[0], "MANE_Select"
        elif appris:
            chosen, source = max(appris, key=lambda t: t["span"]), "appris_principal_1"
        else:
            chosen, source = max(tx, key=lambda t: t["span"]), "longest_transcript"
        tss_1based = chosen["start"] if chosen["strand"] == "+" else chosen["end"]
        rows.append(dict(gene=gene, chr=chosen["chr"], tss=tss_1based - 1, strand=chosen["strand"], source=source,
                         gene_id=chosen["gene_id"], transcript_id=chosen["transcript_id"]))
    return pd.DataFrame(rows)
