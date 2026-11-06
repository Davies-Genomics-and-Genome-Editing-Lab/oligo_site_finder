"""One TSS per gene from a GENCODE GTF (plain or .gz).

Per gene name: the MANE Select transcript, else the longest APPRIS principal_1 transcript, else the longest
transcript ("longest" = genomic span). Only transcript features of the chosen gene/transcript types on the chosen
chromosomes are read. Genes are grouped by NAME, so a symbol with several gene IDs (e.g. readthroughs) gets one TSS.
Pseudoautosomal genes present on both chrX and chrY keep the chrX copy.
"""
from __future__ import annotations

import gzip
import re
from collections import Counter, defaultdict

import pandas as pd

_ATTR = re.compile(r'^(\S+)\s+"([^"]*)"')


def open_text(path):
    """Open a plain or gzipped (.gz) text file for reading."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def explain_empty(gtf: str, chromosomes: list[str], gene_type, transcript_type, n: int = 20000) -> str:
    """Why extract_tss found no transcripts, from the first n feature lines of the GTF."""
    chroms, feats, keys, types = Counter(), set(), set(), Counter()
    with open_text(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chroms[f[0]] += 1
            feats.add(f[2])
            a = _attributes(f[8])
            keys |= set(a)
            types[a.get("gene_type", a.get("gene_biotype", "?"))] += 1
            if sum(chroms.values()) >= n:
                break
    if not chroms:
        return "it has no GTF lines (nine tab-separated columns). Is it a GTF?"
    names = list(chroms)[:3]
    if not set(chroms) & set(chromosomes):
        if not any(c.startswith("chr") for c in chroms) and any(c.startswith("chr") for c in chromosomes):
            return (f"its chromosomes are named {names}, not chr1, chr2, ... (Ensembl style). Use a GENCODE GTF, "
                    f"which uses chr names.")
        return f"none of its chromosomes ({names} ...) is in the chromosomes setting ({chromosomes[:3]} ...)."
    if "transcript" not in feats:
        return f"it has no 'transcript' lines (its feature types include {sorted(feats)[:5]})."
    lacking = [k for k in ("gene_name", "gene_type", "transcript_type") if k not in keys]
    if lacking:
        return (f"its lines lack the {', '.join(lacking)} attribute(s)"
                + (" (an Ensembl GTF, which calls them gene_biotype and transcript_biotype)" if "gene_biotype" in keys else "")
                + ". Use a GENCODE GTF.")
    return (f"no transcript has gene_type {gene_type!r} and transcript_type {transcript_type!r}. Gene types in the file "
            f"include {[t for t, _ in types.most_common(5)]}; set gene_type and transcript_type to one of these, or to null.")


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
    with open_text(gtf) as fh:
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
