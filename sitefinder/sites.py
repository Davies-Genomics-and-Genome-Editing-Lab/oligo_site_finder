"""Choose one capture site per gene from accessibility peaks in one or more samples.

For each sample, a peak counts at a gene if it scores >= min_score, overlaps TSS +/- tss_window and has its summit
within max_summit_distance of the TSS; the best-scoring such peak is kept. The gene's class is the best score over
samples: strong (>= strong_score), weak (min_score to strong_score) or none.

The site is the summit the samples agree on. Among the kept peaks scoring >= strong_score (all kept peaks if none
do), each peak's support is the number of strong peaks overlapping it, itself included when strong; the most
supported peak wins, then the one whose summit is closest to the TSS, then the first sample in label order.
A gene with no peak gets its TSS.
"""
from __future__ import annotations

import bisect

import numpy as np
import pandas as pd

PEAK_COLUMNS = ["chrom", "start", "end", "score", "summit"]
# header names accepted for each column, compared case-insensitively. LanceOtron's callPeaks writes chrom, start, end,
# overall_peak_score, ..., max_coverage_position; newer LanceOtron writes Chromosome, Start, End, Region_id, Summit,
# Overall_peak_score, ...; a plain file can use the five names above.
COLUMN_NAMES = {"chrom": ["chrom", "chromosome", "chr"], "start": ["start"], "end": ["end"],
                "score": ["score", "overall_peak_score"], "summit": ["summit", "max_coverage_position"]}


class InputError(ValueError):
    """A problem with an input file or setting. The message names the file and how to fix it."""


def _is_number(s) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def load_peaks(path: str, min_score: float, chromosomes: list[str]) -> pd.DataFrame:
    """Peaks from a tab-separated file with a header line: LanceOtron output as it comes, or chrom start end score summit.

    Stops with an InputError naming the fix if the file has no header (e.g. MACS2 narrowPeak), its chromosome names
    do not match (e.g. 1 instead of chr1), its summits are not genome positions, or no peak passes min_score.
    """
    p = pd.read_csv(path, sep="\t", low_memory=False)
    cols = [str(c) for c in p.columns]
    if len(cols) >= 3 and _is_number(cols[1]) and _is_number(cols[2]):       # the first line is a peak, not a header
        fix = ""
        if len(cols) == 10:
            fix = ("\nIt looks like MACS2 narrowPeak, whose summit (column 10) is an offset from start. Convert it with\n"
                   f"  (printf 'chrom\\tstart\\tend\\tscore\\tsummit\\n'; awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $7, $2 + $10}}' "
                   f"{path}) > {path}.tsv\nand set min_score and strong_score on the scale of column 7 (signalValue).")
        raise InputError(f"{path}: no header line; the first line is a peak ({' '.join(cols[:3])} ...). The file needs a "
                         f"header naming chrom, start, end, score, summit (LanceOtron output has one).{fix}")
    lower = {c.lower(): c for c in p.columns}
    use, missing = {}, []
    for col, names in COLUMN_NAMES.items():
        found = [lower[n] for n in names if n in lower]
        if len(found) > 1:
            raise InputError(f"{path}: more than one column could be {col}: {found}")
        if found:
            use[found[0]] = col
        else:
            missing.append(col)
    if missing:
        raise InputError(f"{path}: no column for {missing}. Need a header line naming chrom, start, end, score, summit "
                         f"(LanceOtron output works as it is); found {list(p.columns)}")
    p = p.rename(columns=use)[PEAK_COLUMNS]
    for col in ("start", "end", "score", "summit"):
        if not pd.api.types.is_numeric_dtype(p[col]):
            bad = p[col][pd.to_numeric(p[col], errors="coerce").isna()].iloc[0]
            raise InputError(f"{path}: the {col} column is not all numbers (e.g. {bad!r})")
    p["chrom"] = p.chrom.astype(str)
    on = p.chrom.isin(chromosomes)
    if not on.any():
        names = list(p.chrom.value_counts().index[:4])
        fix = ""
        if any(f"chr{n}" in chromosomes for n in names):
            fix = ("\nThey lack the 'chr' prefix (Ensembl style). Add it with\n"
                   f"  awk 'BEGIN{{FS = OFS = \"\\t\"}} NR > 1 {{$1 = \"chr\" $1}} 1' {path} > {path}.chr.tsv")
        raise InputError(f"{path}: no peak is on the chosen chromosomes ({', '.join(chromosomes[:3])}, ...); this file's "
                         f"chromosomes are {names}.{fix}")
    outside = ~p.summit.between(p.start, p.end)
    if outside.mean() > 0.01:
        raise InputError(f"{path}: {outside.mean():.0%} of summits lie outside their own peak, so the summit column is "
                         f"not the summit's genome position. If it is an offset from start (as in MACS2 narrowPeak "
                         f"column 10), use start + offset.")
    kept = p[on & (p.score >= min_score)]
    if kept.empty:
        raise InputError(f"{path}: no peak on the chosen chromosomes scores >= min_score {min_score} (scores here run "
                         f"from {p.score.min():g} to {p.score.max():g}). Set min_score and strong_score for this score.")
    return kept.reset_index(drop=True)


def best_peak_per_gene(tss: pd.DataFrame, peaks: pd.DataFrame, tss_window: int, max_summit_distance: int,
                       chrom_sizes: dict | None = None) -> dict:
    """{gene: dict(score, summit, start, end)} for genes with a qualifying peak in this sample."""
    out = {}
    for ch, genes in tss.groupby("chr"):
        p = peaks[peaks.chrom == ch].sort_values("start", kind="stable")
        if p.empty:
            continue
        st, en, sc, su = (p[c].to_numpy() for c in ("start", "end", "score", "summit"))
        longest = int((en - st).max())
        for gene, t in zip(genes.gene, genes.tss):
            lo = max(0, t - tss_window)
            hi = t + tss_window if chrom_sizes is None else min(t + tss_window, chrom_sizes.get(ch, np.inf))
            k = np.arange(np.searchsorted(st, lo - longest), np.searchsorted(st, hi))   # peaks starting before the window ends ...
            k = k[(en[k] > lo) & (np.abs(su[k] - t) <= max_summit_distance)]   # ... that reach into it
            if len(k):
                b = k[np.argmax(sc[k])]
                out[gene] = dict(score=float(sc[b]), summit=int(su[b]), start=int(st[b]), end=int(en[b]))
    return out


def choose_sites(tss: pd.DataFrame, hits_by_sample: dict, strong_score: float, min_score: float) -> pd.DataFrame:
    """One row per gene: class, per-sample scores, the chosen site and why."""
    samples = sorted(hits_by_sample)
    rows = []
    for g in tss.itertuples(index=False):
        hits = [dict(sample=s, **hits_by_sample[s][g.gene]) for s in samples if g.gene in hits_by_sample[s]]
        row = dict(gene=g.gene, chr=g.chr, tss=g.tss, strand=g.strand, tss_source=g.source, gene_id=g.gene_id,
                   transcript_id=g.transcript_id, **{f"score_{s}": hits_by_sample[s].get(g.gene, {}).get("score", 0.0)
                                                     for s in samples})
        best = max((h["score"] for h in hits), default=0.0)
        row["max_score"] = best
        row["peak_class"] = "strong" if best >= strong_score else "weak" if best >= min_score else "none"
        row["n_samples_with_peak"] = len(hits)
        if hits:
            strong = [h for h in hits if h["score"] >= strong_score]
            for h in hits:
                h["support"] = sum(o["score"] >= strong_score and (o is h or (h["start"] < o["end"] and h["end"] > o["start"]))
                                   for o in hits)
                h["dist"] = abs(h["summit"] - g.tss)
            pool = strong or hits
            top = max(h["support"] for h in pool)
            c = min((h for h in pool if h["support"] == top), key=lambda h: h["dist"])
            row.update(site=c["summit"], site_source="consensus_summit", site_sample=c["sample"],
                       site_score=c["score"], site_support=c["support"])
        else:
            row.update(site=g.tss, site_source="tss", site_sample="", site_score=0.0, site_support=0)
        # offset of the site from the TSS along the gene: negative = upstream
        row["site_offset_from_tss"] = (row["site"] - g.tss) * (1 if g.strand == "+" else -1)
        rows.append(row)
    return pd.DataFrame(rows)


CLASS_RANK = {"strong": 2, "weak": 1, "none": 0}


def assign_regions(sites: pd.DataFrame, min_spacing: int, no_peak: str = "tss") -> pd.DataFrame:
    """Give each gene a region, so that no two regions' sites are closer than min_spacing.

    no_peak: "tss" keeps genes with no peak (their site is the TSS); "skip" leaves them without a region.
    Genes are taken best first: peak class, then support, then score, then the site closest to the gene's own TSS,
    then gene name. A gene gets its own region if no region already taken is within min_spacing of its site;
    otherwise it shares the nearest such region (e.g. two genes at one bidirectional promoter get one oligo).
    min_spacing 0 gives every gene its own region, even where regions overlap.

    Adds three columns: region (the gene the region is named after; "" if skipped), region_site (its centre) and
    region_genes (every gene in that region, comma-separated).
    """
    if no_peak not in ("tss", "skip"):
        raise ValueError(f"no_peak must be 'tss' or 'skip', not {no_peak!r}")
    s = sites.copy()
    use = s.peak_class != "none" if no_peak == "skip" else pd.Series(True, index=s.index)
    order = s[use].assign(_c=s.peak_class.map(CLASS_RANK), _d=(s.site - s.tss).abs()).sort_values(
        ["_c", "site_support", "site_score", "_d", "gene"], ascending=[False, False, False, True, True], kind="stable")
    taken = {}                     # chr -> [sorted sites, region named at each, order in which each was taken]
    region, region_site = {}, {}
    for g in order.itertuples(index=False):
        pos, name, rank = taken.setdefault(g.chr, ([], [], []))
        k = bisect.bisect_left(pos, g.site)
        near = [j for j in (k - 1, k) if 0 <= j < len(pos) and abs(pos[j] - g.site) < min_spacing]
        if near:                   # share the nearest region; on a tie, the one taken first
            j = min(near, key=lambda j: (abs(pos[j] - g.site), rank[j]))
            region[g.gene], region_site[g.gene] = name[j], pos[j]
        else:
            pos.insert(k, g.site)
            name.insert(k, g.gene)
            rank.insert(k, len(rank))
            region[g.gene], region_site[g.gene] = g.gene, g.site
    s["region"] = s.gene.map(region).fillna("")
    s["region_site"] = s.gene.map(region_site).astype("Int64")
    members = s[s.region != ""].groupby("region").gene.agg(lambda x: ",".join(sorted(x)))
    s["region_genes"] = s.region.map(members).fillna("")
    return s


def regions(sites: pd.DataFrame, length: int) -> pd.DataFrame:
    """A `length`-bp region centred on each region's site, in the oligo design pipeline's region_list format.

    Only genes that carry a region (region == gene) are written; the others share one or were skipped.
    """
    sites = sites[sites.region == sites.gene]
    r = pd.DataFrame({"chr": sites.chr, "oligo_start": sites.site - length // 2, "oligo_stop": sites.site - length // 2 + length,
                      "gene": sites.gene})
    if (r.oligo_start < 0).any():
        raise ValueError(f"site within {length // 2} bp of a chromosome start: {list(r.gene[r.oligo_start < 0])}")
    return r
