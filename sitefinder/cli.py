"""Run the site finder from a config file.

  python -m sitefinder config/example_config.yaml

Writes <outdir>/regions.bed (the oligo design pipeline's region_list: chr, oligo_start, oligo_stop, gene) and
<outdir>/sites.tsv (one row per gene: TSS, peak class, per-sample scores, chosen site and why).
A problem with an input stops the run with a message naming the file and the fix (see "Checks" in the README).
"""
from __future__ import annotations

import difflib
import sys
from pathlib import Path

import pandas as pd
import yaml

from .sites import InputError, assign_regions, best_peak_per_gene, choose_sites, load_peaks, regions
from .tss import explain_empty, extract_tss

REQUIRED = ("gtf", "peaks", "outdir")
DEFAULTS = dict(chromosomes=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"], gene_type="protein_coding",
                transcript_type="protein_coding", genes=None, chrom_sizes=None, min_score=0.5, strong_score=0.8,
                tss_window=250, max_summit_distance=1000, region_length=120, min_site_spacing=None, no_peak="tss",
                build_check=True)

# Genome-build check: how much more often a sample's peaks sit at TSSs than the same peaks moved 1 Mb either way.
# Measured on GENCODE v47 TSSs: ATAC peaks on the right build 16-25x (five cell types); ENCODE CTCF ChIP (CD8 T cells)
# 7.7x; the same CD4 and CD14 ATAC peaks lifted over to hg19 2.0x. With fewer than MIN_PEAKS peaks the ratio is too
# noisy to judge (50 random ATAC peaks: 3.6-6.0x), so the check is skipped and the run says so.
MIN_ENRICHMENT = 5
MIN_PEAKS = 500
SHIFT = 1_000_000


def check_config(cfg) -> dict:
    """The settings with defaults filled in. Stops on a missing or misspelt setting or a missing file."""
    if not isinstance(cfg, dict):
        raise InputError("the config is empty or not a list of  key: value  settings (see config/example_config.yaml)")
    missing = [k for k in REQUIRED if not cfg.get(k)]
    if missing:
        raise InputError(f"the config must set {', '.join(missing)} (see config/example_config.yaml)")
    known = [*REQUIRED, *DEFAULTS]
    unknown = [k for k in cfg if k not in known]
    if unknown:
        def hint(k):
            m = difflib.get_close_matches(str(k), known, 1)
            return f"{k!r} (did you mean {m[0]!r}?)" if m else repr(k)
        raise InputError(f"unknown setting {', '.join(map(hint, unknown))}. The settings are: {', '.join(known)}")
    c = {**DEFAULTS, **cfg}
    if not isinstance(c["peaks"], dict):
        raise InputError("peaks must give one file per sample, as  label: path  lines (see config/example_config.yaml)")
    files = [("gtf", c["gtf"]), *((f"peaks {k}", v) for k, v in c["peaks"].items()), ("genes", c["genes"]),
             ("chrom_sizes", c["chrom_sizes"])]
    for what, path in files:
        if path is not None and not Path(str(path)).is_file():
            raise InputError(f"{what}: no such file: {path}")
    if c["no_peak"] not in ("tss", "skip"):
        raise InputError(f"no_peak must be tss or skip, not {c['no_peak']!r}")
    if not c["min_score"] <= c["strong_score"]:
        raise InputError(f"min_score ({c['min_score']}) is above strong_score ({c['strong_score']})")
    return c


def read_genes(path: str, tss: pd.DataFrame, c: dict, out: Path, notes: list) -> pd.DataFrame:
    """The TSS rows of the requested genes. Unmatched genes are written to <outdir>/genes_not_found.txt."""
    wanted = list(dict.fromkeys(g.strip() for g in open(path) if g.strip()))
    if not wanted:
        raise InputError(f"genes: {path} is empty")
    have = set(tss.gene)
    missing = [g for g in wanted if g not in have]
    selection = f"gene_type {c['gene_type']}, transcript_type {c['transcript_type']}"
    if len(missing) == len(wanted):
        why = (" They look like Ensembl IDs; give gene symbols (e.g. GANC)."
               if all(g.upper().startswith("ENS") for g in wanted[:50]) else "")
        raise InputError(f"genes: none of the {len(wanted):,} genes in {path} is in the GTF selection ({selection}).{why}")
    if missing:
        (out / "genes_not_found.txt").write_text("\n".join(missing) + "\n")
        upper = {g.upper(): g for g in have}
        case = [f"{g} -> {upper[g.upper()]}" for g in missing if g.upper() in upper]
        notes.append(f"{len(missing):,} of {len(wanted):,} requested genes are not in the GTF selection ({selection}) and "
                     f"were left out; they are listed in {out / 'genes_not_found.txt'}. "
                     + "Possible causes: "
                     + (f"different letter case ({', '.join(case[:5])}); " if case else "")
                     + "a renamed symbol (e.g. MARCH1 is now MARCHF1; check genenames.org); or not a protein-coding "
                       "gene (set gene_type and transcript_type to null to keep every type).")
    return tss[tss.gene.isin(wanted)].reset_index(drop=True)


def read_chrom_sizes(path: str, tss: pd.DataFrame) -> dict:
    try:
        s = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["chr", "size"])
        ok = not s.empty and pd.api.types.is_integer_dtype(s["size"])
    except ValueError:                                  # fewer than two columns
        ok = False
    if not ok:
        raise InputError(f"chrom_sizes: {path} must have two tab-separated columns, chromosome and length, no header")
    sizes = dict(zip(s.chr, s["size"]))
    past = tss[tss.tss >= tss.chr.map(sizes)]
    if len(past):
        r = past.iloc[0]
        raise InputError(f"chrom_sizes: {len(past):,} TSSs lie past their chromosome's end (e.g. {r.gene} at {r.chr}:"
                         f"{r.tss:,}; {path} gives {r.chr} {sizes[r.chr]:,} bp). The GTF and chrom_sizes are on "
                         f"different genome builds.")
    return sizes


def check_peaks(label: str, p: pd.DataFrame, cfg: dict, sizes: dict | None):
    """Stops if the scores are not on the scale the thresholds assume, or peaks run past a chromosome's end."""
    if p.score.max() > 1 and not {"min_score", "strong_score"} & set(cfg):
        raise InputError(f"peaks {label}: scores go up to {p.score.max():g}, so they are not LanceOtron's 0-1 peak score, "
                         f"and the default min_score (0.5) and strong_score (0.8) mean nothing on this scale. Set "
                         f"min_score and strong_score in the config for this score.")
    if sizes:
        past = p[p.end > p.chrom.map(sizes)]
        if len(past):
            r = past.iloc[0]
            raise InputError(f"peaks {label}: {len(past):,} peaks run past their chromosome's end (e.g. {r.chrom}:{r.start:,}"
                             f"-{r.end:,}; chrom_sizes gives {r.chrom} {sizes[r.chrom]:,} bp). The peaks and chrom_sizes "
                             f"are on different genome builds.")


def tss_enrichment(n_hit: int, tss: pd.DataFrame, p: pd.DataFrame, c: dict, sizes: dict | None) -> float:
    """TSSs with a peak, over TSSs with a peak once every peak is moved SHIFT bp (mean of both directions)."""
    def hits(shift):
        q = p.assign(start=p.start + shift, end=p.end + shift, summit=p.summit + shift)
        return len(best_peak_per_gene(tss, q, c["tss_window"], c["max_summit_distance"], sizes))
    chance = (hits(SHIFT) + hits(-SHIFT)) / 2
    return (n_hit + 1) / (chance + 1)


def run(cfg: dict) -> pd.DataFrame:
    c = check_config(cfg)
    out = Path(c["outdir"])
    out.mkdir(parents=True, exist_ok=True)
    notes = []
    all_tss = extract_tss(c["gtf"], c["chromosomes"], c["gene_type"], c["transcript_type"])
    if all_tss.empty:
        raise InputError(f"gtf: no genes found in {c['gtf']}: "
                         + explain_empty(c["gtf"], c["chromosomes"], c["gene_type"], c["transcript_type"]))
    print(f"TSS: {len(all_tss):,} genes ({all_tss.source.value_counts().to_dict()})")
    tss = read_genes(c["genes"], all_tss, c, out, notes) if c["genes"] else all_tss
    if c["genes"]:
        print(f"genes: {len(tss):,} requested genes found")
    sizes = read_chrom_sizes(c["chrom_sizes"], all_tss) if c["chrom_sizes"] else None
    keep = set(tss.gene)
    hits = {}
    for label, path in c["peaks"].items():
        p = load_peaks(path, c["min_score"], c["chromosomes"])
        check_peaks(label, p, cfg, sizes)
        found = best_peak_per_gene(all_tss, p, c["tss_window"], c["max_summit_distance"], sizes)
        line = f"{label}: {len(p):,} peaks >= {c['min_score']}; peak at {len(found):,} of {len(all_tss):,} TSSs"
        if len(p) >= MIN_PEAKS:
            e = tss_enrichment(len(found), all_tss, p, c, sizes)
            line += f", {e:.0f}x the rate expected by chance"
            if c["build_check"] and e < MIN_ENRICHMENT:
                raise InputError(
                    f"peaks {label}: {path} has peaks at TSSs only {e:.1f}x as often as chance. ATAC peaks on the right "
                    f"genome build give about 20x and CTCF ChIP about 8x; peaks on a different build than the GTF give "
                    f"about 2x. Check that the peaks and the GTF are on the same build (e.g. not hg19 peaks with a "
                    f"GRCh38 GTF), and lift the peaks over if not. If this file is right, set build_check: false.")
        else:
            notes.append(f"{label} has only {len(p):,} peaks, too few to check that it is on the GTF's genome build.")
        print(line)
        hits[label] = {g: v for g, v in found.items() if g in keep}
    sites = choose_sites(tss, hits, c["strong_score"], c["min_score"])
    spacing = c["region_length"] if c["min_site_spacing"] is None else c["min_site_spacing"]
    sites = assign_regions(sites, spacing, c["no_peak"])
    sites.to_csv(out / "sites.tsv", sep="\t", index=False)
    r = regions(sites, c["region_length"])
    r.to_csv(out / "regions.bed", sep="\t", index=False)
    shared = sites[sites.region_genes.str.contains(",")]
    skipped = int((sites.region == "").sum())
    print(f"classes: {sites.peak_class.value_counts().to_dict()}")
    print(f"regions: {len(r):,} for {len(sites) - skipped:,} genes; "
          + (f"{len(shared):,} genes share {shared.region.nunique():,} regions (sites < {spacing} bp apart, e.g. "
             f"bidirectional promoters; see region_genes in sites.tsv)" if spacing else
             "min_site_spacing 0, so every gene has its own region, even where regions overlap")
          + (f"; {skipped:,} genes with no peak skipped" if skipped else
             f"; {int((sites.peak_class == 'none').sum()):,} genes with no peak placed on their TSS (no_peak: skip leaves them out)"))
    for n in notes:
        print(f"NOTE: {n}")
    print(f"-> {out / 'regions.bed'}, {out / 'sites.tsv'}")
    return sites


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        sys.exit(__doc__)
    try:
        with open(argv[0]) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        sys.exit(f"ERROR: no such config file: {argv[0]}")
    except yaml.YAMLError as e:
        sys.exit(f"ERROR: {argv[0]} is not valid YAML: {e}")
    try:
        run(cfg)
    except InputError as e:
        sys.exit(f"ERROR: {e}")


if __name__ == "__main__":
    main()
