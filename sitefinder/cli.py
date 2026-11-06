"""Run the site finder from a config file.

  python -m sitefinder config/example_config.yaml

Writes <outdir>/regions.bed (the oligo design pipeline's region_list: chr, oligo_start, oligo_stop, gene) and
<outdir>/sites.tsv (one row per gene: TSS, peak class, per-sample scores, chosen site and why).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

from .sites import assign_regions, best_peak_per_gene, choose_sites, load_peaks, regions
from .tss import extract_tss

DEFAULTS = dict(chromosomes=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"], gene_type="protein_coding",
                transcript_type="protein_coding", genes=None, chrom_sizes=None, min_score=0.5, strong_score=0.8,
                tss_window=250, max_summit_distance=1000, region_length=120, min_site_spacing=None, no_peak="tss")


def run(cfg: dict) -> pd.DataFrame:
    c = {**DEFAULTS, **cfg}
    out = Path(c["outdir"])
    out.mkdir(parents=True, exist_ok=True)
    tss = extract_tss(c["gtf"], c["chromosomes"], c["gene_type"], c["transcript_type"])
    if tss.empty:
        sys.exit(f"No transcripts found in {c['gtf']} on {c['chromosomes'][:3]}... with gene_type={c['gene_type']} and "
                 f"transcript_type={c['transcript_type']}. Use a GENCODE GTF (chr-prefixed names, gene_type and "
                 f"transcript_type attributes), or change chromosomes / set the two types to null.")
    print(f"TSS: {len(tss):,} genes ({tss.source.value_counts().to_dict()})")
    if c["genes"]:
        wanted = [g.strip() for g in open(c["genes"]) if g.strip()]
        missing = sorted(set(wanted) - set(tss.gene))
        if missing:
            print(f"WARNING: {len(missing)} requested genes not in the GTF selection: {missing[:20]}")
        tss = tss[tss.gene.isin(wanted)].reset_index(drop=True)
    sizes = None
    if c["chrom_sizes"]:
        s = pd.read_csv(c["chrom_sizes"], sep="\t", header=None, names=["chr", "size"])
        sizes = dict(zip(s.chr, s["size"]))
    hits = {}
    for label, path in c["peaks"].items():
        p = load_peaks(path, c["min_score"], c["chromosomes"])
        hits[label] = best_peak_per_gene(tss, p, c["tss_window"], c["max_summit_distance"], sizes)
        print(f"{label}: {len(p):,} peaks >= {c['min_score']}; peak at {len(hits[label]):,} TSSs")
    sites = choose_sites(tss, hits, c["strong_score"], c["min_score"])
    spacing = c["region_length"] if c["min_site_spacing"] is None else c["min_site_spacing"]
    sites = assign_regions(sites, spacing, c["no_peak"])
    sites.to_csv(out / "sites.tsv", sep="\t", index=False)
    r = regions(sites, c["region_length"])
    r.to_csv(out / "regions.bed", sep="\t", index=False)
    shared = sites[sites.region_genes.str.contains(",")]
    skipped = int((sites.region == "").sum())
    print(f"classes: {sites.peak_class.value_counts().to_dict()}")
    print(f"regions: {len(r):,} for {len(sites) - skipped:,} genes; {len(shared):,} genes share {shared.region.nunique():,} "
          f"regions (sites < {spacing} bp apart, e.g. bidirectional promoters; see region_genes in sites.tsv)"
          + (f"; {skipped:,} genes with no peak skipped" if skipped else
             f"; {int((sites.peak_class == 'none').sum()):,} genes with no peak placed on their TSS (no_peak: skip leaves them out)"))
    print(f"-> {out / 'regions.bed'}, {out / 'sites.tsv'}")
    return sites


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        sys.exit(__doc__)
    run(yaml.safe_load(open(argv[0])))


if __name__ == "__main__":
    main()
