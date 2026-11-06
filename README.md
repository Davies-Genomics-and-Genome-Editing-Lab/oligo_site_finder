# oligo_site_finder

Makes the regions file for the lab's [oligo design pipeline](https://github.com/Davies-Genomics-and-Genome-Editing-Lab/oligo_design_pipeline)
when you want capture oligos at **promoters**. For each gene it picks one site, the ATAC peak summit at the gene's
TSS (or the TSS itself if there is no peak), and writes a 120 bp region centred on it. Genes that share a promoter
(e.g. bidirectional pairs) share one region, so no two oligos in your panel overlap.

The oligo design pipeline used to do this step itself (its `promoter` option). That option was retired in August
2026, and the pipeline now needs a ready-made regions file (`region_list`). This repo makes that file. The oligo
design pipeline then designs and filters the oligos.

It was written for the 2026 genome-wide Micro Capture-C promoter panels. Given the same GTF and peak files, it
reproduces the sites of all 19,704 protein-coding promoters in that design.

## Quick start

```bash
git clone https://github.com/Davies-Genomics-and-Genome-Editing-Lab/oligo_site_finder.git
cd oligo_site_finder
conda env create -f environment.yaml        # small: python, pandas, numpy, pyyaml, pytest
conda activate oligo_site_finder
cp config/example_config.yaml my_config.yaml  # set gtf, peaks and outdir
python -m sitefinder my_config.yaml
```

Then, in the oligo design pipeline's `config/config.yaml`:

```yaml
mode: single_best
region_list: /path/to/outdir/regions.bed
```

## Inputs

| Input | What it is |
|---|---|
| `gtf` | A **GENCODE** GTF for your genome build, e.g. `gencode.v47.primary_assembly.annotation.gtf` from [gencodegenes.org](https://www.gencodegenes.org/human/). It needs chr-prefixed names and the `gene_name`, `gene_type` and `transcript_type` attributes; its `MANE_Select` and `appris_principal_1` tags pick the TSS. Ensembl and RefSeq GTFs are not supported (the run stops with a message if no genes are found). |
| `peaks` | One ATAC peak file per sample, as `label: path`. See [Peak files](#peak-files). |
| `genes` | Optional. A text file with one gene symbol per line. Leave it out to get every protein-coding gene. |
| `chrom_sizes` | Optional. Two columns, chromosome and length. It keeps TSS windows inside chromosome ends. |

The GTF and the peaks must be on the same genome build.

### Peak files

**LanceOtron output works as it is.** Both of its layouts are recognised from the header line, in any letter case:

- `chrom start end overall_peak_score ... max_coverage_position ...` (LanceOtron `callPeaks`)
- `Chromosome Start End Region_id Summit Overall_peak_score ...` (newer LanceOtron)

**Other peak callers (e.g. MACS2):** you need a tab-separated file with the header `chrom start end score summit`,
where `summit` is the summit's genome coordinate. MACS2 narrowPeak gives the summit as an offset (column 10), so use
start + offset. Then set `min_score` and `strong_score` to that caller's score scale. The defaults (0.5 and 0.8)
are LanceOtron's 0–1 peak score and mean nothing for a MACS2 score.

## Outputs (in `outdir`)

- **`regions.bed`**: one 120 bp region per site, with the header `chr oligo_start oligo_stop gene`. This is the
  pipeline's `region_list`. A region shared by several genes is named after one of them (see below).
- **`sites.tsv`**: one row per gene, showing how the site was chosen:

| Column | Meaning |
|---|---|
| `tss`, `strand`, `tss_source` | The TSS used, and which rule picked it (`MANE_Select`, `appris_principal_1` or `longest_transcript`) |
| `score_<sample>` | The best qualifying peak score in each sample (0 = none) |
| `peak_class` | `strong` (best score ≥ 0.8), `weak` (0.5–0.8) or `none` |
| `site`, `site_source` | The chosen position, and whether it is a `consensus_summit` or the `tss` |
| `site_sample`, `site_score`, `site_support` | Which sample's summit was chosen, its score, and how many samples' strong peaks overlap it |
| `site_offset_from_tss` | The site relative to the TSS in bp (negative = upstream) |
| `region`, `region_site` | The region this gene's oligo is in (named after a gene) and its centre; empty if skipped |
| `region_genes` | Every gene in that region, comma-separated. More than one name = a shared promoter |

## How a site is chosen

1. **One TSS per gene**: the MANE Select transcript; else the longest APPRIS principal_1 transcript; else the
   longest transcript. Genes are grouped by name. Genes in the pseudoautosomal region, which appear on both chrX and
   chrY, keep the chrX copy.
2. **A peak at the TSS, per sample**: a peak counts if its score is at least `min_score`, it overlaps TSS ± `tss_window`,
   **and** its summit is within `max_summit_distance` of the TSS. If several peaks count, the highest-scoring one is
   kept.
3. **Class**: set by the best score across samples (see `peak_class` above).
4. **Site**: the summit the samples agree on. Only strong peaks vote; if a gene has no strong peak, all its peaks
   vote. A peak's support is the number of strong peaks overlapping it, itself included. The most-supported peak wins.
   Ties go to the summit closest to the TSS, then to the first sample in alphabetical order of label. A gene with no
   peak gets its TSS.
5. **One region per site**: genes are taken best first (peak class, then support, then score, then the site closest
   to the gene's own TSS, then name). A gene gets its own region unless a region already taken is closer than
   `min_site_spacing` (default 120 bp, so regions never overlap). In that case it shares that region, so it is not
   lost, and `region_genes` lists it. Example: GANC and TMEM87A start at the same promoter, so they get one oligo,
   named GANC, with `region_genes` = `GANC,TMEM87A`. Across all protein-coding genes, about 10% share a region this way.

## Settings (all optional; defaults shown)

| Key | Default | Meaning |
|---|---|---|
| `chromosomes` | chr1–22, chrX, chrY | Chromosomes to use |
| `gene_type` / `transcript_type` | `protein_coding` | Set to `null` to keep every type |
| `min_score` | 0.5 | A peak counts if its score is at least this |
| `strong_score` | 0.8 | A strong peak; only strong peaks vote on the summit |
| `tss_window` | 250 | The peak must overlap TSS ± this many bp ... |
| `max_summit_distance` | 1000 | ... and its summit must be within this many bp of the TSS |
| `region_length` | 120 | Must equal the pipeline's oligo length for single best mode |
| `min_site_spacing` | 120 (= `region_length`) | Sites closer than this share one region. Raise it to keep viewpoints further apart; 0 gives every gene its own region even where they overlap (only for genes split across separate captures) |
| `no_peak` | `tss` | Genes with no peak: `tss` puts their oligo on the TSS (e.g. inactive promoters as controls); `skip` leaves them out |

## Notes

- **Several cell types** mean that a promoter open in any of them gets its summit. If you want sites from your own
  cell type only, give only that sample.
- **Genes with no peak** get an oligo on their TSS by default. Set `no_peak: skip` if you want only genes with a peak.
  The run prints how many there are either way.
- **Check the designed oligos too.** This repo works on positions, so it cannot see sequence. Two things only show
  up after design:
  - Redesign in the oligo pipeline can move an oligo up to `pad` bp, so check the final list for overlaps.
  - An oligo whose 120 bp also occurs elsewhere in the genome (paralog families such as SPDYE, CT45A or TSPY) captures
    every copy, and its reads cannot be placed.
- **Tests:** `python -m pytest tests` checks every rule, shared sites, `no_peak` and both LanceOtron layouts on a tiny
  made-up genome.
