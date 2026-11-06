"""Each rule on a tiny made-up genome. Run: python -m pytest tests  (or: python tests/test_sitefinder.py)"""
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sitefinder.cli import run  # noqa: E402
from sitefinder.sites import load_peaks  # noqa: E402


def tx(chrom, start, end, strand, gene, tid, tags=(), gene_type="protein_coding"):
    attrs = (f'gene_id "{gene}_id"; transcript_id "{tid}"; gene_type "{gene_type}"; gene_name "{gene}"; '
             f'transcript_type "{gene_type}";' + "".join(f' tag "{t}";' for t in tags))
    return "\t".join([chrom, "TEST", "transcript", str(start), str(end), ".", strand, ".", attrs])


GTF = [
    tx("chr1", 5001, 30000, "+", "A", "A2"),                                    # longer, but not MANE
    tx("chr1", 10001, 20000, "+", "A", "A1", ["MANE_Select"]),                  # A: TSS 10000
    tx("chr1", 40001, 90000, "-", "B", "B3"),                                   # longest, but not APPRIS
    tx("chr1", 50001, 60000, "-", "B", "B1", ["appris_principal_1"]),
    tx("chr1", 52001, 70000, "-", "B", "B2", ["appris_principal_1"]),           # B: longest APPRIS, TSS 69999
    tx("chr2", 101, 2000, "+", "C", "C1"),
    tx("chr2", 1001, 5000, "+", "C", "C2"),                                     # C: longest, TSS 1000
    tx("chrY", 3001, 4000, "+", "D", "DY", ["MANE_Select"]),                    # D: PAR, chrX copy kept
    tx("chrX", 3001, 4000, "+", "D", "DX", ["MANE_Select"]),
    tx("chr2", 20001, 30000, "+", "E", "E1", gene_type="lncRNA"),               # E: not protein-coding
    tx("chr3", 50001, 60000, "+", "F", "F1", ["MANE_Select"]),                  # F: TSS 50000
    tx("chr4", 1001, 8000, "+", "G", "G1", ["MANE_Select"]),                    # G: TSS 1000  } bidirectional pair,
    tx("chr4", 101, 901, "-", "H", "H1", ["MANE_Select"]),                      # H: TSS 900   } one peak between them
    tx("chr4", 5001, 9000, "+", "J", "J1", ["MANE_Select"]),                    # J: TSS 5000  } nearby promoters with
    tx("chr4", 5401, 9000, "+", "K", "K1", ["MANE_Select"]),                    # K: TSS 5400  } their own peaks
]
PEAKS = {
    "s1": [("chr1", 9900, 10200, 0.95, 10030),     # A: strong, summit +30
           ("chr1", 69900, 70100, 0.60, 70050),     # B: weak only
           ("chr2", 500, 3000, 0.99, 2500),         # C: overlaps the TSS window but summit 1,500 bp away -> ignored
           ("chr3", 49950, 50060, 0.90, 50005),     # F: strong, summit +5, overlaps nothing else
           ("chr4", 850, 1050, 0.95, 960),          # G and H: the same summit (40 bp from G's TSS, 60 from H's)
           ("chr4", 4900, 5100, 0.60, 5010),        # J: weak, summit 5010
           ("chr4", 5350, 5500, 0.90, 5420)],       # K: strong, summit 5420 (410 bp from J's)
    "s2": [("chr1", 9850, 10100, 0.90, 9990),       # A: strong, summit -10 -> closest of two equally supported
           ("chrX", 2900, 3100, 0.40, 3000),        # D: below the 0.5 minimum -> no peak
           ("chr3", 50070, 50300, 0.85, 50100)],    # F: strong, overlaps s3's peak -> support 2
    "s3": [("chr3", 50080, 50400, 0.88, 50120)],    # F: strong, overlaps s2's peak -> support 2
}


def write_peaks(rows, path, layout="plain"):
    """The same peaks in the plain layout or in either LanceOtron layout (callPeaks, or the newer one)."""
    d = pd.DataFrame(rows, columns=["chrom", "start", "end", "score", "summit"])
    if layout == "lanceotron_callpeaks":
        d = pd.DataFrame({"chrom": d.chrom, "start": d.start, "end": d.end, "overall_peak_score": d.score,
                          "shape_score": 0.5, "enrichment_score": 0.5, "max_coverage_position": d.summit,
                          "region_id": range(len(d)), "pvalue_chrom": 100.0})
    elif layout == "lanceotron_new":
        d = pd.DataFrame({"Chromosome": d.chrom, "Start": d.start, "End": d.end, "Region_id": range(len(d)),
                          "Summit": d.summit, "Overall_peak_score": d.score, "Shape_score": 0.5})
    d.to_csv(path, sep="\t", index=False)


def build(tmp, layouts=None, **settings):
    tmp = Path(tmp)
    (tmp / "test.gtf").write_text("\n".join(GTF) + "\n")
    peaks = {}
    for s, rows in PEAKS.items():
        write_peaks(rows, tmp / f"{s}.bed", (layouts or {}).get(s, "plain"))
        peaks[s] = str(tmp / f"{s}.bed")
    return run(dict(gtf=str(tmp / "test.gtf"), peaks=peaks, outdir=str(tmp / "out"), **settings)).set_index("gene"), tmp / "out"


def test_rules():
    with tempfile.TemporaryDirectory() as tmp:
        s, out = build(tmp)
        assert set(s.index) == {"A", "B", "C", "D", "F", "G", "H", "J", "K"}           # lncRNA E excluded
        assert tuple(s.loc["A", ["tss", "tss_source"]]) == (10000, "MANE_Select")
        assert tuple(s.loc["B", ["tss", "tss_source"]]) == (69999, "appris_principal_1")
        assert tuple(s.loc["C", ["tss", "tss_source"]]) == (1000, "longest_transcript")
        assert s.loc["D", "chr"] == "chrX"
        # A: two strong peaks overlap (support 2 each) -> the summit closest to the TSS
        assert tuple(s.loc["A", ["peak_class", "site", "site_sample", "site_support"]]) == ("strong", 9990, "s2", 2)
        assert s.loc["A", "site_offset_from_tss"] == -10
        # B: weak peak only; minus strand, so a summit at a higher coordinate is upstream
        assert tuple(s.loc["B", ["peak_class", "site", "site_offset_from_tss"]]) == ("weak", 70050, -51)
        # C: summit too far; D: score too low -> no peak, site on the TSS
        assert tuple(s.loc["C", ["peak_class", "site", "site_source"]]) == ("none", 1000, "tss")
        assert tuple(s.loc["D", ["peak_class", "site"]]) == ("none", 3000)
        # F: support beats distance: s1 is closest (+5) but alone; s2 and s3 agree, s2 is the closer of those
        assert tuple(s.loc["F", ["site", "site_sample", "site_support"]]) == (50100, "s2", 2)
        r = pd.read_csv(out / "regions.bed", sep="\t")
        assert list(r.columns) == ["chr", "oligo_start", "oligo_stop", "gene"]
        assert ((r.oligo_stop - r.oligo_start) == 120).all()
        assert r.set_index("gene").loc["A", "oligo_start"] == 9930


def regions_by_gene(out):
    return pd.read_csv(out / "regions.bed", sep="\t").set_index("gene")


def test_shared_sites_get_one_region_and_keep_every_gene():
    with tempfile.TemporaryDirectory() as tmp:
        s, out = build(tmp)                                    # default spacing = 120: regions never overlap
        r = regions_by_gene(out)
        # G and H share one summit: ONE region, named after G (summit closer to its TSS), and H is kept, not lost
        assert "G" in r.index and "H" not in r.index
        assert tuple(s.loc["H", ["site", "region", "region_site", "region_genes"]]) == (960, "G", 960, "G,H")
        assert s.loc["G", "region_genes"] == "G,H"
        # J and K are 410 bp apart: separate regions at the default spacing
        assert {"J", "K"} <= set(r.index) and s.loc["J", "region_genes"] == "J"
        for ch, g in r.groupby("chr"):                         # no two regions overlap
            assert (g.oligo_start.sort_values().diff().dropna() >= 120).all(), ch
        assert len(r) == len(s) - 1                            # every gene has a region but H, which shares G's


def test_wider_spacing_merges_nearby_promoters_into_the_stronger_one():
    with tempfile.TemporaryDirectory() as tmp:
        s, out = build(tmp, min_site_spacing=500)
        assert tuple(s.loc["J", ["region", "region_site", "region_genes"]]) == ("K", 5420, "J,K")   # K strong beats J weak
        assert "J" not in regions_by_gene(out).index


def test_spacing_zero_gives_every_gene_its_own_region():
    with tempfile.TemporaryDirectory() as tmp:
        s, out = build(tmp, min_site_spacing=0)
        r = regions_by_gene(out)
        assert set(r.index) == set(s.index) and (s.region == s.index).all()
        assert r.loc["G", "oligo_start"] == r.loc["H", "oligo_start"]                 # the duplicate is back


def test_no_peak_skip_leaves_out_genes_without_a_peak():
    with tempfile.TemporaryDirectory() as tmp:
        s, out = build(tmp, no_peak="skip")
        none = set(s.index[s.peak_class == "none"])
        assert none == {"C", "D"}
        assert not (none & set(regions_by_gene(out).index)) and (s.loc[list(none), "region"] == "").all()
    with tempfile.TemporaryDirectory() as tmp:
        s, out = build(tmp)                                    # default no_peak: tss keeps them, on the TSS
        assert regions_by_gene(out).loc["C", "oligo_start"] == 1000 - 60


def test_lanceotron_files_as_they_come():
    """LanceOtron output in either layout gives exactly the sites of the plain five-column files."""
    with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
        plain, _ = build(a)
        lance, _ = build(b, {"s1": "lanceotron_callpeaks", "s2": "lanceotron_new", "s3": "lanceotron_callpeaks"})
        pd.testing.assert_frame_equal(plain, lance)


def test_missing_summit_is_an_error():
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / "no_summit.bed"
        pd.DataFrame({"chrom": ["chr1"], "start": [1], "end": [100], "score": [0.9]}).to_csv(f, sep="\t", index=False)
        try:
            load_peaks(str(f), 0.5, ["chr1"])
        except ValueError as e:
            assert "summit" in str(e)
        else:
            raise AssertionError("a peak file with no summit column was accepted")


if __name__ == "__main__":
    test_rules()
    test_shared_sites_get_one_region_and_keep_every_gene()
    test_wider_spacing_merges_nearby_promoters_into_the_stronger_one()
    test_spacing_zero_gives_every_gene_its_own_region()
    test_no_peak_skip_leaves_out_genes_without_a_peak()
    test_lanceotron_files_as_they_come()
    test_missing_summit_is_an_error()
    print("all site finder tests passed")
