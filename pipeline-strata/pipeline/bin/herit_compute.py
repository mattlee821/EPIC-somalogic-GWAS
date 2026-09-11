#!/usr/bin/env python3
"""Per-protein pQTL variance-explained component for the heritability decomposition.

Given the independent index SNPs from clumping and their meta-analysis stats,
compute the variance explained by each (2*f*(1-f)*beta^2, valid because the
phenotype is inverse-rank normalised so beta is in SD units), classify each as
cis or trans relative to the gene encoding the protein (within +/-window of the
gene body, per UKB-PPP), and emit:
  * a summary of counts and summed cis/trans variance (the "pQTL component"); and
  * the genomic intervals to exclude (index SNP +/- window) so the LDSC
    "polygenic component" is estimated on the non-pQTL genome.
"""
import argparse
import sys
from pathlib import Path


def norm_chr(c):
    c = str(c).strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assoc", type=Path, required=True,
                    help="deduped assoc table: SNP CHR POS EA OA EAF BETA P (header)")
    ap.add_argument("--index", type=Path, required=True, help="index SNP IDs (chr:pos), one per line")
    ap.add_argument("--protein-map", type=Path, required=True)
    ap.add_argument("--protein-id", required=True)
    ap.add_argument("--cis-window-kb", type=float, default=1000.0)
    ap.add_argument("--out-summary", type=Path, required=True)
    ap.add_argument("--out-exclude", type=Path, required=True)
    args = ap.parse_args()

    window = int(args.cis_window_kb * 1000)

    # gene location for this protein (may be absent -> all index SNPs are 'trans'/unknown)
    gene_chr = gene_start = gene_end = None
    with args.protein_map.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[idx["protein_id"]] == args.protein_id:
                try:
                    gene_chr = norm_chr(f[idx["chr"]])
                    gene_start = int(float(f[idx["gene_start"]]))
                    gene_end = int(float(f[idx["gene_end"]]))
                except (ValueError, KeyError, IndexError):
                    gene_chr = gene_start = gene_end = None
                break

    # assoc stats keyed by SNP id (chr:pos)
    stats = {}
    with args.assoc.open() as fh:
        h = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                stats[f[col["SNP"]]] = (
                    norm_chr(f[col["CHR"]]), int(f[col["POS"]]),
                    float(f[col["EAF"]]), float(f[col["BETA"]]),
                )
            except (ValueError, IndexError):
                continue

    index_ids = [ln.strip() for ln in args.index.read_text().splitlines() if ln.strip()]

    n_cis = n_trans = 0
    h2_cis = h2_trans = 0.0
    excludes = []
    for sid in index_ids:
        s = stats.get(sid)
        if s is None:
            continue
        chrom, pos, eaf, beta = s
        if not (0.0 < eaf < 1.0):
            continue
        var = 2.0 * eaf * (1.0 - eaf) * beta * beta
        is_cis = (
            gene_chr is not None
            and chrom == gene_chr
            and pos >= (gene_start - window)
            and pos <= (gene_end + window)
        )
        if is_cis:
            n_cis += 1
            h2_cis += var
        else:
            n_trans += 1
            h2_trans += var
        excludes.append((chrom, max(0, pos - window), pos + window))

    with args.out_summary.open("w") as out:
        out.write(f"n_index={n_cis + n_trans}\n")
        out.write(f"n_cis={n_cis}\n")
        out.write(f"n_trans={n_trans}\n")
        out.write(f"h2_cis={h2_cis:.6g}\n")
        out.write(f"h2_trans={h2_trans:.6g}\n")
        out.write(f"gene_chr={gene_chr if gene_chr is not None else 'NA'}\n")

    with args.out_exclude.open("w") as out:
        for chrom, a, b in excludes:
            out.write(f"{chrom}\t{a}\t{b}\n")

    print(f"{args.protein_id}: index={n_cis + n_trans} cis={n_cis} trans={n_trans} "
          f"h2_cis={h2_cis:.4g} h2_trans={h2_trans:.4g}", file=sys.stderr)


if __name__ == "__main__":
    main()
