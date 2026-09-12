#!/usr/bin/env python3
"""Build a protein_id -> gene -> hg38 position map for cis/trans pQTL classification.

Input : meta-feature.txt (SomaLogic feature table) with at least columns
        SeqId (e.g. 10000-28, matching the pipeline protein_id) and UniProt.
Output: protein_gene_map.tsv with columns
        protein_id  uniprot  gene_symbol  chr  gene_start  gene_end

Gene coordinates come from MyGene.info (genomic_pos, GRCh38). This must be run on
a node with internet access (login node), once; the resulting static map is then
consumed by the HERITABILITY pipeline module (params.protein_gene_map).

cis/trans (per UKB-PPP, Sun et al. 2023) is defined as a pQTL within 1 Mb of the
gene encoding the protein, so we keep the gene's start and end on a standard
chromosome and let the module apply the +/-1 Mb window.
"""
import argparse
import csv
import sys
import time
from pathlib import Path

import requests

MYGENE_URL = "https://mygene.info/v3/query"


def query_batch(uniprot_ids):
    ids = [u for u in uniprot_ids if u and u.strip()]
    if not ids:
        return []
    q = "uniprot:" + " OR uniprot:".join(ids)
    params = {"q": q, "fields": "symbol,genomic_pos,uniprot", "species": "human", "size": 1000}
    for attempt in range(3):
        try:
            r = requests.get(MYGENE_URL, params=params, timeout=60)
            r.raise_for_status()
            return r.json().get("hits", [])
        except Exception as e:  # noqa: BLE001
            print(f"  MyGene query failed (attempt {attempt+1}): {e}", file=sys.stderr)
            time.sleep(2)
    return []


def pos_on_standard_chr(genomic_pos):
    """Pick the genomic_pos entry on a numbered/standard chromosome (skip HSCHR alt loci)."""
    if isinstance(genomic_pos, list) and genomic_pos:
        return next(
            (p for p in genomic_pos if p.get("chr") and not str(p.get("chr")).startswith("HSCHR")),
            genomic_pos[0],
        )
    return genomic_pos if isinstance(genomic_pos, dict) else None


def hits_to_uniprot(hits):
    out = {}
    for hit in hits:
        u = hit.get("uniprot")
        keys = []
        if isinstance(u, str):
            keys = [u]
        elif isinstance(u, list):
            keys = [x for x in u if isinstance(x, str)]
        elif isinstance(u, dict):
            for v in u.values():
                keys += [v] if isinstance(v, str) else (v if isinstance(v, list) else [])
        for k in keys:
            out.setdefault(k, hit)  # first hit wins
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--feature-file", type=Path, required=True, help="meta-feature.txt (SeqId, UniProt columns)")
    ap.add_argument("--out", type=Path, required=True, help="output protein_gene_map.tsv")
    ap.add_argument("--seqid-col", default="SeqId")
    ap.add_argument("--uniprot-col", default="UniProt")
    ap.add_argument("--symbol-col", default="EntrezGeneSymbol")
    args = ap.parse_args()

    with args.feature_file.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"ERROR: {args.feature_file} is empty")

    # protein_id == SeqId; collect unique UniProt IDs to query.
    proteins = []
    uniprots = set()
    for r in rows:
        seqid = (r.get(args.seqid_col) or "").strip()
        uni = (r.get(args.uniprot_col) or "").strip().split(",")[0].split(" ")[0]
        sym = (r.get(args.symbol_col) or "").strip()
        if not seqid:
            continue
        proteins.append((seqid, uni, sym))
        if uni:
            uniprots.add(uni)

    print(f">>> {len(proteins)} features, {len(uniprots)} unique UniProt IDs to map via MyGene.info")
    mapping = {}
    uids = list(uniprots)
    for i in range(0, len(uids), 100):
        mapping.update(hits_to_uniprot(query_batch(uids[i : i + 100])))
        time.sleep(0.34)
        print(f"    mapped {min(i+100, len(uids))}/{len(uids)}")

    n_ok = 0
    with args.out.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["protein_id", "uniprot", "gene_symbol", "chr", "gene_start", "gene_end"])
        for seqid, uni, sym in proteins:
            chrom = start = end = ""
            hit = mapping.get(uni)
            if hit:
                if not sym:
                    sym = hit.get("symbol", "") or ""
                p = pos_on_standard_chr(hit.get("genomic_pos"))
                if isinstance(p, dict):
                    chrom = str(p.get("chr", "") or "")
                    start = str(p.get("start", "") or "")
                    end = str(p.get("end", "") or "")
            if chrom and start and end:
                n_ok += 1
            w.writerow([seqid, uni, sym, chrom, start, end])

    print(f">>> wrote {args.out}: {n_ok}/{len(proteins)} features with hg38 coordinates")


if __name__ == "__main__":
    main()
