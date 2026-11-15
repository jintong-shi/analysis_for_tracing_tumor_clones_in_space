#!/usr/bin/env python3
# CalicoST input prep: BAF AD/DP matrices -> A/B allele npz + SNP IDs.
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def parse_args():
    ap = argparse.ArgumentParser(
        description="Convert BAF AD/DP matrices into CalicoST synthetic snpinfo (A/B allele npz + SNP IDs)."
    )
    ap.add_argument("--baf_dir", required=True, help="Folder containing AD.mtx, DP.mtx, barcodes.tsv, snps.tsv")
    ap.add_argument("--output_dir", required=True, help="Output folder for CalicoST snpinfo files")
    ap.add_argument("--ad_mtx", default=None, help="Path to AD.mtx (optional, overrides baf_dir)")
    ap.add_argument("--dp_mtx", default=None, help="Path to DP.mtx (optional, overrides baf_dir)")
    ap.add_argument("--barcodes", default=None, help="Path to barcodes.tsv (optional, overrides baf_dir)")
    ap.add_argument("--snps", default=None, help="Path to snps.tsv (optional, overrides baf_dir)")
    ap.add_argument("--snps_no_header", action="store_true", help="snps.tsv has no header")
    ap.add_argument("--keep_chr_prefix", action="store_true", help="Keep 'chr' prefix in unique_snp_ids")
    ap.add_argument("--max_snps", type=int, default=None, help="Downsample to at most this many SNPs")
    ap.add_argument("--snp_fraction", type=float, default=None, help="Downsample to this fraction of SNPs (0-1)")
    ap.add_argument("--max_cells_per_snp", type=int, default=None,
                    help="Cap number of cells with coverage per SNP (keeps all SNPs)")
    ap.add_argument("--seed", type=int, default=13, help="Random seed for SNP downsampling")
    return ap.parse_args()


def resolve_inputs(args):
    args.ad_mtx = args.ad_mtx or os.path.join(args.baf_dir, "AD.mtx")
    args.dp_mtx = args.dp_mtx or os.path.join(args.baf_dir, "DP.mtx")
    args.barcodes = args.barcodes or os.path.join(args.baf_dir, "barcodes.tsv")
    args.snps = args.snps or os.path.join(args.baf_dir, "snps.tsv")

    for path in (args.ad_mtx, args.dp_mtx, args.barcodes, args.snps):
        if not os.path.exists(path):
            raise SystemExit(f"Missing input: {path}")


def load_snps(snps_path, snps_no_header=False, keep_chr_prefix=False):
    df = pd.read_csv(snps_path, sep="\t", header=None if snps_no_header else "infer")

    # Normalize columns
    if df.shape[1] == 1:
        # Expect format chr_pos_ref_alt
        parts = df.iloc[:, 0].astype(str).str.split("_", expand=True)
        if parts.shape[1] < 4:
            raise SystemExit("snps.tsv appears to have a single column, but values are not chr_pos_ref_alt.")
        df = pd.DataFrame({
            "chr": parts[0],
            "pos": parts[1],
            "ref": parts[2],
            "alt": parts[3],
        })
    else:
        df.columns = [str(c).lower() for c in df.columns]
        col_chr = "chr" if "chr" in df.columns else df.columns[0]
        col_pos = "pos" if "pos" in df.columns else df.columns[1]
        col_ref = next((c for c in df.columns if "ref" in c), df.columns[2])
        col_alt = next((c for c in df.columns if "alt" in c), df.columns[3])
        df = df[[col_chr, col_pos, col_ref, col_alt]]
        df.columns = ["chr", "pos", "ref", "alt"]

    # Normalize chr
    def normalize_chr(x):
        s = str(x)
        if s.startswith("chr") and not keep_chr_prefix:
            return s[3:]
        if (not s.startswith("chr")) and keep_chr_prefix:
            return f"chr{s}"
        return s

    df["chr"] = df["chr"].map(normalize_chr)
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)

    unique_ids = (df["chr"].astype(str) + "_" +
                  df["pos"].astype(str) + "_" +
                  df["ref"] + "_" +
                  df["alt"]).to_numpy()

    return unique_ids, df


def load_mtx(path):
    mat = scipy.io.mmread(path)
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    return mat.tocsr()


def orient_matrices(ad, dp, n_cells, n_snps):
    # Ensure matrices are same shape
    if ad.shape != dp.shape:
        raise SystemExit(f"AD/DP shape mismatch: {ad.shape} vs {dp.shape}")

    if ad.shape == (n_snps, n_cells):
        ad = ad.transpose().tocsr()
        dp = dp.transpose().tocsr()
    elif ad.shape == (n_cells, n_snps):
        ad = ad.tocsr()
        dp = dp.tocsr()
    else:
        raise SystemExit(
            f"AD/DP shape {ad.shape} does not match (n_cells, n_snps)=({n_cells}, {n_snps}) "
            f"or (n_snps, n_cells)=({n_snps}, {n_cells})."
        )
    return ad, dp


def main():
    args = parse_args()
    resolve_inputs(args)

    os.makedirs(args.output_dir, exist_ok=True)

    # Load barcodes and snps
    barcodes = pd.read_csv(args.barcodes, header=None)[0].astype(str).tolist()
    unique_ids, snp_df = load_snps(args.snps, snps_no_header=args.snps_no_header,
                                   keep_chr_prefix=args.keep_chr_prefix)

    n_cells = len(barcodes)
    n_snps = len(unique_ids)

    # Load matrices
    ad = load_mtx(args.ad_mtx)
    dp = load_mtx(args.dp_mtx)

    ad, dp = orient_matrices(ad, dp, n_cells, n_snps)

    # Optional SNP downsampling (random)
    if args.max_snps is not None or args.snp_fraction is not None:
        if args.snp_fraction is not None and not (0 < args.snp_fraction <= 1):
            raise SystemExit("--snp_fraction must be in (0, 1].")
        rng = np.random.default_rng(args.seed)
        total = ad.shape[1]
        k = args.max_snps if args.max_snps is not None else int(total * args.snp_fraction)
        k = max(1, min(k, total))
        idx = rng.choice(total, size=k, replace=False)
        idx.sort()
        ad = ad[:, idx]
        dp = dp[:, idx]
        unique_ids = unique_ids[idx]

    # Optional per-SNP thinning (keep all SNPs, reduce density)
    if args.max_cells_per_snp is not None and args.max_cells_per_snp > 0:
        rng = np.random.default_rng(args.seed)
        dp = dp.tocsc()
        ad = ad.tocsc()
        data_dp = []
        indices_dp = []
        indptr_dp = [0]
        data_ad = []
        indices_ad = []
        indptr_ad = [0]
        thinned = 0
        for col in range(dp.shape[1]):
            dp_start, dp_end = dp.indptr[col], dp.indptr[col + 1]
            dp_idx = dp.indices[dp_start:dp_end]
            dp_dat = dp.data[dp_start:dp_end]
            if dp_idx.size > args.max_cells_per_snp:
                keep = rng.choice(dp_idx.size, size=args.max_cells_per_snp, replace=False)
                keep.sort()
                dp_idx = dp_idx[keep]
                dp_dat = dp_dat[keep]
                thinned += (dp_end - dp_start) - args.max_cells_per_snp
            # keep AD only where DP is kept
            ad_start, ad_end = ad.indptr[col], ad.indptr[col + 1]
            ad_idx = ad.indices[ad_start:ad_end]
            ad_dat = ad.data[ad_start:ad_end]
            if ad_idx.size:
                mask = np.in1d(ad_idx, dp_idx, assume_unique=True)
                ad_idx = ad_idx[mask]
                ad_dat = ad_dat[mask]
            indices_dp.extend(dp_idx.tolist())
            data_dp.extend(dp_dat.tolist())
            indptr_dp.append(len(indices_dp))
            indices_ad.extend(ad_idx.tolist())
            data_ad.extend(ad_dat.tolist())
            indptr_ad.append(len(indices_ad))
        dp = sp.csc_matrix((np.array(data_dp), np.array(indices_dp), np.array(indptr_dp)), shape=dp.shape)
        ad = sp.csc_matrix((np.array(data_ad), np.array(indices_ad), np.array(indptr_ad)), shape=ad.shape)
        dp = dp.tocsr()
        ad = ad.tocsr()

    # Clamp AD to DP and compute REF
    ad = ad.minimum(dp)
    ref = dp - ad
    ref.data[ref.data < 0] = 0
    ref.eliminate_zeros()

    # Save outputs in CalicoST synthetic format
    sp.save_npz(os.path.join(args.output_dir, "cell_snp_Aallele.npz"), ref, compressed=True)
    sp.save_npz(os.path.join(args.output_dir, "cell_snp_Ballele.npz"), ad, compressed=True)
    np.save(os.path.join(args.output_dir, "unique_snp_ids.npy"), unique_ids)

    with open(os.path.join(args.output_dir, "barcodes.txt"), "w") as f:
        for bc in barcodes:
            f.write(f"{bc}\n")

    print("[SUCCESS] Wrote CalicoST synthetic snpinfo to:", args.output_dir)
    print(" - cell_snp_Aallele.npz")
    print(" - cell_snp_Ballele.npz")
    print(" - unique_snp_ids.npy")
    print(" - barcodes.txt")
    if args.max_cells_per_snp is not None and args.max_cells_per_snp > 0:
        print(f"[INFO] Thinned {thinned} cell-SNP entries to max_cells_per_snp={args.max_cells_per_snp}")


if __name__ == "__main__":
    main()
