#!/usr/bin/env python3
"""
01_load_merge.py — 加载 R 转换后的样本目录，用 Ensembl ID 作唯一索引合并 6 样本。

输入: data/processed/GSMxxxx_<group>.sN/{matrix.mtx, genes.tsv, barcodes.tsv}
产出: results/gse131882_merged_raw.h5ad   (var_names=Ensembl ID, var['symbol']=gene symbol)
"""
import os, sys, glob, warnings
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix

warnings.filterwarnings("ignore")
ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
PROC = os.path.join(ROOT, "data", "processed")
RESULTS = os.path.join(ROOT, "results")

def load_sample(d):
    mtx = mmread(os.path.join(d, "matrix.mtx"))          # genes x cells
    genes = pd.read_csv(os.path.join(d, "genes.tsv"), sep="\t", header=None, names=["ens", "symbol"])
    barcodes = pd.read_csv(os.path.join(d, "barcodes.tsv"), header=None)[0].astype(str).tolist()
    a = sc.AnnData(X=csr_matrix(mtx).T)                  # cells x genes
    a.var_names = genes["ens"].astype(str).tolist()      # Ensembl 作唯一索引
    a.var["symbol"] = genes["symbol"].astype(str).tolist()
    a.obs_names = barcodes
    return a

def main():
    sample_dirs = sorted(glob.glob(os.path.join(PROC, "GSM*")))
    if not sample_dirs:
        print(f"!! 未找到样本目录于 {PROC}，先跑 convert_rds.R")
        sys.exit(1)
    print(f"[0] 找到 {len(sample_dirs)} 个样本目录")
    ads = []
    for d in sample_dirs:
        sid = os.path.basename(d)
        a = load_sample(d)
        group = "control" if "control" in sid else ("diabetes" if "diabetes" in sid else "unknown")
        a.obs["sample_id"] = sid
        a.obs["group"] = group
        a.obs["replicate"] = sid.split(".")[-1]
        print(f"    {sid}: {a.n_obs} cells x {a.n_vars} genes (group={group})")
        ads.append(a)
    adata = anndata.concat(ads, join="outer", index_unique="-")
    adata.obs_names_make_unique()
    # concat 会丢 var 列，重建 Ensembl -> symbol 映射并回填
    ens2sym = {}
    for d in sample_dirs:
        g = pd.read_csv(os.path.join(d, "genes.tsv"), sep="\t", header=None, names=["ens", "symbol"])
        for ens, sym in zip(g["ens"].astype(str), g["symbol"].astype(str)):
            ens2sym[ens] = sym
    adata.var["symbol"] = [ens2sym.get(e, e) for e in adata.var_names]
    os.makedirs(RESULTS, exist_ok=True)
    out = os.path.join(RESULTS, "gse131882_merged_raw.h5ad")
    adata.write(out)
    print(f"[done] {adata.n_obs} cells x {adata.n_vars} genes -> {out}")
    print(adata.obs.groupby(["group", "sample_id"]).size().to_string())

if __name__ == "__main__":
    main()
