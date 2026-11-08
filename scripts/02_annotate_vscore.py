#!/usr/bin/env python3
"""
02_annotate_vscore.py — QC → 聚类(仅 HVG+scale 副本) → 标注 PT → 拆分 healthy/injured PT → v-score。

主 adata.X 保持 log1p 归一化计数（用于 marker 与 v-score）；聚类在 HVG+scale 的副本上做。

产出:
  results/gse131882_processed.h5ad
  results/vscore_injury_normal2injured.csv    （损伤签名，正向）
  results/vscore_reversal_injured2normal.csv  （逆转签名 = 取负，找保护药）
"""
import os, sys, warnings
import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
RESULTS = os.path.join(ROOT, "results")
IN = os.path.join(RESULTS, "gse131882_merged_raw.h5ad")

PT_MARKERS = ["LRP2", "SLC34A1", "ALDOB", "AQP1", "CUBN", "SLC5A2"]
INJURY_MARKERS = ["HAVCR1", "VCAM1", "LCN2", "KRT8"]

def qc(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    adata = adata[adata.obs.n_genes_by_counts > 200, :].copy()
    adata = adata[adata.obs.total_counts > 500, :].copy()
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"[QC] {adata.n_obs} cells x {adata.n_vars} genes after QC")
    print(adata.obs.groupby("sample_id").size().to_string())
    return adata

def cluster_on_hvg(adata):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    hvg = adata.var.highly_variable.values
    sub = adata[:, hvg].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=50)
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
    sc.tl.umap(sub)
    sc.tl.leiden(sub, resolution=1.0)
    adata.obs["leiden"] = sub.obs["leiden"].values
    adata.obsm["X_umap"] = sub.obsm["X_umap"].copy()
    print(f"[cluster] {sub.obs['leiden'].nunique()} clusters")
    return adata

def _gexpr(adata, g):
    """取单基因在 log 计数(.X)里的表达向量。"""
    i = adata.var_names.get_loc(g)
    col = adata.X[:, i]
    return np.asarray(col.todense()).ravel() if sparse.issparse(col) else np.asarray(col).ravel()

def annotate_pt(adata):
    pt_genes = [g for g in PT_MARKERS if g in adata.var_names]
    inj_genes = [g for g in INJURY_MARKERS if g in adata.var_names]
    print(f"[annot] PT markers: {pt_genes}")
    print(f"[annot] injury markers: {inj_genes}")
    for g in pt_genes:
        adata.obs[f"pt_{g}"] = _gexpr(adata, g)
    for g in inj_genes:
        adata.obs[f"inj_{g}"] = _gexpr(adata, g)
    adata.obs["pt_score"] = np.mean([adata.obs[f"pt_{g}"] for g in pt_genes], axis=0)
    is_pt = adata.obs["pt_score"] > np.median(adata.obs["pt_score"])
    for g in pt_genes:
        is_pt = is_pt | (adata.obs[f"pt_{g}"] > 0)
    adata.obs["is_pt"] = is_pt.astype(bool)
    print(f"[annot] PT cells: {adata.obs['is_pt'].sum()}")
    return adata

def split_pt(adata):
    inj_genes = [g for g in INJURY_MARKERS if g in adata.var_names]
    pt = adata.obs["is_pt"].values
    injured = pt & (adata.obs[[f"inj_{g}" for g in inj_genes]].max(axis=1) > 0)
    adata.obs["pt_state"] = "other"
    adata.obs.loc[pt & ~injured, "pt_state"] = "healthy_PT"
    adata.obs.loc[injured, "pt_state"] = "injured_PT"
    print("[split] healthy_PT:", (adata.obs.pt_state == "healthy_PT").sum(),
          "| injured_PT:", (adata.obs.pt_state == "injured_PT").sum())
    return adata

def compute_vscores(adata):
    from drugreflector import compute_vscores_adata
    sub = adata[adata.obs.pt_state.isin(["healthy_PT", "injured_PT"])].copy()
    vs = compute_vscores_adata(sub, group_col="pt_state",
                               group1_value="healthy_PT", group2_value="injured_PT")
    vs.to_csv(os.path.join(RESULTS, "vscore_injury_normal2injured.csv"), header=["vscore"])
    (-vs).to_csv(os.path.join(RESULTS, "vscore_reversal_injured2normal.csv"), header=["vscore"])
    print("[vscore] 损伤标志物 v-score（正向，应 >0 为损伤上调）:")
    for g in INJURY_MARKERS:
        if g in vs.index:
            print(f"    {g}: {vs[g]:.3f}")
    return vs

def main():
    if not os.path.exists(IN):
        print(f"!! 找不到 {IN}，先跑 01_load_merge.py")
        sys.exit(1)
    adata = sc.read_h5ad(IN)
    print(f"[0] 载入 {adata.n_obs} cells x {adata.n_vars} genes")
    adata.var_names = adata.var["symbol"].astype(str).tolist()
    adata.var_names_make_unique()
    adata = qc(adata)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata = cluster_on_hvg(adata)
    adata = annotate_pt(adata)
    adata = split_pt(adata)
    adata.write(os.path.join(RESULTS, "gse131882_processed.h5ad"))
    print("[done] processed h5ad saved")
    compute_vscores(adata)
    print("[done] v-score signatures saved")

if __name__ == "__main__":
    main()
