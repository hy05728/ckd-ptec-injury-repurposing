"""
18_marker_vs_unsupervised.py — marker-based 损伤定义与无监督状态结构的一致性检验

问题：injured/healthy 用 4 个标志物定义，是简化。它离真实的状态结构有多远？
做法：
  (A) 在 PTEC 子集上做完全无监督聚类（不输入 marker 标签）
  (B) 再做一次，但**剔除 4 个 marker 基因**——排除"用标志物定义损伤、再用标志物聚类"的循环
  (C) 一致性：ARI / NMI，以及 injured 标签在各簇的富集
  (D) 敏感性：只用无监督聚类判定的 injured 簇重算 v-score -> 重跑 DrugReflector，
      看 ROCK / ALK5 抑制剂的排位是否变化

产出： results/phase3/unsup_consistency.txt / vscore_injury_unsup.csv
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.stats import pearsonr, spearmanr

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
RES = os.path.join(ROOT, "results")
P3 = os.path.join(RES, "phase3")
MARKERS = ["HAVCR1", "VCAM1", "LCN2", "KRT8"]

log = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.append(s)


def cluster(adata, tag, drop_markers=False):
    sub = adata[:, ~adata.var_names.isin(MARKERS)].copy() if drop_markers else adata.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=2000)
    sub = sub[:, sub.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(sub, n_pcs=30, n_neighbors=15)
    sc.tl.leiden(sub, resolution=0.6, key_added="cl")
    return sub


def main():
    P("=" * 66)
    P("marker-based injury definition vs unsupervised PTEC state structure")
    P("=" * 66)

    a = sc.read_h5ad(os.path.join(RES, "gse131882_processed.h5ad"))
    pt = a[a.obs.pt_state.isin(["healthy_PT", "injured_PT"])].copy()
    P(f"\nPTEC subset: {pt.n_obs} cells "
      f"({(pt.obs.pt_state=='injured_PT').sum()} injured, "
      f"{(pt.obs.pt_state=='healthy_PT').sum()} healthy)")

    results = {}
    for tag, drop in [("all genes", False), ("marker genes removed", True)]:
        P(f"\n--- clustering: {tag} ---")
        sub = cluster(pt, tag, drop_markers=drop)
        lab = sub.obs["cl"].astype(str).values
        truth = pt.obs.pt_state.values
        ari = adjusted_rand_score(truth, lab)
        nmi = normalized_mutual_info_score(truth, lab)
        ncl = len(np.unique(lab))
        P(f"  clusters: {ncl}")
        P(f"  ARI = {ari:.3f}   NMI = {nmi:.3f}")

        ct = pd.crosstab(sub.obs["cl"].astype(str), truth)
        ct["injured_frac"] = ct.get("injured_PT", 0) / ct.sum(axis=1)
        ct = ct.sort_values("injured_frac", ascending=False)
        P("  injured fraction per cluster (top 6):")
        for idx, r in ct.head(6).iterrows():
            P(f"    cluster {idx:>3s}: n={int(r.sum()-r['injured_frac']):4d}  "
              f"injured_frac={r['injured_frac']:.2f}")
        results[tag] = dict(ari=ari, nmi=nmi, n_clusters=ncl, ct=ct, sub=sub)

    # ---- 敏感性分析：用无监督簇重新定义 injured / healthy ----
    P("\n" + "=" * 66)
    P("sensitivity: redefine injured/healthy from the marker-free clustering")
    P("=" * 66)
    sub = results["marker genes removed"]["sub"]
    ct = results["marker genes removed"]["ct"]
    inj_clusters = ct[ct.injured_frac >= 0.60].index.tolist()
    hlt_clusters = ct[ct.injured_frac <= 0.20].index.tolist()
    lab = sub.obs["cl"].astype(str)
    inj_mask = lab.isin(inj_clusters).values
    hlt_mask = lab.isin(hlt_clusters).values
    P(f"  injured clusters (injured_frac >= 0.60): {inj_clusters}")
    P(f"  healthy clusters (injured_frac <= 0.20): {hlt_clusters}")
    P(f"  cells: injured {inj_mask.sum()}, healthy {hlt_mask.sum()}")
    overlap = (inj_mask & (pt.obs.pt_state.values == "injured_PT")).sum()
    P(f"  overlap with marker-based injured label: {overlap}/{inj_mask.sum()} "
      f"({100*overlap/max(inj_mask.sum(),1):.1f}%)")

    # v-score：unsup injured vs unsup healthy（稀疏矩阵安全实现）
    from scipy.sparse import issparse

    def mean_var(M):
        mu = np.asarray(M.mean(axis=0)).ravel()
        ex2 = np.asarray(M.multiply(M).mean(axis=0)).ravel() if issparse(M) \
            else np.mean(np.asarray(M) ** 2, axis=0)
        return mu, np.maximum(ex2 - mu ** 2, 1e-10)

    X = pt.X
    gi, vi = mean_var(X[inj_mask])
    gh, vh = mean_var(X[hlt_mask])
    vs = (gi - gh) / np.sqrt(vi + vh)
    names = pt.var["symbol"].astype(str).values if "symbol" in pt.var else pt.var_names.values
    out = pd.DataFrame({"gene": names, "vscore_unsup": vs}).dropna()
    out = out[out.gene.str.len() > 0].drop_duplicates("gene").set_index("gene")
    out.to_csv(os.path.join(P3, "vscore_injury_unsup.csv"))

    orig = pd.read_csv(os.path.join(RES, "vscore_injury_normal2injured.csv"))
    orig.columns = ["gene", "vscore_marker"]
    mm = orig.set_index("gene").join(out, how="inner").dropna()
    r1 = pearsonr(mm.vscore_marker, mm.vscore_unsup)[0]
    r2 = spearmanr(mm.vscore_marker, mm.vscore_unsup)[0]
    P(f"  v-score concordance (marker-based vs unsupervised): "
      f"Pearson r = {r1:.3f}, Spearman rho = {r2:.3f}  (n = {len(mm)} genes)")

    txt = "\n".join(log)
    open(os.path.join(P3, "unsup_consistency.txt"), "w", encoding="utf-8").write(txt)
    P(f"\nwritten -> {os.path.join(P3, 'unsup_consistency.txt')}")


if __name__ == "__main__":
    main()
