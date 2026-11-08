#!/usr/bin/env python3
"""08b_lake_disease_validate.py — Lake 2023 独立队列验证：normal PT vs disease PT（CKD/AKI）。"""
import os, numpy as np, pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr
from drugreflector import compute_vscores_adata

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
H5 = os.path.join(ROOT, "data", "lake2023_kidney_snrna.h5ad")
RES = os.path.join(ROOT, "results"); OUT = os.path.join(RES, "phase3")

ad = sc.read_h5ad(H5, backed="r")
print("[load]", ad.shape)
pt_label = "epithelial cell of proximal tubule"
mask = (ad.obs["cell_type"].astype(str) == pt_label).values
idx = np.where(mask)[0]
print("[PT cells]", len(idx))
sub = ad[idx].to_memory()
sub.obs["pt_state"] = np.where(sub.obs["disease"].astype(str) == "normal", "healthy_PT", "injured_PT")
print("[groups]", sub.obs["pt_state"].value_counts().to_dict())
print("[disease breakdown]", sub.obs["disease"].astype(str).value_counts().to_dict())

# gene symbol
if "feature_name" in sub.var.columns:
    sub.var_names = sub.var["feature_name"].astype(str).tolist()
    sub.var_names_make_unique()

vs = compute_vscores_adata(sub, group_col="pt_state", group1_value="healthy_PT", group2_value="injured_PT")
vs.to_csv(os.path.join(OUT, "vscore_injury_lake2023.csv"), header=["vscore"])
print("[vscore] saved; genes:", len(vs))

g = pd.read_csv(os.path.join(RES, "vscore_injury_normal2injured.csv"))
g.columns = ["gene", "vscore"]; g = g[~g.gene.duplicated()]
l = pd.DataFrame({"gene": vs.index, "vscore": vs.values}); l = l[~l.gene.duplicated()]
m = g.merge(l, on="gene", suffixes=("_gse", "_lake"))
print("[overlap genes]", len(m))
if len(m) > 100:
    pr = pearsonr(m.vscore_gse, m.vscore_lake); sr = spearmanr(m.vscore_gse, m.vscore_lake)
    txt = (f"GSE131882 vs Lake2023(independent cohort) injury signature concordance\n"
           f"common genes = {len(m)}\nPearson r = {pr[0]:.3f} (p = {pr[1]:.2e})\n"
           f"Spearman rho = {sr[0]:.3f} (p = {sr[1]:.2e})\n")
    print(txt); open(os.path.join(OUT, "signature_concordance.txt"), "w").write(txt)

    markers = ["HAVCR1","VCAM1","LCN2","KRT8","KRT18","SPP1","SOX9","ALDOB","LRP2","SLC34A1"]
    rows = []
    for gname in markers:
        if (m.gene == gname).any():
            r = m[m.gene == gname].iloc[0]
            rows.append((gname, round(r.vscore_gse,3), round(r.vscore_lake,3),
                         "same" if np.sign(r.vscore_gse) == np.sign(r.vscore_lake) else "opposite"))
    md = pd.DataFrame(rows, columns=["gene","vscore_GSE131882","vscore_Lake2023","direction"])
    print(md.to_string(index=False)); md.to_csv(os.path.join(OUT, "marker_concordance.csv"), index=False)

# ---- 同口径：在 Lake PT 内用同一套标志物定义 injured/healthy ----
MK = ["HAVCR1","VCAM1","LCN2","KRT8"]
mk_present = [g for g in MK if g in sub.var_names]
print("[marker-based] using", mk_present)
X = sub.X
from scipy import sparse
def col(g):
    i = sub.var_names.get_loc(g); c = sub.X[:, i]
    return np.asarray(c.todense()).ravel() if sparse.issparse(c) else np.asarray(c).ravel()
inj = np.zeros(sub.n_obs, dtype=bool)
for _mg in mk_present:
    inj |= (col(_mg) > 0)
sub.obs["pt_state_mk"] = np.where(inj, "injured_PT", "healthy_PT")
print("[marker-based groups]", sub.obs["pt_state_mk"].value_counts().to_dict())
if sub.obs["pt_state_mk"].value_counts().min() > 50:
    vs2 = compute_vscores_adata(sub, group_col="pt_state_mk", group1_value="healthy_PT", group2_value="injured_PT")
    vs2.to_csv(os.path.join(OUT, "vscore_injury_lake2023_markerdef.csv"), header=["vscore"])
    l2 = pd.DataFrame({"gene": vs2.index, "vscore": vs2.values}); l2 = l2[~l2.gene.duplicated()]
    m2 = g.merge(l2, on="gene", suffixes=("_gse","_lake"))
    if len(m2) > 100:
        pr2 = pearsonr(m2.vscore_gse, m2.vscore_lake); sr2 = spearmanr(m2.vscore_gse, m2.vscore_lake)
        t2 = (f"Marker-defined (same definition) concordance\ncommon genes = {len(m2)}\n"
              f"Pearson r = {pr2[0]:.3f} (p={pr2[1]:.2e})\nSpearman rho = {sr2[0]:.3f} (p={sr2[1]:.2e})\n")
        print(t2); open(os.path.join(OUT,"signature_concordance_markerdef.txt"),"w").write(t2)
print("[done]")
