#!/usr/bin/env python3
"""
05_benchmark.py — ALK5 对接富集基准：已知 ALK5 抑制剂(actives) vs 随机 LINCS decoys。

- actives: 5 个已知 ALK5/TGFBR1 抑制剂
- decoys: 从 LINCS pert_info 随机抽 N 个类药分子（MW 250-500，有效 SMILES，非上述已知）
- 复用 04 的受体/盒子/对接流程
产出: results/phase3/benchmark_decoys.csv, results/phase3/enrichment_summary.txt
"""
import os, sys, random
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
sys.path.insert(0, os.path.dirname(__file__))
import importlib.util
spec = importlib.util.spec_from_file_location("dock", os.path.join(os.path.dirname(__file__), "04_dock_alk5.py"))
dock = importlib.util.module_from_spec(spec); spec.loader.exec_module(dock)

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
OUT = os.path.join(ROOT, "results", "phase3")
PERT = os.path.join(ROOT, "data", "GSE92742_pert_info.txt.gz")
ACIVES = {"LY-364947", "SB-431542", "galunisertib", "SB-505124", "RepSox"}

def sample_decoys(n=60, seed=42):
    df = pd.read_csv(PERT, compression="gzip", sep="\t", dtype=str)
    df = df[df["canonical_smiles"].notna() & (df["canonical_smiles"] != "-666")]
    df = df.drop_duplicates(subset="canonical_smiles")
    out = []
    for _, r in df.iterrows():
        smi = r["canonical_smiles"]
        if r["pert_iname"] in ACIVES or r["pert_iname"] in dock.LIGANDS:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None: continue
        mw = Descriptors.MolWt(m)
        if 250 <= mw <= 500:
            out.append((r["pert_id"] + "|" + str(r["pert_iname"]), smi, mw))
    random.seed(seed); random.shuffle(out)
    return out[:n]

def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--n", type=int, default=60); ap.add_argument("--exh", type=int, default=8)
    a = ap.parse_args()
    center, size = dock.split_and_prep_receptor()
    decoys = sample_decoys(a.n)
    print(f"[decoys] sampled {len(decoys)}")
    rows = []
    for i, (name, smi, mw) in enumerate(decoys):
        p = dock.prep_ligand("decoy_" + name.replace("|", "_"), smi)
        if p is None: continue
        aff, rc = dock.dock("decoy_" + name.replace("|", "_"), p, center, size, exh=a.exh)
        rows.append({"name": name, "smiles": smi, "mw": round(mw, 1), "affinity": aff})
        if (i + 1) % 10 == 0: print(f"  {i+1}/{len(decoys)} done")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "benchmark_decoys.csv"), index=False)

    # 合并 actives + decoys 算富集
    act = pd.read_csv(os.path.join(OUT, "docking_results.csv"))
    act = act[act["name"].isin(ACIVES)][["name", "affinity"]].copy(); act["label"] = 1
    dec = pd.DataFrame(rows)[["name", "affinity"]].copy(); dec["label"] = 0
    allc = pd.concat([act, dec]).dropna().sort_values("affinity").reset_index(drop=True)
    allc["rank"] = np.arange(1, len(allc) + 1)
    n_act = int(allc["label"].sum()); n_tot = len(allc)
    lines = [f"total compounds: {n_tot}, actives: {n_act}"]
    for pct in [5, 10, 20]:
        top = int(np.ceil(n_tot * pct / 100))
        hit = int(allc.iloc[:top]["label"].sum())
        ef = (hit / n_act) / (pct / 100) if n_act else float("nan")
        lines.append(f"EF{pct}% = {ef:.2f}  (actives in top {pct}%: {hit}/{n_act})")
    # ROC-AUC（用 -affinity 作分数）
    from sklearn.metrics import roc_auc_score
    try:
        auc = roc_auc_score(allc["label"], -allc["affinity"])
        lines.append(f"ROC-AUC = {auc:.3f}")
    except Exception as e:
        lines.append(f"ROC-AUC failed: {e}")
    # actives 各自排名
    lines.append("\nactives ranks:")
    for _, r in allc[allc["label"] == 1].iterrows():
        lines.append(f"  {r['name']}: rank {int(r['rank'])}/{n_tot}  aff={r['affinity']}")
    txt = "\n".join(lines)
    open(os.path.join(OUT, "enrichment_summary.txt"), "w").write(txt)
    print(txt)

if __name__ == "__main__":
    main()
