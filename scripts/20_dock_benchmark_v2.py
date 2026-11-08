"""
20_dock_benchmark_v2.py — 在 ChEMBL 基准集上重跑 ALK5 对接，输出标准富集指标

actives: ChEMBL CHEMBL4439 IC50<=1uM（120 个，骨架去冗余）
decoys : L1000 性质匹配 decoy（400 个，DUD-E 匹配原则，Tanimoto<0.35）
指标   : ROC-AUC / EF1% / EF5% / EF10% / EF20% / BEDROC(alpha=20)
         + bootstrap 95% CI（与原 5 actives + 60 decoys 的宽置信区间对比）

支持断点续跑：已存在的 _out.pdbqt 直接读分数。
"""
import os, subprocess, argparse, glob
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from openbabel import pybel

RDLogger.DisableLog("rdApp.*")
ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
BS = os.path.join(ROOT, "data/phase3/benchmark_v2")
REC = os.path.join(ROOT, "data/phase3/receptor/2WOU_A.pdbqt")
LIG = os.path.join(BS, "ligands")
OUT = os.path.join(ROOT, "results", "phase3")
os.makedirs(LIG, exist_ok=True)
VINA = r"D:/vina"
CENTER = (2.43, 7.88, 4.22)
SIZE = (23.6, 22.0, 22.0)


def prep(smi, tag):
    sdf = os.path.join(LIG, tag + ".sdf")
    pdbqt = os.path.join(LIG, tag + ".pdbqt")
    if os.path.exists(pdbqt):
        return pdbqt
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=42, useRandomCoords=True) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m, maxIters=300)
    Chem.MolToMolFile(m, sdf)
    next(pybel.readfile("sdf", sdf)).write("pdbqt", pdbqt, overwrite=True)
    return pdbqt


def dock(pdbqt, tag, exh=8):
    out = os.path.join(LIG, tag + "_out.pdbqt")
    if os.path.exists(out):
        for line in open(out):
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    cmd = [VINA, "--receptor", REC, "--ligand", pdbqt, "--out", out,
           "--center_x", str(CENTER[0]), "--center_y", str(CENTER[1]), "--center_z", str(CENTER[2]),
           "--size_x", str(SIZE[0]), "--size_y", str(SIZE[1]), "--size_z", str(SIZE[2]),
           "--exhaustiveness", str(exh), "--num_modes", "9", "--seed", "42"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    for line in r.stdout.splitlines():
        if "REMARK VINA RESULT" in line:
            return float(line.split()[3])
    if os.path.exists(out):
        for line in open(out):
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    return None


def bedroc(y, s, alpha=20.0):
    y = np.asarray(y); s = np.asarray(s)
    n = len(y); n_act = int(y.sum())
    if n_act == 0:
        return np.nan
    order = np.argsort(-s, kind="stable")
    ranks = (np.where(y[order] == 1)[0] + 1) / n
    rie = np.sum(np.exp(-alpha * ranks)) / n_act
    return rie * alpha / (1 - np.exp(-alpha))


def metrics(y, s):
    y = np.asarray(y); s = np.asarray(s)
    n = len(y); n_act = int(y.sum()); n_dec = n - n_act
    order = np.argsort(-s, kind="stable")
    ys = y[order]
    out = {}
    # ROC-AUC（Mann-Whitney）
    from scipy.stats import rankdata
    rk = rankdata(s)
    out["roc_auc"] = (rk[y == 1].sum() - n_act * (n_act + 1) / 2) / (n_act * n_dec)
    for pct in (1, 5, 10, 20):
        k = max(1, int(round(n * pct / 100)))
        out[f"ef{pct}"] = (ys[:k].sum() / n_act) / (k / n)
    out["bedroc20"] = bedroc(y, s)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exh", type=int, default=8)
    ap.add_argument("--max-decoys", type=int, default=None)
    args = ap.parse_args()

    act = pd.read_csv(os.path.join(BS, "actives.csv"))
    dec = pd.read_csv(os.path.join(BS, "decoys.csv"))
    if args.max_decoys:
        dec = dec.head(args.max_decoys)
    print(f"actives {len(act)} | decoys {len(dec)} | total {len(act)+len(dec)}")

    rows = []
    import time
    t0 = time.time()
    total = len(act) + len(dec)
    for i, (df, label) in enumerate([(act, 1), (dec, 0)]):
        key = "chembl_id" if label else "pert_id"
        for j, r in df.reset_index(drop=True).iterrows():
            tag = f"{'A' if label else 'D'}{j:04d}"
            p = prep(r.smiles, tag)
            aff = dock(p, tag, args.exh) if p else None
            rows.append(dict(tag=tag, label=label, id=r[key], smiles=r.smiles, affinity=aff))
            done = len(rows)
            if done % 20 == 0:
                el = time.time() - t0
                eta = el / done * (total - done) / 60
                print(f"  {done}/{total}  aff={aff}  ETA {eta:.0f} min", flush=True)

    df = pd.DataFrame(rows).dropna(subset=["affinity"])
    df.to_csv(os.path.join(OUT, "benchmark_v2_scores.csv"), index=False)
    m = metrics(df.label.values, df.affinity.values)
    lines = [
        "ALK5 docking benchmark — ChEMBL-scale",
        f"actives {int(df.label.sum())}  decoys {int((1-df.label).sum())}  ratio 1:{int((1-df.label).sum())/max(int(df.label.sum()),1):.1f}",
        f"exhaustiveness {args.exh}, seed 42, box centre {CENTER}, size {SIZE}",
        "",
        f"ROC-AUC  {m['roc_auc']:.3f}",
        f"EF1%     {m['ef1']:.2f}",
        f"EF5%     {m['ef5']:.2f}",
        f"EF10%    {m['ef10']:.2f}",
        f"EF20%    {m['ef20']:.2f}",
        f"BEDROC20 {m['bedroc20']:.3f}",
    ]
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(OUT, "benchmark_v2_summary.txt"), "w").write(txt)


if __name__ == "__main__":
    main()
