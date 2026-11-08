#!/usr/bin/env python3
"""
04_dock_alk5.py — ALK5 (TGFBR1) 分子对接：受体准备 + 配体准备 + Vina 对接。

受体: PDB 2WOU chain A（去掉 HOH/EDO/ZZF），盒子取自共晶配体 ZZF 重心。
配体: SMILES -> RDKit 3D -> PDBQT。
产出: results/phase3/docking_results.csv
"""
import os, sys, subprocess, re, json, argparse
import numpy as np
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
P3 = os.path.join(ROOT, "data", "phase3")
REC_DIR = os.path.join(P3, "receptor")
LIG_DIR = os.path.join(P3, "ligands")
OUT_DIR = os.path.join(ROOT, "results", "phase3")
os.makedirs(LIG_DIR, exist_ok=True); os.makedirs(OUT_DIR, exist_ok=True)
VINA = r"D:/vina"
OBABEL = r"C:/Users/Jesse/.workbuddy/binaries/python/envs/drugreflector/Scripts/obabel.exe"
PDB = os.path.join(REC_DIR, "2WOU.pdb")

# 化合物集合：known = 已知 ALK5 抑制剂(阳性对照)；hit = DrugReflector top 命中；control = 非激酶对照
LIGANDS = {
    # --- 已知 ALK5 抑制剂（基准 actives）---
    "LY-364947":   ("known", "c1[nH]nc(c1-c1ccnc2ccccc12)-c1ccccn1"),
    "SB-431542":   ("known", "NC(=O)c1ccc(cc1)-c1nc(c([nH]1)-c1ccccn1)-c1ccc2OCOc2c1"),
    "galunisertib":("known", "Cc1ccc(cc1)C2=NC3=C(C=NN3C2=O)N4CCN(CC4)C(=O)c5ccncc5"),
    "SB-505124":   ("known", "CC(C)(C)c1ccc(cc1)n2c3ccccc3nc2-c2ccc3OCOc3c2"),
    "RepSox":      ("known", "c1ccncc1-c1nnc(-c2ccccn2)c2ccccc12"),
    # --- DrugReflector 逆转签名 top 命中 ---
    "Y-27632":     ("hit", "C[C@@H](N)[C@H]1CC[C@@H](CC1)C(=O)Nc1ccncc1"),
    "Y-39983":     ("hit", "C[C@@H](N)c1ccc(cc1)C(=O)Nc1ccnc2[nH]ccc12"),
    "fasudil":     ("hit", "O=S(=O)(N1CCCNCC1)c1cccc2cnccc12"),
    "YC-1":        ("hit", "OCc1ccc(o1)c2nn(Cc3ccccc3)c4ccccc24"),
    "dipyridamole":("hit", "OCCN(CCO)c1nc(N2CCCCC2)c2nc(nc(N3CCCCC3)c2n1)N(CCO)CCO"),
    # --- 非激酶/非靶点对照（应结合较差）---
    "irbesartan":  ("control", "CCCCC1=NC2(CCCC2)C(=O)N1Cc1ccc(cc1)-c1ccccc1-c1nnn[nH]1"),
    "telmisartan": ("control", "CCCc1nc2c(C)cc(cc2n1Cc1ccc(cc1)-c1ccccc1C(O)=O)-c1nc2ccccc2n1"),
}

def split_and_prep_receptor():
    prot, lig = [], []
    for line in open(PDB):
        if line.startswith("ATOM") and line[21] == "A":
            prot.append(line)
        elif line.startswith("HETATM") and line[17:20].strip() == "ZZF":
            lig.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    prot_pdb = os.path.join(REC_DIR, "2WOU_A.pdb")
    open(prot_pdb, "w").writelines(prot)
    rec_pdbqt = os.path.join(REC_DIR, "2WOU_A.pdbqt")
    # 用 obabel -xr 生成刚性受体 PDBQT（加极性氢、pH 7.4），pybel 默认会写柔性格式
    subprocess.run([OBABEL, prot_pdb, "-O", rec_pdbqt, "-xr", "-p", "7.4"],
                   capture_output=True, timeout=300)
    L = np.array(lig); center = L.mean(0); ext = (L.max(0) - L.min(0)) + 10.0  # 留 10Å 边距
    center = np.round(center, 2); size = np.round(np.maximum(ext, 22.0), 2)   # 至少 22Å
    print(f"[receptor] {len(prot)} atoms; box center={center.tolist()} size={size.tolist()}")
    return center, size

def prep_ligand(name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    if AllChem.EmbedMolecule(m, randomSeed=42, useRandomCoords=True) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    sdf = os.path.join(LIG_DIR, name + ".sdf")
    pdbqt = os.path.join(LIG_DIR, name + ".pdbqt")
    Chem.MolToMolFile(m, sdf)
    mol = next(pybel.readfile("sdf", sdf))
    mol.write("pdbqt", pdbqt, overwrite=True)
    return pdbqt

def dock(name, pdbqt, center, size, exh=8):
    out = os.path.join(LIG_DIR, name + "_out.pdbqt")
    cmd = [VINA, "--receptor", os.path.join(REC_DIR, "2WOU_A.pdbqt"),
           "--ligand", pdbqt, "--out", out,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exh), "--num_modes", "9", "--seed", "42"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    aff = None
    for line in r.stdout.splitlines():
        if "REMARK VINA RESULT" in line:
            aff = float(line.split()[3]); break
    if aff is None:  # fallback: parse out file
        if os.path.exists(out):
            for line in open(out):
                if line.startswith("REMARK VINA RESULT"):
                    aff = float(line.split()[3]); break
    return aff, r.returncode

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="*", default=None, help="只跑指定化合物")
    ap.add_argument("--exh", type=int, default=8)
    args = ap.parse_args()
    center, size = split_and_prep_receptor()
    targets = {k: v for k, v in LIGANDS.items() if args.names is None or k in args.names}
    rows = []
    for name, (cat, smi) in targets.items():
        p = prep_ligand(name, smi)
        if p is None:
            print(f"[skip] {name}: SMILES 解析/3D 失败"); continue
        aff, rc = dock(name, p, center, size, exh=args.exh)
        print(f"[dock] {name:14s} aff={aff}  rc={rc}")
        rows.append({"name": name, "category": cat, "affinity": aff, "smiles": smi})
    import pandas as pd
    df = pd.DataFrame(rows)
    out = os.path.join(OUT_DIR, "docking_results.csv")
    if os.path.exists(out) and args.names:
        old = pd.read_csv(out); df = pd.concat([old, df]).drop_duplicates(subset="name", keep="last")
    df.to_csv(out, index=False)
    print(f"[done] {len(df)} compounds -> {out}")
    print(df.sort_values("affinity")[["name", "category", "affinity"]].to_string(index=False))

if __name__ == "__main__":
    main()
