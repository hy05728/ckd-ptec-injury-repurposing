#!/usr/bin/env python3
"""
15_dock_rock1.py — ROCK1 交叉对接，用于证明「机制分层」：
  ROCK 抑制剂应对 ROCK1 结合强、对 ALK5 弱；ALK5 抑制剂则相反。

受体: PDB 2ETR chain A（ROCK1 激酶域，与 Y-27632 共晶），盒子取自共晶配体 Y27。
化合物: 与 04_dock_alk5.py 相同的 11 个 + telmisartan。
产出: results/phase3/rock1_docking_results.csv
"""
import os, subprocess, argparse
import numpy as np
import pandas as pd
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
P3 = os.path.join(ROOT, "data", "phase3")
REC_DIR = os.path.join(P3, "receptor")
LIG_DIR = os.path.join(P3, "ligands_rock1")
OUT_DIR = os.path.join(ROOT, "results", "phase3")
os.makedirs(LIG_DIR, exist_ok=True)
VINA = r"D:/vina"
OBABEL = r"C:/Users/Jesse/.workbuddy/binaries/python/envs/drugreflector/Scripts/obabel.exe"
PDB_ID, CHAIN, COFACTOR = "2ETR", "A", "Y27"

# 与 ALK5 对接保持完全相同的 SMILES 集合
LIGANDS = {
    "LY-364947":    ("ALK5_known", "c1[nH]nc(c1-c1ccnc2ccccc12)-c1ccccn1"),
    "SB-431542":    ("ALK5_known", "NC(=O)c1ccc(cc1)-c1nc(c([nH]1)-c1ccccn1)-c1ccc2OCOc2c1"),
    "galunisertib": ("ALK5_known", "Cc1ccc(cc1)C2=NC3=C(C=NN3C2=O)N4CCN(CC4)C(=O)c5ccncc5"),
    "SB-505124":    ("ALK5_known", "CC(C)(C)c1ccc(cc1)n2c3ccccc3nc2-c2ccc3OCOc3c2"),
    "RepSox":       ("ALK5_known", "c1ccncc1-c1nnc(-c2ccccn2)c2ccccc12"),
    "Y-27632":      ("ROCK_hit", "C[C@@H](N)[C@H]1CC[C@@H](CC1)C(=O)Nc1ccncc1"),
    "Y-39983":      ("ROCK_hit", "C[C@@H](N)c1ccc(cc1)C(=O)Nc1ccnc2[nH]ccc12"),
    "fasudil":      ("ROCK_hit", "O=S(=O)(N1CCCNCC1)c1cccc2cnccc12"),
    "YC-1":         ("other_hit", "OCc1ccc(o1)c2nn(Cc3ccccc3)c4ccccc24"),
    "dipyridamole": ("other_hit", "OCCN(CCO)c1nc(N2CCCCC2)c2nc(nc(N3CCCCC3)c2n1)N(CCO)CCO"),
    "irbesartan":   ("control", "CCCCC1=NC2(CCCC2)C(=O)N1Cc1ccc(cc1)-c1ccccc1-c1nnn[nH]1"),
    "telmisartan":  ("control", "CCCc1nc2c(C)cc(cc2n1Cc1ccc(cc1)-c1ccccc1C(O)=O)-c1nc2ccccc2n1"),
}


def prep_receptor():
    pdb = os.path.join(REC_DIR, f"{PDB_ID}.pdb")
    prot, lig = [], []
    for line in open(pdb):
        if line.startswith("ATOM") and line[21] == CHAIN:
            prot.append(line)
        elif (line.startswith("HETATM") and line[17:20].strip() == COFACTOR
              and line[21] == CHAIN):          # 必须限链：2ETR 的 A/B 链各有一个 Y27
            lig.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    out_pdb = os.path.join(REC_DIR, f"{PDB_ID}_{CHAIN}.pdb")
    open(out_pdb, "w").writelines(prot)
    out_pdbqt = os.path.join(REC_DIR, f"{PDB_ID}_{CHAIN}.pdbqt")
    subprocess.run([OBABEL, out_pdb, "-O", out_pdbqt, "-xr", "-p", "7.4"],
                   capture_output=True, timeout=300)
    L = np.array(lig)
    center = np.round(L.mean(0), 2)
    size = np.round(np.maximum((L.max(0) - L.min(0)) + 10.0, 22.0), 2)
    print(f"[receptor] {PDB_ID} chain {CHAIN}: {len(prot)} atoms, {len(lig)} cofactor atoms")
    print(f"[box] center={center.tolist()} size={size.tolist()}")
    return out_pdbqt, center, size


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
    next(pybel.readfile("sdf", sdf)).write("pdbqt", pdbqt, overwrite=True)
    return pdbqt


def dock(rec_pdbqt, name, pdbqt, center, size, exh=8):
    out = os.path.join(LIG_DIR, name + "_out.pdbqt")
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", pdbqt, "--out", out,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exh), "--num_modes", "9", "--seed", "42"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    for line in r.stdout.splitlines():
        if "REMARK VINA RESULT" in line:
            return float(line.split()[3]), r.returncode
    if os.path.exists(out):
        for line in open(out):
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3]), r.returncode
    return None, r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exh", type=int, default=8)
    ap.add_argument("--names", nargs="*", default=None)
    args = ap.parse_args()

    rec, center, size = prep_receptor()
    targets = {k: v for k, v in LIGANDS.items()
               if args.names is None or k in args.names}
    rows = []
    for name, (cat, smi) in targets.items():
        p = prep_ligand(name, smi)
        if p is None:
            print(f"[skip] {name}"); continue
        aff, rc = dock(rec, name, p, center, size, args.exh)
        print(f"[dock] {name:14s} {cat:12s} aff={aff}  rc={rc}", flush=True)
        rows.append({"name": name, "category": cat, "affinity": aff, "smiles": smi})

    df = pd.DataFrame(rows)
    out = os.path.join(OUT_DIR, "rock1_docking_results.csv")
    df.to_csv(out, index=False)
    print(f"\n[done] {len(df)} compounds -> {out}")
    print(df.sort_values("affinity")[["name", "category", "affinity"]].to_string(index=False))


if __name__ == "__main__":
    main()
