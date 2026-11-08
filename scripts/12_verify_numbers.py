"""核对正文中所有排名类数字（1-based position，与稿件一致）"""
import os
import numpy as np
import pandas as pd

D = r"E:\workbuddy\2026-11-08-00-03-12"
b2n = pd.read_csv(os.path.join(D, "data/brd_to_name.tsv"), sep="\t",
                  index_col=0, dtype=str)["pert_iname"].to_dict()
df = pd.read_csv(os.path.join(D, "results/ranks_reversalinjured_normal.csv"),
                 header=[0, 1], index_col=0)
rcol = ("rank", df.columns.get_level_values(1)[0])
ranks0 = df[rcol].astype(int)            # 0-based
pos = ranks0 + 1                          # 1-based position
name_of = {brd: b2n.get(brd, brd) for brd in pos.index}
TOT = len(pos)
print(f"total compounds = {TOT}")

def lookup(drug):
    hits = [b for b, nm in name_of.items() if drug in nm.lower()]
    return int(pos[hits].min()) if hits else None

PROTECTIVE = ["empagliflozin", "dapagliflozin", "canagliflozin",
              "losartan", "valsartan", "telmisartan", "irbesartan",
              "ramipril", "lisinopril", "enalapril",
              "acetylcysteine", "bardoxolone", "pirfenidone", "nintedanib",
              "pioglitazone", "metformin", "finerenone", "spironolactone"]
NEPHROTOXIC = ["cisplatin", "gentamicin", "amphotericin", "vancomycin",
               "cyclosporine", "tacrolimus", "ibuprofen", "naproxen",
               "acyclovir", "methotrexate", "acetaminophen"]

print("\n--- 保护药 / 肾毒性药（1-based position）---")
for label, drugs in [("PROTECTIVE", PROTECTIVE), ("NEPHROTOXIC", NEPHROTOXIC)]:
    found = [(d, lookup(d)) for d in drugs if lookup(d) is not None]
    vals = [v for _, v in found]
    print(f"{label}: {len(found)}/{len(drugs)} found | median {np.median(vals):.0f} | "
          f"min {min(vals)} max {max(vals)}")
    for d, v in sorted(found, key=lambda x: x[1]):
        print(f"    {d:16s} {v:5d}  (top {100*v/TOT:.2f}%)")

print("\n--- 正文点名的三个 ---")
for d in ["irbesartan", "telmisartan", "enalapril", "valsartan", "losartan"]:
    print(f"    {d:14s} position {lookup(d)}")

print("\n--- CKD 19 药集合 ---")
CKD = {
    "RASi_ARB": ["losartan", "valsartan", "irbesartan", "telmisartan", "candesartan", "olmesartan"],
    "RASi_ACEi": ["enalapril", "lisinopril", "ramipril", "perindopril", "captopril"],
    "MRA": ["spironolactone", "eplerenone", "finerenone"],
    "other_antifibrotic": ["pirfenidone", "nintedanib", "bardoxolone", "pentoxifylline",
                           "allopurinol", "febuxostat"],
    "metabolic": ["metformin", "pioglitazone", "atorvastatin", "rosuvastatin", "ezetimibe"],
    "SGLT2i": ["empagliflozin", "dapagliflozin", "canagliflozin", "ertugliflozin"],
}
allv = []
for cat, drugs in CKD.items():
    rows = [(d, lookup(d)) for d in drugs if lookup(d) is not None]
    print(f"  {cat:20s} {len(rows)}/{len(drugs)}")
    for d, v in sorted(rows, key=lambda x: x[1]):
        print(f"      {d:16s} {v:5d}")
    allv += [v for _, v in rows]
allv = np.array(allv)
print(f"\n  CKD 命中总数 = {len(allv)} | median = {np.median(allv):.0f} | "
      f"top5% 数量 = {(allv <= 0.05*TOT).sum()} | top1.3% 数量 = {(allv <= 0.013*TOT).sum()}")
