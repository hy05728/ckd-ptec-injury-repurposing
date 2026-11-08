#!/usr/bin/env python3
"""
03_rank.py — 用真实 DrugReflector API 跑排序 + 已知肾脏药物回溯验证。

关键：模型学的是 P(化合物 | 签名)，即"哪个药会诱导该签名"。
  - 输入损伤签名(vscore_injury_normal2injured) → 排在前的 = 会诱导损伤的药
  - 输入逆转签名(vscore_reversal_injured2normal) → 排在前的 = 逆转损伤(保护)的药
回溯验证：保护药(如 SGLT2i/RASi)应在逆转签名里排得靠前；肾毒性药(顺铂/庆大霉素)应靠后。
"""
import os, sys, warnings
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")
ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")
CKPT = [os.path.join(DATA, "checkpoints", f"model_fold_{i}.pt") for i in range(3)]

# 已知肾脏保护 / 肾毒性药物（回溯名单，CMap 里多数存在）
PROTECTIVE = [
    "empagliflozin", "dapagliflozin", "canagliflozin",
    "losartan", "valsartan", "telmisartan", "irbesartan",
    "ramipril", "lisinopril", "enalapril",
    "acetylcysteine", "bardoxolone", "pirfenidone", "nintedanib",
    "pioglitazone", "metformin", "finerenone", "spironolactone",
]
NEPHROTOXIC = [
    "cisplatin", "gentamicin", "amphotericin", "vancomycin",
    "cyclosporine", "tacrolimus", "ibuprofen", "naproxen",
    "acyclovir", "methotrexate", "acetaminophen",
]

def load_signature(path, name):
    s = pd.read_csv(path, index_col=0)["vscore"]
    s.name = name
    return s

def match_compounds(name_of, query):
    """对 BRD->药名 映射后的药名做大小写不敏感子串匹配，返回 BRD ID 列表。"""
    q = query.lower()
    hits = [brd for brd, nm in name_of.items() if q in nm.lower()]
    return hits

def main():
    from drugreflector import DrugReflector
    for p in CKPT:
        assert os.path.exists(p), f"missing checkpoint {p}"
    model = DrugReflector(checkpoint_paths=CKPT)
    names = list(model.compound_names)
    print(f"[model] {len(names)} compounds; sample: {names[:5]}")

    # BRD -> 药名映射
    brd2name = pd.read_csv(os.path.join(DATA, "brd_to_name.tsv"), sep="\t",
                           index_col=0, dtype=str)["pert_iname"].to_dict()
    name_of = {brd: brd2name.get(brd, brd) for brd in names}

    out_rows = []
    for tag, path, sname in [
        ("injury(normal->injured)", "vscore_injury_normal2injured.csv", "injury"),
        ("reversal(injured->normal)", "vscore_reversal_injured2normal.csv", "reversal"),
    ]:
        fp = os.path.join(RESULTS, path)
        if not os.path.exists(fp):
            print(f"!! missing {fp}, skip")
            continue
        sig = load_signature(fp, sname)
        pred = model.predict(sig, n_top=None)
        rank_col = ("rank", sname)
        ranks = pred[rank_col]  # index = compound name, value = rank (1=best)
        print(f"\n=== {tag} : top 20 ===")
        print(ranks.nsmallest(20).to_string())

        # 回溯
        def grp_stats(drugs, label):
            found = {}
            for d in drugs:
                hits = match_compounds(name_of, d)
                if hits:
                    best = min(ranks[h] for h in hits)
                    found[d] = (best, [name_of[h] for h in hits][:2])
            if found:
                med = np.median([v[0] for v in found.values()])
                print(f"  [{label}] {len(found)}/{len(drugs)} 命中, median rank = {med:.0f} (共 {len(ranks)} 化合物)")
                for d, (best, hitnames) in sorted(found.items(), key=lambda kv: kv[1][0]):
                    print(f"      {d}: best rank {best}  (matched: {hitnames})")
            else:
                print(f"  [{label}] 0/{len(drugs)} 命中")
            return found

        grp_stats(PROTECTIVE, "protective")
        grp_stats(NEPHROTOXIC, "nephrotoxic")

        # 存全量排名
        pred.to_csv(os.path.join(RESULTS, f"ranks_{tag.replace('(', '').replace(')', '').replace('->', '_')}.csv"))

    print("\n[done]")

if __name__ == "__main__":
    main()
