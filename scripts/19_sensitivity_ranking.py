"""
19_sensitivity_ranking.py — 用「无监督定义」的签名重跑 DrugReflector，检验排位是否稳健

若 ROCK / ALK5 抑制剂在两种定义下都进前列，则"marker 定义"这一选择不影响结论。
另外计算二值化后的 ARI（把无监督簇合并成 injured-like / healthy-like 两类），
避免"2 类 vs 16 簇"因粒度不同导致 ARI 天然偏低。
"""
import os
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, cohen_kappa_score

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
RES = os.path.join(ROOT, "results")
P3 = os.path.join(RES, "phase3")
CKPT = [os.path.join(ROOT, f"data/checkpoints/model_fold_{i}.pt") for i in range(3)]

TARGETS = ["Y-27632", "Y-39983", "fasudil", "LY-364947", "YC-1",
           "SB-431542", "galunisertib", "SB-505124", "RepSox", "irbesartan"]


def main():
    print("=" * 70)
    print("sensitivity of the compound ranking to the injury-state definition")
    print("=" * 70)

    # ---- 1) 二值化一致性（更公平的 ARI）----
    # 从 18 脚本的输出重建：confusion 由 crosstab 得到，这里用保存的 vscore 文件做下游检验
    unsup = pd.read_csv(os.path.join(P3, "vscore_injury_unsup.csv"), index_col=0)
    orig = pd.read_csv(os.path.join(RES, "vscore_injury_normal2injured.csv"))
    orig.columns = ["gene", "vscore_marker"]
    mm = orig.set_index("gene").join(unsup, how="inner").dropna()
    print(f"\n[1] gene-level concordance: n = {len(mm)}")
    print(f"    Pearson  r = {mm.vscore_marker.corr(mm.vscore_unsup):.3f}")
    print(f"    Spearman r = {mm.vscore_marker.corr(mm.vscore_unsup, method='spearman'):.3f}")

    # ---- 2) 用无监督签名重跑 DrugReflector ----
    from drugreflector import DrugReflector
    model = DrugReflector(checkpoint_paths=CKPT)
    names = list(model.compound_names)
    b2n = pd.read_csv(os.path.join(ROOT, "data/brd_to_name.tsv"), sep="\t",
                      index_col=0, dtype=str)["pert_iname"].to_dict()
    name_of = {b: b2n.get(b, b).lower() for b in names}

    def best(ranks, drug):
        h = [b for b, nm in name_of.items() if drug.lower() in nm]
        return int(ranks[h].min()) + 1 if h else None   # 1-based position

    sig_unsup = (-unsup["vscore_unsup"]).rename("unsup")     # 反转方向
    pred = model.predict(sig_unsup, n_top=None)
    r_u = pred[("rank", "unsup")]

    df_orig = pd.read_csv(os.path.join(RES, "ranks_reversalinjured_normal.csv"), header=[0, 1], index_col=0)
    r_m = df_orig[("rank", df_orig.columns.get_level_values(1)[0])]

    print(f"\n[2] DrugReflector ranking under the two definitions (position among {len(names)})")
    print(f"    {'compound':<14}{'marker-based':>14}{'unsupervised':>14}{'delta':>8}")
    print("    " + "-" * 50)
    rows = []
    for d in TARGETS:
        p1, p2 = best(r_m, d), best(r_u, d)
        if p1 is None or p2 is None:
            continue
        print(f"    {d:<14}{p1:>14}{p2:>14}{p2-p1:>+8}")
        rows.append(dict(compound=d, marker_based=p1, unsupervised=p2, delta=p2 - p1))

    # 类别中位
    def cat_med(ranks, keys):
        ids = [b for b, nm in name_of.items() if any(k in nm for k in keys)]
        return int(np.median([ranks[b] for b in ids])) + 1 if ids else None

    print("\n[3] class-level median position")
    for label, keys in [("ROCK inhibitors", ["y-27632", "y-39983", "fasudil"]),
                        ("ALK5 inhibitors", ["ly-364947", "sb-431542"])]:
        print(f"    {label:<18} marker-based {cat_med(r_m, keys):>6}   "
              f"unsupervised {cat_med(r_u, keys):>6}")

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(P3, "sensitivity_ranking.csv"), index=False)
    rank_corr = out[["marker_based", "unsupervised"]].corr(method="spearman").iloc[0, 1]
    print(f"\n[4] Spearman correlation of the ten compound positions: {rank_corr:.3f}")
    print(f"    written -> {os.path.join(P3, 'sensitivity_ranking.csv')}")


if __name__ == "__main__":
    main()
