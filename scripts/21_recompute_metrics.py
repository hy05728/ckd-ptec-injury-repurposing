"""
21_recompute_metrics.py — 用正确方向重算基准指标（Vina affinity 越负越好）

注意：20_dock_benchmark_v2.py 里的 metrics() 默认「分数越大越好」，
而 Vina affinity 是越负越好，导致 ROC-AUC 被算成 1-AUC、EF 全为 0。
本脚本读取已完成的对接分数（无需重跑对接），按正确方向计算并给出 bootstrap CI。
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
OUT = os.path.join(ROOT, "results", "phase3")


def bedroc(y, score, alpha=20.0):
    """score 越大越好。"""
    y = np.asarray(y); s = np.asarray(score)
    n = len(y); n_act = int(y.sum())
    if n_act == 0:
        return np.nan
    order = np.argsort(-s, kind="stable")
    ranks = (np.where(y[order] == 1)[0] + 1) / n
    rie = np.sum(np.exp(-alpha * ranks)) / n_act
    return rie * alpha / (1 - np.exp(-alpha))


def metrics(y, score):
    """score 越大越好。"""
    y = np.asarray(y).astype(int); s = np.asarray(score, dtype=float)
    n = len(y); n_act = int(y.sum()); n_dec = n - n_act
    order = np.argsort(-s, kind="stable")
    ys = y[order]
    rk = rankdata(s)                      # 分数越大 rank 越大
    auc = (rk[y == 1].sum() - n_act * (n_act + 1) / 2) / (n_act * n_dec)
    out = {"roc_auc": auc}
    for pct in (1, 5, 10, 20):
        k = max(1, int(round(n * pct / 100)))
        out[f"ef{pct}"] = (ys[:k].sum() / n_act) / (k / n)
    out["bedroc20"] = bedroc(y, s)
    return out


def bootstrap(y, score, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y); s = np.asarray(score)
    idx = np.arange(len(y))
    vals = []
    for _ in range(n_boot):
        b = rng.choice(idx, size=len(idx), replace=True)
        if y[b].sum() == 0 or y[b].sum() == len(b):
            continue
        vals.append(metrics(y[b], s[b])["roc_auc"])
    return np.percentile(vals, [2.5, 97.5]), np.std(vals)


def main():
    df = pd.read_csv(os.path.join(OUT, "benchmark_v2_scores.csv")).dropna(subset=["affinity"])
    y = df.label.values.astype(int)
    score = -df.affinity.values          # 关键：affinity 越负越好 -> 取负后越大越好
    m = metrics(y, score)
    ci, sd = bootstrap(y, score)

    # 旧基准（5 actives + 60 decoys）作对比
    old_txt = ""
    try:
        dk = pd.read_csv(os.path.join(OUT, "docking_results.csv"))
        dec = pd.read_csv(os.path.join(OUT, "benchmark_decoys.csv"))
        y2 = np.r_[np.ones(len(dk[dk.category == "known"])), np.zeros(len(dec))]
        s2 = -np.r_[dk[dk.category == "known"].affinity.values, dec.affinity.values]
        m2 = metrics(y2, s2)
        old_txt = (f"\ncomparison with the original benchmark (5 actives + {len(dec)} decoys):\n"
                   f"  ROC-AUC {m2['roc_auc']:.3f}   EF5% {m2['ef5']:.1f}   BEDROC20 {m2['bedroc20']:.3f}")
    except Exception as e:
        old_txt = f"\n(original benchmark not comparable: {e})"

    lines = [
        "ALK5 docking benchmark — ChEMBL-scale (corrected orientation)",
        f"actives {int(y.sum())}   decoys {int((1-y).sum())}   ratio 1:{int((1-y).sum())/max(int(y.sum()),1):.1f}",
        "score = -Vina affinity (more negative = stronger binding = better)",
        "",
        f"ROC-AUC  {m['roc_auc']:.3f}   (95% CI {ci[0]:.3f}-{ci[1]:.3f}, bootstrap s.d. {sd:.3f})",
        f"EF1%     {m['ef1']:.2f}",
        f"EF5%     {m['ef5']:.2f}",
        f"EF10%    {m['ef10']:.2f}",
        f"EF20%    {m['ef20']:.2f}",
        f"BEDROC20 {m['bedroc20']:.3f}",
        old_txt,
    ]
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(OUT, "benchmark_v2_summary_corrected.txt"), "w").write(txt)


if __name__ == "__main__":
    main()
