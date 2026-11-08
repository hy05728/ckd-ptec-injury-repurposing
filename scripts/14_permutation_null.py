"""
14_permutation_null.py — 随机签名阴性对照（置换检验）

零假设：top 命中只是「任何签名都会给出这批药」。
做法：保持 v-score 的数值分布不变，随机打乱「基因 ↔ v-score」对应关系 N 次，
      重跑 DrugReflector，统计目标化合物的排名分布与经验 P 值。

用法： python 14_permutation_null.py --n 100
产出： results/phase3/null_permutation.csv / null_summary.txt
"""
import argparse, os, time
import numpy as np
import pandas as pd

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
RES = os.path.join(ROOT, "results")
CKPT_DIR = os.path.join(ROOT, "data/checkpoints")
CKPT = [os.path.join(CKPT_DIR, f"model_fold_{i}.pt") for i in range(3)]

# 关注的化合物（正文报告过的 top hits + ALK5 抑制剂 + 非靶点对照）
TARGETS = {
    "Y-27632": "ROCK inhibitor (hit)",
    "Y-39983": "ROCK inhibitor (hit)",
    "fasudil": "ROCK inhibitor (hit)",
    "LY-364947": "ALK5 inhibitor (hit)",
    "YC-1": "signature hit",
    "SB-431542": "ALK5 inhibitor (known)",
    "galunisertib": "ALK5 inhibitor (known)",
    "SB-505124": "ALK5 inhibitor (known)",
    "RepSox": "ALK5 inhibitor (known)",
    "irbesartan": "non-target control",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100, help="置换次数")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    from drugreflector import DrugReflector

    # 真实签名
    sig_df = pd.read_csv(os.path.join(RES, "vscore_reversal_injured2normal.csv"))
    gene_col = sig_df.columns[0]
    sig_df = sig_df.drop_duplicates(subset=gene_col).set_index(gene_col)
    real = sig_df.iloc[:, 0].dropna()
    print(f"[sig] {len(real)} genes, range [{real.min():.2f}, {real.max():.2f}]")

    model = DrugReflector(checkpoint_paths=CKPT)
    names = list(model.compound_names)
    b2n = pd.read_csv(os.path.join(ROOT, "data/brd_to_name.tsv"), sep="\t",
                      index_col=0, dtype=str)["pert_iname"].to_dict()
    name_of = {b: b2n.get(b, b) for b in names}

    def best_rank(ranks, drug):
        hits = [b for b, nm in name_of.items() if drug.lower() in nm.lower()]
        return int(min(ranks[h] for h in hits)) if hits else None

    # ---- 真实签名 ----
    print("[real] ranking ...")
    pred = model.predict(real.rename("real"), n_top=None)
    rcol = ("rank", "real")
    real_ranks = pred[rcol]
    real_pos = {d: best_rank(real_ranks, d) for d in TARGETS}
    for d, p in sorted(real_pos.items(), key=lambda kv: (kv[1] is None, kv[1])):
        print(f"    {d:14s} {p}")

    # ---- 置换 ----
    rng = np.random.default_rng(args.seed)
    vals = real.values.copy()
    rows = []
    t0 = time.time()
    for i in range(args.n):
        perm = pd.Series(rng.permutation(vals), index=real.index, name=f"perm{i}")
        pr = model.predict(perm, n_top=None)
        pranks = pr[("rank", f"perm{i}")]
        rec = {"iter": i}
        for d in TARGETS:
            rec[d] = best_rank(pranks, d)
        rows.append(rec)
        if (i + 1) % 10 == 0:
            el = time.time() - t0
            print(f"    perm {i+1}/{args.n}  ({el:.0f}s, ETA {el/(i+1)*(args.n-i-1):.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "phase3/null_permutation.csv"), index=False)

    # ---- 统计 ----
    lines = []
    lines.append(f"Null permutation test: n = {args.n} shuffled signatures")
    lines.append(f"Compounds in library: {len(name_of)}")
    lines.append("")
    lines.append(f"{'compound':<14}{'role':<28}{'real':>7}{'null median':>13}"
                 f"{'null 5%':>9}{'empirical P':>13}")
    lines.append("-" * 86)
    for d, role in TARGETS.items():
        rv = real_pos[d]
        null = df[d].dropna().astype(int).values
        if rv is None or len(null) == 0:
            continue
        p_emp = (np.sum(null <= rv) + 1) / (len(null) + 1)
        lines.append(f"{d:<14}{role:<28}{rv:>7}{np.median(null):>13.0f}"
                     f"{np.percentile(null, 5):>9.0f}{p_emp:>13.3f}")

    # top-N 捕获数
    lines.append("")
    lines.append("Capture of the ten target compounds within top-N:")
    for N in (10, 11, 50, 100):
        real_cap = sum(1 for d, v in real_pos.items() if v is not None and v <= N)
        null_cap = [sum(1 for d in TARGETS if pd.notna(df.loc[i, d]) and df.loc[i, d] <= N)
                    for i in range(len(df))]
        lines.append(f"  top-{N:<4} real = {real_cap:>2} / {len(TARGETS)}"
                     f"   null = {np.mean(null_cap):.2f} +/- {np.std(null_cap):.2f}"
                     f"   (max {max(null_cap)})")

    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(RES, "phase3/null_summary.txt"), "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    main()
