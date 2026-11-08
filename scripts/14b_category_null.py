"""
14b_category_null.py — 类别级置换检验（基因集式富集）

单个化合物的经验 P 值受限于库内重复条目，统计效力弱。
这里把化合物按药类分组，检验「该类药物的中位排名」是否显著优于随机签名下的期望。

零假设：任取一个随机签名，这些药类的中位排名与真实签名无差别。
用法： python 14b_category_null.py --n 200
产出： results/phase3/null_category.csv / null_category_summary.txt
"""
import argparse, os, time
import numpy as np
import pandas as pd

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
RES = os.path.join(ROOT, "results")
CKPT = [os.path.join(ROOT, f"data/checkpoints/model_fold_{i}.pt") for i in range(3)]

GROUPS = {
    "ROCK inhibitors (hits)": ["y-27632", "y-39983", "fasudil"],
    "ALK5 inhibitors (known)": ["ly-364947", "sb-431542", "galunisertib", "sb-505124", "repsx", "repsox"],
    "RAS inhibitors (ARB)": ["losartan", "valsartan", "irbesartan", "telmisartan",
                             "candesartan", "olmesartan"],
    "RAS inhibitors (ACEi)": ["enalapril", "lisinopril", "ramipril", "perindopril", "captopril"],
    "Nephrotoxic agents": ["cisplatin", "gentamicin", "amphotericin", "vancomycin",
                           "cyclosporine", "tacrolimus", "ibuprofen", "naproxen",
                           "acyclovir", "methotrexate", "acetaminophen"],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--sig", default="vscore_reversal_injured2normal.csv")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()

    from drugreflector import DrugReflector

    sig_df = pd.read_csv(os.path.join(RES, args.sig))
    sig_df = sig_df.drop_duplicates(subset=sig_df.columns[0]).set_index(sig_df.columns[0])
    real = sig_df.iloc[:, 0].dropna()

    model = DrugReflector(checkpoint_paths=CKPT)
    names = list(model.compound_names)
    b2n = pd.read_csv(os.path.join(ROOT, "data/brd_to_name.tsv"), sep="\t",
                      index_col=0, dtype=str)["pert_iname"].to_dict()
    name_of = {b: b2n.get(b, b).lower() for b in names}

    # 预先把每个药类映射到库内 BRD id
    idx = {}
    for g, drugs in GROUPS.items():
        ids = [b for b, nm in name_of.items() if any(d in nm for d in drugs)]
        idx[g] = ids
        print(f"[group] {g:26s} {len(ids)} library entries")

    def group_medians(ranks):
        out = {}
        for g, ids in idx.items():
            out[g] = np.median([int(ranks[b]) for b in ids]) + 1 if ids else np.nan  # 1-based
        return out

    print("[real] ranking ...")
    pred = model.predict(real.rename("real"), n_top=None)
    real_med = group_medians(pred[("rank", "real")])
    for g, v in real_med.items():
        print(f"    {g:26s} median position {v:7.0f}")

    rng = np.random.default_rng(args.seed)
    vals = real.values.copy()
    rows, t0 = [], time.time()
    for i in range(args.n):
        perm = pd.Series(rng.permutation(vals), index=real.index, name=f"p{i}")
        pm = group_medians(model.predict(perm, n_top=None)[("rank", f"p{i}")])
        pm["iter"] = i
        rows.append(pm)
        if (i + 1) % 25 == 0:
            el = time.time() - t0
            print(f"    perm {i+1}/{args.n} ({el:.0f}s)", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "phase3/null_category.csv"), index=False)

    lines = [f"Category-level permutation test: n = {args.n} shuffled signatures",
             f"Compounds in library: {len(name_of)}", "",
             f"{'group':<28}{'n_lib':>7}{'real':>8}{'null median':>13}{'null 5%':>9}{'empirical P':>13}",
             "-" * 80]
    for g in GROUPS:
        rv = real_med[g]
        null = df[g].dropna().values
        if np.isnan(rv) or len(null) == 0:
            continue
        p_emp = (np.sum(null <= rv) + 1) / (len(null) + 1)
        lines.append(f"{g:<28}{len(idx[g]):>7}{rv:>8.0f}{np.median(null):>13.0f}"
                     f"{np.percentile(null, 5):>9.0f}{p_emp:>13.3f}")
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(RES, "phase3/null_category_summary.txt"), "w", encoding="utf-8").write(txt)


if __name__ == "__main__":
    main()
