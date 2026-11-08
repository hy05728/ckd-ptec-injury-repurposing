"""
17_build_benchmark_set.py — 构建标准规模 ALK5 对接基准集

actives : ChEMBL CHEMBL4439 (human TGFBR1) IC50 <= 1 uM，按 Bemis-Murcko 骨架去冗余
decoys  : 按 DUD-E 匹配原则（MW / cLogP / rotB / HBD / HBA / net charge）从 L1000 药库抽取，
          并要求与 actives 拓扑不相似（Tanimoto < 0.35）

产出： data/phase3/benchmark_v2/actives.smi / decoys.smi / descriptor_report.txt
"""
import os, json, time, urllib.request
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
OUT = os.path.join(ROOT, "data/phase3/benchmark_v2")
os.makedirs(OUT, exist_ok=True)
TARGET = "CHEMBL4439"          # human TGF-beta receptor type-1
BASE = "https://www.ebi.ac.uk/chembl/api/data"


def fetch_json(url, retries=4):
    for i in range(retries):
        try:
            return json.load(urllib.request.urlopen(url, timeout=60))
        except Exception as e:
            print(f"    retry {i+1}: {type(e).__name__}")
            time.sleep(3)
    return None


def get_actives(max_n=120):
    print("[1/4] ChEMBL activities (IC50 <= 1 uM, human TGFBR1)")
    ids = {}   # molecule_chembl_id -> best IC50
    for st in ("IC50", "Ki"):
        offset = 0
        while True:
            u = (f"{BASE}/activity.json?target_chembl_id={TARGET}&standard_type={st}"
                 f"&standard_units=nM&standard_value__lte=1000&limit=1000&offset={offset}")
            d = fetch_json(u)
            if not d or not d.get("activities"):
                break
            for a in d["activities"]:
                mid = a.get("molecule_chembl_id")
                val = a.get("standard_value")
                if mid and val:
                    try:
                        val = float(val)
                    except Exception:
                        continue
                    if mid not in ids or val < ids[mid]:
                        ids[mid] = val
            total = d["page_meta"]["total_count"]
            offset += 1000
            print(f"      {st}: {min(offset, total)}/{total}")
            if offset >= total:
                break
    print(f"    unique molecules: {len(ids)}")

    print("[2/4] fetching SMILES in batches")
    mids = list(ids)
    smi = {}
    for i in range(0, len(mids), 50):
        batch = ",".join(mids[i:i + 50])
        d = fetch_json(f"{BASE}/molecule.json?molecule_chembl_id__in={batch}&limit=50&only=molecule_chembl_id,molecule_structures")
        if not d:
            continue
        for m in d.get("molecules", []):
            st = m.get("molecule_structures") or {}
            s = st.get("canonical_smiles")
            if s:
                smi[m["molecule_chembl_id"]] = s
        if i % 500 == 0:
            print(f"      {i}/{len(mids)}")
    print(f"    with SMILES: {len(smi)}")

    # 去重 + 骨架聚类
    print("[3/4] dedup + Murcko scaffold clustering")
    seen, rows = set(), []
    for mid, s in smi.items():
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            continue
        can = Chem.MolToSmiles(mol)
        if can in seen:
            continue
        seen.add(can)
        mw = Descriptors.MolWt(mol)
        if not (200 <= mw <= 600):
            continue
        if rdMolDescriptors.CalcNumRotatableBonds(mol) > 15:
            continue
        try:
            sc = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:
            sc = ""
        rows.append(dict(chembl_id=mid, smiles=can, ic50_nM=ids[mid], scaffold=sc, mw=mw))
    df = pd.DataFrame(rows).sort_values("ic50_nM")
    # 每个骨架最多 2 个，保证化学多样性
    df["n_scaf"] = df.groupby("scaffold").cumcount()
    df = df[df.n_scaf < 2].drop(columns="n_scaf")
    print(f"    deduplicated: {len(df)}, after scaffold cap(2): {len(df)}")
    if len(df) > max_n:
        df = df.sort_values("ic50_nM").head(max_n)
    return df.reset_index(drop=True)


def props(mol):
    return dict(
        mw=Descriptors.MolWt(mol),
        logp=Crippen.MolLogP(mol),
        rotb=rdMolDescriptors.CalcNumRotatableBonds(mol),
        hbd=rdMolDescriptors.CalcNumHBD(mol),
        hba=rdMolDescriptors.CalcNumHBA(mol),
        chg=Chem.GetFormalCharge(mol),
    )


def build_decoys(actives_df, n_target=400):
    print("[3/4] building property-matched decoys from L1000 library")
    pi = os.path.join(ROOT, r"data/GSE92742_pert_info.txt.gz")
    lib = pd.read_csv(pi, sep="\t", compression="gzip", dtype=str,
                      usecols=["pert_id", "pert_iname", "canonical_smiles"])
    lib = lib.dropna(subset=["canonical_smiles"]).drop_duplicates("canonical_smiles")
    print(f"    library with SMILES: {len(lib)}")

    names = {str(x).lower() for x in lib.pert_iname.dropna()}
    known = ["ly-364947", "sb-431542", "galunisertib", "sb-505124", "repsox", "ly2157299"]
    lib = lib[~lib.pert_iname.fillna("").str.lower().isin(known)]
    print(f"    after removing known ALK5 inhibitors: {len(lib)}")

    aprops = [props(Chem.MolFromSmiles(s)) for s in actives_df.smiles]
    askel = [Chem.RDKFingerprint(Chem.MolFromSmiles(s)) for s in actives_df.smiles]

    cand = []
    for _, r in lib.iterrows():
        mol = Chem.MolFromSmiles(r.canonical_smiles)
        if mol is None:
            continue
        p = props(mol)
        if not (150 <= p["mw"] <= 650):
            continue
        fp = Chem.RDKFingerprint(mol)
        # 拓扑不相似（DUD-E 要求 decoy 与 active 骨架不同）
        maxsim = max(DataStructs.TanimotoSimilarity(fp, a) for a in askel)
        if maxsim >= 0.35:
            continue
        # 与任一 active 的理化性质匹配
        matched = False
        for ap in aprops:
            if (abs(p["mw"] - ap["mw"]) / ap["mw"] < 0.28
                    and abs(p["logp"] - ap["logp"]) < 1.6
                    and abs(p["rotb"] - ap["rotb"]) <= 2
                    and abs(p["hbd"] - ap["hbd"]) <= 1
                    and abs(p["hba"] - ap["hba"]) <= 2
                    and p["chg"] == ap["chg"]):
                matched = True
                break
        if matched:
            cand.append(dict(pert_id=r.pert_id, name=r.pert_iname, smiles=r.canonical_smiles,
                             **p))
    cand = pd.DataFrame(cand).drop_duplicates("smiles")
    print(f"    property-matched candidates: {len(cand)}")
    if len(cand) > n_target:
        cand = cand.sample(n_target, random_state=42)
    return cand.reset_index(drop=True)


def main():
    act = get_actives(max_n=120)
    act.to_csv(os.path.join(OUT, "actives.csv"), index=False)
    with open(os.path.join(OUT, "actives.smi"), "w") as f:
        for _, r in act.iterrows():
            f.write(f"{r.smiles} {r.chembl_id}\n")

    dec = build_decoys(act, n_target=400)
    dec.to_csv(os.path.join(OUT, "decoys.csv"), index=False)
    with open(os.path.join(OUT, "decoys.smi"), "w") as f:
        for _, r in dec.iterrows():
            f.write(f"{r.smiles} {r.pert_id}\n")

    lines = [
        "ALK5 (TGFBR1) docking benchmark set",
        f"target: ChEMBL {TARGET} (human TGF-beta receptor type-1)",
        f"actives: {len(act)}  (IC50 <= 1 uM, deduplicated, max 2 per Murcko scaffold)",
        f"decoys : {len(dec)}  (L1000 library, DUD-E matching criteria, Tanimoto < 0.35 to actives)",
        f"ratio  : 1 : {len(dec)/max(len(act),1):.1f}",
        "",
        "actives value ranges:",
        f"  IC50 (nM): median {act.ic50_nM.median():.0f}, range {act.ic50_nM.min():.0f}-{act.ic50_nM.max():.0f}",
        f"  MW       : median {act.mw.median():.0f}",
        "",
        "decoys property ranges:",
        f"  MW   : median {dec.mw.median():.0f}, range {dec.mw.min():.0f}-{dec.mw.max():.0f}",
        f"  cLogP: median {dec.logp.median():.1f}",
        f"  rotB : median {dec.rotb.median():.0f}",
    ]
    txt = "\n".join(lines)
    print("\n" + txt)
    open(os.path.join(OUT, "descriptor_report.txt"), "w").write(txt)


if __name__ == "__main__":
    main()
