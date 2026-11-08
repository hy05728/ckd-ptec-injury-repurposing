#!/usr/bin/env python3
"""10_make_figures.py — 生成正文 Figure 1-5（300 dpi）。"""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = r"E:/workbuddy/2026-11-08-00-03-12"
RES = os.path.join(ROOT, "results"); P3 = os.path.join(RES, "phase3")
FIG = os.path.join(ROOT, "figs"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.linewidth":0.8,
                     "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300,
                     "ps.fonttype":42,"pdf.fonttype":42})   # TrueType，Elsevier 要求

def save(fig, name):
    """输出：PNG 预览(300dpi) / TIFF 投稿(600dpi, LZW) / PDF+EPS 矢量"""
    base = os.path.splitext(name)[0]
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, base + ".png"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(FIG, base + ".tif"), bbox_inches="tight", dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    # Elsevier 要求 RGB（不带 alpha）：转白底 RGB 并写 600 dpi
    from PIL import Image
    tp = os.path.join(FIG, base + ".tif")
    Image.open(tp).convert("RGB").save(tp, compression="tiff_lzw", dpi=(600, 600))
    fig.savefig(os.path.join(FIG, base + ".pdf"), bbox_inches="tight")
    try:
        fig.savefig(os.path.join(FIG, base + ".eps"), bbox_inches="tight")
    except Exception as e:
        print("  !! eps failed:", e)
    plt.close(fig)
    print("saved", base, "(png/tif/pdf/eps)")

# ---------- Figure 1: signature construction & validation ----------
sig = pd.read_csv(os.path.join(RES,"vscore_injury_normal2injured.csv"))
sig.columns=["gene","vscore"]; sig=sig.dropna()
fig,ax=plt.subplots(1,4,figsize=(14,3.2))
ax[0].hist(sig.vscore.clip(-3,3),bins=80,color="#4C72B0",alpha=.85)
ax[0].axvline(0,color="k",lw=.8,ls="--"); ax[0].set_xlabel("v-score (healthy → injured PT)"); ax[0].set_ylabel("genes"); ax[0].set_title("A  Signature distribution")
mk=["VCAM1","KRT8","HAVCR1","LCN2"]
vals=[float(sig.loc[sig.gene==g,"vscore"].iloc[0]) if (sig.gene==g).any() else np.nan for g in mk]
ax[1].bar(mk,vals,color=["#DD8452" if v>0 else "#55A868" for v in vals])
ax[1].axhline(0,color="k",lw=.8); ax[1].set_ylabel("v-score"); ax[1].set_title("B  Injury markers (expected ↑)")
up=sig.sort_values("vscore",ascending=False).head(8); dn=sig.sort_values("vscore").head(8)
lab=list(up.gene)+list(dn.gene); v=list(up.vscore)+list(dn.vscore)
ax[2].barh(range(len(lab))[::-1],v,color=["#C44E52" if x>0 else "#4C72B0" for x in v])
ax[2].set_yticks(range(len(lab))[::-1]); ax[2].set_yticklabels(lab,fontsize=6); ax[2].axvline(0,color="k",lw=.8)
ax[2].set_title("C  Top up / down genes"); ax[2].set_xlabel("v-score")
# D: 与 Lake 2023 独立队列的一致性
import os as _os
lk=_os.path.join(P3,"vscore_injury_lake2023_markerdef.csv")
if _os.path.exists(lk):
    L=pd.read_csv(lk); L.columns=["gene","vscore_lake"]; L=L.drop_duplicates("gene")
    mm=sig[["gene","vscore"]].merge(L,on="gene"); mm=mm.dropna()
    ax[3].scatter(mm.vscore,mm.vscore_lake,s=3,alpha=.25,color="#4C72B0",rasterized=True)
    ax[3].axhline(0,color="k",lw=.5); ax[3].axvline(0,color="k",lw=.5)
    from scipy.stats import pearsonr
    r=pearsonr(mm.vscore,mm.vscore_lake)[0]
    ax[3].set_xlabel("v-score (GSE131882)"); ax[3].set_ylabel("v-score (Lake 2023)")
    ax[3].text(.03,.95,f"r = {r:.2f}\nn = {len(mm)}",transform=ax[3].transAxes,va="top",fontsize=8)
    ax[3].set_title("D  Independent-cohort concordance")
save(fig,"fig1_signature.png")

# ---------- Figure 2: reversal screening ----------
PE=["empagliflozin","dapagliflozin","losartan","valsartan","irbesartan","telmisartan","enalapril","lisinopril","ramipril","perindopril","pirfenidone","nintedanib","acetylcysteine","metformin"]
NT=["cisplatin","gentamicin","amphotericin","vancomycin","tacrolimus","ibuprofen","naproxen","acyclovir","methotrexate"]
b2n=pd.read_csv(os.path.join(ROOT,"data/brd_to_name.tsv"),sep="\t",index_col=0,dtype=str)["pert_iname"].to_dict()
df=pd.read_csv(os.path.join(RES,"ranks_reversalinjured_normal.csv"),header=[0,1],index_col=0)
rc=("rank",df.columns.get_level_values(1)[0]); ranks=df[rc].astype(int); N=len(ranks)
name_of={b:b2n.get(b,b) for b in ranks.index}
def best(drug):
    h=[b for b,nm in name_of.items() if drug in nm.lower()]
    return int(min(ranks[b] for b in h)) if h else None
fig,ax=plt.subplots(1,3,figsize=(11,3.4))
pe=[(d,best(d)) for d in PE if best(d) is not None]; nt=[(d,best(d)) for d in NT if best(d) is not None]
ax[0].boxplot([[r for _,r in pe],[r for _,r in nt]],tick_labels=["Protective","Nephrotoxic"],showfliers=True,widths=.5)
ax[0].scatter([1]*len(pe),[r for _,r in pe],alpha=.7,color="#55A868"); ax[0].scatter([2]*len(nt),[r for _,r in nt],alpha=.7,color="#C44E52")
ax[0].set_ylabel("rank (lower = prioritized)"); ax[0].set_title("A  Known-drug recall")
pe_s=sorted(pe,key=lambda x:x[1]); ax[1].barh([d for d,_ in pe_s][::-1],[r for _,r in pe_s][::-1],color="#4C72B0")
ax[1].set_xlabel("rank / %d"%N); ax[1].set_title("B  Renoprotective drugs"); ax[1].tick_params(labelsize=7)
top=ranks.nsmallest(15); ax[2].barh([name_of[b] for b in top.index][::-1],top.values[::-1],color="#8172B3")
ax[2].set_xlabel("rank"); ax[2].set_title("C  Top-15 reversal hits"); ax[2].tick_params(labelsize=6)
save(fig,"fig2_screening.png")

# ---------- Figure 3: docking + enrichment ----------
dk=pd.read_csv(os.path.join(P3,"docking_results.csv")); dec=pd.read_csv(os.path.join(P3,"benchmark_decoys.csv"))
fig,ax=plt.subplots(1,3,figsize=(11,3.4))
colmap={"known":"#C44E52","hit":"#4C72B0","control":"#999999"}
for cat,g in dk.groupby("category"):
    ax[0].scatter(np.full(len(g),{"known":1,"hit":2,"control":3}[cat])+np.random.uniform(-.12,.12,len(g)),
                  g.affinity,color=colmap[cat],s=45,label=cat,alpha=.9)
ax[0].set_xticks([1,2,3]); ax[0].set_xticklabels(["Known\nALK5 inh.","Signature\nhits","Non-target\ncontrol"])
ax[0].set_ylabel("Vina affinity (kcal/mol)"); ax[0].set_title("A  ALK5 docking"); ax[0].legend(fontsize=6,frameon=False)
allc=pd.concat([dk[dk.category=="known"][["name","affinity"]].assign(label=1),dec[["name","affinity"]].assign(label=0)]).dropna().sort_values("affinity").reset_index(drop=True)
allc["rank"]=np.arange(1,len(allc)+1)
ax[1].plot(allc["rank"],-allc["affinity"],color="#999999",lw=1,zorder=1)
act=allc[allc.label==1]; ax[1].scatter(act["rank"],-act["affinity"],color="#C44E52",s=45,zorder=2,label="known ALK5 inh.")
ax[1].set_xlabel("rank among %d compounds"%len(allc)); ax[1].set_ylabel("-affinity"); ax[1].set_title("B  Actives recovery"); ax[1].legend(fontsize=6,frameon=False)
ef=[8.0,4.0,3.0]; ax[2].bar(["EF5%","EF10%","EF20%"],ef,color="#55A868"); ax[2].axhline(1,color="k",ls="--",lw=.8)
ax[2].set_ylabel("Enrichment factor"); ax[2].set_title("C  Enrichment (AUC 0.873)")
save(fig,"fig3_docking.png")

# ---------- Figure 4: pathways ----------
def top_terms(fn,k=10):
    p=os.path.join(P3,fn); d=pd.read_csv(p); d=d.sort_values("p.adjust").head(k)
    return d["Description"].tolist()[::-1], (-np.log10(d["p.adjust"])).tolist()[::-1]
fig,ax=plt.subplots(1,3,figsize=(13,3.8))
terms,v=top_terms("enrich_GO_BP_up.csv"); ax[0].barh(terms,v,color="#C44E52"); ax[0].set_title("A  GO BP up"); ax[0].tick_params(labelsize=6); ax[0].set_xlabel("-log10 padj")
terms,v=top_terms("enrich_KEGG_up.csv"); ax[1].barh(terms,v,color="#DD8452"); ax[1].set_title("B  KEGG up"); ax[1].tick_params(labelsize=6); ax[1].set_xlabel("-log10 padj")
terms,v=top_terms("enrich_GO_BP_down.csv"); ax[2].barh(terms,v,color="#4C72B0"); ax[2].set_title("C  GO BP down"); ax[2].tick_params(labelsize=6); ax[2].set_xlabel("-log10 padj")
save(fig,"fig4_pathways.png")

# ---------- Figure 5: mechanistic model schematic ----------
fig,ax=plt.subplots(figsize=(9,4.2)); ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,4.5)
def box(x,y,w,h,t,c):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=.06",fc=c,ec="k",lw=.8))
    ax.text(x+w/2,y+h/2,t,ha="center",va="center",fontsize=8,wrap=True)
box(.2,3.1,2.6,.9,"Healthy PT\n(transport, metabolism)","#DCE9F5")
box(.2,.5,2.6,.9,"Injured PT\n(dedifferentiation)","#FBE3D6")
box(3.4,1.8,3.0,1.4,"Signature\nFocal adhesion / actin cytoskeleton ↑\nTubular transport / metabolism ↓","#EFE5F5")
box(7.0,2.9,2.7,.9,"ROCK inhibitors\n(Y-27632, Y-39983, fasudil)","#DCE9F5")
box(7.0,1.2,2.7,.9,"ALK5 / TGFBR1 inhibitor\n(LY-364947)","#DCE9F5")
for a,b in [((2.8,3.55),(3.4,2.9)),((2.8,0.95),(3.4,2.1)),((6.4,2.7),(7.0,3.35)),((6.4,2.3),(7.0,1.65))]:
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle="->",mutation_scale=12,lw=1))
ax.text(5,4.35,"Signature-guided repurposing: pathway → target",ha="center",fontsize=10,weight="bold")
ax.text(8.35,.75,"→ renoprotection\n(literature-validated)",ha="center",fontsize=7,style="italic")
save(fig,"fig5_model.png")
print("ALL FIGURES DONE")
