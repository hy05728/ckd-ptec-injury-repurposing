import os, shutil, zipfile

ROOT = r"E:\workbuddy\2026-11-08-00-03-12"
SRC = os.path.join(ROOT, "submission_CBM")
FIG = os.path.join(ROOT, "figs")
OUT = os.path.join(ROOT, "final_submission")
PDF = os.path.join(OUT, "figures_vector_pdf")

if os.path.exists(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT); os.makedirs(PDF)

# 1) 文档类，按上传顺序重命名
docmap = [
    ("Cover_Letter.docx",             "01_Cover_Letter.docx"),
    ("Manuscript.docx",               "02_Manuscript.docx"),
    ("Highlights.docx",               "09_Highlights.docx"),
    ("Declarations.docx",             "10_Declarations.docx"),
    ("Title_Page.docx",               "11_Title_Page.docx"),
    ("Abstract_and_Keywords.docx",    "12_Abstract_Keywords.docx"),
]
for src, dst in docmap:
    p = os.path.join(SRC, src)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(OUT, dst))

# 2) 图件 TIFF（600 dpi RGB）
figs = [("fig1_signature", "03_Figure_1"), ("fig2_screening", "04_Figure_2"),
        ("fig3_docking", "05_Figure_3"), ("fig4_pathways", "06_Figure_4"),
        ("fig5_model", "07_Figure_5"), ("graphical_abstract", "08_Graphical_Abstract")]
for src, dst in figs:
    for ext in ("tif", "pdf"):
        p = os.path.join(FIG, f"{src}.{ext}")
        if os.path.exists(p):
            tgt = OUT if ext == "tif" else PDF
            shutil.copy2(p, os.path.join(tgt, f"{dst}.{ext}"))

# 3) 说明类
for f in ["Suggested_Reviewers.md", "投稿操作指引_CN.md"]:
    p = os.path.join(SRC, f)
    if os.path.exists(p):
        shutil.copy2(p, os.path.join(OUT, "13_" + f if f.startswith("Suggested") else f))

# 4) 清单
readme = """# CBM 投稿资料包 — 上传清单

期刊：**Computers in Biology and Medicine**（Elsevier, Editorial Manager）
标题：A single-cell injury signature of proximal tubules guides drug repurposing in chronic kidney disease

---

## 按 Editorial Manager 的 Item Type 逐项上传

| # | 文件 | Item Type（下拉选这个） | 说明 |
|---|---|---|---|
| 1 | `01_Cover_Letter.docx` | Cover Letter | 必须第一个传 |
| 2 | `02_Manuscript.docx` | Manuscript | 含 Title page、Abstract、Keywords、正文、Tables、Figure legends、References |
| 3-7 | `03_Figure_1.tif` ~ `07_Figure_5.tif` | Figure | 逐张上传，顺序填 1..5 |
| 8 | `08_Graphical_Abstract.tif` | Graphical Abstract | 必填项，尺寸 19.9 x 8.0 cm |
| 9 | `09_Highlights.docx` | Highlights | 5 条，每条 ≤85 字符 |
| 10 | `10_Declarations.docx` | Declaration of Interest | 含 Funding / Ethics / CRediT / Competing interests / Data & Code availability |
| — | `11_Title_Page.docx` | （备用） | 若系统要求单独作者页才传；否则 Manuscript 里已含 |
| — | `12_Abstract_Keywords.docx` | （备用） | 若系统要求单独摘要才传 |
| — | `13_Suggested_Reviewers.md` | Suggested Reviewers | 不是文件上传，是在系统里逐条填 5 位候选人的姓名/单位/邮箱 |
| — | `figures_vector_pdf/` | （备用） | 若编辑要求矢量图，用这里的 PDF 替代 TIFF |

---

## 上传前最后 4 项检查

1. **通读 `02_Manuscript.docx` 一遍**。机器自检能保证数字与引用一致（已跑过：引用 1-24 双向匹配、CJK 字符 0、图表全部被正文引用），但句子读起来顺不顺只能你判断。
2. **确认 GitHub 仓库已建并公开**：`https://github.com/hy05728/ckd-ptec-injury-repurposing`
   文章里写了这个 URL，审稿人会点。仓库没建的话，改投前先把 `repo_github/` 推上去。
3. **Cover Letter 里没有需要你改的地方**：已按扬大第一单位写死，逐字与 Title Page 一致。
4. **确认无基金、单作者**：两份文件里都已写明，与事实一致即可。

---

## 关键数据一览（审稿人若质疑，直接引这些）

| 环节 | 结果 |
|---|---|
| 签名构建 | 121,788 → 20,703 细胞；7,788 PTEC（6,269 healthy / 1,519 injured）|
| 损伤标志物 v-score | VCAM1 1.03、KRT8 0.77、HAVCR1 0.34、LCN2 0.08 |
| 独立队列复现 | Lake 2023，**Pearson r = 0.73 / Spearman 0.67**，21,772 基因 |
| 置换检验（特异性）| ROCK **P = 0.005**、ALK5 **P = 0.020**；RASi 不显著（0.29 / 0.22）|
| 无监督重定义 | 剔除 marker 基因重聚类，ROCK **P = 0.030**、ALK5 **P = 0.010**；签名相关 r = 0.80 |
| 对接基准 | 120 ChEMBL actives + 398 property-matched decoys；**ROC-AUC 0.904（95% CI 0.873-0.932）**；BEDROC 3.80 |
| 无关药类对照 | ROCK 11 / ALK5 14 vs 抗组胺 1713、β阻滞剂 1511、抗生素 2690、SSRI 4999 |

---

## 已知局限（都已写进 Discussion，审稿人问起不要否认）

- 无湿实验验证
- 未做 MD / MM-PBSA（Limitations 已明写"were not followed by molecular dynamics or free-energy calculations"）
- 个体化合物排位对状态定义敏感（Y-39983 从 2 漂到 496），只有类别层面稳健
- L1000 库不含 SGLT2i 与 finerenone
- EF 在小比例基准下上限为 4.32（文中已注明）
"""
open(os.path.join(OUT, "00_上传清单_README.md"), "w", encoding="utf-8").write(readme)

# 5) 打包 zip
zp = os.path.join(ROOT, "final_submission.zip")
with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
    for dp, dn, fn in os.walk(OUT):
        for f in sorted(fn):
            p = os.path.join(dp, f)
            z.write(p, os.path.relpath(p, OUT))

print("=== final_submission/ ===")
total = 0
for dp, dn, fn in os.walk(OUT):
    rel = os.path.relpath(dp, OUT)
    for f in sorted(fn):
        p = os.path.join(dp, f)
        sz = os.path.getsize(p)
        total += sz
        print(f"  {os.path.join(rel, f) if rel != '.' else f:<42} {sz/1048576:8.2f} MB")
print(f"\n  total {total/1048576:.1f} MB  ->  {zp} ({os.path.getsize(zp)/1048576:.1f} MB zipped)")
