"""13_selfcheck.py — 投稿前程序化终检"""
import re
from docx import Document

MP = r"E:\workbuddy\2026-11-08-00-03-12\submission_CBM\Manuscript.docx"
d = Document(MP)
body = "\n".join(p.text for p in d.paragraphs)
tables = "\n".join(c.text for t in d.tables for r in t.rows for c in r.cells)

print("=" * 62)
print("1) 关键数字核对（与真实计算结果对照）")
print("=" * 62)
checks = [
    ("20,703 cells", True), ("7,788", True), ("6,269", True), ("1,519", True),
    ("9,597", True), ("irbesartan (position 54)", True),
    ("telmisartan (position 66)", True), ("enalapril (position 125)", True),
    ("median position was 939", True), ("median 2,091", True),
    ("2.1 × 10^-9", True), ("5.0 × 10^-9", True),
    ("ROC-AUC 0.904", True), ("0.873-0.932", True),
    ("4.32 at 1%", True), ("BEDROC", True),
    ("120 actives", True), ("398 property-matched", True),
    ("P = 0.030", True), ("P = 0.010", True),
    ("r = 0.80", True), ("62%", True),
    ("r = 0.73", True), ("rho = 0.67", True), ("21,772", True),
    ("hy05728", True),
]
for s, should in checks:
    ok = (s in body) == should
    print(f"  {'OK ' if ok else 'MISS'}  {s}")

print("\n  旧错误值必须消失：")
for bad in ["median position was 938", "2,090", "preferentially ranked", "[GitHub URL"]:
    gone = bad not in body and bad not in tables
    print(f"  {'OK ' if gone else 'STILL'}  {bad!r}")

print("\n" + "=" * 62)
print("2) 参考文献引用完整性")
print("=" * 62)
refs_sec = body.split("References")[-1]
listed = sorted({int(m) for m in re.findall(r"^\s*(\d{1,2})\.\s+[A-Z]", refs_sec, re.M)})
print(f"  参考文献列表编号: {listed}")
cite_text = body.split("References")[0]
cited = set()
for m in re.finditer(r"\[([\d,\-\s]+)\]", cite_text):
    for part in m.group(1).split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            if a.strip().isdigit() and b.strip().isdigit():
                cited.update(range(int(a), int(b) + 1))
        elif part.isdigit():
            cited.add(int(part))
cited = sorted(cited)
print(f"  正文引用编号: {cited}")
unused = [n for n in listed if n not in cited]
missing = [n for n in cited if n not in listed]
print(f"  列出但未引用: {unused or 'NONE'}")
print(f"  引用但未列出: {missing or 'NONE'}")

print("\n" + "=" * 62)
print("3) 图表引用")
print("=" * 62)
for name, pat in [("Fig. 1", r"Fig\. 1(?![0-9])"), ("Fig. 2", r"Fig\. 2(?![0-9])"),
                  ("Fig. 3", r"Fig\. 3(?![0-9])"), ("Fig. 4", r"Fig\. 4(?![0-9])"),
                  ("Fig. 5", r"Fig\. 5(?![0-9])"), ("Table 1", r"Table 1(?![0-9])"),
                  ("Table 2", r"Table 2(?![0-9])")]:
    n = len(re.findall(pat, body))
    flag = "OK " if n >= 1 else "MISS"
    print(f"  {flag} {name}: {n} 处")

print("\n" + "=" * 62)
print("4) 语言与单位污染")
print("=" * 62)
cjk = re.findall(r"[\u2e80-\u9fff\u3000-\u303f\uff00-\uffef]", body + tables)
print(f"  CJK 字符: {len(cjk)} {'OK' if not cjk else 'FOUND -> ' + ''.join(cjk[:20])}")
AFF = "Yangzhou University Medical College, 136 Jiangyang Middle Road, Yangzhou 225100, Jiangsu Province, China"
print(f"  单位串逐字一致: {'OK' if AFF in body else 'MISS'}")
for bad in ["Gaoyou", "高邮", "225600", "People's Hospital", "Department of Nephrology"]:
    if bad in body:
        print(f"  废弃串: {bad} FOUND")

print("\n" + "=" * 62)
print("5) 待替换占位符")
print("=" * 62)
ph = re.findall(r"YOUR-GITHUB-USERNAME", body + tables)
print(f"  GitHub 用户名占位: {len(ph)} 处（需替换）")

print("\n" + "=" * 62)
print(f"字数: {len(body.split())}")
