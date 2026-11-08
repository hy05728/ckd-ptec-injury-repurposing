"""audit_prose.py — 按 academic-humanizer 的分类扫描正文 AI 痕迹（只审计，不改）"""
import re
from docx import Document

MP = r"E:\workbuddy\2026-11-08-00-03-12\submission_CBM\Manuscript.docx"
d = Document(MP)
paras = [p.text.strip() for p in d.paragraphs if p.text.strip()]

AI_WORDS = ["delve", "underscore", "intricate", "tapestry", "testament", "pivotal",
            "showcase", "foster", "leverage", "realm", "seamless", "serves as",
            "landscape", "robust framework", "groundbreaking", "cutting-edge",
            "paves the way", "sheds light", "bridges the gap", "opens new avenues",
            "of paramount importance", "it is worth noting", "it should be noted",
            "importantly,", "notably,", "in order to", "a wide range of",
            "extensive experiments", "comprehensive", "novel", "for the first time",
            "to the best of our knowledge", "furthermore", "moreover", "additionally"]
OVERCLAIM = ["demonstrate", "prove", "establish", "confirm", "guarantee", "significantly"]

print("=" * 70)
print("1) em-dash / en-dash（Layer 1 要求全部删除）")
print("=" * 70)
for i, t in enumerate(paras):
    for ch, nm in [("\u2014", "em-dash"), ("\u2013", "en-dash")]:
        if ch in t:
            print(f"  [{nm}] para {i}: ...{t[max(0,t.find(ch)-60):t.find(ch)+60]}...")

print("\n" + "=" * 70)
print("2) 超长句（>30 词 或 3+ 从句）")
print("=" * 70)
sents = []
for i, t in enumerate(paras):
    for s in re.split(r"(?<=[.!?])\s+", t):
        s = s.strip()
        if len(s.split()) >= 3:
            sents.append((i, s, len(s.split())))
long_s = [x for x in sents if x[2] > 30]
print(f"  句子总数 {len(sents)}，>30 词 {len(long_s)} 句"
      f"（占 {100*len(long_s)/max(len(sents),1):.0f}%）")
for i, s, n in sorted(long_s, key=lambda x: -x[2])[:8]:
    subs = len(re.findall(r"\b(which|that|while|whereas|with|although|because)\b", s, re.I))
    print(f"  [{n}w, {subs} sub] {s[:120]}...")

print("\n" + "=" * 70)
print("3) AI 词表命中")
print("=" * 70)
body = " ".join(paras).lower()
hits = {w: body.count(w) for w in AI_WORDS if body.count(w) > 0}
for w, c in sorted(hits.items(), key=lambda x: -x[1]):
    print(f"  {w:<28} {c}")

print("\n" + "=" * 70)
print("4) 过强动词（需核对是否都有数字支撑）")
print("=" * 70)
for w in OVERCLAIM:
    c = body.count(w)
    if c:
        print(f"  {w:<16} {c}")

print("\n" + "=" * 70)
print("5) 正文结构")
print("=" * 70)
print(f"  段落 {len(paras)}，总词数 {len(body.split())}")
