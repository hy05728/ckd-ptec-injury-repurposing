#!/usr/bin/env python3
"""16_make_graphical_abstract.py — CBM 投稿必需的 Graphical Abstract（20 x 8 cm，合规）"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wf import draw_workflow, save_multi

FIG = r"E:/workbuddy/2026-11-08-00-03-12/figs"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})

# 简化文案：GA 在 13 x 5 cm 显示尺寸下必须可读，每行 <= 16 字符
STEPS_GA = [
    dict(title="scRNA-seq", lines=["GSE131882", "7,788 PTEC"], face="#DCE9F5", edge="#2B6CB0"),
    dict(title="Signature", lines=["healthy vs injured", "v-score"], face="#DFF0E3", edge="#2F855A"),
    dict(title="Reversal screen", lines=["9,597 compounds", "ROCK/ALK5 on top"], face="#FBEED7", edge="#B7791F"),
    dict(title="Docking", lines=["ALK5, 60 decoys", "ROC-AUC 0.87"], face="#EBE2F2", edge="#6B46C1"),
    dict(title="Candidates", lines=["ROCK inhibitors", "ALK5 inhibitors"], face="#F6E0DD", edge="#9B2C2C"),
]

fig = plt.figure(figsize=(7.87, 3.15))          # 20 x 8 cm
ax = fig.add_axes([0.015, 0.03, 0.97, 0.78])
draw_workflow(ax, compact=True, steps=STEPS_GA)

fig.text(0.5, 0.93,
         "Injured-tubule signature guides drug repurposing in chronic kidney disease",
         ha="center", va="center", fontsize=9.2, fontweight="bold", color="#1a202c")
fig.text(0.5, 0.855,
         "signature from 7,788 proximal tubule cells  →  reversal screen  →  docking benchmark",
         ha="center", va="center", fontsize=6.0, color="#4a5568")

save_multi(fig, os.path.join(FIG, "graphical_abstract"))
plt.close(fig)

from PIL import Image
im = Image.open(os.path.join(FIG, "graphical_abstract.png"))
print(f"graphical abstract: {im.size[0]}x{im.size[1]} px (requirement: >=1328x531)")
print(f"physical size at 300 dpi: {im.size[0]/300*2.54:.1f} x {im.size[1]/300*2.54:.1f} cm")
