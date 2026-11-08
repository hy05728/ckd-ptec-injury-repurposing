"""wf.py — 工作流示意图绘制（供 Figure 1A 与 Graphical Abstract 复用）"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# 每步一种浅色，边框用同色深版；打印与黑白复印均可区分
STEPS = [
    dict(
        title="Kidney scRNA-seq",
        lines=["GSE131882, 6 donors", "121,788 cells",
               "QC: 20,703 cells", "7,788 proximal tubule cells"],
        face="#DCE9F5", edge="#2B6CB0",
    ),
    dict(
        title="Injury signature",
        lines=["healthy vs injured PTEC", "v-score, per gene",
               "VCAM1 1.03  KRT8 0.77", "HAVCR1 0.34  LCN2 0.08"],
        face="#DFF0E3", edge="#2F855A",
    ),
    dict(
        title="Reversal screen",
        lines=["injured → healthy", "DrugReflector, 9,597 cmpds",
               "RAS inhibitors top 1.3%", "ROCK / ALK5 hits first"],
        face="#FBEED7", edge="#B7791F",
    ),
    dict(
        title="Structure validation",
        lines=["ALK5 (PDB 2WOU) docking", "60 decoys + 5 known ligands",
               "EF5% = 8.0", "ROC-AUC 0.873"],
        face="#EBE2F2", edge="#6B46C1",
    ),
    dict(
        title="Candidates",
        lines=["ROCK: Y-27632, Y-39983,", "           fasudil",
               "ALK5: LY-364947", "pathway: focal adhesion ↑"],
        face="#F6E0DD", edge="#9B2C2C",
    ),
]


def draw_workflow(ax, compact=False, steps=None):
    """在给定 ax 上画 5 步横向流程（ax 坐标 0-100）。"""
    steps = steps or STEPS
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    n = len(steps)
    gap = 3.0
    w = (100 - gap * (n - 1)) / n
    y0, h = 6, 88
    bar = 14                       # 标题条高度
    tsize = 8.2 if compact else 9.6
    bsize = 7.0 if compact else 7.8

    for i, s in enumerate(steps):
        x = i * (w + gap)
        ax.add_patch(FancyBboxPatch(
            (x, y0), w, h, boxstyle="round,pad=0.15,rounding_size=1.8",
            facecolor=s["face"], edgecolor=s["edge"], linewidth=1.1))
        # 标题条
        ax.add_patch(FancyBboxPatch(
            (x, y0 + h - bar), w, bar,
            boxstyle="round,pad=0.15,rounding_size=1.8",
            facecolor=s["edge"], edgecolor=s["edge"], linewidth=1.1, alpha=0.92))
        # 序号贴左，标题居中于剩余空间，避免长标题压住序号
        ax.text(x + 1.7, y0 + h - bar / 2, str(i + 1), ha="left", va="center",
                fontsize=tsize - 0.4, color="white", alpha=0.95, fontweight="bold")
        ax.text(x + w / 2 + 1.6, y0 + h - bar / 2, s["title"], ha="center", va="center",
                fontsize=tsize, color="white", fontweight="bold")
        # 正文
        nline = len(s["lines"])
        top = y0 + h - bar - 12
        step = 12 if nline > 2 else 14
        for j, ln in enumerate(s["lines"]):
            ax.text(x + w / 2, top - j * step, ln, ha="center", va="center",
                    fontsize=bsize, color="#1a202c")
        # 箭头
        if i < n - 1:
            ax.add_patch(FancyArrowPatch(
                (x + w + 0.3, y0 + h / 2), (x + w + gap - 0.3, y0 + h / 2),
                arrowstyle="-|>", mutation_scale=13, linewidth=1.3, color="#718096"))


def save_multi(fig, path_base, dpi_tiff=600):
    """一次输出 png/tif/pdf；TIFF 转 RGB（Elsevier 要求）。"""
    import os
    from PIL import Image
    fig.savefig(path_base + ".png", bbox_inches="tight", dpi=300)
    fig.savefig(path_base + ".tif", bbox_inches="tight", dpi=dpi_tiff,
                pil_kwargs={"compression": "tiff_lzw"})
    Image.open(path_base + ".tif").convert("RGB").save(
        path_base + ".tif", compression="tiff_lzw", dpi=(dpi_tiff, dpi_tiff))
    fig.savefig(path_base + ".pdf", bbox_inches="tight")
