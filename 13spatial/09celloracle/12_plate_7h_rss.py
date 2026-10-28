import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCENIC = "scenic"
OUT = "output"
CORE = f"{SCENIC}/03_core_TFs.csv"
RSS = f"{SCENIC}/04_rss.csv"

KEY_CELLTYPE = "Stromal"
CT_ORDER = [KEY_CELLTYPE, "Epithelial", "Myeloid", "T cell", "B/Plasma", "NK/ILC"]
COL_KEY, COL_OTHER, COL_NS = "#C0392B", "#2166AC", "#BDBDBD"
FS_PANEL, FS_TITLE, FS_TICK, FS_LAB = 11, 7, 6, 5.5

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "sans"

rss = pd.read_csv(RSS, index_col=0)
rss.index = [str(i) for i in rss.index]
rss.columns = [str(c).split("(")[0].strip() for c in rss.columns]
core = pd.read_csv(CORE)
core_tfs = core["TF"].tolist()
rss_core = rss.loc[:, ~rss.columns.duplicated(keep="first")]
present_core = [t for t in core_tfs if t in rss_core.columns]

n_reg = rss_core.shape[1]

def draw_rss_curve(ax, ct, n_label=6, key=False, title=True):
    s = rss_core.loc[ct].sort_values(ascending=False)
    col = COL_KEY if key else COL_OTHER
    ax.plot(np.arange(1, len(s) + 1), s.values, color=col, lw=1.1 if key else 0.9)
    ax.set_xlim(0, len(s) + 1)
    lo = float(s.min()) - 0.02 * (float(s.max()) - float(s.min()))
    hi = float(s.max()) + 0.10 * (float(s.max()) - float(s.min()))
    ax.set_ylim(lo, hi)
    top = s.head(n_label)
    x_lab = 0.46 * len(s)
    y_lab = np.linspace(float(top.max()), lo + 0.30 * (hi - lo), n_label)
    for k, (tf, v) in enumerate(top.items()):
        ax.annotate(tf, xy=(k + 1, v), xytext=(x_lab, y_lab[k]), fontsize=FS_LAB,
                    va="center", ha="left", color="#111111",
                    fontweight="bold" if tf in present_core else "normal",
                    arrowprops=dict(arrowstyle="-", color="#999999", lw=0.4,
                                    shrinkA=0, shrinkB=1.5))
    ax.set_xlabel("Regulons", fontsize=FS_TITLE)
    ax.set_ylabel("Specificity score", fontsize=FS_TITLE)
    if title:
        ax.set_title(ct, fontsize=FS_TITLE)
    ax.tick_params(labelsize=FS_TICK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return s

fig = plt.figure(figsize=(6.1, 5.4))
gs = fig.add_gridspec(2, 5, height_ratios=[1.25, 1.0], hspace=0.55, wspace=0.75)

axA = fig.add_subplot(gs[0, 0:3])
s_key = draw_rss_curve(axA, KEY_CELLTYPE, n_label=8, key=True, title=False)
axA.text(0.0, 1.14, "A", transform=axA.transAxes, fontsize=FS_PANEL, fontweight="bold", va="bottom")
axA.text(0.10, 1.14, f"{KEY_CELLTYPE}: regulons ranked by specificity",
         transform=axA.transAxes, fontsize=FS_TITLE, va="bottom")

axB = fig.add_subplot(gs[0, 3:5])
mat = rss_core.loc[CT_ORDER, present_core].T
im = axB.imshow(mat.values, cmap="YlOrRd", aspect="auto", vmin=float(mat.values.min()),
                vmax=float(mat.values.max()))
axB.set_xticks(range(len(CT_ORDER)))
axB.set_xticklabels([c if len(c) <= 8 else c[:8] for c in CT_ORDER], rotation=45,
                    ha="right", fontsize=FS_LAB)
axB.set_yticks(range(len(present_core)))
axB.set_yticklabels(present_core, fontsize=FS_LAB)
for i in range(len(present_core)):
    for j in range(len(CT_ORDER)):
        axB.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center", fontsize=4.5,
                 color="#111111")
axB.set_title("Core TF RSS by cell type", fontsize=FS_TITLE)
axB.text(0.0, 1.14, "B", transform=axB.transAxes, fontsize=FS_PANEL, fontweight="bold", va="bottom")
cb = fig.colorbar(im, ax=axB, fraction=0.05, pad=0.02)
cb.ax.tick_params(labelsize=FS_TICK)
cb.set_label("RSS", fontsize=FS_TICK)

others = [c for c in CT_ORDER if c != KEY_CELLTYPE]
axesC = []
for i, ct in enumerate(others[:5]):
    ax = fig.add_subplot(gs[1, i])
    draw_rss_curve(ax, ct, n_label=3, key=False)
    axesC.append(ax)
axesC[0].text(0.0, 1.22, "C", transform=axesC[0].transAxes, fontsize=FS_PANEL,
              fontweight="bold", va="bottom")

fig.savefig(f"{OUT}/06_plate_7h.pdf", dpi=600, bbox_inches="tight")
fig.savefig(f"{OUT}/06_plate_7h.png", dpi=600, bbox_inches="tight")
plt.close(fig)

