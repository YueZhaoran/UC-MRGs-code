import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = "output"
EFFECT = f"{OUT}/04_celltype_effect.csv"
PLATE = f"{OUT}/05_plate_7i"

CT_COLORS = {"Stromal": "#C0392B", "Epithelial": "#F4A582", "Myeloid": "#7B3294"}
CT_ORDER = ["Stromal", "Epithelial", "Myeloid"]
KEY_CELLTYPE = "Stromal"
FS_PANEL, FS_TITLE, FS_TICK = 11, 7, 6

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "sans"

eff = pd.read_csv(EFFECT)
eff["self_retention_diff_pp"] = (eff["self_retention_perturbed"] -
                                 eff["self_retention_random"]) * 100

ko = eff[eff["perturbation"] == "KO"].copy()
tfs = sorted(ko["TF"].unique())
cts = [c for c in CT_ORDER if c in set(ko["celltype"])]

panels = [("mean_imputed_expression", "Imputed expression", False),
          ("ratio_vs_random", "Shift / randomized control", True),
          ("self_retention_diff_pp", "Self-retention change (pp)", False)]

fig, axes = plt.subplots(1, 3, figsize=(6.1, 2.35))
x = np.arange(len(tfs))
w = 0.8 / len(cts)
for ax, (metric, ylab, ref1) in zip(axes, panels):
    for j, c in enumerate(cts):
        vals = [ko.loc[(ko["TF"] == t) & (ko["celltype"] == c), metric].mean() for t in tfs]
        ax.bar(x + j * w, vals, width=w, color=CT_COLORS[c], label=c,
               edgecolor="none")
    ax.set_xticks(x + w * (len(cts) - 1) / 2)
    ax.set_xticklabels(tfs, fontsize=FS_TICK, rotation=45, ha="right")
    ax.set_ylabel(ylab, fontsize=FS_TITLE)
    ax.tick_params(labelsize=FS_TICK)
    ax.axhline(1 if ref1 else 0, color="#2166AC", ls="--", lw=0.7)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
for k, ax in enumerate(axes):
    ax.text(0.02, 0.98, chr(65 + k), transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", va="top")
axes[0].legend(fontsize=FS_TICK, frameon=False, loc="upper right", ncol=1,
               handlelength=1.0, borderpad=0.2, labelspacing=0.25)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{PLATE}.{ext}", dpi=600)
plt.close(fig)

