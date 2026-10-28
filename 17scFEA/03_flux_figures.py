import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INFO = "scFEA_src/data/Human_M168_information.symbols.csv"
MEAN_BY_TYPE = "output/03_flux_mean_by_celltype.csv"
DIFF = "output/03_flux_keycell_vs_other.csv"
FLUX = "flux.csv"
META = "01_cell_metadata.csv"
KEY_CELLTYPE = "Stromal"
TOP_N = 20
OUT = "output"

os.makedirs(OUT, exist_ok=True)

info = pd.read_csv(INFO)
info["reaction"] = info["Compound_IN_name"] + " -> " + info["Compound_OUT_name"]
ann = info.set_index("Module_id")[["reaction", "Compound_IN_name", "Compound_OUT_name", "Supermodule_id"]]
ann.index = "M_" + ann.index.astype(str)

diff = pd.read_csv(DIFF)
diff["reaction"] = diff["module"].map(ann["reaction"])
diff["supermodule"] = diff["module"].map(ann["Supermodule_id"])
diff.to_csv(f"{OUT}/04_keycell_module_top.csv", index=False)

mean_by_type = pd.read_csv(MEAN_BY_TYPE, index_col=0).T

top_modules = [m for m in diff.head(TOP_N)["module"].tolist() if m in mean_by_type.columns]
mat = mean_by_type[top_modules]

fig, axes = plt.subplots(1, 2, figsize=(13, 6),
                         gridspec_kw={"width_ratios": [1.25, 1]})

z = (mat - mat.mean(axis=0)) / mat.std(axis=0).replace(0, np.nan)
im = axes[0].imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
axes[0].set_xticks(range(len(top_modules)))
axes[0].set_xticklabels([f"{m}\n{ann.loc[m, 'reaction']}" for m in top_modules],
                        rotation=90, fontsize=5)
axes[0].set_yticks(range(mat.shape[0]))
axes[0].set_yticklabels(mat.index, fontsize=7)
axes[0].set_title("A  Module flux by cell type\n(top 20 Stromal-enriched modules, z-scored per module)",

                  fontsize=8, fontweight="bold")
fig.colorbar(im, ax=axes[0], shrink=0.7, label="z-score")

d_plot = diff.head(TOP_N).iloc[::-1]
axes[1].barh(range(len(d_plot)), d_plot["delta"], color="#C0392B", height=0.7)
axes[1].set_yticks(range(len(d_plot)))
axes[1].set_yticklabels([f"{m}  {r}" for m, r in zip(d_plot["module"], d_plot["reaction"])],
     fontsize=6)
axes[1].set_xlabel("Mean flux: Stromal - other cells", fontsize=8)
axes[1].set_title("B  Stromal-enriched metabolic modules", fontsize=8, fontweight="bold")
axes[1].tick_params(axis="x", labelsize=7)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/01_plate_7g_flux.{ext}", dpi=600, bbox_inches="tight")
plt.close(fig)

