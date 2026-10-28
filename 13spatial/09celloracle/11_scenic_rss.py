import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pyscenic.plotting import plot_rss
from pyscenic.rss import regulon_specificity_scores

AUC = "03_aucell.csv"
CELLS = "../09_scenic_cells.csv"
TOP_N = 6

auc = pd.read_csv(AUC, index_col=0)
cells = pd.read_csv(CELLS)
if not set(auc.index).issubset(set(cells["cell"])):
    auc = auc.T
cells = cells.set_index("cell").loc[auc.index]
labels = cells["celltype_group"]

rss = regulon_specificity_scores(auc, labels)
rss.to_csv("04_rss.csv")

celltypes = list(rss.index)
fig, axes = plt.subplots(1, len(celltypes), figsize=(4 * len(celltypes), 4))
for ax, ct in zip(axes, celltypes):
    plot_rss(rss, ct, top_n=TOP_N, max_n=None, ax=ax)
    ax.set_title(ct)
    ax.set_xlabel("Regulons")
    ax.set_ylabel("Specificity score")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"05_rss_curves.{ext}", dpi=300)
plt.close(fig)

top = rss.loc["Stromal"].sort_values(ascending=False)
top.to_csv("06_rss_stromal_ranked.csv", header=["rss"])
