import os
import pickle
import time

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(HERE, ".numba_cache"))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(HERE, ".mplcache"))

import matplotlib

matplotlib.use("Agg")
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import palantir
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import spearmanr, zscore

MTX = "01_palantir_counts.mtx"
GENES = "01_palantir_genes.txt"
BARCODES = "01_palantir_barcodes.txt"
META = "01_palantir_metadata.csv"
OUT = "output"

SECTIONS = None
KEY_CELLTYPE = "Stromal"
KEY_ONLY = True
KEY_GENES = ["ADAMTS1", "FBN1", "SPARC", "VCAM1"]
REST_MARKERS = ["RGS5", "NOTCH3", "PDGFRB"]
ACT_MARKERS = ["COL1A1", "COL1A2", "MMP1", "MMP3", "IL11", "POSTN", "INHBA"]
N_HVG = 1500
N_PCS = 150
N_NEIGHBORS = 30
N_DIFF_COMP = 5
NUM_WAYPOINTS = 300
N_JOBS = 1
DO_MAGIC = True

RUN_TAG = ("all" if SECTIONS is None else "-".join(SECTIONS)) + ("_key" if KEY_ONLY else "_allspots")
CKPT = os.path.join("checkpoint", RUN_TAG)

COL_UC, COL_NC = "#C0392B", "#2166AC"
PT_CMAP = LinearSegmentedColormap.from_list("pt", ["#2166AC", "#FFFFFF", "#B2182B"])
FATE_CMAP = LinearSegmentedColormap.from_list("fate", ["#FFFFFF", "#7B3294"])
FS_PANEL, FS_TITLE, FS_TICK = 11, 8, 6

def marker_handles():
    start = plt.Line2D([], [], marker="*", color="none", markeredgecolor="black",
                       markersize=8, linestyle="none")
    term = plt.Line2D([], [], marker="^", color="black", markersize=4, linestyle="none")
    return [start, term], ["start spot (resting end)", "terminal state"]

os.makedirs(OUT, exist_ok=True)
os.makedirs(CKPT, exist_ok=True)

class step:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *a):
        pass

def cached(tag, fn):
    path = os.path.join(CKPT, f"{tag}.pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    with step(tag):
        val = fn()
    with open(path, "wb") as f:
        pickle.dump(val, f)
    return val

def build_adata():
    adata = sc.read_mtx(MTX).T
    adata.var_names = [g.strip() for g in open(GENES)]
    adata.obs_names = [b.strip() for b in open(BARCODES)]
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    meta = pd.read_csv(META).set_index("cell").loc[adata.obs_names]
    for c in meta.columns:
        adata.obs[c] = meta[c].values
    return adata

def preprocess(a):
    sc.pp.normalize_per_cell(a)
    sc.pp.log1p(a)
    a.layers["log_norm"] = a.X.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=N_HVG)
    n_pcs = min(N_PCS, a.n_obs - 1)
    sc.pp.pca(a, n_comps=n_pcs)
    sc.pp.neighbors(a, n_neighbors=min(N_NEIGHBORS, a.n_obs - 1))
    return a

def rest_act_axis(a):
    key = "MAGIC_imputed_data" if "MAGIC_imputed_data" in a.layers else "log_norm"
    expr = a.layers[key]
    expr = np.asarray(expr.todense()) if hasattr(expr, "todense") else np.asarray(expr)
    rest = [g for g in REST_MARKERS if g in a.var_names]
    act = [g for g in ACT_MARKERS if g in a.var_names]
    rest_score = expr[:, [a.var_names.get_loc(g) for g in rest]].mean(1)
    act_score = expr[:, [a.var_names.get_loc(g) for g in act]].mean(1)
    axis = pd.Series(zscore(rest_score) - zscore(act_score), index=a.obs_names)
    return axis, rest_score, rest, act, key

def run_section(sec, a):
    a = cached(f"{sec}_02_preprocessed", lambda: preprocess(a))
    dm_res = cached(f"{sec}_03_diffusion", lambda: palantir.utils.run_diffusion_maps(
        pd.DataFrame(a.obsm["X_pca"], index=a.obs_names), n_components=N_DIFF_COMP))
    ms_data = cached(f"{sec}_04_multiscale",
                     lambda: palantir.utils.determine_multiscale_space(dm_res))
    if DO_MAGIC:
        palantir.utils.run_magic_imputation(a, dm_res, n_jobs=N_JOBS)
    axis_score, rest_score, rest_use, act_use, axis_src = rest_act_axis(a)
    start_cell = axis_score.idxmax()

    def _pt():
        return palantir.core.run_palantir(ms_data, start_cell,
                                          num_waypoints=min(NUM_WAYPOINTS, a.n_obs), n_jobs=N_JOBS)

    pr_res = cached(f"{sec}_05_palantir", _pt)
    a.obs["pseudotime"] = pr_res.pseudotime.reindex(a.obs_names).values

    rho_rest = spearmanr(rest_score, a.obs["pseudotime"])[0]
    if rho_rest > 0:
        print(f"[警告] {sec} 的静息评分随拟时序升高，请核对起点标志；保留 Palantir 起点方向，不自动翻转。")

    a.obsm["X_umap"] = cached(f"{sec}_06_umap", lambda: (sc.tl.umap(a), a.obsm["X_umap"])[1])
    return a, pr_res, dm_res, rho_rest, start_cell

def keygene_matrix(sec, a):
    expr_source = "MAGIC_imputed_data" if "MAGIC_imputed_data" in a.layers else "log_norm"
    genes_use = [g for g in KEY_GENES if g in a.var_names]
    mat = a.layers[expr_source]
    mat = np.asarray(mat.todense()) if hasattr(mat, "todense") else np.asarray(mat)
    mat = mat[:, [a.var_names.get_loc(g) for g in genes_use]]
    if expr_source == "MAGIC_imputed_data":
        del a.layers["MAGIC_imputed_data"]
    return mat, genes_use, expr_source

adata_full = cached("01_adata", build_adata)
sections = sorted(adata_full.obs["section"].astype(str).unique(),
                  key=lambda s: (adata_full.obs.loc[adata_full.obs["section"].astype(str) == s, "group"].iloc[0] != "UC",
                                 s))
if SECTIONS is not None:
    sections = [s for s in sections if s in SECTIONS]

adata_sec, pt_sec, kg_sec, rho_sec, rho_depth_sec, start_sec, term_sec = {}, {}, {}, {}, {}, {}, {}
genes_use, expr_source = None, None
for sec in sections:
    a = adata_full[adata_full.obs["section"].astype(str) == sec].copy()
    if KEY_ONLY:
        a = a[a.obs["celltype"].astype(str) == KEY_CELLTYPE].copy()
    a, pr_res, dm_res, rho_rest, start_cell = run_section(sec, a)
    mat, genes_use, expr_source = keygene_matrix(sec, a)
    rho_depth_sec[sec] = spearmanr(a.obs["pseudotime"], a.obs["nCount"])[0]
    adata_sec[sec] = a
    pt_sec[sec] = pr_res
    kg_sec[sec] = pd.DataFrame(mat, index=a.obs_names, columns=genes_use)
    rho_sec[sec] = rho_rest
    start_sec[sec] = start_cell
    term_sec[sec] = list(pr_res.branch_probs.columns)

obs_all = pd.concat([adata_sec[s].obs[["section", "group", "celltype", "x", "y", "nCount", "nFeature", "pseudotime"]]
                     for s in sections])
obs_all.rename_axis("cell").to_csv(f"{OUT}/04_pseudotime.csv")
pd.concat([pt_sec[s].branch_probs.assign(section=s, cell=pt_sec[s].branch_probs.index) for s in sections]) \
    .to_csv(f"{OUT}/05_branch_probs.csv", index=False)
pd.concat([kg_sec[s].assign(section=s, cell=kg_sec[s].index) for s in sections]) \
    .reset_index(drop=True).to_csv(f"{OUT}/06_keygene_magic.csv", index=False)

n_spot = {s: int(adata_sec[s].shape[0]) for s in sections}

n_col, n_row = 3, int(np.ceil(len(sections) / 3))
fig, axes = plt.subplots(n_row, n_col, figsize=(6.1, 2.1 * n_row))
axes = np.atleast_1d(axes).ravel()
for i, sec in enumerate(sections):
    o = adata_sec[sec].obs
    um = adata_sec[sec].obsm["X_umap"]
    sca = axes[i].scatter(um[:, 0], um[:, 1], c=o["pseudotime"], s=2, cmap=PT_CMAP,
                          vmin=0, vmax=1, linewidths=0)
    k = adata_sec[sec].obs_names.get_loc(start_sec[sec])
    axes[i].scatter(um[k, 0], um[k, 1], marker="*", s=55, facecolor="none", edgecolor="black",
                    linewidths=0.9, zorder=5)
    for t in term_sec[sec]:
        j = adata_sec[sec].obs_names.get_loc(t)
        axes[i].scatter(um[j, 0], um[j, 1], marker="^", s=16, color="black", linewidths=0, zorder=5)
    axes[i].set_title(f"{sec} ({o['group'].iloc[0]})", fontsize=FS_TITLE)
    axes[i].text(0.02, 0.98, chr(65 + i), transform=axes[i].transAxes, fontsize=FS_PANEL,
                 fontweight="bold", va="top")
    axes[i].set_xticks([]); axes[i].set_yticks([])
    for sp in axes[i].spines.values():
        sp.set_visible(False)
for j in range(len(sections), len(axes)):
    axes[j].axis("off")
fig.tight_layout(rect=[0, 0.045, 0.80, 1])
fig.legend(*marker_handles(), loc="lower center", ncol=2, frameon=False, fontsize=FS_TICK)
cax = fig.add_axes([0.83, 0.18, 0.013, 0.64])
cb = fig.colorbar(sca, cax=cax)
cb.set_label("Palantir pseudotime\n(0 = resting, 1 = activated)", fontsize=FS_TICK)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/01_umap_pseudotime.{ext}", dpi=600)
plt.close(fig)

fig, axes = plt.subplots(n_row, n_col, figsize=(6.1, 2.2 * n_row))
axes = np.atleast_1d(axes).ravel()
for i, sec in enumerate(sections):
    o = adata_sec[sec].obs
    sca = axes[i].scatter(o["x"],
               o["y"],
               c=o["pseudotime"],
               s=3,
               cmap=PT_CMAP,
               vmin=0,
               vmax=1,
               linewidths=0)
    s0 = o.loc[start_sec[sec]]
    axes[i].scatter(s0["x"], s0["y"], marker="*", s=55, facecolor="none", edgecolor="black",
                    linewidths=0.9, zorder=5)
    for t in term_sec[sec]:
        t0 = o.loc[t]
        axes[i].scatter(t0["x"], t0["y"], marker="^", s=16, color="black", linewidths=0, zorder=5)
    axes[i].set_title(f"{sec} ({o['group'].iloc[0]})", fontsize=FS_TITLE)
    axes[i].text(0.02, 0.98, chr(65 + i), transform=axes[i].transAxes, fontsize=FS_PANEL,
                 fontweight="bold", va="top")
    axes[i].invert_yaxis()
    axes[i].set_aspect("equal")
    axes[i].axis("off")
for j in range(len(sections), len(axes)):
    axes[j].axis("off")
fig.tight_layout(rect=[0, 0.045, 0.80, 1])
fig.legend(*marker_handles(), loc="lower center", ncol=2, frameon=False, fontsize=FS_TICK)
cax = fig.add_axes([0.83, 0.18, 0.013, 0.64])
cb = fig.colorbar(sca, cax=cax)
cb.set_label("Palantir pseudotime\n(0 = resting, 1 = activated)", fontsize=FS_TICK)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/02_spatial_pseudotime.{ext}", dpi=600)
plt.close(fig)

fate_secs = [s for s in sections if len(term_sec[s]) >= 2]
fate_panels = [(s, t) for s in fate_secs for t in term_sec[s]]
if len(fate_panels) > 0:
    f_col, f_row = 2, len(fate_secs)
    fig, axes = plt.subplots(f_row, f_col, figsize=(6.1, 2.0 * f_row))
    axes = np.atleast_1d(axes).ravel()
    for i, (sec, term) in enumerate(fate_panels):
        o = adata_sec[sec].obs
        k = term_sec[sec].index(term) + 1
        sca = axes[i].scatter(o["x"], o["y"], c=pt_sec[sec].branch_probs[term], s=3,
                              cmap=FATE_CMAP, vmin=0, vmax=1, linewidths=0)
        s0 = o.loc[start_sec[sec]]
        axes[i].scatter(s0["x"], s0["y"], marker="*", s=55, facecolor="none", edgecolor="black",
                        linewidths=0.9, zorder=5)
        axes[i].scatter(o.loc[term, "x"], o.loc[term, "y"], marker="^", s=22, color="black",
                        linewidths=0, zorder=5)
        axes[i].set_title(f"{sec} ({o['group'].iloc[0]})  fate {k}", fontsize=FS_TITLE)
        axes[i].text(0.02, 0.98, chr(65 + i), transform=axes[i].transAxes, fontsize=FS_PANEL,
                     fontweight="bold", va="top")
        axes[i].invert_yaxis()
        axes[i].set_aspect("equal")
        axes[i].axis("off")
    for j in range(len(fate_panels), len(axes)):
        axes[j].axis("off")
    fig.tight_layout(rect=[0, 0.05, 0.80, 1])
    fig.legend(*marker_handles(), loc="lower center", ncol=2, frameon=False, fontsize=FS_TICK)
    cax = fig.add_axes([0.83, 0.20, 0.013, 0.60])
    cb = fig.colorbar(sca, cax=cax)
    cb.set_label("Fate probability", fontsize=FS_TICK)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/07_fate_probability.{ext}", dpi=600)
    plt.close(fig)

fig, axes = plt.subplots(1, len(genes_use), figsize=(6.1, 1.9), sharex=True)
axes = np.atleast_1d(axes)
for k, g in enumerate(genes_use):
    for sec in sections:
        o = adata_sec[sec].obs
        order = np.argsort(o["pseudotime"].values)
        y = kg_sec[sec][g].values[order]
        col = COL_UC if o["group"].iloc[0] == "UC" else COL_NC
        axes[k].plot(o["pseudotime"].values[order],
                     pd.Series(y).rolling(50, min_periods=10, center=True).mean(),
                     color=col, linewidth=0.7, alpha=0.8)
    axes[k].set_title(g, fontsize=FS_TITLE)
    axes[k].set_xlabel("pseudotime", fontsize=FS_TITLE)
    axes[k].tick_params(labelsize=FS_TICK)
    axes[k].text(0.03, 0.95, chr(65 + k), transform=axes[k].transAxes, fontsize=FS_PANEL,
                 fontweight="bold", va="top")
axes[0].set_ylabel(f"{expr_source}\nexpression", fontsize=FS_TITLE)
handles = [plt.Line2D([], [], color=COL_UC, linewidth=1),
           plt.Line2D([], [], color=COL_NC, linewidth=1)]
fig.legend(handles, ["UC", "Normal"], loc="lower center", ncol=2, frameon=False, fontsize=FS_TICK)
fig.tight_layout(rect=[0, 0.07, 1, 1])
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/03_keygene_pseudotime.{ext}", dpi=600)
plt.close(fig)


adata_out = ad.concat([adata_sec[s] for s in sections],
                      join="inner",
                      label="section_id",
                      keys=sections)
adata_out.write_h5ad(f"{OUT}/adata_palantir.h5ad")
