import gc
import json
import os
import re
import sys
import warnings
from datetime import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co

BASE_GRN = os.environ.get(
    "BASE_GRN",
    os.path.join(os.path.expanduser("~"), "celloracle_data", "promoter_base_GRN",
                 "hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"))
H5AD = "03_celloracle_input.h5ad"
TF_FILE = "key_transcription_factors.txt"
OUT = "output"

CLUSTER_COL = "celltype_group"
KEY_CELLTYPE = "Stromal"
EMBEDDING = "X_umap"
N_COMPS = 40
N_JOBS = 10
KO_EXPRESSION = 0.0
OE_FACTOR = 1.5
N_PROPAGATION = 3
N_NEIGHBORS = 200
SIGMA_CORR = 0.05
MC_STEPS = 500
MC_DUP = 5
MC_SEED = 123

CT_COLORS = {"Stromal": "#C0392B", "Epithelial": "#F4A582", "Myeloid": "#7B3294"}

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 10
plt.rcParams["savefig.dpi"] = 300
warnings.filterwarnings("ignore")

os.makedirs(OUT, exist_ok=True)

def safe_name(x):
    return re.sub(r"[^\w\-.]+", "_", str(x))

def save_fig(fig, base, dpi=300):
    for ext in ("pdf", "png"):
        fig.savefig(f"{base}.{ext}", bbox_inches="tight", dpi=dpi)
    plt.close(fig)

def plot_tf_expression(oracle, tf, outdir):
    try:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sc.pl.umap(oracle.adata, color=[tf], layer="imputed_count", use_raw=False,
                   cmap="viridis", frameon=False, ax=axes[0], show=False,
                   title=f"{tf} (imputed)")
        sc.pl.umap(oracle.adata, color=CLUSTER_COL, frameon=False, ax=axes[1],
                   show=False, title="Cell types")
        fig.tight_layout()
        save_fig(fig, f"{outdir}/{tf}_expression_umap")

        fig2, ax = plt.subplots(figsize=(6, 5))
        sc.pl.umap(oracle.adata, color=[tf], layer="imputed_count", use_raw=False,
                   cmap="viridis", frameon=False, ax=ax, show=False,
                   title=f"{tf} expression (imputed)")
        fig2.tight_layout()
        save_fig(fig2, f"{outdir}/{tf}_expression_only")
    except Exception:
        pass

def plot_tf_violin(oracle, tf, outdir):
    try:
        fig, ax = plt.subplots(figsize=(6, 4))
        sc.pl.violin(oracle.adata, keys=[tf], groupby=CLUSTER_COL, layer="imputed_count",
                     use_raw=False, ax=ax, show=False, title=f"{tf} by cell type")
        fig.tight_layout()
        save_fig(fig, f"{outdir}/{tf}_violin")
    except Exception:
        pass

def export_targets(links, tf, outdir):
    n = 0
    try:
        for ct in links.cluster:
            if ct not in links.filtered_links:
                continue
            d = links.filtered_links[ct].copy()
            if "source" not in d.columns:
                continue
            d = d[d["source"] == tf]
            if not d.empty:
                d.to_csv(f"{outdir}/targets_{tf}_in_{safe_name(ct)}.csv", index=False)
                n += len(d)
    except Exception:
        pass

def quantify_effect(oracle, tf, label, expr_by_ct):
    ct = oracle.adata.obs[CLUSTER_COL].astype(str).values
    d_emb = np.asarray(oracle.delta_embedding)
    d_rnd = np.asarray(getattr(oracle, "delta_embedding_random", np.full_like(d_emb, np.nan)))
    mag = np.linalg.norm(d_emb, axis=1)
    mag_r = np.linalg.norm(d_rnd, axis=1)

    tr = getattr(oracle, "markvov_transition", None)
    tr_r = getattr(oracle, "markvov_transition_random", None)

    rows = []
    for c in sorted(set(ct)):
        m = ct == c
        row = {"TF": tf, "perturbation": label, "celltype": c, "n_cell": int(m.sum()),
               "mean_imputed_expression": float(expr_by_ct.get(c, np.nan)),
               "mean_shift": float(mag[m].mean()),
               "mean_shift_random": float(mag_r[m].mean()),
               "ratio_vs_random": float(mag[m].mean() / mag_r[m].mean()) if mag_r[m].mean() > 0 else np.nan}
        for tag, tab in (("perturbed", tr), ("random", tr_r)):
            if tab is not None and c in tab.index:
                tot = float(tab.loc[c].sum())
                row[f"self_retention_{tag}"] = float(tab.loc[c, c] / tot) if tot > 0 else np.nan
            else:
                row[f"self_retention_{tag}"] = np.nan
        row["self_retention_diff_pp"] = (row["self_retention_perturbed"] -
                                         row["self_retention_random"]) * 100
        rows.append(row)
    return pd.DataFrame(rows)

def run_perturbation(oracle, tf, value, label, outdir, expr_by_ct, scale=25, scale_sim=15):
    res = {"TF": tf, "perturbation": label, "success": False}
    try:
        oracle.simulate_shift(perturb_condition={tf: value}, n_propagation=N_PROPAGATION)
    except Exception as e:
        return res
    try:
        oracle.estimate_transition_prob(n_neighbors=N_NEIGHBORS, knn_random=True,
                                        sampled_fraction=1)
        oracle.calculate_embedding_shift(sigma_corr=SIGMA_CORR)
    except Exception as e:
        return res

    try:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
        oracle.plot_quiver(scale=scale, ax=ax[0])
        ax[0].set_title(f"{tf} {label}")
        oracle.plot_quiver_random(scale=scale, ax=ax[1])
        ax[1].set_title("Randomized control")
        for a in ax:
            a.set_xlabel("UMAP1")
            a.set_ylabel("UMAP2")
        fig.tight_layout()
        save_fig(fig, f"{outdir}/{tf}_{label}_quiver")
    except Exception:
        pass

    try:
        oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=N_NEIGHBORS)
        oracle.calculate_mass_filter(min_mass=0.01, plot=False)
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
        oracle.plot_simulation_flow_on_grid(scale=scale_sim, ax=ax[0])
        ax[0].set_title(f"{tf} {label}")
        oracle.plot_simulation_flow_random_on_grid(scale=scale_sim, ax=ax[1])
        ax[1].set_title("Randomized control")
        for a in ax:
            a.set_xlabel("UMAP1")
            a.set_ylabel("UMAP2")
        fig.tight_layout()
        save_fig(fig, f"{outdir}/{tf}_{label}_grid")

        fig, ax = plt.subplots(figsize=(5.5, 5.2))
        oracle.plot_cluster_whole(ax=ax, s=10)
        oracle.plot_simulation_flow_on_grid(scale=scale_sim, ax=ax, show_background=False)
        ax.set_title(f"{tf} {label}: simulated state transition")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        fig.tight_layout()
        save_fig(fig, f"{outdir}/{tf}_{label}_cluster_grid")
    except Exception:
        pass

    try:
        oracle.prepare_markov_simulation()
        oracle.run_markov_chain_simulation(n_steps=MC_STEPS, n_duplication=MC_DUP,
                                           seed=MC_SEED, calculate_randomized=True)
        oracle.get_markov_simulation_cell_transition_table(cluster_column_name=CLUSTER_COL,
                                                          return_df=False)
        tr = getattr(oracle, "markvov_transition", None)
        tr_r = getattr(oracle, "markvov_transition_random", None)
        if tr is not None:
            tr.to_csv(f"{outdir}/{tf}_{label}_transition_perturbed.csv")
        if tr_r is not None:
            tr_r.to_csv(f"{outdir}/{tf}_{label}_transition_random.csv")
        eff = quantify_effect(oracle, tf, label, expr_by_ct)
        eff.to_csv(f"{outdir}/{tf}_{label}_celltype_effect.csv", index=False)
        res["effect"] = eff
    except Exception:
        pass

    res["success"] = True
    return res

adata = sc.read_h5ad(H5AD)

if CLUSTER_COL not in adata.obs.columns:
    if "external_annotation" in adata.obs.columns:
        adata.obs[CLUSTER_COL] = adata.obs["external_annotation"].astype(str)
    else:
        sys.exit(f"h5ad 中没有 {CLUSTER_COL} 列")
adata.obs[CLUSTER_COL] = adata.obs[CLUSTER_COL].astype("category")

if EMBEDDING not in adata.obsm:
    for e in ("X_umap", "umap", "X_tsne"):
        if e in adata.obsm:
            EMBEDDING = e
            break

X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
n_int = bool(np.allclose(X, np.round(X)))
if not n_int:
    sys.exit("adata.X 不是原始 counts")

base_grn = pd.read_parquet(BASE_GRN)
base_tfs = set(base_grn.columns[2:].tolist())

key_tf = [l.strip() for l in open(TF_FILE) if l.strip()]
key_tf = [t for t in key_tf if t in adata.var_names and t in base_tfs]

if not key_tf:
    sys.exit("没有可用于扰动的核心 TF")

oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name=CLUSTER_COL,
                                   embedding_name=EMBEDDING)
oracle.import_TF_data(TF_info_matrix=base_grn)

oracle.perform_PCA()
cumvar = np.cumsum(oracle.pca.explained_variance_ratio_)
fig, ax = plt.subplots(figsize=(5.2, 4.2))
ax.plot(np.arange(1, min(101, len(cumvar) + 1)), cumvar[:100], lw=2)
ax.axvline(N_COMPS, color="black", linestyle="--", lw=1.2)
ax.set_xlabel("Principal components")
ax.set_ylabel("Cumulative explained variance")
ax.set_title(f"PCA cutoff = {N_COMPS}")
fig.tight_layout()
save_fig(fig, f"{OUT}/Figure_PCA_cutoff_PC{N_COMPS}")

n_cell = oracle.adata.shape[0]
n_comps = max(min(N_COMPS, oracle.adata.shape[0] - 1, oracle.adata.shape[1] - 1), 2)
k = max(min(int(0.025 * n_cell), n_cell - 1), 20)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True,
                      b_sight=min(k * 8, n_cell - 1),
                      b_maxl=min(k * 4, n_cell - 1), n_jobs=N_JOBS)

fig, ax = plt.subplots(figsize=(6, 5))
sc.pl.umap(oracle.adata, color=CLUSTER_COL, frameon=False, ax=ax, show=False)
ax.set_title("Cell types")
save_fig(fig, f"{OUT}/Figure_celltype_umap")

n_cluster = oracle.adata.obs[CLUSTER_COL].nunique()
links = oracle.get_links(cluster_name_for_GRN_unit=CLUSTER_COL, alpha=10, verbose_level=1)
links.filter_links(p=0.05, weight="coef_abs", threshold_number=10000)

oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

summary = []
effect_rows = []
failed = []

for idx, tf in enumerate(key_tf, 1):
    tf_dir = f"{OUT}/TF_{tf}"
    os.makedirs(tf_dir, exist_ok=True)

    expr = sc.get.obs_df(oracle.adata, keys=[tf], layer="imputed_count")
    mx, mn = float(expr[tf].max()), float(expr[tf].mean())
    oe_value = OE_FACTOR * mx
    ct_vals = oracle.adata.obs[CLUSTER_COL].astype(str).values
    expr_by_ct = pd.Series(expr[tf].values).groupby(ct_vals).mean().to_dict()

    expr.to_csv(f"{tf_dir}/imputed_expression_{tf}.csv")
    with open(f"{tf_dir}/expression_stats.json", "w") as f:
        json.dump({"TF": tf, "max_expression": mx, "mean_expression": mn,
                   "KO_value": KO_EXPRESSION, "OE_value": oe_value,
                   "mean_expression_by_celltype": expr_by_ct}, f, indent=2)

    plot_tf_expression(oracle, tf, tf_dir)
    plot_tf_violin(oracle, tf, tf_dir)
    export_targets(links, tf, tf_dir)

    ko = run_perturbation(oracle, tf, KO_EXPRESSION, "KO", tf_dir, expr_by_ct)
    oe = run_perturbation(oracle, tf, oe_value, "OE", tf_dir, expr_by_ct)

    for r in (ko, oe):
        if "effect" in r:
            effect_rows.append(r["effect"])

    row = {"TF": tf, "max_imputed_expression": mx, "mean_imputed_expression": mn,
           "OE_value": oe_value, "KO_success": ko["success"], "OE_success": oe["success"]}
    summary.append(row)
    pd.DataFrame([row]).to_csv(f"{tf_dir}/summary.csv", index=False)
    if not (ko["success"] and oe["success"]):
        failed.append(tf)

    gc.collect()

pd.DataFrame(summary).to_csv(f"{OUT}/TF_perturbation_summary.csv", index=False)

if effect_rows:
    eff_all = pd.concat(effect_rows, ignore_index=True)
    eff_all["is_key_cell"] = eff_all["celltype"] == KEY_CELLTYPE
    eff_all.to_csv(f"{OUT}/04_celltype_effect.csv", index=False)

    cts = [c for c in ["Stromal", "Epithelial", "Myeloid"] if c in set(eff_all["celltype"])]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    labels = sorted(eff_all["perturbation"].unique())
    x = np.arange(len(cts))
    w = 0.8 / max(len(labels), 1)
    for ax, metric, ylab in ((axes[0], "mean_shift", "Mean embedding shift"),
                             (axes[1], "ratio_vs_random", "Shift / randomized control"),
                             (axes[2], "self_retention_diff_pp", "Self-retention change (pp)")):
        for j, lab in enumerate(labels):
            sub = eff_all[eff_all["perturbation"] == lab]
            vals = [sub.loc[sub["celltype"] == c, metric].mean() for c in cts]
            ax.bar(x + j * w, vals, width=w, label=lab)
        ax.set_xticks(x + w * (len(labels) - 1) / 2)
        ax.set_xticklabels(cts)
        ax.set_ylabel(ylab)
        ax.axhline(0, color="grey", lw=0.6)
        if metric == "ratio_vs_random":
            ax.axhline(1, color="#2166AC", ls="--", lw=0.8)
        ax.legend(frameon=False, fontsize=7)
    fig.suptitle("Core TF perturbation effect by cell type (UC)")
    fig.tight_layout()
    save_fig(fig, f"{OUT}/Figure_7I_celltype_effect")

    for lab in labels:
        sub = eff_all[eff_all["perturbation"] == lab]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
        tfs = sorted(sub["TF"].unique())
        x = np.arange(len(tfs))
        w = 0.8 / max(len(cts), 1)
        for ax, metric, ylab in ((axes[0], "mean_shift", "Mean embedding shift"),
                                 (axes[1], "ratio_vs_random", "Shift / randomized control"),
                                 (axes[2], "self_retention_diff_pp", "Self-retention change (pp)")):
            for j, c in enumerate(cts):
                vals = [sub.loc[(sub["celltype"] == c) & (sub["TF"] == t), metric].mean() for t in tfs]
                ax.bar(x + j * w, vals, width=w, label=c, color=CT_COLORS.get(c))
            ax.set_xticks(x + w * (len(cts) - 1) / 2)
            ax.set_xticklabels(tfs)
            ax.set_ylabel(ylab)
            ax.axhline(0, color="grey", lw=0.6)
            if metric == "ratio_vs_random":
                ax.axhline(1, color="#2166AC", ls="--", lw=0.8)
            ax.legend(frameon=False, fontsize=7)
        fig.suptitle(f"{lab}: core TF perturbation effect by cell type (UC)")
        fig.tight_layout()
        save_fig(fig, f"{OUT}/Figure_7I_celltype_effect_{lab}")

    ko_tab = eff_all[eff_all["perturbation"] == "KO"]