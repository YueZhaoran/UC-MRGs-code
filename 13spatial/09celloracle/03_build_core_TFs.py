import os
import sys

import anndata
import numpy as np
import pandas as pd

BASE_GRN = os.environ.get(
    "BASE_GRN",
    os.path.join(os.path.expanduser("~"), "celloracle_data", "promoter_base_GRN",
                 "hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet"))
H5AD = "03_celloracle_input.h5ad"
RSS = os.environ.get("SCENIC_RSS", "scenic/04_rss.csv")
OUT = "output"

KEY_CELLTYPE = "Stromal"
CELLTYPE_COL = "celltype_group"
N_TOP = 5
NOT_TF = ["CD59"]

os.makedirs(OUT, exist_ok=True)

adata = anndata.read_h5ad(H5AD)
genes = set(str(g).upper() for g in adata.var_names)

base_grn = pd.read_parquet(BASE_GRN)
base_tfs = set(t.upper() for t in base_grn.columns if t not in ("peak_id", "gene_short_name"))

if not os.path.exists(RSS):
    sys.exit(f"缺少 SCENIC RSS 结果：{RSS}（先跑 10_scenic_chain.sh 与 11_scenic_rss.py）")
rss = pd.read_csv(RSS, index_col=0)
if KEY_CELLTYPE not in rss.index:
    sys.exit(f"RSS 矩阵里没有 {KEY_CELLTYPE}")
ranked = rss.loc[KEY_CELLTYPE].sort_values(ascending=False)
ranked.index = [str(i).split("(")[0].strip() for i in ranked.index]
ranked = ranked[~ranked.index.duplicated(keep="first")]

rows = []
for rank, (tf, val) in enumerate(ranked.items(), 1):
    tf_up = tf.upper()
    rows.append({"TF": tf_up, "rss": float(val), "rss_rank": rank,
                 "in_base_grn": tf_up in base_tfs,
                 "in_expr": tf_up in genes,
                 "is_tf": tf_up not in NOT_TF})
ranked_df = pd.DataFrame(rows)
ranked_df.to_csv(f"{OUT}/03_tf_rss_ranked.csv", index=False)

usable = ranked_df[ranked_df["in_base_grn"] & ranked_df["in_expr"] & ranked_df["is_tf"]]

if usable.empty:
    sys.exit("没有可用的核心 TF")

core = usable.head(N_TOP).copy()
core_tfs = core["TF"].tolist()

core.to_csv(f"{OUT}/03_core_TFs.csv", index=False)
with open(f"{OUT}/03_key_transcription_factors.txt", "w") as f:
    f.write("\n".join(core_tfs) + "\n")
