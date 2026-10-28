import os
import subprocess
import sys

import pandas as pd

SRC = "scFEA_src"
DATA = os.path.abspath(os.path.join(SRC, "data"))
RES = "output"
INPUT_CSV = "01_input_for_scfea.csv"
META_CSV = "01_cell_metadata.csv"
KEY_CELLTYPE = "Stromal"
TRAIN_EPOCH = 100
SC_IMPUTATION = "False"
PY = sys.executable

os.makedirs(RES, exist_ok=True)
os.makedirs(os.path.join(RES, "input"), exist_ok=True)

df = pd.read_csv(INPUT_CSV, index_col=0)
df.to_csv(os.path.join(RES, "input", "expr.csv"))

cmd = [
    PY, os.path.join(SRC, "src", "scFEA.py"),
    "--data_dir", DATA,
    "--input_dir", os.path.join(RES, "input"),
    "--res_dir", RES,
    "--test_file", "expr.csv",
    "--moduleGene_file", "module_gene_m168.csv",
    "--stoichiometry_matrix", "cmMat_c70_m168.csv",
    "--cName_file", "cName_c70_m168.csv",
    "--sc_imputation", SC_IMPUTATION,
    "--train_epoch", str(TRAIN_EPOCH),
    "--output_flux_file", "flux.csv",
    "--output_balance_file", "balance.csv",
]
subprocess.run(cmd, check=True)

flux_path = "flux.csv"
if not os.path.exists(flux_path):
    sys.exit("未生成 flux.csv")

flux = pd.read_csv(flux_path, index_col=0)

meta = pd.read_csv(META_CSV).set_index("cell")
common = [c for c in flux.index if c in meta.index]
flux = flux.loc[common]
meta = meta.loc[common]

mean_by_celltype = flux.groupby(meta["celltype"]).mean()
mean_by_celltype.T.to_csv(os.path.join(RES, "03_flux_mean_by_celltype.csv"))

key_cells = meta.index[meta["celltype"] == KEY_CELLTYPE]
other_cells = meta.index[meta["celltype"] != KEY_CELLTYPE]
diff = pd.DataFrame({
    "module": flux.columns,
    "mean_keycell": flux.loc[key_cells].mean(axis=0).values,
    "mean_other": flux.loc[other_cells].mean(axis=0).values,
})
diff["delta"] = diff["mean_keycell"] - diff["mean_other"]
diff = diff.sort_values("delta", ascending=False)
diff.to_csv(os.path.join(RES, "03_flux_keycell_vs_other.csv"), index=False)
