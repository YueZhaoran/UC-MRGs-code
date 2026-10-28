#!/usr/bin/env bash
# 方案第七节(5)：CellOracle 链式流程（服务器执行）
#
# 用法：nohup bash 06_celloracle_chain.sh > celloracle_chain.log 2>&1 &
# 可覆盖的环境变量：BASE_GRN、CONDA_ENV、CELLORACLE_PYTHON（直接指定解释器）
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p output

BASE_GRN="${BASE_GRN:-$HOME/celloracle_data/promoter_base_GRN/hg19_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10_20210630.parquet}"
CONDA_ENV="${CONDA_ENV:-celloracle}"
CONDA_SH="${CONDA_SH:-$HOME/anaconda3/etc/profile.d/conda.sh}"

step() { echo; echo "############ [$(date '+%F %T')] $* ############"; }

step "CellOracle 链式流程开始"
echo "工作目录：$ROOT"
echo "base GRN：$BASE_GRN"
[ -f "$BASE_GRN" ] || { echo "[错误] base GRN 不存在：$BASE_GRN"; exit 1; }

# ===== 输入检查 =====
for f in 05_make_h5ad.R 01_tf_candidates_all.csv 03_build_core_TFs.py 04_celloracle_perturbation.py; do
  [ -f "$f" ] || { echo "[错误] 缺少文件：$f"; exit 1; }
done
if [ ! -f 05_scData_annotated.rds ] && [ ! -f 03_celloracle_input.h5ad ]; then
  echo "[错误] 既没有 05_scData_annotated.rds 也没有 03_celloracle_input.h5ad"; exit 1
fi

# ===== Python 环境 =====
if [ -n "${CELLORACLE_PYTHON:-}" ]; then
  PY="$CELLORACLE_PYTHON"
else
  source "$CONDA_SH"
  conda activate "$CONDA_ENV"
  PY="$(which python)"
fi
echo "python：$PY"
"$PY" -c "import celloracle, anndata, pandas; print('celloracle', celloracle.__version__, '| anndata', anndata.__version__)"
export BASE_GRN

# ===== Seurat → h5ad =====
if [ -f 03_celloracle_input.h5ad ]; then
  step "STEP1 h5ad 已存在，跳过"
else
  step "STEP1 Seurat -> h5ad"
  Rscript 05_make_h5ad.R
fi

# ===== 核心 TF =====
step "STEP2 核心 TF 选择"
"$PY" 03_build_core_TFs.py
cp output/03_key_transcription_factors.txt key_transcription_factors.txt
cat key_transcription_factors.txt

# ===== 建 GRN + KO/OE 扰动 =====
step "STEP3 建 GRN + KO/OE 虚拟扰动"
"$PY" 04_celloracle_perturbation.py

step "完成"
echo "结果目录：$ROOT/output"
ls -la "$ROOT/output" | head -30
