#!/usr/bin/env bash
# 方案第七节(5) Figure 7H：pySCENIC 三步（GRN → motif 富集 → AUCell）
# 用法：nohup bash 10_scenic_chain.sh > scenic_chain.log 2>&1 &
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p scenic
cd scenic

RANKINGS="${RANKINGS:-$HOME/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather}"
MOTIFS="${MOTIFS:-$HOME/motifs-v9-nr.hgnc-m0.001-o0.0.tbl}"
PYSCENIC="${PYSCENIC:-$HOME/anaconda3/envs/pyscenic/bin/pyscenic}"
N_WORKERS="${N_WORKERS:-8}"

step() { echo; echo "############ [$(date '+%F %T')] $* ############"; }

step "pySCENIC 链式流程开始"
echo "工作目录：$(pwd)"
for f in "$RANKINGS" "$MOTIFS" ../09_scenic_counts.csv ../09_scenic_cells.csv; do
  [ -f "$f" ] || { echo "[错误] 缺少文件：$f"; exit 1; }
done
"$PYSCENIC" --help > /dev/null

# TF 列表：从 motif 注释表的 gene_name 列导出（服务器上无现成 allTFs_hg38.txt）
step "STEP0 导出 TF 列表"
awk -F'\t' 'NR>1 && $6!="" {print $6}' "$MOTIFS" | sort -u > allTFs_hgnc.txt
echo "TF 数：$(wc -l < allTFs_hgnc.txt)"

# 表达矩阵为「基因 × 细胞」，pySCENIC 的 CSV 约定是「细胞 × 基因」，故全部加 -t
step "STEP1 GRN 推断（GRNBoost2）"
"$PYSCENIC" grn ../09_scenic_counts.csv allTFs_hgnc.txt -t \
  -o 01_grn.csv --num_workers "$N_WORKERS" --seed 123
echo "GRN 边数：$(($(wc -l < 01_grn.csv) - 1))"

step "STEP2 motif 富集（ctx）"
"$PYSCENIC" ctx 01_grn.csv "$RANKINGS" -t \
  --annotations_fname "$MOTIFS" \
  --expression_mtx_fname ../09_scenic_counts.csv \
  -o 02_regulons.csv --num_workers "$N_WORKERS"
echo "调控子数：$(($(wc -l < 02_regulons.csv) - 1))"

step "STEP3 AUCell"
"$PYSCENIC" aucell ../09_scenic_counts.csv 02_regulons.csv -t \
  -o 03_aucell.csv --num_workers "$N_WORKERS" --seed 123

step "完成"
ls -la
