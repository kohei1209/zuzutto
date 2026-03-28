#!/usr/bin/env Rscript
# =============================================================================
# 03. DEG解析 (発現変動遺伝子解析) — RStudio用スクリプト
# =============================================================================
#
# Jupyter Notebook (03_DEG_Analysis.ipynb) の代替スクリプトです。
# Rカーネルが不安定な場合にRStudioまたはコマンドラインから実行できます。
#
# 使い方:
#   RStudio: このファイルを開いて Source ボタンで実行
#   コマンドライン: Rscript 03_DEG_Analysis.R
#
# 入力:
#   - results/count_matrix.csv  (02_Mapping_and_Counting.ipynb の出力)
#   - sample_metadata.csv
#
# 出力:
#   - results/deg_*_vs_*.csv    (各ペアワイズ比較結果)
#   - results/all_deg_results.csv (全比較統合)
#   - results/pca_plot.pdf
#
# 注意: 作業ディレクトリを RNA-seq/ フォルダに設定してから実行してください
#   setwd("/path/to/RNA-seq")
# =============================================================================

# ===========================
# 設定 (必要に応じて変更)
# ===========================
LFC_THRESHOLD <- 1.0
PADJ_THRESHOLD <- 0.05
DEG_TOOL <- "deseq2"  # "deseq2" or "edger"

# 前フィルタリングパラメータ
MIN_COUNT_THRESHOLD <- 10   # 最小カウント数

# ===========================
# ライブラリ読み込み
# ===========================
required_packages <- c("DESeq2", "edgeR", "ggplot2", "pheatmap",
                       "dplyr", "readr", "tidyr", "RColorBrewer")

missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste0(
    "以下のRパッケージがインストールされていません:\n  ",
    paste(missing, collapse = ", "),
    "\n\nconda環境が有効か確認してください:\n  conda activate rnaseq_env\n",
    "または BiocManager でインストール:\n",
    "  BiocManager::install(c('DESeq2', 'edgeR'))"
  ))
}

library(DESeq2)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(tidyr)
library(RColorBrewer)

cat("全パッケージの読み込み完了\n")
cat("DEG解析ツール:", DEG_TOOL, "\n\n")

# ===========================
# 出力ディレクトリ作成
# ===========================
dir.create("results", showWarnings = FALSE)

# ===========================
# データ読み込み
# ===========================
if (!file.exists("results/count_matrix.csv")) {
  stop("results/count_matrix.csv が見つかりません。\n先に 02_Mapping_and_Counting.ipynb を実行してください。")
}
if (!file.exists("sample_metadata.csv")) {
  stop("sample_metadata.csv が見つかりません。\nREADME.md のセットアップ手順を確認してください。")
}

count_data <- read.csv("results/count_matrix.csv", row.names = 1, check.names = FALSE)

sample_info <- read.csv("sample_metadata.csv")
rownames(sample_info) <- sample_info$sample_id

# --- 条件名のバリデーション ---
conditions <- unique(sample_info$condition)
has_space <- grepl("\\s", conditions)
if (any(has_space)) {
  bad <- conditions[has_space]
  stop(paste0(
    "条件名にスペースが含まれています: ", paste(bad, collapse = ", "),
    "\nスペースをアンダースコア(_)に置換してください。",
    "\n例: 'compound A 0.5 nM' -> 'CompA_0.5nM'"
  ))
}
cat("条件名チェック: OK\n")

# 順序を合わせる
count_data <- count_data[, sample_info$sample_id]

cat("カウント行列:", nrow(count_data), "genes x", ncol(count_data), "samples\n")
cat("条件:", paste(conditions, collapse = ", "), "\n")
cat("各条件のサンプル数:\n")
print(table(sample_info$condition))

n_comp <- choose(length(conditions), 2)
cat("\n全ペアワイズ比較数:", n_comp, "組\n\n")

# ===========================
# 前フィルタリング
# ===========================
min_group_size <- min(table(sample_info$condition))
keep <- rowSums(count_data >= MIN_COUNT_THRESHOLD) >= min_group_size
count_data <- count_data[keep, ]
cat("フィルタ後:", nrow(count_data), "genes\n")
cat("(閾値: カウント >=", MIN_COUNT_THRESHOLD, "のサンプルが", min_group_size, "以上)\n\n")

# ===========================
# DEG解析関数
# ===========================

# --- DESeq2 ---
run_deseq2 <- function(counts, coldata, pair) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts),
    colData = coldata,
    design = ~ condition
  )
  dds$condition <- relevel(dds$condition, ref = pair[2])
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", pair[1], pair[2]), alpha = 0.05)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$comparison <- paste0(pair[1], "_vs_", pair[2])

  return(list(dds = dds, results = res_df))
}

# --- edgeR ---
run_edger <- function(counts, coldata, pair) {
  idx <- coldata$condition %in% pair
  sub_counts <- counts[, idx]
  sub_coldata <- coldata[idx, ]

  group <- factor(sub_coldata$condition, levels = c(pair[2], pair[1]))
  y <- DGEList(counts = sub_counts, group = group)
  y <- calcNormFactors(y)

  design <- model.matrix(~ group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)

  res <- topTags(qlf, n = Inf, sort.by = "none")$table
  res$gene_id <- rownames(res)
  res$comparison <- paste0(pair[1], "_vs_", pair[2])

  colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
  colnames(res)[colnames(res) == "PValue"] <- "pvalue"
  colnames(res)[colnames(res) == "FDR"] <- "padj"

  return(list(fit = fit, results = res))
}

# ===========================
# 全ペアワイズ比較の実行
# ===========================
comparisons <- combn(conditions, 2, simplify = FALSE)
all_results <- list()

for (pair in comparisons) {
  comp_name <- paste0(pair[1], "_vs_", pair[2])
  cat("比較実行中:", comp_name, "\n")

  if (DEG_TOOL == "deseq2") {
    result <- run_deseq2(count_data, sample_info, pair)
  } else {
    result <- run_edger(count_data, sample_info, pair)
  }

  all_results[[comp_name]] <- result$results

  # 個別CSV出力
  write.csv(result$results,
            file = paste0("results/deg_", comp_name, ".csv"),
            row.names = FALSE)

  # DEG数を表示
  n_up <- sum(result$results$padj < PADJ_THRESHOLD &
              result$results$log2FoldChange > LFC_THRESHOLD, na.rm = TRUE)
  n_down <- sum(result$results$padj < PADJ_THRESHOLD &
                result$results$log2FoldChange < -LFC_THRESHOLD, na.rm = TRUE)
  cat("  Up:", n_up, " Down:", n_down, "\n")
}

# 全結果を統合
all_deg <- do.call(rbind, all_results)
write.csv(all_deg, "results/all_deg_results.csv", row.names = FALSE)
cat("\n全比較結果を results/all_deg_results.csv に出力しました\n\n")

# ===========================
# QC: PCAプロット
# ===========================
cat("PCAプロットを生成中...\n")

dds_all <- DESeqDataSetFromMatrix(
  countData = round(count_data),
  colData = sample_info,
  design = ~ condition
)
vsd <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggtitle("PCA - 全サンプル")

ggsave("results/pca_plot.pdf", p, width = 8, height = 6)
cat("出力: results/pca_plot.pdf\n")

# ===========================
# 完了
# ===========================
cat("\n")
cat("========================================\n")
cat("  DEG解析が完了しました\n")
cat("========================================\n")
cat("\n")
cat("出力ファイル:\n")
for (pair in comparisons) {
  comp_name <- paste0(pair[1], "_vs_", pair[2])
  cat("  - results/deg_", comp_name, ".csv\n", sep = "")
}
cat("  - results/all_deg_results.csv (全比較統合)\n")
cat("  - results/pca_plot.pdf\n")
cat("\n次のステップ:\n")
cat("  04_Visualization.ipynb (Pythonカーネル) に進んでください\n")
