# RNA-seq 解析パイプライン (Jupyter Notebook版)

FASTQ 生データから DEG（発現変動遺伝子）のインタラクティブ可視化までを  
**4つの Jupyter Notebook** で実行する RNA-seq 解析パイプラインです。

---

## 動作環境

| 要件 | 推奨 |
|------|------|
| OS | Linux (Ubuntu 20.04+) |
| Miniconda | インストール済み |
| RAM | **32GB以上**（STAR マッピング時） |
| ストレージ | 100GB以上の空き |
| Jupyter | Notebook または JupyterLab |

---

## セットアップ手順

### 1. conda 環境の構築

```bash
# 環境作成（全ツール・ライブラリを一括インストール）
conda env create -f environment.yml

# アクティベート
conda activate rnaseq_env
```

> **高速化ヒント**: `conda` が遅い場合は `mamba` を使ってください  
> ```bash
> conda install -n base -c conda-forge mamba
> mamba env create -f environment.yml
> ```

### 2. Jupyter カーネルの登録

```bash
conda activate rnaseq_env

# Python カーネル
python -m ipykernel install --user --name rnaseq_env --display-name "RNA-seq (Python)"

# R カーネル
R -e "IRkernel::installspec(name='rnaseq_r', displayname='RNA-seq (R)')"
```

### 3. サンプルメタデータの準備

```bash
cp sample_metadata_template.csv sample_metadata.csv
```

エディタで `sample_metadata.csv` を開き、自分のサンプル情報を記入してください。

| カラム | 説明 | 例 |
|--------|------|-----|
| sample_id | ユニークなサンプル名 | Sample1 |
| fastq_r1 | R1 FASTQファイルパス | raw_data/Sample1_R1.fastq.gz |
| fastq_r2 | R2 FASTQファイルパス | raw_data/Sample1_R2.fastq.gz |
| condition | 実験条件 | Control / Treatment |
| replicate | レプリケート番号 | 1, 2, 3... |

### 4. リファレンスゲノムの準備 (hg38)

```bash
mkdir -p reference/hg38 && cd reference/hg38

# GENCODE v44 からダウンロード
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
gunzip *.gz

# STAR インデックス構築（30〜60分、32GB RAM 必要）
STAR --runMode genomeGenerate \
     --genomeDir star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v44.primary_assembly.annotation.gtf \
     --sjdbOverhang 149 \
     --runThreadN 8
```

### 5. FASTQファイルの配置

```bash
mkdir -p raw_data
# FASTQファイルを raw_data/ に配置（またはシンボリックリンク）
```

---

## ノートブック実行順序

**上から順に実行してください。** 各ノートブックは前のステップの出力に依存します。

| 順番 | ファイル | カーネル | 内容 | 所要時間目安 |
|:----:|---------|:-------:|------|:----------:|
| 1 | `01_QC_and_Trimming.ipynb` | Python | FastQC → Trim Galore | 10〜30分 |
| 2 | `02_Mapping_and_Counting.ipynb` | Python | STAR → featureCounts | 1〜3時間 |
| 3 | `03_DEG_Analysis.ipynb` | **R** | DESeq2 / edgeR | 5〜15分 |
| 4 | `04_Visualization.ipynb` | Python | Volcano Plot, ヒートマップ | 1〜5分 |

> **重要**: ノートブック03は **Rカーネル (RNA-seq (R))** に切り替えてから実行してください。
> ノートブック04は Python に戻してください。

---

## 各ノートブックの詳細

### 01_QC_and_Trimming.ipynb
- **FastQC**: 生データの品質チェック
- **Trim Galore**: アダプター除去・低品質リードトリミング (Paired-end)
- **MultiQC**: 全サンプルの QC レポートを統合
- 出力: `trimmed/`, `qc_reports/multiqc_report.html`

### 02_Mapping_and_Counting.ipynb
- **STAR**: hg38 へのマッピング (Paired-end)
- **featureCounts**: BAM → 遺伝子カウント行列
- マッピング統計のサマリー表示
- 出力: `mapped/`, `results/count_matrix.csv`

### 03_DEG_Analysis.ipynb (Rカーネル)
- **DESeq2 / edgeR** から選択可能（冒頭の `DEG_TOOL` 変数で切替）
- 全条件間のペアワイズ比較を自動実行
- PCA プロット / MDS プロット
- 出力: `results/deg_*_vs_*.csv`, `results/all_deg_results.csv`

### 04_Visualization.ipynb
- **Volcano Plot**: ドロップダウンで比較条件を切替可能（Plotly）
- **全条件比較ヒートマップ**: 赤=増加 / グレー=変化なし / 青=減少
- DEG 数バープロット
- 出力: `results/volcano_plot_interactive.html`, `results/deg_heatmap_interactive.html`

---

## カスタマイズ可能なパラメータ

各ノートブック冒頭の **設定セル** で変更できます。

| パラメータ | デフォルト | 説明 |
|-----------|:---------:|------|
| `THREADS` | 8 | 並列処理スレッド数 |
| `LFC_THRESHOLD` | 1.0 | log2FoldChange の閾値 |
| `PADJ_THRESHOLD` | 0.05 | 調整済 P 値の閾値 |
| `GENOME_DIR` | reference/hg38/star_index | STAR インデックスのパス |
| `GTF_FILE` | (GENCODE v44 パス) | GTF アノテーションファイル |
| `DEG_TOOL` | "deseq2" | DEG ツール選択 ("deseq2" or "edger") |

---

## ディレクトリ構造

```
rnaseq_project/
├── README.md                        ← このファイル
├── environment.yml                  ← conda 環境定義
├── requirements.txt                 ← pip 追加パッケージ
├── sample_metadata_template.csv     ← テンプレート（コピーして編集）
├── sample_metadata.csv              ← ユーザーが編集
│
├── 01_QC_and_Trimming.ipynb         ← Step 1-2: QC, トリミング
├── 02_Mapping_and_Counting.ipynb    ← Step 3-5: STAR, featureCounts
├── 03_DEG_Analysis.ipynb            ← Step 6: DEG 解析 (R カーネル)
├── 04_Visualization.ipynb           ← Step 7: 可視化
│
├── raw_data/                        ← FASTQ ファイルを配置
│   ├── Sample1_R1.fastq.gz
│   ├── Sample1_R2.fastq.gz
│   └── ...
│
├── reference/hg38/                  ← リファレンスゲノム
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── gencode.v44.primary_assembly.annotation.gtf
│   └── star_index/                  ← STAR インデックス
│
├── trimmed/                         ← (自動生成) トリミング後 FASTQ
├── qc_reports/                      ← (自動生成) QC レポート
├── mapped/                          ← (自動生成) BAM ファイル
└── results/                         ← (自動生成) 解析結果
    ├── count_matrix.csv
    ├── featureCounts_output.txt
    ├── deg_CondA_vs_CondB.csv       ← 各比較の DEG 結果
    ├── all_deg_results.csv          ← 全比較統合
    ├── pca_plot.pdf
    ├── volcano_plot_interactive.html
    └── deg_heatmap_interactive.html
```

---

## 可視化の操作方法

### Volcano Plot (volcano_plot_interactive.html)
- **ドロップダウンメニュー** で比較条件を切替
- マウスホバーで遺伝子名・log2FC・P 値を表示
- マウスドラッグでズーム、ダブルクリックでリセット
- 色: **赤** = Up (増加), **青** = Down (減少), **グレー** = NS (変化なし)

### DEG ヒートマップ (deg_heatmap_interactive.html)
- 全比較条件を一覧で表示
- **赤** = 増加, **グレー** = 変化なし, **青** = 減少
- 遺伝子方向に Ward 法でクラスタリング済み
- マウスホバーで遺伝子名・比較条件・log2FC・DEG ステータスを表示

---

## トラブルシューティング

### conda が遅い / 依存解決に失敗する
```bash
conda install -n base -c conda-forge mamba
mamba env create -f environment.yml
```

### STAR でメモリエラー (cannot allocate memory)
hg38 では STAR に約 32GB RAM が必要です。メモリ不足の場合は HISAT2 への切替を検討してください。

### featureCounts で "Unassigned" が多すぎる
GTF ファイルとゲノム FASTA が **同じ GENCODE リリース** (v44) であることを確認してください。

### R カーネルが Jupyter に表示されない
```bash
conda activate rnaseq_env
R -e "install.packages('IRkernel', repos='https://cloud.r-project.org')"
R -e "IRkernel::installspec(name='rnaseq_r', displayname='RNA-seq (R)')"
```
Jupyter を再起動してください。

### Plotly が JupyterLab で表示されない
```bash
pip install jupyterlab-plotly
```

### MultiQC でエラーが出る
```bash
pip install --upgrade multiqc
```

---

## 使用ツール / ライセンス

| ツール | バージョン | 用途 |
|-------|:---------:|------|
| FastQC | 0.12.1 | 品質管理 |
| MultiQC | 1.21 | QC レポート統合 |
| Trim Galore | 0.6.10 | アダプター除去 |
| STAR | 2.7.11b | マッピング |
| Subread (featureCounts) | 2.0.6 | リードカウント |
| DESeq2 | Bioconductor | DEG 解析 |
| edgeR | Bioconductor | DEG 解析 |
| Plotly | 5.18+ | インタラクティブ可視化 |
