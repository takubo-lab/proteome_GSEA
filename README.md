# Proteome pseudo-GSEA workflow

Proteome の定量テーブルから差次解析と fgsea を実行するための公開向けワークフローです。count ベースの RNA-seq 解析ではなく、proteome の連続値データを前提に、log2 変換、中央値センタリング、limma による線形モデリング、fgsea による pathway enrichment を行います。

このワークフローは次の点を前提に設計しています。

- サンプル数は固定ではありません。
- グループ数は 2 群に限定されません。
- 群間比較は、入力 metadata に含まれるすべての condition の総当たりで自動生成されます。
- batch 情報がある場合は線形モデルに組み込み、なければ condition のみで解析します。

## Files

- [04_proteome_limma.R](04_proteome_limma.R)
  - proteome 定量値の前処理と差次解析を行います。

- [05_proteome_fgsea.R](05_proteome_fgsea.R)
  - 差次解析結果をもとに fgsea を行います。
  - MsigDB の H, C2, C5, C7 に加えて [Gene_set/Human_old_HSC_set2.gmt](Gene_set/Human_old_HSC_set2.gmt) を解析対象に含めます。

## Input format

### 1. Proteome quantification table

入力は CSV 形式を想定しています。最低限、次のどちらかの条件を満たしてください。

- `Genes` 列がある
- `Protein.Group` 列があり、gene symbol の代替として使える

推奨カラム構成は次の通りです。

| Column | Required | Description |
| --- | --- | --- |
| Protein.Group | optional | protein ID または protein group ID |
| Protein.Names | optional | protein name |
| Genes | recommended | gene symbol。重複 protein を gene 単位にまとめる際に使用 |
| First.Protein.Description | optional | protein description |
| sample columns | required | 各サンプルの定量値 |

サンプル列は数値である必要があります。0 以下の値は欠損として扱われます。

例:

```csv
Protein.Group,Protein.Names,Genes,First.Protein.Description,Ctrl_1,Ctrl_2,Treat_1,Treat_2
P00001,Protein A,GeneA,Example protein A,1250000,1180000,1580000,1490000
P00002,Protein B,GeneB,Example protein B,240000,255000,198000,205000
P00003,Protein C,GeneC,Example protein C,53000,61000,88000,91000
```

### 2. Sample metadata table

公開向けの利用では、sample metadata を別ファイルで渡すことを推奨します。これにより、任意のサンプル名、任意のグループ数、任意の batch 設計に対応できます。

推奨形式は TSV または CSV で、必須列は次の 2 つです。

| Column | Required | Description |
| --- | --- | --- |
| sample | required | proteome table のサンプル列名と完全一致する名前 |
| condition | required | 群名 |
| batch | optional | バッチ名。省略時は単一 batch として扱う |
| replicate | optional | 任意の replicate ID |

例:

```tsv
sample	condition	batch	replicate
Ctrl_1	Ctrl	Batch1	1
Ctrl_2	Ctrl	Batch1	2
Treat_1	Treat	Batch2	1
Treat_2	Treat	Batch2	2
```

metadata を指定しない場合は、サンプル列名が `condition-batchRep` 形式であると仮定して自動推定します。これは後方互換のための fallback であり、公開用途では metadata 指定を推奨します。

CSV 形式でも問題ありません。例えば次のような内容です。

```csv
sample,condition,batch,replicate
Ctrl_1,Ctrl,Batch1,1
Ctrl_2,Ctrl,Batch1,2
Treat_1,Treat,Batch2,1
Treat_2,Treat,Batch2,2
```

## Analysis strategy

### Differential analysis

[04_proteome_limma.R](04_proteome_limma.R) は次の流れで差次解析を行います。

1. 定量値を読み込む
2. 正の値の最小値をもとに pseudocount を設定して log2 変換する
3. サンプルごとの中央値センタリングを行う
4. condition ごとに最低限の non-missing 値を満たす protein のみ残す
5. 同一 gene symbol に対応する重複 protein は、欠損が少なく平均強度が高いものを優先して 1 行に集約する
6. limma で全 condition 間の pairwise comparison を実行する
7. fgsea 用 ranking metric として moderated t statistic を出力する

### Gene set enrichment analysis

[05_proteome_fgsea.R](05_proteome_fgsea.R) は次の gene set を読み込みます。

- MsigDB H
- MsigDB C2
- MsigDB C5
- MsigDB C7
- custom GMT

結果の pathway 名には collection の接頭辞が付きます。

- `H::`
- `C2::`
- `C5::`
- `C7::`
- `CUSTOM::`

## Why this is proteome-specific

proteome の GSEA では、RNA-seq と異なるバイアスに注意が必要です。このワークフローでは以下を明示的に考慮しています。

1. count モデルは使わない
   proteome は count ではなく連続値なので、DESeq2 ではなく limma を用います。

2. 全体強度のずれを抑える
   サンプル間のロード量や測定スケール差を減らすため、log2 変換後に中央値センタリングを行います。

3. 欠損の過剰補完を避ける
   差次解析本体では全体 imputation を行いません。PCA 可視化のみ行単位中央値で簡易補完します。

4. batch を任意で取り込む
   batch 列が与えられ、かつ設計行列が full rank になる場合のみ batch を線形モデルに入れます。

5. ranking metric を moderated t にする
   raw intensity や raw logFC のみで ranking すると高発現 protein に偏りやすいため、moderated t statistic を使います。

6. species 差を分けて扱う
   MsigDB は `msigdbr` から species 指定で取得します。一方で custom GMT はファイル内容をそのまま使うため、human 系の GMT を mouse proteome に使う場合は別途 ortholog 変換を検討してください。

7. gene symbol 表記を一貫化する
   proteome 側の gene symbol、MsigDB、custom GMT はすべて `trimws + toupper` で正規化してから突き合わせます。human と mouse で大文字小文字の流儀が異なっていても一致しやすくするためです。ただし、これは表記ゆれ吸収であって ortholog 変換ではありません。

## Requirements

R で以下の package が必要です。

- data.table
- limma
- fgsea
- msigdbr

スクリプト内で不足 package の install を自動的に試みます。

## Usage

### RStudio

metadata を使わない最小実行例:

```r
source("04_proteome_limma.R")
source("05_proteome_fgsea.R")
```

metadata を使う推奨実行例:

```r
source("04_proteome_limma.R", local = list(args = c(
  "input_proteome.csv",
  "Proteome_DE",
  "sample_metadata.tsv"
)))

source("05_proteome_fgsea.R", local = list(args = c(
  "Proteome_DE",
  "Gene_set/Human_old_HSC_set2.gmt",
  "Proteome_fGSEA",
  "Mus musculus"
)))
```

### Command line

metadata を使う推奨実行例:

```bash
Rscript 04_proteome_limma.R input_proteome.csv Proteome_DE sample_metadata.tsv
Rscript 05_proteome_fgsea.R Proteome_DE Gene_set/Human_old_HSC_set2.gmt Proteome_fGSEA "Mus musculus"
```

## Output

### Proteome_DE

- `sample_metadata.tsv`
- `normalized_proteome_matrix.tsv`
- `sample_boxplot_before_normalization.pdf`
- `sample_boxplot_after_normalization.pdf`
- `sample_pca_after_normalization.pdf`
- `analysis_summary.txt`
- `*_DE_Results.tsv`

### Proteome_fGSEA

- comparison ごとの `fgsea_Results_*.tsv`
- comparison ごとの `fgsea_barplot_*.pdf`
- `fgsea_summary.txt`

## Notes for interpretation

- condition 数が増えると pairwise comparison の数は自動的に増えます。
- replicate 数が少ない condition では、欠損フィルタ後に利用できる protein 数が減ることがあります。
- custom GMT が human ベースの場合でも、gene symbol は大文字化して照合します。ただし species 間の 1 対 1 対応を保証するものではありません。
- `overlap_n` が小さい pathway は過剰解釈しない方が安全です。