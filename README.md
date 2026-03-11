# Proteome pseudo-GSEA

このディレクトリは、proteome データに対して RNA-seq の [04_deseq2.R](04_deseq2.R) と [05_fgsea.R](05_fgsea.R) の流れに準拠した擬似 GSEA を行うための最小構成です。

## 追加したファイル

- [04_proteome_limma.R](04_proteome_limma.R)
  - proteome_data.csv を読み込みます。
  - サンプル名から condition と batch を抽出します。
  - log2 変換、中央値センタリング正規化、batch を含む limma による差次解析を行います。
  - GSEA 用に moderated t statistic を rank_metric として出力します。

- [05_proteome_fgsea.R](05_proteome_fgsea.R)
   - MsigDB の H, C2, C5, C7 と [Gene_set/Human_old_HSC_set2.gmt](Gene_set/Human_old_HSC_set2.gmt) を読み込みます。
  - [04_proteome_limma.R](04_proteome_limma.R) が出力した差次解析結果を使って fgsea を実行します。

## 実行方法

必要な R パッケージ:

- data.table
- limma
- fgsea
- msigdbr

現在のスクリプトは、不足している package があれば自動で install を試みます。

RStudio での実行例:

```r
setwd("C:/Users/hikob/Dropbox/Research/Proteome_GSEA")
source("04_proteome_limma.R")
source("05_proteome_fgsea.R")
```

出力先を変えたい場合:

```r
source("04_proteome_limma.R", local = list(args = c("proteome_data.csv", "Proteome_DE")))
source("05_proteome_fgsea.R", local = list(args = c("Proteome_DE", "Gene_set/Human_old_HSC_set2.gmt", "Proteome_fGSEA", "Mus musculus")))
```

Rscript で実行する場合:

```r
Rscript 04_proteome_limma.R
Rscript 05_proteome_fgsea.R
```

引数を明示する場合:

```r
Rscript 04_proteome_limma.R proteome_data.csv Proteome_DE
Rscript 05_proteome_fgsea.R Proteome_DE Gene_set/Human_old_HSC_set2.gmt Proteome_fGSEA "Mus musculus"
```

## 出力

- Proteome_DE
  - sample_metadata.tsv
  - normalized_proteome_matrix.tsv
  - sample_boxplot_before_normalization.pdf
  - sample_boxplot_after_normalization.pdf
  - sample_pca_after_normalization.pdf
  - 各比較ごとの *_DE_Results.tsv

- Proteome_fGSEA
  - 各比較ごとの fgsea_Results_*.tsv
  - 各比較ごとの barplot PDF
   - fgsea_summary.txt

結果の pathway 名は `H::`, `C2::`, `C5::`, `C7::`, `CUSTOM::` の接頭辞付きで出力されます。

## proteome に GSEA を使う際の注意

RNA-seq と異なり、proteome の定量値は count ではなく連続値で、欠損や測定バッチの影響も強く受けます。そのため DESeq2 をそのまま流用せず、以下の点でバイアスを抑えています。

1. 正規化
   - サンプルごとの log2 強度を中央値センタリングし、ロード量や全体強度のずれを抑えます。

2. batch 補正
   - このデータでは列名に EXP と REQ が含まれているため、これを batch として線形モデルに入れます。

3. 欠損値の扱い
   - 差次解析そのものには一括 imputation を入れていません。
   - 欠損補完は fold change を過大にしやすいため、解析本体では観測値ベースを優先しています。
   - PCA 図だけは可視化のために行単位中央値で簡易補完しています。

4. GSEA の順位指標
   - raw intensity や raw logFC ではなく、limma の moderated t statistic を rank に使います。
   - これにより高発現タンパク質だけが過度に優位になる偏りを減らします。

5. species 差
   - MsigDB の H, C2, C5, C7 は `msigdbr` から `Mus musculus` 指定で取得するため、proteome 側との整合性は custom GMT より高いです。
   - 一方で [Gene_set/Human_old_HSC_set2.gmt](Gene_set/Human_old_HSC_set2.gmt) は human 名称のままなので、現在のスクリプトでは大文字化による簡便マッチングを使っています。
   - 厳密な解釈が必要なら、custom GMT 側には正式な ortholog 変換を別途入れてください。

## GitHub について

このディレクトリは git 初期化済みで、origin は takubo-lab/proteome_GSEA に設定しています。通常は以下で更新できます。

```bash
git add .
git commit -m "Add proteome pseudo-GSEA workflow"
git push -u origin main
```