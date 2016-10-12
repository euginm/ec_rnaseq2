# Table of contents

- [Dependencies](#dependencies)
- [Quantification with `kallisto`](#quantification-with-kallisto)
 * [Download and index transcripts](#download-and-index-transcripts)
 * [Make `tx2gene` file](#make-tx2gene-file)
 * [Run `kallisto`](#run-kallisto)
- [Fit model in `R`](#fit-model-in-r)
- [Basic visualisation](#basic-visualisation)
 * [Dispersions](#dispersions)
 * [PCA](#pca)
 * [Overlap in expression](#overlap-in-expression)

# Dependencies

* `kallisto` 0.43.0
* `R`: tximport, DESeq2, dplyr, readr, tidyr, tibble, ggplot2, gplots, GenomicFeatures

# Directory structure

```
├───ec_rnaseq2/
│       README.md
│       plots/

...

├───$KALL/
│
└───$ADAT/
        C1_S10_L001_R1_001.fastq.gz
        C1_S10_L001_R2_001.fastq.gz
        C2_S11_L001_R1_001.fastq.gz
        C2_S11_L001_R2_001.fastq.gz
        C3_S12_L001_R1_001.fastq.gz
        C3_S12_L001_R2_001.fastq.gz
        Cip1_S1_L001_R1_001.fastq.gz
        Cip1_S1_L001_R2_001.fastq.gz
        Cip2_S2_L001_R1_001.fastq.gz
        Cip2_S2_L001_R2_001.fastq.gz
        Cip3_S3_L001_R1_001.fastq.gz
        Cip3_S3_L001_R2_001.fastq.gz
        Mel1_S7_L001_R1_001.fastq.gz
        Mel1_S7_L001_R2_001.fastq.gz
        Mel2_S8_L001_R1_001.fastq.gz
        Mel2_S8_L001_R2_001.fastq.gz
        Mel3_S9_L001_R1_001.fastq.gz
        Mel3_S9_L001_R2_001.fastq.gz
        Pex1_S4_L001_R1_001.fastq.gz
        Pex1_S4_L001_R2_001.fastq.gz
        Pex2_S5_L001_R1_001.fastq.gz
        Pex2_S5_L001_R2_001.fastq.gz
        Pex3_S6_L001_R1_001.fastq.gz
        Pex3_S6_L001_R2_001.fastq.gz

```


# Quantification with `kallisto`

```sh
source vars
```

```
cat vars
export ADAT="/path/to/fastqs/"
export KALL="/path/to/kallisto/output/"
export R1="L001_R1_001.fastq.gz"
export R2="L001_R2_001.fastq.gz"
```

```sh
# get sample prefixes.
cd $ADAT &&
for i in *fastq.gz; do
    echo $i
done | cut -f 1-2 -d "_" | sort -u > $KALL/samples
```

## Download and index transcripts

```sh
cd $KALL &&
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-28/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.cdna.all.fa.gz
gunzip Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.cdna.all.fa.gz
kallisto index -i cdna Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.cdna.all.fa
```

### Make `tx2gene` file

```sh
grep ">" cdna.fasta | sed 's/>//g' | cut -f1,4 -d" " | sed 's/ gene:/\t/g' | sed '1 i\TXNAME\tGENEID' > tx2gene.tsv
```

### Run `kallisto`

```sh
for i in `cat samples`; do
    kallisto quant -i cdna -t 4 -o $i --rf-stranded $ADAT/${i}_$R1 $ADAT/${i}_$R2
done
```

## Fit model in `R`

```r
library(readr)
library(tximport)
library(tidyr)
library(dplyr)
library(DESeq2)
library(cowplot) # ggplot2 themes
library(tibble) # as_tibble()
library(gplots) # euler diagram
library(GenomicFeatures)
dir <- "./"
run <- readLines("samples")
files <- file.path(dir, run, "abundance.tsv")
names(files) <- run
all(file.exists(files))
tx2gene <- read_tsv("tx2gene.tsv")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)
sampleTable <- as.data.frame(colnames(txi$counts))
colnames(sampleTable) <- "library"
sampleTable <- separate(sampleTable, library, into = "treatment", sep = "_", remove = FALSE, extra = "drop")
sampleTable$treatment <- tolower(sampleTable$treatment)
sampleTable$treatment <- gsub("c[1-3]", "con", sampleTable$treatment)
sampleTable$treatment <- gsub("[1-3]", "", sampleTable$treatment) %>% as.factor %>% relevel(ref = "con")
rownames(sampleTable) <- sampleTable$library
sampleTable <- dplyr::select(sampleTable, -library)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ treatment)
dds <- DESeq(dds, parallel = TRUE, fitType = "local")
ddp <- DESeq(dds, parallel = TRUE, fitType = "parametric")
```

## Basic visualisation

### Dispersions

```r
par(mfrow = c(1, 2))
plotDispEsts(dds)
title("local")
plotDispEsts(ddp)
title("parametric")
```

![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/disp.png)

We go with `fitType = "local"`.

### PCA

```r
rld <- rlog(dds, fitType = "local")
rldf <- assay(rld) %>% as.data.frame
pca <- prcomp(t(rldf))
plot(pca)
```

![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/pcs.png)

Let's keep the first 3 PCs.

```r
pcax <- as.data.frame(pca$x)
pcax$treatment <- sampleTable$treatment
# this is super hacky
prop <- summary(pca) %>% as.list %>% .$importance %>% t %>% as_tibble %>%
  .[1:3,] %>% dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
pc12 <- ggplot(pcax, aes(x=PC1, y=PC2, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc13 <- ggplot(pcax, aes(x=PC1, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC1: ", round(prop[1]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  scale_color_discrete(guide=FALSE)
pc23 <- ggplot(pcax, aes(x=PC2, y=PC3, colour=treatment)) + geom_point(size=5) +
  xlab(paste("PC2: ", round(prop[2]), "% variance", sep = "")) +
  ylab(paste("PC3: ", round(prop[3]), "% variance", sep = "")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0))
plot_grid(pc12, pc13, pc23, nrow = 1)
```

![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/pca.png)

```r
rot <- as.data.frame(pca$rotation) %>% rownames_to_column
par(mfrow = c(1, 3))
plot(rot$PC1, xlab = "", xaxt = "n", ylab = "PC1 loading")
abline(h = 0)
title("A. PC1")
plot(rot$PC2, xlab = "", xaxt = "n", ylab = "PC2 loading")
abline(h = 0)
title("B. PC2")
plot(rot$PC3, xlab = "", xaxt = "n", ylab = "PC3 loading")
title("C. PC3")
abline(h = 0)
```

![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/rot.png)

### Overlap in expression

```r
cip_res <- results(dds, alpha = 0.05, addMLE = TRUE, contrast = c("treatment", "cip", "con")) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1) %>% as_tibble %>% rownames_to_column
pex_res <- results(dds, alpha = 0.05, addMLE = TRUE, contrast = c("treatment", "pex", "con")) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1) %>% as_tibble %>% rownames_to_column
mel_res <- results(dds, alpha = 0.05, addMLE = TRUE, contrast = c("treatment", "mel", "con")) %>% subset(padj < 0.05 & abs(log2FoldChange) > 1) %>% as_tibble %>% rownames_to_column
```

We compare the overlap in gene expression between the three treatments. Plot A shows all differentially-expressed genes regardless of sign. Plots B and C show up- and down-regulated genes respectively.

```r
eul1 <- list(
  cip = cip_res$rowname,
  pex = pex_res$rowname,
  mel = mel_res$rowname
)
eul2 <- list(
  cip = filter(cip_res, sign(log2FoldChange) == 1) %>% .$rowname,
  pex = filter(pex_res, sign(log2FoldChange) == 1) %>% .$rowname,
  mel = filter(mel_res, sign(log2FoldChange) == 1) %>% .$rowname
)
eul3 <- list(
  cip = filter(cip_res, sign(log2FoldChange) != 1) %>% .$rowname,
  pex = filter(pex_res, sign(log2FoldChange) != 1) %>% .$rowname,
  mel = filter(mel_res, sign(log2FoldChange) != 1) %>% .$rowname
)
par(mfrow = c(2, 2), oma = c(0, 0, 0, 0), mar = c(0, 0, 4, 0), cex = 1.05)
venn(eul1)
title("A. All DE genes")
venn(eul2)
title("B. Up-regulated genes")
venn(eul3)
title("C. Down-regulated genes")
```

![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/eul.png)

### Annotated results

```r
aDB <- makeTxDbFromGFF("Escherichia_coli_str_k_12_substr_mg1655.GCA_000005845.2.28.gtf")
gDB <- genes(aDB, columns = columns(aDB)) %>% as.data.frame %>% rownames_to_column
cip_res_db <- inner_join(dplyr::filter(gDB, rowname %in% cip_res$rowname), cip_res)
pex_res_db <- inner_join(dplyr::filter(gDB, rowname %in% pex_res$rowname), pex_res)
mel_res_db <- inner_join(dplyr::filter(gDB, rowname %in% mel_res$rowname), mel_res)
```

## To do

- [ ] Annotated results
- [ ] use rtracklayer::import() to get additional annotation fields (gene_name, protein_name etc.)
- [ ] GOs
