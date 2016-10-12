# Table of contents

- [Dependencies](#dependencies)
- [Directory structure](#directory-structure)
- [Different quantification methods overview](#different-quantification-methods-overview)
- [Download transcripts and annotations](#download-transcripts-and-annotations)
- [QC](#qc)
- [Check rRNA contamination](#check-rrna-contamination)
 * [Index ncrna](##index-ncrna)
 * [Run kallisto with ncrna](##run-kallisto-with-ncrna)
 * [Summarize](##summarize)
- [Index and quantify reads with `kallisto` and `salmon`](#index-and-quantify-reads-with-kallisto-and-salmon)
- [Analysis](#analysis)
 * [Fit model in `R`](##fit-model-in-r)
 * [Basic visualisation](##basic-visualisation)
 * [Overlap in gene expression between groups](##overlap-in-gene-expression-between-groups)
 * [Compare `kallisto` and `salmon` results](##compare-kallisto-and-salmon-results)
 * [Export DE genes tables](##export-de-genes-tables)
- [GO enrichment](#go-enrichment)

# Dependencies

* `kallisto` 0.43.0
* `salmon` 0.7.2
* `python` annotation.py
* `R`: Rqc, tximport, readr, DESeq2, tidyr, dplyr, pheatmap, RColorBrewer, magrittr, tibble, gplots, ggplot2, GenomicFeatures, cowplot, plotly, goseq

# Directory structure

```
├───ec_rnaseq2/
│       README.md
│       ec2_qc.html
│       plots/
│       python_scripts/

...

├───ref/
│
├───$KALL/
│       ncrna/
│
├───$SALM/
│
└───$READS/
        C1_S10_L001_R1_001.fastq
        C1_S10_L001_R2_001.fastq
        C2_S11_L001_R1_001.fastq
        C2_S11_L001_R2_001.fastq
        C3_S12_L001_R1_001.fastq
        C3_S12_L001_R2_001.fastq
        Cip1_S1_L001_R1_001.fastq
        Cip1_S1_L001_R2_001.fastq
        Cip2_S2_L001_R1_001.fastq
        Cip2_S2_L001_R2_001.fastq
        Cip3_S3_L001_R1_001.fastq
        Cip3_S3_L001_R2_001.fastq
        Mel1_S7_L001_R1_001.fastq
        Mel1_S7_L001_R2_001.fastq
        Mel2_S8_L001_R1_001.fastq
        Mel2_S8_L001_R2_001.fastq
        Mel3_S9_L001_R1_001.fastq
        Mel3_S9_L001_R2_001.fastq
        Pex1_S4_L001_R1_001.fastq
        Pex1_S4_L001_R2_001.fastq
        Pex2_S5_L001_R1_001.fastq
        Pex2_S5_L001_R2_001.fastq
        Pex3_S6_L001_R1_001.fastq
        Pex3_S6_L001_R2_001.fastq
```
# Different quantification methods overview
Libraries was quantified with:
* `kallisto` with cdna
* `salmon` with cdna
* `RSEM` with cdna, aligner - `bowtie`
* `HTSeq-count` with cdna, aligner - `bowtie`
* `RSEM` with full genome, aligner - `bowtie`
* `HTSeq-count` with full genome, aligner - `bowtie`

Overview tables:

Number of reads aligned (pseudoaligned):
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/compare_hits.png)
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/compare_hits_bars.png)
Note as Pex1_S4 library has more than 60% rRNA in it (rRNA-checking highlighted later), results vary massive depending on choosed reference (cdna/genome)

DE genes number:
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/compare_de.png)

Outliers number:
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/compare_outliers.png)

Independent filter treshold, genes under treshold:
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/compare_treshold.png)

Kallisto and salmon methods were choosed for final DE genes table building and GO enrichment.

# Download transcripts and annotations
Download and extract cdna, ncrna and gtf to ref folder
```sh
mkdir ref && cd ref
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/current/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cdna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa.gz
gunzip Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/current/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/ncrna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.ncrna.fa.gz
gunzip Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.ncrna.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/bacteria/current/gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.32.gtf.gz
gunzip Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.32.gtf.gz
```

# QC
```r
library(Rqc)
qc <- rqc(path = '~/path/to/reads', pattern = '.fastq', pair = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12), group = c(rep('con',6), rep('cip',6), rep('mel',6), rep('pex',6)), outdir = 'path/to/qc/output', file = 'ec2_qc.html')
```

# Check rRNA contamination
## Index ncrna
```sh
cd .. && mkdir kallisto
cd kallisto && mkdir ncrna
cd ..

cat > vars
export READS="path/to/fastqs"
export KALL="path/to/kallisto/output"
export REF="path/to/transcripts"
export R1="L001_R1_001.fastq"
export R2="L001_R2_001.fastq"
```
```sh
source vars

# get sample prefixes
for i in *fastq; do
    echo $i
done | cut -f 1-2 -d "_" | sort -u > $KALL/samples

cd $KALL

kallisto index -i ncrna/index $REF/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.ncrna.fa
```

## Run kallisto with ncrna
```sh
for i in `cat $READS/samples`; do
	kallisto quant \
	-i ncrna/index \
	-t 4 \
	-o ncrna/$i \
	--rf-stranded \
	$READS/${i}_$R1 $READS/${i}_$R2
done 
```

Tximport requires tx2gene file, we do it with python script
```sh
python tx2gene_file_creator.py $REF/ Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.ncrna.fa
mv $REF/tx2gene.csv $REF/tx2gene_rrna.csv
```

Import required libraries:
```r
library(tximport)
library(readr)
library(DESeq2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
library(tibble)
library(gplots)
library(ggplot2)
library(GenomicFeatures)
library(cowplot)
library(plotly)
library(goseq)
```

# Summarize
Make a table with total counts and rrna counts per sample (rrna counts obtained by pseudoalignment (`kallisto`) of reads with ncrna of E.Coli).
```r
table <- perFileInformation(qc)[c('filename', 'total.reads')]
table$filename <- gsub(pattern = '_L001_R[12]_001.fastq', replacement = '', x = table$filename)
table <- unique(table)
colnames(table)[names(table) == 'filename'] <- 'rowname'

# import kallisto rrna-counts
dir <- "path/to/kallisto/output/ncrna/"
run <- readLines("path/to/reads/samples")
files <- file.path(dir, run, "abundance.tsv")
names(files) <- run
all(file.exists(files))  # check if all files exist
tx2gene <- read.csv("path/to/transcripts/tx2gene_rrna.csv")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)

# add total counts for every sample to the table
rrna_counts <- colSums(txi$counts) %>% as_tibble %>% rownames_to_column
table <- merge(table, rrna_counts, by='rowname') %>% column_to_rownames(var = 'rowname')
table
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/rRNA_table.png)
```r
colnames(table)[names(table) == 'value'] <- 'rrna.counts'
table$rrna.perc <- round((table$rrna.counts / table$total.reads)*100, digits = 2) %>% paste0("%")
plot_ly(x = rownames(table), y=table$total.reads, name="total reads", type="bar") %>% add_trace(x=rownames(table), y=table$rrna.counts, name="rrna counts", type="bar") %>% layout(title = "rRNA in samples",yaxis = list(title=""), xaxis = list(title = ""))
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/rRNA.png)

# Index and quantify reads with kallisto and salmon
```sh
cd $KALL

kallisto index -i index $REF/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa

for i in `cat $READS/samples`; do
	kallisto quant \
	-i index \
	-t 4 \
	-o $i \
	--rf-stranded \
	$READS/${i}_$R1 $READS/${i}_$R2
done

cd $READS && mkdir salmon
echo export SALM="path/to/salmon/output" >> vars

source vars

cd $SALM

salmon index -t $REF/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa -i index --type quasi -k 31

for i in `cat $READS/samples`; do
	salmon quant \
	-i index \
	-l ISR \
	-1 $READS/${i}_$R1 -2 $READS/${i}_$R2 \
	-o ${i}
done 

python tx2gene_file_creator.py $REF Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cdna.all.fa
```

# Analysis:
## Fit model in `R`
```r
# Locate the directories containing files:
methods <- c("salmon", "kallisto")
for (method in methods){
  assign(x = paste0("dir_", method), value = paste0("path/to/reads", method))
}

# Create named vector pointing to the files:
run <- readLines(file.path(dir_kallisto, "samples"))
files_salmon <- file.path(dir_salmon, run, "quant.sf")
files_kallisto <- file.path(dir_kallisto, run, "abundance.tsv")
names(files_kallisto) <- run
names(files_salmon) <- run

# Check if all the files exist:
all(file.exists(files_kallisto))
all(file.exists(files_salmon))

# Associate transcripts with gene ID's for gene-level summarisation:
tx2gene <- read.csv(file.path(dir_kallisto, "tx2gene.csv"))

# Importing transcript-level estimates:
txi.salmon <- tximport(files_salmon, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
txi.kallisto <- tximport(files_kallisto, type = "kallisto", tx2gene = tx2gene, reader = read_tsv)

# C reating a table with column and row names:
sampleTable <- data.frame(treatment = factor(rep(c("con","cip","mel","pex"), each = 3)))
rownames(sampleTable) <- colnames(txi.kallisto$counts)

# Creating a DESeqDataSet for use with DESeq2:
dds.salmon <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~ treatment)
dds.kallisto <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~ treatment)
```

Perform a default analysis through the steps:
  1. Estimation of size factors: estimateSizeFactors
  2. Estimation of dispersion: estimateDispersions
  3. Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
We use both "parametric" and "local" types of fitting of dispersions to the mean intensity.
```r
dds.salmon_local <- DESeq(dds.salmon, fitType = "local")
dds.salmon_parametric <- DESeq(dds.salmon, fitType = "parametric")
dds.kallisto_local <- DESeq(dds.kallisto, fitType = "local")
dds.kallisto_parametric <- DESeq(dds.kallisto, fitType = "parametric")
```

## Basic visualisation
### Dispersions
Compare fitting types for `salmon` using dispersion plots:
```r
par(mfrow = c(1, 2))
plotDispEsts(dds.salmon_local)
title("Salmon - local")
plotDispEsts(dds.salmon_parametric)
title("Salmon - parametric")
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/disp_salm.png)

We choose local.
```r
dds.salmon <- dds.salmon_local
dds_list <- list()
dds_list[["salmon"]] <- dds.salmon_local
comment(dds_list[["salmon"]]) <- "salmon"
```

Compare fitting types for `kallisto` using dispersion plots:
```r
par(mfrow = c(1, 2))
plotDispEsts(dds.kallisto_local)
title("Kallisto - local")
plotDispEsts(dds.kallisto_parametric)
title("Kallisto - parametric")
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/disp_kal.png)

We choose local again.
```r
dds.kallisto <- dds.kallisto_local
dds_list[["kallisto"]] <- dds.kallisto_local
comment(dds_list[["kallisto"]]) <- "kallisto"
```

### PCA plots
Make PCA plots for `salmon` and `kalisto`:
```r
pca_list <- list()
for (dds in dds_list){
  pca_list[[comment(dds)]] <- rlog(object = dds, fitType = "local") %>% assay %>% as.data.frame %>% t %>% prcomp
  comment(pca_list[[comment(dds)]]) <- paste0("PCA ", comment(dds))
}

par(mfcol=c(1,2))
for (pca in pca_list) {
  plot(pca, main = comment(pca))
}
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PCA.png)

We keep the first 3 PCs and plot the samples in the 2D plane spanned by two principal components:
```r
PC <- 3
par(mfrow = c(2,3))
for (pca in pca_list){
  pca_x <- as.data.frame(pca$x)
  pca_x$treatment <- sampleTable$treatment
  prop <- summary(pca) %>% 
    as.list %>% .$importance %>% t %>% as_tibble %>% .[1:3,] %>% 
    dplyr::select(prop = starts_with("Proportion")) %>% .$prop *100
  for (i in 1:PC){
    pair <- combn(1:PC, m = 2)[,i]
    plot <- ggplot(pca_x, aes(x = pca_x[,pair[1]], y =  pca_x[,pair[2]], color = treatment)) +
      geom_point(size=5) +
      xlab(paste0("PC", pair[1], " :" ,round(prop[pair[1]]), "% variance", sep = "")) +
      ylab(paste0("PC", pair[2], " :" ,round(prop[pair[2]]), "% variance", sep = "")) +
      theme(legend.justification=c(1,0), legend.position=c(1,0)) +
      ggtitle(c(comment(pca), " PC", pair[1], " vs. PC", pair[2]) %>% paste(collapse = ""))
    print(plot)
  }
}
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PC12_salm.png) ![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PC13_salm.png) ![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PC23_salm.png)
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PC12_kal.png) ![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PC13_kal.png) ![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/PC23_kal.png)

### Heatmap and clustering
Make pheatmap showing the expression data of the 50 most expressed genes. The data is of log2 normalized counts.
```r
select_salm <- order(rowMeans(counts(dds.salmon, normalized=T)),decreasing = T)[1:50]
nt_salm <- normTransform(dds.salmon)
log2.norm.counts_salm <- assay(nt_salm)[select_salm,]

select_kal <- order(rowMeans(counts(dds.kallisto, normalized=T)),decreasing = T)[1:50]
nt_kal <- normTransform(dds.kallisto)
log2.norm.counts_kal <- assay(nt_kal)[select_kal,]

pheatmap(log2.norm.counts_salm, cluster_rows = F, show_rownames = T, cluster_cols = T, legend = T, main = "Pheatmap 50 most highly expressed genes Salmon")
pheatmap(log2.norm.counts_kal, cluster_rows = F, show_rownames = T, cluster_cols = T, legend = T, main = "Pheatmap 50 most highly expressed genes Kallisto")
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/pheat_salm.png)
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/pheat_kal.png)

### Euclidian distances
Heatmap of the Euclidian distances between the samples as calculated from the regularized log transformation:
```r
rld_salm <- rlog(dds.salmon, fitType = "local")
sampleDists_salm <- dist(t(assay(rld_salm))) #apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances
sampleDistMatrix_salm <- as.matrix(sampleDists_salm)
pheatmap(sampleDistMatrix_salm, clustering_distance_rows = sampleDists_salm, clustering_distance_cols = sampleDists_salm, main = "Euclidian distances between the samples, Salmon")


rld_kal <- rlog(dds.kallisto, fitType = "local")
sampleDists_kal <- dist(t(assay(rld_kal))) 
sampleDistMatrix_kal <- as.matrix(sampleDists_kal)
pheatmap(sampleDistMatrix_kal, clustering_distance_rows = sampleDists_kal, clustering_distance_cols = sampleDists_kal, main = "Euclidian distances between the samples, Kallisto")
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/dist_salm.png) ![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/dist_kal.png)

### MA-Plots
MA-Plot visualising the distribution of differentially expressed genes by plotting their normalized mean expressions against the Log Fold Change
```r
methods <- c("salmon", "kallisto")
treatments <- c("con", "cip", "mel", "pex")
contrast <- "treatment"
params <- c("alpha=0.05", "addMLE=TRUE")

# this function makes a list with results for evary quantification method and for every treatment in it
make_deseq2_results <- function(experiments=NULL, 
                                contrast=FALSE, 
                                factor_name=NULL, 
                                treatments=NULL, 
                                dds_file_prefix="dds",
                                params=NULL){
  results_list <- list()
  params <- paste(params, collapse=",")
  
  # if contrast=TRUE but no factor name or treatments provided
  if(contrast & (is.null(treatments) | is.null(factor_name))){
    stop("Please provide factor name and contrast treatments or set contrast=FALSE")
  }
  
  # if experiments specified:
  if(!is.null(experiments)){
    
    for (experiment in experiments){
      
      results_list[[experiment]] <- list()
      comment(results_list[[experiment]]) <- experiment
      
      # if contrast and treatments specified:
      if(contrast){
        for (i in 2:length(treatments)){
          res_name <- treatments[i]
          fun_string <- c("results(",
                   paste(c("object=", dds_file_prefix, experiment), collapse=""), ",",
                   paste(c("contrast=c('", factor_name, "','", treatments[1], "','", treatments[i], "')"), collapse=""), ",",
                   params,
                   ")") %>% paste(collapse="")
          results_list[[experiment]][[res_name]] <- eval(parse(text=fun_string))
          comment(results_list[[experiment]][[res_name]]) <- paste(c(treatments[1], " vs. ", treatments[i], " (", experiment, ")"), collapse="")
        }
      }
      
      # if no contrast specified
      else{
        res_name <- paste(c("res", experiment, collapse="_"))
        fun_string <- fun_string <- c("results(",
                   paste(c("object=", dds_file_prefix, experiment), collapse=""), ",",
                   paste(params, collapse=","), ",",
                   ")") %>% paste(collapse="")
        results_list[[res_name]] <- eval(parse(text=fun_string))
        comment(results_list[[res_name]]) <- experiment
      }
    }
  }
  
  # if no experiments:
  else{
    # if contrast and treatments specified:
    if(contrast){
        for (i in 2:length(treatments)){
          res_name <- paste(c("res", treatments[i]), collapse="_")
          fun_string <- c("results(",
                   paste(c("object=", dds_file_prefix), collapse=""), ",",
                   paste(c("contrast=c('", factor_name, "','", treatments[1], "','", treatments[i], "')"), collapse=""), ",",
                   params,
                   ")") %>% paste(collapse="")
          results_list[[res_name]] <- eval(parse(text=fun_string))
          comment(results_list[[res_name]]) <- paste(c(treatments[1], " vs. ", treatments[i]), collapse="")
        }
      }
    
    # if no experiments and no contrast
    else{
      fun_string <- c("results(",
                      dds_file_prefix, ",",
                      params, ")") %>% paste(collapse="")
      results_list[["res"]] <- eval(parse(text=fun_string))
    }
  }
  
  return(results_list)
}

res <- make_deseq2_results(experiments = methods,
                           contrast = TRUE,
                           factor_name = contrast,
                           treatments = treatments,
                           dds_file_prefix = "dds.",
                           params = params)

names(res)
for (method in res){
  print(names(method))
}

par(mfrow=c(2,3))
for (experiment in experiments){
  for (treatment in res[[experiment]]){
    plotMA(treatment, main=comment(treatment), ylim=c(-6,6))
  }
}
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/ma.png)

### Vulcano plots
Vulcano plot visualising DE genes by plotting the log2FoldChange vs. the adjusted p-value
```r
for (method in res){
  for (results in method){
    plot(data=results, padj*(-1) ~ log2FoldChange, 
         main = comment(results),
         xlab="log2FoldChange", ylab="-padj",
         subset = (results$padj>=0.05),
         pch=20, cex=0.2, col = "navy", 
         xlim = c(-6,6), ylim=c(-1,0))
    points(data=results, padj*(-1) ~ log2FoldChange,
           subset=(results$padj<0.05 & results$padj>=0.01), 
           pch=20, cex=0.2, col = "darkmagenta")
    points(data=results, padj*(-1) ~ log2FoldChange,
           subset=(results$padj<0.01), 
           pch=20, cex=0.2, col = "darkorange")
    legend(x="bottomleft",
           c("padj < 0.01", "padj < 0.05", "padj >= 0.05"), 
           pch = c(20, 20, 20), 
           col = c("darkorange", "darkmagenta", "navy"), 
           bty = "n", cex=1)
    grid()
  }
}
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/vulcano.png)

## Overlap in gene expression between groups
Compare the overlap in gene expression between the three treatments groups. We take only the genes with padj < 0.05 and |log2FoldChange| > 1. $\alpha =  0.05$

```r
overlap_allsalm <- list(cip = subset(res$salmon$cip, padj < 0.05 & abs(log2FoldChange) > 1) %>% rownames,
                        mel = subset(res$salmon$mel, padj < 0.05 & abs(log2FoldChange) > 1) %>% rownames, 
                        pex = subset(res$salmon$pex, padj < 0.05 & abs(log2FoldChange) > 1) %>% rownames)
overlap_upRegsalm <- list(cip = subset(res$salmon$cip, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) == 1) %>% rownames,
                      mel = subset(res$salmon$mel, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) == 1) %>% rownames,
                      pex = subset(res$salmon$pex, padj < 0.05 & abs(log2FoldChange) > 1 &sign(log2FoldChange) == 1) %>% rownames)
overlap_downRegsalm <- list(cip = subset(res$salmon$cip, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) != 1) %>% rownames,
                      mel = subset(res$salmon$mel, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) != 1) %>% rownames,
                      pex = subset(res$salmon$pex, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) != 1) %>% rownames)

overlap_allkal <- list(cip = subset(res$kallisto$cip, padj < 0.05 & abs(log2FoldChange) > 1) %>% rownames,
                        mel = subset(res$kallisto$mel, padj < 0.05 & abs(log2FoldChange) > 1) %>% rownames, 
                        pex = subset(res$kallisto$pex, padj < 0.05 & abs(log2FoldChange) > 1) %>% rownames)
overlap_upRegkal <- list(cip = subset(res$kallisto$cip, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) == 1) %>% rownames,
                      mel = subset(res$kallisto$mel, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) == 1) %>% rownames,
                      pex = subset(res$kallisto$pex, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) == 1) %>% rownames)
overlap_downRegkal <- list(cip = subset(res$kallisto$cip, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) != 1) %>% rownames,
                      mel = subset(res$kallisto$mel, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) != 1) %>% rownames,
                      pex = subset(res$kallisto$pex, padj < 0.05 & abs(log2FoldChange) > 1 & sign(log2FoldChange) != 1) %>% rownames)

venn(overlap_allsalm)
title("A. All DE genes (Salmon)")
venn(overlap_upRegsalm)
title("B. Up-regulated genes (Salmon)")
venn(overlap_downRegsalm)
title("C. Down-regulated genes (Salmon)")
venn(overlap_allkal)
title("A. All DE genes (Kallisto")
venn(overlap_upRegkal)
title("B. Up-regulated genes (Kallisto)")
venn(overlap_downRegkal)
title("C. Down-regulated genes (Kallisto)")
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/overlap.png)

## Compare kallisto and salmon results
What is the percentage of the same DE genes found by both? Padj < 0.05 and |log2FoldChange| > 1
For cip:
```r
res$salmon$cip <- res$salmon$cip %>% subset(padj < .05 & abs(log2FoldChange) > 1) # filter not DE genes
res$kallisto$cip <- res$kallisto$cip %>% subset(padj < .05 & abs(log2FoldChange) > 1) 
compare_cip <- 100 * sum(rownames(res$salmon$cip) %in% rownames(res$kallisto$cip)) / max(length(rownames(res$salmon$cip)), length(rownames(res$kallisto$cip)))
paste(round(compare_cip, 2), "%")
```
`92.69 %`

For mel:
```r
res$salmon$mel <- res$salmon$mel %>% subset(padj < .05 & abs(log2FoldChange) > 1)
res$kallisto$mel <- res$kallisto$mel %>% subset(padj < .05 & abs(log2FoldChange) > 1) 
compare_mel <- 100 * sum(rownames(res$salmon$mel) %in% rownames(res$kallisto$mel)) / max(length(rownames(res$salmon$mel)), length(rownames(res$kallisto$mel)))
paste(round(compare_mel, 2), "%")
```
`91.77 %`

For pex:
```r
res$salmon$pex <- res$salmon$pex %>% subset(padj < .05 & abs(log2FoldChange) > 1)
res$kallisto$pex <- res$kallisto$pex %>% subset(padj < .05 & abs(log2FoldChange) > 1) 
compare_pex <- 100 * sum(rownames(res$salmon$pex) %in% rownames(res$kallisto$pex)) / max(length(rownames(res$salmon$pex)), length(rownames(res$kallisto$pex)))
paste(round(compare_pex, 2), "%")
```
`92.45 %`

MA-plots for Salmon and Kallisto overlaping and not overlaping DE genes:
```r
# super puzzling code
kal.salm.overlap_cip <- subset(res$kallisto$cip, rownames(res$kallisto$cip) %in% rownames(res$salmon$cip))
kal.salm.not_overlap_cip <- subset(res$kallisto$cip, !(rownames(res$kallisto$cip) %in% rownames(res$salmon$cip)))
salm.kal.overlap_cip <- subset(res$salmon$cip, rownames(res$salmon$cip) %in% rownames(res$kallisto$cip))
salm.kal.not_overlap_cip <- subset(res$salmon$cip, !(rownames(res$salmon$cip) %in% rownames(res$kallisto$cip)))

par(mfrow = c(2,2))
plot(data = kal.salm.overlap_cip, log2FoldChange ~ baseMean, log = "x", pch=20, cex=0.3, col = "navy", main="A. Salmon vs. Kallisto (cip)", xlim=c(1, max(kal.salm.overlap_cip$baseMean)), ylim=c(-7, 7))
points(data = salm.kal.overlap_cip, log2FoldChange ~ baseMean, pch=20, cex=0.3, col = "indianred2")
points(data = kal.salm.not_overlap_cip, log2FoldChange ~ baseMean, pch=0, cex=1, col = "blue")
points(data = salm.kal.not_overlap_cip, log2FoldChange ~ baseMean, pch=5, cex=1, col = "black")

kal.salm.overlap_mel <- subset(res$kallisto$mel, rownames(res$kallisto$mel) %in% rownames(res$salmon$mel))
kal.salm.not_overlap_mel <- subset(res$kallisto$mel, !(rownames(res$kallisto$mel) %in% rownames(res$salmon$mel)))
salm.kal.overlap_mel <- subset(res$salmon$mel, rownames(res$salmon$mel) %in% rownames(res$kallisto$mel))
salm.kal.not_overlap_mel <- subset(res$salmon$mel, !(rownames(res$salmon$mel) %in% rownames(res$kallisto$mel)))

plot(data = kal.salm.overlap_mel, log2FoldChange ~ baseMean, log = "x", pch=20, cex=0.3, col = "navy", main="B. Salmon vs. Kallisto (mel)", xlim=c(1, max(kal.salm.overlap_mel$baseMean)), ylim=c(-7, 7))
points(data = salm.kal.overlap_mel, log2FoldChange ~ baseMean, pch=20, cex=0.3, col = "indianred2")
points(data = kal.salm.not_overlap_mel, log2FoldChange ~ baseMean, pch=0, cex=1, col = "blue")
points(data = salm.kal.not_overlap_mel, log2FoldChange ~ baseMean, pch=5, cex=1, col = "black")

kal.salm.overlap_pex <- subset(res$kallisto$pex, rownames(res$kallisto$pex) %in% rownames(res$salmon$pex))
kal.salm.not_overlap_pex <- subset(res$kallisto$pex, !(rownames(res$kallisto$pex) %in% rownames(res$salmon$pex)))
salm.kal.overlap_pex <- subset(res$salmon$pex, rownames(res$salmon$pex) %in% rownames(res$kallisto$pex))
salm.kal.not_overlap_pex <- subset(res$salmon$pex, !(rownames(res$salmon$pex) %in% rownames(res$kallisto$pex)))

plot(data = kal.salm.overlap_pex, log2FoldChange ~ baseMean, log = "x", pch=20, cex=0.3, col = "navy", main="C. Salmon vs. Kallisto (pex)", xlim=c(1, max(kal.salm.overlap_pex$baseMean)), ylim=c(-7, 7))
points(data = salm.kal.overlap_pex, log2FoldChange ~ baseMean, pch=20, cex=0.3, col = "indianred2")
points(data = kal.salm.not_overlap_pex, log2FoldChange ~ baseMean, pch=0, cex=1, col = "blue")
points(data = salm.kal.not_overlap_pex, log2FoldChange ~ baseMean, pch=5, cex=1, col = "black")
par(fig = c(0, 1, 0, 1), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(x=0, y=0,
       c("Shared genes - Salmon", "Shared genes - Kallisto", "IN Salmon, NOT IN Kallisto", "IN Kallisto, NOT IN Salmon"),
       pch = c(20, 20, 5, 0), col = c("indianred2", "navy", "black", "blue"), cex=1.5)
```
![](https://github.com/Perugolate/ec_rnaseq2/blob/master/plots/de_overlap.png)

## Export DE genes tables
Export tables with up- and down-regulated genes for every treatment. We take only significant for both salmon and kallisto genes:
```r
res$salmon$cip <- res$salmon$cip %>% as_tibble %>% rownames_to_column(var="GENEID")
res$kallisto$cip <- res$kallisto$cip %>% as_tibble %>% rownames_to_column(var="GENEID")

res_cip <- subset(res$kallisto$cip, res$kallisto$cip$GENEID %in% res$salmon$cip$GENEID)  # take overlaping genes (salmon and kallisto)

res$salmon$mel <- res$salmon$mel %>% as_tibble %>% rownames_to_column(var="GENEID")
res$kallisto$mel <- res$kallisto$mel %>% as_tibble %>% rownames_to_column(var="GENEID")

res_mel <- subset(res$kallisto$mel, res$kallisto$mel$GENEID %in% res$salmon$mel$GENEID)

res$salmon$pex <- res$salmon$pex %>% as_tibble %>% rownames_to_column(var="GENEID")
res$kallisto$pex <- res$kallisto$pex %>% as_tibble %>% rownames_to_column(var="GENEID")

res_pex <- subset(res$kallisto$pex, res$kallisto$pex$GENEID %in% res$salmon$pex$GENEID)

# add TXNAME column, sort by lfcMLE and export tables
tx2gene <- read.csv(file = "path/to/transcripts/tx2gene.csv")
res.cip <- merge(res_cip, tx2gene, by="GENEID") %>% arrange(desc(abs(lfcMLE))) %>% write.table(file="path/to/reads/cip.tsv", sep="\t", row.names=F)
res.mel <- merge(res_mel, tx2gene, by="GENEID") %>% arrange(desc(abs(lfcMLE))) %>% write.table(file="path/to/reads/mel.tsv", sep="\t", row.names=F)
res.pex <- merge(res_pex, tx2gene, by="GENEID") %>% arrange(desc(abs(lfcMLE))) %>% write.table(file="path/to/reads/pex.tsv", sep="\t", row.names=F)
```

Now we extend this tables with some data from UniProt (function, GO etc.) using python function annotation
```sh
python annotation.py path/to/reads/ cip.tsv cip_extended.tsv TXNAME ENSEMBLGENOME_TRS_ID --alt_key GENEID --alt_category GENENAME
python annotation.py path/to/reads/ mel.tsv mel_extended.tsv TXNAME ENSEMBLGENOME_TRS_ID --alt_key GENEID --alt_category GENENAME
python annotation.py path/to/reads/ pex.tsv pex_extended.tsv TXNAME ENSEMBLGENOME_TRS_ID --alt_key GENEID --alt_category GENENAME
```

# GO enrichment:
Prepare GO table from goa file:
```sh
cd ref
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa
python goa_to_goseq.py
```
Do hypergeometric test
```r
# make a data frame with all measured genes and DE genes (0 or 1)
measured_genes <- counts(dds.kallisto) %>% rownames

genes_cip_up <- as.integer(measured_genes %in% subset(res_cip, sign(log2FoldChange) == 1)$GENEID)
names(genes_cip_up) <- measured_genes

genes_cip_down <- as.integer(measured_genes %in% subset(res_cip, sign(log2FoldChange) != 1)$GENEID)
names(genes_cip_down) <- measured_genes

genes_mel_up <- as.integer(measured_genes %in% subset(res_mel, sign(log2FoldChange) == 1)$GENEID)
names(genes_mel_up) <- measured_genes

genes_mel_down <- as.integer(measured_genes %in% subset(res_mel, sign(log2FoldChange) != 1)$GENEID)
names(genes_mel_down) <- measured_genes

genes_pex_up <- as.integer(measured_genes %in% subset(res_pex, sign(log2FoldChange) == 1)$GENEID)
names(genes_pex_up) <- measured_genes

genes_pex_down <- as.integer(measured_genes %in% subset(res_pex, sign(log2FoldChange) != 1)$GENEID)
names(genes_pex_down) <- measured_genes

# make length data
txdb <- makeTxDbFromGFF(file="path/to/transcripts/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.32.gtf")
txsByGene <- transcriptsBy(txdb, "gene")
lengthData <- median(width(txsByGene)) 
lengthData <- subset(lengthData, names(lengthData) %in% measured_genes)

# make category mapping data frame
go <- read.table("path/to/transcripts/18.E_coli_MG1655_extracted.goa", sep = "\t") %>% arrange(V1)
go <- data.frame(go$V1, go$V2)
colnames(go) <- c("gene_id", "GO")

# fitting the probability weighting function (pwf)
pwf_cip_up <- nullp(DEgenes = genes_cip_up, bias.data = lengthData)
pwf_cip_down <- nullp(DEgenes = genes_cip_down, bias.data = lengthData)
pwf_mel_up <- nullp(DEgenes = genes_mel_up, bias.data = lengthData)
pwf_mel_down <- nullp(DEgenes = genes_mel_down, bias.data = lengthData)
pwf_pex_up <- nullp(DEgenes = genes_pex_up, bias.data = lengthData)
pwf_pex_down <- nullp(DEgenes = genes_pex_down, bias.data = lengthData)


go_cip_up_nobias <- goseq(pwf = pwf_cip_up, gene2cat = go, method = "Hypergeometric", test.cats = c("GO:BP"))
go_cip_down_nobias <- goseq(pwf = pwf_cip_down, gene2cat = go, method = "Hypergeometric", test.cats = c("GO:BP"))

go_mel_up_nobias <- goseq(pwf = pwf_mel_up, gene2cat = go, method = "Hypergeometric", test.cats = c("GO:BP"))
go_mel_down_nobias <- goseq(pwf = pwf_mel_down, gene2cat = go, method = "Hypergeometric", test.cats = c("GO:BP"))

go_pex_up_nobias <- goseq(pwf = pwf_pex_up, gene2cat = go, method = "Hypergeometric", test.cats = c("GO:BP"))
go_pex_down_nobias <- goseq(pwf = pwf_pex_down, gene2cat = go, method = "Hypergeometric", test.cats = c("GO:BP"))

write.table(subset(go_cip_up_nobias, go_cip_up_nobias$over_represented_pvalue < 0.05), file="path/to/reads/go_cip_up.tsv", sep="\t", row.names=F)
write.table(subset(go_cip_down_nobias, go_cip_down_nobias$over_represented_pvalue < 0.05), file="path/to/reads/go_cip_down.tsv", sep="\t", row.names=F)

write.table(subset(go_mel_up_nobias, go_mel_up_nobias$over_represented_pvalue < 0.05), file="path/to/reads/go_mel_up.tsv", sep="\t", row.names=F)
write.table(subset(go_mel_down_nobias, go_mel_down_nobias$over_represented_pvalue < 0.05), file="path/to/reads/go_mel_down.tsv", sep="\t", row.names=F)

write.table(subset(go_pex_up_nobias, go_pex_up_nobias$over_represented_pvalue < 0.05), file="path/to/reads/go_pex_up.tsv", sep="\t", row.names=F)
write.table(subset(go_pex_down_nobias, go_pex_down_nobias$over_represented_pvalue < 0.05), file="path/to/reads/go_pex_down.tsv", sep="\t", row.names=F)
```