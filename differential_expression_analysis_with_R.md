---
layout: default_mod
title: <span class="title-underline">How to run your own RNA-seq differential expression analysis from a Counts table with R</span>
---

Pau Berenguer-Molins, Júlia Perera-Bel (MARData-BU, Hospital del Mar
Research Institute)
January 23, 2025

**WARNING: THIS TUTORIAL IS STILL UNDER CONSTRUCTION.**

# 1. Introduction to Differential Expression Analyses in RNA-seq

Differential gene expression (DGE) analysis is a technique used in genomics to understand how genes behave under different conditions, such as treatment versus control.
This type of analysis identifies which genes are more or less expressed in response to specific conditions, such as a disease, drug treatment, or environmental change. By comparing gene activity between groups, scientists can uncover important biological insights, such as how diseases develop or how treatments work at the genetic level.

The results of a differential gene expression analysis provide a list of genes, each showing whether their activity has increased, decreased, or remained the same between the two groups. Understanding these changes can guide further research and contribute to the development of better medical treatments or enhanced biological understanding.

# 2. Before you start: what files do you need and what is required in your computer?

As stated in the title of this tutorial, you should already have a Counts table (in txt, csv format, or similar), where **genes** are the rows and **samples** are the columns. This tutorial has been developed for R/RStudio. If you do not have them installed, you can download the latest versions [here](https://posit.co/download/rstudio-desktop/).

Once you have R and RStudio in your computer, you should install the following libraries:

-   For plot generation: gplots, RColorBrewer, ggplot2, ggrepel, [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html).
-   Data analysis: [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), compareGroups, corrplot, [limma](https://bioconductor.org/packages/release/bioc/html/limma.html), [sva](https://www.bioconductor.org/packages/release/bioc/html/sva.html).
-   MARData-BU generated packages: [QualityGraphs](https://github.com/MARData-BU/QualityGraphs), [VennPlots](https://github.com/MARData-BU/VennPlots), [AnalysisFunctions](https://github.com/MARData-BU/AnalysisFunctions).

Those without a link can be easily installed running the following command on your console:

```
install.packages("gplots") # change "gplots" for any other package
```

Or else you can manually install them by clicking in "Packages" on your Environments section (usually on the bottom-right of your RStudio), going to "Install" and looking for the package you wish to install.

![package_install](https://github.com/MARData-BU/Tutorials/raw/main/Images/package_install.png)

Those package with a link on them (from Bioconductor) can be installed running the following command:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sva") # change "sva" for any other package
```

Lastly, to install MARData-BU packages, you will be required to have "devtools" package installed and loaded in your R/RStudio session, and then run the following command:

```
install.packages("devtools") # run this only the 1st time; once installed, you will just need to load it with the next command
library(devtools)
install_github("MARData-BU/QualityGraphs") # change "QualityGraphs" for any other MARData-BU package
```

Once you have all packages installed, simply load them running the *library()* function:

```
library("gplots") # change "gplots" for any other package
```

# 3. First steps

## 3.1 Reading the table of counts

How to read your counts table will depend on its file type. If you have an xlsx file, you can use *read.xlsx(path_to_file)* from openxlsx package. If you have a csv or tsv file, you can simply read the Counts table using *read.table(path_to_file, sep = specify separation in file)*. This table of counts should have samples as columns and genes as rows. If that is not your case, you can always transpose the table by using *t(table_of_counts)*.

Your table of counts may contain relevant information for each of the genes, such as the chromosome, the start and end position, the length and the strand. If that is the case, save this information in another R object, as you may want to retrieve it after.

```
# Read table of Counts
counts <- read.table(file="path_to_file/Counts_example.txt", sep="\t", header=T, dec=".", stringsAsFactors = F)
rownames(counts) <- counts$Geneid # ensure gene names are set as rownames

# Print the first 2 rows of Counts table
head(counts, 2)
                    Geneid  Chr   Start     End Strand Length KO_A1 KO_A2 KO_A3 KO_B1 KO_B2 KO_B3 WT_A1 WT_A2 WT_A3 WT_B1 WT_B2
4933401J01Rik 4933401J01Rik chr1 3143476 3144545      +   1070     0     0     0     0     0    0     0     0     0     0     0
Gm26206             Gm26206 chr1 3172239 3172348      +    110     0     0     0     0     0    0     0     0     0     0     0

# Split the gene information and counts information into two separate objects
Annot.RNAseq <- counts[,c(1:6)] # gene information
counts.m <- counts[,7:length(counts)] # Counts
```

## 3.2 Initial quality check

Before running any analysis, it is advisable to perform some checks over the data.

-   Check the total counts per sample.
-   Check sequencing depth consistency by comparing the maximum and minimum values of total counts.

If the sequencing depth is reasonably consistent across the RNA samples, then the simplest and most robust approach to differential expression is to use limma-trend approach, where counts are converted to logCPM (log2-counts-per-million). This approach will usually work well if the ratio of the largest library size to the smallest is not more than about 3-fold. In this case, the tutorial will focus on limma package ([Ritchie et al., 1995](#limma)).

```
sample.totals <- apply(counts.m, 2, sum) # calculate the sum (total counts) for each sample (each column)
range(sample.totals) # retrieve the range of total counts in samples
[1] 80730737 95226607
max(sample.totals)/min(sample.totals)
[1] 1.179558
```

## 3.3 Counts table filtering

In order remove the effect those genes that are lowly expressed (either genes with a low expression, or genes that are expressed in few samples) it is advisable to filter out those genes that have 10 or less counts in the N of the smallest group in your analysis. In our case, we have 3 KO_A samples, 3 KO_B samples, 3 WT_A samples and 2 WT_B samples, so it would be advisable to remove those genes with 10 or less counts in less than 2 samples. Again considering our case, we move from 56,791 genes to 17,423 genes after filtering.

```
dim(counts.m) # check the dimensions of the counts table
[1] 56791    11
keep <- rowSums(counts.m>10) >= 2 # keep those genes with more than 10 counts in 2 or more samples
counts.m.f <- counts.m[keep,] # filter the counts table
dim(counts.m.f) # check the dimensions of the new counts table
[1] 17432    11
annot.m.f <- Annot.RNAseq[match(rownames(counts.m.f), Annot.RNAseq$Geneid),] # filter the annotations table
all.equal(rownames(counts.m.f), annot.m.f$Geneid) # check that the gene names are in the same order in the counts and annotation tables
[1] TRUE
```

## 3.4 Data normalization

For this next step, you will need to have package edgeR ([Robinson, McCarthy and Smyth, 2010](#edgeR1), [McCarthy, Chen and Smyth,2012](#edgeR2), [Chen, Lun and Smyth, 2016](#edgeR3)) installed and loaded. Data normalization involves:

-     **1**: transforming the table of counts in a DGEList object (edgeR object).
-     **2**: calculating scaling factors for each sample to convert raw library sizes (total number of counts) into effective library sizes with *calcNormFactors()* function. In our case, we will use the trimmed mean of M values (TMM) method ([Robinson and Oshlack, 2010](#TMM)), which computes normalization factors that represent sample-specific biases given their total counts. These factors are then multiplied by the library size to generate the effective library size.
-     **3**: compute the counts per million (CPM) and log2CPM with *cpm()* function. The *prior.count* argument is required when log-transforming the CPM. This value is the average count to be added to each observation (each gene per each sample) to avoid taking log of 0.

```
library(edgeR)
d <- DGEList(counts=counts.m.f) # transform table of counts into DGEList
Norm.Factor <- calcNormFactors(d, method="TMM") # compute normalization factors
cpm.matrx <- cpm(Norm.Factor, log=T, prior.count=3) # calculate log2CPM
cpm.matrx.nonlog <- cpm(Norm.Factor, log=F) # calculate CPM
```

## 3.5 QC plots

In order to assess the overall behavior of the samples and check whether there is any batch effect - or else the biological effect can already be seen - several plots can be generated to assess the samples given the expression of all the genes (already filtered). It is advisable to create a dataframe in which you have the samples linked to their known conditions (both the condition that wants to be assessed and other conditions that may be affecting the results, such as sex, age, BMI, batch, technician, etc.). In this example we will only focus in our experimental condition. Another possible option rather than creating a dataframe is to directly create a vector with the conditions. This is less advisable, as it is easier to mismatch the samples with their conditions.

Please note that the order in which the conditions are displayed must match the sample order in the *cpm.matrx* object.

```
targets
   Sample Condition
1   KO_A1      KO_A
2   KO_A2      KO_A
3   KO_A3      KO_A
4   KO_B1      KO_B
5   KO_B2      KO_B
6   KO_B3      KO_B
7   WT_A1      WT_A
8   WT_A2      WT_A
9   WT_A3      WT_A
10  WT_B1      WT_B
11  WT_B2      WT_B

all.equal(targets$Sample ,colnames(cpm.matrx)) # check that sample names in cpm.matrx and targets objects are in the same order
[1] TRUE

c(rep("KO_A", 3), rep ("KO_B", 3), rep ("WT_A", 3), rep("WT_B", 2))
 [1] "KO_A" "KO_A" "KO_A" "KO_B" "KO_B" "KO_B" "WT_A" "WT_A" "WT_A" "WT_B" "WT_B"
```

The recommended plots for gaining an overview of the data are dendrograms, PCAs, and MDS plots. Despite being generated using different methodologies, all three types of plots are interpreted similarly: the closer two samples are, the more similar they are in terms of gene expression. Conversely, the farther apart they are, the less similar they are. Samples are color-coded according to a specific condition, providing an initial view of how the samples behave. This also helps identify batch effects or the influence of other variables apart from the variable of interest.

```
library(QualityGraphs)
library(RColorBrewer)

# DENDROGRAMS
clusterdend(est_noctrls=cpm.matrx, conditions=targets$Condition, picname="Condition", resDir=file.path(resultsDir, "QC_Plots"))

# PCAs
makePCA(est_noctrls=cpm.matrx, conditions=t![dendrograms](https://github.com/MARData-BU/Tutorials/raw/main/Images/pca_dendrogram_mds.png)
argets$Condition, picname="Condition", resDir=file.path(resultsDir, "QC_Plots"), dist = 2)

# 2D PCAs
makePCA.2D(est_noctrls=cpm.matrx, conditions=targets$Condition, picname="Condition", resDir=file.path(resultsDir, "QC_Plots"), dist = 0)

# MDS  Plots
COLOR <- c('KO_A'='#1b9e77', 'KO_B' = '#d95f02', 'WT_A' = '#67639f', 'WT_B' = '#e7298a')
png(file.path(resultsDir, "MDS_Condition.png"), units="in",  width=12, height=8, res=200)
conditions <- factor(targets$Condition)  #to be used as block-variable
plotMDS(cpm.matrx, label=targets$Sample,
        col=COLOR[conditions])
legend("topleft", legend=levels(conditions), pch=15, col=COLOR, ncol=1)
dev.off()

```

As you can see on the plots below, sample *KO_A2* appears as an outlier, which indicates that this sample displays a gene expression quite different than the rest of the samples. In some situations, this could lead to the elimination of this sample from the analysis. You can also see that samples from the WT_B group tend to cluster together, as well as WT_A samples and KO_B samples.

![dendrograms](https://github.com/MARData-BU/Tutorials/raw/main/Images/pca_dendrogram_mds_example.png)

# 4. Differential expression analysis

```
library(limma)
library(edgeR)



```

# References

- <a id="limma"></a>Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015;43(7):e47. doi: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007).
- <a id="edgeR1"></a>Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*. 2010;26(1):139–140. doi: [10.1093/bioinformatics/btp616](https://doi.org/10.1093/bioinformatics/btp616).
- <a id="edgeR2"></a>McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Research*. 2012;40(10):4288–4297. doi: [10.1093/nar/gks042](https://doi.org/10.1093/nar/gks042).
- <a id="edgeR3"></a>Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Research*. 201
- <a id="TMM"></a>Robinson MD, Oshlack A. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biology*. 2010;11(3):R25. doi: [10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25).




**WARNING: THIS TUTORIAL IS STILL UNDER CONSTRUCTION.**
