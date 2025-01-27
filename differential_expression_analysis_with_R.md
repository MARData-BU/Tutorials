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

Once you have all packages installed, simply load them running the "library" function:

```
library("gplots") # change "gplots" for any other package
```

# 3. First steps: reading and filtering your counts table

How to read your counts table will depend on its file type. If you have an xlsx file, you can use read.xlsx(*path_to_file*) from openxlsx package. If you have a csv or tsv file, you can simply read the Counts table using read.table(*path_to_file*, sep = *specify separation in file*). This table of counts should have samples as columns and genes as rows. If that is not your case, you can always transpose the table by using t(*table_of_counts*).

Your table of counts may contain relevant information for each of the genes, such as the chromosome, the start and end position, the length and the strand. If that is the case, save this information in another R object, as you may want to retrieve it after.

```
# Read table of Counts
counts <- read.table(file="path_to_file/Counts_example.txt"), sep="\t", header=T, dec=".", stringsAsFactors = F)

# Print the first 2 rows of Counts table
head(counts, 2)
         Geneid  Chr   Start     End Strand Length KO_A1 KO_A2 KO_A3 KO_B1 KO_B2 KOB3 WT_A1 WT_A2 WT_A3 WT_B1 WT_B2
1 4933401J01Rik chr1 3143476 3144545      +   1070     0     0     0     0     0    0     0     0     0     0     0
2       Gm26206 chr1 3172239 3172348      +    110     0     0     0     0     0    0     0     0     0     0     0

# Split the gene information and counts information into two separate objects
Annot.RNAseq <- counts[,c(1:6)] # gene information
counts.m <- counts[,7:length(counts)] # Counts
```

Before running any analysis, it is advisable to perform some checks over the data.

-   Check the total counts per sample.
-   Check sequencing depth consistency by comparing the maximum and minimum values of total counts.

If the sequencing depth is reasonably consistent across the RNA samples, then the simplest and most robust approach to differential expression is to use limma-trend approach, where counts are converted to logCPM (log2-counts-per-million). This approach will usually work well if the ratio of the largest library size to the smallest is not more than about 3-fold.

```
sample.totals <- apply(counts.m, 2, sum) # calculate the sum (total counts) for each sample (each column)
range(sample.totals) # retrieve the range of total counts in samples
max(sample.totals)/min(sample.totals)
```

**WARNING: THIS TUTORIAL IS STILL UNDER CONSTRUCTION.**
