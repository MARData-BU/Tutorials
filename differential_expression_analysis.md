---
layout: default_mod
title: <span class="title-underline">How to interpret your RNA-seq differential gene expression analysis results</span>
---

Pau Berenguer-Molins, Júlia Perera-Bel (MARData-BU, Hospital del Mar
Research Institute)
January 14, 2025

# 1. Introduction to Differential Expression Analyses in RNA-seq

Differential gene expression (DGE) analysis is a technique used in genomics to understand how genes behave under different conditions, such as treatment versus control.
This type of analysis identifies which genes are more or less expressed in response to specific conditions, such as a disease, drug treatment, or environmental change. By comparing gene activity between groups, scientists can uncover important biological insights, such as how diseases develop or how treatments work at the genetic level.

The results of a differential gene expression analysis provide a list of genes, each showing whether their activity has increased, decreased, or remained the same between the two groups. Understanding these changes can guide further research and contribute to the development of better medical treatments or enhanced biological understanding.

# 2. First approach to the data: dendrograms, PCAs and MDS plots

An initial approach to gaining an overview of the data is to interpret dendrograms, PCA, and MDS plots. Despite being generated using different methodologies, all three types of plots are interpreted similarly: the closer two samples are, the more similar they are in terms of gene expression. Conversely, the farther apart they are, the less similar they are.

Samples are typically color-coded according to a specific condition, providing an initial view of how the samples behave. This also helps identify batch effects or the influence of other variables apart from the variable of interest.

![dendrograms](https://github.com/MARData-BU/Tutorials/raw/main/Images/pca_dendrogram_mds.png)

# 3. Interpreting the first results: differentially expressed genes and filters

The results from a DGE analysis include a list of genes with corresponding p-values, fold-change (FC), log2(FC), and adjusted p-values. The number of differentially expressed genes depends on the thresholds you decide to apply. Naturally, the stricter the filter, the more reliable the results.

-   **p-value**: p-value obtained with the moderated t-test ([Ritchie et al., 2015](#Ritchie_2015)). It indicates the statistical significance of the result.
-   **adjusted p-value**: adjusted p-value for multiple comparisons ([Benjamini and Hochberg, 1995](#Benjamini_Hochberg_1995))(Benjamini and Hochberg 1995).
-   **fold-change**: the value represents the 2^ transformation of the log2 fold-change field.
-   **log2 fold-change**: the log2 difference between the mean expression levels of the two groups (i.e., the mean of the first group minus the mean of the second group).

# 4. Visual results and plots

Results from a DGE analysis can be viewed in several plots.

## 4.1 Volcanoplots

Volcano plots display the results of a DGE analysis by combining statistical significance and the magnitude of gene expression changes.

-   **X-axis**: represents the log2 fold-change (log2(FC)). Genes on the right are upregulated (positive log2(FC)), while genes on the left are downregulated (negative log2(FC)).
-   **Y-axis**: represents the negative logarithm of the p-value or adjusted p-value (-log(p-value)). This transformation ensures all values are positive. Genes higher on the plot are more statistically significant.

Typically, genes are represented as dots, and their colors indicate whether they meet certain thresholds for log2(FC) and p-value (or adjusted p-value). This provides a clear visual summary of significant and biologically relevant genes.

![volcanoplot](https://github.com/MARData-BU/Tutorials/raw/main/Images/volcanoplot.png)

## 4.2 Heatmaps

Heatmaps visualize the expression levels of a set of genes using a color-coded matrix.

-   **Rows**: represent individual genes.
-   **Columns**: represent samples.
-   **Clustering**: both rows and columns are typically clustered based on expression patterns, allowing groups of genes and samples with similar behaviors to be easily identified.

Heatmaps are useful for discovering expression patterns across groups of samples or for highlighting differences between conditions. Genes displayed on the heatmap can either be manually selected or filtered based on specific thresholds (log2(FC) and/or (adjusted) p-value).

![heatmap](https://github.com/MARData-BU/Tutorials/raw/main/Images/heatmap.png)

# 5. References

- <a id="Ritchie_2015"></a>Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015 Apr;43(7):e47. doi: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007).
- <a id="Benjamini_Hochberg_1995"></a>Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *Journal of the Royal Statistical Society: Series B (Methodological)*. 1995;57:289–300. doi: [10.1186/gb-2010-11-3-r25](https://doi.org/10.1186/gb-2010-11-3-r25).
