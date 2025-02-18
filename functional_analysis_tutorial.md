---
layout: default_mod
title: <span class="title-underline">Understand and interpret functional analyses</span>
---

Pau Berenguer-Molins, Júlia Perera-Bel (MARData-BU, Hospital del Mar
Research Institute)
January 10, 2025

# 1. Introduction to Functional Analysis in Gene Expression Studies

In gene expression studies, functional analysis is a critical approach
for understanding the biological meaning behind differential gene
expression (DGE) results. After identifying genes that show significant
changes in expression levels under different conditions (e.g: case vs
control, treatment vs vehicle), one might wonder what these changes
imply for cellular function, biological pathways, and/or disease
mechanisms. This is where functional analysis comes in: it enables us to
contextualize gene expression patterns within biological processes,
pathways, and gene sets, helping to connect the molecular data to
potential biological or clinical insights.

Functional analysis allows to interpret patterns in gene sets, moving
beyond individual gene analysis to capture how groups of related genes
contribute to specific functions. These analyses provide a broader
picture, as genes rarely act in isolation; instead, they work together
within pathways and network to drive diverse biological processes.
Through functional analysis, we can uncover which pathways, biological
processes, or molecular functions are most affected by the changes
observed in gene expression, offering a more holistic view of the
underlying biology.

In this tutorial, we’ll focus on two primary methods for functional
analysis in the context of DGE results:

-   **Pre-Ranked Gene Set Enrichment Analysis (GSEA)**: the pre-ranked
    gene set enrichment analysis (GSEA) is a method used to determine
    whether predefined sets of genes show statistically significant,
    concordant differences between two biological states. Unlike many
    other techniques, GSEA does **not** require a cutoff for
    differentially expressed genes, which allows it to capture broader,
    more subtle patterns. Here, genes resulting from a DGE analysis are
    ranked (in our case, using *-*(p.val) \* signFC) and tested for
    functional enrichment against predefined gene sets belonging to a
    gene set collection.

-   **Over-Representation Analysis (ORA)**: in this case, the genes do
    **not** necessarily come from a DGE analysis (although they usually
    do). This analysis tests whether a specific category of genes is
    statistically overrepresented over a specific group of genes (e.g:
    genes fulfilling a specific statistical threshold from a DGE
    analysis). ORA is widely used for its straightforward interpretation
    and ability to highlight the most significantly impacted pathways or
    processes, as it usually focuses on genes that already fulfill
    statistical thresholds.

By following this tutorial, you’ll learn not only how to implement these
methods programmatically but also how to interpret the results through
numerical metrics and graphical representations. Each analysis technique
offers unique insights, and by the end of this tutorial, you’ll be able
to apply these methods effectively to gain a deeper understanding of the
functional implications of your gene expression data.

# 2. Databases

In order to run a functional analysis, one needs to provide information
regarding the pathways (and the genes within it) that will be tested.
These pathways need to be in “.gmt” format. These files can be either
manually created or else downloaded. One common webpage in which gene
sets can be browsed, viewed and downloaded (in gmt format, as well as
other formats) is the [Molecular Signatures
Database](https://www.gsea-msigdb.org/gsea/msigdb), which contains human
and murine gene sets from several databases. Some examples are:

-   **c5.bp**: gene sets derived from the Biological Process Gene
    Ontology (GO) ([Köhler et al., 2021](#HPO_2021)).
-   **c2.cp.kegg**: canonical Pathways gene sets derived from the KEGG
    pathway database ([Kanehisa, 2019](#KEGG_2019), [Kanehisa and Goto, 2000](#KEGG_2000), [Kanehisa et al., 2022](#KEGG_2022)).
-   **h.all**: hallmark gene sets summarize and represent specific
    well-defined biological states or processes and display coherent
    expression. These gene sets were generated by a computational
    methodology based on identifying overlaps between gene sets in other
    MSigDB collections and retaining genes that display coordinate
    expression ([Liberzon et al. 2015](#Hallmark_2015)).

![MSigDb](https://github.com/MARData-BU/Tutorials/raw/main/Images/msigdb.png)

# 3. Statistics behind the functional analyses

## 3.1 GSEA

GSEA is based in a permutation approach to assess whether a predefined
set of ranked genes shows a statistically significant, coordinated
expression pattern.

1.  **Gene ranking**: genes are ranked by a metric that reflects their
    association with the condition of interest. We will rank the genes
    by **-**(p.val) \* signFC.
2.  **Enrichment Score (ES)**: GSEA calculates an enrichment score (ES)
    that reflects the degree to which genes in the gene set are
    overrepresented at the top (upregulated) or bottom (downregulated)
    of the ranked gene list. GSEA calculates the ES by walking down the
    ranked list of genes, increasing a running-sum statistic when a gene
    is in the gene set and decreasing it when it is not. The magnitude
    of the increment depends on the correlation of the gene with the
    phenotype. The ES is the maximum deviation from zero encountered in
    walking the list. A positive ES indicates gene set enrichment at the
    top of the ranked list; a negative ES indicates gene set enrichment
    at the bottom of the ranked list. The leading-edge subset are those
    genes in the gene set that appear in the ranked list at, or before,
    the point where the running sum reaches its maximum deviation from
    zero.
3.  **Significant testing**: GSEA permutates the sample labels (e.g:
    case vs control, treated vs vehicle) multiple times to generate a
    null distribution of ES values for each gene set, allowing for a
    comparison with the observed ES.
4.  **Normalized Enrichment Score (NES)**: the normalized enrichment
    score (NES) is the primary statistic for examining gene set
    enrichment results. The NES accounts for differences in gene set
    size and in correlations between gene sets and the expression
    dataset. Therefore, the NES can be used to compare analysis results
    across gene sets. As this statistic is based on the gene set
    enrichment scores for all dataset permutations, changing the
    permutation method, the number of permutations, or the size of the
    expression dataset will affect the NES. NES=(actual ES)/mean(ES
    against all permutations of the dataset)
5.  **False Discovery Rate (FDR)**: the false discovery rate (FDR) is
    the estimated probability that a gene set with a given NES
    represents a false positive finding. An FDR of 25% indicates that
    the result is likely to be valid 3 out of 4 times (1 out of 4 will
    be a false positive). The FDR is a ratio of two distributions: (1)
    the actual enrichment score versus the enrichment scores for all
    gene sets against all permutations of the dataset and (2) the actual
    enrichment score versus the enrichment scores of all gene sets
    against the actual dataset. For example, if you analyze four gene
    sets and run 1000 permutations, the first distribution contains 4000
    data points and the second contains 4.
6.  **Nominal P value**: the nominal p value estimates the statistical
    significance of the enrichment score for a single gene set. However,
    when you are evaluating multiple gene sets, you must correct for
    gene set size and multiple hypothesis testing. Because the p value
    is not adjusted for either, it is of limited value when comparing
    gene sets.
7.  **Multiple Hypothesis Testing**: after obtaining p-values for each
    gene set, GSEA applies methods to control for multiple hypothesis
    testing, typically using the false discovery rate (FDR).

## 3.2 ORA

ORA is based on the Fisher’s exact test to determine if a predefined
gene set is overrepresented among a specific set of genes (usually
up/down regulated genes given a specific threshold).

1.  **Gene selection**: ORA requires a list of genes, which are usually
    defined using a threshold of p-value and/or log2(FC) in a
    differential expression analysis. Only these genes are tested for
    enrichment.
2.  **Contingency table construction**: for each gene set, ORA builds a
    2x2 table to compare the number of genes in the gene set that are
    differentially expresed versus the number of genes not in the set,
    where **A** is the number of differentially expressed (DE) genes in
    the gene set, **B** is the number of DE genes **not** in the gene
    set, **C** is the number of non-DE genes in the gene set, and **D**
    is the number of non-DE genes **not** in the gene set.

<!-- -->

    ##                              In GeneSet Not in GeneSet
    ## Differentially Expressed              A              B
    ## Not Differentially Expressed          C              D

3.  **Fisher’s Exact Test**: ORA uses Fisher’s exact test to determine
    if the overlap between the DE genes and the gene set is greater than
    the expected by chance.
4.  **Multiple hypothesis testing**: after obtaining p-values for each
    gene set, ORA applies methods to control for multiple hypothesis
    testing, typically using the false discovery rate (FDR).

# 4. How do I interpret my functional analysis results?

As a summary, results can be divided in two groups: images and Excel
files. While images are the most user-friendly and visual way to see
your results, Excel files will contain all your results (regardless the
statistical significance).

## 4.1 GSEA

### 4.1.1 Excel file

As the GSEA is derived from a DGE analysis, there will be one Excel file
for each comparison performed (GSEA.Database.Group1.vs.Group2.xlsx),
being one sheet for the upregulated pathways, and another for the
downregulated. Each sheet contains the following columns:

-   **ID**: gene set name. For MSigDB gene sets, the description can be
    browsed in the gene set page on the [MSigDb
    website](https://www.gsea-msigdb.org/gsea/msigdb).
-   **setSize**: number of genes in the gene set after filtering out
    those genes not in the expression dataset.
-   **enrichmentScore**: enrichment score for the gene set; that is, the
    degree to which this gene set is overrepresented at the top or
    bottom of the ranked list of genes in the expression dataset. This
    score **cannot be used to compare analysis resutls across gene
    sets**.
-   **NES**: normalized enrichment score; that is, the enrichment score
    for the gene set after it has been normalized across analyzed gene
    sets. It accounts for differences in gene set size and in
    correlations between gene sets and the expression dataset;
    therefore, **the normalized enrichment scores (NES) can be used to
    compare analysis results across gene sets**.
-   **pvalue**: nominal p value; that is, the statistical significance
    of the enrichment score. The nominal p value is **not** adjusted for
    gene set size or multiple hypothesis testing; therefore, it is of
    limited use in comparing gene sets. A lower bound of 1e-10 is used
    for estimating P-values.
-   **p.adjust**: false discovery rate-adjusted p value; that is, the
    estimated probability that a gene set with a given NES represents a
    false positive finding.
-   **rank**: the position in the ranked list at which the maximum
    enrichment score occurred. The more interesting gene sets achieve
    the maximum enrichment score near the top or bottom of the ranked
    list; that is, the rank at max is either very small or very large.
-   **leading_edge**: displays the three statistics used to define the
    leading edge subset:
    -   Tags: the percentage of gene hits before (for positive ES) or
        after (for negative ES) the peak in the running enrichment
        score. This gives an indication of the percentage of genes
        contributing to the enrichment score.
    -   List: the percentage of genes in the ranked gene list before
        (for positive ES) or after (for negative ES) the peak in the
        running enrichment score. This gives an indication of where in
        the list the enrichment score is attained.
    -   Signal: the enrichment signal strength that combines the two
        previous statistics. If the gene set is entirely within the
        first positions in the list, then the signal strength is maximal
        or 100%. If the gene set is spread throughout the list, then the
        signal strength decreases towards 0%.
-   **core_enrichment**: genes that contribute to the leading-edge
    subset within the gene set. This is the subset of genes that
    contributes most to the enrichment result ordered by
    -log(p.val)\*signFC. Hence, in positively enriched gene sets (NES
    \> 0) the genes are ordered in decreasing order of importance, and
    in negatively enriched gene sets (NES \< 0) in increasing order of
    importance.

![Enrichment_score](https://github.com/MARData-BU/Tutorials/raw/main/Images/enrichment_score.jpg)
\### Images

The most common plots for the GSEA are the following:

-   **Barplots**: enriched functional pathways are displayed in the Y
    axis on the left. Bars represent the normalized enrichment score
    (NES) and are colored according to the adjusted p-value. NES
    indicates the distribution of functional pathway genes across a list
    of ranked genes (in this case, ranked by -log(p.val)\*signFC). A
    positive NES indicates an enrichment in the first condition of the
    comparison (cases), whether a negative value corresponds to an
    enrichment in the second condition of the comparison (controls).

![Barplot](https://github.com/MARData-BU/Tutorials/raw/main/Images/GSEA_barplot_example.png)

-   **Dotplots**: enriched functional pathways are displayed in the Y
    axis on the left. The dot’s color shows the FDR of each term
    involved in the analysis, and the dot’s size indicates how many
    genes from that pathway were identified in the analysis. The
    horizontal axis (gene ratio) quantifies the proportion of the
    pathway’s genes in the analysis relative to its total gene count.

![Dotplot](https://github.com/MARData-BU/Tutorials/raw/main/Images/GSEA_dotplot_example.png)

-   **Enrichment maps**: nodes represent gene-sets and edges represents
    mutual overlap. Highly redundant gene-sets are grouped together as
    clusters. The nodes’ color corresponds to the adjusted p-value,
    whereas its size represents the gene count that contributes to the
    enrichment of that pathway.

![Enrichment_map](https://github.com/MARData-BU/Tutorials/raw/main/Images/GSEA_enrichment_map_example.png)

## 4.2 ORA

### 4.2.1 Excel file

In this case, results do not necessarily come from a DGE analysis, but
are a list of genes. As mentioned, it is common to filter genes
resulting from a DGE analysis via the (adjusted) p-value and/or the
log2(FC) to obtain these lists. Thus, there will be one Excel file per
gene list provided. Each excel contains the following columns:

-   **ID**: Gene set name. For MSigDB gene sets, the description can be
    browsed in the gene set page on the [MSigDb
    website](https://www.gsea-msigdb.org/gsea/msigdb).
-   **GeneRatio**: number of genes in the gene set as a fraction of
    genes in the input gene list.
-   **BgRatio**: background ratio. Size of the gene set as a fraction of
    the total number of background genes (ie. universe).
-   **pvalue**: nominal p value; that is, the statistical that the
    overlap between the DE genes and the gene set is greater than the
    expected by chance. In other words, it tests wether the GeneRatio is
    greater than expected when compared to the BgRatio. The nominal p
    value is **not** adjusted for gene set size or multiple hypothesis
    testing; therefore, it is of limited use in comparing gene sets.
-   **p.adjust**: false discovery rate-adjusted p value. It represents
    the probability that the observed overlap between the DE genes and
    the gene set is a false positive result.
-   **qvalue**: alternative FDR estimation by qvalue R package.
-   **geneID**: genes in the input gene list present in the gene set.

### 4.2.2 Images

As ORA does not directly com from a DGE analysis, the Enrichment Map
(which provides graphical information for the up and downregulated
pathways) is not a possibility. The most common plots for the ORA are
the following:

-   **Dotplots**: enriched functional pathways are displayed in the Y
    axis on the left. The dot’s color shows the FDR of each term
    involved in the analysis, and the dot’s size indicates how many
    genes from that pathway were identified in the analysis. The
    horizontal axis (gene ratio) quantifies the proportion of the
    pathway’s genes in the analysis relative to its total gene count.

![Dotplot](https://github.com/MARData-BU/Tutorials/raw/main/Images/GSEA_dotplot_example.png)

-   **Enrichment maps**: nodes represent gene-sets and edges represents
    mutual overlap. Highly redundant gene-sets are grouped together as
    clusters. The nodes’ color corresponds to the adjusted p-value,
    whereas its size represents the gene count that contributes to the
    enrichment of that pathway.

![Enrichment_map](https://github.com/MARData-BU/Tutorials/raw/main/Images/GSEA_enrichment_map_example.png)

# 5. References

- <a id="KEGG_2019"></a>Kanehisa, Minoru. Toward understanding the origin and evolution of cellular organisms. *Protein Science*. 2019;28(11):1947-1951. doi: [10.1002/pro.3715](https://doi.org/10.1002/pro.3715).

- <a id="KEGG_2022"></a>Kanehisa, Minoru, Miho Furumichi, Yoko Sato, Masayuki Kawashima, and Mari Ishiguro-Watanabe. KEGG for taxonomy-based analysis of pathways and genomes. *Nucleic Acids Research*. 2022;51(D1):D587-D592. doi: [10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963).

- <a id="KEGG_2000"></a>Kanehisa, Minoru, and Susumu Goto. KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research*. 2000;28(1):27-30. doi: [10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27).

- <a id="Hallmark_2015"></a>Liberzon, Arthur, Chet Birger, Helga Thorvaldsdóttir, Mahmoud Ghandi, Jill P. Mesirov, and Pablo Tamayo. The Molecular Signatures Database Hallmark Gene Set Collection. *Cell Systems*. 2015;1(6):417-425. doi: [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004).

- <a id="HPO_2021"></a>Köhler, Sebastian, Nicolas Matentzoglu, Michael Gargano. The Human Phenotype Ontology in 2021. *Nucleic Acids Research*. 2021;49(50):D687-D692. doi: [10.1093/nar/gkab1028](https://doi.org/10.1093/nar/gkab1028).
