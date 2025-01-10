---
layout: default_mod
title: Online functional analysis tutorial with g:Profiler - MARData-BU
---

Pau Berenguer-Molins, Júlia Perera-Bel (MARData-BU, Hospital del Mar
Research Institute)
January 7, 2025

# 1 Introduction to Functional Analysis in Gene Expression Studies

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

By following this tutorial, you’ll learn how to run and interpret a functional
analysis given your set of results in g:Profiler.

# 2 g:Profiler online tool

g:Profiler is an online tool that allows users to perform functional enrichment
analyses integrating many databases, organisms and analysis options. It was first
released in 2007 [Reimand et al., 2007](#gprofiler_2007) and has received many updates
[Kolberg et al., 2023](#gprofiler_2023). It provides a user-friendly interface
that facilitates the analysis without the requirement of a prior bioinformatic knowledge.

The user needs to provide a whitespaced-separated gene list (you can directly
copy the genes from an Excel file, for example) and set the organism that will be
used for the analysis (Human as default), as well as whether the query is ordered
or not, which will direct the analysis towards a GSEA or an ORA, respectively.
In the former case, elements are to be arranged in **decreasing** order of importance
(e.g: down-regulated to up-regulated). Results of ordered queries in g:Profiler
should not be treated as p-values. Instead, users should only infer whether genes
belonging to a term are evenly distributed across the query or primarily located at the top.

![gprofiler_overview](https://github.com/MARData-BU/Tutorials/raw/main/Images/gprofiler_overview.png)

Users can also define several advanced options, such as:

-   **All results**: whether to show all results or only the statistically significant ones (default).
-   **Measure underrepresentation**: measure underrepresentation instead of the default enrichment analysis (default).
-   **No evidence codes**: whether to show (default) or not the evidence codes. If not shown, querying and rendering is faster, but interactions will not be available in the result.
-   **Statistical domain scope**: whether to run the analysis on annotated genes, all known genes or custom genes.
-   **Significance threshold**: method for multiple testing correction, either gSCS (default), Bonferroni or FDR.
-   **User threshold**: significance threshold (0.05 by default).
-   **Numeric IDs treated as**: default prefix for fully numeric IDs.

The data sources that g:Profiler uses for the analysis can be (un)selected as well, or else provide a customized gmt file.

![gprofiler_options](https://github.com/MARData-BU/Tutorials/raw/main/Images/gprofiler_options.png)

## 2.1 Random example results

This random example was generated with the "random example" option in the "Query" section and using the default parameters of g:Profiler version e111_eg58_p18_f463989d on 1/7/2025, 11:23:07 AM. Results are displayed in both an image and a table. Image displays the statistically significant results (Y axis; -log10(p-adjusted)) for each database analysed (X axis). The table below shows the same information as the image, but in a table format and displaying only 8 pathways.  

![gprofiler_example_result](https://github.com/MARData-BU/Tutorials/raw/main/Images/gprofiler_example_result.png)

The detailed results section displays all pathways for each database. It allows the user to (un)select the pathways that are highlighted and appear on the table in the "Overview" section. Pathways are ordered by significance. A heatmap on the right shows the genes that are part of each pathway and are coloured according to the evidence codes.

![gprofiler_gomf_example](https://github.com/MARData-BU/Tutorials/raw/main/Images/gprofiler_gomf_example.png)

The "GO Context" section displays the relationships between the pathways in those databases that are hierarchical. Lastly, "Query Info" section displays all parameters set for the specific query.

# 3 References

- <a id="gprofiler_2007"></a>Reimand J, Kull M, Peterson H, Hansen J, Vilo J. g:Profiler--a web-based toolset for functional profiling of gene lists from large-scale experiments. *Nucleic Acids Research*. 2007 Jul;35(Web Server issue):W193-200. doi: [10.1093/nar/gkm226](https://doi.org/10.1093/nar/gkm226). Epub 2007 May 3. PMID: [17478515](https://pubmed.ncbi.nlm.nih.gov/17478515); PMCID: [PMC1933153](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1933153/).

- <a id="gprofiler_2023"></a>Kolberg L, Raudvere U, Kuzmin I, Adler P, Vilo J, Peterson H. g:Profiler-interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update). *Nucleic Acids Research*. 2023 Jul 5;51(W1):W207-W212. doi: [10.1093/nar/gkad347](https://doi.org/10.1093/nar/gkad347). PMID: [37144459](https://pubmed.ncbi.nlm.nih.gov/37144459); PMCID: [PMC10320099](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10320099/).
