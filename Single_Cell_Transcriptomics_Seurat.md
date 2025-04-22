Single Cell Transcriptomics with Seurat
================
Júlia Perera-Bel (MARData-BU, Hospital del Mar Research Institute)

-   <a href="#introduction-to-single-cell-transcriptomics"
    id="toc-introduction-to-single-cell-transcriptomics">Introduction to
    Single Cell Transcriptomics</a>
    -   <a href="#preprocessing-with-cell-ranger"
        id="toc-preprocessing-with-cell-ranger">Preprocessing with Cell
        Ranger</a>
        -   <a href="#qc" id="toc-qc">QC</a>
    -   <a href="#seurat" id="toc-seurat">Seurat</a>
-   <a href="#loading-data-and-initial-qc"
    id="toc-loading-data-and-initial-qc">Loading data and initial QC</a>
    -   <a href="#visualizing-qc" id="toc-visualizing-qc">Visualizing QC</a>
    -   <a href="#cell-filtering" id="toc-cell-filtering">Cell filtering</a>
-   <a href="#normalization-and-scaling"
    id="toc-normalization-and-scaling">Normalization and Scaling</a>
-   <a href="#dimentionality-reduction"
    id="toc-dimentionality-reduction">Dimentionality reduction</a>
-   <a href="#integration" id="toc-integration">Integration</a>
-   <a href="#cell-clustering-and-annotation"
    id="toc-cell-clustering-and-annotation">Cell clustering and
    annotation</a>
-   <a href="#differential-expression"
    id="toc-differential-expression">Differential expression</a>

**WARNING: THIS TUTORIAL IS STILL UNDER CONSTRUCTION.**

# Introduction to Single Cell Transcriptomics

scRNA-seq measures the abundance of mRNA molecules per cell. Extracted
biological tissue samples constitute the input for single-cell
experiments. Tissues are digested during single-cell dissociation,
followed by single-cell isolation to profile the mRNA per cell
separately. Plate-based protocols isolate cells into wells on a plate,
whereas droplet-based methods capture cells in microfluidic droplets.

The obtained mRNA sequence reads are mapped to genes and cells of origin
in raw data processing pipelines that use either cellular barcodes or
unique molecular identifiers (UMIs) and a reference genome to produce a
count matrix of cells by genes (Fig. 2a). You can find a detailed
overview of single cell sequencing in the [Single-Cell Best Practices
online
book](https://www.sc-best-practices.org/introduction/scrna_seq.html).
Here we will give a brief description of the preprocessing and consider
count matrices as the starting point for our analysis workflow of
scRNA-seq data.

## Preprocessing with Cell Ranger

Cell Ranger is the commercial software from 10x Genomics. It comrpises a
set of analysis pipelines that perform sample demultiplexing, barcode
processing, single cell 3’ and 5’ gene counting, V(D)J transcript
sequence assembly and annotation, and Feature Barcode analysis from
single cell data.

For this tutorial we will be using data generated with the [cellranger
count](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)
pipeline. The input you need to run cellranger count are the sequence
reads and a reference. The reference has a specific format. You can
download precomputed human and mouse references from the [10X
website](https://www.10xgenomics.com/support/software/cell-ranger/downloads).

``` r
cellranger count --id=$name --fastqs=$FASTQDIR --transcriptome=$REFDIR --sample=$base_name --create-bam=true --include-introns=true
```

The cellranger count pipeline outputs are in the pipestance directory in
the outs folder. List the contents of this directory with ls -1.

    ├── analysis
    ├── cloupe.cloupe
    ├── filtered_feature_bc_matrix
    ├── filtered_feature_bc_matrix.h5
    ├── metrics_summary.csv
    ├── molecule_info.h5
    ├── possorted_genome_bam.bam
    ├── possorted_genome_bam.bam.bai
    ├── raw_feature_bc_matrix
    ├── raw_feature_bc_matrix.h5
    └── web_summary.html

The **web_summary.html** contains a summary of the QC and results of the
experiment. You can also load the cloupe.cloupe file into the Loupe
Browser and start an analysis. The **outs/** directory contains the
outputs that can be used as input for software tools developed outside
of 10x Genomics, such as the Seurat R package.

### QC

It is common in droplet-based protocols that certain barcodes are
associated with ambient RNA instead of the RNA of a captured cell. This
happens when droplets fail to capture a cell. Many approaches exist to
assess whether a barcode likely corresponds to an empty droplet or not.

One simple method implemented in CellRanger is to examine the cumulative
frequency plot of the barcodes, in which barcodes are sorted in
descending order of the number of distinct UMIs with which they are
associated. This plot often contains a “knee” that can be identified as
a likely point of discrimination between properly captured cells and
empty droplet. While this “knee” method is intuitive and can often
estimate a reasonable threshold, not all cumulative histograms display
an obvious knee. Finally, the total UMI count associated with a barcode
may not, alone, be the best signal to determine if the barcode was
associated with an empty or damaged cell.

![web_summary](https://github.com/MARData-BU/Tutorials/raw/main/Images/web_summary.png)

There are several tools specifically designed to detect empty or damaged
droplets, or cells generally deemed to be of “low quality”. We recommend
importing all cells in R and filter empty droplets during the first
quality check, in which several filters and heuristics are applied to
obtain robust cells and genes. We also recommend to use
[SoupX](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831?login=true)
for quantifying the extent of the contamination and estimating
“background-corrected” cell expression profiles.

## Seurat

The next step after the generation of the count matrices with cellranger
count, is the data analysis. The R package Seurat is currently the most
popular software to do this. To start working with Seurat you can load
the data in two ways:

-   Using the barcodes/features/matrix files:

<!-- -->


    filtered_feature_bc_matrix
      ├── barcodes.tsv.gz
      ├── features.tsv.gz
      └── matrix.mtx.gz
      
    raw_feature_bc_matrix
      ├── barcodes.tsv.gz
      ├── features.tsv.gz
      └── matrix.mtx.gz

-   Using h5 files:

<!-- -->

    raw_feature_bc_matrix.h5
    filtered_feature_bc_matrix.h5

``` r
library(Seurat)

sc <- Read10X(file.path(dir,"raw_feature_bc_matrix")) #folder with the three barcodes, features and matrix files
sc <- Read10X_h5(file.path(dir,"raw_feature_bc_matrix.h5"))
```

# Loading data and initial QC

In this tutorial we will be using scRNASeq data from lamina propia from
the article by [Martin et
al](https://pmc.ncbi.nlm.nih.gov/articles/PMC7060942/). The data is
available under the accession GEO
[GSE134809](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134809).

``` r
# Download dataset
data_dir="GSE134809_RAW"
samples=list.files(data_dir)
seurat_list=list()

for (name in samples){
  cat(paste("\n Readding sample", name, "\n", sep = " "))
  temp_dir <- file.path(data_dir, name)
  # Now Read10X
  counts <- Read10X(data.dir = temp_dir)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = name, min.cells = 3, min.features = 200)
  # Add sample info
  seurat_obj$name <- name
  # calculate % of mithocondrial RNA
  seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") 
  # calculate % of ribosomal RNA
  seurat_obj$percent.rp <- PercentageFeatureSet(seurat_obj, pattern = "^RPL|^RPS") 
  
  # Store it
  seurat_list[[name]] <- seurat_obj
  
}

# Merge all Seurat objects into one
merged_seurat <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)],add.cell.ids = names(seurat_list))
```

Now we will need to add metadata

### Visualizing QC

While generating the `Seurat` object, there were already some quality
measures calculated for each cell, namely the total UMI counts per cell
(`nCount_RNA`) and the total number of detected features per cell
(`nFeature_RNA`).

-   nFeature: The number of unique genes detected in each cell.
    Low-quality cells or empty droplets will often have very few genes.
    Cell doublets or multiplets may exhibit an aberrantly high gene
    count.

-   nCount: Similarly, the total number of molecules detected within a
    cell (correlates strongly with unique genes)

We can also calculate other useful metrics and use for filtering cells
include:

-   Mitochondrial genes: If a cell membrane is damaged, it looses free
    RNA quicker compared to mitochondrial RNA, because the latter is
    part of the mitochondrion. A high relative amount of mitochondrial
    counts can therefore point to damaged cells.

-   Ribosomal genes: Are not rRNA (ribosomal RNA) but is mRNA that code
    for ribosomal proteins. They do not point to specific issues, but it
    can be good to have a look at their relative abundance. They can
    have biological relevance.

All of them can be plot in a violin plot and evaluate their distribution
per sample:

``` r
#head(seurat_obj@meta.data)

VlnPlot(merged_seurat, features=c("percent.mt"), pt.size = 0, group.by = "orig.ident",) + 
  geom_hline(yintercept = 10, size = 1, linetype = "dashed")+
  theme(legend.position = 'none', axis.title = element_text(size = 16), axis.text = element_text(size = 12))
```

![](Single_Cell_Transcriptomics_Seurat_files/figure-gfm/QC-1.png)<!-- -->

``` r
VlnPlot(merged_seurat, features=c("percent.rp"), pt.size = 0, group.by = "orig.ident",) + 
  theme(legend.position = 'none', axis.title = element_text(size = 16), axis.text = element_text(size = 12))
```

![](Single_Cell_Transcriptomics_Seurat_files/figure-gfm/QC-2.png)<!-- -->

``` r
VlnPlot(merged_seurat, features=c("nFeature_RNA"), pt.size = 0, group.by = "orig.ident",) + 
  geom_hline(yintercept = 300, size = 1, linetype = "dashed") + geom_hline(yintercept = 3000, size = 1, linetype = "dashed")+
  theme(legend.position = 'none', axis.title = element_text(size = 16), axis.text = element_text(size = 12))
```

![](Single_Cell_Transcriptomics_Seurat_files/figure-gfm/QC-3.png)<!-- -->

``` r
VlnPlot(merged_seurat, features=c("nCount_RNA"), pt.size = 0, group.by = "orig.ident",) + ylim(0,25000)+
  theme(legend.position = 'none', axis.title = element_text(size = 16), axis.text = element_text(size = 12))
```

![](Single_Cell_Transcriptomics_Seurat_files/figure-gfm/QC-4.png)<!-- -->

### Cell filtering

Based on the QC process we went through we can come to the following
conclusions:

-   There are no cells with very high mitochondrial gene counts.
-   There are some cells with a hemoglobin and low ribosomal counts, and
    these are probably erythrocytes.
-   There are some cells with a very low and very high number of
    features. These might point to non-informative cells and doublets
    respectively.
-   The ‘usual suspect’ MALAT1 sometimes makes up a large part of the
    counts per cell. As we do not see any other suggestions of
    dying/stressed cells, we leave it in.

In the M&M of the
[publication](https://www.nature.com/articles/s41598-020-64929-x#Sec7),
the authors describe that they have used a threshold of \< 8%
mitochondrial counts and \> 200 features per cell. To filter against
possible doublets, here, we also filter out cells with \> 5000 detected
features/cell. Filtering `Seurat` objects can be done with the `subset`
method for class `SeuratObject`:

``` r
seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 8)
```

To evaluate this did the trick we can visualize those parameters again
in a violin plot:

``` r
Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))
```

# Normalization and Scaling

# Dimentionality reduction

# Integration

# Cell clustering and annotation

# Differential expression
