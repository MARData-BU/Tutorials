---
layout: default_mod
title: <span class="title-underline">Single Cell Transcriptomics with Seurat</span>
---

Júlia Perera-Bel (MARData-BU, Hospital del Mar
Research Institute)
April 15, 2025

**WARNING: THIS TUTORIAL IS STILL UNDER CONSTRUCTION.**

# Introduction to Single Cell Transcriptomics

## 

## Preprocessing with Cell Ranger

## Seurat

# Loading data and initial QC

### Visualizing QC per cell and gene

While generating the `Seurat` object, there were already some quality measures calculated for each cell, namely the total UMI counts per cell (`nCount_RNA`) and the total number of detected features per cell (`nFeature_RNA`). We can plot those in a violin plot and evaluate their distribution per sample:

```
Seurat::VlnPlot(seu, features = c("nCount_RNA",
                                  "nFeature_RNA"))
```


### Cell filtering

Based on the QC process we went through we can come to the following conclusions:

- There are no cells with very high mitochondrial gene counts.
- There are some cells with a hemoglobin and low ribosomal counts, and these are probably erythrocytes.
- There are some cells with a very low and very high number of features. These might point to non-informative cells and doublets respectively. 
- The 'usual suspect' MALAT1 sometimes makes up a large part of the counts per cell. As we do not see any other suggestions of dying/stressed cells, we leave it in. 

In the M&M of the [publication](https://www.nature.com/articles/s41598-020-64929-x#Sec7), the authors describe that they have used a threshold of < 8% mitochondrial counts and > 200 features per cell. To filter against possible doublets, here, we also filter out cells with > 5000 detected features/cell. Filtering `Seurat` objects can be done with the `subset` method for class `SeuratObject`:

```
seu <- subset(seu, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 8)
```

To evaluate this did the trick we can visualize those parameters again in a violin plot:

```
Seurat::VlnPlot(seu, features = c("nFeature_RNA",
                                  "percent.mito"))
```

```
#| echo: false
saveRDS(seu, "seu_day1-2.rds")
```

# Normalization and Scaling

# Dimentionality reduction

# Integration

# Cell clustering and annotation

# Differential expression








