# Single-cell RNA-seq and Analysis of Drosophila Olfactory Projection Neurons

This repository holds code for analysis of single-cell RNA-sequencing data from Drosophila olfactory projection neurons.

# Quick Start for ICIM

ICIM code is in the Python module sct.

A simple example of usage is analysis/ICIM_Example.ipynb.

It boils down to:

```
import sct

myICIM = sct.ICIM(X, df)

myICIM.calc()

marker_genes = myICIM.get_all_markers()
```

# Organization

Code for preprocessing sequence data (mapping reads to the genome and counting reads mapping to genes) is contained in the pipeline directory.

Code for visualization and analysis, including reproducing figures shown in the paper, is in the analysis directory.

Preprocessed sequence data and various intermediate files used during analysis are in the data directory.

Raw sequence data can be obtained from the Sequence Read Archive (accession GSE100058).

# Iterative Clustering for Identifying Markers (ICIM)

We introduced an unsupervised machine-learning algorithm for identifying informative genes for separating cell types in single-cell RNA-seq data, which we call ICIM.

Input: reduced count matrix X (gene x cell) and full count matrix df (gene x cell). The reduced count matrix may optionally be prefiltered to remove genes from consideration by ICIM. In practice, we use log-transformed counts.

Output: list of genes that distinguish populations.

You can use the list of genes for further dimensionality reduction and clustering.

Adjustable parameters are:

Related to filtering for informative genes for cell type identification:

* N = number of overdispersed genes considered at each step
* correlation_cutoff = minimum Pearson correlation for identifying correlated genes
* min_hits = minimum number of correlated genes required to keep a gene
* exclude_max = number of top expressing cells that are excluded when calculating robust correlation (for robustness to outliers)
* dropout_rate_low = minimum fraction of cells that a gene must be absent from to keep it
* dropout_rate_high = maximum fraction of cells that a gene must be absent from to keep it

Related to termination:

* stop_condition = {"linkage_dist", "num_cells"} = termination condition
* N_stop = number of cells in a subpopulation below which iteration stops, ignored unless stop_condition = "num_cells"
* linkage_dist_stop = linkage distance above which iteration stops, ignored unless stop_condition = "linkage_dist"
