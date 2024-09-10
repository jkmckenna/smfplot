# smfplot - This repo has been merged into [smftools](https://github.com/jkmckenna/smftools) and is accessible through the preprocessing, tools, and plotting modules
A python package to visualize single locus single molecule footprinting data contained within an anndata object produced by [smftools](https://github.com/jkmckenna/smftools)

**Current biologic questions include, but are not limited to:**
1) Site specific transcription factor binding kinetics, both in vitro and in permeabilized nuclei.
2) Locus specific nucleosome remodeling kinetics in permeabilized nuclei.
3) Thermodynamic modeling of transcriptional states as a function of transcription factor binding occupancy in permeabilized nuclei.
4) Long range cis-regulatory interactions throughout differentiation and during maintainance of cell state in permeabilized nuclei.
5) How enzymatic activities of transcriptional coactivators/corepressors kinetically alter chromatin microstates to regulate transcription in permeabilized nuclei.

## Table of Contents
- [User Guide](https://github.com/jkmckenna/smfplot/edit/main/README.md#user-guide)
  - [Intended usage](https://github.com/jkmckenna/smfplot/edit/main/README.md#intended-usage)
  - [Basic data structures](https://github.com/jkmckenna/smfplot/edit/main/README.md#basic-data-structures)
  - [Types of analyses](https://github.com/jkmckenna/smfplot/edit/main/README.md#types-of-analyses)
  - [Basic Usage](https://github.com/jkmckenna/smfplot/edit/main/README.md#basic-usage)
    - [From anndata objects to visualization using Plot_SMF_h5ad.py](https://github.com/jkmckenna/smfplot/edit/main/README.md#from-anndata-objects-to-visualization-using-plot_smf_h5adpy)
      - [Plotting basic heatmaps](https://github.com/jkmckenna/smfplot/edit/main/README.md#plotting-basic-heatmaps)
      - [Plotting pairwise correlation heatmaps](https://github.com/jkmckenna/smfplot/edit/main/README.md#plotting-pairwise-correlation-heatmaps)
      - [Plotting clustered heatmaps](https://github.com/jkmckenna/smfplot/edit/main/README.md#plotting-clustered-heatmaps)
      - [Higher dimensional analyses](https://github.com/jkmckenna/smfplot/edit/main/README.md#higher-dimensional-analyses)
  - [Performance Notes](https://github.com/jkmckenna/smfplot/edit/main/README.md#performance-notes)
  - [In the works](https://github.com/jkmckenna/smfplot/edit/main/README.md#in-the-works)

## User Guide
This repository currently contains python scripts and jupyter notebooks demonstrating basic analyses and visualizations centered around the anndata object. These scripts will eventually be consolidated into a python package for ease of use.

## Basic data structure
The core data structure is the [anndata](https://github.com/scverse/anndata) object. This object handles complex data nicely and is being actively developed and maintained by a solid community of developers in the single-cell genomics field. We have abstracted the usage of [scanpy](https://github.com/scverse/scanpy) and anndata objects to contain a Read X Position matrix analagous to the typical Cell X Gene matrix used in single cell RNA sequencing analyses. While many operations can be directly run on the anndata object, [anndata can be interfaced with pytorch](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/annloader.html) to enable more powerful modeling. This enables the use of CUDA if a NVIDIA GPU is available.

## Basic analyses
1) Detection of pairwise methylation states through linkage maps. When this is applied to large windows (low to moderate kb range), this detects pairwise accessibility of regions within the window (ex enhancer-promoter coaccessibility). When this is applied to smaller windows (sub kb), this can look at protein binding cooperativity and general protein binding footprints. Analyses can be visualized with an all versus all positional probability matrix (Commonly used in MicroC/HiC analyses) or as a one versus all positional plot (Commonly used in 4C analyses) depending on your question.
2) Clustering of methylation states (kmeans, hierarchical, dbscan). When applied to small windows (sub 500bp windows), this can look at the proportion of DNA molecules in a population that are bound at a given position. When applied to moderate size windows (500bp - 2000bp), prevalence of higher complexity binding microstates can be assayed.
3) Scanpy derived analyses including: principle component anlysis, dimensionality reduction, computing neighborhood graphs, leiden/louvain clustering, and UMAP. A means of analyzing loci with high complexity microstates. Facilitates the projection of metadata into PCA space and UMAP space to enable high dimensional analyses.

## Basic Usage
### General Dataflow:
**Input requirements ->** 
1) An anndata object derived from [smftools](https://github.com/jkmckenna/smftools)

### Plotting basic heatmaps

### Plotting pairwise correlation heatmaps
**80bp pairwise unmethylated data of a dCas9 bound to DNA as profiled by Hia5 m6A SMF**
![](https://github.com/jkmckenna/smfplot/blob/main/images/15min_Hia5_dCas9_bound_dsDNA_m6A_plus_strand_100bp_window_pairwise_heatmap.png?raw=true)

**1500bp pairwise unmethylated data of a dCas9 bound to DNA as profiled by Hia5 m6A SMF**
![](https://github.com/jkmckenna/smfplot/blob/main/images/dCas9_invitro_m6A_plus_strand_1500bp_window_pairwise_unmethylated_heatmap.png?raw=true)

**3000bp pairwise all versus all unmethylated linkage data of a an enhancer-promoter pair at a single allelic locus in F1 hybrid B6 NK cells**
![](https://github.com/jkmckenna/smfplot/blob/main/images/NKG2A_F1_B6_hybrid_enhancer_promoter_positional_all_vs_all_pairwise_methylation_heatmap.png)

**3000bp pairwise one versus all unmethylated linkage data of a an enhancer-promoter pair at a single allelic locus in F1 hybrid B6 NK cells**
![](https://github.com/jkmckenna/smfplot/blob/main/images/NKG2A_F1_B6_hybrid_enhancer_promoter_positional_one_vs_all_pairwise_methylation_linkage_plot.png)

### Plotting clustered heatmaps
**200bp kmeans clustered methylation data of a dCas9 bound to DNA as profiled by Hia5 m6A SMF**
![](https://github.com/jkmckenna/smfplot/blob/main/images/15min_Hia5_dCas9_bound_dsDNA_m6A_plus_strand_200bp_window_kmeans_clustered_heatmap.png)

**3000bp kmeans and hierarchical clustering of unmethylated reads by position data of a an enhancer-promoter pair at a single allelic locus in F1 hybrid B6 NK cells**
![](https://github.com/jkmckenna/smfplot/blob/main/images/NKG2A_F1_B6_hybrid_enhancer_promoter_position_by_read_kmeans_heatmap.png)
![](https://github.com/jkmckenna/smfplot/blob/main/images/NKG2A_F1_B6_hybrid_enhancer_promoter_position_by_read_hierarchical_heatmap.png)

### Higher dimensional analyses

## In the works
**1) Consolidating code into a python package and providing accompanying notebook tutorials to plot SMF data contained in anndata object**
I am currently composing scripts that will be compiled into a python package to be distributed through PyPI. This package will be imported into a Jupyter notebook to be used for tutorials on data visualization with the python package.

The package takes an anndata object produced from the CLI application as input and performs steady state and kinetic analyses/visualizations. These analyses include, but are not limited to:
 1) Inferring identity of proteins bound to a locus through an integration of structural and motif-based modeling derived from footprints.
 2) Calculating the fraction of a binding sites occupied across a population of single molecules of DNA.
 3) Single molecule cooccupancy analyses across proximal and distal binding sites to assess cooperativity of protein-binding.
 4) Long-range (ones to tens of kb) single molecule coaccessibility analyses to assess functional linkage of cis-regulatory elements.
 5) Bulk occupancy analyses over populations of single molecules (for visualization in genome track browsers, such as [IGV](https://igv.org/)). Essentially provides an interpretation analogous to high-coverage ATACseq, but at a single locus.

