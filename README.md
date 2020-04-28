# sarcoma_shiny
This repository contains Shiny web apps for visualizaling RNAseq data analysis results

* Shiny RNAseq viewer: This app displays 3 various plots from an RNAseq analysis of human and mouse sarcoma cells. These plots include, MA plot, Volcano plot and a count plot (of a selected gene, for each sample). These plots are generated dynamically plotted by user click input. There are also an html tables displaying the results of the analysis which can be queried, which interact with the display plots. There are additional tabs which display experimental design, PCA of samples, GO terms and KEGG pathways of the GSEA analysis of the data, and genes relating to those GO terms and KEGG pathways appears in corresponding volcano plots.

* GeneTonic viewer: this app produces numerous plots, including GO enrichment and gene networks, but can only display the results from one DE comparison.
