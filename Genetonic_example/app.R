## genetonic app

# BiocManager::install("mikelove/macrophage")
# BiocManager::install("GeneTonic")

library("GeneTonic")
library("macrophage")
library("DESeq2")
library("org.Hs.eg.db")
library("AnnotationDbi")

# example("GeneTonic")
## which will in the end run

## dds object
data("gse", package = "macrophage")
dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
rownames(dds_macrophage) <- substr(rownames(dds_macrophage), 1, 15)
dds_macrophage <- estimateSizeFactors(dds_macrophage)

## annotation object
anno_df <- data.frame(
    gene_id = rownames(dds_macrophage),
    gene_name = mapIds(org.Hs.eg.db,
                       keys = rownames(dds_macrophage),
                       column = "SYMBOL",
                       keytype = "ENSEMBL"),
    stringsAsFactors = FALSE,
    row.names = rownames(dds_macrophage)
)

## res object
data(res_de_macrophage, package = "GeneTonic")
res_de <- res_macrophage_IFNg_vs_naive

## res_enrich object
data(res_enrich_macrophage, package = "GeneTonic")
res_enrich <- shake_topGOtableResult(topgoDE_macrophage_IFNg_vs_naive)

GeneTonic(dds = dds_macrophage,
          res_de = res_de,
          res_enrich = res_enrich,
          annotation_obj = anno_df,
          project_id = "my_first_genetonic")
