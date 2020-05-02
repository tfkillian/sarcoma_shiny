## genetonic_conversions

# BiocManager::install("GeneTonic")

library("GeneTonic")
library("DESeq2")
# library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("SummarizedExperiment")
library("clusterProfiler")
library("GO.db")

# define your directory
myDirectory <- "~/tmp/sarcoma_shiny"

## dds object (from raw mouse data)
# countdata <- read.csv(paste0(myDirectory, "/data/rms_sk_combo.csv"),
#                       header = TRUE, row.names = 1)
# rownames(countdata) <- gsub("\\..*","", rownames(countdata))
# condition <- factor(c(rep("T0", 3), rep("T3", 4), rep("T5", 2), rep("KMR", 2)))
# coldata <- data.frame(row.names = colnames(countdata), condition)
# dds_1 <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata,
#                                 design = ~ condition)
# dds_1 <- estimateSizeFactors(dds_1)
# 
# ## annotation object
# anno_df <- data.frame(
#     gene_id = rownames(dds_1),
#     gene_name = mapIds(org.Mm.eg.db,
#                        keys = rownames(dds_1),
#                        column = "SYMBOL",
#                        keytype = "ENSEMBL"),
#     stringsAsFactors = FALSE,
#     row.names = rownames(dds_1)
# )
# 
# ## res object
# res_de <- readRDS(paste0(myDirectory, "/results/mouse_KMR_TO_res.rds"))
# rownames(res_de) <- gsub("\\..*","", rownames(res_de))
# 
# ## convert res object to enrich_object
# res_de_df <- as.data.frame(res_de)
# res_de_df$Gene <- rownames(res_de_df)
# res_de_df <- res_de_df %>%
#     mutate(ENTREZID = mapIds(org.Mm.eg.db, Gene, "ENTREZID", "ENSEMBL") %>% unname())
# 
# # Use simulated data set
# dataset <- res_de_df %>% 
#     filter(!is.na(padj)) %>% 
#     filter(!is.na(ENTREZID)) %>% 
#     filter(!duplicated(ENTREZID))
# 
# padj_thr <- 0.05
# log2FC_thr <- 0.05
# 
# geneList <- dataset$log2FoldChange
# names(geneList)<- dataset$ENTREZID
# genes <- dataset$ENTREZID[dataset$padj <= padj_thr & 
#                               abs(dataset$log2FoldChange) >= log2FC_thr]
# 
# res_GO <- enrichGO(gene          = genes,
#                    universe      = names(geneList),
#                    OrgDb         = org.Mm.eg.db,
#                    ont           = "BP", ### you can change this to another ontology
#                    pAdjustMethod = "BH",
#                    minGSSize     = 5,
#                    maxGSSize     = 500,
#                    pvalueCutoff  = 1,
#                    qvalueCutoff  = 1,
#                    readable      = TRUE)
# 
# res_enrich <- shake_enrichResult(res_GO)
# 
# saveRDS(dds_1, paste0(myDirectory, "/results/dds_1.rds"))
# saveRDS(anno_df, paste0(myDirectory, "/results/anno_df.rds"))
# saveRDS(res_de, paste0(myDirectory, "/results/res_de.rds"))
# saveRDS(res_enrich, paste0(myDirectory, "/results/res_enrich.rds"))