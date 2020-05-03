## genetonic_conversions

# BiocManager::install("GeneTonic")

library("GeneTonic")
library("DESeq2")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("AnnotationDbi")
library("SummarizedExperiment")
library("clusterProfiler")
library("GO.db")
library("dplyr")

# define your directory
myDirectory <- "~/tmp/sarcoma_shiny"

########################### convert raw mouse data #############################

## dds object 
countdata <- read.csv(paste0(myDirectory, "/data/rms_sk_combo.csv"),
                      header = TRUE, row.names = 1)
rownames(countdata) <- gsub("\\..*","", rownames(countdata))
condition <- factor(c(rep("T0", 3), rep("T3", 4), rep("T5", 2), rep("KMR", 2)))
coldata <- data.frame(row.names = colnames(countdata), condition)
dds_1 <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata,
                                design = ~ condition)
dds_1 <- estimateSizeFactors(dds_1)

## annotation object
anno_df <- data.frame(
    gene_id = rownames(dds_1),
    gene_name = mapIds(org.Mm.eg.db,
                       keys = rownames(dds_1),
                       column = "SYMBOL",
                       keytype = "ENSEMBL"),
    stringsAsFactors = FALSE,
    row.names = rownames(dds_1)
)

## res object
res_de <- readRDS(paste0(myDirectory, "/results/mouse_KMR_TO_res.rds"))
rownames(res_de) <- gsub("\\..*","", rownames(res_de))

## convert res object to enrich_object
res_de_df <- as.data.frame(res_de)
res_de_df$Gene <- rownames(res_de_df)
res_de_df <- res_de_df %>%
    mutate(ENTREZID = mapIds(org.Mm.eg.db, Gene, "ENTREZID", "ENSEMBL") %>% unname())

# Use simulated data set
dataset <- res_de_df %>%
    filter(!is.na(padj)) %>%
    filter(!is.na(ENTREZID)) %>%
    filter(!duplicated(ENTREZID))

padj_thr <- 0.05
log2FC_thr <- 0.05

geneList <- dataset$log2FoldChange
names(geneList)<- dataset$ENTREZID
genes <- dataset$ENTREZID[dataset$padj <= padj_thr &
                              abs(dataset$log2FoldChange) >= log2FC_thr]

res_GO <- enrichGO(gene          = genes,
                   universe      = names(geneList),
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP", ### you can change this to another ontology
                   pAdjustMethod = "BH",
                   minGSSize     = 5,
                   maxGSSize     = 500,
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   readable      = TRUE)

res_enrich <- shake_enrichResult(res_GO)

saveRDS(dds_1, paste0(myDirectory, "/results/dds_2.rds"))
saveRDS(anno_df, paste0(myDirectory, "/results/anno_df_2.rds"))
saveRDS(res_de, paste0(myDirectory, "/results/res_de_2.rds"))
saveRDS(res_enrich, paste0(myDirectory, "/results/res_enrich_2.rds"))

############################## human data conversion ###########################
countdata_no_row <-  read.csv(paste0(
    myDirectory, "/data/human-RMS-RNAseq_selection.csv"), header = TRUE)
countdata_no_row$gene <- as.character(countdata_no_row$gene)
countdata_no_row <- countdata_no_row %>%
    mutate(gene_name = mapIds(org.Hs.eg.db, gene, "ENSEMBL", "SYMBOL") %>%
               unname()) %>% dplyr::select(gene, gene_name, everything())

# Remove duplicates based on gene columns
countdata_1 <- countdata_no_row[!duplicated(countdata_no_row$gene_name), ]
countdata_no_na <- countdata_1[!is.na(countdata_1$gene_name), ]
countdata <- countdata_no_na %>% dplyr::select(-gene)
rownames(countdata) <- countdata$gene_name
countdata <- countdata %>% dplyr::select(-gene_name)
names(countdata)[2] <- "RMS2028"

################ need to convert gene row names to ENS

## Assign condition 'coldata' 
condition <- factor(c(rep("RMS", 16), rep("ctrl", 5)))
coldata <- data.frame(row.names=colnames(countdata), condition)
dds_1 <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata,
                                design = ~ condition)
dds_1 <- estimateSizeFactors(dds_1)

## annotation object
anno_df <- data.frame(
    gene_id = rownames(dds_1),
    gene_name = mapIds(org.Hs.eg.db,
                       keys = rownames(dds_1),
                       column = "SYMBOL",
                       keytype = "ENSEMBL"),
    stringsAsFactors = FALSE,
    row.names = rownames(dds_1)
)

## res object - has to be a DESeq2 results object
#res_de <- readRDS(paste0(myDirectory, "/results/mouse_KMR_TO_res.rds"))

## convert res object to enrich_object
res_de_df <- as.data.frame(res1)
res_de_df$Gene <- rownames(res_de_df)
res_de_df <- res_de_df %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene, "ENTREZID", "ENSEMBL") %>%
               unname())

# Use simulated data set
dataset <- res_de_df %>%
    filter(!is.na(padj)) %>%
    filter(!is.na(ENTREZID)) %>%
    filter(!duplicated(ENTREZID))

padj_thr <- 0.05
log2FC_thr <- 0.05

geneList <- dataset$log2FoldChange
names(geneList)<- dataset$ENTREZID
genes <- dataset$ENTREZID[dataset$padj <= padj_thr &
                              abs(dataset$log2FoldChange) >= log2FC_thr]

res_GO <- enrichGO(gene          = genes,
                   universe      = names(geneList),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP", ### you can change this to another ontology
                   pAdjustMethod = "BH",
                   minGSSize     = 5,
                   maxGSSize     = 500,
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 1,
                   readable      = TRUE)

res_enrich <- shake_enrichResult(res_GO)

saveRDS(dds_1, paste0(myDirectory, "/results/dds_1.rds"))
saveRDS(anno_df, paste0(myDirectory, "/results/anno_df_1.rds"))
saveRDS(res1, paste0(myDirectory, "/results/res_de_1.rds"))
saveRDS(res_enrich, paste0(myDirectory, "/results/res_enrich_1.rds"))
