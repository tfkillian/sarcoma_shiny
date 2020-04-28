# conversions to get correct results for shiny app

setwd("~/tmp/sarcoma_shiny/scripts/")

library("tidyverse")
library("dplyr")
library("limma")
library("GO.db")
library("org.Hs.eg.db")
library("org.Mm.eg.db")

## save human 'coldata'  
countdata_1 <- c("RMS2021", "RMS2028", "RMS2051", "RMS2052", "RMS2056", "RMS2058",
                 "RMS2059", "RMS2075", "RMS2080", "RMS2094", "RMS2103", "RMS2110",
                 "RMS2114", "RMS216", "RMS231", "NCI0075", "NS127muscle",
                 "NS128muscle", "NS129muscle", "NS132muscle", "NS134muscle")
condition_1 <- factor(c(rep("RMS", 16), rep("ctrl", 5)))
col_1 <- data.frame(row.names = countdata_1, condition_1)

## save mouse 'coldata' 
countdata_2 <- c("T01", "T02", "T03", "T31", "T32", "T33", "T34", "T5A5", "T5A6", "KMR19", "KMR46")
condition_2 <- factor(c(rep("T0", 3), rep("T3", 4), rep("T5", 2), rep("KMR", 2)))
col_2 <- data.frame(row.names = countdata_2, condition_2)

saveRDS(as.data.frame(col_1), file = "../results/col_1.rds")
saveRDS(as.data.frame(col_2), file = "../results/col_2.rds")


### dds results
dds_df_1 <- read.csv("../results/human_results_1.csv", header = TRUE, row.names = 1)
dds_df_2 <- read.csv("../results/mouse_results_1.csv", header = TRUE, row.names = 1)

dds_df_1$Gene <- as.character(dds_df_1$Gene)
dds_df_2$Gene <- as.character(dds_df_2$Gene)

dds_df_1 <- dds_df_1 %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname()) %>%
    dplyr::select(Gene, geneName, ENTREZID, baseMean,   
                  log2FoldChange, lfcSE, stat, pvalue, padj,
                  everything())

dds_df_2 <- dds_df_2 %>%
    mutate(ENTREZID = mapIds(org.Mm.eg.db, Gene,"ENTREZID", "ENSEMBL") %>% unname()) %>%
    dplyr::select(Gene, geneName, ENTREZID, baseMean,   
                  log2FoldChange, lfcSE, stat, pvalue, padj,
                  everything())

saveRDS(as.data.frame(dds_df_1), file = "../results/res_1.rds")
saveRDS(as.data.frame(dds_df_2), file = "../results/res_2.rds")

## we need to save the results with ENTREZ IDs from now on

## get all GO IDs, and all ENTREZ IDs associated with each ID
# go_list <- mapIds(org.Hs.eg.db, keys(org.Hs.eg.db, "GO"),
#                   "ENTREZID", "GO", multiVals = "list")
# go_vector <- lapply(go_list, as.vector)
# ezs <- sapply(go_vector, paste0, collapse = ";")
# go_df <- data.frame(GOID = names(go_vector), ENTREZID = ezs)

## this also misses some genes
# go_list <- as.list(org.Hs.egGO2ALLEGS)
# go_vector <- lapply(go_list, as.vector)
# ezs <- sapply(go_vector, paste0, collapse = ";")
# go_df <- data.frame(GOID = names(go_vector), ENTREZID = ezs)

## okay, so apparently you should use this method from goana
obj <- paste0("org.Hs.egGO2ALLEGS")
orgPkg <- paste0("org.Hs.eg.db")
egGO2ALLEGS <- getFromNamespace(obj,orgPkg)
GeneID_PathID <- AnnotationDbi::toTable(egGO2ALLEGS)
go_df <- GeneID_PathID %>%
         dplyr::select(go_id, gene_id) %>%
         dplyr::rename(GOID = go_id) %>% 
         group_by(GOID) %>%
         unique() %>%
         summarise(ENTREZID_in_term = paste(gene_id, collapse = ";"))

# go_32501 <- go_set %>% dplyr::filter(GOID == "GO:0032501") 
# 
# go_ez <- go_32501 %>%
#          group_by(GOID) %>%
#          dplyr::summarise(ENTREZIDs = paste0(gene_id, collapse = ";"))
# 
# go_list <- go_set %>%
#            group_by(GOID) %>%
#            summarise(ENTREZIDs = paste(gene_id, collapse = ";"))

### then map GO terms to genes between tables via ENTREZID
go_1 <- go_1 %>% dplyr::select(-ENTREZID_in_term)
go_2 <- go_2 %>% dplyr::select(-ENTREZID_in_term)

go_1 <- go_1 %>% dplyr::left_join(go_df, by ="GOID")
go_2 <- go_2 %>% dplyr::left_join(go_df, by ="GOID")

###################### KEGG conversion ###################################
# https://www.researchgate.net/post/How_i_can_get_a_list_of_KEGG_pathways_and_its_list_of_genes
# obj <- paste0("org.Hs.egPATH2EG")
# orgPkg <- paste0("org.Hs.eg.db")
# egPATH2EG <- getFromNamespace(obj,orgPkg)
# GeneID_KEGGID <- AnnotationDbi::toTable(egPATH2EG)

####################### we need to construct exactly what limma goana does ################################
library("limma")
GeneID_PathID <- getGeneKEGGLinks("hsa", convert=FALSE)
# https://rdrr.io/bioc/limma/src/R/kegga.R
gene_list <- GeneID_PathID %>%
    dplyr::select(PathwayID, GeneID) %>%
    group_by(PathwayID) %>%
    unique() %>%
    summarise(ENTREZID_in_path = paste(GeneID, collapse = ";"))
    
PathID_PathName <- getKEGGPathwayNames("hsa", remove.qualifier=TRUE)
path_id_kegg <- dplyr::full_join(PathID_PathName, gene_list, by = "PathwayID")
names(path_id_kegg) <- c("PathwayID", "Pathway", "ENTREZID_in_path")

# library("org.Hs.eg.db")
# mapped <- mappedkeys(org.Hs.egPATH2EG)
# L <- as.list(org.Hs.egPATH2EG[mapped])
# Kegg_ID <- names(L)
# Gene_IDs <- sapply(L, paste, collapse = ";")
# kegg_genes <- cbind(Kegg_ID, Gene_IDs)
# kegg_genes <- as.data.frame(kegg_genes)
# names(kegg_genes) <- c("PathwayID", "ENTREZID_in_path")
# rownames(kegg_genes) <- NULL
# path_id_kegg <- kegg_genes
# write.table(cbind(Kegg_ID, Gene_IDs), file="KEGG to Genes.txt", sep="\t", row.names=FALSE, col.names=FALSE)

### this method seems to provide erroneous ENRTEZIDs 
# kegg_organism <- "mmu"
# GK <- getGeneKEGGLinks(species.KEGG = kegg_organism)
# GK_list <- as.list(GK)
# ezs <- sapply(GK_list, paste0, collapse = ";")

## get all KEGG PathwayIDs, and all ENTREZ IDs associated with each PathwayID
# GK_desc <- getKEGGPathwayNames(species.KEGG = kegg_organism, remove = TRUE)
# GK_col <- GK %>% group_by(PathwayID) %>% summarise(GeneIDs = paste(GeneID, collapse = ";"))
# path_id_kegg <- dplyr::full_join(GK_col, GK_desc, by = "PathwayID")
# names(path_id_kegg) <- c("PathwayID", "ENTREZ_GeneIDs", "Pathway")
# saveRDS(as.data.frame(path_id_kegg), file = "./app/data/kegg_path_2_id_mouse.rds")

# kegg_1$PathwayID <- gsub("path:mmu", "", kegg_1$PathwayID)
# kegg_2$PathwayID <- gsub("path:mmu", "", kegg_2$PathwayID)

kegg_1 <- kegg_1 %>% dplyr::select(-ENTREZID_in_path) %>% dplyr::select(-PathwayID)
kegg_2 <- kegg_2 %>% dplyr::select(-ENTREZID_in_path) %>% dplyr::select(-PathwayID)

kegg_1 <- kegg_1 %>% left_join(path_id_kegg, by = "Pathway")
kegg_2 <- kegg_2 %>% left_join(path_id_kegg, by = "Pathway")

# kegg_1$ENTREZID_in_path <- as.character(kegg_1$ENTREZID_in_path)
# kegg_2$ENTREZID_in_path <- as.character(kegg_2$ENTREZID_in_path)
# kegg_3$ENTREZID_in_path <- as.character(kegg_3$ENTREZID_in_path)
# kegg_4$ENTREZID_in_path <- as.character(kegg_4$ENTREZID_in_path)

# names(kegg_1) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")
# names(kegg_2) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")
# names(kegg_3) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")
# names(kegg_4) <- c("Pathway", "N" , "DE", "P.DE", "PathwayID", "ENTREZID_in_path")

############################# save ######################################

saveRDS(dds_df_1, file = "../results/res1.rds")
saveRDS(dds_df_2, file = "../results/res2.rds")

saveRDS(go_1, file = "../results/go_1.rds")
saveRDS(go_2, file = "../results/go_2.rds")

saveRDS(kegg_1, file = "../results/kegg_1.rds")
saveRDS(kegg_2, file = "../results/kegg_2.rds")
