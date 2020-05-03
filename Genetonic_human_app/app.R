## genetonic app with human data

# BiocManager::install("GeneTonic")

library("GeneTonic")
library("DESeq2")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("SummarizedExperiment")
library("clusterProfiler")
library("GO.db")

# define your directory
myDirectory <- "~/tmp/sarcoma_shiny"

## load data
dds_1 <- readRDS(paste0(myDirectory, "/results/dds_1.rds"))
anno_df <- readRDS(paste0(myDirectory, "/results/anno_df_1.rds"))
res_de <- readRDS(paste0(myDirectory, "/results/res_de_1.rds"))
res_enrich <- readRDS(paste0(myDirectory, "/results/res_enrich_1.rds"))

## this part runs the app
GeneTonic(dds = dds_1, ## dds object (SummarizedExperiment)
          res_de = res_de, ## results object
          res_enrich = res_enrich, ## enriched GO object
          annotation_obj = anno_df, ## annotation object
          project_id = "Human Sarcoma Dashboard")
