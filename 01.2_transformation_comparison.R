source("00_helper_functions.R")
library(tidyverse)
library(DESeq2)
library(edgeR)
library(VennDiagram)
library(caret)

# How different are vst and rlog ?? ==================
dds <- readRDS("OutputTables/TCGA-BEAT_raw_DESeqDF_filterByExpr_byDatabase.rds")

# perform rlog
start <- Sys.time()
ntd <- normTransform(dds)
end <- Sys.time()
end - start
saveRDS(ntd, file = "TCGA-BEAT_raw_DESeq2_ntd_filterByExpr_byDatabase.rds")
# write.table(as.data.frame(assay(rld_rlog)), "InputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_rlog.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

assay <- as.data.frame(assay(ntd))









