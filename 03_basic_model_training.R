library(caret)

# Load data and labels ========================
TCGA.BEAT_data <- read.table("OutputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt",
                             header=T, sep="\t")
TCGA.BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels.txt", header=T, sep="\t")

# Transpose to prepare for input 
TCGA.BEAT_data <- t(TCGA.BEAT_data)



























