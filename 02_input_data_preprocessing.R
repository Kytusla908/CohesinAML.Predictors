library(tidyverse)
library(DESeq2)
library(caret)


# Load all the data ============================================================
# TCGA-BEAT data
TCGA_BEAT_mut_table <- read.table("InputTables/TCGA_BEAT_ALL_mutation_table.txt", sep="\t", header=T)
TCGA_BEAT_df <- t(read.table("InputTables/TCGA-BEAT_raw_counts.txt", sep="\t", header=T))

# Fischer data
Fischer_df <- read.table("InputTables/Fischer_counts_table_raw.txt", header=T)
Fischer_metadata <- read.table("InputTables/Metadata_RNAseq_cohesin_AML.inclQC.txt",
                               sep = "\t", header = T)

# Get AML samples and set same order for both df
TCGA_BEAT_mut_table <- subset(TCGA_BEAT_mut_table, TCGA_BEAT_mut_table$Disease == "AML")
TCGA_BEAT_df <- TCGA_BEAT_df[row.names(TCGA_BEAT_mut_table),] # Select and order samples in the same way

# Fix Fischer metadata
Fischer_metadata$Sample_Name <- sapply(strsplit(Fischer_metadata$Sample_Name, "_"), function(parts) {
  paste(tail(parts, 2), collapse = "_")})
row.names(Fischer_metadata) <- Fischer_metadata$Sample_Name





