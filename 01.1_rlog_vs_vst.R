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
rld_rlog <- rlog(dds)
end <- Sys.time()
end - start
# saveRDS(rld_rlog, file = "TCGA-BEAT_raw_DESeq2rld_rlog_filterByExpr_byDatabase.rds")
# write.table(as.data.frame(assay(rld_rlog)), "InputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_rlog.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Check PCA
cols <- c("TCGA" = "royalblue", "BEAT" = "gold")
plotPCA.DESeqTransform(rld_rlog, intgroup="Database") +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Database:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_rlog_raw.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)

# Read DEGenes
res.vol <- read.table("OutputTables/DE_genes_filterByExpr_by_Database_TCGA_BEAT.txt", header=T, sep = "\t")

# Filter out DEGenes
non_DEGenes <- res.vol %>% 
  filter(condition == "UNCHANGED") %>% 
  dplyr::select(Gene_ID) %>% 
  rownames(.)
keep <- rownames(rld_rlog) %in% non_DEGenes
rld_rlog <- rld_rlog[keep,]

# Save filterByExpr dataset
Normalized_counts_rlog <- as.data.frame(assay(rld_rlog))
# write.table(Normalized_counts_rlog, "InputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_rlog_DEGfiltered.txt", 
#             row.names = T, col.names = T, quote = F, sep = "\t")

cols <- c("TCGA" = "royalblue", "BEAT" = "gold")
plotPCA.DESeqTransform(rld_rlog, intgroup="Database") +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Database:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_rlog_removed_DEG.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)


# rlog Feature Selection =====================================
TCGA_BEAT_mut_table <- read.table("InputTables/TCGA_BEAT_ALL_mutation_table.txt", sep="\t", header=T)
TCGA_BEAT_mut_table <- subset(TCGA_BEAT_mut_table, TCGA_BEAT_mut_table$Disease == "AML")
TCGA_BEAT_df <- read.table("InputTables/TCGA-BEAT_raw_counts.txt", sep="\t", header=T)
TCGA_BEAT_df <- TCGA_BEAT_df[, row.names(TCGA_BEAT_mut_table)]
genes_keep <- scan("OutputTables/Genes_input_filterByExpr_rlog_DEGfiltered.txt", what = "", sep = "\n")

# Correlation based
TCGA_BEAT_DEGfiltered <- TCGA_BEAT_df[genes_keep,]
dds <- DESeqDataSetFromMatrix(countData = TCGA_BEAT_DEGfiltered,
                              colData = TCGA_BEAT_mut_table,
                              design = ~ 1)
dds <- DESeq(dds)
rld_rlog <- rlog(dds, blind=T)
Normalized_counts_rlog <- as.data.frame(assay(rld_rlog))

cor_matrix_rlog <- cor(t(Normalized_counts_rlog))

start <- Sys.time()
to_drop_rlog <- findCorrelation(cor_matrix_rlog, cutoff = 0.8, exact  = TRUE,
                           names  = TRUE, verbose=TRUE)
end <- Sys.time()
end - start
write.table(to_drop_rlog, "OutputTables/Genes_to_drop_by_correlation_filterByExpr_rlog.txt",
            row.names = F, col.names = F, sep = "\t")

## Eliminating highly correlated variables
Normalized_counts_rlog_corr_filtered <- Normalized_counts_rlog[, !(colnames(Normalized_counts_rlog) %in% to_drop_rlog)]
# write.table(t(Normalized_counts_rlog_corr_filtered), "InputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_rlog_corrFiltered.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

## Evaluate
cat("No of variables after setting 0.8 as correlation cut-off: ", ncol(Normalized_counts_rlog_corr_filtered),
    "\nEliminated variables:\n", to_drop_rlog)

# Variance based Feature Selection
variances <- data.frame(t(apply(Normalized_counts_vst_corr_filtered, 2, var)))
table(variances > 0)


# Check similarities/differences =========================
# Load vst Normalized data
Normalized_counts_vst_corr_filtered <- read.table("InputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt",
                                                  header=T, sep="\t")

# Genes left in the dataframes
input <- list(vst_transformation=row.names(Normalized_counts_vst_corr_filtered),
              rlog_transformation=row.names(Normalized_counts_rlog_corr_filtered))
venn.diagram(input, filename ='plots/VennDiagram_vst_vs_rlog_transformations.png', height = 2500,
             width = 2500, resolution = 500, imagetype = "tiff", units = "px", 
             compression ="lzw", na = "stop",
             main = "Genes left in each transformation after removing DDBB DEGenes and highly correlated genes", 
             main.fontface = "plain", main.fontfamily = "sans",
             fontfamily = "sans", main.col = "black", main.pos = c(0.42, 1.1),
             main.cex = 1, main.just = c(0.5, 1), category.names = c("vst_transformation", "rlog_transformation"),
             cat.pos = c(0, 0), cat.dist = c(0.03, -0.02), cat.cex = 1, 
             force.unique =TRUE, print.mode = "raw", sigdigs = 2, direct.area = F,
             hyper.test = FALSE, total.population = NULL, 
             lower.tail = TRUE,  alpha = 0.5, lty="blank",
             fill = c("darkgreen", "lightgreen"), cex = 0.5)
file.remove(list.files(path = "./plots", pattern = ".log$", full.names = TRUE))

# Genes removed by correlation
to_drop_vst <- read.table("OutputTables/Genes_to_drop_by_correlation_filterByExpr_vst.txt", sep="\t", header=F)
input <- list(vst_transformation=to_drop_vst,
              rlog_transformation=to_drop_rlog)
venn.diagram(input, filename ='plots/VennDiagram_vst_vs_rlog_transformations.png', height = 2500,
             width = 2500, resolution = 500, imagetype = "tiff", units = "px", 
             compression ="lzw", na = "stop",
             main = "Genes highly corelated in each transformation", 
             main.fontface = "plain", main.fontfamily = "sans",
             fontfamily = "sans", main.col = "black", main.pos = c(0.42, 1.1),
             main.cex = 1, main.just = c(0.5, 1), category.names = c("vst_transformation", "rlog_transformation"),
             cat.pos = c(0, 0), cat.dist = c(0.03, -0.02), cat.cex = 1, 
             force.unique =TRUE, print.mode = "raw", sigdigs = 2, direct.area = F,
             hyper.test = FALSE, total.population = NULL, 
             lower.tail = TRUE,  alpha = 0.5, lty="blank",
             fill = c("darkgreen", "lightgreen"), cex = 0.5)
file.remove(list.files(path = "./plots", pattern = ".log$", full.names = TRUE))










