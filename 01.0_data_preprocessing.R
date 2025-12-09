source("00_helper_functions.R")
library(tidyverse)
library(DESeq2)
library(edgeR)
library(VennDiagram)
library(gridExtra)
library(caret)
  

# Read data ===============
# TCGA-BEAT data
TCGA_BEAT_mut_table <- read.table("InputTables/TCGA_BEAT_ALL_mutation_table.txt", sep="\t", header=T)
TCGA_BEAT_df <- t(read.table("InputTables/TCGA-BEAT_raw_counts.txt", sep="\t", header=T))

# Get AML samples and set same order for both df
TCGA_BEAT_mut_table <- subset(TCGA_BEAT_mut_table, TCGA_BEAT_mut_table$Disease == "AML")
TCGA_BEAT_df <- TCGA_BEAT_df[row.names(TCGA_BEAT_mut_table),] # Select and order samples in the same way

# Read Fischer data
Fischer_df <- read.table("InputTables/Fischer_counts_table.txt", skip=1, header=T)
row.names(Fischer_df) <- Fischer_df$Geneid
Fischer_df <- t(Fischer_df[, 7:dim(Fischer_df)[2]])

# Load Fischer metadata
Fischer_metadata <- read.table("InputTables/Metadata_RNAseq_cohesin_AML.inclQC.txt",
                                 sep = "\t", header = T)
Fischer_metadata$Sample_Name <- sapply(strsplit(Fischer_metadata$Sample_Name, "_"), function(parts) {
  paste(tail(parts, 2), collapse = "_")})
row.names(Fischer_metadata) <- Fischer_metadata$Sample_Name

# Fix Fischer_df sample names
row.names(Fischer_df) <- sapply(strsplit(row.names(Fischer_df), "\\."), function(parts) {
  id <- sub("_.*", "", parts[2])
  hits <- grep(id, Fischer_metadata$Sample_Name)
  return(Fischer_metadata$Sample_Name[hits])})

# Set Fischer_df to the metadata order
Fischer_df <- Fischer_df[Fischer_metadata$Sample_Name,]
# write.table(t(Fischer_df), "InputTables/Fischer_counts_table_raw.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Get labels for the samples ======================================
# Set labels for TCGA-BEAT data
TCGA_BEAT_labels_df <- data.frame(sample=row.names(TCGA_BEAT_mut_table),
                                  label=c(ifelse(TCGA_BEAT_mut_table$Cohesin=="Mutant", "cohesinAML", "wtAML")))
table(TCGA_BEAT_labels_df$label)
# write.table(TCGA_BEAT_labels_df, "InputTables/TCGA_BEAT_labels.txt",
#             row.names = F, quote = F, sep = "\t")

# Set labels for Fischer data
Fischer_labels_df <- data.frame(sample=row.names(Fischer_metadata),
                                label=c(ifelse(Fischer_metadata$cohesin_mut_type=="", "wtAML", "cohesinAML")))
table(Fischer_labels_df$label)
# write.table(Fischer_labels_df, "InputTables/Fischer_labels.txt",
#             row.names = F, quote = F, sep = "\t")


# Select Genes in all datasets ===========================
common_genes_table <- as.data.frame(table(colnames(TCGA_BEAT_df) %in% colnames(Fischer_df)))
colnames(common_genes_table) <- c("Common", "Count")
ggplot(common_genes_table, aes(x = "", y = Count, fill = Common)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black") +
  coord_polar(theta = "y") +
  labs(title="Common genes in both TCGA-BEAT and Fischer datasets") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
# ggsave("plots/common_genes_in_datasets.png", device="png", units="cm", dpi=500,
#        width=15, height=12)

keep <- colnames(TCGA_BEAT_df) %in% colnames(Fischer_df)
TCGA_BEAT_df <- TCGA_BEAT_df[, keep]


# Plot mutants proportions =====================================
## TCGA-BEAT
cohesin_tab <- as.data.frame(table(TCGA_BEAT_mut_table$Cohesin))
colnames(cohesin_tab) <- c("Cohesin", "Count")
ggplot(cohesin_tab, aes(x = "", y = Count, fill = Cohesin)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black") +
  coord_polar(theta = "y") +
  labs(title="Cohesin mutant and WT AML in TCGA-BEAT samples") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
# ggsave("plots/Cohesin_mut_proportions_TCGA_BEAT.png", device="png",
#        units="cm", dpi=500, width=15, height=12)

## Fischer
Fischer_coh_mut <- ifelse(Fischer_metadata$cohesin_mut_type=="", "WT", "Mutant")
cohesin_tab <- as.data.frame(table(Fischer_coh_mut))
colnames(cohesin_tab) <- c("Cohesin", "Count")
ggplot(cohesin_tab, aes(x = "", y = Count, fill = Cohesin)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black") +
  coord_polar(theta = "y") +
  labs(title="Cohesin mutant and WT AML in Fischer samples") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
# ggsave("plots/Cohesin_mut_proportions_Fischer.png", device="png",
#        units="cm", dpi=500, width=15, height=12)


# Remove DEGenes between TCGA and BEAT Databases =======================
dds <- DESeqDataSetFromMatrix(countData = t(TCGA_BEAT_df),
                              colData = TCGA_BEAT_mut_table,
                              design = ~ Database)
dds <- DESeq(dds)

# Remove low expressed genes, check after and before removing
dim(dds)
keep <-filterByExpr(dds)
dds <- dds[keep,]
dim(dds)
# saveRDS(dds, file = "OutputTables/TCGA-BEAT_raw_DESeqDF_filterByExpr_byDatabase.rds")

dds <- readRDS("OutputTables/TCGA-BEAT_raw_DESeqDF_filterByExpr_byDatabase.rds")

# Normalize
vsd <- vst(dds)
# saveRDS(vsd, file = "OutputTables/TCGA-BEAT_raw_DESeq2vst_filterByExpr_byDatabase.rds")
# write.table(as.data.frame(assay(vsd)), "InputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst.txt", 
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Check PCA
cols <- c("TCGA" = "royalblue", "BEAT" = "gold")
plotPCA.DESeqTransform(vsd, intgroup="Database") +
  # geom_text(aes(label = name), vjust = -0.5, size = 3) +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Database:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_vst_raw.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)
# ggsave("plots/PCA_TCGA_BEAT_vst_raw_with_labels.png", device = "png", width = 40, height = 27,
#        units = "cm", pointsize = 10, dpi = 500)

# Get the DEGenes
res.vol <- as.data.frame(results(dds, contrast=c("Database", "TCGA", "BEAT")))
res.vol$Gene_ID <- rownames(res.vol)
res.vol <- res.vol %>% 
  mutate(condition = if_else(log2FoldChange>2 & padj < 0.0005, 'UP',
                             if_else(log2FoldChange < -2 & padj < 0.0005, 'DOWN', 'UNCHANGED')))
summary(is.na(res.vol))
res.vol$condition <- as.factor(res.vol$condition)
summary(res.vol$condition)
# write.table(res.vol, "OutputTables/DE_genes_filterByExpr_by_Database_TCGA_BEAT.txt",
#             col.names = T, row.names = T, sep = "\t", quote = F)

# Volcano plot of the DEGenes
cols <- c("UP" = "indianred3", "DOWN" = "royalblue3", "UNCHANGED" = "grey80")
ggplot(res.vol, aes(y = -log10(padj), x = log2FoldChange, col = as.factor(condition))) + 
  geom_point(size = 0.5) + 
  labs(title="TCGA vs BEAT DEGenes",
       y = "-log10(Padj)", x = "log2FC", color = "Gene classification") +
  scale_colour_manual(values = cols) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, size=16, face="bold"), 
        axis.text = element_text(size=14, color = "black"),
        axis.title = element_text(size=14, face="bold", color = "black"), 
        legend.position = "none") + 
  geom_hline(yintercept=0, linetype="dashed") + 
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_hline(yintercept= -log10(0.0005), linetype="dotted", color = "grey40") +
  geom_vline(xintercept= c(-2,2), linetype="dotted", color = "grey40") +
  annotate("text", x = -6, y = 200, label = summary(res.vol$condition)[1],
           size = 5, hjust = 0, color = "royalblue3", fontface="bold") +
  annotate("text", x = 5, y = 200, label = summary(res.vol$condition)[3],
           size = 5, hjust = 0, color = "indianred3", fontface="bold")
# ggsave("plots/Volcano_plot_vst_DEGenes_by_DATABASE.png",device="png", units = "cm",
#        width = 10, height = 10, pointsize = 10, dpi = 500)

# Filter out these DEGenes
non_DEGenes <- res.vol %>% 
  filter(condition == "UNCHANGED") %>% 
  dplyr::select(Gene_ID) %>% 
  rownames(.)
keep <- rownames(vsd) %in% non_DEGenes
vsd <- vsd[keep,]

# Save filterByExpr dataset
Normalized_counts <- as.data.frame(assay(vsd))
# write.table(Normalized_counts, "OutputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst_DEGfiltered.txt", 
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Also save set of genes to be used as input
# write.table(row.names(Normalized_counts), "OutputTables/Genes_input_filterByExpr_DEGfiltered.txt",
#             row.names = F, col.names = F, quote = F, sep = "\t")

cols <- c("TCGA" = "royalblue", "BEAT" = "gold")
plotPCA.DESeqTransform(vsd, intgroup="Database") +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Database:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_vst_removed_DEG.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)


# How different are filterByExpr and rowSums(counts(dds)) >= 10 ?? ==================

## FilterByExpr is done above

## rowSums(counts(dds)) >= 10
dds_rowSum <- DESeqDataSetFromMatrix(countData = t(TCGA_BEAT_df),
                              colData = TCGA_BEAT_mut_table,
                              design = ~ Database)
dds_rowSum <- DESeq(dds_rowSum)
keep <- rowSums(counts(dds_rowSum)) >= 10
dds_rowSum <- dds_rowSum[keep,]
res.vol_rowSum <- as.data.frame(results(dds_rowSum, contrast=c("Database", "TCGA", "BEAT")))
res.vol_rowSum$Gene_ID <- rownames(res.vol_rowSum)
res.vol_rowSum <- res.vol_rowSum %>% 
  mutate(condition = if_else(log2FoldChange>2 & padj < 0.0005, 'UP',
                             if_else(log2FoldChange < -2 & padj < 0.0005, 'DOWN', 'UNCHANGED')))
summary(is.na(res.vol_rowSum))
res.vol_rowSum$condition <- as.factor(res.vol_rowSum$condition)

## Check similarities
input <- list(raw_data=colnames(TCGA_BEAT_df), filterByExpr=res.vol$Gene_ID, rowSum=res.vol_rowSum$Gene_ID)
venn.diagram(input, filename ='plots/VennDiagram_genes_filterByExpr_vs_rowSum.png', height = 2500,
             width = 2500, resolution = 500, imagetype = "tiff", units = "px", 
             compression ="lzw", na = "stop",
             main = "Comparison of methods to filter genes out by expression", 
             main.fontface = "plain", main.fontfamily = "sans",
             fontfamily = "sans", main.col = "black", main.pos = c(0.42, 1.1),
             main.cex = 1, main.just = c(0.5, 1), category.names = c("raw_data", "filterByExpr", "rowSum"),
             cat.pos = c(0, 0, 0), cat.dist = c(0.03, -0.05, -0.02), cat.cex = 1, 
             force.unique =TRUE, print.mode = "raw", sigdigs = 2, direct.area = F,
             hyper.test = FALSE, total.population = NULL, 
             lower.tail = TRUE,  alpha = 0.5, lty="blank",
             fill = c("darkgreen", "seagreen", "lightgreen"), cex = 0.5)
file.remove(list.files(path = "./plots", pattern = ".log$", full.names = TRUE))

genes_df <- data.frame(data_set=c("raw_data", "rowSum", "filterByExpr"),
                       genes=c(dim(TCGA_BEAT_df)[2], length(res.vol_rowSum$Gene_ID), length(res.vol$Gene_ID)))
rownames(genes_df) <- c(""," ","  ")
print(genes_df)
grid.table(genes_df)


# Feature Selection ==========================================
# Correlation based 
genes_keep <- scan("OutputTables/Genes_input_filterByExpr_DEGfiltered.txt", what = "", sep = "\n")
TCGA_BEAT_DEGfiltered <- TCGA_BEAT_df[, genes_keep]
dds <- DESeqDataSetFromMatrix(countData = t(TCGA_BEAT_DEGfiltered),
                              colData = TCGA_BEAT_mut_table,
                              design = ~ 1)
dds <- DESeq(dds)
vsd <- vst(dds, blind=T)
# saveRDS(vsd, file = "OutputTables/Input_TCGA-BEAT_raw_DESeq2_vst_filterByExpr_by1.rds")
Normalized_counts_vst <- as.data.frame(assay(vsd))
# write.table(Normalized_counts_vst, "InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_DEGfiltered.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Correlation based Feature Selection
cor_matrix_vst <- cor(t(Normalized_counts_vst))
to_drop_vst <- findCorrelation(cor_matrix_vst, cutoff = 0.8, exact  = TRUE,
                               names  = TRUE, verbose=TRUE)
# write.table(to_drop_vst, "OutputTables/Genes_to_drop_by_correlation_filterByExpr_vst.txt",
#             row.names = F, col.names = F, sep = "\t")

# Eliminating highly correlated variables
Normalized_counts_vst_corr_filtered <- Normalized_counts_vst[!(row.names(Normalized_counts_vst) %in% to_drop_vst),]
# write.table(Normalized_counts_vst_corr_filtered, "InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Evaluate
cat("No of variables after setting 0.8 as correlation cut-off: ", ncol(Normalized_counts_corr_filtered),
    "\nEliminated variables:\n", to_drop_vst)


# Variance based Feature Selection
variances <- data.frame(t(apply(Normalized_counts_vst_corr_filtered, 1, var)))
table(variances > 0)




