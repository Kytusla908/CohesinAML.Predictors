source("helper_functions.R")
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


# Plot mutants proportions =====================================
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
# ggsave("plots/proportions_TCGA_BEAT.png", device="png", units="cm", dpi=500,
#        width=15, height=12)

# Set labels
TCGA_BEAT_labels_df <- data.frame(sample=row.names(TCGA_BEAT_mut_table),
                                  label=c(ifelse(TCGA_BEAT_mut_table$Cohesin=="Mutant", "cohesinAML", "wtAML")))
table(TCGA_BEAT_labels_df$label)
# write.table(TCGA_BEAT_labels_df, "InputTables/TCGA_BEAT_labels.txt",
#             row.names = F, quote = F, sep = "\t")


# Change Fischer sample names =====================================
sample_guide <- read.csv('InputTables/SampleGuide.csv', header=F)
sample_guide.names <- sapply(strsplit(sample_guide$V2, "_"), function(parts) {
  paste(tail(parts, 2), collapse = "_")
})

Fischer.row.names <- sapply(strsplit(row.names(Fischer_df), "\\."), function(parts) {
  parts[2]})
Fischer.row.names <- sub("_.*", "", Fischer.row.names)

newcol.names <- sapply(Fischer.row.names, function(name) {
  match <- sample_guide.names[grep(name, sample_guide.names)]
  if (length(match) > 0) {
    match[1]
  }
})
row.names(Fischer_df) <- newcol.names


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


# Remove DEGenes between Databases =======================
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
rld <- vst(dds)

# Check PCA
cols <- c("TCGA" = "royalblue", "BEAT" = "gold")
plotPCA.DESeqTransform(rld, intgroup="Database") +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Database:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_raw.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)

# Get the DEGenes
res.vol<-as.data.frame(results(dds, contrast=c("Database", "TCGA", "BEAT")))
res.vol$Gene_ID<-rownames(res.vol)
res.vol <- res.vol %>% 
  mutate(condition = if_else(log2FoldChange>2 & padj < 0.0005, 'UP',
                             if_else(log2FoldChange < -2 & padj < 0.0005, 'DOWN', 'UNCHANGED')))
summary(is.na(res.vol))
res.vol$condition <- as.factor(res.vol$condition)
summary(res.vol$condition)
# write.table(res.vol, "OutputTables/DE_genes_filterByExpr_by_Database_TCGA_BEAT.txt", col.names = T, row.names = T, sep = "\t", quote = F)

# Volcano plot of the DEGenes
cols <- c("UP" = "indianred3", "DOWN" = "royalblue3", "UNCHANGED" = "grey80")
ggplot(res.vol, aes(y = -log10(padj), x = log2FoldChange, col = as.factor(condition))) + 
  geom_point(size = 0.5) + 
  ggtitle("BEAT vs TCGA DEGenes") + 
  labs(y = "-log10(Padj)", x = "log2FC", color = "Gene classification") +
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
# ggsave("plots/Volcano_plot_DEGenes_by_DATABASE.png",device="png", units = "cm",
#        width = 10, height = 10, pointsize = 10, dpi = 500)

# Filter out these DEGenes
non_DEGenes <- res.vol %>% 
  filter(condition == "UNCHANGED") %>% 
  dplyr::select(Gene_ID) %>% 
  rownames(.)
keep <- rownames(rld) %in% non_DEGenes
rld <- rld[keep,]

cols <- c("TCGA" = "royalblue", "BEAT" = "gold")
plotPCA.DESeqTransform(rld, intgroup="Database") +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Database:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_removed_DEG.png", device = "png", width = 12, height = 12,
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
common_genes <- intersect(colnames(TCGA_BEAT_df), intersect(res.vol$Gene_ID, res.vol_rowSum$Gene_ID))
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

## Save filterByExpr dataset
Normalized_counts <- as.data.frame(assay(rld))
# write.table(Normalized_counts, "OutputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst.txt", 
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Correlation based Feature Selection ===================
Normalized_counts <- t(read.table("OutputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst.txt",
                                sep = "\t"))
cor_matrix <- cor(Normalized_counts)
to_drop <- findCorrelation(cor_matrix, cutoff = 0.8, exact  = TRUE,
                           names  = TRUE, verbose=TRUE)
# write.table(to_drop, "OutputTables/Genes_to_drop_by_correlation.txt",
#             row.names = F, col.names = F, sep = "\t")

# Eliminating highly correlated variables
Normalized_counts_corr_filtered <- Normalized_counts[, !(colnames(Normalized_counts) %in% to_drop)]
# write.table(t(Normalized_counts_corr_filtered), "OutputTables/NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")

# Evaluate
cat("No of variables after setting 0.8 as correlation cut-off: ", ncol(Normalized_counts_corr_filtered),
    "\nEliminated variables:\n", to_drop)


# Variance based Feature Selection ==========
variances <- data.frame(t(apply(Normalized_counts_corr_filtered, 2, var)))
table(variances > 0)

# Check presence in Fischer =============
table(colnames(Normalized_counts_corr_filtered) %in% colnames(Fischer_df))












