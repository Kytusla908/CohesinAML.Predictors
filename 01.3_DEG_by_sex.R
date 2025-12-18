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
TCGA_BEAT_df <- t(read.table("InputTables/Input_TCGA-BEAT_raw_counts_common_genes.txt", sep="\t", header=T))

# Get AML samples and set same order for both df
TCGA_BEAT_mut_table <- subset(TCGA_BEAT_mut_table, TCGA_BEAT_mut_table$Disease == "AML")

# Check sample order are correct
table(colnames(TCGA_BEAT_df) == row.names(TCGA_BEAT_mut_table))


# Get sex info =============================
tcga_sex <- read.table("InputTables/TCGA_laml_tcga_pub_clinical_data.tsv", header = T, sep = "\t")
tcga_sex$Sample.ID <- gsub("-", ".", tcga_sex$Sample.ID)
beat_sex <- read_tsv("InputTables/BEAT_aml_ohsu_2022_clinical_data.tsv")
beat_sex$`Sample ID` <- paste0(substr(beat_sex$`Sample ID`, 20, 40), "R")

#Merge them
beat_sex <- beat_sex %>% 
  select(`Sample ID`, Sex) %>% 
  filter(`Sample ID` %in% rownames(TCGA_BEAT_mut_table))
colnames(beat_sex) <- c("sample_ID", "Sex")

tcga_sex <- tcga_sex %>%
  select(Sample.ID, Sex) %>% 
  filter(Sample.ID %in% rownames(TCGA_BEAT_mut_table))
colnames(tcga_sex) <- c("sample_ID", "Sex")

sex_table <- rbind(tcga_sex, beat_sex)
TCGA_BEAT_mut_table$sample_ID <- rownames(TCGA_BEAT_mut_table)
TCGA_BEAT_mut_table <- merge(TCGA_BEAT_mut_table, sex_table, by = "sample_ID")
row.names(TCGA_BEAT_mut_table) <- TCGA_BEAT_mut_table$sample_ID
TCGA_BEAT_mut_table <- TCGA_BEAT_mut_table[colnames(TCGA_BEAT_df),]


# Remove DEGenes between TCGA and BEAT Databases =======================
dds <- DESeqDataSetFromMatrix(countData = TCGA_BEAT_df,
                              colData = TCGA_BEAT_mut_table,
                              design = ~ Sex)
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
cols <- c("Male" = "indianred3", "Female" = "forestgreen")
plotPCA.DESeqTransform(vsd, intgroup="Sex") +
  # geom_text(aes(label = name), vjust = -0.5, size = 3) +
  theme(panel.grid.major=element_line(colour="white"), panel.grid.minor=element_line(colour="white")) +
  labs(color="Sex:") +
  scale_colour_manual(values = cols) +
  theme(legend.position = "right")
# ggsave("plots/PCA_TCGA_BEAT_vst_raw_by_sex.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)

# Get the DEGenes
res.vol <- as.data.frame(results(dds, contrast=c("Sex", "Male", "Female")))
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
  labs(title="Male vs Female DEGenes",
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
  annotate("text", x = 3, y = 200, label = summary(res.vol$condition)["UP"],
           size = 5, hjust = 0, color = "indianred3", fontface="bold")
# ggsave("plots/Volcano_plot_vst_DEGenes_by_sex.png",device="png", units = "cm",
#        width = 10, height = 10, pointsize = 10, dpi = 500)

DEGene <- res.vol %>% filter(condition == "UP") %>% 
  dplyr::select(Gene_ID) %>% rownames(.)

# write.table(res.vol, "OutputTables/DE_genes_filterByExpr_by_Sex_TCGA_BEAT.txt",
#             col.names = T, row.names = T, sep = "\t", quote = F)

