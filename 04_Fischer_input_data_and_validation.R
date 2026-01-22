source("00_helper_functions.R")
library(tidyverse)
library(DESeq2)
library(preprocessCore)
library(caret)

# Load Fischer data ============================================================
Fischer_df <- read.table("InputTables/Fischer_counts_table_raw.txt", header=T)
Fischer_metadata <- read.table("InputTables/Metadata_RNAseq_cohesin_AML.inclQC.txt",
                               sep = "\t", header = T)
Fischer_metadata$Sample_Name <- sapply(strsplit(Fischer_metadata$Sample_Name, "_"), function(parts) {
  paste(tail(parts, 2), collapse = "_")})
row.names(Fischer_metadata) <- Fischer_metadata$Sample_Name

# Cohesin labels
Fischer_labels_df <- data.frame(sample=row.names(Fischer_metadata),
                                label=c(ifelse(Fischer_metadata$cohesin_mut_type=="", "wtAML", "cohesinAML")))
table(Fischer_labels_df$label)
Fischer_labels <- Fischer_labels_df$label
Fischer_labels <- as.factor(Fischer_labels)
# write.table(Fischer_labels_df, "InputTables/Fischer_labels_cohesin.txt",
#             row.names = F, quote = F, sep = "\t")

# NPM1 labels 
Fischer_labels_df <- data.frame(sample=row.names(Fischer_metadata),
                                label=c(ifelse(Fischer_metadata$NPM1=="pos", "NPM1.Mutant", "NPM1.WT")))
table(Fischer_labels_df$label)
# write.table(Fischer_labels_df, "InputTables/Fischer_labels_NPM1.txt",
#             row.names = F, quote = F, sep = "\t")

# FLT3 labels
Fischer_labels_df <- data.frame(sample=row.names(Fischer_metadata),
                                label=c(ifelse(Fischer_metadata$FLT3_ITD=="pos", "FLT3.Mutant", 
                                               ifelse(Fischer_metadata$FLT3_TKD=="pos", "FLT3.Mutant", "FLT3.WT"))))
table(Fischer_labels_df$label)
# write.table(Fischer_labels_df, "InputTables/Fischer_labels_FLT3.txt",
#             row.names = F, quote = F, sep = "\t")

# FLT3-ITD labels
Fischer_labels_df <- data.frame(sample=row.names(Fischer_metadata),
                                label=c(ifelse(Fischer_metadata$FLT3_ITD=="pos", "FLT3_ITD.Mutant", "FLT3_ITD.WT")))
table(Fischer_labels_df$label)
# write.table(Fischer_labels_df, "InputTables/Fischer_labels_FLT3_ITD.txt",
#             row.names = F, quote = F, sep = "\t")

# FLT3-TKD labels
Fischer_labels_df <- data.frame(sample=row.names(Fischer_metadata),
                                label=c(ifelse(Fischer_metadata$FLT3_TKD=="pos", "FLT3_TKD.Mutant","FLT3_TKD.WT")))
table(Fischer_labels_df$label)
# write.table(Fischer_labels_df, "InputTables/Fischer_labels_FLT3_TKD.txt",
#             row.names = F, quote = F, sep = "\t")


# Load TCGA-BEAT transformation ==============================
TCGA_BEAT_vst_min_max <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_min_max.txt", check.names = FALSE)
TCGA_BEAT_vst_zScore <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_zScore.txt", check.names = FALSE)
TCGA_BEAT_vst_QN <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_QN.txt", check.names = FALSE)
TCGA_BEAT_rlog <- read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_rlog_varFiltered.txt", check.names = FALSE)
TCGA_BEAT_rlog_min_max <- read.table("InputTables/Input_TCGA-BEAT_top3000_rlog_min_max.txt", check.names = FALSE)
TCGA_BEAT_rlog_zScore <- read.table("InputTables/Input_TCGA-BEAT_top3000_rlog_zScore.txt", check.names = FALSE)
TCGA_BEAT_rlog_QN <- read.table("InputTables/Input_TCGA-BEAT_top3000_rlog_QN.txt", check.names = FALSE)


# Load best varFiltered CV+ROSE model ============================
gbm_vst_min_max <- readRDS("OutputTables/basic_gbm_vst_min_max_3000_CV_SMOTE_model.rds")


# Select same genes =====================================
# Run DESeq2 
fischer_dds <- DESeq(DESeqDataSetFromMatrix(countData = Fischer_df,
                                            colData = Fischer_metadata,
                                            design = ~ 1))

# Apply vst
fischer_vsd <- as.data.frame(assay(vst(fischer_dds, blind=T)))
fischer_rld <- as.data.frame(assay(rlog(fischer_dds, blind=T)))

# Check genes are the same
table(colnames(TCGA_BEAT_vst_min_max) %in% row.names(fischer_vsd))
table(colnames(gbm_vst_min_max[["trainingData"]]) %in% row.names(fischer_vsd))

table(row.names(TCGA_BEAT_rlog) %in% row.names(fischer_rld))

# Select same genes as training and testing
fischer_vst <- t(fischer_vsd[colnames(TCGA_BEAT_vst_min_max), ])
dim(fischer_vst)
table(colnames(gbm_vst_min_max[["trainingData"]]) %in% colnames(fischer_vst))

fischer_rlog <- t(fischer_rld[row.names(TCGA_BEAT_rlog), ])
dim(fischer_rlog)

# Save Fischer Input data
# write.table(fischer_vst, "InputTables/Input_NormalizedCounts_Fischer_top3000_vst.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")
# write.table(fischer_rlog, "InputTables/Input_NormalizedCounts_Fischer_top3000_rlog.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Perform min-max transformation (vst) ======================
fischer_vst_min_max <- apply(fischer_vst, 2, min_max)
table(is.na(fischer_vst_min_max))

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(TCGA_BEAT_vst_min_max))
FISCHER_values <- as.vector(as.matrix(fischer_vst_min_max))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of min-max scaled values",
       x = "Scaled value", y = "Density")
# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_vst_min_max_top3000.png", device = "png",
#        width = 12, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_vst_min_max, "InputTables/Input_Fischer_top3000_vst_min_max.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Perform zScore transformation (vst) ======================
fischer_vst_zScore <- scale(fischer_vst)
table(is.na(fischer_vst_zScore))
fischer_vst_zScore[, apply(fischer_vst, 2, sd, na.rm = TRUE) == 0] <- 0
table(is.na(fischer_vst_zScore))

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(TCGA_BEAT_vst_zScore))
FISCHER_values <- as.vector(as.matrix(fischer_vst_zScore))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of zScored values",
       x = "Scaled value", y = "Density")
# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_vst_zScore_top3000.png", device = "png",
#        width = 12, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_vst_zScore, "InputTables/Input_Fischer_top3000_vst_zScore.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Perform QN transformation (vst)======================
fischer_vst_QN <- data.frame(t(normalize.quantiles(t(fischer_vst))))
table(is.na(fischer_vst_QN))
rownames(fischer_vst_QN) <- rownames(fischer_vst)
colnames(fischer_vst_QN) <- colnames(fischer_vst)

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(TCGA_BEAT_vst_QN))
FISCHER_values <- as.vector(as.matrix(fischer_vst_QN))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of QN values",
       x = "Scaled value", y = "Density")

# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_vst_QN_top3000.png", device = "png",
#        width = 12, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_vst_QN, "InputTables/Input_Fischer_top3000_vst_QN.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Perform min-max transformation (rlog) ======================
fischer_rlog_min_max <- apply(fischer_rlog, 2, min_max)
table(is.na(fischer_rlog_min_max))

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(TCGA_BEAT_rlog_min_max))
FISCHER_values <- as.vector(as.matrix(fischer_rlog_min_max))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of min-max scaled values",
       x = "Scaled value", y = "Density")
# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_rlog_min_max_top3000.png", device = "png",
#        width = 12, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_rlog_min_max, "InputTables/Input_Fischer_top3000_rlog_min_max.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Perform zScore transformation (rlog) ======================
fischer_rlog_zScore <- scale(fischer_rlog)
table(is.na(fischer_rlog_zScore))
fischer_rlog_zScore[, apply(fischer_rlog, 2, sd, na.rm = TRUE) == 0] <- 0
table(is.na(fischer_rlog_zScore))

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(TCGA_BEAT_rlog_zScore))
FISCHER_values <- as.vector(as.matrix(fischer_rlog_zScore))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of zScored values",
       x = "Scaled value", y = "Density") +
  coord_cartesian(xlim = c(-5, 5))
# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_rlog_zScore_top3000.png", device = "png",
#        width = 12, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_rlog_zScore, "InputTables/Input_Fischer_top3000_rlog_zScore.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Perform QN transformation (rlog)======================
fischer_rlog_QN <- data.frame(t(normalize.quantiles(t(fischer_rlog))))
table(is.na(fischer_rlog_QN))
rownames(fischer_rlog_QN) <- rownames(fischer_rlog)
colnames(fischer_rlog_QN) <- colnames(fischer_rlog)

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(TCGA_BEAT_rlog_QN))
FISCHER_values <- as.vector(as.matrix(fischer_rlog_QN))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of QN values",
       x = "Scaled value", y = "Density")

# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_rlog_QN_top3000.png", device = "png",
#        width = 12, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_rlog_QN, "InputTables/Input_Fischer_top3000_rlog_QN.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Check performance of the basic model using Fischer data (validation) ============
get_performance(gbm_vst_min_max, fischer_vst_min_max, Fischer_labels,
                classes = c("cohesinAML", "wtAML"),
                plot_title = "ROC curve: GBM on Fischer vst_min_max+CV+SMOTE")
# ggsave("plots/ROC_Fischer_basic_gbm_vst_min_max_3000_CV_SMOTE.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)

