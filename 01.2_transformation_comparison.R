source("00_helper_functions.R")
library(tidyverse)
library(preprocessCore)
library(ggplot2)
library(cowplot)

# Load all Input data ==================
# Load transformed and/or filtered data
vst_transformed <- t(read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt", header=T))
rlog_transformed <- t(read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_rlog_corrFiltered.txt", header=T))

# Load and filter raw_counts
raw_counts <- read.table("InputTables/Input_TCGA-BEAT_raw_counts_common_genes.txt", sep="\t", header=T)

# check cohesin genes presence
cohesin_genes <- c("STAG2", "STAG1", "SMC3", "SMC1A", "RAD21", "WAPL", "PDS5B", "PDS5A")
cohesin_genes %in% colnames(rlog_transformed)
cohesin_genes %in% colnames(vst_transformed)
"WAPL" %in% colnames(raw_counts)

# Min-Max scaling ==============================================================
vst_min_max <- apply(vst_transformed, 2, min_max)
table(is.na(vst_min_max))

rlog_min_max <- apply(rlog_transformed, 2, min_max)
table(is.na(rlog_min_max))


# Z-score per gene =====================================
vst_zScore <- scale(vst_transformed)
table(is.na(vst_zScore))

rlog_zScore <- scale(rlog_transformed)
table(is.na(rlog_zScore))


# Quantile Normalization ===============================
vst_QN <- data.frame(t(normalize.quantiles(t(vst_transformed))))
rownames(vst_QN) <- rownames(vst_transformed)
colnames(vst_QN) <- colnames(vst_transformed)

rlog_QN <- data.frame(t(normalize.quantiles(t(rlog_transformed))))
rownames(rlog_QN) <- rownames(rlog_transformed)
colnames(rlog_QN) <- colnames(rlog_transformed)


# Let's see the distribution ====================================
df_transformation <- bind_rows(make_long(vst_transformed, "VST"),
                               make_long(rlog_transformed, "rlog"))
df_min_max <- bind_rows(make_long(vst_min_max, "VST Min-Max"),
                        make_long(rlog_min_max, "rlog Min-Max"))
df_zscore <- bind_rows(make_long(vst_zScore, "VST Z-score"),
                       make_long(rlog_zScore, "rlog Z-score"))
df_QN <- bind_rows(make_long(vst_QN, "VST QN"),
                   make_long(rlog_QN, "rlog QN"))

p1 <- ggplot(df_transformation, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +
  coord_cartesian(ylim = c(0,50)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  ylab("Expression values") +
  ggtitle("Global Distribution of vst/rlog transformation")
p2 <- ggplot(df_min_max, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  ylab("Expression values") +
  ggtitle("Global Distribution of vst/rlog + min-max scaling")
p3 <- ggplot(df_zscore, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  ylab("Expression values") +
  ggtitle("Global Distribution of vst/rlog + per gene z-score")
p4 <- ggplot(df_QN, aes(x = Transformation, y = Value, fill = Transformation)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, alpha = 0.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none") +
  ylab("Expression values") +
  ggtitle("Global Distribution of vst/rlog + Quantile Normalization")

plot_grid(p1, p2, p3, p4, ncol = 4)
# ggsave("plots/Global_distribution_transformations.png", device = "png", width = 80, height = 20,
#        units = "cm", pointsize = 10, dpi = 500)


# Save the tables ===================================
write.table(rlog_min_max, "InputTables/Input_TCGA-BEAT_rlog_min_max.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")
write.table(rlog_zScore, "InputTables/Input_TCGA-BEAT_rlog_zScore.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")
write.table(rlog_QN, "InputTables/Input_TCGA-BEAT_rlog_QN.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")

write.table(vst_min_max, "InputTables/Input_TCGA-BEAT_vst_min_max.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")
write.table(vst_zScore, "InputTables/Input_TCGA-BEAT_vst_zScore.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")
write.table(vst_QN, "InputTables/Input_TCGA-BEAT_vst_QN.txt",
            row.names = T, col.names = T, quote = F, sep = "\t")

