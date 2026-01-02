source("00_helper_functions.R")
library(tidyverse)
library(sva)
library(caret)
library(gbm)


# Load data =================================
# X data
vst_min_max <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_min_max.txt", header=T)
Fischer_vst_min_max <- read.table("InputTables/Input_Fischer_top3000_vst_min_max.txt", header=T)
stopifnot(identical(colnames(Fischer_vst_min_max), colnames(vst_min_max)))

# y data 
# Load labels 
TCGA_BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels_cohesin.txt", header=T)
TCGA_BEAT_labels <- c(TCGA_BEAT_labels$label)
Fischer_labels <- read.table("InputTables/Fischer_labels_cohesin.txt", header=T)
Fischer_labels <- c(Fischer_labels$label)
Fischer_labels <- as.factor(Fischer_labels)

# Load partition indexes
train_index <- scan("InputTables/Input_train_indexes.txt", sep="\n")

# Splitting data
y_train <- as.factor(TCGA_BEAT_labels[train_index])
y_test <- as.factor(TCGA_BEAT_labels[-train_index])
x_train_vst_min_max <- vst_min_max[train_index, ]
x_test_vst_min_max <- vst_min_max[-train_index, ]


# Set random seed =================================
set.seed(12345)


# Read best tuned gbm model =======================
gbm_tuned <- readRDS("OutputTables/cohesin_tuned_gbm_vst_min_max.rds")
performance <- get_performance(gbm_tuned, Fischer_vst_min_max, Fischer_labels,
                               classes = c("cohesinAML", "wtAML"))

# Optimized threshold for Fischer data ===========================
coordinates <- coords(performance$roc_obj, "best", best.method = "youden",
                      ret = c("threshold", "sensitivity", "specificity"))
threshold <- as.numeric(coordinates["threshold"])
optimal_preds <- factor(ifelse(performance$predicted_probabilities > threshold,
                               "cohesinAML", "wtAML"),
                        levels = levels(Fischer_labels))
confusionMatrix(optimal_preds, Fischer_labels, positive = "cohesinAML")
write_confusion_matrix(confusionMatrix(optimal_preds, Fischer_labels, positive = "cohesinAML"),
                       file="OutputTables/cohesin_validation_Fischer_tuned_gbm_opt_threshold_cm.txt")


# Remove batch effects ==============================
# Load original data
TCGA_BEAT_vst <- read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_varFiltered.txt")
row.names(TCGA_BEAT_vst) <- gsub("-", ".", row.names(TCGA_BEAT_vst))
Fischer_vst <- t(read.table("InputTables/Input_NormalizedCounts_Fischer_top3000_vst.txt"))

# Get matrix
tcga_matrix <- as.matrix(TCGA_BEAT_vst)
fischer_matrix <- as.matrix(Fischer_vst)
dim(tcga_matrix)
dim(fischer_matrix)

# Combine data and set batch labels
combined_matrix <- cbind(tcga_matrix, fischer_matrix)
dim(combined_matrix)
batch <- c(rep("TCGA-BEAT", ncol(tcga_matrix)),
           rep("Fischer", ncol(fischer_matrix)))


# Run ComBat
combat_matrix <- ComBat(dat = combined_matrix, batch = batch,
                        par.prior = TRUE, prior.plots = FALSE)
dim(combat_matrix)

# Separate the matrix
tcga_combat <- t(combat_matrix[, batch == "TCGA-BEAT"])
fischer_combat <- t(combat_matrix[, batch == "Fischer"])
dim(tcga_combat)
dim(fischer_combat)

# Apply min-max scaling
min <- apply(tcga_combat, 2, min)
max <- apply(tcga_combat, 2, max)

scale_with_train <- function(matrix, min, max) {
  sweep(sweep(matrix, 2, min, "-"), 2, max - min, "/")
}

tcga_combat_min_max <- scale_with_train(tcga_combat, min, max)
fischer_combat_min_max <- scale_with_train(fischer_combat, min, max)

# Check distribution
TCGA_BEAT_values <- as.vector(as.matrix(tcga_combat_min_max))
FISCHER_values <- as.vector(as.matrix(fischer_combat_min_max))
distribution_data <- data.frame(value = c(TCGA_BEAT_values, FISCHER_values),
                                dataset = rep(c("TCGA-BEAT", "FISCHER"), 
                                              times = c(length(TCGA_BEAT_values), length(FISCHER_values))))
ggplot(distribution_data, aes(x = value, color = dataset)) +
  geom_density(size = 1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA)) +
  labs(title = "Global distribution of batch corrected min-max values",
       x = "Scaled value", y = "Density")

# ggsave("plots/Global_distribution_TCGA_BEAT_Fischer_vst_min_max_top3000_ComBat.png", device = "png",
#        width = 13, height = 12, units = "cm", pointsize = 10, dpi = 500)
# write.table(fischer_combat_min_max, "InputTables/Input_Fischer_top3000_vst_min_max_ComBat.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")


# Check performance after batch adjusting =====================
fischer_combat_min_max <- read.table("InputTables/Input_Fischer_top3000_vst_min_max_ComBat.txt")
# Best model
performance <- get_performance(gbm_tuned, fischer_combat_min_max, Fischer_labels,
                               classes = c("cohesinAML", "wtAML"),
                               plot_title = "ROC curve of gbm with top Sensitivity (on ComBat Fischer data)",
                               plot_subtitle=paste0("n.trees=", gbm_tuned[["bestTune"]][["n.trees"]],
                                                    " interaction.depth=", gbm_tuned[["bestTune"]][["interaction.depth"]],
                                                    " shrinkage=", gbm_tuned[["bestTune"]][["shrinkage"]],
                                                    " n.minobsinnode=", gbm_tuned[["bestTune"]][["n.minobsinnode"]]))
# ggsave("plots/ROC_Fischer_ComBat_tuned_gbm_vst_min_max_3000.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)
# write_confusion_matrix(performance$confusion_matrix,
#                        "OutputTables/cohesin_validation_Fischer_ComBat_tuned_gbm_cm.txt")

# Top 10 best performing gbms
top10_gbm <- readRDS("OutputTables/cohesin_tuned_top10_gbm_models.rds")
fischer_top5_roc <- vector("list", 5)
fischer_top5_auc <- c()
for (i in 1:5) {
  model <- top10_gbm[[i]]
  performance_df <- get_performance(top10_gbm[[i]], fischer_combat_min_max, Fischer_labels,
                                    classes = c("cohesinAML", "wtAML"))
  fischer_top5_auc <- append(fischer_top5_auc, performance_df$auc)
  fischer_top5_roc[[i]] <- performance_df$roc_df
}

roc_5_df <- bind_rows(lapply(1:5, function(i) {
  df <- fischer_top5_roc[[i]]
  df$model <- paste0("Model ", i)
  df
}))
auc_labels <- paste0("Model ", 1:5,
                     " (AUC = ",
                     round(fischer_top5_auc[1:5], 3), ")")
roc_5_df$model <- factor(roc_5_df$model,
                         levels = paste0("Model ", 1:5),
                         labels = auc_labels)
ggplot(roc_5_df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves for Top 5 gbm models on ComBat Fischer data",
       x = "False Positive Rate", y = "True Positive Rate",
       color = "Model") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
# ggsave("plots/ROC_Fischer_ComBat_top5_gbm_vst_min_max_3000.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)


# Variable Importance ====================
varImp_df <- data.frame(variable=row.names(varImp(gbm_tuned)$importance),
                        importance=varImp(gbm_tuned)$importance)
varImp_top20 <- varImp_df %>%
  arrange(desc(Overall)) %>%
  head(20)
ggplot(varImp_top20, aes(x = reorder(variable, Overall), y = Overall)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Variable", y = "Importance",
       title = "Top 20 most important genes for cohesin mutation predicition (tuned gbm)")
# ggsave("plots/VariableImportance_plot_gbm_tuned.png", device = "png",
#        width = 18, height = 12, units = "cm", pointsize = 10, dpi = 500)

varImp100 <- varImp_df %>% arrange(desc(Overall)) %>% head(100) %>% select(variable)
# write.table(varImp100, "OutputTables/cohesin_varImp100.txt", quote = FALSE,
#             row.names = FALSE, col.names = FALSE)


# Which samples are wrongly predicted? =========================
Fischer_metadata <- read.table("InputTables/Metadata_RNAseq_cohesin_AML.inclQC.txt", sep = "\t", header = T)
wrong_samples <- which(Fischer_labels == "cohesinAML" & Fischer_labels != performance$predicted_classes)
Fischer_metadata_wrong <- Fischer_metadata[wrong_samples,]
common_value <- sapply(Fischer_metadata_wrong, function(x) length(unique(x)) == 1)
Fischer_metadata_wrong[,common_value == "TRUE"]











