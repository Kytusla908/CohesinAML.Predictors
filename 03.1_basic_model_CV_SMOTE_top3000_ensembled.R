source("00_helper_functions.R")
library(tidyverse)
library(preprocessCore)
library(ggplot2)
library(cowplot)
library(caret)
library(pROC)


load("03_CV_SMOTE_3000.RData")
all_models <- readRDS("OutputTables/basic_all_models_CV_SMOTE_3000.rds")


# Data partition  ================================
raw_data <- t(read.table("InputTables/Input_TCGA-BEAT_raw_counts_common_genes.txt", sep="\t", header=T))
vst_transform <- t(Normalized_counts_vst_var_filtered)
rlog_transform <- t(Normalized_counts_rlog_var_filtered)

# Variance select raw data
variances <- apply(raw_data, 1, var)
table(variances > 0)
top_variances <- variances %>%
  enframe(name = "gene", value = "variance") %>%
  filter(variance > 0) %>%
  arrange(desc(variance)) %>%
  slice(1:3000)
raw_data <- t(raw_data[top_variances$gene, ])

# Load labels 
TCGA_BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels_cohesin.txt", header=T)
TCGA_BEAT_labels <- c(TCGA_BEAT_labels$label)

# Load partition indexes
train_index <- scan("InputTables/Input_train_indexes.txt", sep="\n")

# Splitting Labels
y_train <- as.factor(TCGA_BEAT_labels[train_index])
y_test <- as.factor(TCGA_BEAT_labels[-train_index])

# Subset train data
x_train_raw <- raw_data[train_index, ]
x_train_vst <- vst_transform[train_index, ]
x_train_vst_min_max <- vst_min_max[train_index, ]
x_train_vst_zScore <- vst_zScore[train_index, ]
x_train_vst_QN <- vst_QN[train_index, ]
x_train_rlog <- rlog_transform[train_index, ]
x_train_rlog_min_max <- rlog_min_max[train_index, ]
x_train_rlog_zScore <- rlog_zScore[train_index, ]
x_train_rlog_QN <- rlog_QN[train_index, ]

# Subset test data
x_test_raw <- raw_data[-train_index, ]
x_test_vst <- vst_transform[-train_index, ]
x_test_vst_min_max <- vst_min_max[-train_index, ]
x_test_vst_zScore <- vst_zScore[-train_index, ]
x_test_vst_QN <- vst_QN[-train_index, ]
x_test_rlog <- rlog_transform[-train_index, ]
x_test_rlog_min_max <- rlog_min_max[-train_index, ]
x_test_rlog_zScore <- rlog_zScore[-train_index, ]
x_test_rlog_QN <- rlog_QN[-train_index, ]


# Map test_data ================================
test_data_map <- list(
  kNN_vst_min_max = x_test_vst_min_max,
  kNN_vst_zScore = x_test_vst_zScore,
  kNN_rlog_min_max = x_test_rlog_min_max,
  kNN_rlog_zScore = x_test_rlog_zScore,
  
  nb_raw = x_test_raw, nb_vst = x_test_vst, nb_rlog = x_test_rlog,
  nb_vst_min_max = x_test_vst_min_max, nb_vst_zScore = x_test_vst_zScore,
  nb_vst_QN = x_test_vst_QN, nb_rlog_min_max = x_test_rlog_min_max,
  nb_rlog_zScore = x_test_rlog_zScore, nb_rlog_QN = x_test_rlog_QN,
  
  svm_vst_min_max = x_test_vst_min_max, svm_vst_zScore = x_test_vst_zScore,
  svm_rlog_min_max = x_test_rlog_min_max, svm_rlog_zScore = x_test_rlog_zScore,
  
  glm_raw = x_test_raw, glm_vst = x_test_vst, glm_rlog = x_test_rlog,
  glm_vst_min_max = x_test_vst_min_max, glm_vst_zScore = x_test_vst_zScore,
  glm_vst_QN = x_test_vst_QN, glm_rlog_min_max = x_test_rlog_min_max,
  glm_rlog_zScore = x_test_rlog_zScore, glm_rlog_QN = x_test_rlog_QN,
  
  rf_raw = x_test_raw, rf_vst = x_test_vst, rf_rlog = x_test_rlog,
  rf_vst_min_max = x_test_vst_min_max, rf_vst_zScore = x_test_vst_zScore,
  rf_vst_QN = x_test_vst_QN, rf_rlog_min_max = x_test_rlog_min_max,
  rf_rlog_zScore = x_test_rlog_zScore, rf_rlog_QN = x_test_rlog_QN,
  
  gbm_raw = x_test_raw, gbm_vst = x_test_vst, gbm_rlog = x_test_rlog,
  gbm_vst_min_max = x_test_vst_min_max, gbm_vst_zScore = x_test_vst_zScore,
  gbm_vst_QN = x_test_vst_QN, gbm_rlog_min_max = x_test_rlog_min_max,
  gbm_rlog_zScore = x_test_rlog_zScore, gbm_rlog_QN = x_test_rlog_QN
)


# Measure performance ==========================
performance_df <- data.frame(model = character(), sensitivity = numeric(),
                             specificity = numeric(), precision = numeric(),
                             accuracy = numeric(), kappa = numeric(),
                             auc = numeric(), stringsAsFactors = FALSE)
preds_list <- vector("list", length(all_models))
probs_list <- vector("list", length(all_models))
names(preds_list) <- names(all_models)
names(probs_list) <- names(all_models)

pb <- txtProgressBar(min = 0, max = length(names(all_models)), style = 3)
i <- 0
for (model_name in names(all_models)){
  i <- i + 1
  cat("Model: ", model_name)
  model <- all_models[[model_name]]
  test_data <- test_data_map[[model_name]]
  performance <- get_performance(model, test_data, y_test,
                                 classes = c("cohesinAML", "wtAML"),
                                 plot_title = paste0(model_name, " model ROC curve"))
  
  preds_list[[model_name]] <- performance$predicted_classes
  probs_list[[model_name]] <- performance$predicted_probabilities
  
  # Extract metrics from confusion matrix
  cm <- performance$confusion_matrix
  sensitivity <- cm$byClass["Sensitivity"]
  specificity <- cm$byClass["Specificity"]
  precision <- cm$byClass["Pos Pred Value"]
  accuracy <- cm$overall["Accuracy"]
  kappa <- cm$overall["Kappa"]
  auc <- performance$auc
  
  performance_df <- rbind(performance_df, data.frame(model = model_name, sensitivity = sensitivity,
                                                     specificity = specificity, precision = precision,
                                                     accuracy = accuracy, kappa = kappa,
                                                     auc = auc, stringsAsFactors = FALSE))
}

# write.table(performance_df, "OutputTables/basic_all_models_CV_SMOTE_3000_performance.txt",
#             col.names=TRUE, sep="\t")
# saveRDS(preds_list, "OutputTables/basic_all_models_CV_SMOTE_3000_preds.rds")
# saveRDS(probs_list, "OutputTables/basic_all_models_CV_SMOTE_3000_probs.rds")


# Check missclassified samples =================================
TCGA_BEAT_coh_labels <- read.table("InputTables/TCGA_BEAT_labels_cohesin.txt", header=T)
performance_df <- read.table("OutputTables/basic_all_models_CV_SMOTE_3000_performance.txt", header=T)
performance_df <- performance_df %>%
  separate(model, into = c("model", "transformation"), sep = "_", extra = "merge") 
preds_list <- readRDS("OutputTables/basic_all_models_CV_SMOTE_3000_preds.rds")
probs_list <- readRDS("OutputTables/basic_all_models_CV_SMOTE_3000_probs.rds")

selected_models <- c("kNN_vst_min_max", "nb_rlog_QN",
                     "svm_vst_min_max", "glm_rlog_QN",
                     "rf_vst_zScore", "gbm_vst_min_max")

positive_samples <- which(y_test == "cohesinAML")


miss_df <- data.frame(sampleID = TCGA_BEAT_coh_labels[positive_samples,]$sample,
                      kNN_vst_min_max = preds_list$kNN_vst_min_max[positive_samples],
                      nb_rlog_QN = preds_list$nb_rlog_QN[positive_samples],
                      svm_vst_min_max = preds_list$svm_vst_min_max[positive_samples],
                      glm_rlog_QN = preds_list$glm_rlog_QN[positive_samples],
                      rf_vst_zScore = preds_list$rf_vst_zScore[positive_samples],
                      gbm_vst_min_max = preds_list$gbm_vst_min_max[positive_samples])

# Check any sample classified as WT by all models
only_wtAML <- apply(miss_df[ , -1], 1,
  function(x) all(x == "wtAML"))
any(only_wtAML)


# Get weights for each contributor model ================================
for (model_name in selected_models){
  cat("Model: ", model_name, "\n")
  model <- all_models[[model_name]]
  test_data <- test_data_map[[model_name]]
  performance <- get_performance(model, test_data, y_test,
                                 classes = c("cohesinAML", "wtAML"),
                                 plot_title = paste0(model_name, " model ROC curve"))
}

# Set weights
selected_models_weights <- c(kNN_vst_min_max = 0.18,
                             nb_rlog_QN = 0.12,
                             svm_vst_min_max = 0.11,
                             glm_rlog_QN = 0.14,
                             rf_vst_zScore = 0.10,
                             gbm_vst_min_max = 0.35)


# Predict with ensemble model ===========================
# Test data
ensemble_predictions <- ensemble_predict(all_models = all_models,
                                         selected_models = selected_models,
                                         test_data_map = test_data_map,
                                         positive_class = "cohesinAML",
                                         negative_class = "wtAML",
                                         threshold = 0.5,
                                         weights = selected_models_weights)
ensemble_cm <- confusionMatrix(ensemble_predictions$classes, y_test)
# write_confusion_matrix(ensemble_cm, "OutputTables/cohesin_validation_test_ensemble_cm.txt")


# Test on Fischer ===============================
# Load Fischer data
Fischer_rlog_QN <- read.table("InputTables/Input_Fischer_top3000_rlog_QN.txt", check.names = F)
Fischer_vst_min_max <- read.table("InputTables/Input_Fischer_top3000_vst_min_max.txt", check.names = F)
Fischer_vst_zScore <- read.table("InputTables/Input_Fischer_top3000_vst_zScore.txt", check.names = F)
Fischer_labels <- read.table("InputTables/Fischer_labels_cohesin.txt", header = T)
Fischer_labels <- c(Fischer_labels$label)
Fischer_labels <- as.factor(Fischer_labels)

# Set up data map
fischer_data_map <- list(kNN_vst_min_max = Fischer_vst_min_max,
                         nb_rlog_QN = Fischer_rlog_QN,
                         svm_vst_min_max = Fischer_vst_min_max,
                         glm_rlog_QN = Fischer_rlog_QN,
                         rf_vst_zScore = Fischer_vst_zScore,
                         gbm_vst_min_max = Fischer_vst_min_max)

# Get predictions
ensemble_predictions <- ensemble_predict(all_models = all_models,
                                         selected_models = selected_models,
                                         test_data_map = fischer_data_map,
                                         positive_class = "cohesinAML",
                                         negative_class = "wtAML",
                                         threshold = 0.5,
                                         weights = selected_models_weights)
ensemble_cm <- confusionMatrix(ensemble_predictions$classes, Fischer_labels)
# write_confusion_matrix(ensemble_cm, "OutputTables/cohesin_validation_Fischer_ensemble_cm.txt")













