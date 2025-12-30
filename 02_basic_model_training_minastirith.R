source("00_helper_functions.R")
library(tidyverse)
library(caret)
library(pROC)


# Load data and labels ========================
raw_data <- read.table("InputTables/Input_TCGA-BEAT_raw_counts_common_genes.txt", sep="\t", header=T)

vst_transform <- t(read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt", header=T))
vst_min_max <- read.table("InputTables/Input_TCGA-BEAT_vst_min_max.txt", header=T)
vst_zScore <- read.table("InputTables/Input_TCGA-BEAT_vst_zScore.txt", header=T)
vst_QN <- read.table("InputTables/Input_TCGA-BEAT_vst_QN.txt", header=T)

rlog_transform <- t(read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_corrFiltered.txt", header=T))
rlog_min_max <- read.table("InputTables/Input_TCGA-BEAT_rlog_min_max.txt", header=T)
rlog_zScore <- read.table("InputTables/Input_TCGA-BEAT_rlog_zScore.txt", header=T)
rlog_QN <- read.table("InputTables/Input_TCGA-BEAT_rlog_QN.txt", header=T)

train_index <- scan("InputTables/Input_train_indexes.txt", sep="\n")

# Load labels =================================
TCGA_BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels_cohesin.txt", header=T)
TCGA_BEAT_labels <- c(TCGA_BEAT_labels$label)

# Data splitting ==============================
# Labels
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

# Try on some models ============================
# Set train control obj
trControl <- trainControl(method = "none", classProbs = TRUE)

# k-NN
k_val <- round(sqrt(length(y_train)), 0)
print("Fitting model kNN_vst_min_max")
kNN_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="knn", trControl=trControl,
                         tuneGrid = data.frame(k = k_val))
print("Fitting model kNN_vst_zScore")
kNN_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="knn", trControl=trControl,
                        tuneGrid = data.frame(k = k_val))
print("Fitting model kNN_rlog_min_max")
kNN_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="knn", trControl=trControl,
                          tuneGrid = data.frame(k = k_val))
print("Fitting model kNN_rlog_zScore")
kNN_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="knn", trControl=trControl,
                         tuneGrid = data.frame(k = k_val))

# Naive Bayes
print("Fitting model nb_raw")
nb_raw <- train(x=x_train_raw, y=y_train, method="naive_bayes", trControl=trControl,
                tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_vst")
nb_vst <- train(x=x_train_vst, y=y_train, method="naive_bayes", trControl=trControl,
                tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_rlog")
nb_rlog <- train(x=x_train_rlog, y=y_train, method="naive_bayes", trControl=trControl,
                 tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_vst_min_max")
nb_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="naive_bayes", trControl=trControl,
                        tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_vst_zScore")
nb_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="naive_bayes", trControl=trControl,
                       tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_vst_QN")
nb_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="naive_bayes", trControl=trControl,
                   tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_rlog_min_max")
nb_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="naive_bayes", trControl=trControl,
                         tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_rlog_zScore")
nb_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="naive_bayes", trControl=trControl,
                        tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))
print("Fitting model nb_rlog_QN")
nb_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="naive_bayes", trControl=trControl,
                    tuneGrid = data.frame(usekernel = TRUE, laplace = 0, adjust = 1))

# SVM
print("Fitting model svm_vst_min_max")
svm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="svmLinear2", trControl=trControl,
                         tuneGrid = data.frame(cost = 1))
print("Fitting model svm_vst_zScore")
svm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="svmLinear2", trControl=trControl,
                        tuneGrid = data.frame(cost = 1))
print("Fitting model svm_rlog_min_max")
svm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="svmLinear2", trControl=trControl,
                          tuneGrid = data.frame(cost = 1))
print("Fitting model svm_rlog_zScore")
svm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="svmLinear2", trControl=trControl,
                         tuneGrid = data.frame(cost = 1))

# Logistic Regression
print("Fitting model glm_vst_min_max")
glm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="glm",
                         family="binomial", trControl=trControl)
print("Fitting model glm_vst_zScore")
glm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="glm",
                        family="binomial", trControl=trControl)
print("Fitting model glm_rlog_min_max")
glm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="glm",
                          family="binomial", trControl=trControl)
print("Fitting model glm_rlog_zScore")
glm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="glm",
                         family="binomial", trControl=trControl)

# Random Forest
mtry_val <- floor(sqrt(ncol(x_train_raw)))

print("Fitting model rf_raw")
rf_raw <- train(x=x_train_raw, y=y_train, method="rf", trControl=trControl,
                tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_vst")
rf_vst <- train(x=x_train_vst, y=y_train, method="rf", trControl=trControl,
                tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_rlog")
rf_rlog <- train(x=x_train_rlog, y=y_train, method="rf", trControl=trControl,
                 tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_vst_min_max")
rf_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="rf", trControl=trControl,
                        tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_vst_zScore")
rf_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="rf", trControl=trControl,
                       tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_vst_QN")
rf_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="rf", trControl=trControl,
                   tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_rlog_min_max")
rf_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="rf", trControl=trControl,
                         tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_rlog_zScore")
rf_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="rf", trControl=trControl,
                        tuneGrid = data.frame(mtry = mtry_val))
print("Fitting model rf_rlog_QN")
rf_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="rf", trControl=trControl,
                    tuneGrid = data.frame(mtry = mtry_val))

# Boosted Decision Trees
print("Fitting model gbm_raw")
gbm_raw <- train(x=x_train_raw, y=y_train, method="gbm", trControl=trControl,
                       tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                             n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_vst")
gbm_vst <- train(x=x_train_vst, y=y_train, method="gbm", trControl=trControl,
                 tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                       n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_rlog")
gbm_rlog <- train(x=x_train_rlog, y=y_train, method="gbm", trControl=trControl,
                  tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                        n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_vst_min_max")
gbm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="gbm", trControl=trControl,
                         tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                               n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_vst_zScore")
gbm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="gbm", trControl=trControl,
                        tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                              n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_vst_QN")
gbm_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="gbm", trControl=trControl,
                    tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                          n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_rlog_min_max")
gbm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="gbm", trControl=trControl,
                          tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                                n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_rlog_zScore")
gbm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="gbm", trControl=trControl,
                         tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                               n.minobsinnode = 10), verbose = FALSE)
print("Fitting model gbm_rlog_QN")
gbm_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="gbm", trControl=trControl,
                     tuneGrid = data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                                           n.minobsinnode = 10), verbose = FALSE)


# Save models =================================
# List of all models
all_models <- list(
  kNN_vst_min_max = kNN_vst_min_max,
  kNN_vst_zScore = kNN_vst_zScore,
  kNN_rlog_min_max = kNN_rlog_min_max,
  kNN_rlog_zScore = kNN_rlog_zScore,
  
  nb_raw = nb_raw,
  nb_vst = nb_vst,
  nb_rlog = nb_rlog,
  nb_vst_min_max = nb_vst_min_max,
  nb_vst_zScore = nb_vst_zScore,
  nb_vst_QN = nb_vst_QN,
  nb_rlog_min_max = nb_rlog_min_max,
  nb_rlog_zScore = nb_rlog_zScore,
  nb_rlog_QN = nb_rlog_QN,
  
  svm_vst_min_max = svm_vst_min_max,
  svm_vst_zScore = svm_vst_zScore,
  svm_rlog_min_max = svm_rlog_min_max,
  svm_rlog_zScore = svm_rlog_zScore,
  
  glm_vst_min_max = glm_vst_min_max,
  glm_vst_zScore = glm_vst_zScore,
  glm_rlog_min_max = glm_rlog_min_max,
  glm_rlog_zScore = glm_rlog_zScore,
  
  rf_raw = rf_raw,
  rf_vst = rf_vst,
  rf_rlog = rf_rlog,
  rf_vst_min_max = rf_vst_min_max,
  rf_vst_zScore = rf_vst_zScore,
  rf_vst_QN = rf_vst_QN,
  rf_rlog_min_max = rf_rlog_min_max,
  rf_rlog_zScore = rf_rlog_zScore,
  rf_rlog_QN = rf_rlog_QN,
  
  gbm_raw = gbm_raw,
  gbm_vst = gbm_vst,
  gbm_rlog = gbm_rlog,
  gbm_vst_min_max = gbm_vst_min_max,
  gbm_vst_zScore = gbm_vst_zScore,
  gbm_vst_QN = gbm_vst_QN,
  gbm_rlog_min_max = gbm_rlog_min_max,
  gbm_rlog_zScore = gbm_rlog_zScore,
  gbm_rlog_QN = gbm_rlog_QN
)

# Save each model as an RDS file
# saveRDS(all_models, file = "OutputTables/basic_all_models.rds")
all_models <- readRDS("OutputTables/basic_all_models.rds")
load("02_basic_model_training.RData")

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
  
  glm_vst_min_max = x_test_vst_min_max, glm_vst_zScore = x_test_vst_zScore,
  glm_rlog_min_max = x_test_rlog_min_max, glm_rlog_zScore = x_test_rlog_zScore,
  
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

pb <- txtProgressBar(min = 0, max = length(names(all_models)), style = 3)
i <- 0
for (model_name in names(all_models)){
  i <- i + 1
  setTxtProgressBar(pb, i)
  
  model <- all_models[[model_name]]
  test_data <- test_data_map[[model_name]]
  performance <- get_performance(model, test_data, y_test,
                                 plot_title = paste0(model_name, " model ROC curve"))
  
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
close(pb)

write.table(performance_df, "OutputTables/basic_all_models_performance.txt",
            col.names=TRUE, sep="\t")



