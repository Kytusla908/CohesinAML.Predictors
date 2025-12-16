source("00_helper_functions.R")
library(tidyverse)
library(caret)
library(pROC)


# Load data and labels ========================
load("02_train_test_split.RData")


# Try on some models ============================
# Set ROSE resampling and rCV
trControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, sampling = "smote")

# k-NN
k_val <- round(sqrt(length(y_train)), 0)
tuneGrid <- data.frame(k = k_val)
print("Fitting model kNN_vst_min_max")
kNN_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="knn", trControl=trControl,
                         tuneGrid = tuneGrid)
print("Fitting model kNN_vst_zScore")
kNN_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="knn", trControl=trControl,
                        tuneGrid = tuneGrid)
print("Fitting model kNN_rlog_min_max")
kNN_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="knn", trControl=trControl,
                          tuneGrid = tuneGrid)
print("Fitting model kNN_rlog_zScore")
kNN_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="knn", trControl=trControl,
                         tuneGrid = tuneGrid)

# Naive Bayes
tuneGrid <- data.frame(usekernel = TRUE, laplace = 0, adjust = 1)
print("Fitting model nb_vst")
nb_vst <- train(x=x_train_vst, y=y_train, method="naive_bayes", trControl=trControl,
                tuneGrid = tuneGrid)
print("Fitting model nb_rlog")
nb_rlog <- train(x=x_train_rlog, y=y_train, method="naive_bayes", trControl=trControl,
                 tuneGrid = tuneGrid)
print("Fitting model nb_vst_min_max")
nb_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="naive_bayes", trControl=trControl,
                        tuneGrid = tuneGrid)
print("Fitting model nb_vst_zScore")
nb_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="naive_bayes", trControl=trControl,
                       tuneGrid = tuneGrid)
print("Fitting model nb_vst_QN")
nb_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="naive_bayes", trControl=trControl,
                   tuneGrid = tuneGrid)
print("Fitting model nb_rlog_min_max")
nb_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="naive_bayes", trControl=trControl,
                         tuneGrid = tuneGrid)
print("Fitting model nb_rlog_zScore")
nb_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="naive_bayes", trControl=trControl,
                        tuneGrid = tuneGrid)
print("Fitting model nb_rlog_QN")
nb_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="naive_bayes", trControl=trControl,
                    tuneGrid = tuneGrid)

# SVM
tuneGrid <- data.frame(cost = 1)
print("Fitting model svm_vst_min_max")
svm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="svmLinear2", trControl=trControl,
                         tuneGrid = tuneGrid)
print("Fitting model svm_vst_zScore")
svm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="svmLinear2", trControl=trControl,
                        tuneGrid = tuneGrid)
print("Fitting model svm_rlog_min_max")
svm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="svmLinear2", trControl=trControl,
                          tuneGrid = tuneGrid)
print("Fitting model svm_rlog_zScore")
svm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="svmLinear2", trControl=trControl,
                         tuneGrid = tuneGrid)

# Logistic Regression
print("Fitting model glm_vst")
glm_vst <- train(x=x_train_vst, y=y_train, method="glm",
                 family="binomial", trControl=trControl)
print("Fitting model glm_rlog")
glm_rlog <- train(x=x_train_rlog, y=y_train, method="glm",
                  family="binomial", trControl=trControl)
print("Fitting model glm_vst_min_max")
glm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="glm",
                         family="binomial", trControl=trControl)
print("Fitting model glm_vst_zScore")
glm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="glm",
                        family="binomial", trControl=trControl)
print("Fitting model glm_vst_QN")
glm_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="glm",
                        family="binomial", trControl=trControl)
print("Fitting model glm_rlog_min_max")
glm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="glm",
                          family="binomial", trControl=trControl)
print("Fitting model glm_rlog_zScore")
glm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="glm",
                         family="binomial", trControl=trControl)
print("Fitting model glm_rlog_QN")
glm_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="glm",
                    family="binomial", trControl=trControl)

# Random Forest
mtry_val <- floor(sqrt(ncol(x_train_raw)))
tuneGrid <- data.frame(mtry = mtry_val)

print("Fitting model rf_vst")
rf_vst <- train(x=x_train_vst, y=y_train, method="rf", trControl=trControl,
                tuneGrid = tuneGrid)
print("Fitting model rf_rlog")
rf_rlog <- train(x=x_train_rlog, y=y_train, method="rf", trControl=trControl,
                 tuneGrid = tuneGrid)
print("Fitting model rf_vst_min_max")
rf_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="rf", trControl=trControl,
                        tuneGrid = tuneGrid)
print("Fitting model rf_vst_zScore")
rf_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="rf", trControl=trControl,
                       tuneGrid = tuneGrid)
print("Fitting model rf_vst_QN")
rf_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="rf", trControl=trControl,
                   tuneGrid = tuneGrid)
print("Fitting model rf_rlog_min_max")
rf_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="rf", trControl=trControl,
                         tuneGrid = tuneGrid)
print("Fitting model rf_rlog_zScore")
rf_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="rf", trControl=trControl,
                        tuneGrid = tuneGrid)
print("Fitting model rf_rlog_QN")
rf_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="rf", trControl=trControl,
                    tuneGrid = tuneGrid)

# Boosted Decision Trees
tuneGrid <- data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                      n.minobsinnode = 10)

print("Fitting model gbm_vst")
gbm_vst <- train(x=x_train_vst, y=y_train, method="gbm", trControl=trControl,
                 tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_rlog")
gbm_rlog <- train(x=x_train_rlog, y=y_train, method="gbm", trControl=trControl,
                  tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_vst_min_max")
gbm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="gbm", trControl=trControl,
                         tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_vst_zScore")
gbm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="gbm", trControl=trControl,
                        tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_vst_QN")
gbm_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="gbm", trControl=trControl,
                    tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_rlog_min_max")
gbm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="gbm", trControl=trControl,
                          tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_rlog_zScore")
gbm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="gbm", trControl=trControl,
                         tuneGrid = tuneGrid, verbose = FALSE)
print("Fitting model gbm_rlog_QN")
gbm_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="gbm", trControl=trControl,
                     tuneGrid = tuneGrid, verbose = FALSE)


# Save models =================================
# List of all models
all_models <- list(
  kNN_vst_min_max = kNN_vst_min_max,
  kNN_vst_zScore = kNN_vst_zScore,
  kNN_rlog_min_max = kNN_rlog_min_max,
  kNN_rlog_zScore = kNN_rlog_zScore,
  
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
  
  glm_vst = glm_vst,
  glm_rlog = glm_rlog,
  glm_vst_min_max = glm_vst_min_max,
  glm_vst_zScore = glm_vst_zScore,
  glm_vst_QN = glm_vst_QN,
  glm_rlog_min_max = glm_rlog_min_max,
  glm_rlog_zScore = glm_rlog_zScore,
  glm_rlog_QN = glm_rlog_QN,
  
  rf_vst = rf_vst,
  rf_rlog = rf_rlog,
  rf_vst_min_max = rf_vst_min_max,
  rf_vst_zScore = rf_vst_zScore,
  rf_vst_QN = rf_vst_QN,
  rf_rlog_min_max = rf_rlog_min_max,
  rf_rlog_zScore = rf_rlog_zScore,
  rf_rlog_QN = rf_rlog_QN,
  
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
# saveRDS(all_models, file = "OutputTables/basic_all_models_rCV_ROSE.rds")

# Map test_data ================================
test_data_map <- list(
  kNN_vst_min_max = x_test_vst_min_max,
  kNN_vst_zScore = x_test_vst_zScore,
  kNN_rlog_min_max = x_test_rlog_min_max,
  kNN_rlog_zScore = x_test_rlog_zScore,
  
  nb_vst = x_test_vst, nb_rlog = x_test_rlog,
  nb_vst_min_max = x_test_vst_min_max, nb_vst_zScore = x_test_vst_zScore,
  nb_vst_QN = x_test_vst_QN, nb_rlog_min_max = x_test_rlog_min_max,
  nb_rlog_zScore = x_test_rlog_zScore, nb_rlog_QN = x_test_rlog_QN,
  
  svm_vst_min_max = x_test_vst_min_max, svm_vst_zScore = x_test_vst_zScore,
  svm_rlog_min_max = x_test_rlog_min_max, svm_rlog_zScore = x_test_rlog_zScore,
  
  glm_vst = x_test_vst, glm_rlog = x_test_rlog,
  glm_vst_min_max = x_test_vst_min_max, glm_vst_zScore = x_test_vst_zScore,
  glm_vst_QN = x_test_vst_QN, glm_rlog_min_max = x_test_rlog_min_max,
  glm_rlog_zScore = x_test_rlog_zScore, glm_rlog_QN = x_test_rlog_QN,
  
  rf_vst = x_test_vst, rf_rlog = x_test_rlog,
  rf_vst_min_max = x_test_vst_min_max, rf_vst_zScore = x_test_vst_zScore,
  rf_vst_QN = x_test_vst_QN, rf_rlog_min_max = x_test_rlog_min_max,
  rf_rlog_zScore = x_test_rlog_zScore, rf_rlog_QN = x_test_rlog_QN,
  
  gbm_vst = x_test_vst, gbm_rlog = x_test_rlog,
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

# write.table(performance_df, "OutputTables/basic_all_models_rCV_ROSE_performance.txt",
#             col.names=TRUE, sep="\t")

performance_df <- read.table("OutputTables/basic_all_models_performance.txt", header=T)
performance_df <- performance_df %>%
  separate(model, into = c("model", "transformation"), sep = "_", extra = "merge") 

# Plot Sensitivity
ggplot(performance_df, aes(x = transformation, y = sensitivity, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Sensitivity", color = "Model")
# ggsave("plots/basic_model_performance_sensitivity.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)

# Plot Specificity
ggplot(performance_df, aes(x = transformation, y = specificity, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Specificity", color = "Model")
# ggsave("plots/basic_model_performance_specificity.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)

# Plot kappa
ggplot(performance_df, aes(x = transformation, y = kappa, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Kappa", color = "Model")
# ggsave("plots/basic_model_performance_kappa.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)

# Plot AUC
ggplot(performance_df, aes(x = transformation, y = auc, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Kappa", color = "Model")
# ggsave("plots/basic_model_performance_AUC.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)



