source("00_helper_functions.R")
library(tidyverse)
library(caret)


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


# Define custom function to evaluate models =================
customSummary <- function(data, lev = NULL, model = NULL) {
  c(twoClassSummary(data, lev, model),
    defaultSummary(data, lev, model))
}


# Set up parallelization ==========================
# Detect cores allocated by Slurm
cores <- as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
if (is.na(cores) || cores <= 0) {
  cores <- parallel::detectCores() - 1
}

# Create cluster
cl <- makeCluster(cores)
registerDoParallel(cl)

# Make sure parallelization stop when exit
on.exit({
  stopCluster(cl)
  registerDoSEQ()
}, add = TRUE)


# Set resampling methods ==========================
trControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                          classProbs = TRUE, summaryFunction = customSummary,
                          sampling = "smote")


# Set tune grid ===================================
gbmGrid <- expand.grid(n.trees = c(100, 500, 1000, 3000),
                       interaction.depth = c(1, 2, 3, 4),
                       shrinkage = c(0.1, 0.05, 0.01, 0.005),
                       n.minobsinnode = c(10, 30, 50))


# Fit models ======================================
start_time <- Sys.time()
print("Fitting model gbm_vst_min_max")
gbm_tuned <- train(x=x_train_vst_min_max, y=y_train, method="gbm",
                   trControl=trControl, tuneGrid = gbmGrid,
                   metric = "Sens", verbose = T)
print("Time took to compute models:")
end_time <- Sys.time()
end_time - start_time

# saveRDS(gbm_tuned, "OutputTables/cohesin_tuned_gbm_vst_min_max.rds")
gbm_tuned <- readRDS("OutputTables/cohesin_tuned_gbm_vst_min_max.rds")


# Obtain top 10 best performing models ====================
top10_grid <- gbm_tuned$results %>%
  arrange(desc(Sens)) %>%
  head(10) %>%
  select(n.trees, interaction.depth, shrinkage, n.minobsinnode)

# List to save the models
top10_gbm <- vector("list", nrow(top10_grid))
for (i in seq_len(nrow(top10_grid))) {
  
  cat("Training model", i, "of", nrow(top10_grid), "\n")
  
  top10_gbm[[i]] <- train(
    x = x_train_vst_min_max,
    y = y_train,
    method = "gbm",
    trControl = trControl,
    tuneGrid = top10_grid[i, , drop = FALSE],
    verbose = FALSE
  )
}
# saveRDS(top10_gbm, file = "OutputTables/cohesin_top10_gbm_models.rds")
top10_gbm <- readRDS("OutputTables/cohesin_tuned_top10_gbm_models.rds")


# Assess performance on test data ===================================
# Single best tuned gbm
  ## Check ROC curve
performance <- get_performance(gbm_tuned, x_test_vst_min_max, y_test,
                               classes = c("cohesinAML", "wtAML"),
                               plot_title = "ROC curve for gbm with highest Sensitivity (on test data)",
                               plot_subtitle=paste0("n.trees=", gbm_tuned[["bestTune"]][["n.trees"]],
                                                    " interaction.depth=", gbm_tuned[["bestTune"]][["interaction.depth"]],
                                                    " shrinkage=", gbm_tuned[["bestTune"]][["shrinkage"]],
                                                    " n.minobsinnode=", gbm_tuned[["bestTune"]][["n.minobsinnode"]]))
# ggsave("plots/ROC_test_tuned_gbm_vst_min_max_3000.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)
# write_confusion_matrix(performance$confusion_matrix,
#                        "OutputTables/cohesin_validation_test_tuned_gbm_vst_min_max_3000_cm.txt")

  ## Get class probabilities
probs <- predict(gbm_tuned, Fischer_vst_min_max, type = "prob")[, "cohesinAML"]
quantile(probs[Fischer_labels == "cohesinAML"])  # cohesinAML distribution
quantile(probs[Fischer_labels == "wtAML"])       # wtAML distribution

# Top 10 best performing gbms
test_performance_grid <- top10_grid
test_performance_grid$Accuracy <- NA_real_
test_performance_grid$Sens <- NA_real_
test_performance_grid$Spec <- NA_real_
test_performance_grid$Prec <- NA_real_
test_performance_grid$Kappa <- NA_real_
test_performance_grid$AUC <- NA_real_
test_performance_top10_roc <- vector("list", nrow(test_performance_grid))

for (i in seq_len(nrow(test_performance_grid))) {
  
  cat("Assessing performance of model", i, "of", nrow(test_performance_grid), "\n")
  
  model <- top10_gbm[[i]]
  plot_title <- paste0("ROC curve for model ", i)
  plot_subtitle <- ""
  for (j in 1:4){
    parameter <- colnames(top10_gbm[[i]][["results"]])[j]
    value <- top10_gbm[[i]][["results"]][[parameter]]
    plot_subtitle <- paste0(plot_subtitle, parameter, "=", value, " ")
  }
  performance_df <- get_performance(top10_gbm[[i]], x_test_vst_min_max, y_test,
                                    classes = c("cohesinAML", "wtAML"),
                                    plot_title = plot_title, plot_subtitle = plot_subtitle)
  
  # Extract metrics from confusion matrix
  cm <- performance_df$confusion_matrix
  test_performance_grid$Sens[i] <- cm$byClass["Sensitivity"]
  test_performance_grid$Spec[i] <- cm$byClass["Specificity"]
  test_performance_grid$Prec[i] <- cm$byClass["Pos Pred Value"]
  test_performance_grid$Accuracy[i] <- cm$overall["Accuracy"]
  test_performance_grid$Kappa[i] <- cm$overall["Kappa"]
  test_performance_grid$AUC[i] <- performance_df$auc
  
  test_performance_top10_roc[[i]] <- performance_df$roc_df
}

names(test_performance_top10_roc) <- paste0("model ", seq_along(test_performance_top10_roc))
test_performance_grid <- cbind(test_performance_grid, data.frame(model=paste0("model ", 1:10)))
# write.table(test_performance_grid, "OutputTables/cohesin_tuned_top10_gbm_models_validation_test_performance_table.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")
# saveRDS(test_performance_top10_roc, "OutputTables/cohesin_tuned_top10_gbm_models_validation_test_performance_roc.rds")


# Plot test ROC curves ==============================
roc_5_df <- bind_rows(lapply(1:5, function(i) {
  df <- test_performance_top10_roc[[i]]
  df$model <- paste0("Model ", i)
  df
  }))
auc_labels <- paste0("Model ", 1:5,
                     " (AUC = ",
                     round(test_performance_grid$AUC[1:5], 3), ")")
roc_5_df$model <- factor(roc_5_df$model,
                         levels = paste0("Model ", 1:5),
                         labels = auc_labels)
ggplot(roc_5_df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves for Top 5 gbm models on test data",
       x = "False Positive Rate", y = "True Positive Rate",
       color = "Model") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
# ggsave("plots/ROC_test_top5_gbm_vst_min_max_3000.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)


# Assess performance on Fischer data ===================================
# Single best tuned gbm
  ## Check ROC curve
performance <- get_performance(gbm_tuned, Fischer_vst_min_max, Fischer_labels,
                               classes = c("cohesinAML", "wtAML"),
                               plot_title = "ROC curve for gbm with highest Sensitivity (on Fischer data)",
                               plot_subtitle=paste0("n.trees=", gbm_tuned[["bestTune"]][["n.trees"]],
                                                    " interaction.depth=", gbm_tuned[["bestTune"]][["interaction.depth"]],
                                                    " shrinkage=", gbm_tuned[["bestTune"]][["shrinkage"]],
                                                    " n.minobsinnode=", gbm_tuned[["bestTune"]][["n.minobsinnode"]]))
# ggsave("plots/ROC_Fischer_tuned_gbm_vst_min_max_3000.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)
# write_confusion_matrix(performance$confusion_matrix,
#                        "OutputTables/cohesin_validation_Fischer_tuned_gbm_vst_min_max_3000_cm.txt")

  ## Get class probabilities
probs <- predict(gbm_tuned, Fischer_vst_min_max, type = "prob")[, "cohesinAML"]
quantile(probs[Fischer_labels == "cohesinAML"])  # cohesinAML distribution
quantile(probs[Fischer_labels == "wtAML"])       # wtAML distribution

# Top 10 best performing gbms
fischer_performance_grid <- top10_grid
fischer_performance_grid$Accuracy <- NA_real_
fischer_performance_grid$Sens <- NA_real_
fischer_performance_grid$Spec <- NA_real_
fischer_performance_grid$Prec <- NA_real_
fischer_performance_grid$Kappa <- NA_real_
fischer_performance_grid$AUC <- NA_real_
fischer_top10_roc <- vector("list", nrow(fischer_performance_grid))

for (i in seq_len(nrow(fischer_performance_grid))) {
  
  cat("Assessing performance of model", i, "of", nrow(fischer_performance_grid), "\n")
  
  model <- top10_gbm[[i]]
  plot_title <- paste0("ROC curve for model ", i)
  plot_subtitle <- ""
  for (j in 1:4){
    parameter <- colnames(top10_gbm[[i]][["results"]])[j]
    value <- top10_gbm[[i]][["results"]][[parameter]]
    plot_subtitle <- paste0(plot_subtitle, parameter, "=", value, " ")
  }
  performance_df <- get_performance(top10_gbm[[i]], Fischer_vst_min_max, Fischer_labels,
                                    classes = c("cohesinAML", "wtAML"),
                                    plot_title = plot_title, plot_subtitle = plot_subtitle)
  
  # Extract metrics from confusion matrix
  cm <- performance_df$confusion_matrix
  fischer_performance_grid$Sens[i] <- cm$byClass["Sensitivity"]
  fischer_performance_grid$Spec[i] <- cm$byClass["Specificity"]
  fischer_performance_grid$Prec[i] <- cm$byClass["Pos Pred Value"]
  fischer_performance_grid$Accuracy[i] <- cm$overall["Accuracy"]
  fischer_performance_grid$Kappa[i] <- cm$overall["Kappa"]
  fischer_performance_grid$AUC[i] <- performance_df$auc
  
  fischer_top10_roc[[i]] <- performance_df$roc_df
}

names(fischer_top10_roc) <- paste0("model ", seq_along(fischer_top10_roc))
fischer_performance_grid <- cbind(fischer_performance_grid, data.frame(model=paste0("model ", 1:10)))
# write.table(fischer_performance_grid, "OutputTables/cohesin_tuned_top10_gbm_models_validation_Fischer_performance_table.txt",
#             row.names = T, col.names = T, quote = F, sep = "\t")
# saveRDS(fischer_top10_roc, "OutputTables/cohesin_tuned_top10_gbm_models_validation_Fischer_performance_roc.rds")


# Plot Fischer ROC curves ==============================
roc_5_df <- bind_rows(lapply(1:5, function(i) {
  df <- fischer_top10_roc[[i]]
  df$model <- paste0("Model ", i)
  df
}))
auc_labels <- paste0("Model ", 1:5,
                     " (AUC = ",
                     round(fischer_performance_grid$AUC[1:5], 3), ")")
roc_5_df$model <- factor(roc_5_df$model,
                         levels = paste0("Model ", 1:5),
                         labels = auc_labels)
ggplot(roc_5_df, aes(x = fpr, y = tpr, color = model)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves for Top 5 gbm models on Fischer data",
       x = "False Positive Rate", y = "True Positive Rate",
       color = "Model") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
# ggsave("plots/ROC_Fischer_top5_gbm_vst_min_max_3000.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)


