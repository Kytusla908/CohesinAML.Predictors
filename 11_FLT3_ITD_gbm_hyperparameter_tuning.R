source("00_helper_functions.R")
library(tidyverse)
library(caret)
library(doParallel)


# Load data =================================
# X data
vst_min_max <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_min_max.txt", header=T)
Fischer_vst_min_max <- read.table("InputTables/Input_Fischer_top3000_vst_min_max_ComBat.txt")
stopifnot(identical(colnames(Fischer_vst_min_max), colnames(vst_min_max)))

# y data 
# Load labels 
TCGA_BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels_FLT3_ITD.txt", header=T)
TCGA_BEAT_labels <- as.factor(c(TCGA_BEAT_labels$x))
Fischer_labels <- read.table("InputTables/Fischer_labels_FLT3_ITD.txt", header=T)
Fischer_labels <- as.factor(c(Fischer_labels$label))

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
    prSummary(data, lev, model),
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
                          classProbs = TRUE, sampling = "smote", savePredictions = "final",
                          summaryFunction = customSummary, allowParallel = FALSE)


# Set grid of parameters ===================================
# Small version
gbmGrid <- expand.grid(n.trees = c(600, 1200, 2500),
                       interaction.depth = c(3, 5, 7),
                       shrinkage = c(0.05, 0.01, 0.005),
                       n.minobsinnode = c(10, 30, 50))
# Not-so-Small version
gbmGrid <- expand.grid(n.trees = c(3000, 5000),
                       interaction.depth = c(3, 5, 7),
                       shrinkage = c(0.05, 0.01, 0.005),
                       n.minobsinnode = c(10, 30, 50))


# Fit models ======================================
start_time <- Sys.time()
print("Fitting model gbm_vst_min_max")
gbm_tuned <- train(x=x_train_vst_min_max, y=y_train, method="gbm",
                   trControl=trControl, tuneGrid = gbmGrid,
                   metric = "Sens", verbose = T)
print("Time took to compute models:")
print(Sys.time()-start_time)

# saveRDS(gbm_tuned, "OutputTables/FLT3_ITD_tuned_gbm_vst_min_max_smaller.rds")
# saveRDS(gbm_tuned, "OutputTables/FLT3_ITD_tuned_gbm_vst_min_max_not_so_small.rds")
# gbm_tuned_small <- readRDS("OutputTables/FLT3_ITD_tuned_gbm_vst_min_max_smaller.rds")
# gbm_tuned_big <- readRDS("OutputTables/FLT3_ITD_tuned_gbm_vst_min_max_not_so_small.rds")


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
# saveRDS(top10_gbm, file = "OutputTables/FLT3_ITD_top10_gbm_models_smaller.rds")
# saveRDS(top10_gbm, file = "OutputTables/FLT3_ITD_top10_gbm_models_not_so_small.rds")
top10_gbm_small <- readRDS("OutputTables/FLT3_ITD_top10_gbm_models_smaller.rds")
top10_gbm_big <- readRDS("OutputTables/FLT3_ITD_top10_gbm_models_not_so_small.rds")


# Check performance ====================================
performance <- get_performance(gbm_tuned_small, x_test_vst_min_max, y_test,
                               classes = c("FLT3.Mutant", "FLT3.WT"))
# write_confusion_matrix(performance$confusion_matrix,
#                        "OutputTables/FLT3_validation_test_gbm_cm.txt")
performance <- get_performance(gbm_tuned_small, Fischer_vst_min_max, Fischer_labels,
                               classes = c("FLT3.Mutant", "FLT3.WT"))
# write_confusion_matrix(performance$confusion_matrix,
#                        "OutputTables/FLT3_validation_Fischer_gbm_cm.txt")


# Top 10 best performing gbms ========================
fischer_top5_roc <- vector("list", 5)
fischer_top5_auc <- c()
for (i in 1:5) {
  model <- top10_gbm_small[[i]]
  performance_df <- get_performance(top10_gbm_small[[i]], Fischer_vst_min_max, Fischer_labels,
                                    classes = c("FLT3_ITD.Mutant", "FLT3_ITD.WT"))
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
  labs(title = "ROC Curves for Top 5 gbm models for FLT3_ITD mutation on Fischer data",
       x = "False Positive Rate", y = "True Positive Rate",
       color = "Model") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
# ggsave("plots/FLT3/ROC_Fischer_top5_gbm_vst_min_max_3000.png", device = "png", width = 18, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)

performance <- get_performance(top10_gbm_small[[5]], Fischer_vst_min_max, Fischer_labels,
                               classes = c("FLT3_ITD.Mutant", "FLT3_ITD.WT"),
                               plot_title = "Performance of gbm model 5 for FLT3 mutation on Fischer data",
                               plot_subtitle=paste0("n.trees=", top10_gbm_small[[5]][["bestTune"]][["n.trees"]],
                                                    " interaction.depth=", top10_gbm_small[[5]][["bestTune"]][["interaction.depth"]],
                                                    " shrinkage=", top10_gbm_small[[5]][["bestTune"]][["shrinkage"]],
                                                    " n.minobsinnode=", top10_gbm_small[[5]][["bestTune"]][["n.minobsinnode"]]))


# Check Variable importance ===========================
varImp_df <- data.frame(variable=row.names(varImp(svm_tuned)$importance),
                        importance=varImp(svm_tuned)$importance)
varImp_top20 <- varImp_df %>%
  arrange(desc(importance.NPM1.Mutant)) %>%
  head(20)
ggplot(varImp_top20, aes(x = reorder(variable, importance.NPM1.Mutant), y = importance.NPM1.Mutant)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Variable", y = "Importance",
       title = "Top 20 most important genes for NPM1 mutation predicition (tuned SVM)")
# ggsave("plots/NPM1/VariableImportance_plot_svm_tuned.png", device = "png",
#        width = 18, height = 12, units = "cm", pointsize = 10, dpi = 500)










