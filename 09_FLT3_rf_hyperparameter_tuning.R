source("00_helper_functions.R")
library(tidyverse)
library(caret)
library(doParallel)


# Load data =================================
# X data
vst_QN <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_QN.txt", header=T)
Fischer_vst_QN <- read.table("InputTables/Input_Fischer_top3000_vst_QN.txt", header=T)
stopifnot(identical(colnames(Fischer_vst_QN), colnames(vst_QN)))

# y data 
# Load labels 
TCGA_BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels_FLT3.txt", header=T)
TCGA_BEAT_labels <- c(TCGA_BEAT_labels$x)
Fischer_labels <- read.table("InputTables/Fischer_labels_FLT3.txt", header=T)
Fischer_labels <- c(Fischer_labels$label)
Fischer_labels <- as.factor(Fischer_labels)

# Load partition indexes
train_index <- scan("InputTables/Input_train_indexes.txt", sep="\n")

# Splitting data
y_train <- as.factor(TCGA_BEAT_labels[train_index])
y_test <- as.factor(TCGA_BEAT_labels[-train_index])
x_train_vst_QN <- vst_QN[train_index, ]
x_test_vst_QN <- vst_QN[-train_index, ]


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
                          classProbs = TRUE, sampling = "smote", savePredictions = "final",
                          summaryFunction = customSummary, allowParallel = FALSE)


# Set grid of parameters ===================================
mtry_vals <- c(40, 55, 70)
ntree_vals <- c(500, 1000, 2000)
nodesize_vals <- c(3, 5, 10)
maxnodes_vals <- c(30, 50, 80)

param_grid <- expand.grid(mtry = mtry_vals,
                          ntree = ntree_vals,
                          nodesize = nodesize_vals,
                          maxnodes = maxnodes_vals)

# Fit models ======================================
rf_results <- vector("list", nrow(param_grid))

for (i in seq_len(nrow(param_grid))) {
  params <- param_grid[i, ]
  
  cat("\nTraining model", i, "of", nrow(param_grid), "\n")
  cat("mtry =", params$mtry,
      "\nntree =", params$ntree,
      "\nnodesize =", params$nodesize,
      "\nmaxnodes =", params$maxnodes, "\n")
  
  start_time <- Sys.time()
  
  rf_results[[i]] <- train(x = x_train_vst_QN, y = y_train,
                           method = "rf", metric = "ROC", trControl = trControl,
                           tuneGrid = expand.grid(mtry = params$mtry),
                           ntree = params$ntree, nodesize = params$nodesize,
                           maxnodes = params$maxnodes, importance = TRUE, replace = TRUE)
  
  cat("Time to train model:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")
}


# saveRDS(rf_results, "OutputTables/FLT3_tuned_rf_vst_QN.rds")
rf_results <- readRDS("OutputTables/FLT3_tuned_rf_vst_QN.rds")

# Store performance metrics to the grid
test_performance_grid <- param_grid
test_performance_grid$Accuracy <- NA_real_
test_performance_grid$Sens <- NA_real_
test_performance_grid$Spec <- NA_real_
test_performance_grid$Prec <- NA_real_
test_performance_grid$Kappa <- NA_real_
test_performance_grid$AUC <- NA_real_
test_performance_roc <- vector("list", nrow(test_performance_grid))

for (i in seq_len(nrow(test_performance_grid))) {
  model <- rf_results[[i]]
  plot_title <- paste0("ROC curve for model ", i)
  plot_subtitle <- ""
  for (j in 1:4){
    parameter <- colnames(test_performance_grid)[j]
    value <- test_performance_grid[i,j]
    plot_subtitle <- paste0(plot_subtitle, parameter, "=", value, " ")
  }
  performance_df <- get_performance(rf_results[[i]], x_test_vst_QN, y_test,
                                    classes = c(levels(y_test)[1], levels(y_test)[2]),
                                    plot_title = plot_title, plot_subtitle = plot_subtitle)
  
  # Extract metrics from confusion matrix
  cm <- performance_df$confusion_matrix
  test_performance_grid$Sens[i] <- cm$byClass["Sensitivity"]
  test_performance_grid$Spec[i] <- cm$byClass["Specificity"]
  test_performance_grid$Prec[i] <- cm$byClass["Pos Pred Value"]
  test_performance_grid$Accuracy[i] <- cm$overall["Accuracy"]
  test_performance_grid$Kappa[i] <- cm$overall["Kappa"]
  test_performance_grid$AUC[i] <- performance_df$auc
  
  test_performance_roc[[i]] <- performance_df$roc_df
}

names(test_performance_roc) <- paste0("model ", seq_along(test_performance_roc))
test_performance_grid <- cbind(test_performance_grid,
                               data.frame(model=paste0("model ",1:length(test_performance_grid$mtry))))





