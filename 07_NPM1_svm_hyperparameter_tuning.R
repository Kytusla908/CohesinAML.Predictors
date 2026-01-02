source("00_helper_functions.R")
library(tidyverse)
library(caret)
library(doParallel)


# Load data =================================
# X data
vst_min_max <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_min_max.txt", header=T)
Fischer_vst_min_max <- read.table("InputTables/Input_Fischer_top3000_vst_min_max.txt", header=T)
stopifnot(identical(colnames(Fischer_vst_min_max), colnames(vst_min_max)))

# y data 
# Load labels 
TCGA_BEAT_labels <- read.table("InputTables/TCGA_BEAT_labels_NPM1.txt", header=T)
TCGA_BEAT_labels <- c(TCGA_BEAT_labels$x)
Fischer_labels <- read.table("InputTables/Fischer_labels_NPM1.txt", header=T)
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


# Set Recursive Feature Elimination control ========
ctrl_rfe <- rfeControl(functions = caretFuncs, method = "repeatedcv",
                       number = 5, repeats = 3, verbose = TRUE, allowParallel = TRUE)


# Set resampling methods ==========================
trControl <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                          classProbs = TRUE, sampling = "smote",
                          summaryFunction = customSummary,
                          allowParallel = FALSE)


# Set tune grid ===================================
svmGrid <- expand.grid(cost = c(0.01, 0.1, 1, 10, 100))


# Fit models ======================================
start_time <- Sys.time()
print("Fitting model svm_vst_min_max")

svm_rfe <- rfe(x = x_train_vst_min_max, y = y_train,
               sizes = c(10, 50, 100, 300, 1000, 2000),
               method = "svmLinear2", metric = "Sens",
               tuneGrid = svmGrid, trControl = trControl, rfeControl = ctrl_rfe)

print("Time took to compute models:")
print(Sys.time() - start_time)

# saveRDS(svm_rfe, "OutputTables/NPM1_rfe_results_svm_vst_min_max.rds")
svm_rfe <- readRDS("OutputTables/NPM1_rfe_results_svm_vst_min_max.rds")


# Train final model on RFE-selected features ======================
selected_features <- predictors(svm_rfe)
x_train_rfe <- x_train_vst_min_max[, selected_features]
x_test_rfe <- x_test_vst_min_max[, selected_features]

start_time <- Sys.time()
print("Training final SVM on RFE-selected features...")
svm_tuned <- train(x = x_train_rfe, y = y_train,
                   method = "svmLinear2", trControl = trControl,
                   tuneGrid = svmGrid, metric = "Sens")

print("Final model training finished. Time took:")
print(Sys.time() - start_time)

# Save final model
# saveRDS(svm_tuned, "OutputTables/NPM1_rfe_svm_vst_min_max_model.rds")
svm_tuned <- readRDS("OutputTables/NPM1_rfe_svm_vst_min_max_model.rds")


# Test model on test data ======================
test_performance <- get_performance(svm_tuned, x_test_vst_min_max[, selected_features], y_test,
                                    classes = c("NPM1.Mutant", "NPM1.WT"),
                                    plot_title = "SVM best tuned model for NPM1 prediction")

# write_confusion_matrix(test_performance$confusion_matrix,
#                        "OutputTables/NPM1_validation_test_svm_vst_min_max_3000_CV_SMOTE_cm.txt")
# ggsave("plots/NPM1/Validation_xtest_svm_vst_min_max_3000_CV_SMOTE_performance.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)


# Test model on Fischer data ================
Fischer_performance <- get_performance(svm_tuned, Fischer_vst_min_max[, selected_features], Fischer_labels,
                                       classes = c("NPM1.Mutant", "NPM1.WT"),
                                       plot_title = "Testing SVM best tuned model for NPM1 mutation on Fischer data")
# write_confusion_matrix(Fischer_performance$confusion_matrix,
#                        "OutputTables/NPM1_validation_Fischer_svm_vst_min_max_3000_CV_SMOTE_cm.txt")
# ggsave("plots/NPM1/Validation_Fischer_svm_vst_min_max_3000_CV_SMOTE_performance.png", device = "png", width = 15, height = 15,
#        units = "cm", pointsize = 10, dpi = 500)


# Check Variable importance
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

varImp100 <- varImp_df %>% arrange(desc(importance.NPM1.Mutant)) %>% head(100) %>% select(variable)
write.table(varImp100, "OutputTables/NPM1_varImp100.txt", quote = FALSE,
            row.names = FALSE, col.names = FALSE)









