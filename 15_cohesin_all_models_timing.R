load("03_CV_SMOTE_3000.RData")
source("00_helper_functions.R")
library(tidyverse)
library(caret)
library(lme4)
library(lmerTest)


# Data partition (from here, run on cluster) ================================
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


# Time record table =====================
timing_log <- data.frame(model = character(),
                         start = as.POSIXct(character()),
                         end = as.POSIXct(character()),
                         elapsed_sec = numeric(),
                         stringsAsFactors = FALSE)


# Try on some models ============================
# Set train control obj
trControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, sampling = "smote")

# k-NN
k_val <- round(sqrt(length(y_train)), 0)

t0 <- Sys.time()
kNN_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="knn",
                         trControl=trControl, tuneGrid=data.frame(k=k_val))
timing_log <- log_time("kNN_vst_min_max", t0, timing_log)

t0 <- Sys.time()
kNN_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="knn",
                        trControl=trControl, tuneGrid=data.frame(k=k_val))
timing_log <- log_time("kNN_vst_zScore", t0, timing_log)

t0 <- Sys.time()
kNN_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="knn",
                          trControl=trControl, tuneGrid=data.frame(k=k_val))
timing_log <- log_time("kNN_rlog_min_max", t0, timing_log)

t0 <- Sys.time()
kNN_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="knn",
                         trControl=trControl, tuneGrid=data.frame(k=k_val))
timing_log <- log_time("kNN_rlog_zScore", t0, timing_log)

# Naive Bayes
nb_grid <- data.frame(usekernel=TRUE, laplace=0, adjust=1)

t0 <- Sys.time()
nb_raw <- train(x=x_train_raw, y=y_train, method="naive_bayes",
                trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_raw", t0, timing_log)

t0 <- Sys.time()
nb_vst <- train(x=x_train_vst, y=y_train, method="naive_bayes",
                trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_vst", t0, timing_log)

t0 <- Sys.time()
nb_rlog <- train(x=x_train_rlog, y=y_train, method="naive_bayes",
                 trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_rlog", t0, timing_log)

t0 <- Sys.time()
nb_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="naive_bayes",
                        trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_vst_min_max", t0, timing_log)

t0 <- Sys.time()
nb_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="naive_bayes",
                       trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_vst_zScore", t0, timing_log)

t0 <- Sys.time()
nb_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="naive_bayes",
                   trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_vst_QN", t0, timing_log)

t0 <- Sys.time()
nb_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="naive_bayes",
                         trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_rlog_min_max", t0, timing_log)

t0 <- Sys.time()
nb_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="naive_bayes",
                        trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_rlog_zScore", t0, timing_log)

t0 <- Sys.time()
nb_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="naive_bayes",
                    trControl=trControl, tuneGrid=nb_grid)
timing_log <- log_time("nb_rlog_QN", t0, timing_log)

# SVM
svm_grid <- data.frame(cost=1)

t0 <- Sys.time()
svm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="svmLinear2",
                         trControl=trControl, tuneGrid=svm_grid)
timing_log <- log_time("svm_vst_min_max", t0, timing_log)

t0 <- Sys.time()
svm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="svmLinear2",
                        trControl=trControl, tuneGrid=svm_grid)
timing_log <- log_time("svm_vst_zScore", t0, timing_log)

t0 <- Sys.time()
svm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="svmLinear2",
                          trControl=trControl, tuneGrid=svm_grid)
timing_log <- log_time("svm_rlog_min_max", t0, timing_log)

t0 <- Sys.time()
svm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="svmLinear2",
                         trControl=trControl, tuneGrid=svm_grid)
timing_log <- log_time("svm_rlog_zScore", t0, timing_log)

# Logistic Regression
t0 <- Sys.time()
glm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="glm",
                         family="binomial", trControl=trControl)
timing_log <- log_time("glm_vst_min_max", t0, timing_log)

t0 <- Sys.time()
glm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="glm",
                        family="binomial", trControl=trControl)
timing_log <- log_time("glm_vst_zScore", t0, timing_log)

t0 <- Sys.time()
glm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="glm",
                          family="binomial", trControl=trControl)
timing_log <- log_time("glm_rlog_min_max", t0, timing_log)

t0 <- Sys.time()
glm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="glm",
                         family="binomial", trControl=trControl)
timing_log <- log_time("glm_rlog_zScore", t0, timing_log)

# Random Forest
mtry_val <- floor(sqrt(ncol(x_train_raw)))

t0 <- Sys.time()
rf_raw <- train(x=x_train_raw, y=y_train, method="rf",
                trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_raw", t0, timing_log)

t0 <- Sys.time()
rf_vst <- train(x=x_train_vst, y=y_train, method="rf",
                trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_vst", t0, timing_log)

t0 <- Sys.time()
rf_rlog <- train(x=x_train_rlog, y=y_train, method="rf",
                 trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_rlog", t0, timing_log)

t0 <- Sys.time()
rf_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="rf",
                        trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_vst_min_max", t0, timing_log)

t0 <- Sys.time()
rf_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="rf",
                       trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_vst_zScore", t0, timing_log)

t0 <- Sys.time()
rf_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="rf",
                   trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_vst_QN", t0, timing_log)

t0 <- Sys.time()
rf_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="rf",
                         trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_rlog_min_max", t0, timing_log)

t0 <- Sys.time()
rf_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="rf",
                        trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_rlog_zScore", t0, timing_log)

t0 <- Sys.time()
rf_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="rf",
                    trControl=trControl, tuneGrid=data.frame(mtry=mtry_val))
timing_log <- log_time("rf_rlog_QN", t0, timing_log)

# Boosted Decision Trees
gbm_grid <- data.frame(n.trees=100, interaction.depth=1,
                       shrinkage=0.1, n.minobsinnode=10)

t0 <- Sys.time()
gbm_raw <- train(x=x_train_raw, y=y_train, method="gbm",
                 trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_raw", t0, timing_log)

t0 <- Sys.time()
gbm_vst <- train(x=x_train_vst, y=y_train, method="gbm",
                 trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_vst", t0, timing_log)

t0 <- Sys.time()
gbm_rlog <- train(x=x_train_rlog, y=y_train, method="gbm",
                  trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_rlog", t0, timing_log)

t0 <- Sys.time()
gbm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="gbm",
                         trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_vst_min_max", t0, timing_log)

t0 <- Sys.time()
gbm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="gbm",
                        trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_vst_zScore", t0, timing_log)

t0 <- Sys.time()
gbm_vst_QN <- train(x=x_train_vst_QN, y=y_train, method="gbm",
                    trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_vst_QN", t0, timing_log)

t0 <- Sys.time()
gbm_rlog_min_max <- train(x=x_train_rlog_min_max, y=y_train, method="gbm",
                          trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_rlog_min_max", t0, timing_log)

t0 <- Sys.time()
gbm_rlog_zScore <- train(x=x_train_rlog_zScore, y=y_train, method="gbm",
                         trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_rlog_zScore", t0, timing_log)

t0 <- Sys.time()
gbm_rlog_QN <- train(x=x_train_rlog_QN, y=y_train, method="gbm",
                     trControl=trControl, tuneGrid=gbm_grid, verbose=FALSE)
timing_log <- log_time("gbm_rlog_QN", t0, timing_log)


# Final timing table ===================
timing_log <- timing_log[order(-timing_log$elapsed), ]
print(timing_log)
write.table(timing_log, "OutputTables/cohesin_basic_all_model_training_times.txt", col.names=TRUE, sep="\t")

timing_log <- read.table("OutputTables/cohesin_basic_all_model_training_times.txt", sep="\t", row.names = 1)
parts <- strsplit(timing_log$model, "_")

timing_log$model_type <- as.factor(sapply(parts, `[`, 1))
timing_log$stabilizing_transformation <- factor(sapply(parts, `[`, 2), levels = c("raw", "vst", "rlog"))
timing_log$transformation <- sapply(parts, function(x) {
  if (length(x) >= 3) {
    paste(x[3:length(x)], collapse = "_")
  } else {
    "none"
  }
})
timing_log$transformation <- factor(timing_log$transformation,
                                                levels = c("none", "min_max", "zScore", "QN"))

order <- c("model", "model_type", "stabilizing_transformation","transformation", "start", "end", "elapsed")
timing_log <- timing_log[, order]
timing_log$log_elapsed <- log10(timing_log$elapsed)


# Test vst vs rlog vs transformations =============================
lmm <- lmer(log_elapsed ~ stabilizing_transformation + transformation + (1 | model_type),
            data = timing_log)
summary(lmm)


# vst vs rlog FC
10^(lmm@beta[3] - lmm@beta[2])


# Plot distribution of timings in stabilizing transformations =================
ggplot(timing_log,
       aes(x = stabilizing_transformation, y = log_elapsed, fill = stabilizing_transformation)) +
  geom_boxplot() +
  labs(x = "Stabilizing transformation",
       y = "Training time (log10(s))") + guides(fill = "none")
# ggsave("plots/timing_distribution_stabilizing_transformation_boxplot.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)


# Plot distribution of timings in transformations =================
ggplot(timing_log,
       aes(x = transformation, y = log_elapsed, fill = transformation)) +
  geom_boxplot() +
  labs(x = "Data Transformation",
       y = "Training time (log10(s))") + guides(fill = "none")
# ggsave("plots/timing_distribution_data_transformations_boxplot.png", device = "png", width = 12, height = 12,
#        units = "cm", pointsize = 10, dpi = 500)

































