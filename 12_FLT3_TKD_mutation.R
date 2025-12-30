source("00_helper_functions.R")
library(tidyverse)
library(patchwork)
library(caret)
library(PRROC)

dir.create("plots/FLT3_TKD", recursive = TRUE, showWarnings = FALSE)


# Load data =================================
# TCGA-BEAT data
TCGA_BEAT_mut_table <- read.table("InputTables/TCGA_BEAT_ALL_mutation_table.txt", sep="\t", header=T)
TCGA_BEAT_mut_table <- subset(TCGA_BEAT_mut_table, TCGA_BEAT_mut_table$Disease == "AML")
vst_transform <- t(read.table("InputTables/Input_NormalizedCounts_TCGA-BEAT_filterByExpr_vst_varFiltered.txt", header=T))
vst_min_max <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_min_max.txt", header=T)
vst_zScore <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_zScore.txt", header=T)
vst_QN <- read.table("InputTables/Input_TCGA-BEAT_top3000_vst_QN.txt", header=T)

# Check samples are in same order
table(row.names(vst_transform) == row.names(TCGA_BEAT_mut_table))
table(row.names(vst_min_max) == row.names(TCGA_BEAT_mut_table))
table(row.names(vst_zScore) == row.names(TCGA_BEAT_mut_table))
table(row.names(vst_QN) == row.names(TCGA_BEAT_mut_table))


# Read mutation type ==========================
# TCGA
TCGA_FLT3_mut_type <- read.delim("InputTables/TCGA_FLT3_mutations.tsv", header = TRUE,
                            stringsAsFactors = FALSE, check.names = FALSE)
TCGA_FLT3_mut_type$`Sample ID` <- gsub("-", ".", TCGA_FLT3_mut_type$`Sample ID`)
TCGA_FLT3_mut_type <- TCGA_FLT3_mut_type %>% select(`Sample ID`, `Protein Change`, `Mutation Type`, `Variant Type`)
coords <- regmatches(TCGA_FLT3_mut_type$`Protein Change`,
                     gregexpr("[0-9]+", TCGA_FLT3_mut_type$`Protein Change`))
TCGA_FLT3_mut_type$coord_start <- sapply(coords, function(x) as.integer(x[1]))
TCGA_FLT3_mut_type$coord_end <- sapply(coords, function(x) {
  if (length(x) >= 2) as.integer(x[2]) else NA_integer_
})
TCGA_FLT3_mut_type$FLT3_TKD <- ifelse(
  TCGA_FLT3_mut_type$`Mutation Type` == "Missense_Mutation" &
    TCGA_FLT3_mut_type$coord_start >= 835 & TCGA_FLT3_mut_type$coord_start <= 842,
  "Mutant", "WT")

# BEAT
BEAT_FLT3_mut_type <- read.delim("InputTables/BEAT_FLT3_mutations.tsv", header = TRUE,
                                 stringsAsFactors = FALSE, check.names = FALSE)
BEAT_FLT3_mut_type$`Sample ID` <- paste0(sub(".*_", "", BEAT_FLT3_mut_type$`Sample ID`), "R")
BEAT_FLT3_mut_type <- BEAT_FLT3_mut_type %>% select(`Sample ID`, `Protein Change`, `Mutation Type`, `Variant Type`)
coords <- regmatches(BEAT_FLT3_mut_type$`Protein Change`,
                     gregexpr("[0-9]+", BEAT_FLT3_mut_type$`Protein Change`))
BEAT_FLT3_mut_type$coord_start <- sapply(coords, function(x) as.integer(x[1]))
BEAT_FLT3_mut_type$coord_end <- sapply(coords, function(x) {
  if (length(x) >= 2) as.integer(x[2]) else NA_integer_
})
BEAT_FLT3_mut_type$FLT3_TKD <- ifelse(
  BEAT_FLT3_mut_type$`Mutation Type` == "Missense_Mutation" &
    BEAT_FLT3_mut_type$coord_start >= 835 & BEAT_FLT3_mut_type$coord_start <= 842,
  "Mutant", "WT")


# Ensure no repeated samples (some has both ITD and TKD) ==============
TCGA_BEAT_mut_table$Sample_ID <- row.names(TCGA_BEAT_mut_table)
FLT3_all <- rbind(TCGA_FLT3_mut_type[, c("Sample ID", "FLT3_TKD")],
                  BEAT_FLT3_mut_type[, c("Sample ID", "FLT3_TKD")])
FLT3_all <- FLT3_all %>%
  group_by(`Sample ID`) %>%
  summarise(FLT3_TKD = ifelse(any(FLT3_TKD == "Mutant"), "Mutant", "WT"),
            .groups = "drop")


# Append to full mutation table =========================
TCGA_BEAT_mut_table <- TCGA_BEAT_mut_table %>%
  left_join(FLT3_all, by = c("Sample_ID" = "Sample ID"))
TCGA_BEAT_mut_table$FLT3_TKD[is.na(TCGA_BEAT_mut_table$FLT3_TKD)] <- "WT"
row.names(TCGA_BEAT_mut_table) <- TCGA_BEAT_mut_table$Sample_ID
TCGA_BEAT_mut_table$Sample_ID <- NULL

# Get labels for FLT3 mutations =================
TCGA_BEAT_FLT3_TKD_labels <- paste0("FLT3_TKD.", TCGA_BEAT_mut_table$FLT3_TKD)
TCGA_BEAT_FLT3_TKD_labels <- as.factor(TCGA_BEAT_FLT3_TKD_labels)
table(TCGA_BEAT_FLT3_TKD_labels)

# Plot mutants proportions
flt3_tab <- as.data.frame(table(TCGA_BEAT_FLT3_TKD_labels))
colnames(flt3_tab) <- c("flt3", "Count")
ggplot(flt3_tab, aes(x = "", y = Count, fill = flt3)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black") +
  coord_polar(theta = "y") +
  labs(title="FLT3 mutant proportions in TCGA-BEAT samples") +
  theme_void() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
# ggsave("plots/FLT3_TKD/FLT3_TKD_mutants_proportions_TCGA_BEAT.png", device="png",
#        units="cm", dpi=500, width=15, height=12)
# write.table(TCGA_BEAT_FLT3_TKD_labels, "InputTables/TCGA_BEAT_labels_FLT3_TKD.txt",
#             row.names = F, quote = F, sep = "\t")


# Data partition ===============================
# Load partition indexes
train_index <- scan("InputTables/Input_train_indexes.txt", sep="\n")

# Splitting Labels
y_train <- as.factor(TCGA_BEAT_FLT3_TKD_labels[train_index])
y_test <- as.factor(TCGA_BEAT_FLT3_TKD_labels[-train_index])
train_df <- data.frame(Class=y_train) %>% group_by(Class) %>%
  summarise(Count = n())
test_df <- data.frame(Class=y_test) %>% group_by(Class) %>%
  summarise(Count = n())
# Plot proportions
p_train <- ggplot(train_df, aes(x = "", y = Count, fill = Class)) +
  geom_col(width = 1) + coord_polar(theta = "y") + theme_void() + 
  labs(title = "Train Split") +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
p_test <- ggplot(test_df, aes(x = "", y = Count, fill = Class)) +
  geom_col(width = 1) + coord_polar(theta = "y") + theme_void() +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black") +
  labs(title = "Test Split")+
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
p_train + p_test + plot_layout(guides = "collect") & theme(legend.position = "right")
# ggsave("plots/FLT3_TKD/train_test_split_proportions.png", device = "png", width = 15, height = 8,
#        units = "cm", pointsize = 10, dpi = 500)

# Subset train data
x_train_vst <- vst_transform[train_index, ]
x_train_vst_min_max <- vst_min_max[train_index, ]
x_train_vst_zScore <- vst_zScore[train_index, ]
x_train_vst_QN <- vst_QN[train_index, ]

# Subset test data
x_test_vst <- vst_transform[-train_index, ]
x_test_vst_min_max <- vst_min_max[-train_index, ]
x_test_vst_zScore <- vst_zScore[-train_index, ]
x_test_vst_QN <- vst_QN[-train_index, ]


# Try on some models ============================
set.seed(12345)

# Set ROSE resampling and CV
trControl <- trainControl(method = "cv", number = 5, classProbs = TRUE, sampling = "smote")

# k-NN
start <- Sys.time()

k_val <- round(sqrt(length(y_train)), 0)
tuneGrid <- data.frame(k = k_val)
print("Fitting model kNN_vst_min_max")
kNN_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="knn", trControl=trControl,
                         tuneGrid = tuneGrid)
print("Fitting model kNN_vst_zScore")
kNN_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="knn", trControl=trControl,
                        tuneGrid = tuneGrid)
stop <- Sys.time()
cat("\nkNN Execution time: ", stop-start)

# Naive Bayes
start <- Sys.time()

tuneGrid <- data.frame(usekernel = TRUE, laplace = 0, adjust = 1)
print("Fitting model nb_vst")
nb_vst <- train(x=x_train_vst, y=y_train, method="naive_bayes", trControl=trControl,
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
stop <- Sys.time()
cat("\nNaive Bayes Execution time: ", stop-start)

# SVM
start <- Sys.time()

tuneGrid <- data.frame(cost = 1)
print("Fitting model svm_vst_min_max")
svm_vst_min_max <- train(x=x_train_vst_min_max, y=y_train, method="svmLinear2", trControl=trControl,
                         tuneGrid = tuneGrid)
print("Fitting model svm_vst_zScore")
svm_vst_zScore <- train(x=x_train_vst_zScore, y=y_train, method="svmLinear2", trControl=trControl,
                        tuneGrid = tuneGrid)
stop <- Sys.time()
cat("\nSVM Execution time: ", stop-start)

# Logistic Regression
start <- Sys.time()

print("Fitting model glm_vst")
glm_vst <- train(x=x_train_vst, y=y_train, method="glm",
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
stop <- Sys.time()
cat("\nGLM Execution time: ", stop-start)

# Random Forest
start <- Sys.time()

mtry_val <- floor(sqrt(ncol(x_train_vst)))
tuneGrid <- data.frame(mtry = mtry_val)

print("Fitting model rf_vst")
rf_vst <- train(x=x_train_vst, y=y_train, method="rf", trControl=trControl,
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
stop <- Sys.time()
cat("\nRandom Forest Execution time: ", stop-start)

# Boosted Decision Trees
start <- Sys.time()

tuneGrid <- data.frame(n.trees = 100, interaction.depth = 1, shrinkage = 0.1,
                       n.minobsinnode = 10)

print("Fitting model gbm_vst")
gbm_vst <- train(x=x_train_vst, y=y_train, method="gbm", trControl=trControl,
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
stop <- Sys.time()
cat("\nGBM Execution time: ", stop-start)

# Save models =================================
# List of all models
all_models <- list(
  kNN_vst_min_max = kNN_vst_min_max,
  kNN_vst_zScore = kNN_vst_zScore,
  
  nb_vst = nb_vst,
  nb_vst_min_max = nb_vst_min_max,
  nb_vst_zScore = nb_vst_zScore,
  nb_vst_QN = nb_vst_QN,
  
  svm_vst_min_max = svm_vst_min_max,
  svm_vst_zScore = svm_vst_zScore,
  
  glm_vst = glm_vst,
  glm_vst_min_max = glm_vst_min_max,
  glm_vst_zScore = glm_vst_zScore,
  glm_vst_QN = glm_vst_QN,
  
  rf_vst = rf_vst,
  rf_vst_min_max = rf_vst_min_max,
  rf_vst_zScore = rf_vst_zScore,
  rf_vst_QN = rf_vst_QN,
  
  gbm_vst = gbm_vst,
  gbm_vst_min_max = gbm_vst_min_max,
  gbm_vst_zScore = gbm_vst_zScore,
  gbm_vst_QN = gbm_vst_QN
)

# Save each model as an RDS file
# saveRDS(all_models, file = "OutputTables/FLT3_TKD_basic_all_models_CV_SMOTE.rds")
all_models <- readRDS("OutputTables/FLT3_TKD_basic_all_models_CV_SMOTE.rds")
list2env(all_models, envir = .GlobalEnv)


# Map test_data ================================
test_data_map <- list(
  kNN_vst_min_max = x_test_vst_min_max, kNN_vst_zScore = x_test_vst_zScore,
  
  nb_vst = x_test_vst, nb_vst_min_max = x_test_vst_min_max,
  nb_vst_zScore = x_test_vst_zScore, nb_vst_QN = x_test_vst_QN,
  
  svm_vst_min_max = x_test_vst_min_max, svm_vst_zScore = x_test_vst_zScore,
  
  glm_vst = x_test_vst, glm_vst_min_max = x_test_vst_min_max,
  glm_vst_zScore = x_test_vst_zScore, glm_vst_QN = x_test_vst_QN,
  
  rf_vst = x_test_vst, rf_vst_min_max = x_test_vst_min_max, 
  rf_vst_zScore = x_test_vst_zScore, rf_vst_QN = x_test_vst_QN,
  
  gbm_vst = x_test_vst, gbm_vst_min_max = x_test_vst_min_max,
  gbm_vst_zScore = x_test_vst_zScore, gbm_vst_QN = x_test_vst_QN
)


# Measure performance ==========================
performance_df <- data.frame(model = character(), sensitivity = numeric(),
                             specificity = numeric(), precision = numeric(),
                             accuracy = numeric(), kappa = numeric(),
                             auc = numeric(), pr_auc = numeric(),
                             stringsAsFactors = FALSE)
roc_df_list <- list()

pb <- txtProgressBar(min = 0, max = length(names(all_models)), style = 3)
i <- 0
for (model_name in names(all_models)){
  i <- i + 1
  setTxtProgressBar(pb, i)
  
  model <- all_models[[model_name]]
  test_data <- test_data_map[[model_name]]
  performance <- get_performance(model, test_data, y_test,
                                 classes = c("FLT3_TKD.Mutant", "FLT3_TKD.WT"),
                                 plot_title = paste0(model_name, " model ROC curve"))
  # Get PRROC
  pr <- pr.curve(
    scores.class0 = probs[y_test == "FLT3_TKD.Mutant"],
    scores.class1 = probs[y_test == "FLT3_TKD.WT"],
    curve = TRUE)
  pr_auc <- pr$auc.integral
  
  # Save ROC and PR-ROC curve info
  roc_df_list[[model_name]] <- performance$roc_df
  prroc_df_list[[model_name]] <- performance$roc_df
  
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
                                                     auc = auc, pr_auc = pr_auc, stringsAsFactors = FALSE))
}
close(pb)

# saveRDS(roc_df_list, "OutputTables/FLT3_TKD_basic_models_ROC_data.rds")
# write.table(performance_df, "OutputTables/FLT3_TKD_basic_all_models_CV_SMOTE_performance.txt",
#             col.names=TRUE, sep="\t")

roc_df_list <- readRDS("OutputTables/FLT3_TKD_basic_models_ROC_data.rds")
performance_df <- read.table("OutputTables/FLT3_TKD_basic_all_models_CV_SMOTE_performance.txt", header=T)
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
ggsave("plots/FLT3_TKD/FLT3_basic_model_performance_CV_SMOTE_3000_sensitivity.png", device = "png", width = 15, height = 15,
       units = "cm", pointsize = 10, dpi = 500)

# Plot Specificity
ggplot(performance_df, aes(x = transformation, y = specificity, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Specificity", color = "Model")
ggsave("plots/FLT3_TKD/FLT3_basic_model_performance_CV_SMOTE_3000_specificity.png", device = "png", width = 15, height = 15,
       units = "cm", pointsize = 10, dpi = 500)

# Plot Precision
ggplot(performance_df, aes(x = transformation, y = precision, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Precision", color = "Model")
ggsave("plots/FLT3_TKD/FLT3_basic_model_performance_CV_SMOTE_3000_precision.png", device = "png", width = 15, height = 15,
       units = "cm", pointsize = 10, dpi = 500)

# Plot kappa
ggplot(performance_df, aes(x = transformation, y = kappa, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "Kappa", color = "Model")
ggsave("plots/FLT3_TKD/FLT3_basic_model_performance_CV_SMOTE_3000_kappa.png", device = "png", width = 15, height = 15,
       units = "cm", pointsize = 10, dpi = 500)

# Plot AUC
ggplot(performance_df, aes(x = transformation, y = auc, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Performance Across Transformations",
       x = "Transformation", y = "AUC", color = "Model")
ggsave("plots/FLT3_TKD/FLT3_basic_model_performance_CV_SMOTE_3000_AUC.png", device = "png", width = 15, height = 15,
       units = "cm", pointsize = 10, dpi = 500)


# Save best performing model =====================================









