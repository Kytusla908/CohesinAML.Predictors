library(ggplot2)
library(tidyverse)
library(pROC)

plotPCA.DESeqTransform = function(object, intgroup="condition", ntop=500, returnData=F)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  Samples <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=Samples, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="Samples"))  + geom_point(size = 1.5) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) + 
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() + theme_bw()
}


# Min-Max scaling function ============================
min_max <- function(x) {
  rng <- range(x)
  if (rng[1] == rng[2]) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}


# Convert matrices to long format with a label ============================
make_long <- function(mat, label){
  as.data.frame(mat) %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Value") %>%
    mutate(Transformation = label)
}


# Helper function to check performance ============================
get_performance <- function(model, test_data, test_labels,
                            classes=c("positive_class", "negative_class"),
                            plot_title="ROC Curve", plot_subtitle = NULL) {
  
  # Check train and test contain same variables
  missing_vars <- setdiff(colnames(model$trainingData)[-1], colnames(test_data))
  missing_vars <- setdiff(missing_vars, ".outcome")
  if (length(missing_vars) > 0) {
    stop("Test data is missing predictors: ", paste(missing_vars, collapse = ", "))
  }
  
  # Predictions
  preds <- predict(model, newdata = test_data)
  probs <- predict(model, newdata = test_data, type = "prob")
  
  # Extract positive-class probability
  prob_pos <- probs[, classes[1]]
  
  # ROC and AUC
  roc_obj <- pROC::roc(response = test_labels, predictor = prob_pos,
                       levels = c(classes[2], classes[1]), direction = "<")
  auc_value <- pROC::auc(roc_obj)
  
  # ROC dataframe
  roc_df <- data.frame(fpr = 1 - roc_obj$specificities,
                       tpr = roc_obj$sensitivities,
                       thresholds = roc_obj$thresholds)
  
  # ROC curve
  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(size = 1.2, color = "blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = plot_title, subtitle = plot_subtitle,
         x = "False Positive Rate", y = "True Positive Rate") +
    annotate("text", x = 0.9, y = 0.05, label = paste0("AUC = ", round(auc_value, 3)),
             size = 5, color = "black") +
    theme_minimal() + 
    theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))
  print(roc_plot)
  
  # Confusion matrix
  cm <- confusionMatrix(preds, test_labels)
  print(cm)
  
  # Return everything
  return(list(
    predicted_classes = preds,
    predicted_probabilities = prob_pos,
    auc = auc_value,
    roc_df = roc_df,
    roc_obj = roc_obj,
    roc_plot = roc_plot,
    confusion_matrix = cm))
}


# Helper function to write out Confusion Matrix =======================
write_confusion_matrix <- function(cm, file="confusion_matrix.txt") {
  txt <- capture.output(print(cm), type = "output")
  con <- file(file, open = "wb")
  writeBin(paste0(paste(txt, collapse = "\n"), "\n"), con)
  close(con)
}


# Helper function to get ensemble predictions =======================
ensemble_predict <- function(all_models, selected_models, test_data_map,
                             positive_class = "cohesinAML", negative_class = "wtAML",
                             threshold = 0.5, weights = NULL) {
  
  # checks 
  if (!all(selected_models %in% names(all_models))) {
    stop("Some selected_models are not present in all_models")
  }
  
  if (!all(selected_models %in% names(test_data_map))) {
    stop("Some selected_models are not present in test_data_map")
  }
  if (!is.null(weights)) {
    if (!all(names(weights) %in% selected_models)) {
      stop("All weight names must match selected_models")
    }
    if (any(weights < 0)) {
      stop("All weights must be non-negative")
    }
    # all weights must sum 1
    weights <- weights / sum(weights)
  }
  
  # select models
  ensemble_models <- all_models[selected_models]
  
  # predict probabilities
  ensemble_preds <- lapply(selected_models, function(m) {
    preds <- predict(ensemble_models[[m]],
                     newdata = test_data_map[[m]],
                     type = "prob")
    
    # handle matrix / data.frame output
    if (is.matrix(preds) || is.data.frame(preds)) {
      preds[, positive_class]
    } else { preds }
    })
  
  names(ensemble_preds) <- selected_models
  
  # average probabilities
  preds_df <- as.data.frame(ensemble_preds)
  if (is.null(weights)) {
    avg_pred <- rowMeans(preds_df)
  } else {
    avg_pred <- rowSums(sweep(preds_df, 2, weights[names(preds_df)], `*`))
  }
  
  # final class
  ensemble_class <- ifelse(avg_pred > threshold,
                           positive_class,
                           negative_class)
  
  return(list(probabilities = avg_pred,
              classes = as.factor(ensemble_class),
              per_model_probabilities = preds_df))
}


# Helper function to register execution times for each model =================
log_time <- function(model_name, start_time, timing_log = NULL) {
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  cat(sprintf("[%s] %s | elapsed: %.2f sec\n",
              format(end_time, "%Y-%m-%d %H:%M:%S"), model_name, elapsed))
  
  if (!is.null(timing_log)){
    timing_log <- rbind(timing_log,
                        data.frame(model = model_name,
                                   start = start_time,
                                   end = end_time,
                                   elapsed = elapsed,
                                   stringsAsFactors = FALSE))
    return(timing_log)
  }
  invisible(NULL)
}

