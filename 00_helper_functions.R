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

# Min-Max scaling function
min_max <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# Convert matrices to long format with a label
make_long <- function(mat, label){
  as.data.frame(mat) %>%
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Value") %>%
    mutate(Transformation = label)
}

# Helper function to check performance
get_performance <- function(model, test_data, test_labels, plot_title="ROC Curve") {
  # Predictions
  preds <- predict(model, newdata = test_data)
  probs <- predict(model, newdata = test_data, type = "prob")
  
  # Extract positive-class probability
  prob_pos <- probs[, "cohesinAML"]
  
  # ROC and AUC
  roc_obj <- roc(response = test_labels, predictor = prob_pos,
                 levels = c("wtAML", "cohesinAML"), direction = "<")
  auc_value <- auc(roc_obj)
  roc_df <- data.frame(fpr = 1 - roc_obj$specificities,
                       tpr = roc_obj$sensitivities)
  
  # ROC curve
  roc_plot <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(size = 1.2, color = "blue") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = plot_title, x = "False Positive Rate", y = "True Positive Rate") +
    annotate("text", x = 0.9, y = 0.05, label = paste0("AUC = ", round(auc_value, 3)),
             size = 5, color = "black") +
    theme_minimal()
  print(roc_plot)
  
  # Confusion matrix
  cm <- confusionMatrix(preds, test_labels)
  print(cm)
  
  # Return everything
  return(list(
    predicted_classes = preds,
    predicted_probabilities = prob_pos,
    auc = auc_value,
    roc_plot = roc_plot,
    confusion_matrix = cm))
}
