
#' @title Heatmaps of signature genes in low- and high-risk group
#'
#' @description This function compute the heatmap of genes selected by the combination of BMD-screening and network methods by dividing the patiens in high- and low-risk survival group.
#' @param x1 input training matrix \eqn{n1xp}.
#' @param x2 input training matrix \eqn{n2xp}.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param beta regression coefficients.
#' @param PI.train prognostic index on training set \eqn{T}.
#' @param PI.test prognostic index on testing set \eqn{D}.
#' @param th threshold for the \eqn{z}-score. Default is \code{th=3.5}.
#'
#' @return heatmap on training set {T} and testing set {D} divided for high-and low-groups.
#' @export
HeatmapSurvGroup <- function(x1, x2, screenVars, beta, PI.train, opt.cutoff, PI.test, th=3.5){
  
  indexScreen <- match(screenVars,colnames(x1))
  beta <- as.matrix(beta)
  
  PI.train.ordered <- as.matrix(PI.train[order(PI.train, decreasing = TRUE),])
  x1.ordered <- x1[order(PI.train, decreasing = TRUE),indexScreen]
  # dim(x1.ordered)

  # Divide patients in high and low risk group in the training set
  group.risk.train=vector()
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:nrow(PI.train.ordered))){
    if (PI.train.ordered[i] >= opt.cutoff) {group.risk.train[i]="High Risk"; gr1 = 1;} # 1=high-low
    else if (PI.train.ordered[i] < opt.cutoff) {group.risk.train[i]="Low Risk"; gr2 = 1;} # 2=low-risk
    #print(i)
  }

  # beta non zero
  ind.non.zero <- which(beta!=0)

  data.train <- data.frame(colnames(x1.ordered[,ind.non.zero]),t(x1.ordered[,ind.non.zero]))
  colnames(data.train) <- c("genes", rownames(x1.ordered[,ind.non.zero]))
  rownames(data.train) <- c(1:dim(data.train)[1])
  # dim(data.train)
  
  s.train <- length(which(group.risk.train == "High Risk"))
  
  HL.group.train <- data.frame(GroupRisk=group.risk.train)
  rownames(HL.group.train) <- rownames(PI.train.ordered)

  zscore.train <- as.matrix(t(scale(t(as.matrix(data.train[,2:dim(data.train)[2]])))))
  zscore.train[zscore.train < -th] <- -th
  zscore.train[zscore.train > th] <- th
  rownames(zscore.train) <- data.train$genes

  library(pheatmap)

  my_colour = list(
    GroupRisk = c('High Risk' = "#FC4E07", 'Low Risk' = "#00AFBB")
    )

  heatmap.train <- pheatmap(zscore.train, color = colorRampPalette(c("green", "black", "red"))(100), cluster_row = TRUE, cluster_cols = FALSE,
                                     gaps_col = s.train, annotation_col = HL.group.train, cutree_col = 2,
                                     annotation_colors = my_colour, #cutree_rows = 6,
                                     show_colnames = F, fontsize = 6.5, fontsize_row=5, filename = "zscoreTrainingData.pdf")

  clust.genes <- heatmap.train$tree_row[["order"]]
  coeff.non.zero <- beta[ind.non.zero,]
  genes.test <- as.character(coeff.non.zero[clust.genes])
  order.clust <- cutree(heatmap.train$tree_row,6)[heatmap.train$tree_row[["order"]]]

  # Testing set
  PI.test.ordered <- as.matrix(PI.test[order(PI.test, decreasing = TRUE),])
  x2.ordered <- x2[order(PI.test, decreasing = TRUE),indexScreen]

  # Divide patients in high and low risk group in the testing set
  group.risk.test=vector()
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:nrow(PI.test.ordered))){
    if (PI.test.ordered[i] >= opt.cutoff) {group.risk.test[i]="High Risk"; gr1 = 1;}
    else if (PI.test.ordered[i] < opt.cutoff) {group.risk.test[i]="Low Risk"; gr2 = 1;}
    #print(i)
  }

  data.test <- data.frame(colnames(x2.ordered[,ind.non.zero]),t(x2.ordered[,ind.non.zero]))
  colnames(data.test) <- c("genes", rownames(x2.ordered[,ind.non.zero]))
  rownames(data.test) <- c(1:dim(data.test)[1])

  HL.group.test <- data.frame(GroupRisk=group.risk.test)
  rownames(HL.group.test) <-rownames(PI.test.ordered)

  s.test <- length(which(group.risk.test == "High Risk"))
  
  zscore.test <- as.matrix(t(scale(t(as.matrix(data.test[clust.genes,2:dim(data.test)[2]])))))
  zscore.test[zscore.test < -th] <- -th
  zscore.test[zscore.test > th] <- th
  rownames(zscore.test) <- data.test[,1][clust.genes]
  
  heatmap.test <- pheatmap(zscore.test+1, color = colorRampPalette(c("green", "black", "red"))(100), cluster_row = FALSE, cluster_cols = FALSE,
                                    gaps_col = s.test, annotation_col = HL.group.test, cutree_col = 2,
                                    annotation_colors = my_colour,
                                    show_colnames = F, fontsize = 6.5, fontsize_row=5, filename = "zscoreTestingData.pdf")
}