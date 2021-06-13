
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
#' @return Heatmap on training set \eqn{T} and testing set \eqn{D} divided in high-and low-risk survival groups.
#' @export
HeatmapSurvGroup <- function(x1, x2, screenVars, beta, PI.train, opt.cutoff, PI.test, th=3.5){
  
  library(ggplot2)
  library(ggplotify)
  library(pheatmap)
  library(patchwork)
  
  if(is.list(screenVars)==TRUE){
    index <- match(unlist(screenVars),colnames(x1))
  } else {index <- match(screenVars,colnames(x1))}
  
  if(length(which(is.na(index)))==0){indexScreen=index
  } else {indexScreen <- index[-which(is.na(index))]}
  beta <- as.matrix(beta)
  
  PI.train.ordered <- as.numeric(PI.train[order(PI.train, decreasing = TRUE)])
  x1.ordered <- x1[order(PI.train, decreasing = TRUE),indexScreen]
  # dim(x1.ordered)

  # Divide patients in high and low risk group in the training set
  group.risk.train=vector()
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:length(PI.train.ordered))){
    if (PI.train.ordered[i] >= opt.cutoff) {group.risk.train[i]="High Risk"; gr1 = 1;} # 1=high-low
    else if (PI.train.ordered[i] < opt.cutoff) {group.risk.train[i]="Low Risk"; gr2 = 1;} # 2=low-risk
    #print(i)
  }

  # beta non zero
  ind.non.zero1 <- which(rownames(beta)[which(beta!=0)]%in%colnames(x1)==TRUE)

  data.train <- data.frame(colnames(x1.ordered[,ind.non.zero1]),t(x1.ordered[,ind.non.zero1]))
  colnames(data.train) <- c("genes", rownames(x1.ordered[,ind.non.zero1]))
  rownames(data.train) <- c(1:dim(data.train)[1])
  # dim(data.train)
  
  s.train <- length(which(group.risk.train == "High Risk"))
  
  HL.group.train <- data.frame(GroupRisk=group.risk.train)
  colnames(HL.group.train) <- "Group Risk"
  rownames(HL.group.train) <- rownames(x1.ordered)

  zscore.train <- as.matrix(t(scale(t(as.matrix(data.train[,2:dim(data.train)[2]])))))
  zscore.train[zscore.train < -th] <- -th
  zscore.train[zscore.train > th] <- th
  rownames(zscore.train) <- data.train$genes
 
  my_colour = list(
    'Group Risk' = c('High Risk' = "#FC4E07", 'Low Risk' = "#00AFBB")
    )

  heatmap.train <- pheatmap(zscore.train, color = colorRampPalette(c("green", "black", "red"))(100), cluster_row = TRUE, cluster_cols = FALSE,
                                     gaps_col = s.train, annotation_col = HL.group.train, cutree_col = 2,
                                     annotation_colors = my_colour,
                                     show_colnames = F, fontsize = 6.5, fontsize_row=5)

  clust.genes <- heatmap.train$tree_row[["order"]]
  coeff.non.zero <- beta[ind.non.zero1,]
  genes.test <- as.character(coeff.non.zero[clust.genes])
  order.clust <- cutree(heatmap.train$tree_row,6)[heatmap.train$tree_row[["order"]]]

  # Testing set
  PI.test.ordered <- as.numeric(PI.test[order(PI.test, decreasing = TRUE)])
  x2.ordered <- x2[order(PI.test, decreasing = TRUE),indexScreen]

  # Divide patients in high and low risk group in the testing set
  group.risk.test=vector()
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:length(PI.test.ordered))){
    if (PI.test.ordered[i] >= opt.cutoff) {group.risk.test[i]="High Risk"; gr1 = 1;}
    else if (PI.test.ordered[i] < opt.cutoff) {group.risk.test[i]="Low Risk"; gr2 = 1;}
    #print(i)
  }
  
  ind.non.zero2 <- which(rownames(beta)[which(beta!=0)]%in%colnames(x2)==TRUE)
  data.test <- data.frame(colnames(x2.ordered[,ind.non.zero2]),t(x2.ordered[,ind.non.zero2]))
  colnames(data.test) <- c("genes", rownames(x2.ordered[,ind.non.zero2]))
  rownames(data.test) <- c(1:dim(data.test)[1])

  HL.group.test <- data.frame(GroupRisk=group.risk.test)
  colnames(HL.group.test) <- "Group Risk"
  rownames(HL.group.test) <-rownames(x2.ordered)

  s.test <- length(which(group.risk.test == "High Risk"))
  
  zscore.test <- as.matrix(t(scale(t(as.matrix(data.test[clust.genes,2:dim(data.test)[2]])))))
  zscore.test[zscore.test < -th] <- -th
  zscore.test[zscore.test > th] <- th
  rownames(zscore.test) <- data.test[,1][clust.genes]
  
  heatmap.test <- pheatmap(zscore.test+1, color = colorRampPalette(c("green", "black", "red"))(100), cluster_row = FALSE, cluster_cols = FALSE,
                                    gaps_col = s.test, annotation_col = HL.group.test, cutree_col = 2,
                                    annotation_colors = my_colour,
                                    show_colnames = F, fontsize = 6.5, fontsize_row=5)
  
  # h <- ggpubr::ggarrange(heatmap.train[[4]], heatmap.test[[4]], ncol = 2, nrow = 1)
  # print(h)
  # ggsave("heatmap.pdf", plot = h)
  
  # save plots
  p1 <- as.ggplot(heatmap.train)
  p2 <- as.ggplot(heatmap.test)
  # use patchwork to arrange them together
  print(p1 + p2)
}