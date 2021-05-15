
#' @title Normalization between training and testing sets
#'
#' @description This function performs the normalization between training and testing set. We normalize the data so that the two sets have the same shape (or distribution).
#' @param x1 a matrix of intensities, rows are genes and columns are samples \eqn{pxn1}. 
#' @param x2 a matrix of intensities, rows are genes and columns are samples \eqn{pxn2}.
#' @param method character string specifying the \code{quantile} normalization method.
#' @param plot logical flag for plot the density distribution. Default is \code{plot=FALSE}.
#'
#' @return The following objects are returned:
#' \item{x1.norm}{a normalized matrix \eqn{n1xp}}
#' \item{x2.norm}{a normalized matrix \eqn{n2xp}}
#' 
#' @details This function uses a quantile normalization in order to make the distributions of training and testing set the same across samples.
#' The normalization approach used in our package consist of adding by column each sample of the test set to the train set and normalize the new dataset. Then, we take the test column normalized and build the normalized testing set column by column.
#' This improves the performance and stability of the models and make the two datasets comparable between them. 
#' 
#' @export
NormalizeBetweenData <- function(x1, x2, norm.method="quantile", plot=FALSE){
  
  x2.comb <- NULL
  for(i in 1:dim(x2)[2]){
    # i=1
    comb <- cbind(x1,x2[,i])
    library(preprocessCore)
    normalized <- normalize.quantiles(as.matrix(comb),copy=FALSE)
    x2.col <- normalized[,dim(normalized)[2]]
    x2.comb <- cbind(x2.comb, x2.col)
  }
  
  x1.norm <- as.matrix(t(x1))
  x2.norm <- as.matrix(t(x2.comb))
  rownames(x2.norm) <- colnames(x2)
  
  if(plot){ 
    # par(mfrow=c(1,2))
    # boxplot(x1.norm, las=2, outline=F, main="Training data", xaxt="n")
    # boxplot(x2.norm, las=2, outline=F, main="Testing data", xaxt="n")
    
    par(mfrow = c(1,2))
    colramp = colorRampPalette(c(3,"white",2))(20)
    plot(density(x1.norm[1,]),col=colramp[1],lwd=3,ylim=c(0,.25),main="Training")
    for(i in 2:20){lines(density(x1.norm[i,]),lwd=3,col=colramp[i])}
    plot(density(x2.norm[1,]),col=colramp[1],lwd=3,ylim=c(0,.25),main="Testing")
    for(i in 2:20){lines(density(x2.norm[i,]),lwd=3,col=colramp[i])}
  }
  
  return(list(x1.norm=x1.norm,x2.norm=x2.norm))
}