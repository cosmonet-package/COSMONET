
#' @title Functional linkage matrix and Co-expression network matrix
#'
#' @description This function calculate the functional linkage matrix (\code{FL}) and the co-expression network matrix (\code{adjm}). 
#' Each element of the \code{FL} matrix represents the probability that two adjacent genes are correlated between them. 
#' The \code{adjm} matrix is a binary matrix 0/1 obtained from the \code{FL} matrix depending upon a user-selected threshold.
#' The \code{adjm} matrix has zero on the diagonal and it is used used for \code{penalty = Net} to calculate Laplacian matrix. 
#' @param x input matrix \eqn{nxp}
#' @param th user-selected threshold 
#'
#' @return The following objects are returned as matrices \eqn{pxp}:
#' \item{FL}{functional linkage matrix (probability matrix)}
#' \item{adjm}{co-expression network matrix (0/1 matrix)}
#' 
#' @references 
#'  Zhang, W., Ota, T., Shridhar, V., Chien, J., Wu, B., and Kuang, R. (2013),
#'  Network-based survival analysis reveals subnetwork signatures for predicting outcomes of ovarian cancer treatment. \cr
#'  \emph{PLoS computational biology, 9(3), e1002975}\cr \cr
#' @export
NetworkMatrix <- function(x, th){

  # n: number of samples
  # p: number of genes

  n <- dim(x)[1]
  p <- dim(x)[2]

  vec <- matrix(1,n)
  mean <- as.matrix(t(apply(x,2,mean)))
  x <- x - vec%*%mean

  W <- matrix(1,p,p)
  for(i in 1:p){
    for(j in (i+1):p){
      if(i+1<p){
      c <- sum(x[,i]*x[,j])/(sqrt(sum(x[,i]^2))*sqrt(sum(x[,j]^2)))
      W[i,j] <- abs(c)
      W[j,i] <- W[i,j]
      }
    }
  }

  IX <- apply(W, 2, function(x)(order(x,decreasing = TRUE)))
  IXI <- apply(IX, 2, function(x)(order(x,decreasing = FALSE)))

  W <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){W[i,j] <- 0} else {W[i,j] <- 1/(IXI[i,j])/(IXI[j,i])}
    }
  }

  Sum_R <- rowSums(W)
  Sum_C <- colSums(W)
  S <- matrix(0,p,p)

  for(i in 1:p){
    for(j in 1:p){
      S[i,j] <- W[i,j]/sqrt(Sum_R[i])/sqrt(Sum_C[j])
    }
  }

  # Functional linkage network
  FL <- S

  # 0-1 matrix
  adjm <- mat.or.vec(p,p)
  for (i in 1:p){
    for (j in 1:p){
      if(i==j){FL[i,j]=0}
      else if(FL[i,j]>=th){adjm[i,j]=1}
      else if(FL[i,j]<th){adjm[i,j]=0}
    }
  }

  colnames(FL) <- colnames(x)
  rownames(FL) <- colnames(x)

  colnames(adjm) <- colnames(x)
  rownames(adjm) <- colnames(x)
    
 return(list(FL=FL, adjm=adjm))
}