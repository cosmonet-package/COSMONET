
#' @title Marginal Cox ranking function
#'
#' @description This function computes the ranking of genes that are highly-correlated with the patient survival according to their corresponding marginal cox utilities from the largest to the smallest. 
#' @param x1 input training matrix \eqn{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#'
#' @return A data frame composed by the ranking of genes that are highly-correlated with the patient survival, the corresponding regression coefficients and p.values.
#' @references 
#' Fan, J., Feng, Y., & Wu, Y. (2010). 
#' High-dimensional variable selection for Cox’s proportional hazards model.\cr
#' \emph{In Borrowing strength: Theory powering applications–a Festschrift for Lawrence D.} Brown (pp. 70-86). Institute of Mathematical Statistics.\cr \cr
#' @export
MarginalCoxRanking <- function(x1, y1){

  s <- NULL
  res <- 0
  beta <- NULL
  pvals <- NULL

  x1 <- data.frame(x1)
  n <- dim(x1)[1]
  p <- dim(x1)[2]

  times <- y1[,1]
  status <- y1[,2]

  library(survival)
  for(m in 1:p){
    beta[m] <- coxph(Surv(times,status)~x1[,m], x1)$coefficients
    pvals[m] <- broom::tidy(coxph(Surv(times,status)~x1[,m], x1))$p.value
  }

  # R(t_i):= {X_j:t_j>=t_i, j=1,...,n}

  ind <- which(status==1)
  # print(length(ind))
  for(t in 1:length(ind)){
    i <- ind[t]
    ex <- exp(x1[which(times >= times[i]),]*beta)
    if(length(which(times >= times[i]))>1){
      s <- colSums(ex)
    }
    else{
      s <- ex
    }
    res <- res + status[i]*log(s)
    # print(t)
  }

  x1 <- as.matrix(x1)
  u_m <- status[ind]%*%(x1[ind,]*beta) - res # marginal utilities

  index <-order(u_m,decreasing = T)
  rankingData <- data.frame(ranking=index, symbol=colnames(x1)[index], coef= beta[index], p.value=pvals[index])
  #value <- (u_m[index])
  
  return(rankingData=rankingData)
  
}