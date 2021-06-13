
#' @title Marginal Cox ranking function
#'
#' @description This function ranks the genes that are highly-correlated with the patient survival according to their corresponding marginal cox utilities. The genes are ordered from the largest to the smallest. 
#' @param x1 input matrix \eqn{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#'
#' @return A data frame composed by the ranking of genes (original ndices and gene symbols) that are highly-correlated with the patient survival, the marginal regression coefficients, the coefficients and the relative p.values.
#' @references 
#' Fan, J., Feng, Y., & Wu, Y. (2010). 
#' High-dimensional variable selection for Cox’s proportional hazards model.\cr
#' \emph{IMS Collections}, \bold{6}, 70-86.\cr \cr
#' Saldana, D. F. and Feng, Y. (2018).
#' SIS: An R package for Sure Independence Screening in Ultrahigh Dimensional Statistical Models.\cr
#'  \emph{Journal of Statistical Software}, \bold{83}, 2, 1-25.\cr\cr
#' @export
MarginalCoxRanking <- function(x, y){

  library(plyr)
  library(survival)
  
  x.stand <- data.frame(standardize(x))
  n <- dim(x.stand)[1]
  p <- dim(x.stand)[2]

  # survival data
  times <- y[,1]
  status <- y[,2]

  beta <- NULL
  pvals <- NULL
  
  library(plyr)
  pbar <- create_progress_bar('text')
  pbar$init(p)
  
  for(m in 1:p){
    beta[m] <- coxph(Surv(times,status)~x.stand[,m], x.stand)$coefficients
    pvals[m] <- broom::tidy(coxph(Surv(times,status)~x.stand[,m], x.stand))$p.value
    pbar$step()
    }
  
  # Use function from SIS library
  y_surv <- Surv(times,status)
  x.stand <- standardize(x)
  margcoef <- margcoef(x.stand, y_surv, family = "cox", null.model = FALSE, iterind = 100)
  
  # Sort by margcoef
  rankcoef <- sort(margcoef, decreasing = TRUE, index.return = TRUE)
  rankingData <- data.frame(rankcoef$ix, colnames(x)[rankcoef$ix], margcoef[rankcoef$ix], beta[rankcoef$ix], pvals[rankcoef$ix])
  colnames(rankingData) <- c("original.index", "symbol", "marg.coeff", "beta", "pvals")
  # head(rankingData)
  
  # par(mfrow = c(2, 2))
  # plot(rankingData$marg.coeff, rankingData$pvals, col = "blue", pch = 20)
  # abline(h=0.05, col="red",lwd=2, lty=2)
  # plot(rankingData$pvals, 1:p, col = "blue", pch = 20)
  # abline(v=0.05, col="red",lwd=2, lty=2)
  # abline(h=4000, col="red",lwd=2, lty=2)
  # plot(rankingData$beta, 1:length(rankingData$beta),col = "blue", pch = 20)
  # plot(rankingData$marg.coeff,abs(rankingData$beta),col = "blue", pch = 20)

  return(rankingData=rankingData)
  
}

# Standardizes the columns of a high-dimensional design matrix to mean zero and unit Euclidean norm.
standardize <- function(x){
  center <- colMeans(x)
  x.c <- sweep(x, 2, center)
  unit.var <- sqrt(apply(x.c, 2, crossprod))
  val <- sweep(x.c, 2, unit.var, "/")
  return(val)
}

# Marginal utility function from SIS package according to Fan et al. 2010 page 75
margcoef <- function(x, y, condind = NULL, family, null.model = FALSE, iterind){
  n <- dim(x)[1]
  p <- dim(x)[2]
  ones <- rep(1, n)
  candind <- setdiff(1:p, condind)
  if (iterind == 0) {
    if (family == "cox") 
      margcoef <- abs(cor(x, y[, 1])) else margcoef <- abs(cor(x, y))
  } else {
    if (null.model == TRUE) {
      if (is.null(condind) == TRUE) {
        x <- x[sample(1:n), ]
      }
      if (is.null(condind) == FALSE) {
        x[, candind] <- x[sample(1:n), candind]
      }
    }
    margcoef <- abs(sapply(candind, mg, x, y, ones, family, condind))
  }
  return(margcoef)
}

mg <- function(index, x = x, y = y, ones = ones, family = family, condind = condind) {
  margfit <- switch(family, gaussian = coef(glm.fit(cbind(ones, x[, index], x[,condind]), y, family = gaussian()))[2], binomial = coef(glm.fit(cbind(ones, x[, index], x[, condind]), y, family = binomial()))[2], poisson = coef(glm.fit(cbind(ones, x[, index], x[, condind]), y, family = poisson()))[2], cox = coef(coxph(y ~ cbind(x[, index], x[,condind])))[1])
}

# # For each covariate define its marginal utility as the maximum of the partial likelihood of the single covariate.
# utility <- NULL
# x1 <- as.matrix(x1)
# beta <- as.matrix(beta) # [px1]
# 
# # Compute risk set: R(t_i):= {X_j:t_j>=t_i, j=1,...,n}
# index <- which(status==1) 
# # print(length(ind))
# for(m in 1:p){
#   # m=1
#   for(t in 1:length(index)){
#     expFun <- exp(x1[which(times >= times[index[t]]),m]*beta[m]) # px1
#     if(length(which(times >= times[index[t]]))>1){
#       s <- sum(expFun)
#     } else {s <- expFun}
#     res <- status[index[t]]*log(s)
#     u <- max(status[index[t]]*(x1[index[t],m]*beta[m]) - res)
#   }
#   # marginal utilities [px1]
#   utility <- rbind(utility, u)
# }
# rankcoef <- sort(utility, decreasing = TRUE, index.return = TRUE)
# ix0 <- rankcoef$ix
# rankingData <- data.frame(ranking=ix0, symbol=colnames(x1)[ix0], 
#                           utility=rankcoef$x[ix0], coef= beta[ix0], p.value=pvals[ix0])
