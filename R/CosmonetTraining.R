
##' @title Fit network-regularized Cox regression models on training set
#'
#' @description This function fits penalized Cox regression methods in order to incorporate gene regulatory relationships and to select signature genes using the training set \eqn{T}.
#' @param k times to loop through cross validation.
#' @param x input training matrix \eqn{nxp}. Each row is an observation vector.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param family Cox proportional hazards regression model.
#' @param penalty penalty type. Can choose "Net" where \code{Omega} matrix is requested.
#' For \code{penalty = Net}, the penalty is defined as
#' \eqn{\lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}}, 
#' where \eqn{L} is a Laplacian matrix calculated from \code{Omega}.
#' @param alpha ratio between \code{L_1} and Laplacian for \code{Net}. Default is \code{alpha = 0.5}.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param nfolds number of cross-validation performed for tuning optimal parameters over runs. Default is \code{nfolds = 5}.
#' @param plot plot the cross-validation curve as a function of the lambda values used.
#'
#' @return The following objects are returned:
#' \item{beta}{a sparse Matrix of coefficients, stored in class\code{dgCMatrix}.}
#' \item{opt.lambda}{value of \eqn{lambda} based on minimum \eqn{cvm} over runs.}
#' \item{df}{data frame of sample, prognostic index, time, status and group risk for each quantile \eqn{q_{gamma}}, with \eqn{\gamma=0.20, ..., 0.80}.}
#' \item{summary}{summary information about number of patients at risk, cutoff and p.value for each quantile.}
#' \item{p.value}{\eqn{p}-value resulting from the log-rank test (the significance level is \eqn{p}-value < 0.05).}
#' \item{plots}{survival curves and distribution plot of prognostic index \eqn{PI^{T}}.}
#' @export
CosmonetTraining <- function(k,x,y,screenVars,family="Cox",penalty="Net",Omega,alpha=0.5,nfolds=5,plot=TRUE){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x))
  } else {index <- match(screenVars,colnames(x))}
  
  if(length(which(is.na(index)))==0){indexScreen=index
  } else {indexScreen <- index[-which(is.na(index))]}
  
  penalty <- match.arg(penalty)
  
  if(family=="Cox"){
    fitTrain <- switch(penalty,
                       "Net"=NetworkCox(k,x[,indexScreen],y,Omega[indexScreen,indexScreen],alpha,nfolds,plot))
    fitTrain$family <- "cox"
  }
  
  beta <- fitTrain$beta
  opt.lambda <- fitTrain$opt.lambda

  ## Compute the optimal cutoff on $T$
  x.screened <- x[,indexScreen]
  index.non.zero.beta <- which(fitTrain$beta!=0)
  select.cutoff <- SelectOptimalCutoff(x.screened[,index.non.zero.beta],y,fitTrain$beta[index.non.zero.beta])
  df <- select.cutoff$df
  summary <- select.cutoff$summary
  opt.cutoff <- select.cutoff$opt.cutoff
  p.value <- select.cutoff$p.value
  
  return(list(beta=beta,opt.lambda=opt.lambda,df=df,summary=summary,opt.cutoff=opt.cutoff,p.value=p.value))
}

Lambdas <- function(x,y,Omega=Omega,alpha=alpha,nfolds=nfolds){
  # Cross Validation
  # set.seed(4321)
  cv <- APML0(as.matrix(x),as.matrix(y),family="cox",penalty="Net",Omega=Omega,alpha=alpha,
              nlambda=100,nfolds=nfolds,ifast=TRUE,isd = FALSE,ifastr = TRUE)
  require(data.table)
  return(data.table(cvm=cv$fit$cvm, lambda=cv$fit$lambda))
}

OptimLambda <- function(k, plot=FALSE, ...){
  # Returns optimal lambda for APML0.
  #
  # Args:
  #   k: # times to loop through cross validation
  #   ...: Other args passed to cross validation
  #
  # Returns:
  #   Lambda associated with minimum average CV error over runs.
  require(parallel)
  require(plyr)
  MSEs <- data.frame(rbind.fill(mclapply(seq(k), function(dummy) Lambdas(...))))
  res <- ddply(MSEs, .(lambda), summarise, mean=mean(cvm))
  
  if(plot){ 
    # Plot Cross-validation curve
    x <- log(as.matrix(res$lambda))
    y <- as.matrix(res$mean)
    plot(x,y,"p",xlab = "Log(lambda)",ylab ="cvm",col = "red", pch = 20)
    abline(v=log(res[order(res$mean),][1,1]), col="black",lwd=3, lty=2)
  }
  
  return(res[order(res$mean),][1,1])
}

NetworkCox <- function(k,x,y,Omega=Omega,alpha=alpha,nfolds=nfolds,plot=TRUE){
  
  require(APML0)
  set.seed(1234)
  count <- 0
  colnames(y) <- c("time","status")
  # Returns optimal lambda for APML0.
  opt.lambda <- OptimLambda(k,x,y,Omega,alpha,nfolds,plot)
  # Fit optimal model
  fit <- APML0(as.matrix(x),as.matrix(y),family="cox",penalty="Net",Omega=Omega,lambda=opt.lambda,alpha)
  beta <- as.matrix(fit$Beta)
  rownames(beta) <- colnames(x)
  indexBeta <- which(beta!=0)
  nonzerobeta <- beta[indexBeta,]
  
  if(length(indexBeta)==0){count = count + 1;
  if(count > 5) {print(sprintf("More than 5 beta, file are null"));}
  }
  
  return(list(beta=beta,opt.lambda=opt.lambda,nonzerobeta=nonzerobeta))
}
