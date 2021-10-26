
##' @title Fit network-regularized Cox regression models on training set
#'
#' @description This function fits penalized Cox regression methods in order to incorporate gene regulatory relationships and to select signature genes using the training set \code{T}.
#' @param k times to loop through cross validation.
#' @param x input training matrix \code{nxp}. Each row is an observation vector.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening or by the user. A list or a character string can be used.
#' @param family Cox proportional hazards regression model. `Family=Cox`
#' @param penalty penalty type. Can choose `Net` where `Omega` matrix is requested. For \code{penalty = Net}, the penalty is defined as \eqn{\lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}},
#' where `L` is a Laplacian matrix calculated from `Omega`.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param alpha ratio between `L_1` and Laplacian for `Net`. Default is `alpha = 0.5`.
#' @param lambda a user supplied decreasing sequence. If `lambda = NULL`, a sequence of lambda is generated based on `nlambda` and `rlambda` (for more details,  see `APML0` package). Supplying a value of lambda overrides this.
#' @param nlambda number of lambda values. Default is 50.
#' @param nfolds number of folds performed for tuning optimal parameters over runs. Default is `nfolds = 5`.
#' @param foldid an optional vector of values between 1 and nfolds specifying which fold each observation is in. 
#' @param selOptLambda a character string for selecting the lambda parameter. Options are `min` which uses the regularisation procedure implemented in the APML0 package or `1se"` to select the lambda parameter within one standard error from the optimal value.
#' @param optCutpoint a character string for choosing the optimal cutpoint on training set \code{T} based on prognostic index \code{PI^{T}}. Can choose `minPValue`, `median` and `survCutpoint`.
#'
#' @return The following objects are returned:
#' \item{beta}{a sparse Matrix of coefficients, stored in class\code{dgCMatrix}.}
#' \item{opt.lambdas}{\eqn{lambda} values based on minimum \code{cvm} over runs and on \code{1se} to select the lambda parameter within one standard error from the optimal value.}
#' \item{df}{data frame composed by samples, relative prognostic indices, times, status and group risk for each quantile \eqn{q_{gamma}}, with \eqn{\gamma=0.20, ..., 0.80}.}
#' \item{summary}{summary table on number of patients at risk, cutoff and p.value for each quantile.}
#' \item{opt.cutoff}{optimal cutoff selected on the training.}
#' \item{p.value}{resulting from the log-rank test (the significance level is \code{p}-value < 0.05).}
#' @export
CosmonetTraining <- function(k,x,y,screenVars,family="Cox",penalty="Net",Omega,alpha=0.5,lambda=NULL,nlambda=50,nfolds=5,foldid=NULL,selOptLambda=min("min","1se"),optCutpoint=c("minPValue","median","survCutpoint")){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x))
  } else {index <- match(screenVars,colnames(x))}
  
  if(length(which(is.na(index)))==0){indexScreen=index
  } else {indexScreen <- index[-which(is.na(index))]}
  
  penalty <- match.arg(penalty)
  
  if(family=="Cox"){
    fitTrain <- switch(penalty,
                       "Net"=NetworkCox(k,x[,indexScreen],y,Omega[indexScreen,indexScreen],alpha,lambda,nlambda,nfolds,foldid,selOptLambda))
    fitTrain$family <- "cox"
  }
  
  beta <- fitTrain$beta
  opt.lambda <- fitTrain$opt.lambda

  ## Compute the optimal cutoff on $T$
  x.screened <- x[,indexScreen]
  index.non.zero.beta <- which(fitTrain$beta!=0)
  select.cutoff <- SelectOptimalCutoff(x.screened[,index.non.zero.beta],y,fitTrain$beta[index.non.zero.beta],optCutpoint)
  df <- select.cutoff$df
  summary <- select.cutoff$summary
  opt.cutoff <- select.cutoff$opt.cutoff
  p.value <- select.cutoff$p.value
  
  return(list(beta=beta,opt.lambda=opt.lambda,df=df,summary=summary,opt.cutoff=opt.cutoff))
}

getLambda.min <- function(x,y,Omega=Omega,alpha=alpha,lambda=NULL,nlambda=50,nfolds=5,foldid=NULL){
  # Cross Validation
  library(APML0)
  cv <- APML0(as.matrix(x),as.matrix(y),family="cox",penalty="Net",Omega=Omega,alpha=alpha,
              lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,ifast=TRUE,isd = FALSE,ifastr = TRUE)
  require(data.table)
  return(data.table(cvm=cv$fit$cvm, lambda=cv$fit$lambda))
}

getLambda.1se <- function(x,y,Omega=Omega,alpha=alpha,lambda=NULL,nlambda=50,nfolds=5,foldid=NULL){
  # Cross Validation
  library(APML0)
  cv <- APML0(as.matrix(x),as.matrix(y),family="cox",penalty="Net",Omega=Omega,alpha=alpha,
              lambda=lambda,nlambda=nlambda,nfolds=nfolds,foldid=foldid,ifast=TRUE,isd = FALSE,ifastr = TRUE)
  lambda <- cv$fit$lambda
  cvm <- cv$fit$cvm    
  cvsd <- cv$fit0$cvse         # cvse computed with l0 regularisation
  cvmin <- cv$fit0$cvm         # cvmin computed with l0 regularisation
  return(list(lambda=lambda,cvm=cvm,cvsd=cvsd,cvmin=cvmin))
}

OptimLambda <- function(k,x,y,Omega=Omega,alpha=alpha,lambda=NULL,nlambda=50,nfolds=5,foldid=NULL,selOptLambda=c("min","1se")){
  # Returns optimal lambda for APML0.
  #
  # Args:
  #   k: # times to loop through cross validation
  #   ...: Other args passed to cross validation
  #
  # Returns:
  # min: lambda associated with minimum average CV error over runs.
  # 1se: largest value of lambda such that error is within 1 standard error of the minimum.
  require(parallel)
  require(plyr)
  require(pbmcapply)
  require(APML0)
  
  MSEs <- data.frame(rbind.fill(pbmclapply(seq(k), function(dummy) getLambda.min(x,y,Omega,alpha,lambda,nlambda,nfolds,foldid))))
  res <- ddply(MSEs, .(lambda), summarise, mean=mean(cvm))
  
  if(selOptLambda=="min"){
  opt.lambda <- res[order(res$mean),][1,1]
  } else if(selOptLambda=="1se"){
    lambda.vec <- NULL
    for(i in seq(k)){
      pars <- getLambda.1se(x,y,Omega,alpha,lambda,nlambda,nfolds,foldid)
      semin <- (pars$cvmin + pars$cvsd)
      id1se <- pars$cvm <= semin
      lambda.1se <- max(pars$lambda[id1se], na.rm = TRUE) # change this to min to change direction of lambda selection
      lambda.vec <- rbind(lambda.vec,lambda.1se)
    }
    opt.lambda <- lambda.vec[order(lambda.vec)][1]
  }
  
  # Plot Cross-validation curve
  x <- log(as.matrix(res$lambda))
  y <- as.matrix(res$mean)
  plot(x,y,"p",xlab = "Log(lambda)",ylab ="cvm",col = "red", pch = 20)
  abline(v=log(opt.lambda), col="black",lwd=3, lty=2)
  
  return(list(opt.lambda=opt.lambda))
}

NetworkCox <- function(k,x,y,Omega=Omega,alpha=alpha,lambda=NULL,nlambda=50,nfolds=5,foldid=NULL,selOptLambda=c("min","1se")){
  
  require(APML0)
  colnames(y) <- c("time","status")
  
  # Get Optimal Lambda
  opt.lambda <- OptimLambda(k,x,y,Omega,alpha,lambda,nlambda,nfolds,foldid,selOptLambda)$opt.lambda
  
  # Fit optimal model
  fit <- APML0(as.matrix(x),as.matrix(y),family="cox",penalty="Net",Omega=Omega,alpha,lambda=opt.lambda)
  beta <- as.matrix(fit$Beta)
  rownames(beta) <- colnames(x)
  nonzerobeta <- beta[which(beta!=0),]

  return(list(beta=beta,opt.lambda=opt.lambda,nonzerobeta=nonzerobeta))
}

