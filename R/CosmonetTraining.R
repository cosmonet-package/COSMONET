
##' @title Fit network-regularized Cox regression models on training set
#'
#' @description This function fits penalized Cox regression methods in order to incorporate gene regulatory relationships and to select a subset
#' of potential biomarkers by using the training set \eqn{T}.
#' @param x1 input training matrix \eqn{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param family Cox proportional hazards regression model.
#' @param penalty penalty type. Can choose "Net" where \code{Omega} matrix is requested.
#' For \code{penalty = Net}, the penalty is defined as
#' \eqn{\lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}}, 
#' where \eqn{L} is a Laplacian matrix calculated from \code{Omega}.
#' @param alpha ratio between \code{L_1} and Laplacian for \code{Net}. Default is \code{alpha = 0.5}.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param ncv number of cross-validation performed for tuning optimal parameters (number of folds by default is 5).
#'
#' @return The following objects are returned:
#' \item{fitTrain}{values of the \code{NetworkCox()} function.}
#' \item{PI.train}{prognostic index \eqn{PI} on the training set. Compute the quantile \eqn{q_{\gamma}} of \eqn{PI} with \eqn{\gamma = \{0.20, 0.25, 0.30,... 0.80\}}. Each patient \eqn{i} in \eqn{T} is assigned to the high-risk (or low-risk) group if its prognostic index \eqn{PI_i} is above (or below) the \eqn{q_{\gamma}}-quantile.}
#' \item{cutoff.opt}{optimal cutoff selected adaptively on the training set \eqn{T}. It corresponds to the best separation in high-and-low risk group with respect to the log-rank test.}
#' \item{p.value}{\eqn{p}-value on the training set \eqn{T} (the significance level is \eqn{p}-value < 0.05).}
#' @export
CosmonetTraining <- function(x1, y1, screenVars, family="Cox", penalty="Net", Omega=NULL, alpha=0.5, ncv=NULL){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x1))
  } else {index <- match(screenVars,colnames(x1))}
  
  if(length(which(is.na(index)))==0){
    indexScreen=index} 
  else {
    indexScreen <- index[-which(is.na(index))]
    }
  
  penalty <- match.arg(penalty)
  
  if(family=="Cox") {
    fitTrain <- switch(penalty,
                       "Net"=NetworkCox(x1[,indexScreen], y1, penalty="Net", Omega[indexScreen,indexScreen], alpha, ncv))

    fitTrain$family <- "cox"
  }
  
  ## Compute the optimal cutoff on $T$
  select.cutoff <- SelectOptimalCutoff(x1[,indexScreen], y1, fitTrain$beta)
  opt.cutoff <- select.cutoff$opt.cutoff
  p.value <- select.cutoff$p.value
  PI.train <- select.cutoff$PI.train
  
  return(list(fitTrain=fitTrain, PI.train=PI.train, opt.cutoff=opt.cutoff, p.value=p.value))
}