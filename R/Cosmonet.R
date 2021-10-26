
#' @title Survival prediction via screening-network Cox methods
#'
#' @description This function performs two main step: (i) penalized Cox regression methods to select a subset of potential biomarkers by using the training set \code{T}; (ii) the validation test by using the testing set \code{D}.
#' @param k times to loop through cross validation.
#' @param x1 input training matrix \code{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param x2 input testing matrix \code{n2xp}. Each row represents an observation, while each column a variable.
#' @param y2 response variable, \code{y2} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param family Cox proportional hazards regression model. `Family=Cox`.
#' @param penalty penalty type. Can choose `Net` where `Omega` matrix is requested.
#' For penalty = `Net`, the penalty is defined as
#' \eqn{\lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}}, 
#' where `L` is a Laplacian matrix calculated from `Omega`.
#' @param alpha ratio between `L1` and Laplacian for `Net`. Default is alpha = 0.5`.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param lambda a user supplied decreasing sequence. If lambda = NULL, a random sequence of lambda is generated. For more details for instance `APML0` package
#' @param nlambda number of lambda values. Default is 50.
#' @param nfolds number of folds performed for tuning optimal parameters over runs. Default is `nfolds = 5`.
#' @param foldid an optional vector of values between 1 and nfolds specifying which fold each observation is in. 
#' @param selOptLambda a character string for selecting the lambda parameter. Options are `min` which uses the regularisation procedure implemented in the APML0 package or `1se"` to select the lambda parameter within one standard error from the optimal value.
#' @param optCutpoint a character string for choosing the optimal cutpoint on training set \code{T} based on prognostic index \code{PI^{T}}. Can choose `minPValue`, `median` and `survCutpoint`.
#' 
#' @return An object of class `COSMONET` is returned composed by:
#' \item{fitTrain}{see `CosmonetTraining` function outputs.}
#' \item{fitTest}{see `CosmonetTesting` function outputs.}
#' 
#' @details The first step is the variable screening of the data which aimed to reduce the number of variables for a 
#' large to a moderate scale. To this purpose, we assume that only a small number of these \code{p} variables is affecting 
#' the survival outcome. Therefore, we filter out variables that are considered not relevant for the disease under 
#' investigation. To this purpose, we consider three different types of variable screenings: biomedical 
#' screening (BMD-screening), data-driven screening (DAD-screening) and the fusion of biomedical and 
#' data-driven screening (BMD+DAD-screening). The second step is the application of penalized methods using 
#' the subset of screened variables \eqn{\{x_j,j \in \mathcal{I}\}} (where \eqn{\mathcal{I}} depends on the type screening performed) 
#' as new feature space to further remove not significant variables from the model. To assess the stability of 
#' the survival prediction we performed the k-fold cross-validation different times and we take as estimate 
#' the average value of \eqn{\lambda} and the corresponding \eqn{\alpha}. These two parameters are used to fit 
#' the corresponding penalized Cox model and obtain the parameter estimate of \eqn{\beta_\mathcal{I}}. 
#' Survival analysis is performed using the Kaplan Meier curves after dividing the patients in 
#' two risk groups (high-and-low risk group) on the basis of the prognostic index \code{PI} computed with the gene signature.
#' The \code{p-value}, used to test the null hypothesis that the survival curves are identical vs. the alternative that 
#' the two groups have different survival, is calculated by using the \code{log-rank test}.
#' 
#' @references Iuliano, A., Occhipinti, A., Angelini, C., De Feis, I., and Liò, p. (2018).
#' Combining Pathway Identification and Breast Cancer Survival Prediction via Screening-Network Methods. \cr
#' \emph{Frontiers in genetics, 9, 206.}\cr \cr
#' Iuliano, A., Occhipinti, A., Angelini, C., De Feis, I., & Lió, P. (2016). 
#' Cancer markers selection using network-based Cox regression: A methodological and computational practice. \cr
#' \emph{Frontiers in physiology, 7, 208.}\cr \cr
#' @export
Cosmonet <- function(k,x1,y1,x2,y2,screenVars,family="Cox",penalty="Net",Omega,alpha=0.5,lambda=NULL,nlambda=50,nfolds=5,foldid=NULL,selOptLambda=c("min","1se"),optCutpoint=c("minPValue","median","survCutpoint")){
  
  # Training phase
  fitTrain <- CosmonetTraining(k,x1,y1,screenVars,family="Cox",penalty="Net",Omega,alpha,lambda,nlambda,folds,foldid,selOptLambda,optCutpoint)
  
  # Regression coefficients and optimal cutoff computed  on training set
  beta <- fitTrain$beta
  opt.cutoff <- fitTrain$opt.cutoff
    
  # Testing phase
  fitTest <- CosmonetTesting(x2,y2,screenVars,beta,opt.cutoff)
                               
  obj <- c(list(fitTrain=fitTrain, fitTest=fitTest))
  class(obj) = "COSMONET"
  obj
}