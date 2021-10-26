
#' @title Make predictions on data
#'
#' @description This function performes the model validation and prediction using the testing set \code{D}. It generates the Kaplan-Meier curve resulting from the log-rank test between high and low-risk group and the distribution plot of prognostic index \code{PI^{D}} computed using the signature genes and the optimal cutoff \code{PI^{*,T}} obtained from the training set \code{T}.
#' @param x input testing matrix \code{nxp}. Each row is an observation vector.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param beta regression coefficients estimated on the training set.
#' @param opt.cutoff optimal cutoff selected adaptively on the training set \code{T}.
#'
#' @return The following objects are returned:
#' \item{df}{data frame composed by prognostic index \code{PI^{D}}, sample, time and status.}
#' \item{p.value}{resulting from the log-rank test (the significance level is \code{p}-value < 0.05).}
#' @export
CosmonetTesting <- function(x,y,screenVars,beta,opt.cutoff){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x))
  } else {index <- match(screenVars,colnames(x))}
  
  if(length(which(is.na(index)))==0){indexScreen=index
  } else {indexScreen <- index[-which(is.na(index))]}

  x.screened <- x[,indexScreen]
  index.non.zero.beta <- which(beta!=0)
  testingValue <- ValidationTest(x.screened[,index.non.zero.beta],y,beta[index.non.zero.beta],opt.cutoff)
  df <- testingValue$df
  p.value <- testingValue$p.value

  return(list(df=df,p.value=p.value))
}