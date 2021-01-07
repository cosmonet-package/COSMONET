
#' @title Model Validation and Prediction
#'
#' @description This function performes the model validation and prediction using the testing set \eqn{D}.
#' @param x2 input testing matrix \eqn{n2xp}.
#' @param y2 response variable, \code{y2} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param beta regression cofficients estimated from variable selection methods on the training set.
#' @param opt.cutoff optimal cutoff selected adaptively on the training set \eqn{T}.
#' @param screening type of screening applied on the training set \eqn{T}. Can choose BMD, DAD or BMD+DAD. 
#'
#' @return The following objects are returned:
#' \item{PI.test}{prognostic index \eqn{PI} computed on the testing set \eqn{D}.}
#' \item{p.value}{\eqn{p}-value on the testing set \eqn{D} (the significance level is \eqn{p}-value < 0.05).}
#' \item{HL.groups}{high-risk and low-risk prognosis groups.}
#' \item{survival curves}{plot of survival curves.}
#' @export
CosmonetTesting <- function(x2, y2, screenVars, beta, opt.cutoff, screening = c("BMD", "DAD", "BMD+DAD")){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x2))
  } else {index <- match(screenVars,colnames(x2))}
  
  if(length(which(is.na(index)))==0){
    indexScreen=index} 
  else {
    indexScreen <- index[-which(is.na(index))]
  }
  
  beta <- as.matrix(beta)
  testingValue <- ValidationTest(x2[,indexScreen], y2, beta, opt.cutoff, screening)
  p.value <- testingValue$p.value
  HL.groups <- testingValue$HL.groups
  PI.test <- testingValue$PI.test

  return(list(PI.test=PI.test,p.value=p.value, HL.groups=HL.groups))
}