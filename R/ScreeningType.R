#' @title Screening techniques (DAD- and BMD+DAD-screening)
#'
#' @description This function computes the and BMD+DAD-screening considering three different ranking cases 
#' in order to select the optimal threshold `d` and optimize data prediction. The three cases are:
#'  (a) ranking by `regression coefficients`; (b) ranking by `p.values`; (c) ranking by `regression coefficients` and `p.values`.
#' @param screenType three different type of screening based on a, b and c. The alternative are: `regCoef`, `pValue` or `regCoef+pValue`.
#' @param ranking the ranking of genes that are highly-correlated with the patient survival.
#' @param th the thresholds-dependent genes selected by the users.
#' @param genesList the list of genes identified using the BMD-screening. 
#'
#' @return A object contains the following two list of genes that are highly-correlated with the patient survival according to the type of ranking selected (a,b or c):
#'  \item{screenVars.DAD}{list of the screened variables using DAD-screening.}
#'  \item{screenVars.BMD.DAD}{list of the screened variables using BMD+DAD-screening.}
#' @export
ScreeningType <- function(screenType=c("regCoef","pValue", "regCoef+pValue"), ranking, th, genesList){

  # DAD-screening
  if(screenType=="regCoef"){
    screenVars.DAD <- list()
    for(i in 1:length(th)){
      screenVars.DAD[[i]] <- ranking[1:th[i],2]
      }
  }
  
  if(screenType=="pValue"){
    screenVars.DAD <- list(ranking[which(ranking[,4] < 0.05),2])
  }
  
  if(screenType=="regCoef+pValue"){
    screenVars.DAD <- list()
    for(i in 1:length(th)){
      ind <- which(ranking[1:th[i],4] < 0.05)
      screenVars.DAD[[i]] <- ranking[ind,2]
    }
  }
   
  # BMD+DAD-screening
  screenVars.BMD.DAD <- list()
  for(i in 1:length(screenVars.DAD)){
    # i=1
    screenVars.BMD.DAD[[i]] <- union(genesList,screenVars.DAD[[i]])
  }
  
  return(list(screenVars.DAD=screenVars.DAD, screenVars.BMD.DAD=screenVars.BMD.DAD))
}