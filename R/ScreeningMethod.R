#' @title Screening methods (BMD-, DAD- and BMD+DAD-screening)
#'
#' @description This function computes the BMD, DAD and BMD+DAD-screening considering three different ranking cases 
#' in order to select the optimal threshold and optimize data prediction.
#' @param x1 input training matrix \eqn{n1xp}.
#' @param myList list of cancer-specific genes.
#' @param screening type of screening (BMD-, DAD- and BMD+DAD-screening). 
#' @param ranking ranking of genes that are highly-correlated with the patient survival. Default is `ranking=NULL`, i.e., for BMD-screening.  
#' @param thresh number of the ranked genes choosen. Default is `thresh=NULL`, i.e., for BMD-screening.  
#'
#' @return A object contains the following list of genes:
#'  \item{screenVars}{list of the screened genes.}
#' @export
ScreeningMethod <- function(x1,myList,screening=c("bmd","dad","bmd+dad"), ranking=NULL,thresh=NULL){

  # BMD-screening
  if(screening=="bmd"){
    pvars <- colnames(x1)
    screenVars <- list(intersect(myList, pvars))
    print(paste0("Number of BMD-screened genes: ", length(screenVars[[1]])))
  }
    
  # DAD-screening
  if(screening=="dad"){
    screenVars <- list()
    for(i in 1:length(thresh)){
      screenVars[[i]] <- ranking[1:thresh[i],2]
      print(paste0("Number of DAD-screened genes: ", length(screenVars[[i]])))
    }
  }
   
  # BMD+DAD-screening
  if(screening=="bmd+dad"){
    # DAD
    screenVars.DAD <- list()
    for(i in 1:length(thresh)){
      screenVars.DAD[[i]] <- ranking[1:thresh[i],2]
    }
    # Union BMD and DAD
    screenVars <- list()
    for(i in 1:length(screenVars.DAD)){
      # i=1
      screenVars[[i]] <- union(myList,screenVars.DAD[[i]])
      print(paste0("Number of BMD+DAD-screened genes: ", length(screenVars[[i]])))
    }
  }
  
  return(list(screenVars=screenVars))
}