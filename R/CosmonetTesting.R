
#' @title Validation test and Kaplan Meier curve
#'
#' @description This function performes the model validation and prediction using the testing set \eqn{D}.
#' @param x2 input testing matrix \eqn{n2xp}. Each row represents an observation, while each column a variable.
#' @param y2 response variable, \code{y2} should be a three-column matrix with columns named sample, time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param beta regression coefficients estimated from variable selection methods on the training set.
#' @param opt.cutoff optimal cutoff selected adaptively on the training set \eqn{T}.
#' @param screening type of screening applied on the training set \eqn{T}. Can choose BMD, DAD or BMD+DAD. 
#'
#' @return The following objects are returned:
#' \item{PI.test}{prognostic index \eqn{PI} computed on the testing set \eqn{D}.}
#' \item{p.value}{\eqn{p}-value on the testing set \eqn{D} (the significance level is \eqn{p}-value < 0.05).}
#' \item{HL.groups}{high-risk and low-risk prognosis groups.}
#' \item{survival curves}{plot of survival curves.}
#' @export
#*****************************************************************************************************
# Validation set
#********************************************************************************************
ValidationTest <- function(x2, y2, beta, opt.cutoff, screening = c("BMD", "DAD", "BMD+DAD")){
 
  library(survival)
  
  # Compute the PI 
  PI <- x2 %*% beta
  df.PI <- data.frame(sample=rownames(PI),PI.score=PI)
  df.PI <- merge(df.PI,y2,by.x=1,by.y=1)
  df.PI <- df.PI[order(df.PI$PI.score),]
  rownames(df.PI) <- 1:dim(df.PI)[1]
  
  # We stratified patients into high- and low-risk subtypes based on PI_{T.*} 
  group <- vector(mode = "numeric", length = dim(df.PI)[1])
  for (i in c(1:nrow(df.PI))){
      if(df.PI$PI.score[i] >= opt.cutoff){group[i]=1}
      else if(df.PI$PI.score[i] < opt.cutoff){group[i]=2}
    }
  
  if((length(unique(group)) > 1)==TRUE){
    df.PI$id <- 1:dim(df.PI)[1]
    df.PI$risk <- group
    logranktest <- survdiff(Surv(y2[,2], y2[,3]) ~ risk, data = df.PI, rho = 0)
    p.quantile <- 1-pchisq(logranktest$chisq, 1)
    p.value <- signif(p.quantile,3)
    
    # Kaplan-Meier survival curves
    fitTest <- survfit(Surv(y2[,2], y2[,3]) ~ risk, data = df.PI)
    
    filename <- sprintf("survPlot_%s.pdf",paste0(screening))
    #ggsave(filename, print(survp))
    pdf(filename)
    
    library(survminer)
    survp <- ggsurvplot(
      fitTest,                   # survfit object with calculated statistics.
      data = df.PI,                  # data used to fit survival curves.
      conf.int = TRUE,           # show confidence intervals for point estimaes of survival curves.
      pval = p.value,               # show p-value of log-rank test.
      risk.table = TRUE,         # show risk table.
      ggtheme = theme_minimal(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.col = "strata",
      risk.table.y.text = FALSE,
      legend.title = paste0("Optimal Cutoff: ",round(opt.cutoff,4)),
      legend.labs = c("High Risk", "Low Risk"))
    
    print(survp, newpage = FALSE)
    dev.off()
    
    # Plot survival curves
    print(survp)
    
  } else {p.value <- 1
  print("Warning: no splitting!")} # no splitting - only one group
  
  # } 
  
  return(list(df.PI=df.PI, p.value=p.value))
}
#*********************************************************************************************************************
# Run CosmonetTesting() function
#*********************************************************************************************************************
CosmonetTesting <- function(x2, y2, screenVars, beta, opt.cutoff, screening = c("BMD", "DAD", "BMD+DAD")){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x2))
  } else {index <- match(screenVars,colnames(x2))}
  
  if(length(which(is.na(index)))==0){indexScreen=index
  } else {indexScreen <- index[-which(is.na(index))]}
  
  # Select the non-zero regression coeffiecients 
  index.non.zero.beta <- which(beta!=0)
  
  x2.screened <- x2[,indexScreen]
  testingValue <- ValidationTest(x2.screened[,index.non.zero.beta], y2, beta[index.non.zero.beta], opt.cutoff, screening)
  p.value <- testingValue$p.value
  HL.groups <- testingValue$RiskGroup
  PI.test <- testingValue$PI.test

  return(list(PI.test=PI.test,p.value=p.value, RiskGroup=RiskGroup))
}