
#' @title Validation test
#'
#' @description This function performs the validation test using the testing set \eqn{D}, the regression coefficients and the optimal cutoff computed on the training set \eqn{T}.
#' @param x2 input testing matrix \eqn{n2xp}.
#' @param y2 response variable, \code{y2} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param beta regression cofficients estimated from variable selection methods on the screened subset (BMD, DAD or BMD+DAD-screening).
#' @param opt.cutoff optimal cutoff selected adaptively on the training set \eqn{T}.
#' @param screening type of screening applied on the training set \eqn{T}. Can choose BMD, DAD or BMD+DAD. 
#'
#' @return The following objects are returned:
#' \item{PI.test}{prognostic index \eqn{PI} on the testing set \eqn{D}.}
#' \item{p.value}{\eqn{p}-value on the testing set \eqn{D} (the significance level is \eqn{p}-value < 0.05).}
#' \item{HL.groups}{high-and low-risk groups.}
#' \item{SurvPlot}{plot of survival curves.}
#' @export
ValidationTest <- function(x2, y2, beta, opt.cutoff, screening = c("BMD", "DAD", "BMD+DAD")){

 library(survival)
 # if(sum(beta)!=0){
    PI.test <- x2 %*% beta
    group <- vector(mode = "numeric", length = dim(x2)[1])

    for (i in c(1:nrow(PI.test))){
      if(PI.test[i] >= opt.cutoff){group[i]=1}
      else if(PI.test[i] < opt.cutoff){group[i]=2}
      #print(i)
        }
    ind.high <- which(group == 1)
    # high <- data.frame(high_risk_group=rownames(x2)[ind_high])
    ind.low <- which(group == 2)
    # low <- data.frame(low_risk_group=rownames(x2)[ind_low])

    HL.groups <- data.frame(sample=rownames(x2),group=group)
    HL.groups$risk <- ""
    HL.groups$risk[ind.high] <- "High"
    HL.groups$risk[ind.low] <- "Low"

    if((length(unique(group)) > 1)==TRUE){
      x <- data.frame(y2,groupRisk=group)
      logranktest <- survdiff(Surv(y2[,1], y2[,2]) ~ groupRisk, data = x, rho = 0)
      p.quantile <- 1-pchisq(logranktest$chisq, 1)
      p.value <- signif(p.quantile,3)
      # Kaplan-Meier survival curves test set
      fitTest <- survfit(Surv(y2[,1], y2[,2]) ~ groupRisk, data = x)
     
     filename <- sprintf("survPlot_%s.pdf",paste0(screening))
     #ggsave(filename, print(survp))
     pdf(filename)

     library(survminer)
     survp <- ggsurvplot(
     fitTest,                   # survfit object with calculated statistics.
     data = x,                  # data used to fit survival curves.
     conf.int = TRUE,           # show confidence intervals for point estimaes of survival curves.
     pval = TRUE,               # show p-value of log-rank test.
     risk.table = TRUE,         # show risk table.
     ggtheme = theme_minimal(), # customize plot and risk table with a theme.
     risk.table.y.text.col = T, # colour risk table text annotations.
     risk.table.col = "strata",
     risk.table.y.text = FALSE,
     legend.title = "",
     legend.labs = c("High Risk", "Low Risk")) # show bars instead of names in text annotations in legend of risk table

     print(survp, newpage = FALSE)
     dev.off()
    } else {p.value <- 1
           print("Warning: no splitting!")} # no splitting - only one group
   
   # } 
  
    return(list(PI.test=PI.test, p.value=p.value, HL.groups=HL.groups))
}
