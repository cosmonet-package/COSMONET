
#' @title Prognostic index (PI) and Kaplan-Meier curves on testing set
#'
#' @description This function compute the prognostic index \eqn{PI^{D}} on the testing set \eqn{D} using the regression coefficients and the optimal cutoff computed on the training set \eqn{T}. Then, the kaplan-Meier curves are shown as resulting of the log-rank test between high- and low-risk group.
#' @param x input testing matrix \eqn{nxp}.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param beta regression cofficients estimated on the training set \eqn{T}.
#' @param opt.cutoff optimal cutoff selected adaptively on the training set \eqn{T}.
#'
#' @return The following objects are returned:
#' \item{df}{data frame about sample, prognostic index, time, status and group risk on testing set \eqn{D}.}
#' \item{p.value}{\eqn{p}-value resulting from the log-rank test (the significance level is \eqn{p}-value < 0.05).}
#' \item{plots}{survival curves and distribution plot of prognostic index \eqn{PI^{D}}.}
#' @export
ValidationTest <- function(x,y,beta,opt.cutoff){
  
  # Compute the PI 
  PI <- as.matrix(x) %*% beta # [nxp]*[px1]=[nx1]
  PI <- data.frame(PI)
  PI <- merge(PI,y,by.x=0,by.y=0)
  PI <- PI[order(PI[,2]),]
  colnames(PI)[1] <- "sample"
  
  # We stratified patients into high- and low-risk subtypes based on PI_{T,*} 
  groupRisk <- vector(mode = "numeric", length = dim(PI)[1])
  for (i in c(1:nrow(PI))){
    if(PI[i,2] >= opt.cutoff){groupRisk[i]=1}
    else if(PI[i,2] < opt.cutoff){groupRisk[i]=2}
  }
  PI$groupRisk <- groupRisk
  
  if((length(unique(groupRisk)) > 1)==TRUE){
    library(survival)
    logranktest <- survdiff(Surv(time, status) ~ factor(groupRisk), data = PI, rho = 0)
    p.quantile <- 1-pchisq(logranktest$chisq, 1)
    p.value <- signif(p.quantile,3)
    
    # Kaplan-Meier survival curves
    fitTest <- survfit(Surv(time, status) ~ factor(groupRisk), data = PI)
    
    # filename <- sprintf("survPlot_%s.pdf",paste0(screening))
    # ggsave(filename, print(survp))
    # pdf(filename)
    
    library(ggplot2)
    survp <- survminer::ggsurvplot(
      fitTest,                  
      data = PI,                  
      conf.int = TRUE,           
      pval = p.value,               
      risk.table = TRUE,        
      ggtheme = theme_minimal(), 
      legend.title = paste0("Optimal Cutoff: ", opt.cutoff),
      legend.labs = c("High Risk", "Low Risk"))
    
    # print(survp)
    # print(survp, newpage = FALSE)
    # dev.off()
    
    # Distribution plot of PI 
    PI$id <- 1:dim(PI)[1]
    PI$groupRisk  <- as.factor(PI$groupRisk)
    d <- ggpubr:: ggline(PI,"id","PI",color = "groupRisk", palette = c("#FC4E07","#00AFBB")) +
      ylim(c(min(PI$PI),max(PI$PI))) +
      scale_color_discrete(labels = c("High Risk", "Low Risk"), name = paste0("Optimal cutoff: ",opt.cutoff)) +
      labs(x="Sample", y="Prognostic Index") 
    # print(d)
    
    ggsurv <- ggpubr::ggarrange(survp$plot, survp$table, heights = c(2, 0.7), ncol = 1, nrow = 2)
    print(ggpubr::ggarrange(ggsurv, d, labels = c("A", "B"), ncol = 2, nrow = 1))
    
  } else {p.value <- 1
  print("Warning: no splitting!")} # no splitting - only one group
  
  # } 
  
  df <- PI
  return(list(df=df, p.value=p.value))
}