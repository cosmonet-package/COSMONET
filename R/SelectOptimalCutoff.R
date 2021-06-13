
#' @title Select adaptively the PI-optimal cutoff on the training set 
#'
#' @description This function computes the optimal cutoff \eqn{PI^{*,T}} on the training set \eqn{T}, i.e., the value that corresponds to the best separation in high-and-low risk group with respect to the log-rank test using the prognostic index \eqn{PI^{T}}.
#' @param x input training matrix \eqn{nxp}.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param beta regression coefficients estimated on the training set \eqn{T}.
#'
#' @return The following objects are returned:
#' \item{df}{data frame abount sample, prognostic index, time, status and group risk for each quantile \eqn{q_{gamma}}, with \eqn{\gamma=0.20, ..., 0.80}.}
#' \item{summary}{summary information about number of patients at risk, cutoff and p.value for each quantile.}
#' \item{cutoff.opt}{optimal cutoff selected adaptively.}
#' \item{p.value}{\eqn{p}-value resulting from the log-rank test (the significance level is \eqn{p}-value < 0.05).}
#' \item{plots}{survival curves and distribution plot of prognostic index \eqn{PI^{T}}.}
#' @export
SelectOptimalCutoff <- function(x,y,beta){
  
  # if(sum(beta)!=0){
  PI <- as.matrix(x) %*% beta # [nxp]*[p]=[n]
  PI <- data.frame(PI)
  PI <- merge(PI,y,by.x=0,by.y=0)
  PI.data <- PI[order(PI[,2]),]
  colnames(PI.data)[1] <- "sample"
  
  # Compute the quantile q_gamma, with gamma = 0.2, ..., 0.80
  probs <- seq(0.2 ,0.8, by=0.05)
  perc <- paste(probs*100, "%", sep="")
  q <- quantile(PI.data[,2], probs=probs)
  p <- vector()
  # par(mfrow=c(3,3))
  
  # We stratified patients into high- and low-risk subtypes based on PI
  highRisk <- NULL
  lowRisk <- NULL
  risk.matrix <- NULL
  for (j in 1:length(q)){
    # print(paste0("cutoff: ",q[j]))
    groupRisk <- vector(mode = "numeric", length = dim(PI.data)[1]);
    
    for (i in c(1:nrow(PI.data))){
      if (PI.data[i,2] >= q[j]){groupRisk[i]=1}
      else if (PI.data[i,2]< q[j]){groupRisk[i]=2}
    }
    
    PI.data$groupRisk <- groupRisk
    
    ind.high <- which(groupRisk == 1)
    HR <- length(ind.high)
    highRisk <- cbind(highRisk, HR)
    ind.low <- which(groupRisk == 2)
    LR <- length(ind.low)
    lowRisk <- cbind(lowRisk, LR)
    lowHighRisk <- rbind(lowRisk,highRisk)
    
    if((length(unique(groupRisk)) > 1)==TRUE){
      # print(table(groupRisk)
      library(survival)
      logranktest <- survdiff(Surv(time, status) ~ factor(groupRisk), data = PI.data, rho = 0)
      p.quantile <- 1-pchisq(logranktest$chisq, 1)
      p[j] <- as.vector(signif(p.quantile,3))
      # print(paste0("p.value: ", p[j]))
      fitTrain <- survfit(Surv(time, status) ~ factor(groupRisk), data = PI.data)
      # cutoff <- signif(q[j],3)
    } else {
      p <- rep(1, length(q))  # no splitting - only one group
    }
    
    # Keplan-Meier curve
    library(ggplot2)
    survp <- survminer::ggsurvplot(
      fitTrain,                  
      data = PI.data,                  
      conf.int = TRUE,  
      pval = p[j],              
      risk.table = TRUE,  
      ggtheme = theme_minimal(),  
      legend.title = paste0("Cutoff: ",round(q[j],6)),
      legend.labs = c("High Risk", "Low Risk")
    )
    # print(survp)
    
    # Distribution plot of PI 
    PI.data$id <- 1:dim(PI.data)[1]
    PI.data$groupRisk  <- as.factor(PI.data$groupRisk)
    d <- ggpubr:: ggline(PI.data,"id","PI",color = "groupRisk", palette = c("#FC4E07","#00AFBB")) +
      ylim(c(min(PI.data$PI),max(PI.data$PI))) +
      scale_color_discrete(labels = c("High Risk", "Low Risk"), name = paste0("Cutoff: ",round(q[j],6))) +
      labs(x="Sample", y="Prognostic Index") 
    # print(d)
    
    ggsurv <- ggpubr::ggarrange(survp$plot, survp$table, heights = c(2, 0.7), ncol = 1, nrow = 2)
    print(ggpubr::ggarrange(ggsurv, d, labels = c("A", "B"), ncol = 2, nrow = 1))
    
    risk.matrix <- cbind(risk.matrix, groupRisk)
  }
  
  # Risk matrix
  risk.matrix <- data.frame(risk.matrix)
  rownames(risk.matrix) <- PI.data$sample
  colnames(risk.matrix) <- perc
  
  df <- merge(PI.data,risk.matrix,by.x=1,by.y=0)
  df <- df[order(df$PI,decreasing = FALSE),]
  rownames(df) <- 1:dim(df)[1]
    
  lowHighRisk <- as.data.frame(lowHighRisk)
  rownames(lowHighRisk) <- c("Low Risk", "High Risk")
  colnames(lowHighRisk) <- as.character(perc)
  # print(lowHighRisk)
  
  q.p.values <- as.data.frame(rbind(q,p))
  rownames(q.p.values) <- c("cutoff","p.value")
  # print(q.p.values)
  
  summary <- data.frame(t(lowHighRisk),t(q.p.values))
  # print(summary)
  
  index.non.zero <- which(q.p.values[2,]!=0)
  
  if(length(index.non.zero)==1){
    opt.cutoff <- q.p.values[1,index.non.zero]
    opt.cutoff <- signif(opt.cutoff,6)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,index.non.zero]
    p.value <- signif(p.value,3)
  }
  
  if(length(index.non.zero)!=1 && sum(index.non.zero)!=0){
    p.values <- q.p.values[2,index.non.zero]
    index.p <- which(p.values==min(p.values))
    cutoff.value.min <- min(q.p.values[1,index.p])
    opt.cutoff <- signif(cutoff.value.min,6)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,index.p]
    p.value <- signif(p.value,3)
  }
  
  if(sum(index.non.zero)==0){
    opt.cutoff <- q.p.values[1,1]
    opt.cutoff <- signif(opt.cutoff,6)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,1]
    p.value <- signif(p.value,3)        
  }
  
  if(sum(index.non.zero)==length(q)){
    opt.cutoff <- 0
    p.value <- 1       
  }
  
  #} 
  return(list(df=df, summary=summary, opt.cutoff=opt.cutoff, p.value=p.value))
}