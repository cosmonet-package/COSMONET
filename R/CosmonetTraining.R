
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
#************************************************************************************
# Fit Network-regularized Cox regression models
#************************************************************************************
NetworkCox <- function(x1, y1, penalty="Net", Omega=NULL, alpha=0.5, ncv=NULL){
  
  set.seed(1234)
  count <- 0
  colnames(y1) <- c("time","status")
  
  lambdaOpt <- NULL
  bestlambda <- NULL
  for(j in 1:length(alpha)){
    # j=1
    print(paste0("alpha: ",alpha[j]))
    
    lambdavec <- NULL
    for(i in 1:ncv){
      print(paste("cv",i,sep=""))
      
      cv <- APML0::APML0(x1, y1, family="cox", penalty="Net", Omega=Omega, alpha = alpha[j], nlambda=100, nfolds = 5)
      lambda.opt <- cv$lambda.opt
      #lambda.min <- cv$lambda.min
      lambdavec <- cbind(lambdavec,lambda.opt)
      
      # x <- log(as.matrix(cv$fit[1]))
      # y <- as.matrix(cv$fit[2])
      # plot(x,y,"p",xlab = "Log(lambda)",ylab ="cvm",col = "red", pch = 20)
      # abline(v=log(cv$lambda.opt), col="black",lwd=3, lty=2)
      # title(paste("alpha:",alpha[j]," ", "cv:",i, sep=""))
      # Sys.sleep(0.1)
    }
    
    lambdaOpt <- rbind(lambdaOpt,lambdavec)
    # colnames(lambdaOpt) <- paste0("cv",1:ncv,sep = "")
    
    lambdaOpt.mean <- mean(as.numeric(lambdaOpt[j,]))
    bestlambda <- rbind(bestlambda, lambdaOpt.mean)
  }
  
  opt.pars <- cbind(alpha,bestlambda[,1])
  colnames(opt.pars)[2] <- c("lambdaOpt.mean")
  rownames(opt.pars) <- c(1:length(alpha))
  opt.pars <- data.frame(opt.pars)
  
  ind.lambda.min <- which.min(opt.pars[,2])
  opt.tuningPars <- opt.pars[ind.lambda.min,]
  
  # Fit optimal model
  fit <- APML0::APML0(x = x1, y = y1, family="cox", penalty="Net", Omega=Omega,lambda=opt.tuningPars$lambdaOpt.mean,alpha=opt.tuningPars$alpha)
  
  beta <- as.matrix(fit$Beta)
  rownames(beta) <- colnames(x1)
  indexBeta <- which(beta!=0)
  nonzerobeta <- beta[indexBeta,]
  # length(indexBeta)
  
  if(length(indexBeta)==0){count = count + 1;
  if(count > 5) {print(sprintf("More than 5 beta, file are null"));}
  }
  
  return(list(beta=beta, opt.tuningPars=opt.tuningPars, nonzerobeta=nonzerobeta))
}
#************************************************************************************
# Select optimal cutoff on the training set
#************************************************************************************
SelectOptimalCutoff <- function(x1, y1, beta){
  
  library(survival)
  # if(sum(beta)!=0){
  PI.train <- x1 %*% beta # [nxp]*[px1]=[nx1]
  PI.data <- data.frame(sample=rownames(PI),PI.score=PI)
  PI.train <- PI.data[order(PI.data$PI.score),]
  
  # Compute the quantile q_gamma, with gamma = 0.2, ..., 0.80
  probs <- seq(0.2 ,0.8, by=0.05)
  perc <- paste(probs*100, "%", sep="")
  q <- quantile(PI.train$PI.score, probs=probs)
  p <- vector()
  #par(mfrow=c(3,3))
  
  # We stratified patients into high- and low-risk subtypes based on PI
  highRisk <- NULL
  lowRisk <- NULL
  for (j in 1:length(q)){
    # print(paste0("cutoff: ",q[j]))
    group <- vector(mode = "numeric", length = dim(PI.train)[1]);
    
    for (i in c(1:nrow(PI.train))){
        if (PI.train[i] >= q[j]){group[i]=1}
        else if (PI.train[i]< q[j]){group[i]=2}
    }
    
    ind.high <- which(group == 1)
    HR <- length(ind.high)
    highRisk <- cbind(highRisk, HR)
    ind.low <- which(group == 2)
    LR <- length(ind.low)
    lowRisk <- cbind(lowRisk, LR)
    highLowRisk <- rbind(lowRisk,highRisk)
    
    if((length(unique(group)) > 1)==TRUE){
      #print(table(groupRisk)
      y1 <- merge(PI.train$sample,y1,by.x=1,by.y=0)
      x <- data.frame(y1,Risk=group)
      logranktest <- survdiff(Surv(y1[,1], y1[,2]) ~ Risk, data = x, rho = 0)
      p.quantile <- 1-pchisq(logranktest$chisq, 1)
      p[j] <- as.vector(signif(p.quantile,3))
      # print(paste0("p.value: ", p[j]))
      fitTrain <- survfit(Surv(y1[,1], y1[,2]) ~ Risk, data = x)
      # cutoff <- signif(q[j],3)
    } else {
      p <- rep(1, length(q))  # no splitting - only one group
    }
    
    library(survminer)
    # Keplan-Meier curve
    survp <- ggsurvplot(
      fitTrain,                  # survfit object with calculated statistics.
      data = x,                  # data used to fit survival curves.
      conf.int = TRUE,           # show confidence intervals for point estimaes of survival curves.
      pval = p[j],               # show p-value of log-rank test.
      risk.table = FALSE,       
      ggtheme = theme_minimal(),  # customize plot and risk table with a theme.
      legend.title = paste0("Cutoff: ",round(q[j],4)),
      legend.labs = c("High Risk", "Low Risk")
    )
    print(survp)
  }
  
  highLowRisk <- as.data.frame(highLowRisk)
  rownames(highLowRisk) <- c("High Risk", "Low Risk")
  colnames(highLowRisk) <- as.character(perc)
  # print(highLowRisk)
  
  q.p.values <- as.data.frame(rbind(q,p))
  rownames(q.p.values) <- c("cutoff","p.value")
  # print(q.p.values)
  
  summary <- data.frame(t(highLowRisk),t(q.p.values))
  print(summary)
  
  index.non.zero <- which(q.p.values[2,]!=0)
  
  if(length(index.non.zero)==1){
    opt.cutoff <- q.p.values[1,index.non.zero]
    opt.cutoff <- signif(opt.cutoff,3)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,index.non.zero]
    p.value <- signif(p.value,3)
  }
  
  if(length(index.non.zero)!=1 && sum(index.non.zero)!=0){
    p.values <- q.p.values[2,index.non.zero]
    index.p <- which(p.values==min(p.values))
    cutoff.value.min <- min(q.p.values[1,index.p])
    opt.cutoff <- signif(cutoff.value.min,3)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,index.p]
    p.value <- signif(p.value,3)
  }
  
  if(sum(index.non.zero)==0){
    opt.cutoff <- q.p.values[1,1]
    opt.cutoff <- signif(opt.cutoff,3)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,1]
    p.value <- signif(p.value,3)        
  }
  
  if(sum(index.non.zero)==length(q)){
    opt.cutoff <- 0
    p.value <- 1       
  }
  #} 
  return(list(PI.train=PI.train,opt.cutoff=opt.cutoff,p.value=p.value))
}
#*********************************************************************************************************************
# Run CosmonetTraining() function
#*********************************************************************************************************************
CosmonetTraining <- function(x1, y1, screenVars, family="Cox", penalty="Net", Omega=NULL, alpha=0.5, ncv=NULL){
  
  if(is.list(screenVars)==TRUE){
  index <- match(unlist(screenVars),colnames(x1))
  } else {index <- match(screenVars,colnames(x1))}
  
  if(length(which(is.na(index)))==0){indexScreen=index
  } else {indexScreen <- index[-which(is.na(index))]}
  
  penalty <- match.arg(penalty)
  
  if(family=="Cox"){
    fitTrain <- switch(penalty,
                       "Net"=NetworkCox(x1[,indexScreen], y1, penalty="Net", Omega[indexScreen,indexScreen], alpha, ncv))

    fitTrain$family <- "cox"
  }
  
  ## Compute the optimal cutoff on $T$
  x1.screened <- x1[,indexScreen]
  index.non.zero.beta <- which(fitTrain$beta!=0)
  select.cutoff <- SelectOptimalCutoff(x1.screened[,index.non.zero.beta], y1, fitTrain$beta[index.non.zero.beta])
  opt.cutoff <- select.cutoff$opt.cutoff
  p.value <- select.cutoff$p.value
  PI.train <- select.cutoff$PI.train
  
  return(list(fitTrain=fitTrain, PI.train=PI.train, opt.cutoff=opt.cutoff, p.value=p.value))
}