
#' @title Select optimal cutoff on the training set
#'
#' @description This function computes the optimal cutoff \eqn{PI^{*,T}} on \eqn{T}, i.e., the value that corresponds to the best separation
#' in high-and-low risk group with respect to the log-rank test using the training set \eqn{T}.
#' @param x1 input training matrix \eqn{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param beta regression coefficients.
#'
#' @return The following objects are returned:
#' \item{PI.train}{prognostic index \eqn{PI} on the training set \eqn{T}. Compute the quantile \eqn{q_{\gamma}} of \eqn{PI} with \eqn{\gamma = \{0.20, 0.25, 0.30,... 0.80\}}. Each patient \eqn{i} in \eqn{T} is assigned to the high-risk (or low-risk) group if its prognostic index \eqn{PI_i} is above (or below) the \eqn{q_{\gamma}}-quantile.}
#' \item{cutoff.opt}{optimal cutoff selected adaptively on the training set \eqn{T}. It corresponds to the best separation in high-and-low risk group with respect to the log-rank test.}
#' \item{p.value}{\eqn{p}-value on the training set (the significance level is \eqn{p}-value < 0.05).}
#' @export
SelectOptimalCutoff <- function(x1, y1, beta){

  library(survival)
  # if(sum(beta)!=0){
    PI.train <- x1 %*% beta # [nxp]*[px1]=[nx1]

    # Compute the quantile q_gamma, with gamma = 0.2, ..., 0.80
    probs=seq(0.2 ,0.8, by=0.05)
    perc=paste(probs*100, "%", sep="")
    q <- quantile(PI.train, probs=probs)
    p <- vector()
    #par(mfrow=c(3,3))
    
    for (j in 1:length(q)){
      group <- vector(mode = "numeric", length = dim(x1)[1]);
      for (i in c(1:nrow(PI.train))){
        if (PI.train[i]>=q[j]) {group[i]=1}
        else if (PI.train[i]<q[j]) {group[i]=2}
      }

     if((length(unique(group)) > 1)==TRUE){
        #print(table(groupRisk)
        x <- data.frame(x1,groupRisk=group)
        logranktest <- survdiff(Surv(y1[,1], y1[,2]) ~ groupRisk, data = x, rho = 0)
        p.quantile <- 1-pchisq(logranktest$chisq, 1)
        p[j] <- as.vector(signif(p.quantile,3))
        # fit <- survfit(Surv(y1[,1], y1[,2]) ~ groupRisk, data = x)
        # cutoff <- signif(q[j],3)
      } else {
        p <- rep(1, length(q))  # no splitting - only one group
        }
      }
      
     q.p.values <- rbind(q,p)
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