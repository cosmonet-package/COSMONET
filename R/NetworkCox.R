
#' @title Fit Network-regularized Cox regression models
#'
#' @description This function fits network-regularized Cox regression methods with net (L1 and Laplacian) penalty. See for more details \code{APML0} package.
#' @param x1 input training matrix \eqn{n1xp}
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param penalty penalty type. Can choose \code{Net} where \code{Omega} matrix is requested.
#' For \code{penalty = Net}, the penalty is defined as
#' \eqn{\lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}}, 
#' where L is a Laplacian matrix  calculated from \code{Omega}.
#' @param alpha ratio between \code{L_1} and Laplacian for \code{Net}. Default is \code{alpha = 0.5}.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param ncv number of cross-validation performed for tuning optimal parameters (number of folds by default is 5).
#'
#' @return The following objects are returned:
#' \item{beta}{a matrix of coefficients.}
#' \item{opt.tuningPars}{optimal tuning parameters.}
#' \item{nonzerobeta}{non-zero regression coefficients.}
#' 
#' @references 
#' Li, X., Xie, S., Zeng, D., Wang, Y. (2018). 
#' Efficient l0-norm feature selection based on augmented and penalized minimization. \cr
#' \emph{Statistics in medicine, 37(3), 473-486.}
#' Sun, H., Lin, W., Feng, R., and Li, H. (2014).
#' Network-regularized high-dimensional cox regression for analysis of genomic data.\cr
#' \emph{Statistica Sinica}.
#' Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011). 
#' Distributed optimization and statistical learning via the alternating direction method of multipliers. \cr
#' \emph{Foundations and Trends in Machine Learning, 3(1), 1-122.}
#' Friedman, J., Hastie, T., Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. \cr 
#'  \emph{Journal of Statistical Software, Vol. 33(1), 1.} \cr \cr
#' @export
NetworkCox <- function(x1, y1, penalty="Net", Omega=NULL, alpha=0.5, ncv=NULL){

 set.seed(2020)
 count <- 0;

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
