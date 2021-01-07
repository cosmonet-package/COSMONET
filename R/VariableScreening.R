##' @title Fit network-regularized Cox regression models on screened data 
#'
#' @description This function implemets network-regularized Cox regression models which are iteratively trained and validated on different screened sets.
#' @param x1 input training matrix \eqn{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param x2 input testing matrix \eqn{n2xp}.
#' @param y2 response variable, \code{y2} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from DAD-and BMD+DAD-screening by using `ScreeningType()` function.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param alpha ratio between \code{L_1} and Laplacian for \code{Net}. Default is \code{alpha = 0.5}.
#' @param ncv number of cross-validation performed for tuning optimal parameters (number of folds by default is 5).
#' @param name type of screening used from `ScreeningType()` function.
#'
#' @return The following objects are returned:
#' \item{summary}{a data frame composed by number of screened variables (`nGenes`), optimal tuning parameters (`alpha` and `lambda`), number of selected genes (`nonZeroBeta`) and the \eqn{p}-value on the testing set (`pValue`).}
#' \item{beta}{a data frame of gene symbols and relative regression coefficients which is also saved in a file `txt`.}
#' \item{PI.train and PI.test}{prognostic index \eqn{PI} computed on the training \eqn{T} and testing set \eqn{D}.}
#' \item{survival curves}{plot of survival curves.}
#' \item{summary plots}{two scatter plots: `nGenes vs pValue` and `nGenes vs nonZeroBeta` only if more threshods (`nGenes`) are selected}
#' 
#' @export
VariableScreening <- function(x1,y1,x2,y2,screenVars,Omega,alpha,ncv,name){
    
    betaCoeff <- list()
    opt.tuningPars <- NULL
    nonzerobeta <- NULL
    ngenes <- NULL
    pvals <- NULL

    for (t in 1:length(screenVars)){
      th <- length(screenVars[[t]])
      print(paste0("Top genes: ", th))
      
      # Training phase
      fitTraining <- CosmonetTraining(x1,y1,screenVars[[t]],family="Cox",penalty="Net",alpha=alpha,Omega=Omega,ncv=ncv)
      
      # Regression coefficents
      beta <- fitTraining$fit$beta
      write.table(beta, paste0("BetaCoeff_",name,"_",th,".txt"), sep = "\t", row.names = T, col.names = F)
      betaCoeff[[t]] <- beta
      
      # Optimal tuning parameters
      opt.cutoff <- fitTraining$opt.cutoff
      optPars <- fitTraining$fitTrain$opt.tuningPars
      opt.tuningPars <- rbind(opt.tuningPars, optPars)
      
      # Prognostic index
      PI.train <- fitTraining$PI.train
      
      # Testing phase (Control that the length of the beta vector is not null!)
      if(sum(beta)!=0){
        fitTesting <- CosmonetTesting(x2,y2,as.character(screenVars[[t]]),beta,opt.cutoff,screening = paste0(name,"_",th))
        
        nonzero <- length(which(beta!=0))
        p.value <- fitTesting$p.value
        
        # Survival risk-groups
        survGroups <- fitTesting$HL.groups
        print(paste0("High-risk group: ",length(survGroups$risk[which(survGroups$risk=="High")])))
        print(paste0("Low-risk group: ",length(survGroups$risk[which(survGroups$risk=="Low")])))
        
        # Prognostic index
        PI.test <- fitTesting$PI.test
        
        } else {
          
        print("Warning: No genes have been selected!")
        nonzero <- 0
        p.value <- 1
        
      }

      ngenes <- rbind(ngenes,length(screenVars[[t]]))
      nonzerobeta <- rbind(nonzerobeta, nonzero)
      pvals <- rbind(pvals,p.value)
    }

    summary <- data.frame(nGenes=ngenes, opt.tuningPars, nonZeroBeta=nonzerobeta, pValue=pvals)
    rownames(summary) <- 1:dim(summary)[1]
    write.table(summary, paste0("summary_",name, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)  
    
    if(dim(summary)[1] > 1){
    filename <- sprintf(paste0("summaryPlot_",name, ".pdf"))
    pdf(filename, width = 14, height = 7)
    
    library(ggpubr)
    p1 <- ggscatter(summary, x = "nGenes",y = "nonZeroBeta", color = "#00AFBB", size = 3) +
      geom_line(color="#00AFBB") +
      geom_label(aes(nGenes+1,nonZeroBeta+1,label=paste0(nGenes,",",nonZeroBeta))) +
      coord_cartesian(xlim = c(length(screenVars[[1]]),th), ylim = c(0, max(summary$nonZeroBeta)+1))
    
    p2 <- ggscatter(summary, x = "nGenes", y = "pValue", fill = "pValue", size = 3, shape = 21) +
      geom_hline(yintercept=0.05, linetype="dashed", color = "green") +
      scale_fill_gradientn(colours=topo.colors(7),na.value = "transparent",breaks=c(0,0.5,1),
                           labels=c(0,0.5,1),limits=c(0,1)) + 
      coord_cartesian(xlim = c(length(screenVars[[1]]),th), ylim = c(0, 1))
    
    plot <- ggarrange(p1,p2,labels = c("A", "B"), ncol = 2, nrow = 1)
    
    print(plot, newpage = FALSE)
    dev.off()
    }
    return(list(summary=summary, betaCoeff=betaCoeff, PI.train=PI.train, PI.test=PI.test))
    
}
