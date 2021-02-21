
#' @title Survival prediction via screening-network Cox methods
#'
#' @description This function performs two main step: (i) penalized Cox regression methods to select a subset of potential biomarkers by using the training set \eqn{T}; (ii) the validation test by using the testing set \eqn{D}.
#' @param x1 input training matrix \eqn{n1xp}.
#' @param y1 response variable, \code{y1} should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening.
#' @param family Cox proportional hazards regression model.
#' @param penalty penalty type. Can choose "Net" where \code{Omega} matrix is requested.
#' For penalty = "Net", the penalty is defined as
#' \eqn{\lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}}, 
#' where \eqn{L} is a Laplacian matrix calculated from \code{Omega}.
#' @param alpha ratio between \code{L1} and Laplacian for \code{Net}. Default is \code{alpha = 0.5}.
#' @param Omega adjacency matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix.
#' @param ncv number of cross-validation performed for tuning optimal parameters (number of folds by default is 5).
#'
#' @return An object of class `"COSMONET"` is returned, which is a list of the following outputs:
#' \item{fitTrain}{all values of the fitted model (see \code{NetworkCox()} function)}
#' \item{PI.train}{the prognostic index \eqn{PI} on the training set. Compute the quantile \eqn{q_{\gamma}} of \eqn{PI} with \eqn{\gamma = \{0.20, 0.25, 0.30,... 0.80\}}. Each patient \eqn{i} in \eqn{T} is assigned to the high-risk (or low-risk) group if its prognostic index \eqn{PI_i} is above (or below) the \eqn{q_{\gamma}}-quantile.}
#' \item{cutoff.opt}{the optimal cutoff selected adaptively on the training set \eqn{T}. It corresponds to the best separation in high-and-low risk group with respect to the log-rank test.}
#' \item{p.value}{the \eqn{p}-value computed with respect to the log-rank test on the training set (the significance level is \eqn{p}-value < 0.05).}
#' 
#' @details The first step is the variable screening of the data which aimed to reduce the number of variables for a 
#' large to a moderate scale. To this purpose, we assume that only a small number of these \eqn{p} variables is affecting 
#' the survival outcome. Therefore, we filter out variables that are considered not relevant for the disease under 
#' investigation. To this purpose, we consider three different types of variable screenings: biomedical 
#' screening (BMD-screening), data-driven screening (DAD-screening) and the fusion of biomedical and 
#' data-driven screening (BMD+DAD-screening). The second step is the application of penalized methods using 
#' the subset of screened variables \eqn{\{x_j,j \in \mathcal{I}\}} (where \eqn{\mathcal{I}} depends on the type screening performed) 
#' as new feature space to further remove not significant variables from the model. To assess the stability of 
#' the survival prediction we performed the five-fold cross-validation different times and we take as estimate 
#' the average value of \eqn{\lambda} and the corrisponding \eqn{\alpha}. These two parameters are used to fit 
#' the corresponding penalized Cox model and obtain the parameter estimate of \eqn{\beta_\mathcal{I}}. 
#' Survival analysis is performed using the Kaplan Meier curves after dividing the patients in 
#' two risk groups (high-and-low risk group) on the basis of the prognostic index \eqn{PI} computed with the gene signature.
#' The \eqn{p-value}, used to test the null hypothesis that the survival curves are identical vs. the alternative that 
#' the two groups have different survival, is calculated by using the \eqn{log-rank test}.
#' 
#' @references Iuliano, A., Occhipinti, A., Angelini, C., De Feis, I., and Liò, p. (2018).
#' Combining Pathway Identification and Breast Cancer Survival Prediction via Screening-Network Methods. \cr
#' \emph{Frontiers in genetics, 9, 206.}\cr \cr
#' Iuliano, A., Occhipinti, A., Angelini, C., De Feis, I., & Lió, P. (2016). 
#' Cancer markers selection using network-based Cox regression: A methodological and computational practice. \cr
#' \emph{Frontiers in physiology, 7, 208.}\cr \cr
#' @export
Cosmonet <- function(x1, y1, x2, y2, screenVars, family="Cox", penalty="Net", Omega=NULL, alpha=0.5, ncv=NULL, screening = c("BMD", "DAD", "BMD+DAD")){
    
    training <- CosmonetTraining(x1, y1, screenVars, family="Cox", penalty="Net", Omega=NULL, alpha=0.5, ncv=NULL)
    beta <- as.matrix(fit$fit$beta)
    opt.cutoff <- fit$opt.cutoff
    
    testing <- CosmonetTesting(x2, y2, screenVars, beta, opt.cutoff, screening)
    
    obj <- c(list(training=training, testing=testing))
    class(obj) = "COSMONET"
    obj
}