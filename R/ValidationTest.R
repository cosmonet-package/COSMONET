
#' @title Prognostic index (PI) and Kaplan-Meier curves on testing set
#'
#' @description This function compute the prognostic index \code{PI^{D}} on the testing set \code{D} using the regression coefficients and the optimal cutoff computed on the training set \code{T}. Then, the kaplan-Meier curves are shown as resulting of the log-rank test between high- and low-risk group.
#' @param x input testing matrix \code{nxp}.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param beta regression cofficients estimated on the training set \code{T}.
#' @param opt.cutoff optimal cutoff selected adaptively on the training set \code{T}.
#'
#' @return The following objects are returned:
#' \item{df}{data frame about sample, prognostic index, time, status and group risk on testing set \code{D}.}
#' \item{p.value}{\code{p}-value resulting from the log-rank test (the significance level is \code{p}-value < 0.05).}
#' \item{plots}{survival curves and distribution plot of prognostic index \code{PI^{D}}.}
#' @export
ValidationTest <- function(x,y,beta,opt.cutoff){
  
  library(survival)
  library(survminer)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(ggpubr)
  
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
    p.value <- signif(p.quantile,4)
    
    # Kaplan-Meier survival curves
    fitTest <- survfit(Surv(time, status) ~ factor(groupRisk), data = PI)
    # filename <- sprintf("survPlot_%s.pdf",paste0(screening))
    # ggsave(filename, print(survp))
    # pdf(filename)
    
    survp <- ggsurvplot(
      fitTest,                  
      data = PI,   
      palette = c("red", "blue"),
      conf.int = TRUE,           
      pval = p.value, 
      pval.size = 5,
      risk.table = TRUE,  
      ggtheme = theme_minimal(), 
      legend.title = "",
      legend.labs = c("High Risk", "Low Risk"),
      ) 
    
    # Changing the font size, style and color
    survp <- customize_labels(
      survp,
      font.x        = c(14, "plain"),
      font.y        = c(14, "plain"),
      font.xtickslab = c(14, "plain"),
      font.ytickslab = c(14, "plain"),
      font.legend = c(14, "plain"),
      font.legendlab = c(14, "plain")
    )
    
    # Font for Risk Table
    survp$table <- customize_labels(
      survp$table,
      font.x        = c(14, "plain"),
      font.y        = c(14, "plain"),
      font.xtickslab = c(14, "plain"),
      font.ytickslab = c(14, "plain"),
      font.legend = c(14, "plain"),
      font.legendlab = c(14, "plain")
    )
    
    # print(survp)
    # print(survp, newpage = FALSE)
    # dev.off()
    
    # Distribution plot of PI 
    PI$id <- 1:dim(PI)[1]
    PI$groupRisk[which(PI$groupRisk==1)]  <- "High Risk"
    PI$groupRisk[which(PI$groupRisk==2)]  <- "Low Risk"
    PI$groupRisk <- as.factor(PI$groupRisk)
    d <- ggpubr::ggline(PI,"id","PI", color = "groupRisk", palette = c("red", "blue")) +
      ylim(c(min(PI$PI),max(PI$PI))) + # , name = paste0("Optimal cutoff: ",opt.cutoff)
      labs(x="Sample", y="Prognostic Index", color = "") 
    d <- d + theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14),
                   axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
                   legend.text = element_text(size = 14))
    # print(d)
    
    ggsurv <- ggarrange(survp$plot, survp$table, heights = c(2, 0.7), ncol = 1, nrow = 2, font.label=list(color="black",size=14))
    ggsurvComb <- ggarrange(ggsurv, d, labels = c("A", "B"), ncol = 2, nrow = 1, font.label=list(color="black",size=14))
    # annotate_figure(gg, top = "Testing set")
    print(ggsurvComb)

  } else {p.value <- 1
  print("Warning: no splitting!")} # no splitting - only one group
  
  df <- PI[,-grep("id",  colnames(PI))]
  df$groupRisk <- as.character(df$groupRisk)
  df$groupRisk[which(df$groupRisk=="High Risk")] <- 1
  df$groupRisk[which(df$groupRisk=="Low Risk")]  <- 2
  
  return(list(df=df, p.value=p.value))
}

# Helper function to customize plot labels
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL, 
                              font.legend = NULL, font.legendlab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      if (!is.null(font.legend)) p <- p + theme(legend.title = .set_font(font.legend))
      if (!is.null(font.legendlab)) p <- p + theme(legend.text = .set_font(font.legendlab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}
