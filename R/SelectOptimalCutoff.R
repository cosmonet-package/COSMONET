
#' @title Determine the optimal cutpoint and draw survival curves
#' 
#' @description This function computes the optimal cutoff \code{PI^{*,T}} on the training set \code{T}, i.e., the value that corresponds to the best separation in high-and-low risk group with respect to the log-rank test using the prognostic index \code{PI^{T}}.
#' @param x input training matrix \code{nxp}.
#' @param y response variable, \code{y} should be a two-column data frame with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored. The rownames indicate the sample names ordered as the samples in the input testing matrix.
#' @param beta regression coefficients estimated on training set \code{T}.
#' @param optCutpoint a character string for choosing the optimal cutpoint on training set \code{T} based on prognostic index \code{PI^{T}}. Can choose `minPValue`, `median` and `survCutpoint`.
#'
#' @details We use `surv_cutpoint()` to determine the optimal cutpoint and  `ggsurvplot()` function to plot survival curves from `survminer` package.
#' 
#' @return The following objects are returned:
#' \item{df}{data frame about sample, prognostic index, time, status and group risk for each quantile \eqn{q_{gamma}}, with \eqn{\gamma=0.25, ..., 0.80}.}
#' \item{summary}{summary information about number of patients at risk, cutoff and p.value for each quantile.}
#' \item{opt.cutoff}{optimal cutoff selected.}
#' @export
SelectOptimalCutoff <- function(x,y,beta,optCutpoint=c("minPValue","median","survCutpoint")){
  
  library(survival)
  library(survminer)
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(ggpubr)
  library(miscset)
  
  # if(sum(beta)!=0){
  PI <- as.matrix(x) %*% beta # [nxp]*[p]=[n]
  PI <- data.frame(PI)
  PI <- merge(PI,y,by.x=0,by.y=0)
  PI.data <- PI[order(PI[,2]),]
  colnames(PI.data)[1] <- "sample"
  
  if(optCutpoint == "minPValue"){
    # Compute the quantile q_gamma, with gamma = 0.25, ..., 0.80
    probs <- seq(0.25 ,0.8, by=0.05)
    perc <- paste(probs*100, "%", sep="")
    quantile <- quantileCutoff(PI.data, probs, perc, plot=TRUE)
    df <- quantile$df
    summary <- quantile$summary
    
    if((length(which(summary$p.value!=0)) >= 1)==TRUE){
      p.value.min <- summary[which(summary$p.value!=0),]
      index.p <- which(p.value.min$p.value==min(p.value.min$p.value))
      opt.cutoff <- p.value.min[index.p,3]
      # p.value <- p.value.min[index.p,4]
    } else if((sum(summary$p.value) == 0)==TRUE){
      opt.cutoff <- summary[grep("50%",rownames(summary)),3]
      # p.value <- summary[grep("50%",rownames(summary)),4]
    } else if((sum(summary$p.value) == length(q))==TRUE){ 
      opt.cutoff <- 0
      # p.value <- 1
    }
  }
  
  if(optCutpoint == "median"){
    # Compute the quantile q_gamma, with gamma = 0.50
    probs <- 0.50
    perc <- paste(probs*100, "%", sep="")
    quantile <- quantileCutoff(PI.data, probs, perc, plot=FALSE)
    df <- quantile$df
    summary <- quantile$summary
    opt.cutoff <- summary[,3]
    # p.value <- summary[,4]
    # opt.cutoff <- q.p.values[1,]
    # p.value <- q.p.values[2,]
  }
  
  if(optCutpoint == "survCutpoint"){
    # Compute the quantile q_gamma using surv_cutpoint() function from survminer package
    rownames(PI.data) <- PI.data$sample
    PI.data <- PI.data[,-1]
    colnames(PI.data)[3] <- "event"
    res.cut <- surv_cutpoint(PI.data, time = "time", event = "event", variables = "PI")
    summary <- summary(res.cut)
    opt.cutoff <- res.cut$cutpoint[1,1]
    # print(plot(res.cut, "PI", palette = "npg"))
    res.cat <- surv_categorize(res.cut)
    colnames(res.cat)[3] <- "groupRisk"
    
    PI <- data.frame(PI=PI.data[,1])
    rownames(PI) <- rownames(PI.data)
    df <- merge(PI, res.cat, by.x=0, by.y=0)
    colnames(df)[c(1,4)] <- c("sample","status")
    df$groupRisk[df$groupRisk=="low"] <- 2
    df$groupRisk[df$groupRisk=="high"] <- 1
    df <- df[order(df$PI,decreasing = FALSE),]
    rownames(df) <- 1:dim(df)[1]
    fitTrain <- survfit(Surv(time, status) ~ factor(groupRisk), data = df)
  
    survp <- ggsurvplot(
      fitTrain,                  
      data = df, 
      palette = c("red", "blue"),
      conf.int = TRUE,  
      pval = TRUE,
      pval.size = 5,
      risk.table = TRUE,  
      ggtheme = theme_minimal(),  
      legend.title = paste0("Cutoff: ", signif(opt.cutoff,4)),  
      legend.labs = c("High Risk", "Low Risk")
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
    
    ggsurv <- ggarrange(survp$plot, survp$table, heights = c(2, 0.7), ncol = 1, nrow = 2, font.label=list(color="black",size=14))
    df$id <- 1:dim(PI)[1]
    df$groupRisk[which(df$groupRisk==1)]  <- "High Risk"
    df$groupRisk[which(df$groupRisk==2)]  <- "Low Risk"
    df$groupRisk  <- as.factor(df$groupRisk)
    d <- ggline(df,"id","PI",color = "groupRisk", palette = c("red", "blue")) +
      ylim(c(min(df$PI),max(df$PI))) +
      labs(x="Sample", y="Prognostic Index", color = "", font.legend =  list(size = 14, color = "black"))
    d <- d + theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14),
                   axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
                   legend.text = element_text(size = 14))
    
    ggsurvComb <- ggarrange(ggsurv, d, labels = c("A", "B"), ncol = 2, nrow = 1, font.label=list(color="black",size=14))
    # annotate_figure(gg, top = "Training set")
    print(ggsurvComb)
    
    df$groupRisk <- as.character(df$groupRisk)
    df$groupRisk[which(df$groupRisk=="High Risk")] <- 1
    df$groupRisk[which(df$groupRisk=="Low Risk")]  <- 2
  }
  
  return(list(df=df,summary=summary,opt.cutoff=opt.cutoff))
}

quantileCutoff <- function(PI, probs, perc, plot){
  
  q <- quantile(PI[,2], probs=probs)
  p <- vector()
  
  # We stratified patients into high- and low-risk subtypes based on PI
  highRisk <- NULL
  lowRisk <- NULL
  risk.matrix <- NULL
  plots.list <- list()
  for (j in 1:length(q)){
    # print(paste0("cutoff: ",q[j]))
    groupRisk <- vector(mode = "numeric", length = dim(PI)[1]);
    
    for (i in c(1:nrow(PI))){
      if (PI[i,2] >= q[j]){groupRisk[i]=1}
      else if (PI[i,2]< q[j]){groupRisk[i]=2}
    }
    
    PI$groupRisk <- groupRisk
    ind.high <- which(groupRisk == 1)
    HR <- length(ind.high)
    highRisk <- cbind(highRisk, HR)
    ind.low <- which(groupRisk == 2)
    LR <- length(ind.low)
    lowRisk <- cbind(lowRisk, LR)
    lowHighRisk <- rbind(lowRisk,highRisk)
    
    if((length(unique(groupRisk)) > 1)==TRUE){
      # print(table(groupRisk)
      logranktest <- survdiff(Surv(time, status) ~ factor(groupRisk), data = PI, rho = 0)
      p.quantile <- 1-pchisq(logranktest$chisq, 1)
      p[j] <- as.vector(signif(p.quantile,4))
      # print(paste0("p.value: ", p[j]))
      fitTrain <- survfit(Surv(time, status) ~ factor(groupRisk), data = PI)
      # cutoff <- signif(q[j],3)
    } else {
      p <- rep(1, length(q))  # no splitting - only one group
    }
    
    # Keplan-Meier curve
    survp <- ggsurvplot(
      fitTrain,                  
      data = PI,   
      palette = c("red", "blue"),             
      conf.int = TRUE,  
      pval = p[j],  
      pval.size = 5,
      risk.table = TRUE,  
      ggtheme = theme_minimal(),  
      legend.title = paste0("Cutoff: ", perc[j]),  # round(q[j],4)
      legend.labs = c("High Risk", "Low Risk")
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

    # par(mfrow = c(4, 3)) 
    ggsurv <- ggarrange(survp$plot, survp$table, heights = c(2, 0.7), ncol = 1, nrow = 2, font.label=list(color="black",size=14))
    plots.list[[j]] <- survp$plot
    # print(ggsurv)
    
    # Distribution plot of PI 
    PI$id <- 1:dim(PI)[1]
    PI$groupRisk[which(PI$groupRisk==1)]  <- "High Risk"
    PI$groupRisk[which(PI$groupRisk==2)]  <- "Low Risk"
    PI$groupRisk  <- as.factor(PI$groupRisk)
    d <- ggline(PI,"id","PI",color = "groupRisk", palette = c("red", "blue")) +
      ylim(c(min(PI$PI),max(PI$PI))) + #, name = paste0("Cutoff: ",round(q[j],6))) +
      labs(x="Sample", y="Prognostic Index", color = "")
    d <- d + theme(axis.text.y = element_text(size = 14), axis.title.y = element_text(size = 14),
                   axis.text.x = element_text(size = 14), axis.title.x = element_text(size = 14),
                   legend.text = element_text(size = 14))
    # print(d)
    
    ggsurvComb <- ggarrange(ggsurv, d, labels = c("A", "B"), ncol = 2, nrow = 1, font.label=list(color="black",size=14))
    # annotate_figure(gg, top = "Training set")
    print(ggsurvComb)
    risk.matrix <- cbind(risk.matrix, groupRisk)
  }
  
  if(plot==TRUE){
  # Plot all in one figure
  # ggplotGridA4(plots.list, "plotCutoff.pdf", ncol=3, nrow = 4)
  plots <- grid.arrange(grobs=plots.list, ncol=3)
  # plots <- do.call(marrangeGrob, c(plots.list, list(nrow = 4, ncol = 3)))
  # To save to file, here on A4 paper
  ggsave("multipage_plot.pdf", plots, width = 21, height = 29.7, units = "cm")
  }
  
  # Risk matrix
  risk.matrix <- data.frame(risk.matrix)
  rownames(risk.matrix) <- PI$sample
  colnames(risk.matrix) <- perc
  
  df <- merge(PI[,-5],risk.matrix,by.x=1,by.y=0)
  df <- df[order(df$PI,decreasing = FALSE),]
  rownames(df) <- 1:dim(df)[1]
  df <- df[,-grep("id", colnames(df))]
  lowHighRisk <- as.data.frame(lowHighRisk)
  rownames(lowHighRisk) <- c("Low Risk", "High Risk")
  colnames(lowHighRisk) <- as.character(perc)
  # print(lowHighRisk)
  q.p.values <- as.data.frame(rbind(q,p))
  rownames(q.p.values) <- c("cutoff","p.value")
  q.p.values[2,] <- as.numeric(q.p.values[2,])
  # print(q.p.values)
  summary <- data.frame(t(lowHighRisk),t(q.p.values))
  # print(summary)
  
  return(list(df=df,summary=summary))
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
