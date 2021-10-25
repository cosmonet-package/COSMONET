# This script contains the supplementary function corrPlots() to perform the correlation analysis between PIs.

corrPlot <- function(df1, df2, groupRisk1, groupRisk2, name1, name2, set){
  
  df <- merge(df1, df2, by.x = 1, by.y = 1)
  df$condition <- ""
  df$condition[df[,groupRisk1]==2 & df[,groupRisk2]==2] <- "High Risk"
  df$condition[df[,groupRisk1]==1 & df[,groupRisk2]==1] <- "Low Risk"
  df$condition[df[,groupRisk1]==1 & df[,groupRisk2]==2 | df[,groupRisk1]==2 & df[,groupRisk2]==1] <- "Not  identified"
  df$condition <- factor(df$condition, labels = c("High Risk","Not  identified","Low Risk"),
                         levels = c("High Risk","Not  identified","Low Risk"))
  
  library(ggpubr)
  corrPlot <- ggscatter(df, x = "PI.x", y = "PI.y",
                    color = "condition", shape = 20, size = 3, palette = c("#FC4E07","grey","#00AFBB"), # Points color, shape and size
                    add = "reg.line",                                                                   # Add regression line
                    add.params = list(color = "blue", fill = "lightgray"),                              # Customize reg. line
                    # conf.int = TRUE,                                                                  # Add confidence interval
                    cor.coef = TRUE,                                                                    # Add correlation coefficient
                    cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top", label.sep = "\n"),
                    cor.method = "spearman",
  ) + labs(x=paste0("PI_",name1), y=paste0("PI_",name2), color = set)
  return(corrPlot)
}

