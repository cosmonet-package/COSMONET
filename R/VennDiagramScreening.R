##' @title Venn Diagram screening-network analysis
#'
#' @description This function takes the list of high-risk genes selected by each screening-network method and creates a PNG Venn Diagram.
#' @param genes.bmd BMD-screened variables.
#' @param genes.dad DAD-screened variables.
#' @param genes.bmd.dad BMD+DAD-screened variables.
#' @param output logical; if TRUE a PNG file of venn diagram is saved.

#' @return A list of lists which contain the values assigned to each of the areas of a venn diagram is returned.
#' 
#' @export
VennDiagramScreening <- function(genes.bmd,genes.dad,genes.bmd.dad,output=TRUE){

print(paste0("BMD-genes: ", sum(genes.bmd[,2]!=0)))
print(paste0("DAD-genes number: ", sum(genes.dad[,2]!=0)))
print(paste0("BMD+DAD-genes number: ", sum(genes.bmd.dad[,2]!=0)))

bmd <- genes.bmd[which(genes.bmd[,2]!=0),1] 
dad <- genes.dad[which(genes.dad[,2]!=0),1]
bmd.dad <- genes.bmd.dad[which(genes.bmd.dad[,2]!=0),1]  

# Load library
library(VennDiagram)
library(tidyverse)
# Chart
venn.diagram(
  x = list(bmd, dad, bmd.dad),
  category.names = c("BMD-genes" , "DAD-genes" , "BMD+DAD-genes"),
  filename = paste0("venn_diagramm.png"),
  output=output,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("orange", '#21908dff', '#440154ff'),
  fill = c(alpha("orange",0.3), alpha('#21908dff',0.3), alpha('#440154ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("orange", '#21908dff', '#440154ff'),
  rotation = 1
)

x <- list(bmd, dad, bmd.dad)
overlap <- calculate.overlap(x)

return(list(overlap=overlap))
}