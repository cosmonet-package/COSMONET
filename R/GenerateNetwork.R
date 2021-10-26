##' @title Generate a three-colour network to investigate genes selected by the network-based Cox model
#'
#' @description This function creates a three-colour network using a list of genes as input (.txt file).
#' Different node colors are used to show a correlation measure between the genes and the disease under investigation as follows:
#' Orange nodes represent the genes that are reported by the Functional Map of the Human Genome (Huttenhower, et. al 2008) as strongly correlated to the analysed cancer (p-value <= 0.05);
#' Green nodes represent the genes that are reported by the Functional Map of the Human Genome as weakly correlated to the analysed cancer (p-value > 0.05);
#' Purple nodes represent the genes that are not explored by the Functional Map of the Human Genome.
#' @param ListOfGenes input list of genes.
#' @param ListHeader true if the genes' list has a header, false otherwise
#' @param outputFile file name for the output network file
#' @param userWidth the width of the graphics region in inches. 
#' @param userHeight the height of the graphics region in inches. 
#' @param includeLegend if TRUE the legend will be included in the same pdf file, otherwise a new pdf file with the legend will be generated
#' to run this file: GenerateNetwork()
#' @export
GenerateNetwork <- function(ListOfGenes=listOfGenes, ListHeader = TRUE, outputFile = "NetworkFile.pdf", userWidth = 12, userHeight = 6, includeLegend = TRUE){

genes_Pathways_matrix <- from_data_to_pathways(ListOfGenes, ListHeader)
network <- create_network(genes_Pathways_matrix$connected_genes, genes_Pathways_matrix$tot_paths, repos$gname, repos$gvalue)

pdf(file = outputFile, width=userWidth, height = userHeight)
par(mfrow=c(1,2))
plot_network(network)
if(includeLegend){
  plot_legend(network$legend)
  dev.off()
  print(sprintf("The network has been saved in the file %s", outputFile));
}
else{
  #close connection to NetworkFile.pdf
  dev.off()
  #create a separate LegendFile
  pdf(file = "NetworkLegend.pdf")
  plot_legend(network$legend)
  dev.off()
  print("The legend of the Network has been saved in the file %s NetworkLegend.pdf");
}
}
