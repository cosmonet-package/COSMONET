#' @title Create 0/1 matrix by using the cancer-specific genes
#'
#' @description This function creates the 0/1 network using the files available from https://hb.flatironinstitute.org/download
#' Since this is a cancer-specific study, the tissue-specific network has been downloaded. 
#' The example below is done with the mammary epithelium top-edges Network file.
#' 
#' @param DiseaseGenes set of genes from the DiseaseGenesSelection() function or as input from the user - as entrez gene ids
#' @param tissueSpecificEdge file downloaded from https://hb.flatironinstitute.org/download. Top Edges: the network is filtered to only include edges with evidence supporting a tissue-specific functional interaction
#' [entrez gene id 1][entrez gene id 2][posterior prob.]
#' @param FLth FL adaptive threshold, default is 0.5 - edges greater than 0.5 are kept in the final graph
#' @param output if TRUE the Omega matrix is saved as .txt file in the current directory
#' @param message if TRUE returns a message showing the number of edges in the graph with the selected FLth threshold
#' 
#' @return The following object isreturned:
#' \item{0/1}{adjcent matrix} of dimention equal to the number of genes selected as tissue-specific genes
#' @export
CreateAdjMatrix <- function(DiseaseGenes, tissueSpecificEdge, FLth = 0.5, output = TRUE, message = TRUE){

  # Step 1: find the index of the elements that are in the DiseaseGenes list and in the cancer data list
  # the instruction below returns all the indices of the elements from the DiseaseGenes list in the cancerData table
  idx <- which(tissueSpecificEdge[,1] %in% DiseaseGenes & tissueSpecificEdge[,2] %in% DiseaseGenes)
  networkDataSubset <- tissueSpecificEdge[idx,]
  print(sprintf("Your list contains %d genes. %d out of %d are used in the network.", length(DiseaseGenes), length(DiseaseGenes), 
                length(unique(union(networkDataSubset[,1], union(networkDataSubset[,2], DiseaseGenes))))))
  
  # to generate the adjacency matrix, we use igraph
  library(igraph)
  # generate the graph
  g <- graph.data.frame(networkDataSubset,directed=FALSE)
  g <- set_edge_attr(g, "weight", value= networkDataSubset[,3])
  # create the adjMatrix - this is the FL - symmetrical
  adjMatrix <- get.adjacency(g,attr='weight',sparse=FALSE)
  # translate adjMatrix into a zero-one matrix using the Threshold FL_th
  adjMatrix[adjMatrix >= FLth] <- 1
  adjMatrix[adjMatrix < FLth] <- 0
  
  if(message == TRUE){
    # how many edges are left with this threshold?
    print(sprintf("With the threshold of %f, there are %d edges left.", FLth, sum(adjMatrix == 1)/2))
  }
  
  if(output == TRUE){
    #save adjMatrix 
    write.table(adjMatrix, "adjmSubMatrix.txt", append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)  
  }

  return(adjMatrix)
  
}



