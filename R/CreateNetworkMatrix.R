#' @title Create network matrix by using the cancer-specific genes
#'
#' @description This function creates a network matrix with zero diagonal and non-negative off-diagonal used to calculate Laplacian matrix and of dimention \code{pxp}.
#' Since this is a cancer-specific study, the tissue-specific network has been downloaded. 
#' The example below is done with the top-edges Network file using the files available at https://hb.flatironinstitute.org/download.
#' 
#' @param x input matrix \code{nxp}. Each row is an observation vector.
#' @param repositoryDisease the repository disease is a repository where the prediction score associated to each gene is a posterior probability. The higher the probability the stronger the functional relation between the gene and the cancer.
#' @param diseaseID the disease ID from disease ontology (http://disease-ontology.org).
#' @param topGenes the number of the top functional genes associated with the disease.
#' @param tissueSpecificEdge file downloaded from https://hb.flatironinstitute.org/download. Top Edges: the network is filtered to only include edges with evidence supporting a tissue-specific functional interaction.
#' [entrez gene id 1][entrez gene id 2][posterior prob.]
#' @param output if TRUE the Omega matrix is saved as .txt file in the current directory.
#' @param message if TRUE returns a message showing the number of edges in the graph.
#' 
#' @return The following objects are returned:
#' \item{DiseaseGenes}{top cancer-disease genes and relative posterior probabilities.}
#' \item{Omega}{network matrix with zero diagonal and non-negative off-diagonal \code{pxp}.}
#' @export
CreateNetworkMatrix <- function(x, repositoryDisease, diseaseID=NULL, topGenes=NULL, tissueSpecificEdge, output = TRUE, message = TRUE){
  
  DiseaseGenes <- DiseaseGenesSelection(repositoryDisease, diseaseID, topGenes)
  DiseaseGenesEntrezID <- as.character(DiseaseGenes[,1])
  
  # Find the index of the elements that are in the DiseaseGenesEntrezID list and in the cancer data list the instruction below returns all the indices of the elements from the DiseaseGenes list in the cancerData table
  idx <- which(tissueSpecificEdge[,1] %in% DiseaseGenesEntrezID & tissueSpecificEdge[,2] %in% DiseaseGenesEntrezID)
  networkDataSubset <- tissueSpecificEdge[idx,]
  print(sprintf("Your list contains %d genes. %d out of %d are used in the network.", length(DiseaseGenesEntrezID), length(DiseaseGenesEntrezID), 
                length(unique(union(networkDataSubset[,1], union(networkDataSubset[,2], DiseaseGenesEntrezID))))))
  
  # To generate the adjacency matrix, we use igraph
  library(igraph)
  # generate the graph
  g <- graph.data.frame(networkDataSubset,directed=FALSE)
  g <- set_edge_attr(g, "weight", value= networkDataSubset[,3])
  # create the adjMatrix - this is the FL - symmetrical
  NetworkMatrix <- get.adjacency(g,attr='weight',sparse=FALSE)
  
  # Translate NetworkMatrix into a zero-one matrix using the Threshold FL_th
  # As option: 0/1 matrix with FLth = 0.5
  # NetworkMatrix[NetworkMatrix >= FLth] <- 1
  # NetworkMatrix[NetworkMatrix < FLth] <- 0
  
  if(message == TRUE){
  #   # how many edges are left with this threshold?
     print(sprintf("There are %d edges left.", sum(NetworkMatrix !=0)/2))
   }
  
  symbol <- DiseaseGenes[match(rownames(NetworkMatrix),DiseaseGenesEntrezID),2]
  colnames(NetworkMatrix) <- symbol
  rownames(NetworkMatrix) <- symbol
  NetworkMatrix <- NetworkMatrix[intersect(rownames(NetworkMatrix),colnames(x)),intersect(colnames(NetworkMatrix), colnames(x))]
  # dim(NetworkMatrix)
  # NetworkMatrix[1:10,1:10]
  
  # Complete the network matrix with 0 blocks for genes that are not covered by the functional map
  noncommongenes <- setdiff(colnames(x),colnames(NetworkMatrix))
  zero_matrix <- matrix(0,length(noncommongenes),length(noncommongenes))
  rownames(zero_matrix) <- noncommongenes
  colnames(zero_matrix) <- noncommongenes
  library(magic)
  Omega <- adiag(NetworkMatrix,zero_matrix)
  dim(Omega)
  
  if(output == TRUE){
    # save Network Matrix 
    write.table(Omega, "Omega.txt", append = FALSE, sep = " ", dec = ".",
                row.names = TRUE, col.names = TRUE)  
  }
  return(list(DiseaseGenes=DiseaseGenes,Omega=Omega))
}

DiseaseGenesSelection <- function(repositoryDisease, diseaseID, topGenes=NULL){
  
  data <- repositoryDisease[repositoryDisease[,2]==diseaseID,]
  subdataTop <- data[1:topGenes,]
  
  # Convert entrezID in symbol gene name
  require("biomaRt")
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(
    mart=mart,
    attributes=c("entrezgene_id","external_gene_name"),
    filter = "entrezgene_id",
    values = subdataTop[,1])
  
  data <- merge(annotLookup, subdataTop, by.x=1, by.y=1)
  data <- data[!duplicated(data[,1]),c(1:4)]
  colnames(data) <- c("entrezgene_id","external_gene_name", "DOID accession", "probability score")
  dataOrd <- data[order(data[,4], decreasing = T),]
  
  #Save list as an external file 
  write.table(dataOrd, paste0(diseaseID,"_selectedGenes.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(dataOrd)
  
}



