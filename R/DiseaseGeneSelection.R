#' @title Selection of cancer-related genes
#'
#' @description This function uses the file disease.partaa.gz to select the disease under investigation - available from https://humanbase-dev.s3-us-west-2.amazonaws.com/humanbase-data/data/predictions/disease.partaa.gz. 
#' It returns a list of tissue-specific genes, the prediction scores as posterior probabilities and the disease IDs from disease ontology (http://disease-ontology.org).
#' and the gene IDs are Entrez.
#' @param repositoryDisease the repository disease where the prediction score associated to each gene is a posterior probability. The higher the probability the stronger the functional relation between the gene and the cancer.
#' @param diseaseID the disease ID from disease ontology (http://disease-ontology.org)
#' @param topGenes the number of the top functional genes associated with the disease - by defult is 500
#'
#' @return The following object is returned:
#' a data frame of the disease genes which is saved as .txt file in the current directory.
#' @export
DiseaseGenesSelection <- function(repositoryDisease, diseaseID, topGenes=500){
  
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
  data <- data[!duplicated(data[,1]),]
  dataOrd <- data[order(data[,4], decreasing = T),]
  
  #Save list as an external file
  write.table(dataOrd, paste0(diseaseID,"_selectedGenes.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(dataOrd)
  
}