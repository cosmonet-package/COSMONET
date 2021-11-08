#' @title Generate a gene-network to perform a KEGG pathway analysis on the genes selected by the Cox-network based model
#'
#' @description This function creates an interactive dashboard with two tabs.
#' The first tab "Not-Isolated Genes" shows a three-colour networks of the not-isolated genes, i.e. the genes sharing 
#' at least one KEGG pathway.
#' The second tab "Full Network" shows a three-colour networks of all the genes selected by the method (all the not zero signatures).
#' Different node colors are used to show a probability measure between the genes and the disease under investigation as follows:
# 'Red color is used for the top 500 tissue-specific genes according to the HumanBase tool (https://hb.flatironinstitute.org)  (mapped up genes),
#' Blue color is used for the bottom 501 tissue-specific genes (mapped down genes), 
#' White color is used for genes that are not listed as tissue-specific in the HumanBase database (not mapped genes). 
#' @param ListOfGenes a two-column data frame composed by gene names and beta coefficients
#' @param ListHeader true if the genes' list has a header, false otherwise
#' @param diseaseID the disease ID identified by Disease Ontology (DO) (https://disease-ontology.org)
#' @param nodesCol the colours of the nodes. Default values are orange, greeen and purple for mapped-up, mapped down and not mapped genes, respectively.
#' 
#' @return Create KEGG network and open the dashboard on the browser.
#' @export
GenerateNetwork <- function(ListOfGenes = NULL, header = TRUE, diseaseID = NULL, nodesCols = c('#EFCFD4','#D9E3FC','#FFFFFF')){
  
  # source("SupportingFunctionsGenerateNetwork.R")
  ## call rmarkdown from script with parameters
  # rmarkdown::render("Network_Dashboard/RUNCosmonetDashboard.Rmd", params = list(inputFile = ListOfGenes, #list of genes and beta
  #                                                                               headerFile = header,# does the list have the header?
  #                                                                               diseaseID = diseaseID, # disease ID 
  #                                                                               nodesColours = nodesCols)) #Nodes Colours
  ## open the dashboard on the browser
  # browseURL("Network_Dashboard/RUNCosmonetDashboard.html")   
  
  rmarkdown::render("RUNCosmonetDashboard.Rmd", params = list(inputFile = ListOfGenes,                     # two-column file with gene names and beta coefficients
                                                              headerFile = header,       # does the list have the header?
                                                              diseaseID = diseaseID,     # disease ID 
                                                              nodesColours = nodesCols)) # Nodes Colours
  # open the dashboard on the browser
  browseURL("RUNCosmonetDashboard.html")   
  
  
}

