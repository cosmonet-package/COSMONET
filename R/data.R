#' @exportPattern "."
#' @name breastGEOData
#' @title Microarray data from GEO database
#'
#' @description Affymetrix HG-U133A microarrays of two different breast cancer datasets are downloaded from GEO database:
#' GSE2034 as training set and GSE2990 as testing set. The data are noromalized and annotated.
#' The `training data` is composed by gene expression profiles from total RNA of frozen tumour samples from 286 lymph-node-negative breast cancer patients.
#' The `testing data` is composed by gene expression profiles of of 189 invasive breast carcinomas.
#' The survival information include breast cancer relapse free survival (RFS): `time` (months) and `status` for each patients in the two datasets. `
#'
#' @docType data
#'
#' @usage data(breastGEOData)
#'
#' @format Objects of class `data frame`.
#'
#' @keywords datasets
#'
#' @references
#' Wang Y, Klijn JG, Zhang Y, Sieuwerts AM et al.
#' Gene-expression profiles to predict distant metastasis of lymph-node-negative primary breast cancer.
#' #' \emph{Lancet 2005 Feb 19-25;365(9460):671-9} ([DOI](https://doi.org/10.1016/S0140-6736(05)17947-1)). \cr \cr
#' Sotiriou C, Wirapati P, Loi S, Harris A et al. 
#' Gene expression profiling in breast cancer: understanding the molecular basis of histologic grade to improve prognosis
#' #' \emph{J Natl Cancer Inst 2006 Feb 15;98(4):262-72} ([DOI](https://doi.org/10.1093/jnci/djj052)). \cr \cr
#'
#' @source GEO Database, <https://www.ncbi.nlm.nih.gov/geo/>
#'
#' @examples
#' data(breastGEOData)
NULL
#' @name breastGDCData
#' @title Count data from GDC Data Portal
#'
#' @description
#' Count data composed by 1095 breast cancer patients. The data are filtered and nromalized by using `TCGAbiolinks` package.
#' The survival information include `time` (years) and `status` for each patients.
#'
#' @docType data
#'
#' @usage data(breastGDCData)
#'
#' @format Objects of class `data frame`.
#'
#' @keywords datasets
#'
#' @references Colaprico A, Silva TC, Olsen C, Garofano L, Cava C, Garolini D, et al. (2015).
#' TCGAbiolinks: An R/Bioconductor package for integrative analysis of TCGA data.\cr
#' \emph{Nucleic Acids Research} ([DOI](http://doi.org/10.1093/nar/gkv1507)). \cr \cr
#' Mounir M, Lucchetta M, Silva C T, Olsen C, Bontempi G, Chen X, et al.(2019).
#' New functionalities in the TCGAbiolinks package for the study and integration of cancer data from GDC and GTEx. \cr
#' \emph{PLoS computational biology, 15(3), e1006701.} \cr \cr
#'
#' @source GDC Data Portal, <https://portal.gdc.cancer.gov/>
#'
#' @examples
#' data(breastGDCData)
NULL
#' @name Repository tissue-specific genes
#' @title Repository tissue-specific genes
#'
#' @description A repository of prediction scores or posterior probabilities associated to each gene. The higher the probability the stronger the functional relation between the gene and the tissue cancer.
#'
#' @docType data
#'
#' @usage data(repositoryDisease)
#'
#' @format Objects of class `"data.frame"`.
#'
#' @keywords datasets
#'
#' @source https://hb.flatironinstitute.org
#'
#' @examples
#' data(repositoryDisease)
NULL
#' @name Tissue Specific functional interaction
#' @title Tissue Specific functional interaction
#'
#' @description Top Edges: the network is filtered to only include edges with evidence supporting a tissue-specific functional interaction
#' [entrez gene id 1][entrez gene id 2][posterior prob.]
#'
#' @docType data
#'
#' @usage data(tissueSpecificEdge)
#'
#' @format Objects of class `"data.frame"`.
#'
#' @keywords datasets
#'
#' @source https://hb.flatironinstitute.org/download
#'
#' @examples
#' data(tissueSpecificEdge)
NULL
#' @name KEGG pathways
#' @title KEGG pathways
#'
#' @description A data frame contains pathway maps representing the a-priori biological knowledge including molecular interaction, reaction and networks relations (KEGG, PPI, etc).
#'
#' @docType data
#'
#' @usage data(KEGGrepository)
#'
#' @format Objects of class `"data.frame"`.
#'
#' @keywords datasets
#'
#' @source https://www.genome.jp/kegg/pathway.html
#'
#' @examples
#' data(KEGGrepository)
NULL
