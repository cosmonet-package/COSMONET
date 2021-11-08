# This script contains the auxiliary functions and library for the GenerateNetwork Function
#
## Input: genes must be a two-column array with genes and beta coefficients
## Output: return the list of non-zero genes and the corresponding beta coefficients

library("igraph")
library("graph")
library("randomcoloR")
library("visNetwork")
library("RColorBrewer")

select_non_zero_genes <- function(ListOfGenes, header=FALSE){
  genes <- ListOfGenes # read.table(ListOfGenes, header = header)
  idx <- which(genes[,2]!=0)
  non_zero_genes <-  genes[idx,]
  print(sprintf("There are %d non-zero genes", length(idx)))
  return(non_zero_genes)
}

from_data_to_pathways <- function(genes, ListHeader){
  # column 1 contains the genes, column 2 contains the beta values
  genes <- genes[,1]
  repos <- load_repos()
  KEGG_rep <- repos$KEGG
  genes_Pathways_matrix <- create_genes_Pathways_matrix(genes, KEGG_rep)
}

load_repos <- function(){
  # load('data/KEGGrepository.RData')
  output <- list()
  output$KEGG <- KEGGrepository; # readMat("KEGG_pathway.mat")
  return(output)
}

### Build the Matrix genes x Pathways - tot_paths
create_genes_Pathways_matrix <- function(genes, KEGG_rep){
  tot_paths <- matrix(0,nrow = length(genes), ncol = length(KEGG_rep[1]$pathway.name)) 
  row.names(tot_paths) <- genes
  colnames(tot_paths) <- unlist(KEGG_rep$pathway.name)
  paths_name <- NULL
  for(i in 1:length(genes))
    for(j in 1:length(KEGG_rep[1]$pathway.name)){
      if(is.element(genes[i], unlist(KEGG_rep$pathway.gene.name[j]))){
        # print("i: ");print(i)
        # print("j: ");print(j)
        tot_paths[i,j] <- 1
        paths_name[j] <-  KEGG_rep$pathway.name[j]
      }
    }
  index <- which(colSums(tot_paths) > 1)
  connected_genes <- NULL #genes connected
  for (i in 1:length(index)){
    genes_involved <- which(tot_paths[,index[i]]>0)
    connected_genes <- union(connected_genes, genes[genes_involved])
  }
  print(sprintf("There are %i connected genes", length(connected_genes)))
  
  output <- list()
  output$tot_paths <- tot_paths
  output$connected_genes <- connected_genes
  output$paths_name <- paths_name
  
  return(output)
}

create_network <- function(connected_genes, tot_paths, gname, gvalue){
  # if there are no connected genes
  if(length(connected_genes)==0){
    print(sprintf("There are %i connected genes", length(connected_genes)))
    output <- list()
    output$network <- make_empty_graph(n=0)
    return(output)
  }
  else{
    df <- create_df(connected_genes, tot_paths)
    g <- graph_from_data_frame(as.matrix(df),vertices =connected_genes,directed = F)
    g <- set.vertex.attribute(g, "pvalue", value=0)
    g <- set.vertex.attribute(g, "colour", value="red")
    edge_info <- create_edges_colours(df)
    edge_labels <- create_edges_labels(df, edge_info)
    legend <- create_legend(df, edge_info)
    output <- list()
    output$network <- g
    output$df <- df
    output$edge_colour <- edge_labels$edge_colour
    output$edge_number <- edge_labels$edge_number
    output$edge_path <- edge_labels$edge_path
    output$legend <- legend
    return(output)
  }
  
}

define_node_colour <- function(listOfGenes, repositoryDisease, diseaseID, genesCol){
  n_colours <- rep(0,length(listOfGenes))
  n_borders <- rep("#2B7CE9",length(listOfGenes)) # default colour
  repos <- sorted_cancer_related_genes(repositoryDisease, diseaseID)
  gname <- repos$gene_name
  gvalue <- repos$probability
  # change nodes colour
  for(i in 1:length(listOfGenes)){
    if(is.element(listOfGenes[i],gname)==TRUE){
      m <- match(listOfGenes[i],gname)
      if(m<=500){
        n_colours[i] = genesCol[1];#'#ffd9b3' # is in the top 500 genes - red
        n_borders[i] = genesCol[1];
      }
      else{
        n_colours[i] = genesCol[2];'#d7fcb5' # is in the bottom 501 genes - blue
        n_borders[i] = genesCol[2];
      }
    }
    else{
      n_colours[i] = genesCol[3]; '#e4b5fc' # not in the list - white
      n_borders[i] = "#e6e4df"; #default colours from VisNetwork package is #2B7CE9
    }
  }
  n_format <- list()
  n_format$n_colours <- n_colours
  n_format$n_borders <- n_borders
  return(n_format)
}

create_df <- function(connected_genes, tot_paths){
  edge_added<- NULL
  from <- NULL
  to <- NULL
  
  sub_tot_paths <- tot_paths[connected_genes, which(colSums(tot_paths) > 1)]
  h = 1
  for(k in 1:dim(sub_tot_paths)[2]){
    for(i in 1:dim(sub_tot_paths)[1]){
      j = i+1
      while(j<=dim(sub_tot_paths)[1]){ 
        if(sub_tot_paths[i,k]==1 && sub_tot_paths[j,k]==1){
          from[h] <- connected_genes[i]
          to[h] <- connected_genes[j]
          edge_added[h] <- colnames(sub_tot_paths)[k]
          h = h +1;
        }
        j = j + 1
      }
    }
  }
  df <- data.frame(from, to, edge_added)
  return(df)
}

create_edges_colours <- function(df){
  pathways <- unique(df$edge_added)
  # generate random colours
  n <- length(pathways)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # remove colours that are too bright
  col_vector <- col_vector[-4]
  col_vector <- col_vector[-26]
  colors_random <- col_vector[1:n]
  edge_colours <- NULL
  for(i in 1: length(df$edge_added)){
    edge_colours[i] = colors_random[match(df$edge_added[i], pathways)]
  }
  output <- list()
  output$edge_colours <- edge_colours
  output$list_of_colours <- colors_random
  return(output)
}

create_edges_labels <- function(df, edge_info){
  pathways <- unique(df$edge_added)
  list_of_colours <- edge_info$list_of_colours
  edge_colour <- NULL
  edge_path <- NULL
  edge_number <- NULL
  # to store the number in the legend
  edge_number <- NULL
  for(i in 1: length(df$edge_added)){
    edge_colour[i] = list_of_colours[match(df$edge_added[i], pathways)]
    edge_path[i] = df$edge_added[i]
    edge_number[i] = which(pathways==df$edge_added[i]);
  }
  output <- list()
  # colours
  output$edge_colour <- edge_colour
  # pathways
  output$edge_path <-  edge_path
  # number in legend
  output$edge_number <-  edge_number
  return(output)
}

create_legend <- function(df, edge_info){
  legend_number <- length(unique(df$edge_added))
  legend_name <- unique(df$edge_added)
  legend <- list()
  legend$legend_number <- legend_number
  legend$legend_name <- legend_name
  legend$legend_colour <- edge_info$list_of_colours
  return(legend)
}

plot_legend <- function(legend){
  y_max = 500
  plot(1:40,xaxt='n', yaxt='n',ann=FALSE,axes=FALSE, col = "white",ylim = c(1,y_max), xlim = c(0,10))
  pos = y_max
  # cex_level is used to separate legend entries on different rows
  cex_level = 10
  # this is the font size of each entry
  cex_entry = 1
  # move x to the left
  x_coord = 0
  for(i in 1:legend$legend_number){
    if (i<10){
      text(x_coord,pos, i, col = legend$legend_colour[i], pos = 4, cex = cex_entry)
      text(x_coord+0.5,pos, legend$legend_name[i], col = legend$legend_colour[i], pos = 4, cex = cex_entry)
    }
    if(i>9 && i<100){
      text(x_coord,pos, i, col = legend$legend_colour[i], pos = 4, cex = cex_entry)
      text(x_coord+0.5,pos, legend$legend_name[i], col = legend$legend_colour[i], pos = 4, cex = cex_entry)
    }
    if(i>99){
      text(x_coord,pos, i, col = legend$legend_colour[i], pos = 4, cex = cex_entry)
      text(x_coord+0.5,pos, legend$legend_name[i], col = legend$legend_colour[i], pos = 4, cex = cex_entry)
    }
    pos = pos-cex_level
  }
  par(xpd=TRUE)
}

# Extract genes information from repositoryDisease.rda using a given DOID
# For example, diseaseID = "DOID:1612" for Breast Cancer
sorted_cancer_related_genes <- function(repositoryDisease, diseaseID){
  #load data with given diseaseID  
  # load("data/repositoryDisease.rda")
  data <- repositoryDisease[repositoryDisease[,2]==diseaseID,]
  # Convert entrezID in symbol gene name
  require("biomaRt")
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  annotLookup <- getBM(
    mart=mart,
    attributes=c("entrezgene_id","external_gene_name"),
    filter = "entrezgene_id",
    values = data[,1])
  
  # match genes names, remove duplicates and sort the data
  data <- merge(annotLookup, data, by.x=1, by.y=1)
  data <- data[!duplicated(data[,1]),]
  #extract only gene names and prob and sort them by probability
  dataOrd <- data[order(data[,4], decreasing = T), c(2,4)]
  colnames(dataOrd) <- c('gene_name', 'probability')
  rownames(dataOrd) <- NULL
  
  return(dataOrd)
}

######################################################################################################################################################################
# # takes data from the user
# load_data <-function(list_of_genes_file, header_true_false = TRUE){
#   genes <- read.table(list_of_genes_file, header = header_true_false) #your list goes here
#   genes <- as.matrix(genes)
#   return(genes)
# }
# 
#####################################################################################################################################################################
#
# create_network_data_frame <- function(genes, repository_name){
#   # Download from http://www.gsea-msigdb.org/gsea/msigdb/index.jsp
#   # c2.all.v7.4.symbols.gmt
#   # c2.cp.biocarta.v7.4.symbols.gmt
#   # c2.cp.kegg.v7.4.symbols.gmt
#   # c4.all.v7.4.symbols.gmt
#   # c5.go.v7.4.symbols.gmt
#   # c6.all.v7.4.symbols.gmt
#   # h.all.v7.4.symbols.gmt
#   # msigdb.v7.4.symbols.gmt
#   
#   library(qusage)
#   file <- "c2.all.v7.4.symbols.gmt"
#   list <- read.gmt(file) 
#   
#   my_genes <- genes
#   my_function <- function(x){
#     is.element(my_genes,x)
#   }
#   
#   my_sets <- list
#   # lapply(my_sets, my_function)
#   test <- as.data.frame(lapply(my_sets, my_function))
#   rownames(test) <- my_genes
#   test[test=="TRUE"] <- 1
#   test[test=="FALSE"] <- 0
#   
#   tot_paths <- test
#   #find index of pathways with more than one gene in it
#   pathways <- which(colSums(tot_paths) > 1)
#   connected_genes <- which(rowSums(tot_paths[,pathways])!=0)
#   connected_genes <- rownames(tot_paths)[connected_genes]
#   
#   network_df <- create_df(connected_genes, tot_paths)
#   
#   output <- list()
#   output$network_edges <- network_df
#   output$network_nodes <- connected_genes
#   
#   return(output)
# }
#
#####################################################################################################################################################################
#
# create_adj_matrix <- function(connected_genes, tot_paths){
#   
#   #connected_genes <- genes_Pathways_matrix$connected_genes
#   sub_tot_paths <- tot_paths[connected_genes, which(colSums(tot_paths) > 1)]
#   
#   adj_matrix <- array(0,dim=c(length(connected_genes),length(connected_genes)))
#   colnames(adj_matrix) <- connected_genes
#   row.names(adj_matrix) <- connected_genes
#   edge_added <- NULL
#   h = 0
#   for(k in 1:dim(sub_tot_paths)[2]){
#     for(i in 1:dim(sub_tot_paths)[1]){
#       #for(j in (i+1):dim(sub_tot_paths)[1]){
#       j = i+1
#       while(j<16){ 
#         #print(sprintf("k: %i, i: %i, j: %i", k, i, j))
#         if(sub_tot_paths[i,k]==1 && sub_tot_paths[j,k]==1){
#           adj_matrix[i,j] = 1
#           adj_matrix[j,i] = 1
#           edge_added[h] <- colnames(sub_tot_paths)[k]
#           #print(sprintf("edge_added: %s between node %s and %s", edge_added[h], connected_genes[i], connected_genes[j]))
#           h = h +1;
#         }
#         j = j + 1
#       }
#     }
#   }
#   print.default(adj_matrix)
#   return(adj_matrix)
# }
#
#####################################################################################################################################################################
#
# set_nodes_shapes <- function(g, v_colours){
#   v_shapes <- rep("circle",vcount(g))
#   v_number_of_papers <- rep("", vcount(g))
#   v_label_colours <- rep(rgb(0,0,0,0), vcount(g))
#   v_labels <- V(g)$name
#   v_border <- rep("#2B7CE9", vcount(g))
#   v_width <- rep(0, vcount(g))
#   v_title <- rep("", vcount(g))
#   g <- set.vertex.attribute(g, "number_of_papers", value=v_number_of_papers)
#   # Change shape if they are in papers
#   database <- read.delim("Find Papers related Tool/Database.csv", header = TRUE, sep = ",")
#   for(i in 1:vcount(g)){
#     if(is.element(V(g)$name[i],database[,1])){
#       m <- match(V(g)$name[i],database[,1])
#       v_shapes[i] = "rectangle" #"square"
#       v_number_of_papers[i] <- database[m,5] 
#       v_label_colours[i] <- "black"
#       v_labels[i] <- paste(V(g)$name[i], v_number_of_papers[i], sep="\n", collapse = "") 
#       v_border[i] = 'red'
#       v_width[i] = 2;
#       v_title[i] = sprintf("Node %s has been found in %s papers", V(g)$name[i], v_number_of_papers[i]);
#       #print(sprintf("The number of papers for node %s is: %s", v_labels[i], v_number_of_papers[i]))
#     }
#   }
#   output <- list()
#   output$v_shapes <- v_shapes
#   output$v_number_of_papers <- v_number_of_papers
#   output$v_colours <- v_colours
#   output$v_label_colours <- v_label_colours
#   output$v_labels <- v_labels
#   output$v_border <- v_border
#   output$v_width <- v_width
#   output$v_title <- v_title
#   #print('set nodes shape function - v_width')
#   #print(output$v_width)
#   return(output)
# }
#
# set_nodes_colous <- function(g, gname, gvalue){
#   v_colours <- rep(0,length(V(g)$name))
#   repos <- load_repos()
#   gname <- repos$gname
#   gvalue <- repos$gvalue
#   #change nodes colour
#   for(i in 1:length(V(g))){
#     if(is.element(V(g)$name[i],gname)==TRUE){
#       m <- match(V(g)$name[i],gname)
#       if(as.numeric(gvalue[m])<0.05){
#         v_colours[i] = '#ffd9b3' # is in the list with pvalue < 0.05 - orange
#       }
#       else{
#         v_colours[i] ='#d7fcb5' # in the list with pvalue >= 0.05 - green
#       }
#     }
#     if(is.element(V(g)[i]$name,gname)==FALSE){
#       v_colours[i] = '#e4b5fc' # not in the list - purple
#     }
#   }
#   return(v_colours)
# }
# 
#####################################################################################################################################################################
# 
# plot_network <- function(graph){
#   g <- graph$network
#   v_shapes_labels <- graph$features
#   groups <- cluster_optimal(g)
#   l <- layout_in_circle(g, order =
#                      order(membership(groups)))
#   # pdf("rplot.pdf", width=6, height = 8) 
#   plot(g, layout=l, vertex.shape = v_shapes_labels$v_shapes, vertex.label = v_shapes_labels$v_labels, vertex.size = 18, vertex.label.degree=-pi/2, 
#        vertex.label.dist=0, vertex.color = v_shapes_labels$v_colours, vertex.frame.color = rgb(0,0,0,0), vertex.label.color="black",
#        edge.curved = FALSE, edge.color = graph$edge_colours, edge.label = graph$edge_labels, edge.label.color=rgb(0,0,0,0))
#   # the instruction below adds the name of the pathway as edge_labels
#   plot(g, layout=l, vertex.frame.color=rgb(0,0,0,0), vertex.color=rgb(0,0,0,0), add=TRUE, #vertex.label = v_shapes_labels$v_number_of_papers,
#       vertex.label.color = rgb(0,0,0,0),#vertex.label.color=  v_shapes_labels$v_label_colours, vertex.label.dist=-1.5, vertex.label.degree=-pi/2,
#      edge.color = graph$edge_colours, edge.label = graph$edge_labels, edge.label.color = graph$edge_colours)
#   # the instruction below adds the number of the pathway as edge_labels
#   plot(g,layout=l, vertex.frame.color=rgb(0,0,0,0), vertex.color=rgb(0,0,0,0), add=TRUE, #vertex.label = v_shapes_labels$v_number_of_papers,
#        vertex.label.color = rgb(0,0,0,0),#vertex.label.color=  v_shapes_labels$v_label_colours, vertex.label.dist=-1.5, vertex.label.degree=-pi/2,
#        edge.color = graph$edge_colours, edge.label = graph$edge_numbers, edge.label.color = graph$edge_colours)
#   # par(new=TRUE)
#   # plot_legend(graph$legend)
#   # dev.off()
# }
# 
# 
# create_network_pdf <- function(graph, output = "destination.pdf"){
#   pdf(file = output, width=12, height = 6)
#   par(mfrow=c(1,2))
#   plot_network(network)
#   plot_legend(network$legend)
#   dev.off()
# }




# 
# define_node_border_and_title <- function(listOfGenes, nodesColours){
#   #v_shapes <- rep("circle",vcount(g))
#   #v_number_of_papers <- rep("", vcount(g))
#   #v_label_colours <- rep(rgb(0,0,0,0), vcount(g))
#   #v_labels <- V(g)$name
#   #v_border <- rep("#2B7CE9", length(listOfGenes))
#   #assign the default border colour
#   v_border <- rep("#808080", length(listOfGenes))#nodesColours
#   #v_width <- rep(0, vcount(g))
#   v_title <- rep("", length(listOfGenes))
#   #Change border colours if they are in papers
#   database <- read.delim("Find Papers related Tool/Database.csv", header = TRUE, sep = ",")
#   for(i in 1:length(listOfGenes)){
#     if(is.element(listOfGenes[i],database[,1])){
#       m <- match(listOfGenes[i],database[,1])
#       #v_shapes[i] = "rectangle" #"square"
#       #v_number_of_papers[i] <- database[m,5] 
#       v_border[i] = 'red'
#       v_title[i] = sprintf("Node %s has been found in %s papers", listOfGenes[i], database[m,5]);
#       #print(sprintf("The number of papers for node %s is: %s", v_labels[i], v_number_of_papers[i]))
#     }
#   }
#   output <- list()
#   #output$v_shapes <- v_shapes
#   #output$v_number_of_papers <- v_number_of_papers
#   #output$v_colours <- v_colours
#   #output$v_label_colours <- v_label_colours
#   #output$v_labels <- v_labels
#   output$v_border <- v_border
#   #output$v_width <- v_width
#   output$v_title <- v_title
#   return(output)
# }




