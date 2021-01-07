# This script contains the auxiliary functions and library for the GenerateNetwork() function

library("igraph")
library("graph")
library("randomcoloR")


load_repos <- function(listDiseaseGenes){
  pvalue = read.delim(listDiseaseGenes, header = FALSE)  
  pvalue = as.matrix(pvalue)
  dim(pvalue)
  data(KEGG_Repository)
  gname = pvalue[,3]
  gvalue = pvalue[,1]
  #given a set of genes, this script returns the number of papers associated with each gene
  #source('Find Papers related Tool/script.R', echo=TRUE)
  output <- list()
  output$gname <- pvalue[,3]
  output$gvalue <- pvalue[,1]
  output$KEGG <- KEGG_Repository;#readMat("KEGG_pathway.mat")
  return(output)
}

load_related_papers <- function(genes,listOfPaper){
  database <- read.delim(listOfPaper, header = TRUE, sep = ",")
  database <- as.matrix(database)
  #dim(database)
  aliases <- NULL#matrix(0, nrow = dim(database)[1], ncol = 10)
  max_length <- 0
  
  for(i in 1:dim(database)[1]){
    #print(i)
    temp_vect <- matrix(0,nrow = 1, ncol = 19)
    split <- strsplit(database[i,3], ",")[[1]]
    if(length(split)>0){
      temp_vect[1:length(split)] <- split
      aliases <- rbind(aliases, temp_vect)
      if(length(split)>max_length){max_length <- length(split)}
    }
    if(length(split)==0){aliases <- rbind(aliases, "")}
    
  }
  database_paper_genes <- intersect(genes, database[,1])
  print(sprintf("There are %i papers related to your list of genes", length(database_paper_genes)))
  return(database_paper_genes)
}

# takes data from the user
load_data <-function(list_of_genes_file, header_true_false = TRUE){
  genes <- read.table(list_of_genes_file, header = header_true_false) #your list goes here
  genes <- as.matrix(genes)
  return(genes)
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
        #print("i: ");print(i)
        #print("j: ");print(j)
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
  output$connected_genes<- connected_genes
  output$paths_name <- paths_name
  
  return(output)
}

create_adj_matrix <- function(connected_genes, tot_paths){
  
  #connected_genes <- genes_Pathways_matrix$connected_genes
  sub_tot_paths <- tot_paths[connected_genes, which(colSums(tot_paths) > 1)]
  
  adj_matrix <- array(0,dim=c(length(connected_genes),length(connected_genes)))
  colnames(adj_matrix) <- connected_genes
  row.names(adj_matrix) <- connected_genes
  edge_added <- NULL
  h = 0
  for(k in 1:dim(sub_tot_paths)[2]){
    for(i in 1:dim(sub_tot_paths)[1]){
      #for(j in (i+1):dim(sub_tot_paths)[1]){
      j = i+1
      while(j<16){ 
        #print(sprintf("k: %i, i: %i, j: %i", k, i, j))
        if(sub_tot_paths[i,k]==1 && sub_tot_paths[j,k]==1){
          adj_matrix[i,j] = 1
          adj_matrix[j,i] = 1
          edge_added[h] <- colnames(sub_tot_paths)[k]
          #print(sprintf("edge_added: %s between node %s and %s", edge_added[h], connected_genes[i], connected_genes[j]))
          h = h +1;
        }
        j = j + 1
      }
    }
  }
  print.default(adj_matrix)
  return(adj_matrix)
}

create_df <- function(connected_genes, tot_paths){
  edge_added<- NULL
  from <- NULL
  to <- NULL
  
  sub_tot_paths <- tot_paths[connected_genes, which(colSums(tot_paths) > 1)]
  h = 1
  for(k in 1:dim(sub_tot_paths)[2]){
    for(i in 1:dim(sub_tot_paths)[1]){
      #for(j in (i+1):dim(sub_tot_paths)[1]){
      j = i+1
      while(j<=dim(sub_tot_paths)[1]){ 
        #print(sprintf("k: %i, i: %i, j: %i", k, i, j))
        if(sub_tot_paths[i,k]==1 && sub_tot_paths[j,k]==1){
          #adj_matrix[i,j] = 1
          #adj_matrix[j,i] = 1
          from[h] <- connected_genes[i]
          to[h] <- connected_genes[j]
          edge_added[h] <- colnames(sub_tot_paths)[k]
          #print(sprintf("edge_added: %s between node %s and %s", edge_added[h], connected_genes[i], connected_genes[j]))
          h = h +1;
        }
        j = j + 1
      }
    }
  }
  df <- data.frame(from, to, edge_added)
  return(df)
}



create_network <- function(connected_genes, tot_paths, gname, gvalue){
  #adj_matrix <- create_adj_matrix(connected_genes, tot_paths)
  #g <- graph_from_adjacency_matrix(adj_matrix, mode = c("undirected"))
  df <- create_df(connected_genes, tot_paths)
  g <- graph_from_data_frame(as.matrix(df),vertices =connected_genes,directed = F)
  g <- set.vertex.attribute(g, "pvalue", value=0)
  g <- set.vertex.attribute(g, "colour", value="red")
  v_colours <- set_nodes_colous(g, gname, gvalue)
  v_shapes_labels <- set_nodes_shapes(g, v_colours)
  #g <- add_edges_labels(g, df)
  edge_info <- create_edges_colours(df)
  edge_labels <- create_edges_labels(df, edge_info)
  legend <- create_legend(df, edge_info)
  output <- list()
  output$network <- g
  output$features <- v_shapes_labels
  output$df <- df
  output$edge_colours <- edge_info$edge_colours
  output$edge_labels <- edge_labels
  output$legend <- legend
  #print(output$edge_colours)
  return(output)
}

set_nodes_shapes <- function(g, v_colours){
  v_shapes <- rep("circle",vcount(g))
  v_number_of_papers <- rep("", vcount(g))
  v_label_colours <- rep(rgb(0,0,0,0), vcount(g))
  v_labels <- V(g)$name
  g <- set.vertex.attribute(g, "number_of_papers", value=v_number_of_papers)
  #Change shape if they are in papers
  database <- read.delim("Find Papers related Tool/Database.csv", header = TRUE, sep = ",")
  for(i in 1:vcount(g)){
    if(is.element(V(g)$name[i],database[,1])){
      m <- match(V(g)$name[i],database[,1])
      v_shapes[i] = "square"
      v_number_of_papers[i] <- database[m,5] 
      v_label_colours[i] <- "black"
      v_labels[i] <- paste(V(g)$name[i], v_number_of_papers[i], sep="\n", collapse = "") 
      #print(sprintf("The number of papers for node %s is: %s", v_labels[i], v_number_of_papers[i]))
    }
  }
  output <- list()
  output$v_shapes <- v_shapes
  output$v_number_of_papers <- v_number_of_papers
  output$v_colours <- v_colours
  output$v_label_colours <- v_label_colours
  output$v_labels <- v_labels
  return(output)
  
}

set_nodes_colous <- function(g, gname, gvalue){
  v_colours <- rep(0,length(V(g)$name))
  repos <- load_repos()
  gname <- repos$gname
  gvalue <- repos$gvalue
  #change nodes colour
  for(i in 1:length(V(g))){
    if(is.element(V(g)$name[i],gname)==TRUE){
      m <- match(V(g)$name[i],gname)
      if(as.numeric(gvalue[m])<0.05){
        v_colours[i] = '#ffd9b3' # is in the list with pvalue < 0.05 - orange
      }
      else{
        v_colours[i] ='#d7fcb5' # in the list with pvalue >= 0.05 - green
      }
    }
    if(is.element(V(g)[i]$name,gname)==FALSE){
      v_colours[i] = '#e4b5fc' # not in the list - purple
    }
  }
  return(v_colours)
}


create_edges_colours <- function(df){
  pathways <- unique(df$edge_added)
  colors_random <- distinctColorPalette(length(pathways))
  edge_colours <- NULL
  #rgb(runif(length(paths),runif(length(paths)),runif(length(paths))))
  #table <- cbind(unlist(KEGG_rep[1]$pathway.name), colors_random) 
  #graph_colours <- colors_random
  #pie(rep(1, length(pathways)), col=colors_random)
  #print(colors_random)
  #graph <- set.edge.attribute(graph, "label", value="df[3]")
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
  edge_labels <- NULL
  edge_path <- NULL
  for(i in 1: length(df$edge_added)){
    edge_labels[i] = list_of_colours[match(df$edge_added[i], pathways)]
    edge_path[i] = df$edge_added[i]
    #print(sprintf("edge_label %s for path %s ", edge_labels, edge_path))
  }
  output <- list()
  output$edge_labels <- edge_labels
  output$edge_path <-  edge_path
}

create_legend <- function(df, edge_info){
  legend_number <- length(unique(df$edge_added))
  #print(length(unique(df$edge_added)))
  #print(sprintf("length(legend$legend_number) is %i", length(legend$legend_number)))
  legend_name <- unique(df$edge_added)
  legend <- list()
  legend$legend_number <- legend_number
  legend$legend_name <- legend_name
  legend$legend_colour <- edge_info$list_of_colours
  return(legend)
}

plot_legend <- function(legend){
  plot(1:10, xaxt='n', ann=FALSE, yaxt='n', axes=FALSE, col = "white",ylim = c(1,20), xlim = c(1,10))
  pos = 15
  for(i in 1:legend$legend_number){
    if (i<10){
      text(0,pos, i, col = legend$legend_colour[i], pos = 4, cex = 0.8)
      text(0.3,pos, legend$legend_name[i], col = legend$legend_colour[i], pos = 4, cex = 0.8)
      pos = pos-0.8
    }
    if(i>9 && i<100){
      text(0,pos, i, col = legend$legend_colour[i], pos = 4, cex = 0.8)
      text(0.4,pos, legend$legend_name[i], col = legend$legend_colour[i], pos = 4, cex = 0.8)
      pos = pos-0.8
    }
    if(i>99){
      text(0,pos, i, col = legend$legend_colour[i], pos = 4, cex = 0.8)
      text(0.5,pos, legend$legend_name[i], col = legend$legend_colour[i], pos = 4, cex = 0.8)
      pos = pos-0.8
    }
  }
  par(xpd=TRUE)
}

plot_network <- function(graph){
  g <- graph$network
  v_shapes_labels <- graph$features
  l <- layout_in_circle(g)
  # pdf("rplot.pdf", width=6, height = 8) 
  plot(g, layout=l, vertex.shape = v_shapes_labels$v_shapes, vertex.label = v_shapes_labels$v_labels, vertex.size = 18, vertex.label.degree=-pi/2, 
       vertex.label.dist=0, vertex.color = v_shapes_labels$v_colours, vertex.frame.color = rgb(0,0,0,0), vertex.label.color="black",
       edge.curved = FALSE, edge.color = graph$edge_colours, edge.label = graph$edge_labels, edge.label.color=rgb(0,0,0,0))
  plot(g, layout=l,  vertex.frame.color=rgb(0,0,0,0), vertex.color=rgb(0,0,0,0), add=TRUE, #vertex.label = v_shapes_labels$v_number_of_papers,
       vertex.label.color = rgb(0,0,0,0),#vertex.label.color=  v_shapes_labels$v_label_colours, vertex.label.dist=-1.5, vertex.label.degree=-pi/2,
       edge.color = graph$edge_colours, edge.label = graph$edge_labels, edge.label.color = graph$edge_colours)
  #par(new=TRUE)
  # plot_legend(graph$legend)
  # dev.off()
}


create_network_pdf <- function(graph, output = "destination.pdf"){
  pdf(file = output, width=12, height = 6)
  par(mfrow=c(1,2))
  plot_network(network)
  plot_legend(network$legend)
  dev.off()
}

from_data_to_pathways <- function(GenesList, ListHeader){
  genes <- load_data(GenesList, ListHeader)
  database_paper_genes <- load_related_papers(genes)
  repos <- load_repos()
  KEGG_rep <- repos$KEGG
  genes_Pathways_matrix <- create_genes_Pathways_matrix(genes, KEGG_rep)
}
