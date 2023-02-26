

#'Performs full CellRanker Workflow
#'
#'@param SeuratObj Seurat Object Containing Counts. If Counts are processed, use
#''CellRanker_noprocess' instead
#'@param vis Visualize Tsne Plot After Processing
#'
#'@export
CellRanker <- function(SeuratObj, vis = FALSE){
  Seuratobj <- process_counts(Seuratobj)
  marker.list <- find_markers(Seuratobj)
  pangao.list <- query_pangao(marker.list)
  celltype.list <- get_cell_type_possibilities(pangao.list)
  cluster.celltypes <- celltype_ranking(celltype.list, pangao.list)
  Seuratobj <- label_seurat(cluster.celltypes, Seuratobj)
  if(vis){
    visualize(Seuratobj )
  }

  return(Seuratobj)
}

#'Performs full CellRanker workflow without count preprocessing
#'
#'@param SeuratObj Seurat Object Containing Processed Counts.
#'@param vis Visualize Tsne Plot After Processing
#'
#'@export
CellRanker_noprocess <- function(SeuratObj, vis = FALSE){
  marker.list <- find_markers(Seuratobj)
  pangao.list <- query_pangao(marker.list)
  celltype.list <- get_cell_type_possibilities(pangao.list)
  cluster.celltypes <- celltype_ranking(celltype.list, pangao.list)
  Seuratobj <- label_seurat(cluster.celltypes, Seuratobj)
  visualize(Seuratobj)
  return(Seuratobj)
}

#'Performs full CellRanker Workflow if marker genes are already known
#'
#'@param SeuratObj Seurat Object Containing Processed Counts.
#'@param vis Visualize Tsne Plot After Processing
#'@param marker.list A list of marker genes. Structured as a list of dataframes
#'that is the same length of the number of clusters at SeuratObj.meta.data$cell_type.
#'Each list entry contaings a dataframe with a single column of adjust pvalues labeled
#'"p_val_adj'. The rownames of the dataframe should be the respective gene. Only
#'Only intented for use in debugging purposes.
#'@export
CellRanker_knownmarkers <- function(SeuratObj, marker.list, vis = FALSE){
  pangao.list <- query_pangao(marker.list)
  celltype.list <- get_cell_type_possibilities(pangao.list)
  cluster.celltypes <- celltype_ranking(celltype.list, pangao.list)
  Seuratobj <- label_seurat(cluster.celltypes, Seuratobj)
  visualize(Seuratobj)
  return(Seuratobj)
}


#'Preprocesses counts by follwing Seurat's suggested preprocessing workflow
#'
#'@param Seuratobj Seurat Object containing non-processed counts. Performs clustering
#'
#'@export
process_counts <- function(Seuratobj){
  library(Seurat)
  #Removes mitochondrial genes and low quality genes/cells
  Seuratobj[["percent.mt"]] <- PercentageFeatureSet(Seuratobj, pattern = "^MT-")
  Seuratobj <- subset(Seuratobj,
                      subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

  #normalizes/scales counts, finds variable features
  Seuratobj <- NormalizeData(Seuratobj)
  Seuratobj <- ScaleData(Seuratobj, features = rownames(Seuratobj))
  Seuratobj <- FindVariableFeatures(object = Seuratobj)

  #calculates TSNE coordinates
  Seuratobj <- RunPCA(Seuratobj)
  Seuratobj <- RunTSNE(Seuratobj)

  #Clustering to identify cell type clusters
  Seuratobj <- FindNeighbors(Seuratobj, dims = 1:10)
  Seuratobj <- FindClusters(Seuratobj, resolution = 0.5)

  #returning processed Seuratobj
  return(Seuratobj)
}

#'Finds the marker genes for each cluster
#'
#'@param Seuratobj Seurat Object containing processed counts and cluster annotations
#'
#'@export
find_markers <- function(Seuratobj){

  clusters <- levels(Seuratobj@meta.data$seurat_clusters)
  marker.list <- list()
  for (i in 1:length(clusters)-1){
    print(i)
    marker.list[[(i+1)]] <- FindMarkers(Seuratobj, ident.1 = i, min.pct = 0.25)
  }

  return(marker.list)
}

#'Retreives dataset from Pangao database, which contains information about which
#'marker genes correspond to a given cell type
#'
#'@export
get_ref <- function(){
  download.file("https://www.panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", destfile = "control.tsv")
  return(read.delim("control.tsv"))
}

#'Queries the pangao dataset for marker genes with signficant relations to
#'different cell types
#'
#'@param marker.list a list of dataframes containing marker genes as row names and
#'a column labeled p_val_adj that contains adusted pvalues for each gene. THe length
#'of the list corresponseds to the number of clusters
#'@export
query_pangao <- function(marker.list){
  pangao.list <- list()
  pangao<- get_ref()

  for(i in 1:length(marker.list)){
    #Gets signficiant marker genes per cluster
    markers <- marker.list[[i]]
    markers <- markers[markers$p_val_adj < 0.5,]

    #Limits analysis to top 100 markers
    if(nrow(markers) > 100){
      markers <- rownames(markers[1:100,])
    }else{
      markers <- rownames(markers)
    }

    #Returns pangao entries for each gene
    tempdf <- data.frame()
    for(j in 1:length(markers)){
      genewise.df <- pangao[grep(markers[j], pangao$official.gene.symbol),]
      tempdf <- rbind(tempdf, genewise.df)
    }
    pangao.list[[i]] <- tempdf
  }

  #final output is a list containing pangao dataframes with markers from each cluster
  return(pangao.list)
}

#'Returns all possible cell types for a given cluster based on marker genes
#'
#'@param pangao.list A list the same length as the number of clusters. Each list
#'entry ontains a sliced dataframe with information regarding each relevant marker
#'gene to the respective cluster.
#'
#'@export
get_cell_type_possibilities <- function(pangao.list){

  celltype.list <- list()
  for(i in 1:length(pangao.list)){
    celltype.list[[i]] <- levels(factor(pangao.list[[i]]$cell.type))
  }

  #returns a list of celltypes for each cluster
  return(celltype.list)
}

#'Facilitates the ranking of celltypes
#'
#'@param celltype.list A list the same length as the number of clusters
#'that contains a vector of all possible celltyeps for each respective
#'cluster
#'@param pangao.listA list the same length as the number of clusters. Each list
#'entry ontains a sliced dataframe with information regarding each relevant marker
#'gene to the respective cluster.
#'
#'@export
celltype_ranking <- function(celltype.list, pangao.list){

  ranking.list <- get_initial_ranking(celltype.list, pangao.list)
  ranking.list <- complex_ranking(celltype.list, pangao.list, ranking.list)
  cluster.celltypes <- cluster_assignments(celltype.list, ranking.list)
  return(cluster.celltypes)
}

#'Helper function for celltype_ranking that intially ranks cell types
#'by number of occurances
#'#'@param celltype.list A list the same length as the number of clusters
#'that contains a vector of all possible celltyeps for each respective
#'cluster
#'@param pangao.listA list the same length as the number of clusters. Each list
#'entry ontains a sliced dataframe with information regarding each relevant marker
#'gene to the respective cluster.
get_initial_ranking <- function(celltype.list, pangao.list){
  ranking.list <- list()
  for(i in 1:length(pangao.list)){
    ranking <- vector()
    for(j in 1:length(celltype.list[[i]])){
      ranking<- c(ranking, table(pangao.list[[i]]$cell.type)[[j]])
    }
    ranking.list[[i]] <- ranking
  }
  return(ranking.list)
}

#'Helper function that ranks by the following:
#'Number of genes for each celltype. Although it sounds counterintuitive, too many
#'genes promotes bias towards more general cell types
#'Sensitivity - probability that a gene is expressed in given celltype
#'Specificity - probaility that gene is not expressed in other celltypes
#'Ubiquitousness index - probability that a gene is expressed at random,
#'  i.e a housekeeping gene
#'Top marker gene rank
#'
#'#'Helper function for celltype_ranking that intially ranks cell types
#'by number of occurances
#'#'@param celltype.list A list the same length as the number of clusters
#'that contains a vector of all possible celltyeps for each respective
#'cluster
#'@param pangao.listA list the same length as the number of clusters. Each list
#'entry ontains a sliced dataframe with information regarding each relevant marker
#'gene to the respective cluster.
#'@param ranking.list A list vectors containing numeric rankings for each cluster.
#'@param weights Wieghts to change the importance of sensitivity, specificity, and ui
#'@param NA.penalty A penalty to incur on the ranking if the database does not contain
#'information on a given value. Necessary to prevent genes with little documentation
#'from being wieghted signficantly
complex_ranking <- function(celltype.list, pangao.list, ranking.list, weights = c(0.05,1,1.5,2,3), NA.penalty = 0.0001){
  for(i in 1:length(pangao.list)){
    for(j in 1:length(celltype.list)){
      #isolates rows of j'th celltype for ith cluster
      interest.rows <- pangao.list[[i]][which(pangao.list[[i]]$cell.type %in% celltype.list[[i]][j]),]
      interest.rows.temp <- interest.rows[!is.na(interest.rows$sensitivity_human),]

      #Adds initial penalty if a certain celltype has too many genes to try to encourage more specific cell types
      ngenes <- nrow(interest.rows)
      ranking.list[[i]][j] <- ranking.list[[i]][j]*(1- weights[1]*ngenes)

      #ranking by sensitivity
      if(nrow(interest.rows.temp) == 0 ){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*NA.penalty*weights[2] #rationale: if info is not available, gene cannot be eliminated
      }else{
        ranking.list[[i]][j] <- ranking.list[[i]][j]*prod(interest.rows$sensitivity_human)*weights[2]
      }

      #ranking by specficity
      interest.rows.temp <-interest.rows[!is.na(interest.rows$specificity_human),]
      if(nrow(interest.rows.temp) == 0 ){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*NA.penalty*weights[3] #rationale: if info is not available, gene cannot be eliminated
      }else{
        ranking.list[[i]][j] <- ranking.list[[i]][j]*prod(interest.rows$specificity_human)*weights[3]
      }

      #ranking by ubiquitous index (probability gene is randomly expressed in cell)
      interest.rows.temp <-interest.rows[!is.na(interest.rows$ubiquitousness.index),]
      if(nrow(interest.rows.temp) == 0 ){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*NA.penalty*weights[4] #rationale: if info is not available, gene cannot be eliminated
      }else{
        ranking.list[[i]][j] <- ranking.list[[i]][j]*prod(interest.rows$ubiquitousness.index)*weights[4]
      }

      #Multiplying by the rank of most signficant marker gene, no penalty if top 2 genes
      top.marker <- interest.rows$official.gene.symbol[1]
      rank <- which(top.marker == top.marker)
      if(!(rank <= 2)){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*(1/log(rank))*weights[5]
      }


    }
  }
  return(ranking.list)
}

#'Determines the highest ranked cell type for each cluster. For duplicate clusters,
#'rankings are repeated with the next highest cell type until all cell type classifications
#'are unique.
#'
#'#'#'@param celltype.list A list the same length as the number of clusters
#'that contains a vector of all possible celltyeps for each respective
#'cluster
#'@param ranking.list A list vectors containing numeric rankings for each cluster.
cluster_assignments <- function(celltype.list, ranking.list){
  library(stringi)

  #order celltype assignments based on rankings
  for(i in 1:length(celltype.list)){
    celltype.list[[i]] <-celltype.list[[i]][order(ranking.list[[i]], decreasing = TRUE)]
    ranking.list[[i]] <- ranking.list[[i]][order(ranking.list[[i]], decreasing = TRUE)]
  }

  #determine top label for each cluster
  next.layer.indices = c(1:length(celltype.list))
  layer = 1
  complete = FALSE
  celltypes.of.layer = vector()
  while(!complete){

      #determines most likely cell type based on layer, with layer being
      #number of positions queried on each rank list
      for(i in next.layer.indices){
          celltypes.of.layer[i] = celltype.list[[i]][layer]
      }


    #Deals with duplicate labels. For duplicates, the score with the highest rank is accepted
    if(length(levels(factor(celltypes.of.layer))) != length(celltype.list)){
        next.layer.indices <- vector()
        repeats <- which(stri_duplicated(celltypes.of.layer))
        for(j in repeats){
            fight.to.death <- which(celltypes.of.layer %in% celltypes.of.layer[j])
            if(ranking.list[[fight.to.death[1]]][layer] >= ranking.list[[fight.to.death[2]]][layer]){
                next.layer.indices <- c(next.layer.indices, fight.to.death[2])
            }else{
                next.layer.indices <- c(next.layer.indices, fight.to.death[1])
            }
        }

      }else{
        complete = TRUE
      }
      layer = layer +1
  }

  return(celltypes.of.layer)
}

#'Adds final cell type annotations to Seurat Object
#'
#'@param cluster.celltypes A vector containing cell type annotations for
#'each cluster.
#'@param Seuratobj A Seurat object with processed counts and tsne coordinates
#'
#'@export
label_seurat <- function(cluster.celltypes, Seuratobj){
  cluster.nos <- Seuratobj@meta.data$seurat_clusters
  celltype_vector <- rep("0", length(cluster.nos))
  for(i in 1:length(levels(cluster.nos))){
    celltype_vector[cluster.nos %in% levels(cluster.nos)[i]] = cluster.celltypes[i]
  }
  Seuratobj@meta.data$cell_type <- celltype_vector
  return(Seuratobj)
}


#'Plots TSNE coordinates with appropriate CellRanker assigned labels
#'
#'@param Seuratobj a Seurat object with processed counts, tsne coordinates, and
#'CellRanker assigned labels
#'@export
visualize <- function(Seuratobj){
  DimPlot(Seuratobj, reduction = "tsne", group.by = "cell_type")
}


