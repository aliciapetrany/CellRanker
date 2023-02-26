

#central function that performs workflow from start to finish
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

#Central CellRanker function for already processed seurat object: 
CellRanker_noprocess <- function(SeuratObj){
  marker.list <- find_markers(Seuratobj)
  pangao.list <- query_pangao(marker.list)
  celltype.list <- get_cell_type_possibilities(pangao.list)
  cluster.celltypes <- celltype_ranking(celltype.list, pangao.list)
  Seuratobj <- label_seurat(cluster.celltypes, Seuratobj)
  visualize(Seuratobj)
  return(Seuratobj) 
}

#with know markers
CellRanker_knownmarkers <- function(SeuratObj, marker.list){
  pangao.list <- query_pangao(marker.list)
  celltype.list <- get_cell_type_possibilities(pangao.list)
  cluster.celltypes <- celltype_ranking(celltype.list, pangao.list)
  Seuratobj <- label_seurat(cluster.celltypes, Seuratobj)
  visualize(Seuratobj)
  return(Seuratobj)
}

#Load Heidegger et al dataset
load_example_data <- function(){
  library(Seurat)
  samples <- list.files("example_data/")
  sample.list <- list()
  for (i in 1:length(samples)){
    dir <- paste("example_data/", samples[i], "/", sep = "")
    sample.list[[i]] <- Read10X(data.dir = dir)
    sample.list[[i]] <- CreateSeuratObject(sample.list[[i]])
  }
  
  Seuratobj <- merge(sample.list[[1]], y = c(sample.list[[2]],
                                             sample.list[[3]],
                                             sample.list[[4]],
                                             sample.list[[5]],
                                             sample.list[[6]],
                                             sample.list[[7]],
                                             sample.list[[8]]),
                                            add.cell.ids = samples)
  return(Seuratobj)
}

process_counts <- function(Seuratobj){
  library(Seurat)
  #Removes mitochondrial genes and low quality genes/cells
  Seuratobj[["percent.mt"]] <- PercentageFeatureSet(Seuratobj, pattern = "^MT-")
  Seuratobj <- subset(Seuratobj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
  
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

#finds the marker genes for each cluster
find_markers <- function(Seuratobj){
  
  clusters <- levels(Seuratobj@meta.data$seurat_clusters)
  marker.list <- list()
  for (i in 1:length(clusters)-1){
    print(i)
    marker.list[[(i+1)]] <- FindMarkers(Seuratobj, ident.1 = i, min.pct = 0.25)
  }
  
  return(marker.list)
}

#Retreives Reference Dataset From Pangao DB, which contains information about 
#which marker genes are canonically associated with different cell types
get_ref <- function(){
  download.file("https://www.panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", destfile = "control.tsv")
  return(read.delim("control.tsv"))
}

#identifies potential cell type subsets from pangao for each cluster
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

##returns all possible cell types for each cluster
get_cell_type_possibilities <- function(pangao.list){
  
  celltype.list <- list()
  for(i in 1:length(pangao.list)){
    celltype.list[[i]] <- levels(factor(pangao.list[[i]]$cell.type))
  }
  
  return(celltype.list)
}

#facilitaing function for ranking cell types
celltype_ranking <- function(celltype.list, pangao.list){
  
  ranking.list <- get_initial_ranking(celltype.list, pangao.list)
  ranking.list <- ranking_on_sens_spec_ubiq(celltype.list, pangao.list, ranking.list)
  cluster.celltypes <- cluster_assignments(celltype.list, ranking.list)
  return(cluster.celltypes)
}

#initial ranking. Based on # of cell type occurances
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

##INVESTIGATE NA PRODUCTION HERE
ranking_on_sens_spec_ubiq <- function(celltype.list, pangao.list, ranking.list, weights = c(1,1,1), NA.penalty = 0.0001){
  for(i in 1:length(pangao.list)){
    for(j in 1:length(celltype.list)){
      interest.rows <- pangao.list[[i]][which(pangao.list[[i]]$cell.type %in% celltype.list[[i]][j]),]
      interest.rows.temp <- interest.rows[!is.na(interest.rows$sensitivity_human),]
      
      #ranking by sensitivity
      if(nrow(interest.rows.temp) == 0 ){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*NA.penalty*weights[1] #rationale: if info is not available, gene cannot be eliminated
      }else{
        ranking.list[[i]][j] <- ranking.list[[i]][j]*prod(interest.rows$sensitivity_human)*weights[1]
      }
      
      #ranking by specficity
      interest.rows.temp <-interest.rows[!is.na(interest.rows$specificity_human),]
      if(nrow(interest.rows.temp) == 0 ){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*NA.penalty*weights[2] #rationale: if info is not available, gene cannot be eliminated
      }else{
        ranking.list[[i]][j] <- ranking.list[[i]][j]*prod(interest.rows$specificity_human)*weights[2]
      }
      
      #ranking by ubiquitous index (probability gene is randomly expressed in cell)
      interest.rows.temp <-interest.rows[!is.na(interest.rows$ubiquitousness.index),]
      if(nrow(interest.rows.temp) == 0 ){
        ranking.list[[i]][j] <- ranking.list[[i]][j]*NA.penalty*weights[3] #rationale: if info is not available, gene cannot be eliminated
      }else{
        ranking.list[[i]][j] <- ranking.list[[i]][j]*prod(interest.rows$ubiquitousness.index)*weights[3]
      }
    }
  }
  return(ranking.list)
}

#assigns celltypes to clusters based on rankings
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

#adds celltype annotations to seurat object
label_seurat <- function(cluster.celltypes, Seuratobj){
  cluster.nos <- Seuratobj@meta.data$seurat_clusters
  celltype_vector <- rep("0", length(cluster.nos))
  for(i in 1:length(levels(cluster.nos))){
    celltype_vector[cluster.nos %in% levels(cluster.nos)[i]] = cluster.celltypes[i]
  }
  Seuratobj@meta.data$cell_type <- celltype_vector
  return(Seuratobj)
}

#plots tsne coordinates with approriate cell type labels
visualize <- function(Seuratobj){
  DimPlot(Seuratobj, reduction = "tsne", group.by = "cell_type")
}


##NOTE TABLE DOES ALPHABETICAL ORDER GO FIX OTHER USE CASE

