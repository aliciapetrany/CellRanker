Seuratobj <- load_example_data()

#count preprocessing

starttime <- Sys.time()
Seuratobj <- process_counts(Seuratobj)

marker.list <- find_markers(Seuratobj)
endtime <- Sys.time()
processing.runtime <- endtime - starttime

#saving marker genes for later use:
filenames <- paste("marker_genes/" , 1:17, ".csv", sep = "")
for(i in 1:length(marker.list)){
  write.csv(marker.list[[i]], file = filenames[i])
}

marker.list.files <- list.files(path = "marker_genes/")
marker.list.files <- marker.list.files[c(1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,8,9)]
marker.list <- list()
for (i in 1:length(marker.list.files)){
  marker.list[[i]] <- read.csv(paste("marker_genes/", marker.list.files[i], sep = ""))
  rownames(marker.list[[i]]) <- marker.list[[i]]$X
  marker.list[[i]] <- marker.list[[i]][,-1]
  }

pangao.list <- query_pangao(marker.list)
celltype.list <- get_cell_type_possibilities(pangao.list)
cluster.celltypes <- celltype_ranking(celltype.list, pangao.list)
Seuratobj <- label_seurat(cluster.celltypes, Seuratobj)
visualize(Seuratobj)

tsne_coords <- data.frame(Seuratobj@reductions$tsne@cell.embeddings)
tsne_coords$labels <- Seuratobj@meta.data$cell_type

ReadMtx(mtx = )