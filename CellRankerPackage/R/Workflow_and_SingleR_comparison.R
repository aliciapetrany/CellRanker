#figure generation and comparison with SingleR
library("SingleR")
library("CellRankerPackage")

#analysis with CellRanker
Seuratobj <- load_example_data()
Seuratobj.CellRanker <- CellRanker()
#...and that's it

#Analysis with singleR
#processing counts with CellRanker preprocessing function for brevity
Seuratobj.SingleR <- process_counts(Seuratobj)

#converted seurat object to single cell experiment class, then
#running SingleR with Human Primary Cell Atlas Data as a reference
SingleR.sce <- as.SingleCellExperiment(Seuratobj.processed)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
SingleR.pred <- SingleR(test = SingleR.sce, ref = hpca.se,
                        labels = hpca.se$label.main)

#assigning SingleR predicions to Seurat meta data
Seuratobj.SingleR$celltype <- SingleR.pred$labels

#Visualizing final results
DimPlot(Seuratobj.SingleR, reduction = "tsne", group.by = 'celltype')
DimPlot(Seuratobj.CellRanker, reduction = "tsne", group.by = 'cell_type')





#not for final vignette
Seuratobj <- load_example_data()
Seuratobj <- process_counts(Seuratobj)
marker.gene.files <-list.files("marker_genes")
marker.gene.files <- marker.gene.files[c(1,10,11,12,13,14,15,16,17,2,3,4,5,6,7,8,9)]
marker.list <- list()
for (i in 1:length(marker.gene.files)){
  marker.list[[i]] <- read.csv(paste("marker_genes/", marker.gene.files[i], sep = ""))
  rownames(marker.list[[i]]) <- marker.list[[i]][,1]
  marker.list[[i]] <- marker.list[[i]][-1]
}
Seuratobj.final <- CellRanker_knownmarkers(Seuratobj, marker.list)
