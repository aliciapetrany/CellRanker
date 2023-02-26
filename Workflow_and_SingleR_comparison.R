#This file shows the CellRanker workflow, and compares its output to SingleR, the primary package for cell type annotation. 

#figure generation and comparison with SingleR
library("SingleR")
library("CellRankerPackage")
 library("Seurat")
#Loading in MTX files fromm Heidegger et. al

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


#analysis with CellRanker
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
