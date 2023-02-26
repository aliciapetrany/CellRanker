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
