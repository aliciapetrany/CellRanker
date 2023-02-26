# CellRanker
## An R framework that performs rank-based cell type annotation of scRNA-Seq data

A Seurat-compatible package to perform cell type cluster annotations on scRNA-seq data

### R Shiny Application:
https://av33d2-alicia-petrany.shinyapps.io/cellranker/  
The application is compatible with 10xMtx formatted data from 10X genomics, or you can view interactive example data

### Package Install instructions
To install CellRanker, download the 'CellRankerPackage' folder, then open the project locally on your computer. Once in the folder, copy and paste the following into the R console: 
```
packages.install("devtools")
devtools::install()
```
Then, you can load in the package with the following command:
```
library(CellRankerPackage)
```
Finally, with any seurat object, CellRanker can be run with a single command A built in example dataset sourced from Heidegger et. al used in the following example: 
```
Seuratobj <- load_example_data()
Seuratobj <- CellRanker(Seuratobj, vis = TRUE)
```
It's as easy as that!  
The final output of CellRanker is an annotated Seurat object, as well as an optional tsne plot with the proper celltype labels. Celltype labels can then be applied in downstream analyses:  
![final_iteration](https://user-images.githubusercontent.com/28795694/221399640-7d94b6b6-d26e-41d5-b0bb-1320a4ce7c3f.svg)

#### References:
Heidegger, I., Fotakis, G., Offermann, A. et al. Comprehensive characterization of the prostate tumor microenvironment identifies CXCR4/CXCL12 crosstalk as a novel antiangiogenic therapeutic target in prostate cancer. Mol Cancer 21, 132 (2022). https://doi.org/10.1186/s12943-022-01597-7

