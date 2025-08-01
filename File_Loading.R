library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)

###Loading in data
sample_sheet <- read_csv("/home/participant/Course_Materials/week3_single_cell/project_single_cell/sample_info.csv")
sample_list<- sample_sheet$sra_run
files <- str_c("/home/participant/Course_Materials/week3_single_cell/project_single_cell/preprocessed/cellranger/"
               ,sample_list , "/outs/filtered_feature_bc_matrix")
names(files) <- sample_list

bp.params <- MulticoreParam(workers = 7)

sce <- read10xCounts(files, col.names=TRUE, BPPARAM = bp.params)

dim(sce)
colData(sce)
rowData(sce)
colnames(sample_sheet)[colnames(sample_sheet) == 'sra_run'] <- 'Sample'

sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), sample_sheet, by="Sample", sort=FALSE)

sce$Sample <- NULL
colnames(colData(sce) )[colnames(colData(sce) ) == 'sample'] <- 'Sample'
rownames(colData(sce)) <- sce$Barcode
