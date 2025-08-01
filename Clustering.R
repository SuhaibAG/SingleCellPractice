library(scater)
library(scran)
library(bluster)
library(cluster)
library(igraph)
library(pheatmap)
library(patchwork)
library(tidyverse)


sce
clustering1 <- clusterCells(sce, use.dimred="UMAP", full=TRUE)
table(clustering1$clusters)




###lpuvain clustering
sce$louvain15 <- clusterCells(sce, 
                              use.dimred = "UMAP", 
                              BLUSPARAM = SNNGraphParam(k = 20, 
                              cluster.fun = "louvain"))

plotReducedDim(sce, 
               dimred = "TSNE",
               colour_by = "louvain15", 
               text_by = "louvain15")


sce
sil.approx <- approxSilhouette(reducedDim(sce, "UMAP"),
                               clusters=sce$louvain15)
sil.approx


colLabels(sce) <- sce$louvain15
sce$
plotReducedDim(sce, 
               dimred = "TSNE",
               colour_by = "stage", 
               text_by = "stage") +
  ggtitle("louvain15")



rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

top_10 <- hvgs[1:10]

plotReducedDim(
    sce,
    dimred = "TSNE",
    by_exprs_values = "logcounts",
    colour_by = top_10,
    text_by = "label"
) 


plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c(top_10))



plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c(sce$stage))

table(sce$stage)
top_10
cald1<- plotReducedDim(sce, 
               dimred = "TSNE",
               by_exprs_values = "logcounts",
               colour_by = "Cald1",
               text_by = "stage")
Tagln<- plotReducedDim(sce, 
                       dimred = "TSNE",
                       by_exprs_values = "logcounts",
                       colour_by = "Tagln",
                       text_by = "stage")
plotReducedDim(sce, 
                       dimred = "TSNE",
                       by_exprs_values = "logcounts",
                       colour_by = "Tagln",
                       text_by = "stage")
sce$detected

rowData(sce)
getTopHVGs(gene_var[], prop=0.1)

summstagesummary(hvgs)




###plotting for my analysis
plotter <- function(gene){
  plotReducedDim(sce, 
                 dimred = "TSNE",
                 by_exprs_values = "logcounts",
                 colour_by = gene,
                 text_by = "stage")
}
top_50 <- getTopHVGs(gene_var, prop=0.1)[1:50]
plots <- lapply(top_50[20:40], plotter)

combined_plot <- wrap_plots(plots)
print(combined_plot)

colData(sce)
rowData(sce)
getTopHVGs(gene_var, prop=0.1)


####

