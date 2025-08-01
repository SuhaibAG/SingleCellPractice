library(PCAtools)
library(patchwork)

### making the gene name the identifier instead of the id 
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)
#checking variation
gene_var <- modelGeneVar(sce)



###feature selection
gene_var %>% 
  as.data.frame() %>% 
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")


hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs)

hvgs[1:10]

plotExpression(sce, features = hvgs[1:20], point_alpha = 0.05)


###dimensionality reduction
##PCA with 50 PCAS
sce <- runPCA(sce, subset_row = hvgs)
sce
dim(reducedDim(sce, "PCA"))
reducedDim(sce, "PCA")[1:10, 1:1]


percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
sum(percent.var[1:6])
##only 6 PCA's are important after the elbow test


##plotting PCAs
plotReducedDim(sce, dimred = "PCA", colour_by = "stage")
plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "stage")


ggcells(sce, aes(x = PCA.1, y = PCA.2, colour = stage)) +
  geom_point(size = 0.5) +
  facet_wrap(~ stage) +
  labs(x = "PC1", y = "PC2", colour = "Sample")


explain_pcs <- getExplanatoryPCs(sce,
                                 variables = c("sum",
                                               "detected",
                                               "stage",
                                               "Sample",
                                               "subsets_Mito_percent")
)

plotExplanatoryPCs(explain_pcs/100)


plotExplanatoryVariables(sce,
                         variables = c("sum",
                                       "detected",
                                       "stage",
                                       "Sample",
                                       "subsets_Mito_percent"))

table(percent.var > 1)
chosen_elbow <- findElbowPoint(percent.var)
chosen_elbow
plot(percent.var)
abline(v=chosen_elbow, col="dodgerblue")


sce.denoised <- denoisePCA(sce, technical = gene_var, subset.row = hvgs)

ncol(reducedDim(sce.denoised, "PCA"))



sce$stage
##tsne for visualization
sce <- runTSNE(sce)
plotTSNE(sce, colour_by = "stage")


sce <- runTSNE(sce, 
               name = "TSNE_perplex50",
               perplexity = 50, 
               dimred = "PCA",
               n_dimred = 10)

ggcells(sce, aes(x = TSNE_perplex50.1, 
                 y = TSNE_perplex50.2, 
                 colour = stage)) +
  geom_point()

##UMAP
sce$Sample
sce <- runUMAP(sce)
plotUMAP(sce, colour_by = "Sample")

sce <- runUMAP(sce, 
               name = "UMAP_neighbors500",
               dimred = "PCA",
               n_neighbors = 500)

ggcells(sce, aes(x = UMAP_neighbors500.1, y = UMAP_neighbors500.2, 
                 colour = Sample)) +
  geom_point() +
  labs(title = "Neighbours = 500")


