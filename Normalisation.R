library(scran)


set.seed(333)

##clustering
clust <- quickCluster(sce, BPPARAM=bp.params)
table(clust)

##using clusters to compute the size factors
sce <- computePooledFactors(sce,
                            clusters = clust,
                            min.mean = 0.1,
                            BPPARAM = bp.params)
deconv.sf <- sizeFactors(sce)
summary(deconv.sf)


lib.sf <- librarySizeFactors(sce)
sce$stage

###plottig size factors
data.frame(LibrarySizeFactors = lib.sf, 
           DeconvolutionSizeFactors = deconv.sf,
           stage = sce$stage) %>%
  ggplot(aes(x=LibrarySizeFactors, y=DeconvolutionSizeFactors)) +
  geom_point(aes(col=stage)) +
  geom_abline(slope = 1, intercept = 0)



##applying the size factors
sce <- logNormCounts(sce)
assayNames(sce)

