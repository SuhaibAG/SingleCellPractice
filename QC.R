###Quality Control
##filtering genes with no counts
dim(sce) #27998 25806
detected_genes <- rowSums(counts(sce)) > 0
sce <- sce[detected_genes,]

dim(sce) #21140 25806
mean(detected_genes) ## 75% of genes are undetected

##annotating the genes
ah <- AnnotationHub()
ens.hs.107<- query(ah, c("Mus musculus", "EnsDb", 107))[[1]] 

genes <- rowData(sce)$ID
gene_annot <- AnnotationDbi::select(ens.hs.107, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME", "SYMBOL")) %>%
  set_names(c("ID", "Chromosome", "GeneName"))


head(gene_annot)

newRowData <- merge(rowData(sce), gene_annot, by = "ID", sort=FALSE, all = FALSE)
sce <- sce[newRowData$ID,]
rownames(rowData(sce)) <- rowData(sce)$ID

rowData(sce) <- merge(rowData(sce), gene_annot, by = "ID", sort=FALSE)


###filtering mito
is.mito <- which(rowData(sce)$Chromosome=="MT")
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)


###filtering
cell_qc_filters <- quickPerCellQC(colData(sce),
                                  sub.fields = TRUE,
                                  batch=sce$Sample)

as.data.frame(cell_qc_filters) %>% summarise(across(everything(), sum))
colData(sce) <- cbind(colData(sce), cell_qc_filters)




sce <- sce[, !sce$discard]

dim(sce) ##20800 24780



