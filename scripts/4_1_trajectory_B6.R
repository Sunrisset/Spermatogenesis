# ============== Trajectory analysis. ==============

#getwd()   # "/home/data/SpermatogenicCells"

library(slingshot)
library(Seurat)

load("/home/data/SpermatogenicCells/subgenes_B6_rename1.Rdata")
class(subgenes_B6_rename)

# Clusters.
sce <- as.SingleCellExperiment(subgenes_B6_rename, assay = "RNA")

sce_slingshot1 <- slingshot(sce,
                            reducedDim = 'UMAP',
                            clusterLabels = colData(sce)$ident,
                            approx_points = 150)

SlingshotDataSet(sce_slingshot1)
# 3861          2

levels(Idents(subgenes_B6_rename))
my_cols <-c("Spermatocytes1"="#E86FAA", "Elongating spermatids0"="#737884", "Elongating spermatids1"="#22232B", "Round spermatids1"="#6167AF", 
            "Spermatocytes0"="#F2B2D1", "Elongating spermatids2"="#737870", "Round spermatids2"="#B3A9D3", "Round spermatids3"="#D8DAEC", 
            "Elongating spermatids3"="#22232A", "Round spermatids4"="#2A348B","Elongating spermatids4"="#606B89","Elongating spermatids5"="#444A5E",
            "Somatic cells0"="#985D25","Round spermatids5"="#9595C9","Spermatogonia"="#046735","Elongating spermatids6"="#778866", 
            "Somatic cells1"="#643F18")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)


plot(reducedDims(sce_slingshot1)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot1), type = 'curves', lwd=2, col='yellow')

plot(reducedDims(sce_slingshot1)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot1), type = 'lineages', lwd=2, col='yellow')


sce_slingshot2 <- slingshot(sce,
                            reducedDim = 'UMAP',
                            clusterLabels = colData(sce)$ident,
                            start.clus = 'Spermatogonia',
                            approx_points = 150)

SlingshotDataSet(sce_slingshot2)
# 3861          2

plot(reducedDims(sce_slingshot2)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot2), type = 'curves', lwd=2, col='yellow')

plot(reducedDims(sce_slingshot2)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot2), type = 'lineages', lwd=2, col='yellow')

# With start.clus = 'Spermatogonia', the result of both the curves and the lineages are the same as that of start.clus = NULL.
# But with start.clus = 'Spermatogonia', the plotting was much more fast and saved time.