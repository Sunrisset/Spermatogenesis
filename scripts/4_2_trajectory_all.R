# ============== Trajectory analysis. ==============

#getwd()   # "/home/data/SpermatogenicCells"

library(slingshot)
library(Seurat)

load("/home/data/SpermatogenicCells/subgenes_rename.Rdata")
class(subgenes_rename)
str(subgenes_rename)
# "SeuratObject".
# $ counts:27918 21774.


# Clusters.
sce <- as.SingleCellExperiment(subgenes_rename, assay = "RNA")
str(sce) # Formal class 'SingleCellExperiment' [package "SingleCellExperiment"] with 9 slots.
colData(sce) # DataFrame with 21774 rows and 8 columns.
class(colData(sce)$seurat_clusters)  # "factor".
class(colData(sce)$ident)  # "factor".

sce_slingshot1 <- slingshot(sce,
                            reducedDim = 'UMAP',
                            clusterLabels = colData(sce)$ident,
                            approx_points = 150)

SlingshotDataSet(sce_slingshot1)
# Samples Dimensions
# 21774          2

levels(Idents(subgenes_rename))
my_cols <- c("Spermatocytes2"="#E86FAA", "Elongating spermatids1"="#737884", "Spermatocytes3"="deeppink2", "Round spermatids0"="#B3A9D3", "Somatic cells3"="goldenrod3", "Spermatocytes1"="#F2B2D1",
             "Spermatocytes5"="#9A1C59", "Elongating spermatids0"="#22232B", "Spermatogonia"="#046735", "Round spermatids2"="#282B69", "Round spermatids1"="#6167AF",
             "Somatic cells0"="#643F18","Spermatocytes8"="mediumorchid4","Spermatocytes4"="#D31281","Spermatocytes10"="#502269", "Somatic cells4"="goldenrod2",
             "Spermatocytes7"="plum2", "Spermatocytes0"="#FBDCEA", "Somatic cells5"="goldenrod1", "Somatic cells6"="yellow1", "Somatic cells7"="khaki1",
             "Spermatocytes9"="#7E2678", "Spermatocytes6"="mediumorchid1", "Somatic cells2"="goldenrod4", "Somatic cells1"="#985D25", "Somatic cells8"="lemonchiffon")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)

plot(reducedDims(sce_slingshot1)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot1), lwd=2, col='greenyellow')

plot(reducedDims(sce_slingshot1)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot1), type = 'lineages', lwd=2, col='greenyellow')



sce_slingshot2 <- slingshot(sce,
                            reducedDim = 'UMAP',
                            clusterLabels = colData(sce)$ident,
                            start.clus = 'Spermatogonia',
                            approx_points = 150)

SlingshotDataSet(sce_slingshot2)
# Samples Dimensions
# 21774          2

plot(reducedDims(sce_slingshot2)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot2), lwd=2, col='greenyellow')

plot(reducedDims(sce_slingshot2)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot2), type = 'curves', lwd=2, col='greenyellow')

plot(reducedDims(sce_slingshot2)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot2), type = 'lineages', lwd=2, col='greenyellow')

plot(reducedDims(sce_slingshot2)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot2), type = 'both', lwd=2, col='greenyellow')



sce_slingshot3 <- slingshot(sce,
                            reducedDim = 'UMAP',
                            clusterLabels = colData(sce)$ident,
                            start.clus = 'Somatic cells0',
                            approx_points = 150)

SlingshotDataSet(sce_slingshot3)

plot(reducedDims(sce_slingshot3)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot3), type = 'curves', lwd=2, col='greenyellow')

plot(reducedDims(sce_slingshot3)$UMAP, col = my_cols2[colData(sce)$ident], pch=16, asp = 1)
lines(SlingshotDataSet(sce_slingshot3), type = 'lineages', lwd=2, col='greenyellow')
