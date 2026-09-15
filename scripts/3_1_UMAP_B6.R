R.Version()
# R version 4.3.3 (2024-02-29)
# Creates a SingleCellExperiment.

# ============== UMAP. ==============

getwd()   # "/home/data/SpermatogenicCells"
ibrary(Seurat)
library(dplyr)

seura_do17815_B6 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17815/outs/filtered_feature_bc_matrix"), project = "do17815")
str(seura_do17815_B6)
# 32285 4036.

seura_do17816_B6 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17816/outs/filtered_feature_bc_matrix"), project = "do17816")
str(seura_do17816_B6)
# 32285 5131.

merge_B6 <- merge(x = seura_do17815_B6, y = seura_do17816_B6, add.cell.ids = c("do17815", "do17816"), project = "sperm_B6")

unique(sapply(X = strsplit(colnames(merge_B6), split = "_"), FUN = "[", 1))
table(merge_B6$orig.ident)
save(merge_B6, file = "seura_merge_B6_sperm.Rdata")

merge_B6_new = JoinLayers(object = merge_B6)

rs_B6 <- rowSums(GetAssayData(object = merge_B6_new, assay="RNA",layer ="counts")>0)

fe_B6 <- rs_B6[which(rs_B6>3)]

subgenes_B6 <- subset(merge_B6,features=names(fe_B6))

subgenes_B6[["percent.mt"]] <- PercentageFeatureSet(subgenes_B6, pattern = "^MT-")
dim(subgenes_B6[["percent.mt"]])

# calculate blood cells.
HB.genes_B6 <-c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m_B6 <- match(HB.genes_B6, rownames(subgenes_B6@assays$RNA)) 
HB.genes_B6 <- rownames(subgenes_B6@assays$RNA)[HB_m_B6] 
HB.genes_B6 <- HB.genes_B6[!is.na(HB.genes_B6)] 
subgenes_B6[["percent.HB"]]<-PercentageFeatureSet(subgenes_B6, features=HB.genes_B6) 

VlnPlot(subgenes_B6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4)
Plot1AB_B6 <- FeatureScatter(subgenes_B6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");Plot1AB_B6

subgenes_B6_filtered<- subset(subgenes_B6, subset = nFeature_RNA > 2000 & nCount_RNA > 10000 & percent.mt < 5)

VlnPlot(subgenes_B6_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Plot1AB_filtered <- FeatureScatter(subgenes_B6_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");Plot1AB_filtered


subgenes_B6_filtered_norm <- NormalizeData(subgenes_B6_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
subgenes_B6_varFeatures <- FindVariableFeatures(subgenes_B6_filtered_norm, selection.method = "vst", nfeatures = 2000);subgenes_B6_varFeatures

top10_B6 <- head(VariableFeatures(subgenes_B6_varFeatures), 10);top10_B6
plot03A_B6 <- VariableFeaturePlot(subgenes_B6_varFeatures);plot03A_B6
plot03B_B6 <- LabelPoints(plot = plot03A_B6, points = top10_B6, repel = TRUE, xnudge = 0, ynudge = 0);plot03B_B6

# ---Scale the data.---
all.genes_B6 <- rownames(subgenes_B6_varFeatures);length(all.genes_B6)
subgenes_B6_varFeatures_Scacle <- ScaleData(subgenes_B6_varFeatures, features = all.genes_B6, 
                                            do.scale = TRUE, do.center = TRUE,
                                            model.use = "linear", scale.max = 10)

# ---Perform linear dimensional reduction.---
subgenes_B6_vs_PCA <- RunPCA(subgenes_B6_varFeatures_Scacle, features = VariableFeatures(object = subgenes_B6_varFeatures_Scacle), npcs = 60, rev.pca = FALSE)

# examine and visualize PCA results in different ways.
print(subgenes_B6_vs_PCA[["pca"]], dims = 1:12, nfeatures = 10)
VizDimLoadings(subgenes_B6_vs_PCA, dims = 1:6, reduction = "pca", col = "purple", nfeatures = 10)
DimPlot(subgenes_B6_vs_PCA, reduction = "pca")
DimHeatmap(subgenes_B6_vs_PCA, dims = 1, cells = 1000, nfeatures = 20, balanced = TRUE, slot = "scale.data")
DimHeatmap(subgenes_B6_vs_PCA, dims = 1:15, cells = 1000, balanced = TRUE)
DimHeatmap(subgenes_B6_vs_PCA, dims = 1:6, cells = 1000, balanced = TRUE)

# ---Determine the 'dimensionality' of the dataset.---
subgenes_B6_PCs <- JackStraw(subgenes_B6_vs_PCA, dims = 30, num.replicate = 200, maxit = 2000)
subgenes_B6_PCs_score <- ScoreJackStraw(subgenes_B6_PCs, dims = 1:30)
JackStrawPlot(subgenes_B6_PCs_score, dims = 1:30)
ElbowPlot(subgenes_B6_PCs_score)


# ---Cluster the cells.---
subgenes_B6_SNN <- FindNeighbors(subgenes_B6_PCs_score, dims = 1:30, k.param = 20)
subgenes_B6_SNN_Cluster <- FindClusters(subgenes_B6_SNN, resolution = 0.5)

# ---Run non-linear dimensional reduction (UMAP/tSNE).---
subgenes_B6_UMAP <- RunUMAP(subgenes_B6_SNN_Cluster, reduction = "pca", dims = 1:60,  n.neighbors = 40L)
DimPlot(subgenes_B6_UMAP, reduction = "umap")
DimPlot(subgenes_B6_UMAP, reduction = "umap", label = TRUE, pt.size = 1.5)
save(subgenes_B6_UMAP, file = "subgenes_B6_UMAP.Rdata")

library(scCustomize)
# ref original paper 22 clusters. If here are 22 clusters, then should be something like:
DimPlot_scCustom(subgenes_B6_UMAP, reduction = "umap", label = TRUE, pt.size = 0.5,
                 colors_use = c("#643F18", "#985D25", "#046735",
                                "#FBDCEA", "#F2B2D1", "#E86FAA", "#D31281", "#9A1C59", "#7E2678", "#502269", "#240F2F",
                                "#B3A9D3", "#D8DAEC", "#9595C9", "#6167AF", "#2A348B", "#282B69", "#171447",
                                "#22232B", "#444A5E", "#606B89", "#737884"))

# ---Find differentially expressed features (cluster biomarkers).---
subgenes_B6_UMAP_new = JoinLayers(object = subgenes_B6_UMAP)
save(subgenes_B6_UMAP_new, file = "subgenes_B6_UMAP_new.Rdata")

# find markers for every cluster compared to all remaining cells, report only the positive ones.
subgenes_B6.markers <- FindAllMarkers(subgenes_B6_UMAP_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log2(2));class(subgenes_B6.markers);dim(subgenes_B6.markers)

# to see all the markers with top 10 avg_logFC in all clusters.
top10_B6 <- subgenes_B6.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# 10*17=170.
DoHeatmap(subgenes_B6_UMAP_new, features = top10_B6$gene[1:20])

# feature/count visualization.
# Somatic cells.
VlnPlot(subgenes_B6_UMAP_new, features = c("Cldn11", "Fabp3"))
VlnPlot(subgenes_B6_UMAP_new, features = c("Cldn11", "Fabp3"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_B6_UMAP_new, features = c("Cldn11", "Fabp3"))
FeaturePlot(subgenes_UMAP_new, features = c("Cldn11"), pt.size = 1)
FeaturePlot(subgenes_UMAP_new, features = c("Fabp3"), pt.size = 1)

# Spermatogonia.
VlnPlot(subgenes_B6_UMAP_new, features = c("Dmrt1"))
VlnPlot(subgenes_B6_UMAP_new, features = c("Dmrt1"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_B6_UMAP_new, features = c("Dmrt1"))

# Spermatocytes.
VlnPlot(subgenes_B6_UMAP_new, features = c("Piwil1"))
VlnPlot(subgenes_B6_UMAP_new, features = c("Piwil1"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_B6_UMAP_new, features = c("Piwil1"))

# Round spermatids.
VlnPlot(subgenes_B6_UMAP_new, features = c("Tex21"))
VlnPlot(subgenes_B6_UMAP_new, features = c("Tex21"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_B6_UMAP_new, features = c("Tex21"))

# Elongating spermatids.
VlnPlot(subgenes_B6_UMAP_new, features = c("Tnp1"))
VlnPlot(subgenes_B6_UMAP_new, features = c("Tnp1"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_B6_UMAP_new, features = c("Tnp1"))


# ---Assign cell type identity to clusters.---
new.cluster.ids_B6 <- c("Spermatocytes", "Elongating spermatids", "Elongating spermatids", "Round spermatids", "Spermatocytes", "Elongating spermatids", "Round spermatids", "Round spermatids", "Elongating spermatids", "Round spermatids","Elongating spermatids","Elongating spermatids","Somatic cells","Round spermatids","Spermatogonia","Elongating spermatids", "Somatic cells")

names(new.cluster.ids_B6) <- levels(subgenes_B6_UMAP_new)
subgenes_B6_rename <- RenameIdents(subgenes_B6_UMAP_new, new.cluster.ids_B6)
DimPlot(subgenes_B6_rename, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
save(subgenes_B6_rename, file = "subgenes_B6_rename.Rdata")

# ---Visualization.---
library(scCustomize)

DimPlot_scCustom(subgenes_B6_UMAP_new, reduction = "umap", label = FALSE, pt.size = 1.5,
                 colors_use = c("#E86FAA", "#737884", "#22232B", "#6167AF", "#F2B2D1", "#737870", "#B3A9D3", "#D8DAEC", "#22232A", "#2A348B","#606B89","#444A5E","#985D25","#9595C9","#046735","#778866", "#643F18"))

library(ggplot2)
DimPlot_scCustom(subgenes_B6_UMAP_new, reduction = "umap", label = FALSE, pt.size = 0.5,
                 colors_use = c("#E86FAA", "#737884", "#22232B", "#6167AF", "#F2B2D1", "#737870", "#B3A9D3", "#D8DAEC", "#22232A", "#2A348B","#606B89","#444A5E","#985D25","#9595C9","#046735","#778866", "#643F18")) +
  theme(legend.position = "none")

# ---Assign cell type identity to clusters 2.---
new.cluster.ids_B6 <- c("Spermatocytes1", "Elongating spermatids0", "Elongating spermatids1", "Round spermatids1", "Spermatocytes0", "Elongating spermatids2", "Round spermatids2", "Round spermatids3", "Elongating spermatids3", "Round spermatids4","Elongating spermatids4","Elongating spermatids5","Somatic cells0","Round spermatids5","Spermatogonia","Elongating spermatids6", "Somatic cells1")

names(new.cluster.ids_B6) <- levels(subgenes_B6_UMAP_new)
subgenes_B6_rename <- RenameIdents(subgenes_B6_UMAP_new, new.cluster.ids_B6)
DimPlot(subgenes_B6_rename, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
save(subgenes_B6_rename, file = "subgenes_B6_rename1.Rdata")

# ---Visualization 2.---
library(scCustomize)

DimPlot_scCustom(subgenes_B6_rename, reduction = "umap", label = FALSE, pt.size = 1.5,
                 colors_use = c("#E86FAA", "#737884", "#22232B", "#6167AF", "#F2B2D1", "#737870", "#B3A9D3", "#D8DAEC", "#22232A", "#2A348B","#606B89","#444A5E","#985D25","#9595C9","#046735","#778866", "#643F18"))

library(ggplot2)
DimPlot_scCustom(subgenes_B6_rename, reduction = "umap", label = FALSE, pt.size = 1.5,
                 colors_use = c("#E86FAA", "#737884", "#22232B", "#6167AF", "#F2B2D1", "#737870", "#B3A9D3", "#D8DAEC", "#22232A", "#2A348B","#606B89","#444A5E","#985D25","#9595C9","#046735","#778866", "#643F18")) +
  theme(legend.position = "none")
