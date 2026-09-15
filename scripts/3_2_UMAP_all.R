R.Version()
# R version 4.3.3 (2024-02-29)
# Creates a SingleCellExperiment.

# ============== UMAP. ==============

getwd()   # "/home/data/SpermatogenicCells"

library(Seurat)
library(dplyr)


sessionInfo(package = "Seurat")
# Seurat_5.0.1

seura_do26386_P5 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do26386/outs/filtered_feature_bc_matrix"), project = "do26386")
str(seura_do26386_P5)    # 32285 5796. rows columns. genes x cells/UMIs/barcodes.
# [1:5796] "AAACCTGAGAAACGAG-1" "AAACCTGAGCGATATA-1" "AAACCTGAGGAATTAC-1"

seura_do26387_P5 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do26387/outs/filtered_feature_bc_matrix"), project = "do26387")
str(seura_do26387_P5)  # 32285 8481.
# chr [1:8481] "AAACCTGAGAGCTATA-1" "AAACCTGAGCGATGAC-1" "AAACCTGAGGCGTACA-1"

seura_do17821_P10 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17821/outs/filtered_feature_bc_matrix"), project = "do17821")
str(seura_do17821_P10)  # 32285 4616.

seura_do18195_P15 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do18195/outs/filtered_feature_bc_matrix"), project = "do18195")
str(seura_do18195_P15)  # 32285 16021.

seura_do17824_P20 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17824/outs/filtered_feature_bc_matrix"), project = "do17824")
str(seura_do17824_P20)  # 32285 8768.

seura_do18196_P25 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do18196/outs/filtered_feature_bc_matrix"), project = "do18196")
str(seura_do18196_P25)  # 32285 7652

seura_do17825_P30 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17825/outs/filtered_feature_bc_matrix"), project = "do17825")
str(seura_do17825_P30)  # 32285 5543.

seura_do17827_P35 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17827/outs/filtered_feature_bc_matrix"), project = "do17827")
str(seura_do17827_P35)  # 32285 8451.

seura_do17815_B6 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17815/outs/filtered_feature_bc_matrix"), project = "do17815")
str(seura_do17815_B6)  # 32285 4036.

seura_do17816_B6 <- CreateSeuratObject(Read10X("/home/data/SpermatogenicCells/aligned/do17816/outs/filtered_feature_bc_matrix"), project = "do17816")
str(seura_do17816_B6)  # 32285 5131.

merge_all <- merge(x = seura_do26386_P5, y = c(seura_do26387_P5, seura_do17821_P10, seura_do18195_P15,
                                               seura_do17824_P20, seura_do18196_P25, seura_do17825_P30,
                                               seura_do17827_P35, seura_do17815_B6, seura_do17816_B6),
                   add.cell.ids = c("do26386", "do26387", "do17821", "do18195", "do17824", "do18196",
                                    "do17825", "do17827", "do17815", "do17816"), project = "sperm")
# counts.1: 32285 5796
# counts.2: 32285 8481
# cells chr [1:74495] "do26386_AAACCTGAGAAACGAG-1" "do26386_AAACCTGAGCGATATA-1"
# 5796+8481+4616+16021+8768+7652+5543+8451+4036+5131=74495.

unique(sapply(X = strsplit(colnames(merge_all), split = "_"), FUN = "[", 1))
table(merge_all$orig.ident)
save(merge_all, file = "seura_merge_all_sperm.Rdata")

# load("/home/data/SpermatogenicCells/seura_merge_all_sperm.Rdata")
merge_all_new = JoinLayers(object = merge_all)

# keep genes with at least 3 cells with 1 or more counts.
# to filter genes. got filtered genes. equal to "min.cells = 3".
rs <- rowSums(GetAssayData(object = merge_all_new, assay="RNA",layer ="counts")>0)
length(rs) # 32285.

# fe <- rs[which(rs>10)]
# length(fe) # 26491. min.cells=10.
fe <- rs[which(rs>3)]
length(fe)  # 27918. min.cells=3.

subgenes <- subset(merge_all,features=names(fe))
str(subgenes)
# chr [1:26491] "Xkr4" "Gm1992" "Gm19938" (min.cells=10).  # chr [1:27918] "Xkr4" "Gm1992" "Gm19938" (min.cells=3).
# chr [1:74495] "do26386_AAACCTGAGAAACGAG-1" "do26386_AAACCTGAGCGATATA-1".
# separated 10 count matrix. "$ counts.do26386". "$ counts.do26387".
# @ project.name: chr "sperm".

subgenes1 <- subset(merge_all_new,features=names(fe))
str(subgenes1)
# chr [1:26491] "Xkr4" "Gm1992" (min.cells=10).  # chr [1:27918] "Xkr4" "Gm1992" (min.cells=3).
# chr [1:74495] "do26386_AAACCTGAGAAACGAG-1" "do26386_AAACCTGAGCGATATA-1". 
# integrated full count matrix of 10 samples. "$ counts".
# @ project.name: chr "sperm".


# try subgenes here.
subgenes[["percent.mt"]] <- PercentageFeatureSet(subgenes, pattern = "^MT-")
# dim(subgenes[["percent.mt"]]) #74495     1


# calculate blood cells.
HB.genes<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(subgenes@assays$RNA)) 
HB.genes <- rownames(subgenes@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
subgenes[["percent.HB"]]<-PercentageFeatureSet(subgenes, features=HB.genes) 

col.num <- length(levels(subgenes@active.ident))


VlnPlot(subgenes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB"), ncol = 4)
Plot1AB <- FeatureScatter(subgenes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");Plot1AB

subgenes_filtered<- subset(subgenes, subset = nFeature_RNA > 2000 & nCount_RNA > 10000 & percent.mt < 5)
str(subgenes_filtered)
# chr [1:21774] cells.  # before filter. chr [1:74495] cells.
# chr [1:27918] "Xkr4" "Gm1992" "Gm19938".
# chr [1:21774] "do26386_AAACCTGAGGAATTAC-1" "do26386_AAACCTGAGTACGCGA-1".
table(subgenes_filtered$orig.ident) 

VlnPlot(subgenes_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Plot1AB_filtered <- FeatureScatter(subgenes_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA");Plot1AB_filtered

subgenes_filtered_norm <- NormalizeData(subgenes_filtered, normalization.method = "LogNormalize", scale.factor = 100000)
# here, nCount_RNA = a*10^5. so scale.factor = 10^5 = 100000

subgenes_varFeatures <- FindVariableFeatures(subgenes_filtered_norm, selection.method = "vst", nfeatures = 2000);subgenes_varFeatures
# 26491 features across 21774 samples within 1 assay 
# Active assay: RNA (26491 features, 2000 variable features).
# min.cells=10.

# 27918 features across 21774 samples within 1 assay 
# Active assay: RNA (27918 features, 2000 variable features)
# min.cells=3.
str(subgenes_varFeatures)

top10 <- head(VariableFeatures(subgenes_varFeatures), 10);top10  # Get and set variable feature information.
plot03A <- VariableFeaturePlot(subgenes_varFeatures);plot03A     # plot variable features with and without labels.
plot03B <- LabelPoints(plot = plot03A, points = top10, repel = TRUE, xnudge = 0, ynudge = 0);plot03B

# ---Scale the data.---
all.genes <- rownames(subgenes_varFeatures);length(all.genes) # 27918. min.cells=3.
subgenes_varFeatures_Scacle <- ScaleData(subgenes_varFeatures, features = all.genes, 
                                         do.scale = TRUE, do.center = TRUE,
                                         model.use = "linear", scale.max = 10)   # All default here.

# ---Perform linear dimensional reduction.---
subgenes_vs_PCA <- RunPCA(subgenes_varFeatures_Scacle, features = VariableFeatures(object = subgenes_varFeatures_Scacle), npcs = 60, rev.pca = FALSE)
class(subgenes_vs_PCA) # "SeuratObject".
str(subgenes_vs_PCA)


# examine and visualize PCA results a few different ways.
print(subgenes_vs_PCA[["pca"]], dims = 1:12, nfeatures = 10)
VizDimLoadings(subgenes_vs_PCA, dims = 1:6, reduction = "pca", col = "purple", nfeatures = 10)
DimPlot(subgenes_vs_PCA, reduction = "pca")
DimHeatmap(subgenes_vs_PCA, dims = 1, cells = 1000, nfeatures = 20, balanced = TRUE, slot = "scale.data")
DimHeatmap(subgenes_vs_PCA, dims = 1:15, cells = 1000, balanced = TRUE)
DimHeatmap(subgenes_vs_PCA, dims = 1:6, cells = 1000, balanced = TRUE)


# ---Determine the 'dimensionality' of the dataset.---
subgenes_PCs <- JackStraw(subgenes_vs_PCA, dims = 30, num.replicate = 200, maxit = 2000)
subgenes_PCs_score <- ScoreJackStraw(subgenes_PCs, dims = 1:30)
JackStrawPlot(subgenes_PCs_score, dims = 1:30)
ElbowPlot(subgenes_PCs_score)

# ---Cluster the cells.---
subgenes_SNN <- FindNeighbors(subgenes_PCs_score, dims = 1:30, k.param = 20)
subgenes_SNN_Cluster <- FindClusters(subgenes_SNN, resolution = 0.5)
head(Idents(subgenes_SNN_Cluster), 5)  # Look at cluster IDs of the first 5 cells.
tail(Idents(subgenes_SNN_Cluster), 5)

# ---Run non-linear dimensional reduction (UMAP/tSNE).---
subgenes_UMAP <- RunUMAP(subgenes_SNN_Cluster, reduction = "pca", dims = 1:60,  n.neighbors = 40L)
DimPlot(subgenes_UMAP, reduction = "umap")
DimPlot(subgenes_UMAP, reduction = "umap", label = TRUE, pt.size = 1)
levels(subgenes_UMAP)
# 26 clusters. min.cells=10.
# save(subgenes_UMAP, file = "subgenes_UMAP.Rdata")

# 26 clusters. min.cells=3.
save(subgenes_UMAP, file = "subgenes_UMAP1.Rdata")

library(scCustomize)
DimPlot_scCustom(subgenes_UMAP, reduction = "umap", label = TRUE, pt.size = 1,
                 colors_use = c("#643F18", "#985D25", "#046735",
                                "#FBDCEA", "#F2B2D1", "#E86FAA", "#D31281", "#9A1C59", "#7E2678", "#502269", "#240F2F",
                                "#B3A9D3", "#D8DAEC", "#9595C9", "#6167AF", "#2A348B", "#282B69", "#171447",
                                "#22232B", "#444A5E", "#606B89", "#737884",
                                "wheat1", "lightgoldenrod1", "yellow1", "gold1"))

# ---Find differentially expressed features (cluster biomarkers).---
subgenes_UMAP_new = JoinLayers(object = subgenes_UMAP)
str(subgenes_UMAP_new)  # "SeuratObject".
levels(subgenes_UMAP_new)


# # find all markers of cluster 1  # Marker = differentially expressed genes
# cluster1.markers <- FindMarkers(subgenes_UMAP_new, ident.1 = 1, min.pct = 0.25);head(cluster1.markers, n = 5)   # test.use = "wilcox"  only.pos = FALSE
# cluster1.markers <- FindMarkers(subgenes_UMAP_new, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE);head(cluster1.markers, n = 5)
# 
# # find all markers distinguishing cluster 5 from clusters 0 and 3
# cluster5.markers <- FindMarkers(subgenes_UMAP_new, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25);head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones.
subgenes.markers <- FindAllMarkers(subgenes_UMAP_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log2(2));class(subgenes.markers);dim(subgenes.markers)
# 55483     7.
head(subgenes.markers)
tail(subgenes.markers)
save(subgenes.markers, file = "subgenes.markers.Rdata")

# to see all the markers with top 2 avg_logFC in all clusters.
subgenes.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) # 26*2=52.

top10 <- subgenes.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);str(top10)   # 26*10=260 genes.
DoHeatmap(subgenes_UMAP_new, features = top10$gene[1:20])  #+ NoLegend()


# feature/count visualization.
# Somatic cells.
VlnPlot(subgenes_UMAP_new, features = c("Cldn11", "Fabp3"))
VlnPlot(subgenes_UMAP_new, features = c("Cldn11", "Fabp3"), slot = "counts", log = TRUE)  # you can plot raw counts as well
# FeaturePlot(subgenes_UMAP_new, features = c("Cldn11", "Fabp3"), pt.size = 1)
FeaturePlot(subgenes_UMAP_new, features = c("Cldn11"), pt.size = 1)
FeaturePlot(subgenes_UMAP_new, features = c("Fabp3"), pt.size = 1)

# Spermatogonia.
VlnPlot(subgenes_UMAP_new, features = c("Dmrt1"))
VlnPlot(subgenes_UMAP_new, features = c("Dmrt1"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_UMAP_new, features = c("Dmrt1"), pt.size = 1)

# Spermatocytes.
VlnPlot(subgenes_UMAP_new, features = c("Piwil1"))
VlnPlot(subgenes_UMAP_new, features = c("Piwil1"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_UMAP_new, features = c("Piwil1"), pt.size = 1)

# Round spermatids.
VlnPlot(subgenes_UMAP_new, features = c("Tex21"))
VlnPlot(subgenes_UMAP_new, features = c("Tex21"), slot = "counts", log = TRUE)
FeaturePlot(subgenes_UMAP_new, features = c("Tex21"), pt.size = 1)

# Elongating spermatids.
VlnPlot(subgenes_UMAP_new, features = c("Tnp1"))
VlnPlot(subgenes_UMAP_new, features = c("Tnp1"), slot = "counts", log = TRUE)  # you can plot raw counts as well
FeaturePlot(subgenes_UMAP_new, features = c("Tnp1"), pt.size = 1)


# ---Assign cell type identity to clusters.---
new.cluster.ids <- c("Spermatocytes2", "Elongating spermatids1", "Spermatocytes3", "Round spermatids0", "Somatic cells3", "Spermatocytes1",
                     "Spermatocytes5", "Elongating spermatids0", "Spermatogonia", "Round spermatids2", "Round spermatids1",
                     "Somatic cells0","Spermatocytes8","Spermatocytes4","Spermatocytes10", "Somatic cells4",
                     "Spermatocytes7", "Spermatocytes0", "Somatic cells5", "Somatic cells6", "Somatic cells7",
                     "Spermatocytes9", "Spermatocytes6", "Somatic cells2", "Somatic cells1", "Somatic cells8")
names(new.cluster.ids) <- levels(subgenes_UMAP_new)
subgenes_rename <- RenameIdents(subgenes_UMAP_new, new.cluster.ids)
save(subgenes_rename, file = "subgenes_rename.Rdata")
# which will be used as seuratobject for velocity analysis later.

DimPlot(subgenes_rename, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

DimPlot_scCustom(subgenes_rename, reduction = "umap", label = TRUE, pt.size = 1,
                 colors_use = c("#E86FAA", "#737884", "deeppink2", "#B3A9D3", "goldenrod3", "#F2B2D1",
                                "#9A1C59", "#22232B", "#046735", "#282B69", "#6167AF",
                                "#643F18","mediumorchid4","#D31281","#502269", "goldenrod2",
                                "plum2", "#FBDCEA", "goldenrod1", "yellow1", "khaki1",
                                "#7E2678", "mediumorchid1", "goldenrod4", "#985D25", "lemonchiffon"))

levels(Idents(subgenes_rename))
my_cols <- c("Spermatocytes2"="#E86FAA", "Elongating spermatids1"="#737884", "Spermatocytes3"="deeppink2", "Round spermatids0"="#B3A9D3", "Somatic cells3"="goldenrod3", "Spermatocytes1"="#F2B2D1",
             "Spermatocytes5"="#9A1C59", "Elongating spermatids0"="#22232B", "Spermatogonia"="#046735", "Round spermatids2"="#282B69", "Round spermatids1"="#6167AF",
             "Somatic cells0"="#643F18","Spermatocytes8"="mediumorchid4","Spermatocytes4"="#D31281","Spermatocytes10"="#502269", "Somatic cells4"="goldenrod2",
             "Spermatocytes7"="plum2", "Spermatocytes0"="#FBDCEA", "Somatic cells5"="goldenrod1", "Somatic cells6"="yellow1", "Somatic cells7"="khaki1",
             "Spermatocytes9"="#7E2678", "Spermatocytes6"="mediumorchid1", "Somatic cells2"="goldenrod4", "Somatic cells1"="#985D25", "Somatic cells8"="lemonchiffon")

my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)

DimPlot(subgenes_rename, reduction = "umap", cols = my_cols2, label = TRUE, pt.size = 1) + NoLegend()

UMAPPlot(subgenes_rename, cols = my_cols2, pt.size = 1, label = FALSE) # used this one.



