# ============== Velocity analysis. ==============

getwd()   # "/home/data/SpermatogenicCells"

library(velocyto.R)

load("/home/data/SpermatogenicCells/subgenes_B6_rename1.Rdata")

adat<-read.loom.matrices("/home/data/SpermatogenicCells/B6_merged.loom")
save(adat, file = "adat.Rdata")

# load("/home/data/SpermatogenicCells/adat.Rdata")
adata <- adat

colnames(adata$spliced)<-paste0(colnames(adata$spliced),"-1")
colnames(adata$spliced)<-gsub("possorted_genome_bam_WJX2O:","do17815_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_O2D4B:","do17816_",colnames(adata$spliced))

colnames(adata$unspliced)<-colnames(adata$spliced)
colnames(adata$ambiguous)<-colnames(adata$spliced)

adata$spliced<-adata$spliced[,rownames(subgenes_B6_rename@meta.data)]
adata$unspliced<-adata$unspliced[,rownames(subgenes_B6_rename@meta.data)]
adata$ambiguous<-adata$ambiguous[,rownames(subgenes_B6_rename@meta.data)]

adata$spliced <- adata$spliced[rownames(unlist(subgenes_B6_rename@assays$RNA@features@.Data))[
  c(1:108, 110:1132, 1134:4837, 4839:7337, 7339:10230,
    10232:11551, 11553:13468, 13470:13730, 13732:16056,
    16058:17015, 17017:19307, 19309:19315, 19317:20096,
    20098:22473, 22475:23315, 23317:24157, 24159:24700)],]
adata$unspliced <- adata$unspliced[rownames(unlist(subgenes_B6_rename@assays$RNA@features@.Data))[
  c(1:108, 110:1132, 1134:4837, 4839:7337, 7339:10230,
    10232:11551, 11553:13468, 13470:13730, 13732:16056,
    16058:17015, 17017:19307, 19309:19315, 19317:20096,
    20098:22473, 22475:23315, 23317:24157, 24159:24700)],]
adata$ambiguous <- adata$ambiguous[rownames(unlist(subgenes_B6_rename@assays$RNA@features@.Data))[
  c(1:108, 110:1132, 1134:4837, 4839:7337, 7339:10230,
    10232:11551, 11553:13468, 13470:13730, 13732:16056,
    16058:17015, 17017:19307, 19309:19315, 19317:20096,
    20098:22473, 22475:23315, 23317:24157, 24159:24700)],]

sp<-adata$spliced; dim(sp)
# 24684  3861.
unsp<-adata$unspliced; dim(unsp) 
# 24684  3861.
umap<-subgenes_B6_rename@reductions$umap@cell.embeddings


# calculate the distance between cells.
cell.dist <- as.dist(1-armaCor(t(subgenes_B6_rename@reductions$umap@cell.embeddings)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(sp,unsp,deltaT=2,kCells=10, cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores=48)
save(rvel.cd, file = "rvel_cd1.Rdata")

load("/home/data/SpermatogenicCells/rvel_cd1.Rdata")
library(ggplot2)
library(Seurat)
library(scCustomize)

pdf("/cell_velocity.pdf",height=6,width=8)
levels(Idents(subgenes_B6_rename))
my_cols <-c("Spermatocytes1"="#E86FAA", "Elongating spermatids0"="#737884", "Elongating spermatids1"="#22232B", "Round spermatids1"="#6167AF", 
            "Spermatocytes0"="#F2B2D1", "Elongating spermatids2"="#737870", "Round spermatids2"="#B3A9D3", "Round spermatids3"="#D8DAEC", 
            "Elongating spermatids3"="#22232A", "Round spermatids4"="#2A348B","Elongating spermatids4"="#606B89","Elongating spermatids5"="#444A5E",
            "Somatic cells0"="#985D25","Round spermatids5"="#9595C9","Spermatogonia"="#046735","Elongating spermatids6"="#778866", 
            "Somatic cells1"="#643F18")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)
gg <- UMAPPlot(subgenes_B6_rename, cols = my_cols2, pt.size = 1.5, label = FALSE)

colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(umap)
p1 <- show.velocity.on.embedding.cor(umap,rvel.cd,n=400,scale='sqrt',cell.colors=ac(colors,alpha=0.5),cex=0.8,
                                     arrow.scale=2,show.grid.flow=T,min.grid.cell.mass=1.0,grid.n=50,
                                     do.par=F,cell.border.alpha =0.1,n.cores=48,main="B6 Cell Velocity",
                                     max.grid.arrow.length=1, fixed.arrow.length=FALSE, arrow.lwd=1)
dev.off()
