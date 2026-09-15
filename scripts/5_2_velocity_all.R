# ============== Velocity analysis. ==============

getwd()   # "/home/data/SpermatogenicCells"

library(velocyto.R)

load("/home/data/SpermatogenicCells/subgenes_rename.Rdata")
str(subgenes_rename)
# $ counts: 27918 21774.
# chr [1:27918] "Xkr4" "Gm1992" "Gm19938". genes.
# chr [1:21774] "do26386_AAACCTGAGGAATTAC-1" "do26386_AAACCTGAGTACGCGA-1". cells.

allten_data <- read.loom.matrices("/home/data/SpermatogenicCells/allten_merged.loom")
save(allten_data, file = "allten_data.Rdata")

# load("/home/data/SpermatogenicCells/allten_data.Rdata")
adata <- allten_data
str(adata)


colnames(adata$spliced)<-paste0(colnames(adata$spliced),"-1")
colnames(adata$spliced)<-gsub("possorted_genome_bam_WJX2O:","do17815_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_O2D4B:","do17816_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_4DYG5:","do26386_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_67K09:","do26387_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_X4IXK:","do17821_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_18F66:","do18195_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_DA6KO:","do17824_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_A70S7:","do18196_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_8E32G:","do17825_",colnames(adata$spliced))
colnames(adata$spliced)<-gsub("possorted_genome_bam_IKGDL:","do17827_",colnames(adata$spliced))

colnames(adata$unspliced)<-colnames(adata$spliced)
colnames(adata$ambiguous)<-colnames(adata$spliced)

adata$spliced<-adata$spliced[,rownames(subgenes_rename@meta.data)]
dim(adata$spliced) # 32285 21774.
adata$unspliced<-adata$unspliced[,rownames(subgenes_rename@meta.data)]
dim(adata$unspliced) # 32285 21774.
adata$ambiguous<-adata$ambiguous[,rownames(subgenes_rename@meta.data)]
dim(adata$ambiguous)  # 32285 21774.


adata$spliced <- adata$spliced[rownames(unlist(subgenes_rename@assays$RNA@features@.Data))[
  which(rownames(unlist(subgenes_rename@assays$RNA@features@.Data)) %in% rownames(adata$spliced))],]
dim(adata$spliced)  # 27898 21774.
adata$unspliced <- adata$unspliced[rownames(unlist(subgenes_rename@assays$RNA@features@.Data))[
  which(rownames(unlist(subgenes_rename@assays$RNA@features@.Data)) %in% rownames(adata$unspliced))],]
dim(adata$unspliced)  # 27898 21774.
adata$ambiguous <- adata$ambiguous[rownames(unlist(subgenes_rename@assays$RNA@features@.Data))[
  which(rownames(unlist(subgenes_rename@assays$RNA@features@.Data)) %in% rownames(adata$ambiguous))],]
dim(adata$ambiguous)  # 27898 21774.


sp<-adata$spliced; dim(sp) # 27898 21774.
unsp<-adata$unspliced; dim(unsp) # 27898 21774.
umap<-subgenes_rename@reductions$umap@cell.embeddings  # num [1:21774, 1:2].


# calculate the distance between cells.
cell.dist <- as.dist(1-armaCor(t(subgenes_rename@reductions$umap@cell.embeddings)))
class(cell.dist)  # dist.
str(cell.dist)
head(cell.dist)
length(cell.dist)  # 237042651.

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(sp,unsp,deltaT=2,kCells=10, cell.dist=cell.dist,fit.quantile=fit.quantile,n.cores=48)
save(rvel.cd, file = "allten_rvel_cd.Rdata")

# load("/home/data/SpermatogenicCells/allten_rvel_cd.Rdata")

library(ggplot2)
library(Seurat)
library(scCustomize)

# pdf("/cell_velocity.pdf",height=6,width=8)
levels(Idents(subgenes_rename))
my_cols <- c("Spermatocytes2"="#E86FAA", "Elongating spermatids1"="#737884", "Spermatocytes3"="deeppink2", "Round spermatids0"="#B3A9D3", "Somatic cells3"="goldenrod3", "Spermatocytes1"="#F2B2D1",
             "Spermatocytes5"="#9A1C59", "Elongating spermatids0"="#22232B", "Spermatogonia"="#046735", "Round spermatids2"="#282B69", "Round spermatids1"="#6167AF",
             "Somatic cells0"="#643F18","Spermatocytes8"="mediumorchid4","Spermatocytes4"="#D31281","Spermatocytes10"="#502269", "Somatic cells4"="goldenrod2",
             "Spermatocytes7"="plum2", "Spermatocytes0"="#FBDCEA", "Somatic cells5"="goldenrod1", "Somatic cells6"="yellow1", "Somatic cells7"="khaki1",
             "Spermatocytes9"="#7E2678", "Spermatocytes6"="mediumorchid1", "Somatic cells2"="goldenrod4", "Somatic cells1"="#985D25", "Somatic cells8"="lemonchiffon")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
scales::show_col(my_cols2)

gg <- UMAPPlot(subgenes_rename, cols = my_cols2, pt.size = 1.5, label = FALSE);gg

colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
names(colors) <- rownames(umap)
p1 <- show.velocity.on.embedding.cor(umap,rvel.cd,n=400,scale='sqrt',cell.colors=ac(colors,alpha=0.5),cex=0.8,
                                     arrow.scale=2,show.grid.flow=T,min.grid.cell.mass=1.0,grid.n=50,
                                     do.par=F,cell.border.alpha =0.1,n.cores=48,main="P5-P10-...-B6 Cell Velocity",
                                     max.grid.arrow.length=1, fixed.arrow.length=FALSE, arrow.lwd=1)

p2 <- show.velocity.on.embedding.cor(umap,rvel.cd,n=200,scale='sqrt',cell.colors=ac(colors,alpha=0.5),cex=0.8,
                                     arrow.scale=2,show.grid.flow=T,min.grid.cell.mass=1.0,grid.n=50,
                                     do.par=F,cell.border.alpha =0.1,n.cores=48, main = "Cell Velocity",
                                     max.grid.arrow.length=1, fixed.arrow.length=FALSE, arrow.lwd=1)

# dev.off()
