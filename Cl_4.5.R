library(Seurat)
library(patchwork)
library(dplyr)

Cl1_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl1_4.5")
Cl1_4.5 <- CreateSeuratObject(counts = Cl1_4.5.data, project = "Cl_4.5")
Cl2_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl2_4.5")
Cl2_4.5 <- CreateSeuratObject(counts = Cl2_4.5.data, project = "Cl_4.5")
Cl_4.5 <- merge(Cl1_4.5, y = Cl2_4.5, add.cell.ids = c("1","2"), project = "Cl_4.5")

max = max(Cl_4.5$nFeature_RNA)
gene_percent = (Cl_4.5$nFeature_RNA)/max * 100
Cl_4.5[["gene_percent"]] = gene_percent

Cl_4.5.new <- subset(Cl_4.5, subset = gene_percent > 0.02 & gene_percent < 99.8)
Cl_4.5.new <- NormalizeData(Cl_4.5.new)
Cl_4.5.new <- FindVariableFeatures(Cl_4.5.new, selection.method = "mean.var.plot")
all.genes <- rownames(Cl_4.5.new)

Cl_4.5.new <- ScaleData(Cl_4.5.new, features = all.genes)
Cl_4.5.new <- RunPCA(Cl_4.5.new)
Cl_4.5.new <- RunTSNE(Cl_4.5.new)

DimPlot(Cl_4.5.new, reduction = "tsne")

Cl_4.5.new <- FindNeighbors(Cl_4.5.new)
Cl_4.5.new <- FindClusters(Cl_4.5.new, resolution = 0.1)
DimPlot(Cl_4.5.new, reduction = "tsne")

Cl_4.5.markers <- FindAllMarkers(Cl_4.5.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cl_4.5.markers%>%
  group_by(cluster)%>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(Cl_4.5.new, features = top10$gene)
FeaturePlot(Cl_4.5.new, features = c("Tcf7", "Ccr7", "Gzmb", "Mt1", "Tox"))

cluster0.markers <- FindMarkers(Cl_4.5.new, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(Cl_4.5.new, ident.1 = 1, ident.2 = c(0), min.pct = 0.25)
cluster2.markers <- FindMarkers(Cl_4.5.new, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25)
cluster3.markers <- FindMarkers(Cl_4.5.new, ident.1 = 3, ident.2 = c(0,1,2), min.pct = 0.25)

new.cluster.ids <- c("Tox–/–", "Progenitor-Like CD8+", "2", "3")
names(new.cluster.ids) <- levels(Cl_4.5.new)
Cl_4.5.new <- RenameIdents(Cl_4.5.new, new.cluster.ids)
DimPlot(Cl_4.5.new, reduction = "tsne")