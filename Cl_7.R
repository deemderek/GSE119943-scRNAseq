library(Seurat)
library(patchwork)
library(dplyr)

Cl1_7.data = Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl1_7")
Cl1_7 <- CreateSeuratObject(counts = Cl1_7.data, project = "Cl_7")
Cl2_7.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl2_7")
Cl2_7 <- CreateSeuratObject(counts = Cl2_7.data, project = "Cl_7")

Cl_7 <- merge(Cl1_7, y = Cl2_7, add.cell.ids = c("1","2"), project = "Cl_7")
max = max(Cl_7$nFeature_RNA)
gene_percent = (Cl_7$nFeature_RNA)/max * 100
Cl_7[["gene_percent"]] = gene_percent
Cl_7.new <- subset(Cl_7, subset = gene_percent > 0.02 & gene_percent < 99.8)

Cl_7.new <- NormalizeData(Cl_7.new)
Cl_7.new <- FindVariableFeatures(Cl_7.new, selection.method = "mean.var.plot")
all.genes <- rownames(Cl_7.new)
Cl_7.new <- ScaleData(Cl_7.new, features = all.genes)
Cl_7.new <- RunPCA(Cl_7.new)
Cl_7.new <- RunTSNE(Cl_7.new)

Cl_7.new <- FindNeighbors(Cl_7.new)
Cl_7.new <- FindClusters(Cl_7.new, resolution = 0.1)

DimPlot(Cl_7.new, reduction = "tsne")
FeaturePlot(Cl_7.new, features = c("Tcf7", "Ccr7", "Gzmb", "Mt1", "Tox"))

Cl_7.markers <- FindAllMarkers(Cl_7.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cl_7.markers %>%
  group_by(cluster)%>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(Cl_7.new, features = top10$gene)

cluster0.markers <- FindMarkers(Cl_7.new, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(Cl_7.new, ident.1 = 1, ident.2 = c(0), min.pct = 0.25)
cluster2.markers <- FindMarkers(Cl_7.new, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25)
cluster3.markers <- FindMarkers(Cl_7.new, ident.1 = 3, ident.2 = c(0,1,2), min.pct = 0.25)

new.cluster.ids <- c("0", "Memory T Cells?", "Gamma Delta", "Natural Killer")
names(new.cluster.ids) <- levels(Cl_7.new)
Cl_7.new <- RenameIdents(Cl_7.new, new.cluster.ids)
DimPlot(Cl_7.new, reduction = "tsne")
