library(Seurat)
library(patchwork)
library(dplyr)

Arm1_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm1_4.5")
Arm1_4.5 <- CreateSeuratObject(counts = Arm1_4.5.data, project = "GSE119943")
Arm2_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm2_4.5")
Arm2_4.5 <- CreateSeuratObject(counts = Arm2_4.5.data, project= "GSE119943")
GSE119943 <- merge(Arm1_4.5, y = Arm2_4.5, add.cell.ids = c("1","2"), project = "GSE119943")

Arm1_7.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm1_7")
Arm1_7 <- CreateSeuratObject(counts = Arm1_7.data, project = "Arm_7")
Arm2_7.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm2_7")
Arm2_7 <- CreateSeuratObject(counts = Arm2_7.data, project = "Arm_7")
Arm_7 <- merge(Arm1_7, y= Arm2_7, add.cell.ids = c("1","2"), project = "Arm_7")

Cl1_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl1_4.5")
Cl1_4.5 <- CreateSeuratObject(counts = Cl1_4.5.data, project = "Cl_4.5")
Cl2_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl2_4.5")
Cl2_4.5 <- CreateSeuratObject(counts = Cl2_4.5.data, project = "Cl_4.5")
Cl_4.5 <- merge(Cl1_4.5, y = Cl2_4.5, add.cell.ids = c("1","2"), project = "Cl_4.5")

Cl1_7.data = Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl1_7")
Cl1_7 <- CreateSeuratObject(counts = Cl1_7.data, project = "Cl_7")
Cl2_7.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Cl2_7")
Cl2_7 <- CreateSeuratObject(counts = Cl2_7.data, project = "Cl_7")
Cl_7 <- merge(Cl1_7, y = Cl2_7, add.cell.ids = c("1","2"), project = "Cl_7")

GSE119943 <- merge(GSE119943, y = c(Arm_7, Cl_4.5, Cl_7), add.cell.ids = c("GSE119943", "Arm_7", "Cl_4.5", "Cl_7"), project = "GSE119943")
max = max(GSE119943$nFeature_RNA)
gene_percent = (GSE119943$nFeature_RNA)/max * 100
GSE119943[["gene_percent"]] = gene_percent

GSE119943.new <- subset(GSE119943, subset = gene_percent > 0.02 & gene_percent < 99.8)
GSE119943.new <- NormalizeData(GSE119943.new)
GSE119943.new <- FindVariableFeatures(GSE119943.new, selection.method = "mean.var.plot")
all.genes <- rownames(GSE119943.new)
GSE119943.new <- ScaleData(GSE119943.new, features = all.genes)
GSE119943.new <- RunPCA(GSE119943.new)
GSE119943.new <- RunTSNE(GSE119943.new)

DimPlot(GSE119943.new, reduction = "tsne")

GSE119943.new <- FindNeighbors(GSE119943.new)
GSE119943.new <- FindClusters(GSE119943.new, resolution = 0.35)
DimPlot(GSE119943.new, reduction = "tsne")

GSE119943.markers <- FindAllMarkers(GSE119943.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
GSE119943.markers%>%
  group_by(cluster)%>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(GSE119943.new, features = top10$gene)
FeaturePlot(GSE119943.new, features = c("Tcf7", "Ccr7", "Gzmb", "Mt1", "Tox"))

cluster0.markers <- FindMarkers(GSE119943.new, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(GSE119943.new, ident.1 = 1, ident.2 = c(0), min.pct = 0.25)
cluster2.markers <- FindMarkers(GSE119943.new, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25)
cluster3.markers <- FindMarkers(GSE119943.new, ident.1 = 3, ident.2 = c(0,1,2), min.pct = 0.25)
cluster4.markers <- FindMarkers(GSE119943.new, ident.1 = 4, ident.2 = c(0,1,2,3), min.pct = 0.25)
cluster5.markers <- FindMarkers(GSE119943.new, ident.1 = 5, ident.2 = c(0,1,2,3,4), min.pct = 0.25)
cluster6.markers <- FindMarkers(GSE119943.new, ident.1 = 6, ident.2 = c(0,1,2,3,4,5), min.pct = 0.25)
cluster7.markers <- FindMarkers(GSE119943.new, ident.1 = 7, ident.2 = c(0,1,2,3,4,5,6), min.pct = 0.25)
cluster8.markers <- FindMarkers(GSE119943.new, ident.1 = 8, ident.2 = c(0,1,2,3,4,5,6,7), min.pct = 0.25)
cluster9.markers <- FindMarkers(GSE119943.new, ident.1 = 9, ident.2 = c(0,1,2,3,4,5,6,7,8), min.pct = 0.25)
cluster10.markers <- FindMarkers(GSE119943.new, ident.1 = 10, ident.2 = c(0,1,2,3,4,5,6,7,8,9), min.pct = 0.25)
