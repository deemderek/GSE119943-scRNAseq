library(Seurat)
library(patchwork)
library(dplyr)

Arm1_4.5.data = Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm1_4.5")
Arm1_4.5 <- CreateSeuratObject(counts = Arm1_4.5.data, project = "Arm_4.5")
Arm2_4.5.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm2_4.5")
Arm2_4.5 <- CreateSeuratObject(counts = Arm2_4.5.data, project = "Arm_4.5")

Arm_4.5 <- merge(Arm1_4.5, y = Arm2_4.5, add.cell.ids = c("1","2"), project = "Arm_4.5")
max = max(Arm_4.5$nFeature_RNA)
gene_percent = (Arm_4.5$nFeature_RNA)/max * 100
Arm_4.5[["gene_percent"]] = gene_percent
Arm_4.5.new <- subset(Arm_4.5, subset = gene_percent > 0.02 & gene_percent < 99.8)

Arm_4.5.new <- NormalizeData(Arm_4.5.new)
Arm_4.5.new <- FindVariableFeatures(Arm_4.5.new, selection.method = "mean.var.plot")
all.genes <- rownames(Arm_4.5.new)
Arm_4.5.new <- ScaleData(Arm_4.5.new, features = all.genes)
Arm_4.5.new <- RunPCA(Arm_4.5.new)
Arm_4.5.new <- RunTSNE(Arm_4.5.new)

Arm_4.5.new <- FindNeighbors(Arm_4.5.new)
Arm_4.5.new <- FindClusters(Arm_4.5.new, resolution = 0.12)
DimPlot(Arm_4.5.new, reduction = "tsne")

Arm_4.5.markers <- FindAllMarkers(Arm_4.5.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Arm_4.5.markers%>%
  group_by(cluster)%>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(Arm_4.5.new, features = top10$gene)
FeaturePlot(Arm_4.5.new, features = c("Tcf7", "Ccr7", "Gzmb", "Mt1", "Tox"))

cluster0.markers <- FindMarkers(Arm_4.5.new, ident.1 = 0, min.pct = 0.25)
cluster1.markers <- FindMarkers(Arm_4.5.new, ident.1 = 1, ident.2 = c(0), min.pct = 0.25)
cluster2.markers <- FindMarkers(Arm_4.5.new, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25)
cluster3.markers <- FindMarkers(Arm_4.5.new, ident.1 = 3, ident.2 = c(0,1,2), min.pct = 0.25)

new.cluster.ids <- c("0", "1", "Progenitor-Like CD8+", "3")
names(new.cluster.ids) <- levels(Arm_4.5.new)
Arm_4.5.new <- RenameIdents(Arm_4.5.new, new.cluster.ids)
DimPlot(Arm_4.5.new, reduction = "tsne")
