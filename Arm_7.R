library(Seurat)
library(patchwork)
library(dplyr)

Arm1_7.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm1_7")
Arm1_7 <- CreateSeuratObject(counts = Arm1_7.data, project = "Arm_7")
Arm2_7.data <- Read10X(data.dir = "/Users/derekmartinez/Desktop/GSE119043_unzip/Arm2_7")
Arm2_7 <- CreateSeuratObject(counts = Arm2_7.data, project = "Arm_7")

Arm_7 <- merge(Arm1_7, y= Arm2_7, add.cell.ids = c("1","2"), project = "Arm_7")
max = max(Arm_7$nFeature_RNA)
gene_percent = (Arm_7$nFeature_RNA)/max * 100
Arm_7[["gene_percent"]] = gene_percent
Arm_7.new <- subset(Arm_7, subset = gene_percent > 0.2 & gene_percent < 99.8)

Arm_7.new <- NormalizeData(Arm_7.new)
Arm_7.new <- FindVariableFeatures(Arm_7.new, selection.method = "mean.var.plot")
all.genes <- rownames(Arm_7.new)
Arm_7.new <- ScaleData(Arm_7.new, features = all.genes)
Arm_7.new <- RunPCA(Arm_7.new)
Arm_7.new <- RunTSNE(Arm_7.new)

Arm_7.new <- FindNeighbors(Arm_7.new)
Arm_7.new <- FindClusters(Arm_7.new, resolution = 0.12)
DimPlot(Arm_7.new, reduction = "tsne")

Arm_7.markers <- FindAllMarkers(Arm_7.new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Arm_7.markers%>%
  group_by(cluster)%>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(Arm_7.new, features = top10$gene)
FeaturePlot(Arm_7.new, features = c("Tcf7", "Ccr7", "Gzmb", "Mt1", "Tox"))
