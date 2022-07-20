library(dplyr)
library(Seurat)
library(patchwork)

sessionInfo()

name_list <- c("MI_D1_IR_1", "MI_D1_IR_2", "MI_D1_RR_1", "MI_D1_RR_2",
               "MI_D7_IR_1", "MI_D7_IR_2", "MI_D7_RR_1", "MI_D7_RR_2",
               "MI_D14_IR_1", "MI_D14_IR_2", "MI_D14_RR_1", "MI_D14_RR_2",
               "Sham_1", "Sham_2")

data_1 <- Read10X("WT_MI_D1_IR_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[1], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[1])
data_2 <- Read10X("WT_MI_D1_IR_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[2], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[2])
data_3 <- Read10X("WT_MI_D1_RR_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[3], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[3])
data_4 <- Read10X("WT_MI_D1_RR_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[4], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[4])
data_5 <- Read10X("WT_MI_D7_IR_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[5], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[5])
data_6 <- Read10X("WT_MI_D7_IR_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[6], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[6])
data_7 <- Read10X("WT_MI_D7_RR_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[7], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[7])
data_8 <- Read10X("WT_MI_D7_RR_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[8], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[8])
data_9 <- Read10X("WT_MI_D14_IR_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[9], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[9])
data_10 <- Read10X("WT_MI_D14_IR_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[10], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[10])
data_11 <- Read10X("WT_MI_D14_RR_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[11], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[11])
data_12 <- Read10X("WT_MI_D14_RR_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[12], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[12])
data_13 <- Read10X("WT_Sham_1/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[13], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[13])
data_14 <- Read10X("WT_Sham_2/outs/filtered_feature_bc_matrix") %>% CreateSeuratObject(project = name_list[14], min.cells = 0, min.features = 0) %>% RenameCells(add.cell.id = name_list[14])

all <- merge(data_1, c(data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14))

all[["percent.mt"]] <- PercentageFeatureSet(all, pattern = "^mt-")

mt.index <- grep(pattern = "^mt-", x = rownames(all[["RNA"]]), value = FALSE)
all_matrix <- all[["RNA"]][-mt.index, ]
all <- CreateSeuratObject(counts = all_matrix, meta.data = all@meta.data)

all_list <- SplitObject(all, split.by = "orig.ident")

all_list <- lapply(X = all_list, FUN = function(x) {
    x <- subset(x, subset = nFeature_RNA > 500 & percent.mt < 60)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = all_list)

all_anchors <- FindIntegrationAnchors(object.list = all_list, anchor.features = features)

combined <- IntegrateData(anchorset = all_anchors)
combined <- ScaleData(combined, features = rownames(all)) %>% RunPCA(npcs = 30, features = VariableFeatures(combined)) %>% RunUMAP(reduction = "pca", dims = 1:20, umap.method = 'umap-learn') %>% FindNeighbors(reduction = "pca", dims = 1:20) 

combined <- FindClusters(combined, resolution = 0.25)

DimPlot(combined, reduction = "umap", pt.size = 1.5)

DefaultAssay(combined) <- "RNA"

CM <- subset(combined, idents = c("0"))
CM <- FindVariableFeatures(CM, selection.method = "vst", nfeatures = 500) %>% ScaleData(features = rownames(CM))
CM <- RunPCA(CM, npcs = 30, features = VariableFeatures(CM)) %>% RunUMAP(reduction = "pca", dims = 1:5, umap.method = 'umap-learn') %>% FindNeighbors(reduction = "pca", dims = 1:5)
CM <- FindClusters(CM, resolution = 0.2)

DimPlot(CM, reduction = "umap", pt.size = 1.5)

cluster0_marker <- FindMarkers(CM, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
cluster1_marker <- FindMarkers(CM, ident.1 = 1, only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
cluster2_marker <- FindMarkers(CM, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
