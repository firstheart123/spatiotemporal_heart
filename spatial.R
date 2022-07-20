library(dplyr)
library(Seurat)
library(patchwork)
library(WGCNA)

sessionInfo()

name_list <- c("WT_Sham",
               "WT_MI_day1_1", "WT_MI_day1_2", "WT_MI_day1_3",
               "WT_MI_day7_1", "WT_MI_day7_2", "WT_MI_day7_3",
               "WT_MI_day14_1", "WT_MI_day14_2", "WT_MI_day14_3"
               )
               
data_1 <- Load10X_Spatial(paste0(name_list[1], "/outs"), slice = paste(name_list[1])) %>% RenameCells(add.cell.id = name_list[1]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_2 <- Load10X_Spatial(paste0(name_list[2], "/outs"), slice = paste(name_list[2])) %>% RenameCells(add.cell.id = name_list[2]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_3 <- Load10X_Spatial(paste0(name_list[3], "/outs"), slice = paste(name_list[3])) %>% RenameCells(add.cell.id = name_list[3]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_4 <- Load10X_Spatial(paste0(name_list[4], "/outs"), slice = paste(name_list[4])) %>% RenameCells(add.cell.id = name_list[4]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_5 <- Load10X_Spatial(paste0(name_list[5], "/outs"), slice = paste(name_list[5])) %>% RenameCells(add.cell.id = name_list[5]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_6 <- Load10X_Spatial(paste0(name_list[6], "/outs"), slice = paste(name_list[6])) %>% RenameCells(add.cell.id = name_list[6]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_7 <- Load10X_Spatial(paste0(name_list[7], "/outs"), slice = paste(name_list[7])) %>% RenameCells(add.cell.id = name_list[7]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_8 <- Load10X_Spatial(paste0(name_list[8], "/outs"), slice = paste(name_list[8])) %>% RenameCells(add.cell.id = name_list[8]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_9 <- Load10X_Spatial(paste0(name_list[9], "/outs"), slice = paste(name_list[9])) %>% RenameCells(add.cell.id = name_list[9]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
data_10 <- Load10X_Spatial(paste0(name_list[10], "/outs"), slice = paste(name_list[10])) %>% RenameCells(add.cell.id = name_list[10]) %>% subset(nFeature_Spatial > 300) %>% NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

data_1$orig.ident <- name_list[1]
data_2$orig.ident <- name_list[2]
data_3$orig.ident <- name_list[3]
data_4$orig.ident <- name_list[4]
data_5$orig.ident <- name_list[5]
data_6$orig.ident <- name_list[6]
data_7$orig.ident <- name_list[7]
data_8$orig.ident <- name_list[8]
data_9$orig.ident <- name_list[9]
data_10$orig.ident <- name_list[10]

object_list <- c(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10) 

anchors <- FindIntegrationAnchors(object.list = object_list, dims = 1:20)

merge <- IntegrateData(anchorset = anchors, dims = 1:20) %>% ScaleData() %>% RunPCA() %>% RunUMAP(reduction = "pca", dims = 1:30, umap.method = 'umap-learn') %>% FindNeighbors(reduction = "pca", dims = 1:30) %>% FindClusters(resolution = 0.20)
 
DimPlot(merge, reduction = "umap", label = TRUE, label.size = 15, group.by = "orig.ident", pt.size = 0.8)

marker0 <- FindMarkers(merge, ident.1 = "0", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
marker1 <- FindMarkers(merge, ident.1 = "1", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
marker2 <- FindMarkers(merge, ident.1 = "2", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
marker3 <- FindMarkers(merge, ident.1 = "3", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
marker4 <- FindMarkers(merge, ident.1 = "4", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
marker5 <- FindMarkers(merge, ident.1 = "5", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)
marker6 <- FindMarkers(merge, ident.1 = "6", only.pos = TRUE, logfc.threshold = 0.25) %>% dplyr::filter(p_val_adj < 0.05)

data <- GetAssayData(merge, assay = "Spatial", slot = "counts") %>% NormalizeData(normalization.method = "RC", scale.factor = 1e6,  verbose = FALSE) %>% as.matrix()

nExpression <- apply(data, 1, function(x){
    result <- sum(x > 0)
    return(result)
    }
)

data <- subset(data, subset = nExpression > 1600)

datExpr <- t(data)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=1.5,col="red");
abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=1.5 ,col="red")

net <- blockwiseModules(datExpr, power = 10, networkType="signed", TOMType = "signed", maxBlockSize = 30000,　deepSplit = 4, nPreclusteringCenters = 500,
                        minModuleSize = 15, reassignThreshold = 1e-6, detectCutHeight = 0.998, mergeCutHeight = 0.15, numericLabels = TRUE, 　pamStage = TRUE, pamRespectsDendro = TRUE,                        
                        saveTOMs = FALSE, verbose = 3)
                        
reference <- readRDS("reference.rds")

library(CARD)

sc_count <- GetAssayData(reference, assay = "RNA", slot = "counts")
DefaultAssay(merge) <- "integrated"
Idents(merge) <- "orig.ident"

n <- length(name_list)

for (i in 1:n){
    data <- subset(merge, idents = name_list[i])
    assign(paste0(name_list[i], "_sub"), data)
}

for (i in 1:n){
    data <- eval(parse(text = paste0(name_list[i], "_sub")))
    result <- GetAssayData(data, assay = "Spatial", slot = "counts")
    assign(paste0(name_list[i], "_spatial_counts"), result)
}


spatial_location_WT_Sham <- merge@images$WT_Sham@coordinates[,2:3]
colnames(spatial_location_WT_Sham)[1] <- "x"
colnames(spatial_location_WT_Sham)[2] <- "y"
head(spatial_location_WT_Sham)

spatial_location_WT_MI_day1_1 <- merge@images$WT_MI_day1_1@coordinates[,2:3]
colnames(spatial_location_WT_MI_day1_1)[1] <- "x"
colnames(spatial_location_WT_MI_day1_1)[2] <- "y"
head(spatial_location_WT_MI_day1_1)

spatial_location_WT_MI_day1_2 <- merge@images$WT_MI_day1_2@coordinates[,2:3]
colnames(spatial_location_WT_MI_day1_2)[1] <- "x"
colnames(spatial_location_WT_MI_day1_2)[2] <- "y"
head(spatial_location_WT_MI_day1_2)

spatial_location_WT_MI_day1_3 <- merge@images$WT_MI_day1_3@coordinates[,2:3]
colnames(spatial_location_WT_MI_day1_3)[1] <- "x"
colnames(spatial_location_WT_MI_day1_3)[2] <- "y"
head(spatial_location_WT_MI_day1_3)

spatial_location_WT_MI_day7_1 <- merge@images$WT_MI_day7_1@coordinates[,2:3]
colnames(spatial_location_WT_MI_day7_1)[1] <- "x"
colnames(spatial_location_WT_MI_day7_1)[2] <- "y"
head(spatial_location_WT_MI_day7_1)

spatial_location_WT_MI_day7_2 <- merge@images$WT_MI_day7_2@coordinates[,2:3]
colnames(spatial_location_WT_MI_day7_2)[1] <- "x"
colnames(spatial_location_WT_MI_day7_2)[2] <- "y"
head(spatial_location_WT_MI_day7_2)

spatial_location_WT_MI_day7_3 <- merge@images$WT_MI_day7_3@coordinates[,2:3]
colnames(spatial_location_WT_MI_day7_3)[1] <- "x"
colnames(spatial_location_WT_MI_day7_3)[2] <- "y"
head(spatial_location_WT_MI_day7_3)

spatial_location_WT_MI_day14_1 <- merge@images$WT_MI_day14_1@coordinates[,2:3]
colnames(spatial_location_WT_MI_day14_1)[1] <- "x"
colnames(spatial_location_WT_MI_day14_1)[2] <- "y"
head(spatial_location_WT_MI_day14_1)

spatial_location_WT_MI_day14_2 <- merge@images$WT_MI_day14_2@coordinates[,2:3]
colnames(spatial_location_WT_MI_day14_2)[1] <- "x"
colnames(spatial_location_WT_MI_day14_2)[2] <- "y"
head(spatial_location_WT_MI_day14_2)

spatial_location_WT_MI_day14_3 <- merge@images$WT_MI_day14_3@coordinates[,2:3]
colnames(spatial_location_WT_MI_day14_3)[1] <- "x"
colnames(spatial_location_WT_MI_day14_3)[2] <- "y"
head(spatial_location_WT_MI_day14_3)

for (i in 1:n){
    spatial_count <- eval(parse(text = paste0(name_list[i], "_spatial_counts")))
    spatial_location <- eval(parse(text = paste0("spatial_location_", name_list[i])))
    result <- createCARDObject(
	  sc_count = sc_count,
	  sc_meta = reference@meta.data,
	  spatial_count = spatial_count,
	  spatial_location = spatial_location,
	  ct.varname = "celltype",
    ct.select = NULL,
	  sample.varname = "orig.ident",
	  minCountGene = 100,
	  minCountSpot = 5) 
    assign(paste0("CARD_obj_", name_list[i]), result)
}

for (i in 1:n){
    data <- eval(parse(text = paste0("CARD_obj_", name_list[i])))
    result <- CARD_deconvolution(CARD_object = data)
    assign(paste0("CARD_obj_", name_list[i]), result)
}

head(CARD_obj_WT_Sham@Proportion_CARD)
