library(Seurat)
library(plyr)
library(ggplot2)
library(DropletUtils)
library(celda)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scDblFinder)
library(viridis)
library(MASS)
library(patchwork)
library(harmony)
library(readr)
library(dplyr)
library(clustree)
library(tidyverse)
library(batchelor)

#####################################################################
cosmx_pre <- readRDS('sce.seurat.RData')
View(cosmx_pre@meta.data)

#Harmony Batch Correction
cosmx_pre <- NormalizeData(cosmx_pre)
cosmx_pre <- FindVariableFeatures(cosmx_pre)
cosmx_pre <- ScaleData(cosmx_pre, split.by='sample_id')
cosmx_pre <- RunPCA(object = cosmx_pre,npcs=50)
ElbowPlot(cosmx_pre)


cosmx2 <- cosmx_pre %>%
  RunHarmony(group.by.vars = 'sample_id', plot_convergence = FALSE)

cosmx2.embed <- Embeddings(cosmx2, 'harmony')

cosmx2 <- cosmx2 %>%
  RunUMAP(reduction = 'harmony', dims=1:22) %>%
  FindNeighbors(reduction='harmony', dims = 1:22) %>%
  FindClusters(resolution = 0.95) 

p1 <- DimPlot(cosmx2, reduction = 'umap', group.by = 'sample_id',raster=FALSE)
p2 <- DimPlot(cosmx2, reduction = 'umap', group.by = 'seurat_clusters',label=TRUE,raster=FALSE)

p1|p2

#Dot plot to identify major cell classes
cell_features <- c('EPCAM','AQP8','DUOXA1','LCN2','PIGR','MKI67','PCNA','TOP2A',
                   'BANK1','CD19','DERL3','MS4A1','MZB1',
                   'CD3D','CD3E','CD3G','CD8A','CD4',
                   'CD68','C1QA','ITGAM','CD14','ITGAX','FCGR3B','LYZ','CD274','HLA-DRA',
                   'ACTA2','ADAMDEC1','COL3A1','VWF','TPSAB1','TPSB2')
DotPlot(cosmx2, features = cell_features, dot.scale = 8) + RotatedAxis()

cell.markers <- FindAllMarkers(cosmx2, only.pos = TRUE)
top10<- cell.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
View(top10)

cell.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top20
DotPlot(cosmx2, features = top20$gene,dot.scale = 8) + RotatedAxis()
DoHeatmap(cosmx2, features = top20$gene) + NoLegend()

#Annotating major cell classes
new.cluster.ids <- c('Epithelial cells','T-cells','Epithelial cells Proliferative','Plasma cells','Plasma cells','Fibroblasts',
                    'Smooth Muscles','Epithelial cells','Fibroblasts','Unknown','Myeloid cells','Endothelial cells','Myeloid cells',
                    'Epithelial cells AQP8+','B-cells','Epithelial cells AQP8+','Mast cells','Myeloid cells','T-cells',
                    'Epithelial cells HLADRA+','Epithelial cells AQP8+')
names(new.cluster.ids) <- levels(cosmx2)
cosmx2 <- RenameIdents(cosmx2, new.cluster.ids)
DimPlot(cosmx2, reduction = "umap", label = TRUE, pt.size = 0.5,raster=FALSE)


unique(cosmx2@meta.data$sample_id)
cosmx2@meta.data<-cosmx2@meta.data %>%
  mutate(Status= case_when(
    sample_id =='HC_a' ~ 'HC',
    sample_id=='HC_b' ~'HC',
    sample_id=='HC_c' ~'HC',
    sample_id=='CD_a' ~'CD',
    sample_id=='CD_b' ~'CD',
    sample_id=='CD_c' ~'CD',
    sample_id=='UC_a' ~'UC',
    sample_id=='UC_b' ~'UC',
    sample_id=='UC_c' ~'UC'))
View(cosmx2@meta.data)

saveRDS(cosmx2, 'CosMx_Salas_pre_processed_harmony_corrected_major_cell_classified.rds')
######################################################################################################

