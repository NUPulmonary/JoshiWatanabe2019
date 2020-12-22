#A spatially restricted fibrotic niche in pulmonary fibrosis is sustained by 
#M-CSF/M-CSFR signaling in monocyte-derived alveolar macrophages

#Authors: Nikita Joshi, Satoshi Watanabe, Nikolay S Markov, 
#Rohan Verma, GR Scott Budinger, Alexander V Misharin
#2018-2019

library(Seurat) #Seurat V3.1.0
library(ggplot2)
library(dplyr)
library(sctransform)

#SC15 - asbestos
#SC14 - TiO2

#Combine two libraries: SC14 - TiO2 and SC15 - asbestos
SC15.data <- Read10X(data.dir = "/SC_15/outs/filtered_gene_bc_matrices/mm10/")
SC15 <- CreateSeuratObject(counts = SC15.data, project = "SC15", min.cells = 3, min.features = 200)
SC15 <- RenameCells(object = SC15, add.cell.id = "SC15")
SC15$stim <- "Asbestos"
head(x = colnames(x = SC15))

SC14.data <- Read10X(data.dir = "/SC_14/outs/filtered_gene_bc_matrices/mm10/")
SC14 <- CreateSeuratObject(counts = SC14.data, project = "SC14", min.cells = 3, min.features = 200)
SC14 <- RenameCells(object = SC14, add.cell.id = "SC14")
head(x = colnames(x = SC14))
SC14$stim <- "TiO2"

SC15
#An object of class Seurat 
#17224 features across 7173 samples within 1 assay 
#Active assay: RNA (17224 features)
SC14
#An object of class Seurat 
#16714 features across 8279 samples within 1 assay 
#Active assay: RNA (16714 features)

#Merge two datasets
lung <- merge(SC15, y = c(SC14), project = "asbestos")
head(x=colnames(x=lung))
tail(x=colnames(x=lung))

lung
#An object of class Seurat 
#17560 features across 15452 samples within 1 assay 
#Active assay: RNA (17560 features)

lung <- PercentageFeatureSet(object = lung, pattern = "^mt-", col.name = "percent.mt")
VlnPlot(object = lung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
FeatureScatter(object = lung, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
FeatureScatter(object = lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

#run sctransform
lung <- SCTransform(object = lung, verbose = T)
lung <- RunPCA(object = lung, npcs = 100, verbose = T)
ElbowPlot(object = lung, ndims = 100)
DimPlot(object = lung, reduction = "pca", group.by = "orig.ident")
DimHeatmap(object = lung, dims = 1:9, cells = 500, balanced = TRUE)
DimHeatmap(object = lung, dims = 10:18, cells = 500, balanced = TRUE)
DimHeatmap(object = lung, dims = 19:27, cells = 500, balanced = TRUE)
DimHeatmap(object = lung, dims = 28:36, cells = 500, balanced = TRUE)
DimHeatmap(object = lung, dims = 37:45, cells = 500, balanced = TRUE)
DimHeatmap(object = lung, dims = 46:54, cells = 500, balanced = TRUE)
DimHeatmap(object = lung, dims = 55:63, cells = 500, balanced = TRUE)

lung <- RunUMAP(object = lung, dims = 1:53, verbose = T)
lung <- FindNeighbors(object = lung, dims = 1:53, verbose = FALSE)
lung <- FindClusters(object = lung, verbose = FALSE, resolution = 1.0)
plot1<-DimPlot(object = lung, label = TRUE) + NoLegend()
plot2<-DimPlot(object = lung, label = TRUE, group.by = "orig.ident")
CombinePlots(plots = list(plot1, plot2))
ggplot(lung@meta.data, aes(x=SCT_snn_res.1, fill=orig.ident)) + geom_bar(position = "fill")

DefaultAssay(lung)
DefaultAssay(lung)<-"RNA"
lung<-NormalizeData(lung)
lung<-ScaleData(lung, verbose = FALSE)
lung.markers <- FindAllMarkers(object = lung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(lung.markers,"lung_PC1_53_res1.0.csv")

#save(lung, file = "lung_SeuratV3.Robj")
#load("lung_SeuratV3.Robj")

#C0 - Classical monocytes
#C1 - AT2
#C2 - B cells
#C3 - T cells
#C4 - AT2
#C5 - Non-classical monocytes
#C6 - Alveolar macrophages
#C7 - Neutrophils
#C8 - Alveolar macrophages
#C9 - Interstitial macrophages
#C10 - NK cells
#C11 - DC2
#C12 - AT2
#C13 - Endothelial-1
#C14 - Secretory (club and goblet)
#C15 - DC1
#C16 - Endothelial-2 (alveolar capillaries )
#C17 - T cells
#C18 - Neutrophils
#C19 - Alveolar macrophages (MoAM, Mmp12, Arg1)
#C20 - Ciliated
#C21 - AT2 cells (asbestos only, Retnla+)
#C22 - Neutrophils
#C23 - ILC (Id2, Gata3, Il7r)
#C24 - AM
#C25 - Proliferating AMs
#C26 - pDC
#C27 - Regulatory T cells
#C28 - Fibroblasts
#C29 - Ccr7+Ccl22+ DC
#C30 - Doublets: AT2 and neutrophils
#C31 - Mesothelium
#C32 - AT1 cells
#C33 - Proliferating DCs
#C34 - Doublets: monocytes and B cells
#C35 - Smooth muscle cells
#C36 - Lymphatics
#C37 - Doublets: AT2 and NK cells
#C38 - Mast cells

#Remove doublets
lung <- subset(lung, idents = c(30,34,37), invert = T)
lung
#An object of class Seurat 
#34455 features across 15288 samples within 2 assays 
#Active assay: SCT (16895 features)
#1 other assay present: RNA
#2 dimensional reductions calculated: pca, umap

#Rename clusters
lung<-RenameIdents(lung, '32'='AT1','1'='AT2', '4'='AT2', '12'='AT2', '21'='AT2_Retnla', '14'='Club_cells', '20'='Ciliated_cells', 
                   '28'='Fibroblasts','35'='Smooth_muscle_cells','36'='Lymphatics', '13'='Endo_1', '16'='Endo_2','31'='Mesothelium',
                   '0'='Classical_monocytes', '5'='Non-classical_monocytes', '7'='Neutrophils','18'='Neutrophils', '22'='Neutrophils', '38'='Mast_cells',
                   '6'='TRAM', '8'='TRAM', '24'='TRAM', '25'='Proliferating_TRAM', '19'='MoAM', '9'='IM', 
                   '15'='DC1', '11'='DC2','26'='pDC',  '29'="Ccr7_DC", '33'='Proliferating_DC',
                   '2'='B_cells', '10'='NK_cells','23'='ILC', '3'='T_cells', '17'='T_cells', '27'='Tregs')

lung.markers <- FindAllMarkers(object = lung, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.7)
write.csv(lung.markers,"Supplemental_Table_S1_all_clusters.csv") #Supplemental_Table_S1_all_clusters.csv
DimPlot(object = lung, label = TRUE) + NoLegend() #Figure 4A
DimPlot(object = lung, label = TRUE, group.by = "orig.ident", cols = c("#00BFC4","#F8766D"), pt.size = 0.1) + NoLegend() #Figure S4A
ggplot(FetchData(lung, c("ident", 'orig.ident')), aes(x=ident, fill=orig.ident)) + geom_bar(position = "fill") + coord_flip()+scale_fill_manual(values=c("#00BFC4","#F8766D")) #Figure S4B
FeaturePlot(object = lung, features = c("Mrc1"), min.cutoff = "q10", max.cutoff = "q90", order = T, pt.size = 0.1) + NoLegend() #Figure 4B
save(lung, file = "lung_clean.Robj")

#Macrophages----
#Explore macrophage subsets using integration workflow in Seurat V3 (https://satijalab.org/seurat/v3.0/pancreas_integration_label_transfer.html)
Macrophages <- subset(lung, idents = c("TRAM","IM","MoAM"))
DimPlot(object = Macrophages, label = TRUE) + NoLegend()
Macrophages.list <- SplitObject(Macrophages, split.by = "orig.ident")

for (i in 1:length(Macrophages.list)) {
  Macrophages.list[[i]] <- NormalizeData(Macrophages.list[[i]], verbose = FALSE)
  Macrophages.list[[i]] <- FindVariableFeatures(Macrophages.list[[i]], selection.method = "vst", 
                                                nfeatures = 2000, verbose = FALSE)
}
reference.list <- Macrophages.list[c("SC14", "SC15")]
Macrophages.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
Macrophages.integrated <- IntegrateData(anchorset = Macrophages.anchors, dims = 1:30)
DefaultAssay(Macrophages.integrated)
DefaultAssay(Macrophages.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
Macrophages.integrated <- ScaleData(Macrophages.integrated, verbose = FALSE)
Macrophages.integrated <- RunPCA(Macrophages.integrated, npcs = 30, verbose = FALSE)
Macrophages.integrated <- RunUMAP(Macrophages.integrated, reduction = "pca", dims = 1:30)
Macrophages.integrated <- FindNeighbors(object = Macrophages.integrated, dims = 1:30, verbose = FALSE)
Macrophages.integrated <- FindClusters(object = Macrophages.integrated, verbose = FALSE, resolution = 0.7)
DimPlot(Macrophages.integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(Macrophages.integrated, reduction = "umap", label = T)
ggplot(Macrophages.integrated@meta.data, aes(x=integrated_snn_res.0.7, fill=orig.ident)) + geom_bar(position = "fill")
FeaturePlot(object = Macrophages.integrated, features = c("Wfdc21"), order = T, min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(object = Macrophages.integrated, features = c("Csf1"), order = T, min.cutoff = "q10", max.cutoff = "q90") + NoLegend()

DefaultAssay(Macrophages.integrated)
DefaultAssay(Macrophages.integrated) <-"RNA"
Macrophages.integrated.markers <- FindAllMarkers(object = Macrophages.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Macrophages.integrated.markers,"Macrophages.integrated.markers.csv")
save(Macrophages.integrated, file = "Macropahges_integrated.Robj")
load("Macropahges_integrated.Robj")

#C0 - TRAM
#C1 - TRAM Ear1
#C2 - Peribronchial IM
#C3 - MoAM
#C4 - TRAM
#C5 - TRAM
#C6 - Perivascular IM
#C7 - B cells, etc
#C8 - proliferating
#C9 - Neutrophils
#C10 - NK cells
#C11 - Monocytes
#C12 - Epithelial cells

#Let's remove contaminating cells
Macrophages2 <- subset(Macrophages.integrated, idents = c(7,8,9,10,11,12), invert = T)
DimPlot(Macrophages2, reduction = "umap", group.by = "orig.ident")
DimPlot(Macrophages2, reduction = "umap", label = T)

#C0 - TRAM
#C1 - TRAM
#C2 - Peribronchial IM
#C3 - MoAM
#C4 - TRAM
#C5 - TRAM
#C6 - Perivascular IM

DefaultAssay(Macrophages2) <-"RNA"
Macrophages2<-NormalizeData(Macrophages2)
Macrophages2<-ScaleData(Macrophages2, verbose = FALSE)
Macrophages2.markers <- FindAllMarkers(object = Macrophages2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Macrophages2.markers,"Macrophages2.markers.csv")
#Cluster 0 has almost no genes differentiating it from cluster 1. Clusters 4 and 5 also look quite similar. Let's test this by running cluster tree. 
Macrophages2 <- BuildClusterTree(Macrophages2, dims=1:30)
PlotClusterTree(Macrophages2)

#Let's merge clusters 0 and 1, and clusters 4 and 5. 
Macrophages2<-RenameIdents(Macrophages2, '0'='AM1', '1'='AM1', '4'='AM2', '5'='AM2', '3'='AM3', '2'='IM1', '6'='IM2')
levels(Macrophages2)
p1<-DimPlot(Macrophages2, reduction = "umap", label = T, pt.size = 0.5) #Figure 4C
p2<-ggplot(FetchData(Macrophages2, c("ident", 'stim')), aes(x=ident, fill=stim)) + geom_bar(position = "fill") + scale_fill_manual(values=c("#F8766D","#00BFC4"))#Figure 4D
p3<-DimPlot(Macrophages2, reduction = "umap", group.by = "stim", pt.size = 0.1) #Figure 4E
CombinePlots(plots = list(p1, p4, p3), ncol = 1) #Figure 4C-E
save(Macrophages2, file= "Macropahges2.Robj")
load("Macropahges2.Robj") 

DimPlot(Macrophages2, reduction = "umap", label = T)
plot <- DimPlot(object = Macrophages2) + NoLegend()
select.cells <- CellSelector(plot = plot)
Macrophages2 <- subset(Macrophages2, cells = select.cells)
DimPlot(Macrophages2, reduction = "umap", label = T)
save(Macrophages2, file= "Macrophages_clean.Robj")

plot1<-FeaturePlot(object = Macrophages2, features = c("Cd68"), min.cutoff = "q10", max.cutoff = "q90", order = T) + NoLegend()
plot2<-FeaturePlot(object = Macrophages2, features = c("Car4"), min.cutoff = "q10", max.cutoff = "q90", order = T) + NoLegend()
plot3<-FeaturePlot(object = Macrophages2, features = c("Mmp12"), min.cutoff = "q10", max.cutoff = "q90", order = T) + NoLegend()
plot4<-FeaturePlot(object = Macrophages2, features = c("Cx3cr1"), min.cutoff = "q10", max.cutoff = "q90", order = T) + NoLegend()
plot5<-FeaturePlot(object = Macrophages2, features = c("Lyve1"), min.cutoff = "q10", max.cutoff = "q90", order = T) + NoLegend()
plot6<-FeaturePlot(object = Macrophages2, features = c("Ccr2"), min.cutoff = "q10", max.cutoff = "q90", order = T) + NoLegend()
CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5, plot6), ncol = 3) #Figure 4F

DotPlot(Macrophages2, features = c("Atf3", "Klf9", "Cebpb", 
                                   "Bhlhe40", "Bhlhe41", "Atf4", 
                                   "Jund", "Litaf"), assay = "RNA", plot.legend = T, dot.min = 0.1)+ coord_flip() #Figure 4G

DotPlot(Macrophages2, features = c("F13a1","Lyve1","Ccr2","C1qa","Cx3cr1","Gpnmb","Itgam","Cd36",
                                   "Itgax","Ear1","Car4","Pparg","Marco","Siglecf","Il18",
                                   "Axl","Mertk","Adgre1","Lyz2","Cd68",
                                   "Mrc1"), assay = "RNA", plot.legend = T, dot.min = 0.1) + coord_flip()#Figure S4C

Macrophages2.markers <- FindAllMarkers(object = Macrophages2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Macrophages2.markers,"Supplemental_Table_S2_Asbestos_Macrophages.csv")
top10 <- Macrophages2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Macrophages2, features = top10$gene, assay = "RNA", angle = 0, slot = "scale.data", raster = F) + NoLegend() #Figure S4D

Macrophages3<-subset(Macrophages2, idents = c("AM1", "AM2", "AM3"))
HarmonizomePFlist <- (read.table("HarmonizomePFlist.txt"))
HarmonizomePFlist=as.vector(HarmonizomePFlist[,1]) #convert to vector
M3.average<-AverageExpression(Macrophages3, return.seurat = T)
Macrophages3.markers <- FindAllMarkers(object = Macrophages3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
M3_Harmonizome <-intersect(Macrophages3.markers$gene, HarmonizomePFlist)
genes.to.plot <- Macrophages3.markers %>% group_by(pct.1) %>% filter(gene %in% c(M3_Harmonizome))
DoHeatmap(object = M3.average,features = genes.to.plot$gene, assay = "RNA", angle = 0, slot = "scale.data", draw.lines=F, raster = F) #Figure 4H

M3.average<-AverageExpression(Macrophages3, return.seurat = T, add.ident = "orig.ident")
Macrophages3.markers <- FindAllMarkers(object = Macrophages3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
M3_Harmonizome <-intersect(Macrophages3.markers$gene, HarmonizomePFlist)
genes.to.plot <- Macrophages3.markers %>% group_by(pct.1) %>% filter(gene %in% c(M3_Harmonizome))
DoHeatmap(object = M3.average,features = genes.to.plot$gene, assay = "RNA", angle = 0, slot = "scale.data", draw.lines=F, raster = F) #Figure 4H

M3.average<-AverageExpression(Macrophages3, return.seurat = T, add.ident = "orig.ident")
DoHeatmap(object = M3.average,features = genes.to.plot$gene, assay = "RNA", angle = 0, slot = "scale.data", draw.lines=F, raster = F) #Figure S4E

#Figure 5----
DotPlot(Macrophages3, features = c("Csf1", "Csf1r", "Csf2rb"), 
        assay = "RNA", dot.scale = 5, dot.min = 0.01, split.by = "stim", 
        cols = c("blue", "red")) #Figure 5A
FeaturePlot(object = Macrophages2, features = c("Csf1"), order = T, pt.size = 1.5) + NoLegend() #Figure 5B

DotPlot(lung, features = c("Csf1", "Il34"), 
        assay = "RNA", dot.scale = 6, dot.min = 0.05, split.by = "stim", 
        cols = c("red", "blue")) #Figure S5A

#Figure 6----
FeaturePlot(object = Macrophages2, features = c("Pdgfa"), order = T, pt.size = 1.5) + NoLegend() #Figure 5B

#AT2 cells, Figure 7 ----
AT2 <- subset(lung, idents = c("AT2", "AT2_Retnla"))
DimPlot(object = AT2, label = TRUE) + NoLegend()
DefaultAssay(AT2)

AT2.list <- SplitObject(AT2, split.by = "orig.ident")

for (i in 1:length(AT2.list)) {
  AT2.list[[i]] <- NormalizeData(AT2.list[[i]], verbose = FALSE)
  AT2.list[[i]] <- FindVariableFeatures(AT2.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000, verbose = FALSE)
}
reference.list <- AT2.list[c("SC14", "SC15")]
AT2.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
AT2.integrated <- IntegrateData(anchorset = AT2.anchors, dims = 1:30)
DefaultAssay(AT2.integrated)
DefaultAssay(AT2.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
AT2.integrated <- ScaleData(AT2.integrated,verbose = FALSE)
#AT2.integrated <- ScaleData(AT2.integrated, vars.to.regress = "percent.mt",verbose = FALSE)
AT2.integrated <- RunPCA(AT2.integrated, npcs = 30, verbose = FALSE)
AT2.integrated <- RunUMAP(AT2.integrated, reduction = "pca", dims = 1:30)
AT2.integrated <- FindNeighbors(object = AT2.integrated, dims = 1:30, verbose = FALSE)
AT2.integrated <- FindClusters(object = AT2.integrated, verbose = FALSE, resolution = 0.7)
DimPlot(AT2.integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(AT2.integrated, reduction = "umap", label = T)
ggplot(AT2.integrated@meta.data, aes(x=integrated_snn_res.0.7, fill=orig.ident)) + geom_bar(position = "fill")

DefaultAssay(AT2.integrated)
DefaultAssay(AT2.integrated) <-"RNA"
FeaturePlot(object = AT2.integrated, features = c("Retnla"), order = T, min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(object = AT2.integrated, features = c("Lgals3"), order = T, min.cutoff = "q10", max.cutoff = "q90") + NoLegend()

AT2.integrated<-NormalizeData(AT2.integrated)
AT2.integrated<-ScaleData(AT2.integrated, verbose = FALSE)
AT2.integrated.markers <- FindAllMarkers(object = AT2.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(AT2.integrated.markers,"AT2.integrated.markers.csv")
save(AT2.integrated, file = "AT2_integrated.Robj")
load("AT2_integrated.Robj")
VlnPlot(AT2.integrated, features = "Retnla", split.by = "orig.ident", pt.size = 0)
VlnPlot(AT2.integrated, features = "nFeature_RNA", pt.size = 0)

#C0 - AT2
#C1 - AT2
#C2 - AT2
#C3 - AT2
#C4 - AT2 Retnla+ Slc26a4+
#C5 - AT2 Ndnf+ Lgals3+ Krt18+ Krt8+
#C6 - B cells
#C7 - Monocytes
#C8 - AT2 Mitochondrial genes high, low UMI count
#C9 - T cells
#C10 - Monocytes/Macrophages
#C11 - Alveolar macrophages
#C12 - Endothelial
#C13 - misc

#Let's remove contaminating cells
AT2_clean <- subset(AT2.integrated, idents = c(6,7,8,9,10,11,12,13), invert = T)
DimPlot(AT2_clean, reduction = "umap", group.by = "orig.ident")
DimPlot(AT2_clean, reduction = "umap", label = T)

plot <- DimPlot(object = AT2_clean) + NoLegend()
select.cells <- CellSelector(plot = plot)
AT2_clean <- subset(AT2.integrated, cells = select.cells)
DimPlot(AT2_clean, reduction = "umap", label = T)

#Clusters 0 has almost no genes differentiating it from cluster 2. Clusters 1 and 3 also look quite similar. Let's test this by running cluster tree. 
AT2_clean <- BuildClusterTree(AT2_clean, dims=1:30)
PlotClusterTree(AT2_clean)

#Let's merge clusters 0 and 1, and clusters 4 and 5. 
AT2_clean<-RenameIdents(AT2_clean, '0'='AT2_1', '2'='AT2_1', '1'='AT2_2', '3'='AT2_2', '5'='AT2_3', '4'='AT2_4')
levels(AT2_clean)
p1<-DimPlot(AT2_clean, reduction = "umap", label = T, pt.size = 0.5) #Figure 7A
p2<-DimPlot(AT2_clean, reduction = "umap", group.by = "stim", pt.size = 0.5) #Figure 7B
CombinePlots(plots = list(p1, p2), ncol = 2)
ggplot(FetchData(AT2_clean, c("ident", 'stim')), aes(x=ident, fill=stim)) + geom_bar(position = "fill") + scale_fill_manual(values=c("#F8766D","#00BFC4"))#Figure 7C
AT2_clean.markers <- FindAllMarkers(object = AT2_clean, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25)
write.csv(AT2_clean.markers,"Supplemental_Table_S7_AT2_cells.csv")
save(AT2_clean, file= "AT2_clean.Robj")

FeaturePlot(object = AT2_clean, features = c("Retnla"), order = T, min.cutoff = "q10", max.cutoff = "q90", pt.size = 1.0) + NoLegend() #Figure 7D
FeaturePlot(object = AT2_clean, features = c("Lgals3"), order = T, min.cutoff = "q10", max.cutoff = "q90") + NoLegend()
FeaturePlot(object = AT2_clean, features = c("Krt18"), order = T, min.cutoff = "q10", max.cutoff = "q90") + NoLegend()

top10 <- AT2_clean.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(AT2_clean, features = top10$gene, assay = "RNA", angle = 0, slot = "scale.data", raster=F) + NoLegend() #Figure S7

#The end----