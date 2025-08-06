#Orthotopics TCDD and Vehicle for Brian
#February 2024
#Seurat Version 4.3.0.1
#SeuratObject Version 4.1.3
#R version 4.3.0 (2023-04-21) -- "Already Tomorrow"

library(Seurat)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(scCustomize)
library(pheatmap)
library(msigdbr)
library(fgsea)
library(viridis)

#Sample Info
#BG-1:Vehicle 1
#BG-2:Vehicle 2
#BG-3:TCDD 1
#BG-4:TCDD 2
#Orthotopic Tumor Cell line-7940B

#### FULL OBJECT PREPROCESSING ####

#Load Dataset, make sure that genes are renamed to features
Vehicle1.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-1/") 
Vehicle2.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-2/") 
TCDD1.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-3/") 
TCDD2.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-4/") 

#Create Seurat objects
Vehicle1  <- CreateSeuratObject(counts = Vehicle1.data, project = 'Vehicle1', min.cells = 3, min.features = 100)
Vehicle2 <- CreateSeuratObject(counts = Vehicle2.data, project = 'Vehicle2', min.cells = 3, min.features = 100)
TCDD1 <- CreateSeuratObject(counts = TCDD1.data, project = 'TCDD1', min.cells = 3, min.features = 100)
TCDD2 <- CreateSeuratObject(counts = TCDD2.data, project = 'TCDD2', min.cells = 3, min.features = 100)

#Group MetaData
Vehicle1$Group <-"Vehicle"
Vehicle2$Group <-"Vehicle"
TCDD1$Group<-"TCDD"
TCDD2$Group <-"TCDD"

#Sample MetaData
Vehicle1$Sample <-"Vehicle 1"
Vehicle2$Sample <-"Vehicle 2"
TCDD1$Sample<-"TCDD 1"
TCDD2$Sample <-"TCDD 2"

#Merge samples
OrthotopicTCDD <- merge(Vehicle1, y =c(Vehicle2, TCDD1, TCDD2), 
                             add.cell.ids = c("Vehicle1", "Vehicle2","TCDD1","TCDD2"))

OrthotopicTCDD <- NormalizeData(OrthotopicTCDD)
OrthotopicTCDD[["percent.mt"]] <- PercentageFeatureSet(OrthotopicTCDD, pattern = "^mt-")
OrthotopicTCDD <- subset(OrthotopicTCDD, subset = nCount_RNA > 1200 & nCount_RNA < 100000 & percent.mt < 15)

Idents(OrthotopicTCDD) <- "Group"
table(OrthotopicTCDD@active.ident)
#Vehicle    TCDD 
#  20995   20922 

#Normalization ScaleData and UMAP Generation
OrthotopicTCDD <- FindVariableFeatures(OrthotopicTCDD, selection.method = "vst", nfeatures = 2000)
OrthotopicTCDD <- ScaleData(OrthotopicTCDD, verbose = T, features = row.names(OrthotopicTCDD))
OrthotopicTCDD <- RunPCA(object = OrthotopicTCDD)
stdev <- OrthotopicTCDD@reductions$pca@stdev
var <- stdev^2
sum(var[1:29])/ sum(var) 

PCNum = 29

OrthotopicTCDD <- FindNeighbors(object = OrthotopicTCDD, dims = 1:PCNum)
OrthotopicTCDD <- FindClusters(object = OrthotopicTCDD, resolution = 3)
OrthotopicTCDD<- RunUMAP(object = OrthotopicTCDD, dims = 1:PCNum)
DimPlot(object = OrthotopicTCDD, reduction = "umap", label = T)

DotPlot(OrthotopicTCDD, features = c("Krt19", "Cdh1","Epcam","Msln","Try4", "Amy2a1", "Col1a2", "Acta2", "Clec3b", "Cspg4","Saa3","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                              "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Label Clusters
Idents(OrthotopicTCDD) <- "seurat_clusters"
OrthotopicTCDD <- RenameIdents(OrthotopicTCDD,
                          "0" = "Fibroblast", 
                          "1" = "Ductal", 
                          "2" = "Ductal - Proliferating", 
                          "3" = "Macrophage", 
                          "4" = "Macrophage", 
                          "5" = "Fibroblast", 
                          "6" = "Macrophage - Proliferating", 
                          "7" = "Ductal", 
                          "8" = "Macrophage", 
                          "9" = "Granulocyte", 
                          "10" = "Ductal", 
                          "11" = "Ductal",
                          "12" = "Granulocyte", 
                          "13" = "Macrophage", 
                          "14" = "Macrophage", 
                          "15" = "Macrophage", 
                          "16" = "Fibroblast", 
                          "17" = "Endothelial", 
                          "18" = "Macrophage", 
                          "19" = "Macrophage", 
                          "20" = "Macrophage", 
                          "21" = "Ductal", 
                          "22" = "T Cell", 
                          "23" = "Macrophage",
                          "24" = "Macrophage", 
                          "25" = "Fibroblast - Proliferating", 
                          "26" = "Fibroblast", 
                          "27" = "Fibroblast", 
                          "28" = "Macrophage", 
                          "29" = "B Cell", 
                          "30" = "T Cell", 
                          "31" = "Granulocyte", 
                          "32" = "Ductal - Proliferating", 
                          "33" = "Dendritic Cell", 
                          "34" = "Pericyte",
                          "35" = "Ductal", 
                          "36" = "RBC", 
                          "37" = "Macrophage",
                          "38" = "Macrophage", 
                          "39" = "Acinar", 
                          "40" = "NK Cell",
                          "41" = "Granulocyte", 
                          "42" = "T Cell - Proliferating", 
                          "43" = "Dendritic Cell",
                          "44" = "T Cell", 
                          "45" = "Fibroblast",
                          "46" = "Endothelial", 
                          "47" = "Macrophage",
                          "48" = "RBC", 
                          "49" = "Macrophage", 
                          "50" = "Fibroblast - Proliferating",
                          "51" = "Ductal", 
                          "52" = "Ductal", 
                          "53" = "Macrophage",
                          "54" = "RBC", 
                          "55" = "RBC - Proliferating",
                          "56" = "Macrophage", 
                          "57" = "Dendritic Cell",
                          "58" = "Endothelial", 
                          "59" = "Dendritic Cell - Proliferating", 
                          "60" = "Mast Cell",
                          "61" = "Endothelial")
OrthotopicTCDD[["global_clusters"]] <- OrthotopicTCDD@active.ident

Idents(OrthotopicTCDD) <- "Group"
new_order <- c("Vehicle", "TCDD")
OrthotopicTCDD@active.ident <- factor(OrthotopicTCDD@active.ident, levels = new_order)
OrthotopicTCDD[["Group"]] <- OrthotopicTCDD@active.ident 

new_order <- c("Ductal",
               "Acinar",
               "Endothelial",
               "Fibroblast",
               "Pericyte",
               "Macrophage",
               "Dendritic Cell",
               "Granulocyte",
               "T Cell",
               "NK Cell",
               "Mast Cell",
               "B Cell",
               "RBC",
               "Ductal - Proliferating",
               "Fibroblast - Proliferating",
               "T Cell - Proliferating",
               "Macrophage - Proliferating",
               "Dendritic Cell - Proliferating",
               "RBC - Proliferating")
OrthotopicTCDD@active.ident <- factor(OrthotopicTCDD@active.ident, levels = new_order)
OrthotopicTCDD[["global_clusters"]] <- OrthotopicTCDD@active.ident

save(OrthotopicTCDD,file="OrthotopicTCDD_KLD_Feb_24.RData")

#### MAKE YOUR FULL OBJECT UMAP FOR FIGURE ####
DimPlot(OrthotopicTCDD, split.by = "Group", group.by = "global_clusters", cols = c()) #add colors

#### T CELL PROCESSING ####

#T cell analysis
Tcells <-subset(OrthotopicTCDD_KLD, idents =  c("T Cell", "NK Cell", "T Cell - Proliferating"))

Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tcells)
Tcells <- ScaleData(Tcells, verbose = T, features = all.genes)
Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE)

#Calculate 90% Variiance:
st_dev <- Tcells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:19])/ sum(var)

Tcells <- FindNeighbors(object = Tcells, dims = 1:19)
Tcells <- FindClusters(object = Tcells, resolution = 4.5)
Tcells <- RunUMAP(Tcells, reduction = "pca", dims = 1:19, verbose = F)
DimPlot(Tcells, reduction = "umap", label = T)

DotPlot(Tcells, features = c("Krt19", "Cdh1","Epcam","Try4", "Amy2a1", "Col1a2", "Pdpn", "Pdgfra", 
                                         "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3","Il2ra", "Cd8a", "Nkg7","Arg1", 
                                         "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33",
                                         "H2-Eb1",    "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", 
                                         "Ms4a1","Kit","Il1rl1","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

Tcells <- subset(Tcells, idents = c(0, 9, 11, 12, 15, 18, 21, 26, 27, 28, 30, 31, 32), invert = T) #remove NK and myeloid cells

Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Tcells)
Tcells <- ScaleData(Tcells, verbose = T, features = all.genes)
Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE)

#Calculate 90% Variiance:
st_dev <- Tcells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

Tcells <- FindNeighbors(object = Tcells, dims = 1:20)
Tcells <- FindClusters(object = Tcells, resolution = 4)
Tcells <- RunUMAP(Tcells, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(Tcells, reduction = "umap", label = T)
DimPlot(Tcells, reduction = "umap", label = T, split.by = "Group")

FeaturePlot(Tcells, features = c("Rorc", "Il17b", "Il17c", "Il17f", "Il13"), cols = c("gainsboro", "firebrick1"), order = T)

VlnPlot(Tcells, features = c("Cd3e","Cd4", "Cd8a", "Trdc"), ncol = 2)

DotPlot(Tcells, features = c("Krt19", "Cdh1","Epcam","Try4", "Amy2a1", "Col1a2", "Pdpn", "Pdgfra", 
                             "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3","Il2ra", "Cd8a", "Nkg7","Arg1", 
                             "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33",
                             "H2-Eb1",    "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", 
                             "Ms4a1","Kit","Il1rl1","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

DotPlot(Tcells, features = c("Cd3e","Cd4","Cd8a","Ccr7","Sell","Foxp3","Gata3","Rora","Ifng","Gzmb","Prf1","Pdcd1","Lag3","Tox"), cols = "RdBu", dot.scale = 8) + RotatedAxis()

DotPlot(Tcells, features = c("Cd3e", "Cd3d","Sell", "Ccr7", "Tcf7","Il7r", "Cd28", "Cd4", "Tbx21", "Ifng", "Cxcr3", "Arg1","Gata3", "Il1rl1", "Il4", "Il13", "Rorc", "Il17a", "Il22", "Il2ra", "Foxp3", "Icos", "Cd8a", "Gzmb","Prf1", "Eomes","Lag3","Pdcd1","Trdc","Nkg7", "Klrg1", "Ccna2", "Mki67", "Ly6c2"), cols = "RdYlBu")+ RotatedAxis()

cluster_7 <- FindMarkers(Tcells, ident.1 = "7", only.pos = T)
cluster_8 <- FindMarkers(Tcells, ident.1 = "8", only.pos = T)

Idents(Tcells) <- "seurat_clusters"
Tcells <- RenameIdents(Tcells,
                    "0" = "Naive CD4 & CD8", 
                    "1" = "Memory CD4", 
                    "2" = "Proliferating CD8", 
                    "3" = "CD8", 
                    "4" = "CD8", 
                    "5" = "Treg", 
                    "6" = "CD8", 
                    "7" = "Naive CD4 & CD8", 
                    "8" = "Naive CD4 & CD8",
                    "9" = "Treg", 
                    "10" = "CD8",
                    "11" = "CD8",
                    "12" = "CD8", 
                    "13" = "Proliferating CD8", 
                    "14" = "CD8", 
                    "15" = "Proliferating CD8", 
                    "16" = "Th1", 
                    "17" = "Treg",
                    "18" = "Th1",
                    "19" = "Memory CD4", 
                    "20" = "Treg",
                    "21" = "Proliferating CD8",
                    "22" = "Naive CD4 & CD8", 
                    "23" = "Proliferating CD4", 
                    "24" = "Naive CD4 & CD8")
Tcells[["t_clusters"]] <- Tcells@active.ident

new_order <- c("CD8",
               "Proliferating CD8",
               "Treg",
               "Th1",
               "Proliferating CD4",
               "Memory CD4",
               "Naive CD4 & CD8")
Tcells@active.ident <- factor(Tcells@active.ident, levels = new_order)
Tcells[["t_clusters"]] <- Tcells@active.ident

Idents(Tcells) <- "t_clusters"
Tcells <- RenameIdents(Tcells,
                       "Naive CD4 & CD8" = "Naive CD4 & CD8", 
                       "Memory CD4" = "Memory CD4", 
                       "Proliferating CD8" = "CD8", 
                       "CD8" = "CD8", 
                       "Treg" = "Treg", 
                       "Proliferating CD4" = "Treg", 
                       "Th1" = "Th1")
Tcells[["simple_t_clusters"]] <- Tcells@active.ident

new_order <- c("CD8",
               "Treg",
               "Th1",
               "Memory CD4",
               "Naive CD4 & CD8")
Tcells@active.ident <- factor(Tcells@active.ident, levels = new_order)
Tcells[["simple_t_clusters"]] <- Tcells@active.ident

save(Tcells,file="TCDD_Tcells_KLD_Feb_24.RData")

#### TOTAL T CELL FIGURES ####

plot_colors <- c("red",
               "yellow",
               "orange",
               "green",
               "blue")

#T cell UMAPs
DimPlot(Tcells, reduction = "umap", label = F, split.by = "Group", pt.size = 1.5, cols = plot_colors)

#T cell abundance bar graphs
Idents(object = Tcells) <- 'Group'
samples.list <- unique(Tcells$Group)
clusters <- lapply(samples.list, function(x){
  subset <- subset(Tcells, subset = Group == x)
  dist <- data.frame(table(subset$simple_t_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("simple_t_clusters","variable","value","Group")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("simple_t_clusters","variable","value","Group")

ggplot(clusters_percent_dist, aes(fill=simple_t_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Vehicle", "TCDD")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance") +scale_fill_manual(values = plot_colors)

#AhR gene signature in t cells
DotPlot(Tcells, features = c("Cyp1a1", "Cyp1b1", "Ahrr", "Tiparp"), split.by = "Group", group.by = "simple_t_clusters", cols = "RdYlBu", dot.scale = 12) + RotatedAxis()
DotPlot(OrthotopicTCDD, features = c("Cyp1a1", "Cyp1b1", "Ahrr", "Tiparp"), split.by = "Group", cols = "RdYlBu", dot.scale = 12) + RotatedAxis()

#checking the number of cells in each group
Idents(Tcells) <- "Group"
table(Tcells@active.ident)

#Number of T cells per group
#Vehicle    TCDD 
#  1282     440 

Idents(Tcells) <- "Group"
t_vehicle <- subset(Tcells, idents = "Vehicle")
t_tcdd <- subset(Tcells, idents = "TCDD")

Idents(t_vehicle) <- "simple_t_clusters"
Idents(t_tcdd) <- "simple_t_clusters"

table(t_vehicle@active.ident)
#      CD8            Treg             Th1      Memory CD4 Naive CD4 & CD8 
#.     732             180             101              97             172 

table(t_tcdd@active.ident)
#            CD8            Treg             Th1      Memory CD4 Naive CD4 & CD8 
#.            86             111               5              60             178 

#Gene set signature scoring (not using):
#Ahr_genes <- c("Cyp1a1", "Cyp1b1", "Ahrr", "Tiparp")
#Add module scores for each list:
#Tcells <- AddModuleScore(Tcells,features = list(Ahr_genes),
#                                         name="Ahr_genes_Signature_Score")

#FeaturePlot(Tcells,
#            features = "Ahr_genes_Signature_Score1", label = FALSE, repel = TRUE, pt.size = 2, order = TRUE, split.by = "Group", cols = c("gainsboro", "firebrick1"))

#VlnPlot(Tcells, features = "Ahr_genes_Signature_Score1", split.by = "Group", group.by = "simple_t_clusters")

#### CD8 PREPROCESSING ####
Idents(Tcells) <- "simple_t_clusters" 
cd8s <- subset(Tcells, idents = "CD8")
cd8s <- FindVariableFeatures(cd8s, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(cd8s)
cd8s <- ScaleData(cd8s, verbose = T, features = all.genes)
cd8s <- RunPCA(cd8s, npcs = 30, verbose = FALSE)

#Calculate 90% Variiance:
st_dev <- cd8s@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

cd8s <- FindNeighbors(object = cd8s, dims = 1:22)
cd8s <- FindClusters(object = cd8s, resolution = 4.5)
cd8s <- RunUMAP(cd8s, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(cd8s, reduction = "umap", label = T)

VlnPlot(cd8s, features = c("Cd3e","Cd4", "Cd8a", "Trdc"), ncol = 2)
DotPlot(cd8s, features = c("Cd3e", "Cd3d","Sell", "Ccr7", "Tcf7","Il7r", "Cd28", "Cd4", "Tbx21", "Ifng", "Cxcr3", "Arg1","Gata3", "Il1rl1", "Il4", "Il13", "Rorc", "Il17c",  "Il22", "Il2ra", "Foxp3", "Icos", "Cd8a", "Gzmb","Prf1", "Eomes","Lag3","Pdcd1","Trdc","Nkg7", "Klrg1", "Ccna2", "Mki67", "Ly6c2"), cols = "RdYlBu")+ RotatedAxis()

cd8s <- subset(cd8s, idents = c("15","16"), invert = T)

#### CD8 FIGURES ####

#CD8 DE (nothing significant)
VlnPlot(cd8s, features = c("Ifng", "Prf1", "Gzmb"), group.by = "Group")
cd8_de <- FindMarkers(cd8s, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group")
cd8_de[["genes"]] <- row.names(cd8_de)

#CD8 HEATMAP
cd8_genes <-  c("Gzmb", "Ifng", "Prf1","Tox", "Nkg7", "Lag3", "Eomes", "Timp3", "Pdcd1", "Ctla4", "Tigit","Mki67", "Ccna2")
cd8_data <- FetchData(cd8s, vars = c(cd8_genes, 'Sample', 'Group'))

cell_avg <- data.frame()
n <- 2 # number of metadata columns

for (id in levels(factor(cd8_data$Sample))) {
  data_subset <- cd8_data %>% filter(Sample == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-n)], 2, mean)
  cell_avg <- rbind(cell_avg, data_subset_avg)
}

colnames(cell_avg) <- colnames(cd8_data)[1:(ncol(cd8_data)-n)]
rownames(cell_avg) <- levels(factor(cd8_data$Sample))

metadata <- unique(cd8_data %>% select('Group', 'Sample'))
rownames(metadata) <- metadata$Sample
metadata <-metadata['Group']
my_colour = list(Group = c("Vehicle" = "#f0f0f0", "TCDD" = "#080808"))

pheatmap(as.matrix(cell_avg), fontsize = 14, cluster_cols = F, annotation_row = metadata, color = colorRampPalette(colors = c('#0000FF','#FFFFFF','#FF0000'))(250), scale = 'column', annotation_colors = my_colour)

#padma likes these colors LOL
pheatmap(as.matrix(cell_avg), fontsize = 14, cluster_cols = F, annotation_row = metadata, color = colorRampPalette(colors = c('#a1d76a','#f7f7f7','#e9a3c9'))(250), scale = 'column', annotation_colors = my_colour)


#CD8 Pathway analysis (nothing interesting)
cd8_DE <- FindMarkers(cd8s, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group", logfc.threshold = 0, min.pct = 0)

h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

ranks <- cd8_DE$avg_log2FC
names(ranks) <- rownames(cd8_DE)

fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

result <- apply(fgseaRes,2,as.character)

fgseaRes_sig <- subset(fgseaRes, fgseaRes$padj <= 0.05)

waterfall_plot_sig <- function (fgseaRes, graph_title) {
  (fgseaRes %>% 
     mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
     ggplot( aes(reorder(short_name,NES), NES)) +
     geom_bar(stat= "identity", aes(fill = padj))+
     coord_flip()+
     labs(x = "Hallmark Gene Set", y = "Normalized Enrichment Score", title = graph_title)+
     theme(axis.text.y = element_text(size = 7), 
           plot.title = element_text(hjust = 1)))+scale_fill_viridis(direction = 1, option = "A")
}
waterfall_plot_sig(fgseaRes_sig, "Pathways Enriched in WT vs. KO - cd8s")

plotEnrichment(msigdbr_list[["HALLMARK_E2F_TARGETS"]], ranks, ticksSize = 0.5) + labs(title = "Vehicle VS TCDD", subtitle="HALLMARK_E2F_TARGETS")

save(cd8s,file="TCDD_CD8_KLD_Feb_24.RData")

#### TREG PREPROCESSING ####
Idents(Tcells) <- "simple_t_clusters" 
tregs <- subset(Tcells, idents = "Treg")
tregs <- FindVariableFeatures(tregs, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tregs)
tregs <- ScaleData(tregs, verbose = T, features = all.genes)
tregs <- RunPCA(tregs, npcs = 30, verbose = FALSE)

#Calculate 90% Variiance:
st_dev <- tregs@reductions$pca@stdev
var <- st_dev^2
sum(var[1:24])/ sum(var)

tregs <- FindNeighbors(object = tregs, dims = 1:24)
tregs <- FindClusters(object = tregs, resolution = 4.5)
tregs <- RunUMAP(tregs, reduction = "pca", dims = 1:24, verbose = F)
DimPlot(tregs, reduction = "umap", label = T)

VlnPlot(tregs, features = "Foxp3")

save(tregs,file="TCDD_Tregs_KLD_Feb_24.RData")

#### TREG PATHWAY ANALYSIS ####
treg_DE <- FindMarkers(tregs, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group", logfc.threshold = 0, min.pct = 0)

h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

ranks <- treg_DE$avg_log2FC
names(ranks) <- rownames(treg_DE)

fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

result <- apply(fgseaRes,2,as.character)

plotEnrichment(msigdbr_list[["HALLMARK_G2M_CHECKPOINT"]], ranks, ticksSize = 0.5) + labs(title = "Vehicle VS TCDD", subtitle="HALLMARK_G2M_CHECKPOINT")
plotEnrichment(msigdbr_list[["HALLMARK_IL2_STAT5_SIGNALING"]], ranks, ticksSize = 0.5) + labs(title = "Vehicle VS TCDD", subtitle="HALLMARK_IL2_STAT5_SIGNALING")
plotEnrichment(msigdbr_list[["HALLMARK_E2F_TARGETS"]], ranks, ticksSize = 0.5) + labs(title = "Vehicle VS TCDD", subtitle="HALLMARK_E2F_TARGETS")
plotEnrichment(msigdbr_list[["HALLMARK_MITOTIC_SPINDLE"]], ranks, ticksSize = 0.5) + labs(title = "Vehicle VS TCDD", subtitle="HALLMARK_MITOTIC_SPINDLE")

#### FIBROBLASTS ####
Idents(OrthotopicTCDD_KLD) <-"global_clusters"
levels(OrthotopicTCDD_KLD)

Fibroblast<-subset(OrthotopicTCDD_KLD, idents =  c("Fibroblast", "Fibroblast - Proliferating"))

#Normalization ScaleData and UMAP Generation
Fibroblast <- FindVariableFeatures(Fibroblast, selection.method = "vst", nfeatures = 2000)
Fibroblast <- ScaleData(Fibroblast, verbose = T, features = row.names(Fibroblast))
Fibroblast <- RunPCA(object = Fibroblast)
stdev <- Fibroblast@reductions$pca@stdev
var <- stdev^2
sum(var[1:31])/ sum(var) #include a number that leads to 0.90 var

PCNum = 31

#Find Neighbors + Find Clusters (without harmony batch correction)
Fibroblast <- FindNeighbors(object = Fibroblast, dims = 1:PCNum)
Fibroblast <- FindClusters(object = Fibroblast, resolution = 0.8)

#Run UMAP and get unlabelled cluster UMAP and violin plot
Fibroblast<- RunUMAP(object = Fibroblast, dims = 1:PCNum)
DimPlot(object = Fibroblast, reduction = "umap", label = TRUE)

#Feature plot fibroblast subset
DotPlot(Fibroblast, features = c("Krt19", "Cdh1","Epcam","Try4", "Amy2a1", "Col1a2", "Pdpn", "Pdgfra", 
                             "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3","Il2ra", "Cd8a", "Nkg7","Arg1", 
                             "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33",
                             "H2-Eb1",    "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", 
                             "Ms4a1","Kit","Il1rl1","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

cluster_11 <- FindMarkers(Fibroblast, ident.1 = "11", only.pos = T)

#remove contaminating immune cells and ductal cells:
Fibroblast<-subset(Fibroblast, idents =  c("2", "8", "11", "12", "13"), invert = T)

#Normalization ScaleData and UMAP Generation
Fibroblast <- FindVariableFeatures(Fibroblast, selection.method = "vst", nfeatures = 2000)
Fibroblast <- ScaleData(Fibroblast, verbose = T, features = row.names(Fibroblast))
Fibroblast <- RunPCA(object = Fibroblast)
stdev <- Fibroblast@reductions$pca@stdev
var <- stdev^2
sum(var[1:32])/ sum(var) #include a number that leads to 0.90 var

PCNum = 32

#Find Neighbors + Find Clusters (without harmony batch correction)
Fibroblast <- FindNeighbors(object = Fibroblast, dims = 1:PCNum)
Fibroblast <- FindClusters(object = Fibroblast, resolution = 1.5)

#Run UMAP and get unlabelled cluster UMAP and violin plot
Fibroblast<- RunUMAP(object = Fibroblast, dims = 1:PCNum)
DimPlot(object = Fibroblast, reduction = "umap", label = TRUE, split.by = "Group")

#Feature plot fibroblast subset
DotPlot(Fibroblast, features = c("Krt19", "Cdh1","Epcam","Try4", "Amy2a1", "Col1a2", "Pdpn", "Pdgfra", 
                                 "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3","Il2ra", "Cd8a", "Nkg7","Arg1", 
                                 "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33",
                                 "H2-Eb1",    "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", 
                                 "Ms4a1","Kit","Il1rl1","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Label Clusters
Fibroblast<-subset(Fibroblast, idents =  "15", invert = T)

DotPlot(Fibroblast, features = c("Pdpn", "Pdgfra","Pdgfrb",
                                 "Acta2", "Col8a1", "Ccn2",
                                 "Clec3b", "Ly6c1", "Col14a1",
                                 "Saa3", "H2-Ab1", "Slpi",
                                 "Mki67", "Ccna2", "Cdk1"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

Idents(Fibroblast) <- "seurat_clusters"
Fibroblast <- RenameIdents(Fibroblast,
                           "0" = "myCAF", 
                           "1" = "myCAF", 
                           "2" = "myCAF", 
                           "3" = "Proliferating CAF", 
                           "4" = "myCAF", 
                           "5" = "myCAF", 
                           "6" = "myCAF", 
                           "7" = "iCAF", 
                           "8" = "apCAF", 
                           "9" = "Proliferating CAF",
                           "10" = "apCAF", 
                           "11" = "iCAF", 
                           "12" = "iCAF", 
                           "13" = "myCAF", 
                           "14" = "myCAF")
new_order <- c("myCAF", "iCAF", "apCAF", "Proliferating CAF")
Fibroblast@active.ident <- factor(Fibroblast@active.ident, levels = new_order)
Fibroblast[["caf_clusters"]] <- Fibroblast@active.ident

DimPlot(object = Fibroblast, reduction = "umap", label = FALSE, pt.size = 1, split.by = "Group", group.by = "caf_clusters")

Idents(Fibroblast) <- "Group"
new_order <- c("Vehicle", "TCDD")
Fibroblast@active.ident <- factor(Fibroblast@active.ident, levels = new_order)
Fibroblast[["Group"]] <- Fibroblast@active.ident 

save(Fibroblast,file="TCDD_Fibroblast_KLD.RData")

#### FIBROBLAST FIGURES ####
#Stacked Violin Plots (for markers and genes of interest):
Idents(Fibroblast) <- "caf_clusters"
Stacked_VlnPlot(seurat_object = Fibroblast, features = c("Pdpn", "Pdgfra","Pdgfrb",
                                                          "Clec3b", "Ly6c1", "Col14a1"), colors_use = c("maroon", "lightskyblue", "palegreen3", "grey"),  x_lab_rotate = TRUE)
Stacked_VlnPlot(seurat_object = Fibroblast, features = c("Acta2", "Col8a1", "Tagln",
                                                          "F11r", "H2-Ab1", "Slpi"), colors_use = c("maroon", "lightskyblue", "palegreen3", "grey"),  x_lab_rotate = TRUE) 


#cell population number histogram:
Idents(object = Fibroblast) <- 'Group'

samples.list <- unique(Fibroblast$Group)
clusters <- lapply(samples.list, function(x){
  subset <- subset(Fibroblast, subset = Group == x)
  dist <- data.frame(table(subset$caf_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})

clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("caf_clusters","variable","value","Group")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("caf_clusters","variable","value","Group")

ggplot(clusters_percent_dist, aes(fill=caf_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Vehicle", "TCDD")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance") 


#Fibroblast DE analysis:
Idents(Fibroblast) <- "caf_clusters"
icaf <- subset(Fibroblast, idents = "iCAF")
mycaf <- subset(Fibroblast, idents = "myCAF")
apcaf <- subset(Fibroblast, idents = "apCAF")

fibro_de <- FindMarkers(Fibroblast, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group")
icaf_de <- FindMarkers(icaf, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group")
mycaf_de <- FindMarkers(mycaf, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group")
apcaf_de <- FindMarkers(apcaf, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group")

write.csv(fibro_de, "All CAFs Vehicle vs TCDD DE KLD 2.5.23.csv")
write.csv(icaf_de, "iCAF Vehicle vs TCDD DE KLD 2.5.23.csv")
write.csv(mycaf_de, "myCAF Vehicle vs TCDD DE KLD 2.5.23.csv")
write.csv(apcaf_de, "apCAF Vehicle vs TCDD DE KLD 2.5.23.csv")

caf_de_cytokines <- c("Il19", "Ccl9", "Cxcl9", "Cxcl10", "Ccl24","Cx3cl1","Ccl19",
                      "Apoe", "Igf1", "Il33", "Spp1", "Cxcl16", "Il6",
                    "Cxcl12", "Wnt4")

DotPlot(Fibroblast, features = caf_de_cytokines, 
        cols = "RdYlBu", dot.scale = 8, split.by = "Group") + RotatedAxis()

Stacked_VlnPlot(seurat_object = Fibroblast, features =  c("Il19", "Ccl9", "Cxcl9", "Cxcl10", "Ccl24","Cx3cl1","Ccl19"), group.by = "Group")
                
Stacked_VlnPlot(seurat_object = Fibroblast, features =  c("Apoe", "Igf1", "Il33", "Spp1", "Cxcl16", "Il6", "Cxcl12", "Wnt4"), group.by = "Group")

Stacked_VlnPlot(seurat_object = Fibroblast, features = caf_de_cytokines, group.by = "Group")

#Fibroblast pathway analysis
h_gene_sets = msigdbr(species = "mouse", category = "H")
msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)

fibro_DE <- FindMarkers(Fibroblast, ident.1 = "Vehicle", ident.2 = "TCDD", group.by = "Group", logfc.threshold = 0, min.pct = 0)

ranks <- fibro_DE$avg_log2FC
names(ranks) <- rownames(fibro_DE)

fgseaRes <- fgsea(pathways = msigdbr_list, 
                  stats = ranks)

result <- apply(fgseaRes,2,as.character)

fgseaRes_sig <- subset(fgseaRes, fgseaRes$padj <= 0.05)

waterfall_plot_sig <- function (fgseaRes, graph_title) {
  (fgseaRes %>% 
     mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
     ggplot( aes(reorder(short_name,NES), NES)) +
     geom_bar(stat= "identity", aes(fill = padj))+
     coord_flip()+
     labs(x = "Hallmark Gene Set", y = "Normalized Enrichment Score", title = graph_title)+
     theme(axis.text.y = element_text(size = 7), 
           plot.title = element_text(hjust = 1)))+scale_fill_viridis(direction = 1, option = "G")
}
waterfall_plot_sig(fgseaRes_sig, "Pathways Enriched in WT vs. TCDD - Fibroblasts")

