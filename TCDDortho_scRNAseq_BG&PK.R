library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db) 
library(org.Mm.eg.db) 
library(DOSE)
library(dplyr)
library(ReactomePA)
library(ggplot2)
library(RColorBrewer)
library(colorRamps)

#Sample Info
#BG-1:Vehicle 1
#BG-2:Vehicle 2
#BG-3:TCDD 1
#BG-4:TCDD 2
#Orthotopic Tumor Cell line-7940B

#Load Dataset, make sure that genes are renamed to features
Vehicle1.data <- Read10X("~/Desktop/TCDD Orthotopic scRNA-seq/BG-1") 
Vehicle2.data <- Read10X("~/Desktop/TCDD Orthotopic scRNA-seq/BG-2") 
TCDD1.data <- Read10X("~/Desktop/TCDD Orthotopic scRNA-seq/BG-3") 
TCDD2.data <- Read10X("~/Desktop/TCDD Orthotopic scRNA-seq/BG-4") 

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
levels(OrthotopicTCDD)

# Changing between meta.data for identities- you can change this by altering what you input into your metadata
Idents(object = OrthotopicTCDD) <- 'Group'
#Check active Identity
levels(OrthotopicTCDD)

#Percent Mito
OrthotopicTCDD[["percent.mt"]] <- PercentageFeatureSet(OrthotopicTCDD, pattern = "^mt-")
VlnPlot(OrthotopicTCDD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .01)

Idents(OrthotopicTCDD) <- "Group"
table(OrthotopicTCDD@active.ident)

#Vehicle   #TCDD
#26331    #27717

#Filter genes further
OrthotopicTCDD <- subset(OrthotopicTCDD, subset = nCount_RNA > 1200 & nCount_RNA < 100000 & percent.mt < 15)
VlnPlot(OrthotopicTCDD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = .01)

Idents(OrthotopicTCDD) <- "Group"
table(OrthotopicTCDD@active.ident)

#Vehicle    TCDD 
#  20995   20922 


#Normalization ScaleData and UMAP Generation
OrthotopicTCDD <- NormalizeData(OrthotopicTCDD)
OrthotopicTCDD <- FindVariableFeatures(OrthotopicTCDD, selection.method = "vst", nfeatures = 2000)
OrthotopicTCDD <- ScaleData(OrthotopicTCDD, verbose = T, features = row.names(OrthotopicTCDD))
OrthotopicTCDD <- RunPCA(object = OrthotopicTCDD)
stdev <- OrthotopicTCDD@reductions$pca@stdev
var <- stdev^2
sum(var[1:30])/ sum(var) #include a number that leads to 0.90 var

PCNum = 30 #put the dimensions here

#Find Neighbors + Find Clusters (without harmony batch correction)

OrthotopicTCDD <- FindNeighbors(object = OrthotopicTCDD, dims = 1:PCNum)
OrthotopicTCDD <- FindClusters(object = OrthotopicTCDD, resolution = 1.8)

#Run UMAP and get unlabelled cluster UMAP and violin plot
OrthotopicTCDD<- RunUMAP(object = OrthotopicTCDD, dims = 1:PCNum)
DimPlot(object = OrthotopicTCDD, reduction = "umap", label = FALSE, pt.size = 1, split.by = 'Group')

save(OrthotopicTCDD,file="OrthotopicTCDD.RData")
FeaturePlot(object = OrthotopicTCDD, features = c("Cd3e","Cd8a","Cd4","Nkg7"), cols = c("grey", "red"), reduction = "umap", pt.size = .5)# T cells 
FeaturePlot(object = OrthotopicTCDD, features = c("Mki67"), cols = c("grey", "red"), reduction = "umap", pt.size = .5)#proliferating Cells 
FeaturePlot(object = OrthotopicTCDD, features = c("Col1a2","Pdgfra","Pdgfrb"), cols = c("grey", "red"), reduction = "umap", pt.size = .5)#Fibroblast
FeaturePlot(object = OrthotopicTCDD, features = c("Krt19", "Muc1", "Krt18","Cdh1","Msln"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Epithelial
FeaturePlot(object = OrthotopicTCDD, features = c("Epcam","Krt19", "Muc1", "Krt18","Cdh1","Msln", "Spink1", "Amy1","Prss1"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Epithelial
FeaturePlot(object = OrthotopicTCDD, features = c("S100a9","S100a8","Mmp9","Il1r2"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Neutrophils
FeaturePlot(object = OrthotopicTCDD, features = c("Vwf", "Cdh5",'Erg'), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Endothelial
FeaturePlot(object = OrthotopicTCDD, features = c("Ms4a1", "Cd19", "Cd79a"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #B cells
FeaturePlot(object = OrthotopicTCDD, features = c("C1qb","C1qa","C1qc","Lyz2"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #TAMs
FeaturePlot(object = OrthotopicTCDD, features = c("Try4","Amy2a2"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Acinar 
FeaturePlot(object = OrthotopicTCDD, features = c("Itgae","Xcr1","Clec9a","Cd209a"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Dcs
FeaturePlot(object = OrthotopicTCDD, features = c("Itgam"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Dcs
FeaturePlot(object = OrthotopicTCDD, features = c("Ptprc"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #Immune Cells

#Lable Clusters
Idents(OrthotopicTCDD) <- "Group"
new_order <- c("Vehicle", "TCDD")
OrthotopicTCDD@active.ident <- factor(OrthotopicTCDD@active.ident, levels = new_order)
OrthotopicTCDD[["Group"]] <- OrthotopicTCDD@active.ident 
Idents(OrthotopicTCDD) <- "seurat_clusters"
OrthotopicTCDD <- RenameIdents(OrthotopicTCDD,
                          "0" = "TAMs", 
                          "1" = "TAMs", 
                          "2" = "Fibroblast 1", 
                          "3" = "Fibroblast 1", 
                          "4" = "Neutrophils", 
                          "5" = "TAMs", 
                          "6" = "TAMs", 
                          "7" = "Epithelial", 
                          "8" = "TAMs", 
                          "9" = "Proliferating Epithelial", 
                          "10" = "T cells", 
                          "11" = "Proliferating TAMs",
                          "12" = "TAMs", 
                          "13" = "Epithelial", 
                          "14" = "Epithelial", 
                          "15" = "Epithelial", 
                          "16" = "Endothelial", 
                          "17" = "Proliferating Epithelial", 
                          "18" = "Proliferating Fibroblast 1", 
                          "19" = "Epithelial", 
                          "20" = "Fibroblast 1", 
                          "21" = "Fibroblast 1", 
                          "22" = "B Cells", 
                          "23" = "Neutrophils",
                          "24" = "Neutrophils", 
                          "25" = "DCs", 
                          "26" = "Fibroblast 2", 
                          "27" = "RBC", 
                          "28" = "Acinar", 
                          "29" = "NK cells", 
                          "30" = "Proliferating T cells", 
                          "31" = "Unk Myeloid", 
                          "32" = "Myeloid & T cell mix", 
                          "33" = "TAMs", 
                          "34" = "Fibroblast 1",
                          "35" = "Endothelial", 
                          "36" = "RBC", 
                          "37" = "TAMs",
                          "38" = "Proliferating Fibroblast 1", 
                          "39" = "Duct", 
                          "40" = "TAMs",
                          "41" = "TAMs", 
                          "42" = "RBC", 
                          "43" = "Unk Myeloid",
                          "44" = "Endothelial", 
                          "45" = "Endothelial")
OrthotopicTCDD[["collapsed_clusters"]] <- OrthotopicTCDD@active.ident
DimPlot(object = OrthotopicTCDD_1, reduction = "umap", label = FALSE, pt.size = 1, split.by = "Group")
save(OrthotopicTCDD,file="OrthotopicTCDD.RData")

#For Main Fig UMAP
Idents(OrthotopicTCDD) <-"collapsed_clusters"
levels(OrthotopicTCDD)
OrthotopicTCDD_1 <-subset(OrthotopicTCDD, idents =  c("TAMs","Fibroblast 1","Neutrophils","Epithelial","Proliferating Epithelial","T cells","Proliferating TAMs","Endothelial","Proliferating Fibroblast 1","B Cells","DCs","Fibroblast 2","Acinar","NK cells","Proliferating T cells","Duct"))
order <- c("T cells","Fibroblast 1","Fibroblast 2","Neutrophils","TAMs","B Cells","DCs","NK cells","Duct","Epithelial","Endothelial","Acinar","Proliferating Epithelial","Proliferating Fibroblast 1","Proliferating T cells","Proliferating TAMs")
OrthotopicTCDD_1 @active.ident <- factor(OrthotopicTCDD_1 @active.ident, levels = order)
DimPlot(object = OrthotopicTCDD_1, reduction = "umap", label = F, pt.size = 0.5, cols = c("TAMs"= "gold",
                                                                                 "Fibroblast 1"="lightblue1",
                                                                                 "Neutrophils"="darkorange",
                                                                                 "Epithelial"="rosybrown2",
                                                                                 "Proliferating Epithelial"= "brown",
                                                                                 "T cells"= "royalblue",
                                                                                 "Proliferating TAMs"="indianred2",
                                                                                 "Endothelial"="mediumpurple",
                                                                                 "Proliferating Fibroblast 1"= "orchid1", 
                                                                                 "B Cells"="grey20",
                                                                                 "DCs"="green",
                                                                                 "Fibroblast 2"="plum",
                                                                                 "Acinar"="green4",
                                                                                 "NK cells"="gold3",
                                                                                 "Proliferating T cells"="cyan",
                                                                                 "Duct"="wheat"))


#To generate a Heatmap subset Cells
Idents(OrthotopicTCDD) <-"collapsed_clusters"
levels(OrthotopicTCDD)
Tcells <-subset(OrthotopicTCDD, idents =  c("T cells"))
Fibro <-subset(OrthotopicTCDD, idents =  c("Fibroblast 1", "Fibroblast 2"))
TAMs <-subset(OrthotopicTCDD, idents =  c("TAMs"))
Neutro <-subset(OrthotopicTCDD, idents =  c("Neutrophils"))

#Tcells Heatmap
Idents(Tcells) <-"Group"
levels(Tcells)       
Tcells.markers <- FindAllMarkers(Tcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Tcells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rp_mt_genes <- Tcells.markers$gene[grep("^ro|^mt-", Tcells.markers$gene)]
Tcells.markers<- Tcells.markers %>% filter(!gene %in% rp_mt_genes)
top10Tcells <- Tcells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Tcells, features = top10Tcells$gene) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
VlnPlot(object = Tcells, features = c("Foxp3","Tcf7"))

#Fibro Heatmap
Idents(Fibro) <-"Group"
levels(Fibro)       
Fibro.markers <- FindAllMarkers(Fibro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Fibro.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rp_mt_genes <- Fibro.markers$gene[grep("^ro|^mt-", Fibro.markers$gene)]
Fibro.markers<- Fibro.markers %>% filter(!gene %in% rp_mt_genes)
top10Fibro <- Fibro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Fibro, features = top10Fibro$gene) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))

#TAMs Heatmap
Idents(TAMs) <-"Group"
levels(TAMs)       
TAMs.markers <- FindAllMarkers(TAMs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TAMs.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rp_mt_genes <- TAMs.markers$gene[grep("^ro|^mt-", TAMs.markers$gene)]
TAMs.markers<- TAMs.markers %>% filter(!gene %in% rp_mt_genes)
top10TAMs <- TAMs.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(TAMs, features = top10TAMs$gene) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
VlnPlot(object = TAMs, features = c("Cxcl9","Tph1"))

#Neutro Heatmap
Idents(Neutro) <-"Group"
levels(Neutro)       
Neutro.markers <- FindAllMarkers(Neutro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Neutro.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rp_mt_genes <- Neutro.markers$gene[grep("^ro|^mt-", Neutro.markers$gene)]
Neutro.markers<- Neutro.markers %>% filter(!gene %in% rp_mt_genes)
top10Neutro <- Neutro.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Neutro, features = top10Neutro$gene) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
VlnPlot(object = Neutro, features = c("Cxcl3"))

#Cell Prop
Idents(object = OrthotopicTCDD_1) <- 'collapsed_clusters'
levels(OrthotopicTCDD_1)

#splitting samples
samples.list <- unique(OrthotopicTCDD_1$Group)
clusters <- lapply(samples.list, function(x){
  subset <- subset(OrthotopicTCDD_1, subset = Group == x)
  dist <- data.frame(table(subset$collapsed_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

#calculate relative freq (fractions) of each cell type
clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})
Idents(object = OrthotopicTCDD_1) <- "Group"
levels(OrthotopicTCDD_1)
#making things ggplot-friendly!
clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("collapsed_clusters","variable","value","Group")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("collapsed_clusters","variable","value","Group")


ggplot(clusters_percent_dist, aes(fill=collapsed_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Vehicle","TCDD")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("Group") + ggtitle("Relative Cell Types Abundance") +
  scale_fill_manual(values = c("gold","lightblue1","darkorange","rosybrown2","brown", "royalblue","indianred2","mediumpurple","orchid1","grey20","green","plum","green4","gold3","cyan","wheat",'grey','grey','grey'))

DotPlot(OrthotopicTCDD_1, features = c("Cd3e","Pdgfra","Pdgfrb","Dcn","Col1a2","Itgam","S100a8","S100a9","Mmp9","C1qa","C1qb","Lyz2","Cd79a","Cd19","Itgax","Nkg7","Ceacam1","Krt19","Cdh1","Muc1","Msln","Vwf","Cdh5","Spink1","Amy1","Mki67"),cols = "RdBu", dot.scale = 8) + RotatedAxis()



######Subset fibroblasts
Idents(OrthotopicTCDD) <-"collapsed_clusters"
levels(OrthotopicTCDD)
Fibroblast<-subset(OrthotopicTCDD, idents =  c("Fibroblast 1", "Fibroblast 2"))

#Normalization ScaleData and UMAP Generation
Fibroblast <- NormalizeData(Fibroblast)
Fibroblast <- FindVariableFeatures(Fibroblast, selection.method = "vst", nfeatures = 2000)
Fibroblast <- ScaleData(Fibroblast, verbose = T, features = row.names(Fibroblast))
Fibroblast <- RunPCA(object = Fibroblast)
stdev <- Fibroblast@reductions$pca@stdev
var <- stdev^2
sum(var[1:32])/ sum(var) #include a number that leads to 0.90 var

PCNum = 32 #put the dimensions here

#Find Neighbors + Find Clusters (without harmony batch correction)
Fibroblast <- FindNeighbors(object = Fibroblast, dims = 1:PCNum)
Fibroblast <- FindClusters(object = Fibroblast, resolution = 0.8)

#Run UMAP and get unlabelled cluster UMAP and violin plot
Fibroblast<- RunUMAP(object = Fibroblast, dims = 1:PCNum)
DimPlot(object = Fibroblast, reduction = "umap", label = TRUE, pt.size = 1, split.by = 'Group')

#Feature plot fibrolabst subset
FeaturePlot(object = Fibroblast, features = c("Col1a2","Pdgfra","Pdgfrb", "Mki67"), cols = c("grey", "red"), order = TRUE, reduction = "umap", pt.size = .5)
FeaturePlot(object = Fibroblast, features = c("Clec3b","Il6","Col14a1", "Has1"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #iCAF
FeaturePlot(object = Fibroblast, features = c("Tagln", "Thy1", "Col12a1", "Thbs2", "Acta2"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #myCAF
FeaturePlot(object = Fibroblast, features = c("H2-Ab1", "Cd74", "Saa3", "Slpi"), cols = c("grey", "red"), reduction = "umap", pt.size = .5) #apCAF

VlnPlot(object = Fibroblast, features = c("Acta2"))
VlnPlot(object = Fibroblast, features = c("Tagln"))
VlnPlot(object = Fibroblast, features = c("Col12a1"))
VlnPlot(object = Fibroblast, features = c("Cd74"))

#Fibro Heatmap
Idents(Fibroblast) <- "seurat_clusters"
levels(Fibroblast)       
Fibroblast.markers <- FindAllMarkers(Fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rp_mt_genes <- Fibroblast.markers$gene[grep("^ro|^mt-", Fibroblast.markers$gene)]
Fibroblast.markers<- Fibroblast.markers %>% filter(!gene %in% rp_mt_genes)
top10Fibroblast <- Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Fibroblast, features = top10Fibroblast$gene) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
VlnPlot(object = Fibroblast, features = c("Cxcl3"))

#Save fibroblast object
save(Fibroblast,file="Fibroblast.RData")

#Lable Clusters
Idents(Fibroblast) <- "Group"
new_order <- c("Vehicle", "TCDD")
Fibroblast@active.ident <- factor(Fibroblast@active.ident, levels = new_order)
Fibroblast[["Group"]] <- Fibroblast@active.ident 
Idents(Fibroblast) <- "seurat_clusters"
Fibroblast <- RenameIdents(Fibroblast,
                               "0" = "Fibro1", 
                               "1" = "Fibro2", 
                               "2" = "TCDD fibro", 
                               "3" = "TCDD fibro", 
                               "4" = "Fibro3", 
                               "5" = "Fibro4", 
                               "6" = "Fibro5", 
                               "7" = "Fibro6", 
                               "8" = "Fibro7", 
                               "9" = "Fibro8", 
                               "10" = "Proliferating fibro", 
                               "11" = "Fibro9",
                               "12" = "Fibro10", 
                               "13" = "Proliferating fibro", 
                               "14" = "Fibro11" )
Fibroblast[["collapsed_clusters"]] <- Fibroblast@active.ident
DimPlot(object = Fibroblast, reduction = "umap", label = FALSE, pt.size = 1, split.by = "Group")
save(Fibroblast,file="Fibroblast.RData")

#Rerun heatmap with collapsed fibro labels
Idents(Fibroblast) <- "collapsed_clusters"
levels(Fibroblast)       
Fibroblast.markers <- FindAllMarkers(Fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rp_mt_genes <- Fibroblast.markers$gene[grep("^ro|^mt-", Fibroblast.markers$gene)]
Fibroblast.markers<- Fibroblast.markers %>% filter(!gene %in% rp_mt_genes)
top10Fibroblast <- Fibroblast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Fibroblast, features = top10Fibroblast$gene) +  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
VlnPlot(object = Fibroblast, features = c("Cyp1b1"), split.by = "Group", cols=c("lightblue", "blue"))
VlnPlot(object = Fibroblast, features = c("Cxcl14"), split.by = "Group", cols=c("lightblue", "blue"))
VlnPlot(object = Fibroblast, features = c("Mmp3"), split.by = "Group", cols=c("lightblue", "blue"))

#Find number of T cells
Idents(OrthotopicTCDD) <- "Group"
table(OrthotopicTCDD@active.ident)



