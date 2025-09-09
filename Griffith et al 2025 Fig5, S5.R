#This script recreates the TCDD and Vehicle orthotopic scRNAseq analysis from Figures 5A-F and S5A-F of Griffith et al. 2025: 

#Data is available at NCBI GEO GSE303102

#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html

#Reference: Hao et al., Integrated analysis of multimodal single-cell data. 
#Cell. 2021 Jun 24;184(13):3573-3587.e29. doi: 10.1016/j.cell.2021.04.048. 
#PMID: 34062119; PMCID: PMC8238499. 

#Seurat Version 4.3.0.1
#SeuratObject Version 4.1.3
#R version 4.3.0 (2023-04-21) -- "Already Tomorrow"

#### Load Required Packages ####
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
library(tidyverse)
library(tibble) 
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(AUCell)

#Sample Info:
#BG-1:Vehicle 1
#BG-2:Vehicle 2
#BG-3:TCDD 1
#BG-4:TCDD 2
#Orthotopic Tumor Cell line-7940B

#### FULL OBJECT PREPROCESSING ####

#Load Data: 
Vehicle1.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-1/") 
Vehicle2.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-2/") 
TCDD1.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-3/") 
TCDD2.data <- Read10X("~/Dropbox (University of Michigan)/TCDD Orthotopic scRNA-seq/BG-4/") 

#Create Seurat objects:
Vehicle1  <- CreateSeuratObject(counts = Vehicle1.data, project = 'Vehicle1', min.cells = 3, min.features = 100)
Vehicle2 <- CreateSeuratObject(counts = Vehicle2.data, project = 'Vehicle2', min.cells = 3, min.features = 100)
TCDD1 <- CreateSeuratObject(counts = TCDD1.data, project = 'TCDD1', min.cells = 3, min.features = 100)
TCDD2 <- CreateSeuratObject(counts = TCDD2.data, project = 'TCDD2', min.cells = 3, min.features = 100)

#Group MetaData:
Vehicle1$Group <-"Vehicle"
Vehicle2$Group <-"Vehicle"
TCDD1$Group<-"TCDD"
TCDD2$Group <-"TCDD"

#Sample MetaData:
Vehicle1$Sample <-"Vehicle 1"
Vehicle2$Sample <-"Vehicle 2"
TCDD1$Sample<-"TCDD 1"
TCDD2$Sample <-"TCDD 2"

#Merge Seurat Objects: 
OrthotopicTCDD <- merge(Vehicle1, y =c(Vehicle2, TCDD1, TCDD2), 
                             add.cell.ids = c("Vehicle1", "Vehicle2","TCDD1","TCDD2"))

#Normalize Data: 
OrthotopicTCDD <- NormalizeData(OrthotopicTCDD)

#Apply Unbiased QC Cutoffs: 
OrthotopicTCDD[["percent.mt"]] <- PercentageFeatureSet(OrthotopicTCDD, pattern = "^mt-")
OrthotopicTCDD <- subset(OrthotopicTCDD, subset = nCount_RNA > 1200 & nCount_RNA < 100000 & percent.mt < 15)

#Find Variable Genes and Scale Data: 
OrthotopicTCDD <- FindVariableFeatures(OrthotopicTCDD, selection.method = "vst", nfeatures = 2000)
OrthotopicTCDD <- ScaleData(OrthotopicTCDD, verbose = T, features = row.names(OrthotopicTCDD))

#Run PCA and Dtermine Dimensions for 90% Variance: 
OrthotopicTCDD <- RunPCA(object = OrthotopicTCDD)
stdev <- OrthotopicTCDD@reductions$pca@stdev
var <- stdev^2
sum(var[1:29])/ sum(var) 
PCNum = 29

#Find Neighbors and Cluster Cells:
OrthotopicTCDD <- FindNeighbors(object = OrthotopicTCDD, dims = 1:PCNum)
OrthotopicTCDD <- FindClusters(object = OrthotopicTCDD, resolution = 3)

#Run UMAP: 
OrthotopicTCDD<- RunUMAP(object = OrthotopicTCDD, dims = 1:PCNum)
DimPlot(object = OrthotopicTCDD, reduction = "umap", label = T)

#Identify cell types: 
DotPlot(OrthotopicTCDD, features = c("Krt19", "Cdh1","Epcam","Msln","Try4", "Amy2a1", "Col1a2", "Acta2", "Clec3b", "Cspg4","Saa3","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                              "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Label clusters: 
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

#Organize cell types: 
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

#Save Seurat Object: 
save(OrthotopicTCDD,file="OrthotopicTCDD_KLD_Feb_24.RData")

DimPlot(OrthotopicTCDD, split.by = "Group", group.by = "global_clusters", cols = c()) #add colors

#### T CELL PROCESSING ####
#Subset and Re-cluster T cells: 
Idents(OrthotopicTCDD_KLD_Feb_24) <-"global_clusters"
Tcells <-subset(OrthotopicTCDD_KLD, idents =  c("T Cell", "NK Cell", "T Cell - Proliferating"))

#Find Variable Genes: 
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)

#Scale Data: 
all.genes <- rownames(Tcells)
Tcells <- ScaleData(Tcells, verbose = T, features = all.genes)

#Run PCA and Dtermine Dimensions for 90% Variance: 
Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE)
st_dev <- Tcells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:19])/ sum(var)

#Find neighbours and cluster cells: 
Tcells <- FindNeighbors(object = Tcells, dims = 1:19)
Tcells <- FindClusters(object = Tcells, resolution = 4.5)

#Run UMAP: 
Tcells <- RunUMAP(Tcells, reduction = "pca", dims = 1:19, verbose = F)
DimPlot(Tcells, reduction = "umap", label = T)

#Identify cell types: 
DotPlot(Tcells, features = c("Krt19", "Cdh1","Epcam","Try4", "Amy2a1", "Col1a2", "Pdpn", "Pdgfra", 
                             "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3","Il2ra", "Cd8a", "Nkg7","Arg1", 
                             "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33",
                             "H2-Eb1",    "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", 
                             "Ms4a1","Kit","Il1rl1","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Remove NK and myeloid cells
Tcells <- subset(Tcells, idents = c(0, 9, 11, 12, 15, 18, 21, 26, 27, 28, 30, 31, 32), invert = T)

#Find Variable Genes: 
Tcells <- FindVariableFeatures(Tcells, selection.method = "vst", nfeatures = 2000)

#Scale Data: 
all.genes <- rownames(Tcells)
Tcells <- ScaleData(Tcells, verbose = T, features = all.genes)

#Run PCA and Dtermine Dimensions for 90% Variance: 
Tcells <- RunPCA(Tcells, npcs = 30, verbose = FALSE)
st_dev <- Tcells@reductions$pca@stdev
var <- st_dev^2
sum(var[1:20])/ sum(var)

#Find neighbours and cluster cells: 
Tcells <- FindNeighbors(object = Tcells, dims = 1:20)
Tcells <- FindClusters(object = Tcells, resolution = 4)

#Run UMAP: 
Tcells <- RunUMAP(Tcells, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(Tcells, reduction = "umap", label = T)
DimPlot(Tcells, reduction = "umap", label = T, split.by = "Group")


#Identify cell types: 
DotPlot(Tcells, features = c("Krt19", "Cdh1","Epcam","Try4", "Amy2a1", "Col1a2", "Pdpn", "Pdgfra", 
                             "Cspg4","Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3","Il2ra", "Cd8a", "Nkg7","Arg1", 
                             "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33",
                             "H2-Eb1",    "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd79a", "Cd19", 
                             "Ms4a1","Kit","Il1rl1","Ccna2", "Mki67","Pecam1", "Cdh5", "Hbb-bt"), 
        cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

DotPlot(Tcells, features = c("Cd3e", "Cd3d","Sell", "Ccr7", "Tcf7","Il7r", "Cd28", "Cd4", "Tbx21", "Ifng", "Cxcr3", "Arg1","Gata3", "Il1rl1", "Il4", "Il13", "Rorc", "Il17a", "Il22", "Il2ra", "Foxp3", "Icos", "Cd8a", "Gzmb","Prf1", "Eomes","Lag3","Pdcd1","Trdc","Nkg7", "Klrg1", "Ccna2", "Mki67", "Ly6c2"), cols = "RdYlBu")+ RotatedAxis()

#Label clusters: 
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

#Organize cell types:
new_order <- c("CD8",
               "Proliferating CD8",
               "Treg",
               "Th1",
               "Proliferating CD4",
               "Memory CD4",
               "Naive CD4 & CD8")
Tcells@active.ident <- factor(Tcells@active.ident, levels = new_order)
Tcells[["t_clusters"]] <- Tcells@active.ident

#Simplify labeling:
Idents(Tcells) <- "t_clusters"
Tcells <- RenameIdents(Tcells,
                       "Naive CD4 & CD8" = "Naive CD4 & CD8", 
                       "Memory CD4" = "Memory CD4", 
                       "Proliferating CD8" = "CD8", 
                       "CD8" = "CD8", 
                       "Treg" = "Treg", 
                       "Proliferating CD4" = "Treg", 
                       "Th1" = "CD4")
Tcells[["simple_t_clusters"]] <- Tcells@active.ident

#Organize cell types:
new_order <- c("CD8",
               "Treg",
               "CD4",
               "Memory CD4",
               "Naive CD4 & CD8")
Tcells@active.ident <- factor(Tcells@active.ident, levels = new_order)
Tcells[["simple_t_clusters"]] <- Tcells@active.ident

#Save Seurat Object: 
save(Tcells,file="TCDD_Tcells_KLD_Feb_24.RData")

#T cell UMAP: 
plot_colors <- c("red",
                 "yellow",
                 "orange",
                 "green",
                 "blue")

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

#Plot cell abundance: 
ggplot(clusters_percent_dist, aes(fill=simple_t_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Vehicle", "TCDD")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("orig.ident") + ggtitle("Relative Cell Types Abundance") +scale_fill_manual(values = plot_colors)

#AhR gene signature in T cells: 
DotPlot(Tcells, features = c("Cyp1a1", "Cyp1b1", "Ahrr", "Tiparp"), split.by = "Group", group.by = "simple_t_clusters", cols = "RdYlBu", dot.scale = 12) + RotatedAxis()
FeaturePlot(Tcells, features = c("Cyp1b1"), cols = c("gainsboro", "firebrick1"), order = T)

#Gene set scoring for exhaustion, activation, memory, and naive signatures: 
#Gene list: 
gene_sets <- list(
  exhaustion = c("Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tigit", "Entpd1", "Tox", 
                 "Cd244", "Cd160", "Eomes", "Nr4a1", "Batf", "Timd4"),
  activation = c("Cd69", "Il2ra", "Cd40lg", "Ifng", "Tnf", "Prf1", "Gzmb", 
                 "Nfatc1", "Il2", "Icos", "Junb", "Cd27"),
  memory = c("Il7r", "Ccr7", "Cd27", "Cd28", "Sell", "Bcl2", "Fas", 
             "Eomes", "Tcf7", "Bcl6", "Hic1", "Il15ra"),
  naive = c("Tcf7", "Lef1", "Ccr7", "Sell", "Cd45ra", "Il7r", "Bach2", 
            "Gimap5", "Rag1", "Rag2", "Id3", "Dapl1")
)

#Gene Scoring: 
expr <- GetAssayData(Tcells,
                     assay = "RNA",
                     layer = "counts")
geneRanking <- AUCell_buildRankings(expr, nCores = 36, splitByBlocks = TRUE)
auc_score <- AUCell_calcAUC(rankings = geneRanking, geneSets = c(gene_sets), aucMaxRank = nrow(geneRanking)*0.1)
aucMatrix <- getAUC(auc_score)
auc_df <- as.data.frame(t(aucMatrix))
head(auc_df)
Tcells <- AddMetaData(Tcells, metadata = auc_df)
genesets <- c('exhaustion', 'activation', 'memory', 'naive')

#T cell UMAP for exhaustion, activation, memory, and naive signatures: 
FeaturePlot(Tcells, features = genesets, label = FALSE, repel = TRUE, 
            order = TRUE, pt.size = 1.0) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

#GSEA for exhaustion, activation, memory, and naive signatures:
#Subset CD8s:
Idents(Tcells) <- "simple_t_clusters" 
cd8_subset <- subset(Tcells, subset = simple_t_clusters == "CD8")
cd8_markers <- FindMarkers(cd8_subset, group.by = "Group", ident.1 = "TCDD", ident.2 = "Vehicle",test.use = "wilcox")
table(cd8_subset@meta.data$Group)

#Rank genes by avg_log2FC:
ranks_cd8 <- cd8_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%  
  deframe()

#Run GSEA:
gsea_results_cd8 <- fgseaMultilevel(
  pathways = gene_sets,
  stats = ranks_cd8,
  minSize = 5,
  maxSize = 500
)

gsea_results_cd8 <- as_tibble(gsea_results_cd8)

#Plot GSEA:
ggplot(gsea_results_cd8, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj)) +  
  scale_fill_gradientn(
    colors = brewer.pal(n = 10, name = "RdYlGn"), 
    guide = guide_colorbar(barwidth = 10, barheight = 1) 
  ) +
  coord_flip() +  
  labs(
    y = "Normalized Enrichment Score (NES)",  
    x = "Pathway",  
    title = "TCDD vs. Vehicle-CD8 T cells" 
  ) + 
  theme_minimal() +
  theme(legend.position = "bottom")

#GeneOntologies for CD8s:
significant_genes_tcdd <- rownames(cd8_markers[cd8_markers$avg_log2FC > 0 & cd8_markers$p_val_adj < 0.05, ])
significant_genes_vehicle <- rownames(cd8_markers[cd8_markers$avg_log2FC < 0 & cd8_markers$p_val_adj < 0.05, ])

o_tcdd <- enrichGO(
  gene          = significant_genes_tcdd,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

o_vehicle <- enrichGO(
  gene          = significant_genes_vehicle,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Convert enrichment results to data frames and get top 10 terms for CD8s:
go_results_tcdd <- as.data.frame(o_tcdd) %>%
  arrange(p.adjust) %>%
  head(10)

go_results_vehicle <- as.data.frame(o_vehicle) %>%
  arrange(p.adjust) %>%
  head(10)

go_results_tcdd$Group <- "TCDD"
go_results_vehicle$Group <- "Vehicle"

combined_go_results <- combined_go_results %>%
  mutate(CountDirection = ifelse(Group == "TCDD", Count, -Count))

#Plot GeneOntologies for CD8s:
ggplot(combined_go_results, aes(x = reorder(Description, abs(CountDirection)), y = CountDirection, fill = Group)) +
  geom_col() +  # No position "dodge" needed
  scale_fill_manual(values = c("TCDD" = "orange", "Vehicle" = "skyblue")) +  
  coord_flip() +  # Flip coordinates for horizontal bars
  geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
  labs(
    y = "Gene Count",
    x = "GO Term",
    title = "Top 10 Enriched GO Terms in TCDD vs. Vehicle-CD8 T cells"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

#Subset CD4s
Idents(Tcells) <- "simple_t_clusters" 
cd4_subset <- subset(Tcells, subset = simple_t_clusters == "CD4")
cd4_markers <- FindMarkers(cd4_subset, group.by = "Group", ident.1 = "TCDD", ident.2 = "Vehicle",test.use = "wilcox")
table(cd4_subset@meta.data$Group)

#Rank genes by avg_log2FC:
ranks_cd4 <- cd4_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%
  mutate(avg_log2FC = avg_log2FC + runif(n(), min = -1e-6, max = 1e-6)) %>%
  deframe()

#Run GSEA:
gsea_results_cd4 <- fgseaMultilevel(
  pathways = gene_sets,
  stats = ranks_cd4,
  minSize = 5,
  maxSize = 500
)

#Plot GSEA:
gsea_results_cd4 <- as_tibble(gsea_results_cd4)

ggplot(gsea_results_cd4, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj)) + 
  scale_fill_gradientn(
    colors = brewer.pal(n = 10, name = "RdYlGn"),  
    guide = guide_colorbar(barwidth = 10, barheight = 1)
  ) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",  
    x = "Pathway",  
    title = "TCDD vs. Vehicle-CD4 T cells"  
  ) + 
  theme_minimal() +
  theme(legend.position = "bottom")


#### CD8 PREPROCESSING ####
cd8s <- cd8_subset

#Find Variable Genes: 
cd8s <- FindVariableFeatures(cd8s, selection.method = "vst", nfeatures = 2000)

#Scale Data: 
all.genes <- rownames(cd8s)
cd8s <- ScaleData(cd8s, verbose = T, features = all.genes)

#Run PCA and Determine Dimensions for 90% Variance: 
cd8s <- RunPCA(cd8s, npcs = 30, verbose = FALSE)
st_dev <- cd8s@reductions$pca@stdev
var <- st_dev^2
sum(var[1:22])/ sum(var)

#Find Neighbors and Cluster Cells: 
cd8s <- FindNeighbors(object = cd8s, dims = 1:22)
cd8s <- FindClusters(object = cd8s, resolution = 4.5)

#Run UMAP: 
cd8s <- RunUMAP(cd8s, reduction = "pca", dims = 1:22, verbose = F)
DimPlot(cd8s, reduction = "umap", label = T)

#Identify Cell Populations: 
VlnPlot(cd8s, features = c("Cd3e","Cd4", "Cd8a", "Trdc"), ncol = 2)
DotPlot(cd8s, features = c("Cd3e", "Cd3d","Sell", "Ccr7", "Tcf7","Il7r", "Cd28", "Cd4", "Tbx21", "Ifng", "Cxcr3", "Arg1","Gata3", "Il1rl1", "Il4", "Il13", "Rorc", "Il17c",  "Il22", "Il2ra", "Foxp3", "Icos", "Cd8a", "Gzmb","Prf1", "Eomes","Lag3","Pdcd1","Trdc","Nkg7", "Klrg1", "Ccna2", "Mki67", "Ly6c2"), cols = "RdYlBu")+ RotatedAxis()

cd8s <- subset(cd8s, idents = c("15","16"), invert = T)

#Vln Plots for CD8 activation markers: 
VlnPlot(cd8s, features = c("Ifng", "Prf1", "Gzmb"), group.by = "Group")

#CD8 HEATMAP: 
cd8_genes <-  c("Gzmb", "Ifng", "Prf1","Tox", "Nkg7", "Lag3", "Timp3", "Pdcd1", "Ctla4")
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
