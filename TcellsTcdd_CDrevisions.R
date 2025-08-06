#Griffith et al., Cancer Discovery Revision T cell analysis 

library(Seurat)
library(dplyr)
library(pheatmap)
library(tidyr)
library(ggplot2)
library(DT)
library(RColorBrewer)
library(AUCell)
library(cowplot)
library(msigdbr)
library(fgsea)
library(knitr)
library(stringr)

DimPlot(object = OrthotopicTCDD, reduction = "umap", label = T, pt.size = 0.25) # Obj from Feb24
FeaturePlot(object = OrthotopicTCDD, features = c("Cd3e"), cols = c("grey", "blue"), reduction = "umap", pt.size = .5)
DotPlot(OrthotopicTCDD, features = rev(c("Cd3e","Mki67","Cd4","Cd8a")), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
Idents(OrthotopicTCDD) <-"seurat_clusters"       #Clusters 22,30,44 are nonprolif T cells, proliferating T cells = cluster 42

#Subset T cells 
Tcells <- subset(OrthotopicTCDD, idents = c("22","30","44"))
#Find Variable Genes:
Tcells <- FindVariableFeatures(object = Tcells, selection.method = "vst", nfeatures = 2000)
#Scale Data, generates "Scale.Data" matrix:
Tcells <- ScaleData(object =Tcells, vars.to.regress = c("nCount_RNA"), features = rownames(Tcells))
#Run PCA:
Tcells <- RunPCA(object = Tcells)
stdev <- Tcells@reductions$pca@stdev
var <- stdev^2
sum(var[1:37])/ sum(var)
PCNum = 37 #put the dimensions here
#Find Neighbors + Find Clusters 
Tcells <- FindNeighbors(object = Tcells, dims = 1:PCNum)
Tcells <- FindClusters(object = Tcells, resolution = 1)
Tcells<- RunUMAP(object = Tcells, dims = 1:PCNum)
DimPlot(object = Tcells, reduction = "umap", label = T, pt.size = 1)
DimPlot(object = Tcells, reduction = "umap", label = F, pt.size = 0.5, cols = c("skyblue", "orange"))
Idents(Tcells) <-"seurat_clusters"
levels(Tcells)
save(Tcells,file="TCDDTcells.RData")

#CD4 and CD8 T cell identification 
FeaturePlot(object = Tcells, features = c("Cd4","Cd8a"), cols = c("grey", "purple"), reduction = "umap", pt.size = 1) #
DotPlot(Tcells, features = c("Cd3e","Cd4","Cd8a"), cols = "RdBu", dot.scale = 8) + RotatedAxis()
Idents(Tcells) <-"seurat_clusters"
Tcells <- RenameIdents(Tcells,"0" = "Naive CD8", 
                                "1" = "Treg", 
                                "2" = "Exhausted CD8", 
                                "3" = "Exhausted CD4&CD8 DP", 
                                "4" = "Exhausted CD8", 
                                "5" = "Exhausted CD4&CD8 DP", 
                                "6" = "Naive Cd4&Cd8 DP", 
                                "7" = "Naive CD8", 
                                "8" = "Naive CD4",
                                "9" = "Memory CD4", 
                                "10" = "Memory CD8", 
                                "11" = "Not Tcells")
Tcells[["T_clusters"]] <- Tcells@active.ident
DimPlot(object = Tcells, reduction = "umap", label = T, pt.size = 1)
DotPlot(Tcells, features = c('Cd3e','Cd4',"Cd8a",'Foxp3','Sell','Tcf7','Ccr7',"Bcl6",'Trdc','Gzmb',"Tbx21","Prf1","Ifng","Pdcd1","Eomes") ,cols = "RdYlBu", dot.scale = 8) + RotatedAxis()

#Reordering
Idents(Tcells) <- factor(Idents(Tcells), levels = c("Naive CD4", "Naive CD8", "Naive Cd4&Cd8 DP","Treg","Exhausted Cd4","Exhausted CD8","Exhausted CD4&CD8 DP","Memory CD4","Memory CD8","Not Tcells"))
DimPlot(Tcells, reduction = "umap",label = T) + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

Tcells$T_clusters <- factor(Idents(Tcells), levels = c("Naive CD4", "Naive CD8", "Naive Cd4&Cd8 DP","Treg","Exhausted Cd4","Exhausted CD8","Exhausted CD4&CD8 DP","Memory CD4","Memory CD8","Not Tcells"))
DotPlot(Tcells, features = c('Cd3e','Cd4',"Cd8a",'Foxp3','Sell','Tcf7','Ccr7',"Bcl6",'Trdc','Gzmb',"Tbx21","Prf1","Ifng","Pdcd1","Tox","Eomes") ,cols = "RdYlBu", dot.scale = 8) + RotatedAxis()




#GSEA analysis for Exhaustion, Activation, Memory, and Naive Signatures 

library(tidyverse)
library(dplyr)    
library(tibble)    
library(fgsea)
library(ggplot2)
library(dplyr)
library(tibble)

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


# CD8s
cd8_subset <- subset(Tcells, subset = T_clusters == "Cd8")
cd8_markers <- FindMarkers(cd8_subset, group.by = "Group", ident.1 = "TCDD", ident.2 = "Vehicle",test.use = "wilcox")
table(cd8_subset@meta.data$Group)


# Rank genes by avg_log2FC
ranks_cd8 <- cd8_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%  
  deframe()

# Run GSEA
gsea_results_cd8 <- fgseaMultilevel(
  pathways = gene_sets,
  stats = ranks_cd8,
  minSize = 5,
  maxSize = 500
)

gsea_results_cd8 <- as_tibble(gsea_results_cd8)

ggplot(gsea_results_cd8, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj)) +  # Color bars based on adjusted p-value
  scale_fill_gradientn(
    colors = brewer.pal(n = 10, name = "RdYlGn"),  # Use color palette 'RdYlGn'
    guide = guide_colorbar(barwidth = 10, barheight = 1)  # Customizing the color bar
  ) +
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(
    y = "Normalized Enrichment Score (NES)",  # Label for Y-axis
    x = "Pathway",  # Label for X-axis
    title = "TCDD vs. Vehicle-CD8 T cells"  # Plot title
  ) + 
  theme_minimal() +
  theme(legend.position = "bottom")


# GeneOntologies 

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)

significant_genes_tcdd <- rownames(cd8_markers[cd8_markers$avg_log2FC > 0 & cd8_markers$p_val_adj < 0.05, ])
significant_genes_vehicle <- rownames(cd8_markers[cd8_markers$avg_log2FC < 0 & cd8_markers$p_val_adj < 0.05, ])

ego_tcdd <- enrichGO(
  gene          = significant_genes_tcdd,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# GO enrichment for Vehicle DEGs
ego_vehicle <- enrichGO(
  gene          = significant_genes_vehicle,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Convert enrichment results to data frames and get top 10 terms
go_results_tcdd <- as.data.frame(ego_tcdd) %>%
  arrange(p.adjust) %>%
  head(10)

go_results_vehicle <- as.data.frame(ego_vehicle) %>%
  arrange(p.adjust) %>%
  head(10)

# Add group information for plotting
go_results_tcdd$Group <- "TCDD"
go_results_vehicle$Group <- "Vehicle"

# Combine results
combined_go_results <- combined_go_results %>%
  mutate(CountDirection = ifelse(Group == "TCDD", Count, -Count))

ggplot(combined_go_results, aes(x = reorder(Description, abs(CountDirection)), y = CountDirection, fill = Group)) +
  geom_col() +  # No position "dodge" needed
  scale_fill_manual(values = c("TCDD" = "orange", "Vehicle" = "skyblue")) +  # Set colors for groups
  coord_flip() +  # Flip coordinates for horizontal bars
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Line to separate groups
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


# CD4s
cd4_subset <- subset(Tcells, subset = T_clusters == "Cd4")
cd4_markers <- FindMarkers(cd4_subset, group.by = "Group", ident.1 = "TCDD", ident.2 = "Vehicle",test.use = "wilcox")
table(cd4_subset@meta.data$Group)

set.seed(42)  # Ensure reproducibility
ranks_cd4 <- cd4_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%
  mutate(avg_log2FC = avg_log2FC + runif(n(), min = -1e-6, max = 1e-6)) %>%
  deframe()

#GSEA with adjusted ranks
gsea_results_cd4 <- fgseaMultilevel(
  pathways = gene_sets,
  stats = ranks_cd4,
  minSize = 5,
  maxSize = 500
)

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
    x = "Pathway",  # Label for X-axis
    title = "TCDD vs. Vehicle-CD4 T cells"  
  ) + 
  theme_minimal() +
  theme(legend.position = "bottom")


# GeneOntologies 

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(RColorBrewer)

significant_genes_tcdd <- rownames(cd4_markers[cd4_markers$avg_log2FC > 0 & cd4_markers$p_val_adj < 0.05, ])
significant_genes_vehicle <- rownames(cd4_markers[cd4_markers$avg_log2FC < 0 & cd4_markers$p_val_adj < 0.05, ])

ego_tcdd <- enrichGO(
  gene          = significant_genes_tcdd,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# GO enrichment for Vehicle DEGs
ego_vehicle <- enrichGO(
  gene          = significant_genes_vehicle,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Convert enrichment results to data frames and get top 10 terms
go_results_tcdd <- as.data.frame(ego_tcdd) %>%
  arrange(p.adjust) %>%
  head(10)

go_results_vehicle <- as.data.frame(ego_vehicle) %>%
  arrange(p.adjust) %>%
  head(10)

# Add group information for plotting
go_results_tcdd$Group <- "TCDD"
go_results_vehicle$Group <- "Vehicle"

# Combine results
combined_go_results <- combined_go_results %>%
  mutate(CountDirection = ifelse(Group == "TCDD", Count, -Count))

ggplot(combined_go_results, aes(x = reorder(Description, abs(CountDirection)), y = CountDirection, fill = Group)) +
  geom_col() +  # No position "dodge" needed
  scale_fill_manual(values = c("TCDD" = "orange", "Vehicle" = "skyblue")) +  # Set colors for groups
  coord_flip() +  # Flip coordinates for horizontal bars
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +  # Line to separate groups
  labs(
    y = "Gene Count",
    x = "GO Term",
    title = "Top 10 Enriched GO Terms in TCDD vs. Vehicle-CD4 T cells"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )


#AUCell
library(AUCell)
expr <- GetAssayData(Tcells,
                     assay = "RNA",
                     layer = "counts")
geneRanking <- AUCell_buildRankings(expr, nCores = 36, splitByBlocks = TRUE)
auc_score <- AUCell_calcAUC(rankings = geneRanking, geneSets = c(gene_sets), aucMaxRank = nrow(geneRanking)*0.1)
aucMatrix <- getAUC(auc_score)
auc_df <- as.data.frame(t(aucMatrix))
head(auc_df)
Tcells <- AddMetaData(Tcells, metadata = auc_df)
VlnPlot(Tcells, features = genesets, split.by = 'Group', pt.size = 1) +
  scale_fill_brewer(palette = "RdYlBu", direction = -1)

genesets <- c('exhaustion', 'activation', 'memory', 'naive')

FeaturePlot(Tcells, features = genesets, label = FALSE, repel = TRUE, 
            order = TRUE, pt.size = 1.0) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))






#2 Gsea for pathways 
#cd8s
gene_sets_db <- msigdbr(species = "Mus musculus", category = "C2")
keywords <- "T_CELL|T lymphocyte|ACTIVATION|EXHAUSTION|MEMORY|NAIVE"
tcell_gene_sets <- gene_sets_db %>%
  filter(grepl(keywords, gs_name, ignore.case = TRUE))

# Create a list of gene sets for the filtered pathways
tcell_gene_sets_list <- split(x = tcell_gene_sets$gene_symbol, f = tcell_gene_sets$gs_name)
ranks_cd8 <- cd8_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%
  deframe()

# Run GSEA with the processed pathways
set.seed(123)
gsea_results_cd8 <- fgsea(pathways = tcell_gene_sets_list, stats = ranks_cd8, minSize = 15, maxSize = 500, nPermSimple = 10000)

# Select and print top pathway results for inspection
top_tcell_pathways <- gsea_results_cd8 %>%
  arrange(padj) %>%
  head(10)  # Adjust as needed for deeper insights

print(top_tcell_pathways)

# Plot the GSEA results
ggplot(top_tcell_pathways, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
  geom_col() +
  scale_fill_gradientn(
    colors = brewer.pal(n = 10, name = "RdYlGn"),
    guide = guide_colorbar(barwidth = 10, barheight = 1),
    name = "Adjusted P-Value"
  ) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",
    x = "Pathway",
    title = "Top 10 cell Functional Pathways for CD8 T Cells (TCDD vs. Vehicle)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )



#cd4s
# Load C2 gene sets for mouse
gene_sets_db <- msigdbr(species = "Mus musculus", category = "C2")
keywords <- "T_CELL|T lymphocyte|ACTIVATION|EXHAUSTION|MEMORY|NAIVE"
tcell_gene_sets <- gene_sets_db %>%
  filter(grepl(keywords, gs_name, ignore.case = TRUE))

# Create a list of gene sets for the filtered pathways
tcell_gene_sets_list <- split(x = tcell_gene_sets$gene_symbol, f = tcell_gene_sets$gs_name)

# Assuming `ranks_cd4` is prepared from your differential expression results
# Handle ties by adding small noise, if necessary
ranks_cd4 <- cd4_markers %>%
  rownames_to_column(var = "gene") %>%
  filter(!is.na(avg_log2FC)) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC) %>%
  mutate(avg_log2FC = avg_log2FC + runif(n(), min = -1e-5, max = 1e-5)) %>%
  deframe()

# Run GSEA with the processed pathways
set.seed(123)
gsea_results_cd4 <- fgsea(pathways = tcell_gene_sets_list, stats = ranks_cd4, minSize = 15, maxSize = 500, nPermSimple = 10000)

# Select and print top pathway results for inspection
top_tcell_pathways <- gsea_results_cd4 %>%
  arrange(padj) %>%
  head(10)  # Adjust as needed for deeper insights

print(top_tcell_pathways)

# Plot the GSEA results
ggplot(top_tcell_pathways, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
  geom_col() +
  scale_fill_gradientn(
    colors = brewer.pal(n = 10, name = "RdYlGn"),
    guide = guide_colorbar(barwidth = 10, barheight = 1),
    name = "Adjusted P-Value"
  ) +
  coord_flip() +
  labs(
    y = "Normalized Enrichment Score (NES)",
    x = "Pathway",
    title = "Top T cell Functional Pathways for CD4 T Cells (TCDD vs. Vehicle)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )


#Histogram
Idents(Tcells) <-"Group"   
samples.list <- unique(Tcells$Group)
clusters <- lapply(samples.list, function(x){
  subset <- subset(Tcells, subset = Group == x)
  dist <- data.frame(table(subset$T_clusters))
  
  return(dist)
})

names(clusters) <- samples.list

#calculate relative freq (fractions) of each cell type
clusters_percent <- lapply(clusters, FUN = function(x){
  summ <- sum(x$Freq)
  x$Freq <- (x$Freq/summ)
  return(x)
})
Idents(object = Tcells) <- "Group"
levels(Tcells)
#making things ggplot-friendly!
clusters_dist <- reshape2::melt(clusters, id.var = "Var1")
colnames(clusters_dist) <- c("T_clusters","variable","value","Group")
clusters_percent_dist <- reshape2::melt(clusters_percent, id.var = "Var1")
colnames(clusters_percent_dist) <- c("T_clusters","variable","value","Group")


ggplot(clusters_percent_dist, aes(fill=T_clusters, y = value, x = Group)) + 
  scale_x_discrete(limits = c("Vehicle","TCDD")) + 
  geom_bar(position = "stack", stat = "identity") + 
  theme_bw()+
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, size = 8, face = "bold", hjust = 1))+
  ylab("count") + xlab("Group") + ggtitle("Relative Cell Types Abundance") +
  scale_fill_manual(values = c("Naive CD4"= "blue",
                               "Naive CD8"="hotpink", 
                               "Naive Cd4&Cd8 DP"="plum",
                               "Treg"="chartreuse1",
                               "Exhausted CD8" ="burlywood",
                               "Exhausted CD4&CD8 DP"="aquamarine",
                               "Memory CD4"="indianred",
                               "Memory CD8"="gold",
                               "Not Tcells" ="thistle"))
#Heatmap