library(SingleCellExperiment)
library(strex)
library(purrr)
library(visNetwork)
library(dplyr)
library(scater)
library(scran)
library(batchelor)
library(bluster)
library(tidyr)
library(here)
library(dittoSeq)
library(bluster)

set.seed(1000)

#### load the data
data_bs <- readRDS(here("Data_AML_mouse/Data_Basel/transfer_184642_files_3dc8eea7/sce_withTSNE_29092023.rds"))
data_bs

table(data_bs$final_annotation)

#### subset the data for the interactions analysis
data_bs_sub <- data_bs[,data_bs$final_annotation %in% c("MPPStemcells")]

saveRDS(data_bs_sub, here("Data_AML_mouse/Data_Basel/MPPStemcells_Basel.rds"))
data_mpp <- readRDS(here("Data_AML_mouse/Data_Basel/MPPStemcells_Basel.rds"))


#### explore cell annotation
table(data_bs$final_annotation)
table(data_bs$Type)
table(data_bs$Batch)
head(colData(data_bs))
table(data_bs$Day)
table(data_bs$Mouse)

plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "final_annotation")
plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Batch")
plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Day")
plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Mouse")

table(cell_type=data_bs$final_annotation, mouse=data_bs$Mouse)

#### check the MLL-AF9 expression and its target genes
fusion_gene_targets <- c("Meis1", "Hoxa9", "Eya1")

plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "MLL-AF9")

multi_dittoPlot(data_bs, vars = "MLL-AF9", 
                group.by = "Mouse", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 1,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

multi_dittoPlot(data_bs, vars = fusion_gene_targets, 
                group.by = "Day", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### prepare data for interactions analysis
# Re-process the data
dec <- modelGeneVar(data_mpp)

# Order by most interesting genes for inspection
dec[order(dec$bio, decreasing=TRUE),] 

# Select highly variable genes based on variance modeling
hvgs <- getTopHVGs(dec, prop=0.2)

length(hvgs)

# Dimensionality reduction
data_mpp <-denoisePCA(data_mpp, technical=dec, subset.row=hvgs)
data_mpp <- runUMAP(data_mpp, dimred = 'PCA', external_neighbors=TRUE)

plotUMAP(data_mpp, colour_by="Mouse")
plotUMAP(data_mpp, colour_by="Day")
plotUMAP(data_mpp, colour_by="Batch")


#harmonization needed
# harmonized_pcs <- HarmonyMatrix(
#   data_mat  = reducedDims(data_mpp)[["PCA"]],
#   meta_data = as.data.frame(colData(data_mpp)), # Dataframe with information for each cell (row)
#   vars_use  = "exp_nb", # Column in meta_data that defines dataset for each cell
#   do_pca    = FALSE      # Since we are providing PCs, do not run PCA
# )
# # rerun dimred and clusering
# reducedDims(data_mpp)[["harmony_pca"]] <- harmonized_pcs
# data_mpp <- runUMAP(data_mpp, dimred = 'harmony_pca', external_neighbors=TRUE, name = "UMAP_harmony")
# 
# # Visualization.
# plotReducedDim(data_mpp, dimred = "UMAP_harmony", colour_by="cond", text_by="label")


#### run  clustering 
kgraph.clusters <- clusterCells(data_mpp, use.dimred="UMAP",
                                BLUSPARAM=TwoStepParam(
                                  first=KmeansParam(centers=500),
                                  second=NNGraphParam(k=70)
                                )
)
table(kgraph.clusters)
colLabels(data_mpp) <- factor(I(kgraph.clusters))
plotReducedDim(data_mpp, dimred = "UMAP", colour_by="label", text_by="label")


# kmeans better for ECs not used in the end, graph clusters were done
# set.seed(150)
# clust.kmeans <- clusterCells(data_mpp, use.dimred="harmony_pca", 
#                              BLUSPARAM=KmeansParam(centers=4))
# table(clust.kmeans)
# 
# colLabels(data_mpp) <- factor(I(clust.kmeans))
# plotReducedDim(data_mpp, dimred = "UMAP_harmony", colour_by="label", text_by="label")


# inspect the clusters in more detail
tab <- table(cluster=data_mpp$label, cond=data_mpp$Mouse)
tab
write.csv(tab, here("Results_AML_mouse/Results_AML_mouse_EDA/Harmonized_Clusters_vs_TimePoint_AML_ECs_mouse.csv"))


tab2 <- table(cluster=data_mpp$label, day=data_mpp$Day)
tab2
write.csv(tab2, here("Results_AML_mouse/Results_AML_mouse_EDA/Harmonized_Clusters_vs_Exp_TimePoint_AML_ECs_mouse.csv"))

tab3 <- table(Mouse=data_mpp$Mouse, day=data_mpp$Day)
tab3

#### run markers
# Find markers for each cluster
markers <- scoreMarkers(data_mpp, lfc=0.25)

# Save markers for all clusters
for (i in names(markers)) {
  markers_ordered <- markers[[i]][order(markers[[i]]$mean.AUC, 
                                        decreasing=TRUE),]
  write.csv(markers_ordered[markers_ordered$mean.AUC>0.5,],
            paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_MPPStemCells/"), 'Markers_MPPStemCells_C', i, '.csv',
                  sep = ""))
  
}
