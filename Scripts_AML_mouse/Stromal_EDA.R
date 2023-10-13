library(SingleCellExperiment)
library(strex)
library(purrr)
library(visNetwork)
library(dplyr)
library(scater)
library(scran)
library(bluster)
library(tidyr)
library(here)
library(harmony)
library(dittoSeq)

set.seed(1000)

# load the harmonized data and subset stromal cells 
data_str <- readRDS(here("R_objects_AML_mouse/AML_fibros_Exp1_Exp2.rds"))

plotUMAP(data_str, colour_by="label")
table(cluster=data_str$label, batch=data_str$exp_time_point)

dec <- modelGeneVar(data_str)

# Order by most interesting genes for inspection
dec[order(dec$bio, decreasing=TRUE),] 

# Select highly variable genes based on variance modeling
hvgs <- getTopHVGs(dec, prop=0.2)

length(hvgs)

# Dimensionality reduction
#data_str <- runPCA(data_str, ncomponents=25, subset_row=hvgs)
data_str <-denoisePCA(data_str, technical=dec, subset.row=hvgs)
#data_str <- runPCA(data_str, ncomponents=25, subset_row=hvgs)
data_str <- runUMAP(data_str, dimred = 'PCA', external_neighbors=TRUE)

plotUMAP(data_str, colour_by="time_point")
plotUMAP(data_str, colour_by="exp_time_point")

# harmonization not needed
harmonized_pcs <- HarmonyMatrix(
  data_mat  = reducedDims(data_str)[["PCA"]],
  meta_data = as.data.frame(colData(data_str)), # Dataframe with information for each cell (row)
  vars_use  = "exp_nb", # Column in meta_data that defines dataset for each cell
  do_pca    = FALSE      # Since we are providing PCs, do not run PCA
)
# rerun dimred and clusering
reducedDims(data_str)[["harmony_pca"]] <- harmonized_pcs
data_str <- runUMAP(data_str, dimred = 'harmony_pca', external_neighbors=TRUE, name = "UMAP_harmony")

# Visualization.
plotReducedDim(data_str, dimred = "UMAP_harmony", colour_by="exp_time_point", text_by="label")


#### run  clustering 
kgraph.clusters <- clusterCells(data_str, use.dimred="harmony_pca",
                                BLUSPARAM=TwoStepParam(
                                  first=KmeansParam(centers=1000),
                                  second=NNGraphParam(k=15)
                                )
)

table(kgraph.clusters)
colLabels(data_str) <- factor(I(kgraph.clusters))

plotReducedDim(data_str, dimred = "UMAP_harmony", colour_by="label", text_by="label")
#plotUMAP(data_str, colour_by="label", text_by="label")

# inspect the clusters in more detail
tab <- table(cluster=data_str$label, cond=data_str$time_point)
tab
write.csv(tab, here("Results_AML_mouse/Results_AML_mouse_EDA/Clusters_vs_TimePoint_AML_Mes_mouse.csv"))


tab2 <- table(cluster=data_str$label, cond=data_str$exp_time_point)
tab2
write.csv(tab2, here("Results_AML_mouse/Results_AML_mouse_EDA/Clusters_vs_Exp_TimePoint_AML_Mes_mouse.csv"))


# plot distribution of cell among clusters
ggplot(data.frame(tab), aes(fill=cond, y=Freq, x=cluster)) + 
  geom_bar(position="stack", stat="identity")



# run markers
markers <- scoreMarkers(data_str, lfc=0.25)

# Save markers for all clusters
for (i in names(markers)) {
  markers_ordered <- markers[[i]][order(markers[[i]]$mean.AUC, 
                                        decreasing=TRUE),]
  write.csv(markers_ordered[markers_ordered$mean.AUC>0.5,],
            paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/Markers_AML_Stromal/"), 
                  'Markers_Stromal_Harmonized_C', i, '.csv',
                  sep = ""))
  }

plotReducedDim(data_str, dimred = "UMAP_harmony", 
               colour_by=rownames(data_str[grep("Col3a1$",rownames(data_str))]))

#### plot selecteed markers
# plots with dittoSeq
markers_mes <- c("Bmp6", "Pdzrn4", "Negr1", "Tnc", "Alpl", "Lrp4", "Wif1", "Mef2c",
                 "Apoe", "Col8a1", "Chrdl1","Angpt1", "Adipoq", "Grem1", "Clu", 
                 "Fos", "Jun", "Egr1", "Fosb", "Zfp36", "Sox5", "Rgs5", "Col3a1",
                 "Sgip1", "Notch3", "Gsn", "Dcn", "Fn1", "Igfbp6", "Ly6a", "Tnxb",
                 "Ltbp4", "Cxcl12", "Kitl", "Col1a1", "Col1a2", "Igfbp5", "Pdgfrb")
# weird duplicated rows:(, remove
data_sce_filt <- data_str[,!duplicated(colnames(data_str))]
dittoDotPlot(data_sce_filt, vars = markers_mes, group.by = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 
# save object
saveRDS(data_str,  here("R_objects_AML_mouse/AML_Exp1_Exp2/", "Mes_AML_Exp1_2_harmonized_7clusters.rds"))
