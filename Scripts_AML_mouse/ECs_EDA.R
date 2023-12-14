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
library(forcats)
library(viridis)
library(dittoSeq)

set.seed(1000)

# load the EC object
data_ec <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_ECs_Exp1_Exp2.rds"))
# This was run after initial ECs prosessing
# load the processed ECs object for markers between healthy and AML ECs
data_ec <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5Clusters.rds"))
# add meta data
# data_ec$state <- fct_collapse(data_ec$label,
#                                     ECs_healthy_d10 = c("1", "3"),
#                                     ECs_d60 = c("2", "4"))

#### quick visual check
plotUMAP(data_ec, colour_by="label")

#### Re-process the data
dec <- modelGeneVar(data_ec)

# Order by most interesting genes for inspection
dec[order(dec$bio, decreasing=TRUE),] 

# Select highly variable genes based on variance modeling
hvgs <- getTopHVGs(dec, prop=0.2)

length(hvgs)

# Dimensionality reduction
data_ec <-denoisePCA(data_ec, technical=dec, subset.row=hvgs)
data_ec <- runUMAP(data_ec, dimred = 'PCA', external_neighbors=TRUE)

plotUMAP(data_ec, colour_by="time_point")
plotUMAP(data_ec, colour_by="exp_nb")
plotUMAP(data_ec, colour_by="Cav1")


#harmonization needed
harmonized_pcs <- HarmonyMatrix(
  data_mat  = reducedDims(data_ec)[["PCA"]],
  meta_data = as.data.frame(colData(data_ec)), # Dataframe with information for each cell (row)
  vars_use  = "exp_nb", # Column in meta_data that defines dataset for each cell
  do_pca    = FALSE      # Since we are providing PCs, do not run PCA
)
# rerun dimred and clusering
reducedDims(data_ec)[["harmony_pca"]] <- harmonized_pcs
data_ec <- runUMAP(data_ec, dimred = 'harmony_pca', external_neighbors=TRUE, name = "UMAP_harmony")

# Visualization.
plotReducedDim(data_ec, dimred = "UMAP_harmony", colour_by="cond", text_by="label")


#### run  clustering 
kgraph.clusters <- clusterCells(data_ec, use.dimred="harmony_pca",
                                BLUSPARAM=TwoStepParam(
                                  first=KmeansParam(centers=1000),
                                  second=NNGraphParam(k=20)
                                )
)
table(kgraph.clusters)
colLabels(data_ec) <- factor(I(kgraph.clusters))
plotReducedDim(data_ec, dimred = "UMAP_harmony", colour_by="label", text_by="label")


# kmeans better for ECs not used in the end, graph clusters were done
# set.seed(150)
# clust.kmeans <- clusterCells(data_ec, use.dimred="harmony_pca", 
#                              BLUSPARAM=KmeansParam(centers=4))
# table(clust.kmeans)
# 
# colLabels(data_ec) <- factor(I(clust.kmeans))
# plotReducedDim(data_ec, dimred = "UMAP_harmony", colour_by="label", text_by="label")


# inspect the clusters in more detail
tab <- table(cluster=data_ec$label, cond=data_ec$time_point)
tab
write.csv(tab, here("Results_AML_mouse/Results_AML_mouse_EDA/Harmonized_Clusters_vs_TimePoint_AML_ECs_mouse.csv"))


tab2 <- table(cluster=data_ec$label, cond=data_ec$exp_time_point)
tab2
write.csv(tab2, here("Results_AML_mouse/Results_AML_mouse_EDA/Harmonized_Clusters_vs_Exp_TimePoint_AML_ECs_mouse.csv"))

# plot distribution of cell among clusters
ggplot(data.frame(tab), aes(fill=cond, y=Freq, x=cluster)) + 
  geom_bar(position="stack", stat="identity")


#### run markers
# split into sinusoidal and arterial ECs did not work
#data_ec_sub <- data_ec[,data_ec$label %in% c(2,4)]
# additional analysis was done betwen healthy and diseased ECs
# par block was not used before:(
data_ec$state <- fct_collapse(data_ec$label,
                                    ECs_healthy_d10 = c("1", "5"),
                                    ECs_d60 = c("2", "3", "4"))

markers <- scoreMarkers(data_ec, lfc=0.25, block=data_ec$exp_nb)
markers <- scoreMarkers(data_ec, groups=data_ec$state, lfc=0.25, block=data_ec$exp_nb)

# Save markers for all clusters
for (i in names(markers)) {
  markers_ordered <- markers[[i]][order(markers[[i]]$mean.AUC, 
                                        decreasing=TRUE),]
  write.csv(markers_ordered[markers_ordered$mean.AUC>0.5,],
            paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/Markers_AML_ECs//"), 'Markers_ECs_AML_Con_State_', i, '.csv',
                  sep = ""))
  }

plotReducedDim(data_ec, dimred = "UMAP_harmony", 
               colour_by=rownames(data_ec[grep("S1p$",rownames(data_ec))]))

#### plot selecteed markers
# plots with dittoSeq
markers_ec <- c("S100a6", "Crip1", "Dach1", "Ccdc85a", "Meox2", "Bmp6", "Mrc1",
                 "Stab2", "Lrg1", "Tfpi", "Mafb", "Igfbp4", 
                "Plvap", "Cav1", "Cdh13")

markers_ec2 <- c("Tfpi", "Mafb", "Adamts5", "Stab2", "Lrg1", "Sparcl1", "Kitl",
                 "Clic5", "Dach1", "Mecom", "S100a6", "Mrc1", "Igfbp4", "Abcc9",
                 "Il6st", "Flt4", "Nr2f2", "Sema3a", "Adam12", "Plvap", "Cav1",
                 "Cdh13", "Cxcl12")

markers_ec3 <- c("Kdr", "Myh9", "Igfbp7", "Col4a1", "Myo10", "Timp3", "Pdgfb",
                 "Macf1", "Esam", "Itga1", "Zeb1", "Epas1", "Shank3", "Tcf12",
                 "Asap1", "Ctnna1", "Ptprb", "Ctnnb1", "Mef2a", "Tgfbr2", "Tgfbr3")

markers_conD10_vs_d60 <-c("Hsp90ab1", "Hspa8", "Shank3", "Myh9", "Eng",
                          "Elmo1", "Dock1", "Nfe2l2", "Slfn5", "Asap1")

markers_bar_int <- c("Dock1", "Elmo1", "Nfe2l2", "Eng", "Kdr", "Cdh5", "Tjp1")

markers_tum_sup <- c("Hsp90ab1", "Myh9", "Slfn5")

fusion_gene_targets <- c("Meis1", "Hoxa9", "Eya1")

# weird duplicated rows:(, remove
data_ec_filt <- data_ec[,!duplicated(colnames(data_ec))]
dittoDotPlot(data_ec_filt, vars = markers_bar_int, group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 

multi_dittoPlot(data_ec_filt, vars = markers_tum_sup, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

multi_dittoPlot(data_ec_filt, vars = fusion_gene_targets, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

# plot fusion genes
for (i in 1:length(fusion_gene_targets)) {
  plotReducedDim(data_ec, dimred = "UMAP_harmony", colour_by=fusion_gene_targets[i]) 
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(fusion_gene_targets[i], "_ECs_AML.pdf", sep = "")), width = 7,  height = 7)
}

# split by time_point
for (i in 1:length(markers_conD10_vs_d60)) {
  dittoDimPlot(data_ec_filt, var=markers_conD10_vs_d60[i], reduction.use = "UMAP_harmony", 
               split.by = c("time_point"),
               main = paste(markers_conD10_vs_d60[i], " Expression in ECs at con, D10 and D60")) + 
    scale_color_viridis(option = "C")
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(markers_conD10_vs_d60[i], "_ECs_TimePoint_AML.pdf", sep = "")), width = 14,  height = 7)
  
}

for (i in 1:length(markers_conD10_vs_d60)) {
  dittoDimPlot(data_ec_filt, var=markers_conD10_vs_d60[i], reduction.use = "UMAP_harmony", 
               split.by = c("time_point"),
               main = paste(markers_conD10_vs_d60[i], " Expression in ECs at con, D10 and D60")) + 
    scale_color_viridis(option = "C")
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(markers_conD10_vs_d60[i], "_ECs_TimePoint_AML.pdf", sep = "")), width = 14,  height = 7)
  
}


# save object
saveRDS(data_ec,  here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/", "ECs_AML_mouse_Con_Adult1_2_harmonized_5Clusters.rds"))
