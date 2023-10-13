library(gridExtra)
library(dittoSeq)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(here)

# load the data
data_sce <- readRDS(here("R_objects_AML_mouse/AML_mouse_Exp1_sce_clusters.rds"))
data_mes <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
data_ecs <- readRDS(here("R_objects_AML_mouse/AML_mouse_Exp1_ECs_sce.rds"))

# load and sort markersfor healthy fibros to inspect
markers_fibros_C4 <- read.csv(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/Markers_AML_Stromal/Markers_Stromal_Harmonized_C4.csv"))

markers_ordered_C4 <- markers_fibros_C4[order(markers_fibros_C4$other.detected, 
                                      decreasing=FALSE),]
markers_fibros_C4_filt <- markers_fibros_C4[markers_fibros_C4$self.detected > 0.7 
                                            & markers_fibros_C4$other.detected < 0.3,]


write.csv(markers_fibros_C4_filt,
          paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/Markers_AML_Stromal///"), 'Markers_Stromal_Harmonized_C4_filtered', '.csv',
                sep = ""))

# select markers per cluster to plot
markers <- c("Cav1", "Cd36", "Cdh13", "Cdh5", "Lmo2", "Hbb-bs", "Nr4a1",
                    "Col3a1", "Gsn", "Lrp1", "Dcn", "Fn1", "Cd79a",
                    "Cd79b", "Igfbp5", "Cxcl14", "Serping1", "Serpine2",
             "Lepr", "Adipoq", "Kitl", "Mgp", "Plvap", "Cldn5", "Cxcl12")

markers_mes <- c("Alpl", "Wif1", "Lrp4", "Ncam1", "Spp1", "Kcnk2", "Mef2c",
                 "Olfml3", "Jun", "Junb", "Fos", "Adipoq", "Esm1", "Ccn1",
                 "Zfp36", "Gsn", "Col3a1", "Fn1", "Igfbp6", "Kalrn", "Thsd4",
                 "Bmp6")

markers_ecs <- c("Cav1", "Cxcl12", "Cdh13", "Kitl", "Dach1", "Sparcl1", "Mecom",
                 "S100a6","Tfpi", "Selenop", "Igfbp4", "Mafb", "Plvap", "Ctsb",
                 "Maf", "Lrg1", "Stab2")
 
# plots with dittoSeq
# weird duplicated rows:(, remove
data_sce_filt <- data_ecs[,!duplicated(colnames(data_ecs))]
dittoDotPlot(data_sce_filt, vars = markers_ecs, group.by = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 

ggsave(here("Plots_AML_mouse/", "Markers_ECs_Exp1_Clusters1_6.pdf"), width = 15, height = 6)

# visualise markers per time point
plotUMAP(data_mes, colour_by="label") + facet_wrap(~data_mes$time_point)
ggsave(here("Plots_AML_mouse/", "Clusters_Mesenchymal_TimePoint_Exp1.pdf"), width = 18, height = 6)



dittoPlotVarsAcrossGroups(data, s_common, group.by = "label2",
                          main = "S phase genes")

multi_dittoDimPlot(data_fibros, prc)
