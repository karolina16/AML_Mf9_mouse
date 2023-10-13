library(gridExtra)
library(dittoSeq)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(here)

#### load data
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
data_ecs<- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_4clusters.rds"))

#### quick visual check
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_fibros, dimred = "UMAP_harmony", 
               colour_by=rownames(data_fibros[grep("Col3a1$",rownames(data_fibros))]))


plotReducedDim(data_ecs, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_ecs, dimred = "UMAP_harmony", 
               colour_by="Cspg4")


#### split UMAP by time point
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="time_point") + facet_wrap(~data_fibros$time_point)
ggsave(here("Plots_AML_mouse/", "Clusters_Mesenchymal_TimePoint_Exp1.pdf"), width = 18, height = 6)
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="Lepr") + facet_wrap(~data_fibros$time_point)

# with dittoSeq we get the gray area of the dataset
data_fibros_filt <- data_fibros[,!duplicated(colnames(data_fibros))]
dittoDimPlot(data_fibros_filt, var="time_point", reduction.use = "UMAP_harmony", 
             split.by = c("time_point"), main = "Mesenchymal cells by time point")
dittoDimPlot(data_fibros_filt, var="label", reduction.use = "UMAP_harmony", 
             split.by = c("time_point"),
             main = "Mesenchymal clusters at con, D10 and D60")

data_ecs_filt <- data_ecs[,!duplicated(colnames(data_ecs))]
dittoDimPlot(data_ecs_filt, var="time_point", reduction.use = "UMAP_harmony", 
             split.by = c("time_point"), main = "Endothelial cells by time point")
dittoDimPlot(data_ecs_filt, var="label", reduction.use = "UMAP_harmony", 
             split.by = c("time_point"),
             main = "Endothelial clusters at con, D10 and D60")

#### violin plots for genes from DM
multi_dittoPlot(data_fibros_filt, vars = c("Kitl", "Cxcl12", "Adipoq", "Col3a1",
                                           "Grem1", "Sox5", "Notch3", "Fn1", 
                                           "Dcn", "Gsn", "Vegfc", "Rgs5"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### violin plots for markers from Tikhonova paper
multi_dittoPlot(data_fibros_filt, vars = c("Lepr", "Col1a1", "Wif1", "Lpl",
                                           "Spp1", "Ibsp", "Fbn1", "Igf1", 
                                           "Bglap", "Car3", "Tnn", "Cxcl12"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### violin plots  CAR cells markers from Omatsu paper
multi_dittoPlot(data_fibros_filt, vars = c("Foxc1", "Runx1", "Runx2", "Ebf3","Ebf1","Kitl", "Cxcl12", "Lepr"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### markers for clusters 6 and 7 based on papers
multi_dittoPlot(data_fibros_filt, vars = c("Col3a1", "Rgs5", "Notch3", "Sgip1",
                                           "Mast4","Fbln2", "Cald1", "Myh11",
                                           "Pdgfra", "Ly6a", "Fn1", "Gsn", "Dcn"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### additional CAR cells markers from Omatsu Immunity 2010 paper
multi_dittoPlot(data_fibros_filt, vars = c("Itgav", "Vcam1", "Pdgfra", "Pdgfrb"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 2,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#  Adipo-osteogenic markers from Omatsu Immunity
multi_dittoPlot(data_fibros_filt, vars = c("Cebpa", "Pparg", "Runx2", ), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))


#### markers that change within Lepr+ cells
markers_lepr_pos <- c("Mef2c", "Kcnk2", "Alpl", "Tnc", "Kitl", "Lepr",
                      "Ebf3", "Runx1", "Runx2", "Ebf1", "Col1a1", "Col1a2",
                      "Adipoq", "Wif1", "Lrp4") 
for (i in 1:length(markers_lepr_pos)) {
  gg <- multi_dittoDimPlot(data_fibros_filt, var=markers_lepr_pos[i], reduction.use="UMAP_harmony",
                           ncol = 1, min.color="#878787", max.color="red",
                           split.by="time_point")
  ggsave(filename=paste0(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_DM/Plots_AML_Mes_DM/"), markers_lepr_pos[i], "_Markers_Mes_LeprPos.pdf"), 
         plot=gg, width = 15, height = 6)
  
}

#### markers that change within clusters 6 and 7
markers_6_7 <- c("Notch3", "Rgs5", "Aspn", "Ebf2", "Sox5", "Lum", "Serpinf1",
                 "Clec3b", "Ltbp4", "Ly6a", "Dcn", "Fn1", "Igfbp6", "Mustn1", 
                 "Myh11", "Myo1b", "Myl9", "Acta2") 
for (i in 1:length(markers_6_7)) {
  gg <- multi_dittoDimPlot(data_fibros_filt, var=markers_6_7[i], reduction.use="UMAP_harmony",
                           ncol = 1, min.color="#878787", max.color="red",
                           split.by="time_point")
  ggsave(filename=paste0(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_DM/Plots_AML_Mes_DM/"), markers_6_7[i], "_Markers_Mes_Clusters_6_7.pdf"), 
         plot=gg, width = 15, height = 6)
  
}

#### additional markers for clusters 1:5
multi_dittoPlot(data_fibros_filt, vars = c("Pdzrn4", "Col8a1", "Tnc", "Alpl", "Mef2c", "Spp1",
                                           "Kcnk2"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))


