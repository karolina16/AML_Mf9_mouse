library(gridExtra)
library(dittoSeq)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(here)
library(viridis)

#### load data
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
data_ecs<- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5clusters.rds")) # old initial object with 4, 5 clusters were used later

data_fibros_filt <- data_fibros[,!duplicated(colnames(data_fibros))]
data_ecs_filt <- data_ecs[,!duplicated(colnames(data_ecs))]


#### quick visual check
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_fibros, dimred = "UMAP_harmony", 
               colour_by="Cdh2")


plotReducedDim(data_ecs, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_ecs, dimred = "UMAP_harmony", 
               colour_by="Myh9")


#### split UMAP by time point
# this was not used
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="time_point") + facet_wrap(~data_fibros$time_point)
ggsave(here("Plots_AML_mouse/", "Clusters_Mesenchymal_TimePoint_Exp1.pdf"), width = 18, height = 6)
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="Lepr") + facet_wrap(~data_fibros$time_point)

# with dittoSeq we get the gray area of the dataset and this was used
dittoDimPlot(data_fibros_filt, var="time_point", reduction.use = "UMAP_harmony", 
             split.by = c("time_point"), main = "Mesenchymal cells by time point")
dittoDimPlot(data_fibros_filt, var="label", reduction.use = "UMAP_harmony", 
             split.by = c("time_point"),
             main = "Mesenchymal clusters at con, D10 and D60")

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
multi_dittoPlot(data_fibros_filt, vars = c("Cebpa", "Pparg", "Runx2", "Sp7",
                                           "Col1a1", "Bglap"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
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

#### plot fusion gene target genes in fibros
fusion_gene_targets <- c("Meis1", "Hoxa9", "Eya1")
# plot fusion genes
for (i in 1:length(fusion_gene_targets)) {
  plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by=fusion_gene_targets[i]) 
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA/", paste(fusion_gene_targets[i], "_Mes_AML.pdf", sep = "")), width = 7,  height = 7)
}

multi_dittoPlot(data_fibros_filt, vars = fusion_gene_targets, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

multi_dittoPlot(data_ecs_filt, vars = fusion_gene_targets, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

for (i in 1:length(fusion_gene_targets)) {
  dittoPlot(data_fibros_filt, fusion_gene_targets[i], group.by = "time_point", 
            split.by = "label",
            split.adjust = list(scales = "free_y"), max = NA)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA//", paste(fusion_gene_targets[i], "_Mes_FusionTargets_Label_AML_Violin.pdf", sep = "")), width = 10,  height = 14)
  
}

for (i in 1:length(fusion_gene_targets)) {
  dittoPlot(data_fibros_filt, fusion_gene_targets[i], group.by = "label", 
            split.by = "time_point",
            split.adjust = list(scales = "free_y"), max = NA)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA//", paste(fusion_gene_targets[i], "_Mes_FusionTargets_TimePoint_AML_Violin.pdf", sep = "")), width = 10,  height = 7)
  
}


# split by time_point
for (i in 1:length(fusion_gene_targets)) {
  dittoDimPlot(data_ecs_filt, var=fusion_gene_targets[i], reduction.use = "UMAP_harmony", 
               split.by = c("time_point"),
               main = paste(fusion_gene_targets[i], " Expression in ECs at con, D10 and D60")) + 
    scale_color_viridis(option = "C")
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(fusion_gene_targets[i], "_ECs_TimePoint_AML.pdf", sep = "")), width = 14,  height = 7)
  
}

for (i in 1:length(fusion_gene_targets)) {
  dittoDimPlot(data_fibros_filt, var=fusion_gene_targets[i], reduction.use = "UMAP_harmony", 
               split.by = c("time_point"),
               main = paste(fusion_gene_targets[i], " Expression in Mes cells at con, D10 and D60")) + 
    scale_color_viridis(option = "C")
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA/", paste(fusion_gene_targets[i], "_Mes_TimePoint_AML.pdf", sep = "")), width = 14,  height = 7)
  
}


#### plot adipo and osteogenic markers
# adipogenic markers
markers_adipo <- c("Mgp", "Gpx3", "Tmem176b", "Scp2", "Lpl", "Fstl1", 
                   "Pdzrn4", "Slc5a3", "Angpt1")

# Bglap and Car3 marekers of mature osteoblasts
markers_osteo <-c("Wif1", "Col8a1", "Kcnk2", "Limch1", "Palld", "Tnfrsf19", 
                  "Spp1", "Col1a1", "Col13a1", "Mmp13", "Ifitm5", "Serpine2", 
                  "Mef2c", "Aqp1", "Igfbp5", "Bglap", "Car3", "Col11a2")

multi_dittoPlot(data_fibros_filt, vars = markers_adipo, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### plot CAF markers
markers_CAF <- c("Fap", "S100a4", "Pdgfra", "Pdgfrb", "Vim", "Cav1", "Postn", 
                 "Tagln", "Itga11", "Col11a1", "Mfap5", "Cd34", "Ly6a", "Cd44", "Thy1")

multi_dittoPlot(data_fibros_filt, vars = markers_CAF, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

#### plot barrier integrity markers in ECs
markers_bar_int <- c("Dock1", "Elmo1", "Nfe2l2", "Eng", "Kdr", "Cdh5", "Tjp1", 
                     "Myl2")

nrf2_targets <- c("Ulk1", "Atg5", "Cul3", "Gclc", "Gclm", "Keap1", "S1pr1",
                     "Hmox1","Cat", "P4hb", "Sema6a", "Il6")

dittoDotPlot(data_ecs_filt, vars = markers_bar_int, group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 

dittoDotPlot(data_ecs_filt, vars = nrf2_targets, group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 


multi_dittoPlot(data_ecs_filt, vars = markers_bar_int, 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

