library(gridExtra)
library(dittoSeq)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(Seurat)
library(scater)
library(here)
library(viridis)


#### load data
data_ecs<- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5clusters.rds")) # old initial object with 4, 5 clusters were used later

data_ecs_filt <- data_ecs[,!duplicated(colnames(data_ecs))]
# remove 3 cells from D60 in cluster 1 for plotting because it affected ditstributions
d60_ecs1 <- colnames(data_ecs_filt[,data_ecs_filt$label == "1" & data_ecs_filt$time_point == "D60"])
data_ecs_filt2 <- data_ecs_filt[,!colnames(data_ecs_filt) %in% d60_ecs1]

#### quick visual check
plotReducedDim(data_ecs, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_ecs, dimred = "UMAP_harmony", 
               colour_by="Spp1")


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

#### plot TFs
tfs_mes <- c("Id3", "Twist2", "Sp7", "Irf1", "Runx2", "Smad4", "Foxg1", "Dlx5")

dittoDotPlot(data_fibros_filt, vars =tfs_mes ,group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 

tfs_ecs <- c("Snai2", "Stat5b", "Klf6", "Id3", "Rara", "Clock", "Mta2")
dittoDotPlot(data_ecs_filt, vars =tfs_ecs ,group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 

#### plot barrier genes by time point
# for ECs
ecs_barrier <- c("Cdh5", "Jam2", "F11r", "Pecam1", "Tjp1", "Esam", "Chd4", 
                 "Igf1r", "Tek", "Ptprm", "Dock1", "Elmo1")
for (i in 1:length(ecs_barrier)) {
  dittoDimPlot(data_ecs_filt, var=ecs_barrier[i], reduction.use = "UMAP_harmony", 
               split.by = c("time_point"),
               main = paste(ecs_barrier[i], " Expression in ECs at con, D10 and D60")) + 
    scale_color_viridis(option = "C")
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(ecs_barrier[i], "_ECs_TimePoint_AML.pdf", sep = "")), width = 14,  height = 7)
  
}

# split by cluster   
for (i in 1:length(ecs_barrier)) {
  dittoPlot(data_ecs_filt2, ecs_barrier[i], group.by = "time_point", 
            split.by = "label",
            split.adjust = list(scales = "free_y"), max = NA)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(ecs_barrier[i], "_ECs_Label_AML.pdf", sep = "")), width = 10,  height = 14)
  
}

# ridge plots
for (i in 1:length(ecs_barrier)) {
  dittoRidgePlot(data_ecs_filt2, ecs_barrier[i], group.by = "time_point", 
                 split.by = "label",
                 split.ncol = 1)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_ECs_EDA/", paste(ecs_barrier[i], "_ECs_Barrier_Ridge_AML.pdf", sep = "")), width = 10,  height = 14)
}

# dot plot with ECs barrier genes
dittoDotPlot(data_ecs_filt, vars =ecs_barrier ,group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 
