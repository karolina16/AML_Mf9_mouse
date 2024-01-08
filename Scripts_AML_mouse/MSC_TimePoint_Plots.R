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

data_fibros_filt <- data_fibros[,!duplicated(colnames(data_fibros))]

#### quick visual check
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_fibros, dimred = "UMAP_harmony", 
               colour_by="Spp1")


#### plot niche markers over time points and clusters
niche_markers <- c("Cxcl12", "Adipoq", "Cdh2", "Kitl", "Angpt1", "Angptl1")
com_markers <- c("Fgf7", "Fgfr1", "Fgfr2", "Gja1", "Igf1", "Itgav", "Ptprm")

#### split violins by cluster
# plot niche markers
for (i in 1:length(niche_markers)) {
  dittoPlot(data_fibros_filt, niche_markers[i], group.by = "time_point", 
            split.by = "label",
            split.adjust = list(scales = "free_y"), max = NA)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA/", paste(niche_markers[i], "_Mes_Label_AML.pdf", sep = "")), width = 10,  height = 14)
  
}

# plot cell-cell communications genes
for (i in 1:length(com_markers)) {
  dittoPlot(data_fibros_filt, com_markers[i], group.by = "time_point", 
            split.by = "label",
            split.adjust = list(scales = "free_y"), max = NA)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA/", paste(com_markers[i], "_Mes_Label_AML.pdf", sep = "")), width = 10,  height = 14)
  
}


#### ridge plots
# plot niche markers
for (i in 1:length(niche_markers)) {
  dittoRidgePlot(data_fibros_filt, niche_markers[i], group.by = "time_point", 
                 split.by = "label",
                 split.ncol = 1)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA/", paste(niche_markers[i], "_Mes_Barrier_Ridge_AML.pdf", sep = "")), width = 10,  height = 14)
}

# plot cell-cell communications genes
for (i in 1:length(com_markers)) {
  dittoRidgePlot(data_fibros_filt, com_markers[i], group.by = "time_point", 
                 split.by = "label",
                 split.ncol = 1)
  ggsave(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Ex1_Exp2_Fibros_EDA/", paste(com_markers[i], "_Mes_Barrier_Ridge_AML.pdf", sep = "")), width = 10,  height = 14)
}


# dot plot with ECs barrier genes
dittoDotPlot(data_fibros_filt, vars = com_markers, group.by  = "label",
             scale = FALSE, size = 12) + theme(text = element_text(size = 25)) 
