library(destiny)
library(ggbeeswarm)
library(ggthemes)
library(SingleCellExperiment)
library(scran)
library(scater)
library(grid)
library(gridExtra)
library(ggplot2)
library(here)


# diffusion map will be run only for fibros

#### load data
data <- readRDS("./R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds")

# subset to use only Lepr+ cells or clusters 6 and 7
data <- data[,data$label %in% c(6:7)]
data$label <- droplevels(as.factor(data$label))
#### explore markers
plotReducedDim(data, dimred = "UMAP_harmony", 
               colour_by=rownames(data[grep("Egr1$",rownames(data))]))

#### run DM and DPT
set.seed(1234)
dm <- DiffusionMap(data, n_pcs = 50, rotate = TRUE)
#saveRDS(dm, "./R_objects_AML_mouse/AML_Exp1_Exp2_DM/Mes_DM_AML_Exp1_Exp2.rds")
dpt <- DPT(dm)

### Plot DM 
plot.DiffusionMap(dm)
plot.DiffusionMap(dm, pch=20, dims = c(2,3), col_by="label")
plot.DiffusionMap(dm, pch=20, dims = c(3,2,1), col_by="label")
plot.DiffusionMap(dm, pch=20, dims = c(3,2,1), col=as.factor(data$time_point),
                  draw_legend = T)
# c(3,1,2) for clusters 6 and 7
# c(2,3,1) used for LeprPos clusters
#plot.DiffusionMap(dm, pch=20, dims = c(4,2,1), col_by="cond") does not work dims used for all clusters

# colour by gene expression
genes_to_plot <- c("Col3a1", "Rgs5", "Sox5", "Notch3", "Mast4", "Gsn", "Dcn", 
                   "Fn1", "Igfbp6", "Ly6a", "Cxcl12", "Igfbp5", "Lrp4", "Wif1",
                   "Kitl", "Apoe", "Cxcl14", "Tgfbr2", "Tgfbr3", "Lepr", "Bmp6",
                   "Clec2d", "Chrdl1", "Angpt1", "Adipoq", "Vegfc", "Grem1", 
                   "Col4a2", "Stat5b", "S100a6", "Col1a1", "Col1a2")

genes_to_plot_lpr <- c("Cxcl12", "Igfbp5", "Lrp4", "Wif1","Kitl", "Apoe", "Cxcl14",
                       "Tgfbr2", "Tgfbr3", "Lepr", "Bmp6", "Clec2d", "Chrdl1", 
                       "Angpt1", "Adipoq", "Vegfc", "Grem1", "Col4a2", "Stat5b", 
                       "S100a6", "Col1a1", "Col1a2", "Ebf1", "Runx1", "Runx2",
                       "Ebf3", "Pdzrn4", "Col8a1", "Tnc", "Alpl", "Mef2c", "Spp1",
                       "Kcnk2")

genes_to_plot_c6_7 <- c("Fbln2", "Col3a1", "Gsn", "Fn1", "Igfbp6", "Dcn", "Ly6a",
                        "Ltbp4", "Clec3b", "Lum", "Serpinf1", "Sox5", "Ebf2", 
                        "Aspn", "Rgs5", "Notch3", "Mustin1", "Myh11", "Tagln")
  
for (i in 1:length(genes_to_plot_c6_7)) {
  pdf(file=paste0(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_DM/Plots_AML_Mes_DM/"), 
               genes_to_plot_c6_7[i], "_Mes_C6_7_DM_AML.pdf"), width = 10, height = 10)
  plot.DiffusionMap(dm, pch=20, dims = c(3,2,1), col_by=genes_to_plot_c6_7[i])
  dev.off()
  }

### Plot pseudotime
plot(eigenvalues(dm)[1:100])
plot.DPT(dpt, dcs = c(1,4))
plot.DPT(dpt, root = 2, paths_to = c(1,3), col_by = 'branch', dcs = c(1,4))
plot.DPT(dpt, root = 2, paths_to = c(1,3), col_by = 'label', dcs = c(1,4))



# not sure if this is useful
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 3],
                  Label = data$label)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Label)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

# visualise DC1 as proxy of pseudotime
data$diffusion_component1 <- rank(eigenvectors(dm)[,1])   # rank cells by their DC
ggplot(as.data.frame(colData(data)), 
       aes(x = diffusion_component1, 
           y = label, colour = label)) + geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +  scale_color_tableau() +
  xlab("Diffusion component 1 (DC1)") + ylab("Cluster nb") +
  ggtitle("Cells ordered by DC1")

# order cells according to pseudotime
data$pseudotime_dpt <- rank(dpt$dpt) 
ggplot(as.data.frame(colData(data)), 
       aes(x = pseudotime_dpt, 
           y = label, colour = cond)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion map pseudotime (dpt)") +
  ylab("Cluster nb") +
  ggtitle("Mesenchymal Cells ordered by diffusion map pseudotime")