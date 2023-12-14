library(Rmagic)
library(ggplot2)
library(here)
# library(reticulate)
library(scran)
library(scater)
library(viridis)
library(cowplot)
library(rgl)

# import_from_path("magic", path = "/Users/KZ/Library/r-miniconda-arm64/envs/r-reticulate/lib/python3.9/site-packages")
# import_from_path("numpy", path = "/Users/KZ/Library/r-miniconda-arm64/envs/r-reticulate/lib/python3.9/site-packages")

#### load data
data_ecs <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5Clusters.rds"))
data_ecs
plotReducedDim(data_ecs, colour_by="Dock1", dimred = "UMAP_harmony")
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))


# run the analysis for ECs
#### prepare data for magic
# keep genes expressed in more than 10 cells
keep_rows <- rowSums(counts(data_ecs) > 0) > 10
table(keep_rows)
data_ecs <- data_ecs[keep_rows,]
data_ecs
data_ecs <- logNormCounts(data_ecs)
# subset the healthy or AML clusters
data_ecs_aml <- data_ecs[,data_ecs$label %in% c(2:4)]
data_ecs_healthy <- data_ecs[,data_ecs$label %in% c(1,5)]

# w/o magic
ggplot(as.data.frame(t(logcounts(data_ecs_healthy)))) +
  geom_point(aes(Dock1, Elmo1, color=Nfe2l2)) +
  scale_color_viridis(option="B")

#### run magic
MAGIC_data_healthy <- magic(t(logcounts(data_ecs_healthy)), genes=c("Dock1", "Tjp1", "Nfe2l2"), knn = 5)
MAGIC_data_aml <- magic(t(logcounts(data_ecs_aml)), genes=c("Dock1", "Tjp1", "Nfe2l2"), knn = 5)

#### plot the results
p_healthy <- ggplot(MAGIC_data_healthy) +
  geom_point(aes(x=Dock1, y=Tjp1, color=Nfe2l2)) + scale_color_viridis(option="B", limits = c(min(MAGIC_data_aml$result$Nfe2l2), max(MAGIC_data_healthy$result$Nfe2l2))) +
  coord_cartesian(xlim = c(min(MAGIC_data_aml$result$Dock1), max(MAGIC_data_healthy$result$Dock1)), 
                  ylim = c(min(MAGIC_data_aml$result$Tjp1), max(MAGIC_data_healthy$result$Tjp1))) +
  ggtitle("Clusters 1 and 5 healthy ECs")
  

p_aml <- ggplot(MAGIC_data_aml) +
  geom_point(aes(x=Dock1, y=Tjp1, color=Nfe2l2)) + scale_color_viridis(option="B", limits = c(min(MAGIC_data_aml$result$Nfe2l2), max(MAGIC_data_healthy$result$Nfe2l2))) +
  coord_cartesian(xlim = c(min(MAGIC_data_aml$result$Dock1), max(MAGIC_data_healthy$result$Dock1)), 
                  ylim = c(min(MAGIC_data_aml$result$Tjp1), max(MAGIC_data_healthy$result$Tjp1))) +
  ggtitle("Clusters 2:4 D10 and D60")

plot_grid(p_healthy,p_aml) 

# run the analysis for mesenchymal cells
# not much interesting found, look for good genes to plot
#### prepare data for magic
# keep genes expressed in more than 10 cells
keep_rows_mes <- rowSums(counts(data_fibros) > 0) > 10
table(keep_rows_mes)
data_fibros <- data_fibros[keep_rows_mes,]
data_fibros
data_fibros <- logNormCounts(data_fibros)
# subset the healthy or AML clusters
data_fibros_aml <- data_fibros[,data_fibros$label %in% c(1:3,5)]
data_fibros_healthy <- data_fibros[,data_fibros$label %in% c(4)]

# w/o magic
ggplot(as.data.frame(t(logcounts(data_fibros_aml)))) +
  geom_point(aes(Col1a1, Postn, color=Kitl)) +
  scale_color_viridis(option="B")

#### run magic
MAGIC_fibros_healthy <- magic(t(logcounts(data_fibros_healthy)), genes=c("Col1a1", "Postn", "Kitl"), knn = 5)
MAGIC_fibros_aml <- magic(t(logcounts(data_fibros_aml)), genes=c("Tnc", "Spp1", "Kitl"), knn = 10, t = 10)

#### plot the results
p_fibros_healthy <- ggplot(MAGIC_fibros_healthy) +
  geom_point(aes(x=Col1a1, y=Postn, color=Kitl)) + scale_color_viridis(option="B") +
  coord_cartesian(xlim = c(min(MAGIC_data_aml$result$Dock1), max(MAGIC_data_healthy$result$Dock1)), 
                  ylim = c(min(MAGIC_data_aml$result$Elmo1), max(MAGIC_data_healthy$result$Elmo1))) +
  ggtitle("Clusters 1 and 5 healthy ECs")


p_fibros_aml <- ggplot(MAGIC_fibros_aml) +
  geom_point(aes(x=Tnc, y=Spp1, color=Kitl)) + scale_color_viridis(option="B") 
# +
#   coord_cartesian(xlim = c(min(MAGIC_data_aml$result$Dock1), max(MAGIC_data_healthy$result$Dock1)), 
#                   ylim = c(min(MAGIC_data_aml$result$Elmo1), max(MAGIC_data_healthy$result$Elmo1))) +
#   ggtitle("Clusters 2:4 D10 and D60")
# 
plot_grid(p_healthy,p_aml) 

# 3D plots not very good
plot3d( 
  x=MAGIC_fibros_aml$result$Tnc, y=MAGIC_fibros_aml$result$Spp1, z=MAGIC_fibros_aml$result$Kitl, 
  col = viridis(option = "B", n=3), 
  type = 's', 
  radius = .1,
  xlab="Sepal Length", ylab="Sepal Width", zlab="Petal Length")
