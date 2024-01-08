library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(SingleCellExperiment)
library(scater)
library(decoupleR)
library(biomaRt)
library(OmnipathR)
library(here)
library(dittoSeq)

### load the data and filter duplicate cells
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
data_filt <- data_fibros[,!duplicated(colnames(data_fibros))]
data_ecs <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5clusters.rds"))
data_filt <- data_ecs[,!duplicated(colnames(data_ecs))]


#### read in the TF regulatory interactions
## We read CollectTRI Regulons for mouse
## new collection is used instead of Dorothea CollectTRI
net <- get_collectri(organism='mouse', split_complexes=FALSE)
net


#### infer TF activities
# Extract the normalized log-transformed counts
mat <- as.matrix(assay(data_filt, "logcounts"))

# Run wmean
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
acts

#prepare df for visualisation
activities <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source')

dim(activities)
# activities are somehow smaller than initial matrix
# subset the intial object accordingly to use cluster assignment
cells <- colnames(activities)
length(cells)
data_sub <- data_filt[,cells]
data_sub

n_tfs <- 50
# Extract activities from object as a long dataframe
df <- t(as.matrix(activities)) %>%
  as.data.frame() %>%
  mutate(cluster = data_sub$label) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks) 


# add activities to metadata 
activities_t <- t(activities)

new_meta <- cbind(colData(data_sub), activities_t)
colData(data_sub) <- new_meta

plotReducedDim(data_sub, dimred = "UMAP_harmony", colour_by="Stat3", text_by="label")
saveRDS(data_sub, here("R_objects_AML_mouse/AML_Exp1_Exp2/CollecTRI_Mes_AML.rds"))

# save regulons for top50 TFs
regs <- net[net$source %in% colnames(top_acts_mat),]  %>% arrange(source)
head(regs)
dim(regs)
regs2 <- regs %>% arrange(source)
regs2

write.csv(regs, here("Results_AML_mouse/Results_AML_CollecTRI/Regs_Top50TFs_ECs_AML.csv"))

regs2 <- net[net$target == "Rgs5",]
regs2
