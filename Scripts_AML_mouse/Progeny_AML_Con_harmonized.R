library(OmnipathR)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(progeny)
library(dplyr)
library(decoupleR)
library(dittoSeq)
library(here)

# read in data
data <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_4clusters.rds"))

data <- data[,!duplicated(colnames(data))]

data <- logNormCounts(data)

# analyse pathways with decouplR
net <- get_progeny(organism = 'mouse', top = 100)
net


# infer pathway activities with wmean
# Extract the normalized log-transformed counts
mat <- as.matrix(assay(data, "logcounts"))

# Run wmean
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', times = 100, minsize = 5)
acts

activities <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source')


# add activities to metadata 
activities_t <- t(activities)

new_meta <- cbind(colData(data), activities_t)
colData(data) <- new_meta

plotReducedDim(data, dimred = "UMAP_harmony", colour_by="Hypoxia", text_by="label")


# create new sce with pathways for visualisation
data_progeny <- SingleCellExperiment(assays = list(counts=activities))
colData(data_progeny) <- colData(data)

dittoHeatmap(data_progeny, annot.by = c("label"))
