library(SingleCellExperiment)
library(strex)
library(purrr)
library(visNetwork)
library(dplyr)
library(scater)
library(scran)
library(batchelor)
library(bluster)
library(tidyr)
library(here)
library(dittoSeq)

#### load the data
data_bs <- readRDS(here("Data_AML_mouse/Data_Basel/transfer_184642_files_3dc8eea7/sce_withTSNE_29092023.rds"))
data_bs

#### explore cell annotation
table(data_bs$final_annotation)
table(data_bs$Type)
table(data_bs$Batch)
head(colData(data_bs))
table(data_bs$Day)
table(data_bs$Mouse)

plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "final_annotation")
plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Batch")
plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Day")
plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Mouse")

table(cell_type=data_bs$final_annotation, mouse=data_bs$Mouse)

#### check the MLL-AF9 expression and its target genes
fusion_gene_targets <- c("Meis1", "Hoxa9", "Eya1")

plotReducedDim(data_bs, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "MLL-AF9")

multi_dittoPlot(data_bs, vars = "MLL-AF9", 
                group.by = "Mouse", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 1,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

multi_dittoPlot(data_bs, vars = fusion_gene_targets, 
                group.by = "Day", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Lognormalized counts", ncol = 3,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))
