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


### load the data and filter duplicate cells
data <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
data_filt <- data[,!duplicated(colnames(data))]
data_ecs <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_4clusters.rds"))

#### read in the TF regulatory interactions
## We read CollectTRI Regulons for mouse
## new collection is used instead of Dorothea CollectTRI
net <- get_collectri(organism='mouse', split_complexes=FALSE)
net

### convert targets and symbols to mouse ID
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

humGenes <- c("MYC", "SMAD3", "STAT5A")
conv <- convertHumanGeneList(humGenes)
conv
regs_mice <- get(data("dorothea_mm", package = "dorothea"))

regulons <- regs_mice %>%
  dplyr::filter(confidence %in% c("A","B","C"))

data_test <- run_viper(mat, regulons,.source='tf', .target='target',
                      options = list(method = "scale", minsize = 4, 
                                     eset.filter = FALSE, cores = 1, 
                                     verbose = FALSE))

#### infer TF activities
# Extract the normalized log-transformed counts
mat <- as.matrix(assay(data_filt, "logcounts"))

# Run wmean
acts <- run_wmean(mat=mat, net=regulons, .source='tf', .target='target',
                  .mor='mor', minsize = 5)

acts

# filter results and add as alternative assay
activities <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source')

dim(activities)


# extract alternative experiment for further analysis
data_doro <- swapAltExp(data_sub, "dorothea")


# Rename the label column before reclustering
data_doro$label_rna <- data_doro$label


# cluster cells based on TF activity
# did not yet work for sce, seurat used instead
# Normalisation again probably not necessary
data_doro <- logNormCounts(data_doro, assay.type="tf_activities")
# Feature selection.
dec <- modelGeneVar(data_sub, density.weights=FALSE)

# Visualizing the fit:
fit.data_sce_tp <- metadata(dec)
plot(fit.data_sce_tp$mean, fit.data_sce_tp$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.data_sce_tp$trend(x), col="dodgerblue", add=TRUE, lwd=2)


# Ordering by most interesting genes for inspection.
dec[order(dec$bio, decreasing=TRUE),] 

# Select HVGs
hvg <- getTopHVGs(dec, prop=0.15)

# Dimensionality reduction
set.seed(1234)
data_sub <- runPCA(data_doro, ncomponents=25)
data_sub <- runUMAP(data_sub, dimred = 'PCA', external_neighbors=TRUE)
#data_sub <- runMDS(data_sub)

plotPCA(data_sub, colour_by="label", text_by="label")

plotUMAP(data_sub, colour_by="label", text_by="label") 

# Reclustering, done only if necessary
set.seed(0101010)
kgraph.clusters <- clusterRows(reducedDim(data_sub, "PCA"),
                               TwoStepParam(
                                 first=KmeansParam(centers=1200),
                                 second=NNGraphParam(k=30)
                               )
)
table(kgraph.clusters)
colLabels(data_sub) <- factor(I(kgraph.clusters))


## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(180, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "SMC + CAF5", angle_col = "45",
                       treeheight_col = 0,  border_color = NA,
                       treeheight_row = 0) 

# Save seurat as sce
sce <- as.SingleCellExperiment(seurat)
plotUMAP(sce, colour_by="dorothea_snn_res.0.25") +  
  scale_color_manual(values=c("#598cbc", 
                              "#becfe8", "#5fa54e")) 

"#598cbc" 
"#becfe8"
"#f18e42"
"#5fa54e"
"#a6ce95"
plotUMAP(sce, text_by = "label", 
         colour_by=rownames(sce[grep("Stat3",rownames(sce))]))

ccl19 <- rownames(sce
                  [grepl("Myc", rownames(sce)), ])

# swap to plot expression data
sce_swap <- swapAltExp(sce, "RNA", saved = "dorothea", withColData = TRUE)

# save specific regulons
write.csv(regulons[regulons$tf=="Trp53",], 
          "./Results/Results_Merge4/Results_Dorothea_Merge4/Trp53_Regulons.csv")

saveRDS(sce, "./Data_Merges/Merge4_Dorothea_sce.rds")


# Plot regulons

# C2 #f18e42
# C5 #a6ce95
# C0 #598cbc
# C1 #becfe8
# C4 #5fa54e

tfea <- read.csv("./Results/Results_Merge4/Results_TFEA_Merge4/TFEA_Merge4_Markers_Sub_Cluster_4.csv")
tfea_motif <- strsplit(tfea[4,]$enrichedGenes, ";")

reg <- read.csv("./Results/Results_Merge4/Results_Dorothea_Merge4/Sp1_Regulons.csv")
genes_reg <- intersect(rownames(data_sub), tfea_motif[[1]])

genes_reg2 <- c("Col18a1", "Timp2")

dittoPlotVarsAcrossGroups(data_sub, genes_reg, group.by = "label",
                          main = "Expression of genes with Etv3-binding motif",
                          color.panel = c("#598cbc", "#becfe8", "#f18e42", "#5fa54e", "#a6ce95"),
                          plots = c("vlnplot","boxplot"))

dittoPlot(data_sub, "Sulf1", group.by = "label",plots = c("vlnplot","jitter"))

dittoRidgePlot(data_sub, "Sulf1", group.by = "label")
