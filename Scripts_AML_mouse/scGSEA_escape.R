library(escape)
library(dittoSeq)
library(SingleCellExperiment)
library(here)
library(scater)
library(forcats)

#### Load the data
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
# before old object with 4 clusters for ECs were used
data_ecs <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5clusters.rds"))

plotReducedDim(data_ecs, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_ecs, dimred = "UMAP_harmony", 
               colour_by=rownames(data_ecs[grep("Sele$",rownames(data_ecs))]))

###  Get the gene sets from MsigDB
GS_MF <- getGeneSets(library = "C5", species = "Mus musculus", subcategory = "GO:MF")
GS_BP <- getGeneSets(library = "C5", species = "Mus musculus", subcategory = "GO:BP")
GS_TFs <- getGeneSets(library = "C3", species = "Mus musculus", subcategory = "TFT:GTRD")
GS.hallmark <- getGeneSets(library = "H", species = "Mus musculus")
GS_Nfe2l2 <- getGeneSets(library = "C3", species = "Mus musculus", 
                         gene.sets = c("NFE2L1_TARGET_GENES", "NFE2L3_TARGET_GENES",
                                       "NRF2_01", "NRF2_Q4"))


GS_GO <- getGeneSets(library = "C5", species = "Mus musculus", 
                     gene.sets =  c("GOMF_CELL_ADHESION_MOLECULE_BINDING",
                                    "GOMF_ACTIN_FILAMENT_BINDING",
                                    "GOMF_INTEGRIN_BINDING"))


#### Calculate NES per cell
# fileter duplicate cells
data_fibros_filt <- data_fibros[,!duplicated(colnames(data_fibros))]
data_ecs_filt <- data_ecs[,!duplicated(colnames(data_ecs))]
# add metadata about time point 
data_ecs_filt$state <- fct_collapse(data_ecs_filt$label,
                                   ECs_healthy_d10 = c("1", "3"),
                                   ECs_d60 = c("2", "4"))


ES.sce <- enrichIt(obj = data_ecs_filt, 
                   gene.sets = GS_Nfe2l2,
                   cores = 2, 
                   min.size = 5,
                   maxRank = 2000,
                   method = "UCell")

# add NES in the sce object
met.data <- merge(colData(data_ecs_filt), ES.sce, by = "row.names", all=TRUE)
row.names(met.data) <- met.data$Row.names
met.data$Row.names <- NULL
colData(data_ecs_filt) <- met.data
dim(colData(data_ecs_filt))

# visualize the NES with heatmap
dittoHeatmap(data_ecs_filt, genes = NULL, metas = names(ES.sce), 
             annot.by = "label", 
             fontsize = 7,
             cluster_cols=FALSE)



# add label info in the enrich result
ES.sce$label <- data_fibros_filt$label
#ES.sce$state <- data_ecs_filt$state

# run the significance tests on the enrichIt result
# not very successful:(
output <- getSignificance(ES.sce, group = "label", fit = "KW")
head(output)

output_ord <- output[order(output[,6], 
                                 decreasing=FALSE),]
output_ord[1:50,]

# select specific terms for violin plot
multi_dittoPlot(data_ecs_filt, vars = c("NFE2L1_TARGET_GENES", "NFE2L3_TARGET_GENES",
                                        "NRF2_01", "NRF2_Q4"), 
                group.by = "label", plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Normalized Enrichment Scores", ncol = 2,
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))

