library(CellChat)
library(svglite)
library(extrafont)
library(NMF)
library(ggalluvial)
library(here)
library(scran)
library(scater)
library(forcats)

#### read in the data 
data_mpp <- readRDS(here("Data_AML_mouse/Data_Basel/MPPStemcells_Basel.rds"))
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))
data_ecs<- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/ECs_AML_mouse_Con_Adult1_2_harmonized_5clusters.rds")) # old initial object with 4, 5 clusters were used later


data_fibros <- data_fibros[,!duplicated(colnames(data_fibros))]
data_ecs <- data_ecs[,!duplicated(colnames(data_ecs))]


plotReducedDim(data_ecs, dimred = "UMAP_harmony",colour_by="label", 
               text_by="label")

plotReducedDim(data_ecs, dimred = "UMAP_harmony",colour_by="Plvap")

#### explore the data
data_mpp
head(colData(data_mpp))
plotReducedDim(data_mpp, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "final_annotation")
plotReducedDim(data_mpp, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Day")
plotReducedDim(data_mpp, dimred = "UMAP_MNN_Sample_k50_cellcycle", colour_by = "Thy1")



#### subset the data and add additional annotation
head(colData(data_fibros))
# annotate MSCs
data_fibros$cluster_annotation <- fct_collapse(data_fibros$label,
                                   Mes_AML_1 = c("1"),
                                   Mes_AML_2 = c("2"),
                                   Mes_AML_3 = c("3"),
                                   Mes_healthy_4 = c("4"),
                                   Mes_AML_5 = c("5"),
                                   Mes_AML_6 = c("6"),
                                   Mes_AML_7 = c("7"))


data_fibros$cluster_annotation_state <- fct_collapse(data_fibros$label,
                                                Mes_AML = c("1", "2", "3", "5", "6", "7"),
                                                Mes_healthy = c("4"))

data_fibros$state <- fct_collapse(data_fibros$label,
                                                     Mes_AML = c("1", "2", "3", "5", "6", "7"),
                                                     Mes_healthy = c("4"))

# annotate ECs
data_ecs$cluster_annotation <- fct_collapse(data_ecs$label,
                                               ECs_Cav1_healthy = c("1"),
                                               ECs_Cav1_AML = c("2", "4"),
                                               ECs_Plvap_AML = c("3"),
                                               ECs_Plvap_healthy = c("5"))

data_ecs$cluster_annotation_state <- fct_collapse(data_ecs$label,
                                                     ECs_AML = c("2", "3", "4"),
                                                     ECs_healthy = c("1", "5"))


#### merge ECs and MSCs
data_ecs_msc <- cbind(data_ecs, data_fibros)

col_names <- paste(data_ecs_msc$exp_nb, data_ecs_msc$Barcode, sep="_")
colnames(data_ecs_msc) <- col_names
data_ecs_msc <- data_ecs_msc[,!duplicated(colnames(data_ecs_msc))]

# remove clusters 6 and 7
data_ecs_msc <- data_ecs_msc[,!data_ecs_msc$cluster_annotation %in% 
                               c("Mes_AML_6", "Mes_AML_7")]

data_ecs_msc$cluster_annotation <- droplevels(data_ecs_msc$cluster_annotation)

table(data_ecs_msc$cluster_annotation)

# Create cellChat object
cellchat <- createCellChat(data_ecs_msc, group.by = "cluster_annotation")

# Set the cellChat DB
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# set the used database in the object
cellchat@DB <- CellChatDB

# Show the structure of the database
dplyr::glimpse(CellChatDB$complex)

# Subset the data to save computational cost and preprocess the data 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat) 
df.netP <- subsetCommunication(cellchat, slot.name = "netP") 
df.netP_ord <- df.netP[order(df.netP$prob, decreasing=TRUE),]
write.csv(df.netP_ord, here("Results_AML_mouse/Results_AML_mouse_CellChat/netP_ord_wo6_7.csv"))
#write.csv(cellchat@netP$pathways, here("Results_AML_mouse/Results_AML_mouse_CellChat/Pathways_AML_mouse_Mes_ECs.csv"))

# Check the average expression of ligand-receptor
#computeAveExpr(cellchat, features = c("Vegfa","Vegfr1"))

# Compute the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

# group cell by cell type 

# Visualise the aggregated cell-cell network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Visualise signaling from each cluster
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot 
#shows signaling to fibroblast and the right portion shows signaling to immune cells 
pathways.show <- c("CDH", "CDH5", "VEGF", "FN1", "JAM", "ESAM")
#vertex.receiver = c(1,3) # a numeric vector. 

# Circle plot
# to save as svglite
#svglite(file = "pathway2.svg", width=10, height=10)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show[6], layout = "circle", small.gap = .1)

# Chord plot with specific groups
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = "CDH5", vertex.receiver = vertex.receiver, layout = "circle")


pathways.show.all <- cellchat@netP$pathways
pathways.show.all

for (i in 1:length(pathways.show.all)) {
  gg1 <- netAnalysis_contribution(cellchat, signaling =  pathways.show.all[i])
  ggsave(filename=paste0(here("Plots_AML_mouse/Plots_AML_mouse_CellChat/CellChat_Plots_AML_wo_6_7/"), pathways.show.all[i], "_L-R_contribution_AML_Mouse_Mes_wo_67_EC.pdf"), 
         plot=gg1, height = 20, width = 20)
  gg3 <- plotGeneExpression(cellchat, signaling = pathways.show.all[i]) + geom_boxplot()
  ggsave(filename=paste0(here("Plots_AML_mouse/Plots_AML_mouse_CellChat/CellChat_Plots_AML_wo_6_7/"), pathways.show.all[i], "_GeneExpression_AML_Mouse_Mes_wo_67_EC.pdf"), 
         plot=gg3, width = 20, height = 20)
}

# to plot expression
# plotExpression(data_ecs_msc, c("Lgals9", "Cd44"), 
#                x="cell_type", ncol = 1) + geom_boxplot() + theme_classic(base_size = 25) +
#   scale_x_discrete(guide=guide_axis(n.dodge=2))

# to make chord plots for all pathways
# for (i in 1:length(pathways.show.all)) {
#   pdf(file = paste0("./Plots_IFNAR/Plots_IFNAR_CellChat/", pathways.show.all[i], "_Pathway_Network.pdf"), width=10, height=10)
#   par(mfrow=c(1,1))
#   netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "chord")
#   dev.off()
# }

#netVisual(cellchat, signaling = pathways.show.all[1])

# Visualise specific pathways with a heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show[2], color.heatmap = "Reds")
dev.off()
# visualisation with bubble plot

pdf(file = paste0("./Plots_WT_Infected/Plots_WT_Infected_CellChat/", pathways.show, "_Bubble.pdf"), width=20, height=10)
netVisual_bubble(cellchat, signaling = pathways.show[2], remove.isolate = FALSE)
dev.off()

#### look at specific L-R pairs
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show.all, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[60,]
netVisual_individual(cellchat, signaling = pathways.show.all, pairLR.use = LR.show, layout = "circle")
# save all L-R pairs
write.csv(pairLR.CXCL, "./Results/Results_CellChat/Ligand_Receptor_Pairs_6Clusters.csv")
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL", width = 8, height = 2.5, font.size = 10)

#### aggregated heatmap
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, signaling=pathways.show.all, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, signaling=pathways.show.all, pattern = "incoming", height = 15)
ht1 + ht2

saveRDS(cellchat, "./Data/Fibroblasts_Tcells_6Clusters_CellChat.rds")


#### find patterns
# Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
pat <- c("incoming", "outgoing")
nPatterns <- 4

pdf(file = "./Plots/Plots_CellChat/SelectK_Outgoing_IFNAR.pdf")
selectK(cellchat, pattern = pat[1])
dev.off()

pdf(file = "./Plots_IFNAR//Plots_IFNAR_CellChat/Patterns_Outgoing_Ifnar.pdf", 
    width=30, height=25)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = pat[1], k = nPatterns,
                                          width = 10, height = 20)
dev.off()

netAnalysis_river(cellchat, pattern = pat[1])

# dot plot
netAnalysis_dot(cellchat, pattern = pat[2])
