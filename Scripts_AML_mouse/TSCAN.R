library(scater)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(bluster)
library(pheatmap)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pathview)
library(TSCAN)
library(slingshot)

# This analysis was done only for LeprPos clusters

#### load data
data_fibros <- readRDS(here("R_objects_AML_mouse/AML_Exp1_Exp2/AML_Exp1_Exp2_Harmonized/Mes_AML_Exp1_2_harmonized_7clusters.rds"))

#### quick visual check
plotReducedDim(data_fibros, dimred = "UMAP_harmony", colour_by="label", text_by="label")
plotReducedDim(data_fibros, dimred = "UMAP_harmony", 
               colour_by=rownames(data_fibros[grep("Dlx5$",rownames(data_fibros))]))


# Subset the data to use only Lepr+ clusters
data_sub <- data_fibros[, data_fibros$label %in% c(1:5)]
data_sub$label <- droplevels(as.factor(data_sub$label))

#### run TSCAN
# Obtain centroids
by.cluster <- aggregateAcrossCells(data_sub, ids=colLabels(data_sub))
centroids <- reducedDim(by.cluster, "harmony_pca")

# Set clusters=NULL as we have already aggregated above.
mst <- TSCAN::createClusterMST(centroids, clusters=NULL)
mst

plot(mst)

# plot trajectory
line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="UMAP_harmony")

plotReducedDim(data_sub, dimred = "UMAP_harmony", colour_by="label", text_by="label") + 
  geom_line(data=line.data, mapping=aes(x=UMAP1, y=UMAP2, group=edge),
            size=1.2, linetype="solid") 
  
# Obtain Pseudotime ordering
map.tscan <- mapCellsToEdges(data_sub, mst=mst, use.dimred="harmony_pca")
tscan.pseudo <- orderCells(map.tscan, mst, start = c(4))
head(tscan.pseudo)

# Assign pseudotime values in the meta data
data_sub$TSCAN.first <- pathStat(tscan.pseudo)[,1]

# Plot pseudotime
common.pseudo <- averagePseudotime(tscan.pseudo) 
plotReducedDim(data_sub, dimred = "UMAP_harmony", colour_by=I(common.pseudo), 
  text_by="label", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=UMAP1, y=UMAP2, group=edge),
            size=1.2, linetype="solid")

#### Run DE on the trajectory
pseudo <- TSCAN::testPseudotime(data_sub, pseudotime=tscan.pseudo[,1])[[1]]
pseudo$SYMBOL <- rownames(data_sub)
sorted <- pseudo[order(pseudo$logFC),]

up.left <- sorted[sorted$logFC < 0,]
up.right <- sorted[sorted$logFC > 0,]

write.csv(up.right, here("Results_AML_mouse/Results_TSCAN/LeprPos_DE_TSCAN_Up.csv"))
write.csv(up.left, here("Results_AML_mouse/Results_TSCAN/LeprPos_DE_TSCAN_Down.csv"))

# select the best DEGs to plot
best_up <- tail(up.right$SYMBOL, 30)
best_up
best_down <- head(up.left$SYMBOL, 30)
best_down

# plot epression changes along trajectory
# upregulated genes
plotExpression(data_sub, features=c("Spp1", "Col1a1", "Bglap", "Fat3", "Alpl", 
                                    "Bglap2", "Mef2c", "S100a6", "Lrp4", "Postn",
                                    "Cadm1", "Wif1", "Tnc", "Ncam1", "Kcnk2", 
                                    "Cfh", "Fap", "Lum"),
               x="TSCAN.first", colour_by="label", ncol = 3, show_smooth=TRUE)
ggsave(here("Plots_AML_mouse/Plots_AML_MES_TSCAN/Top15_Up_TSCAN_LeprPos_AML.pdf"),
       width = 8, height = 12)

# downregulated genes
plotExpression(data_sub, features=c("Cxcl12", "Pdzrn4", "Igfbp5", "Ebf3", "Lpl", 
                                    "Kitl", "Adcy2", "Lepr", "Bmp6", "Negr1",
                                    "Kalrn", "Gpc6", "Adipoq", "Chrdl1", "Thsd4", 
                                    "Fbn1", "Runx1", "Tgfbr3"),
               x="TSCAN.first", colour_by="label", ncol = 3, show_smooth=TRUE)
ggsave(here("Plots_AML_mouse/Plots_AML_MES_TSCAN/Top15_Down_TSCAN_LeprPos_AML.pdf"),
       width = 8, height = 12)

# other niche markers
plotExpression(data_sub, features=c("Il34", "Csf1", "Il7", "Tpo", "Epo"),
               x="TSCAN.first", colour_by="label", ncol = 3, show_smooth=TRUE)
ggsave(here("Plots_AML_mouse/Plots_AML_MES_TSCAN/Top15_NicheMarkers_TSCAN_LeprPos_AML.pdf"),
       width = 8, height = 12)

plotExpression(data_sub, features=c("Pparg"),
               x="TSCAN.first", colour_by="label", ncol = 3, show_smooth=TRUE)


# plotExpression(data_sub, features=c("Ccl11", "Tnc"),
#                x="TSCAN.first", colour_by="label", ncol = 1, show_smooth=TRUE) +
#   scale_color_manual(values=c("#598cbc", "#becfe8", "#5fa54e"))
# ggsave("./Plots/Plots_Merge4/Tscan_Selected5_Clusters_0-4_Merge4.pdf",
#        width = 6, height = 8)


