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
library(ggplot2)
library(here)



files_list <- list.files(paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/Markers_AML_ECs/")), pattern = "_ECs_d60")
files_list


for (i in 1:length(files_list)) {
  # Read files with markers per cluster
  markers <- read.csv(paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/Markers_AML_ECs//"), 
                            sep = "", files_list[i]))
  
  for (j in c("BP", "MF", "CC")) {
    GO <- enrichGO(markers$X, 
                   "org.Mm.eg.db", 
                   ont=j, 
                   keyType = "SYMBOL",
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.01, 
                   qvalueCutoff = 0.05)
    # Write results in csv files and visualize top30 in a barplot
  write.csv(GO, paste(here("Results_AML_mouse/Results_AML_mouse_GO/Results_AML_GO_ECs///"), "GO_", j,  files_list[i], sep=""))
  }
}

genes <- read.csv(paste("./Results/Results_Merge4/Markers_Merge4/", files_list[[1]], sep = ""))
genes_ord <-  genes[order(genes$mean.logFC.cohen, 
                                 decreasing=TRUE),]


# run GSEA analysis not sure if good for scRNAseq
for (i in 1:length(files_list)) {
  # Read files with markers per cluster
  markers <- read.csv  (paste("./Results/Results_", condition, "/Markers_", condition,
                              "/", sep = "", files_list[i]))
  gene_list <- markers$mean.AUC
  
  # edit rownames
  rn_markers <- c((markers[,1]))
  sp <- sapply(strsplit(rn_markers, split='.', fixed=TRUE), function(x) (x[1]))
  names(gene_list) <- sp
  gene_list <- sort(gene_list, decreasing = TRUE)
  
    gse <- gseGO(geneList=gene_list, 
                 ont ="ALL", 
                 keyType = "ENSEMBL",
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = "org.Mm.eg.db", 
                 pAdjustMethod = "none")
    # Write results in csv files and visualize top30 in a barplot
    write.csv(gse, paste("./Results/Results_", 
                        condition, "/", 
                        "Results_GO_", condition, "/", 
                        "gseGO_meanAUC_", "_", files_list[i], sep=""))
    dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) +
    ggtitle(paste("Top30", " gseaGO CAFs", files_list[i], sep = "")) + 
      xlab("gene count")
    
    ggsave(paste("./Plots/Plots_", 
                 condition, "/", 
                 "Plots_GO_", condition, "/", 
                 "gseGO_meanAUC_", files_list[i], device=".pdf", sep=""),
           width = 15, height =10) 
    
  
}



# comparing clusters
onto <- "MF"
clusters <- list()
for (i in c(1:8)) {
  markers <- read.csv(paste("./Results_UnaffectedLung_D23/Markers_Fibros/Unaffected_Markers_Fibro_C", i,  ".csv",
                            sep=""))  
  sp <- sapply(strsplit(markers[,1], split='.', fixed=TRUE), function(x) (x[1]))
  names_l <- paste("Fibros UnaffectedLung_D23", "_", "C", i, sep = "")
  clusters[[names_l]] <- sp  
  
  comp_cl <- compareCluster(geneCluster=clusters, 
                            fun = "enrichGO", 
                            pvalueCutoff = 0.01,  
                            qvalueCutoff = 0.05, 
                            OrgDb='org.Mm.eg.db',
                            ont=onto,
                            keyType = "ENSEMBL")
  
  dotplot(comp_cl, showCategory = 15) + ggtitle(paste("GO ", onto, "Fibros Unaffected Lung D23", sep = ""))
  ggsave(paste("./Plots_UnaffectedLung_Fibros_D23/Plots_GO_Fibros/", "Fibros_UnaffectedLung_D23_", onto, ".pdf", sep=""), height = 30, width = 40) 
  
}

