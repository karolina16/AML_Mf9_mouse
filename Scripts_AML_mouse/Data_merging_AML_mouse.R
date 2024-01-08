library(BiocSingular)
library(scDblFinder)
library(DropletUtils)
library(readxl)
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
library(gridExtra)

set.seed(1000)

# This code was used to merge the 2 AML experiments
# Rough annotation was also performed and ECs and fibros were saved

#### Data loading
files <- list.files(here("Data_AML_mouse"), pattern = "0_exp2|0_exp1|MF9_Ctrl")
files

### sample list
dir_meta <- here("Data_AML_mouse/Metada_AML/Metadata_AML_mouse.xlsx")
exp_nb <- read_excel(dir_meta, range = cell_cols("B:B")) %>% .$Exp_nb
time_point <- read_excel(dir_meta, range = cell_cols("C:C")) %>% .$Time_Point
cond <- read_excel(dir_meta, range = cell_cols("D:D")) %>% .$Condition
exp_timePoint <- read_excel(dir_meta, range = cell_cols("E:E")) %>% .$Exp_TimePoint
#project_name <- read_excel(dir_meta, range = cell_cols("F:F")) %>% .$Project_name

### import cellranger files from different data sets and create SCEs
for (i in seq_along(files)){
  assign(paste0("data_", files[i]), read10xCounts(paste0(here("Data_AML_mouse", files[i], "/", "filtered_feature_bc_matrix"))))
}

## create SCE list, add meta data and change rownames to Symbol
sce_list <- list()

for (i in seq_along(files)) {
  data_new <- paste("data_", files[i], sep = "")
  data_new <- get(data_new)
  data_new$exp_nb <- exp_nb[i]
  data_new$time_point <- time_point[i]
  data_new$cond <- cond[i]
  data_new$exp_time_point <- exp_timePoint[i]
  #data_new$project_name <- project_name[i]
  rownames(data_new) <- rowData(data_new)[,2]
  sce_list[[length(sce_list) + 1]] <- data_new
  
} 

sce_list <- sce_list %>% set_names(exp_timePoint)

#### run QC for exp1 and exp2 ####
sce_list <- lapply(sce_list, function(sce){
  is.mito <- grepl("^MT-", rownames(sce), ignore.case = TRUE)
  sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))
})

sce_list

#### plot the QC metrics distributions before filtering
for (i in 1:length(sce_list)) {
  qc_plot <- gridExtra::grid.arrange(
    plotColData(sce_list[[i]], x="exp_time_point", y="sum") +  
      scale_y_log10() + ggtitle("Total count") +
      geom_hline(yintercept=800, color="red"),
    plotColData(sce_list[[i]], x="exp_time_point", y="detected") +  
      scale_y_log10() + ggtitle("Detected features") +
      geom_hline(yintercept=300, color="red"),
    plotColData(sce_list[[i]], x="exp_time_point", y="subsets_Mito_percent") + 
      ggtitle("Mito percent") +
      geom_hline(yintercept=10, color="red"),
    ncol=1
  )
  ggsave(paste(here("Plots_AML_mouse/Plots_Exp1_Exp2_Con_EDA/Plots_Exp1_Exp2_QC/"), 
               names(sce_list)[[i]], "_QC.pdf", sep = ""), qc_plot)
}


# use this to assess thresholds
summary(sce_list[[1]]$detected)

# filter according to the fixed threholds
sce_list <- lapply(sce_list, function(sce){
  sce <- sce[,sce$detected > 300 & sce$detected < 10000]
  sce <- sce[,sce$sum > 800 & sce$sum < 100000]
  sce <- sce[,sce$subsets_Mito_percent < 10]
})

#### run mean-variance modeling ####
# run lognormalization
sce_list <- lapply(sce_list, function(sce){
  sce <- logNormCounts(sce)
})


sce_list <- sce_list %>% set_names(exp_timePoint)

#determine variance components for each sce separately
dec_list <- lapply(sce_list, function(sce){
  dec <- modelGeneVar(sce)
})

#common genes
universe <- sce_list %>% map(.,rownames) %>% Reduce(intersect, .)
length(universe)

# Subsetting the SingleCellExperiment objects
sce_list <- lapply(sce_list, function(sce){
  sce <- sce[universe,]
})
dec_list <- lapply(dec_list, function(dec){
  dec <- dec[universe,]
})

combined.dec <- do.call(combineVar, dec_list)
chosen.hvgs <- combined.dec$bio > 0
sum(chosen.hvgs)

hvgs <- getTopHVGs(combined.dec, prop=0.2)

sce_list <- lapply(names(sce_list), function(name){
  sce <- sce_list[[name]]
  mcols(sce) <- mcols(sce_list[[1]])
  sce
})



#### merge the objects
sce_com <- do.call(SingleCellExperiment::cbind, sce_list)
dim(sce_com)

# add colnames
col_names <- paste(sce_com$batch, sce_com$Barcode, sep="_")
colnames(sce_com) <- col_names

# compute normalized logcounts
sce_com <- multiBatchNorm(sce_com, batch = as.factor(sce_com$exp_nb))


#run dimension reduction
sce_com <- denoisePCA(sce_com, technical=combined.dec, subset.row=hvgs)
sce_com <- runUMAP(sce_com, dimred = 'PCA', external_neighbors=TRUE)

plotUMAP(sce_com, colour_by="cond")
plotUMAP(sce_com, colour_by="time_point")

# run doublet detection
dbl.dens <- computeDoubletDensity(sce_com, subset.row=hvgs, 
                                  d=ncol(reducedDim(sce_com)))
summary(dbl.dens)

sce_com$DoubletScore <- dbl.dens
plotUMAP(sce_com, colour_by="DoubletScore")

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call")
summary(dbl.calls)

# assign doubled classification in the meta data
sce_com$Dbl_classification <- dbl.calls

# remove doublets
sce_com <- sce_com[,!sce_com$Dbl_classification == "doublet"]
dim(sce_com)

plotUMAP(sce_com, colour_by="cond")

# run  clustering
kgraph.clusters <- clusterRows(reducedDim(sce_com, "PCA"),
                               TwoStepParam(
                                 first=KmeansParam(centers=500),
                                 second=NNGraphParam(k=10)
                               )
)
table(kgraph.clusters)
colLabels(sce_com) <- factor(I(kgraph.clusters))


# inspect the clusters in more detail

tab <- table(cluster=sce_com$label, batch=sce_com$exp_nb)
tab
write.csv(tab, here("Results_AML_mouse/Clusters_Exp_nb_AML_Exp1_Exp2.csv"))

tab2 <- table(cluster=sce_com$label, batch=sce_com$time_point)
tab2
write.csv(tab2, here("Results_AML_mouse/Clusters_TimePoint_AML_Exp1_Exp2.csv"))

tab3 <- table(cluster=sce_com$label, batch=sce_com$exp_time_point)
tab3
write.csv(tab3, here("Results_AML_mouse/Clusters_Exp_TimePoint_AML_Exp1_Exp2.csv"))

plotUMAP(sce_com, colour_by="label", text_by="label")
plotUMAP(sce_com, colour_by="time_point", text_by="label")

# do a quick check where the stromal cell are with Col1a1, Col1a2, Pdgfra, Pdgfrb
plotUMAP(sce_com, colour_by=rownames(sce_com[grep("Cd79b$",rownames(sce_com))]),
         theme_size=25)



#### run markers to optimize clusters
# Find markers for each cluster
markers <- scoreMarkers(sce_com, lfc=0.25)

# Save markers for all clusters
for (i in names(markers)) {
  markers_ordered <- markers[[i]][order(markers[[i]]$mean.AUC, 
                                        decreasing=TRUE),]
  write.csv(markers_ordered[markers_ordered$mean.AUC>0.5,],
            paste(here("Results_AML_mouse/Results_AML_mouse_markers/Results_AML_Exp1_2/"), 'Markers_C', i, '.csv',
                  sep = ""))
  
}


# save the merged exps
saveRDS(sce_com, here("R_objects_AML_mouse/", "AML_Exp1_Exp2_clusters.rds"))
#### subset the ECs and Fibros
ECs <- sce_com[,sce_com$label %in% c(1,6:8)]
ECs$label <- droplevels(ECs$label)
saveRDS(ECs, here("R_objects_AML_mouse/", "AML_ECs_Exp1_Exp2.rds"))

fibros <- sce_com[,sce_com$label %in% c(2,4,5)]
fibros$label <- droplevels(fibros$label)
saveRDS(fibros, here("R_objects_AML_mouse/", "AML_fibros_Exp1_Exp2.rds"))

