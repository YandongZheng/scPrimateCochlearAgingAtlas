library(Seurat)
library(hdf5r)
library(DoubletFinder)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)
library(magrittr)
library(harmony)

library(future)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 1000*1024^3)

work_wd <- paste0("/object/08_monkey_cochlea/snRNA_seq/Downstream_analysis/")
save_wd <- paste0("/object/08_monkey_cochlea/snRNA_seq/Downstream_analysis/01_Aging/")

sample <- c("YF1", "YF2", "YF3", "YM1", "YM2", "YM3", 
            "OF1",  "OF2", "OF3", "OM1", "OM2", "OM3", "OM4")

MT_gene <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6","COX1","COX2","COX3","ATP6","ATP8","CYTB")

# ------Quality control -------
scRNAlist <- list()
for (i in seq_len(length(sample))) {
  scRNAlist[[i]] <- readRDS(paste0(work_wd, "Remove_doublecell/", sample[i], 
                                   "/DoublrFider_output/", sample[i], "_remove_doublecell.rds"))
  mt_gene <- rownames(scRNAlist[[i]])[which(rownames(scRNAlist[[i]]) %in% MT_gene)]
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], features = mt_gene)
  DefaultAssay(scRNAlist[[i]]) <- "SCT"
  row.num <- data.frame(Row.num = ncol(scRNAlist[[i]]), Sample = sample[i])
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 500 & 
                             nFeature_RNA < 4000 & percent.mt < 5 & 
                             Double %in% "Singlet")
  row.num$Filter.num <- ncol(scRNAlist[[i]])
  if (i == 1) {
    row.num.all <- row.num
  }else{
    row.num.all <- rbind(row.num.all, row.num)
  }
}

row.num.all$Delete.num <- row.num.all$Row.num - row.num.all$Filter.num
row.num.all %<>% .[, c("Sample", "Row.num", "Filter.num", "Delete.num")]
write.csv(row.num.all, paste0(save_wd, "Monkey_cochlea_filter.num.csv"), row.names = F)


library(qs)
qsave(scRNAlist, paste0(save_wd, "Monkey_cochlea_Aging_noMet_scRNAlist.qs"))



#-----harmony integration-----
library(Seurat)
library(dplyr)
library(magrittr)
library(harmony)
library(qs)
library(ggplot2)


scRNA_harmony <- merge(x = scRNAlist[[1]], y = scRNAlist[2:length(scRNAlist)], merge.data = TRUE)
qsave(scRNA_harmony, file = paste0(save_wd, "Monkey_cochlea_Aging_scRNAlist_harmony_merge.qs"))
gc()

var.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
scRNA_harmony <- RunPCA(object = scRNA_harmony, assay = "SCT", features = var.features,
                        npcs = 100, verbose = T)

system.time({
  scRNA_harmony <- RunHarmony(scRNA_harmony, assay.use = "SCT", group.by.vars = "orig.ident", 
                              plot_convergence = TRUE, dims.use = 1:50)
})

qsave(scRNA_harmony, file = paste0(save_wd, "Monkey_cochlea_Aging_RunHarmony.qs"))