### dimensionality reduction and cell clustering -----

library(Seurat)
library(dplyr)
library(magrittr)
library(harmony)
library(qs)
library(ggplot2)
library(RColorBrewer)

library(future)
future::plan(strategy = 'multicore', workers = 20)
options(future.globals.maxSize = 1000*1024^3)

save.wd <- paste0("/object/08_monkey_cochlea/snRNA_seq/Downstream_analysis/01_Aging/")

scRNA_harmony <- qread(file = paste0(save.wd, "Monkey_cochlea_Aging_RunHarmony.qs"))

npcs <- 50

plot.wd <- paste0(save.wd, "Plot/")
dir.create(plot.wd)
umap.wd <- paste0(plot.wd, "UMAP.parameter/")
dir.create(umap.wd)


#select_para
select_para <- function(combination){
  combination0 <- combination
  #DefaultAssay(combination0) <- "integrated"
  for (n.neighbors in c(10, 30, 50)) {
    for (min.dist in c(0.1, 0.3, 0.5)) {
      for (spread in c(1, 1.5)) {
        combination0 <- RunUMAP(combination0, reduction = "harmony",
                                dims = 1:npcs,
                                n.neighbors = n.neighbors,
                                min.dist = min.dist,
                                spread = spread)
        p <- DimPlot(combination0, reduction = "umap", label = TRUE, 
                     pt.size = 0.1, raster = F) +
          theme(aspect.ratio = 1) +
          NoLegend()
        ggsave(paste0(umap.wd, "Monkey_cochlea_Aging_umap_", n.neighbors, "_", 
                      min.dist, "_", spread , ".png"), 
               plot = p, width = 8, height = 8)
      }
    }
  }
}


select_para(scRNA_harmony)


# n.neighbors = 10, min.dist = 0.5, spread = 1.5 ------
seurat.object <- scRNA_harmony %>% 
  FindNeighbors(reduction = "harmony", dims = 1:npcs) %>%
  FindClusters(resolution = c(2, 2.5, 3, 3.5, 4, 5))

seurat.object <- RunUMAP(seurat.object, reduction = "harmony", 
                         dims = 1:npcs,
                         n.neighbors = 10,
                         min.dist = 0.5,
                         spread = 1.5) 

qsave(seurat.object, paste0(save.wd, "Monkey_cochlea_Aging_runUMAP.qs"))



for (i in c(2, 2.5, 3, 3.5, 4, 5)) {
  loc.id <- which(colnames(seurat.object@meta.data) == paste0("SCT_snn_res.", i))
  colourCount = length(unique(seurat.object@meta.data[,loc.id]))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  p <- DimPlot(seurat.object, reduction = "umap",label = T,group.by = paste0("SCT_snn_res.", i),
               pt.size = 0.2, label.color = "black", label.size = 5, 
               cols = getPalette(colourCount), raster = FALSE) +  
    scale_size(guide = guide_legend(override.aes = list(size = 5))) +
    theme(aspect.ratio = 1)
  ggsave(paste0(umap.wd, "Monkey_cochlea_Aging_umap_cluster_res.", i, ".png"), 
         plot = p, width = 12, height = 8)
}



# resolution = 2 -----
seurat.object <- scRNA_harmony %>% 
  FindNeighbors(reduction = "harmony", dims = 1:npcs) %>%
  FindClusters(resolution = 2) %>%
  RunUMAP(reduction = "harmony", dims = 1:npcs, n.neighbors = 10, min.dist = 0.5, spread = 1.5) 

qsave(seurat.object, paste0(save.wd, "Monkey_cochlea_Aging_runUMAP.qs"))



# FindAllMarkers -----
cluster.marker <- FindAllMarkers(seurat.object, only.pos = T,min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster.marker, paste0(save.wd,"Macaque_cochlea_NoMet_allcluster_marker_20231102.csv"),row.names = F)

