# extract smFISH quantification according to the 3D segmentation

set.seed(2020)

sce = readRDS("../analysis_output/full_processed_data/seqFISH_spe_all.Rds")

library(SingleCellExperiment)
library(scater)
library(scran)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(limma)
library(patchwork)
library(reshape)
library(gtools)
library(scattermore)
library(batchelor)
library(igraph)
library(cowplot)
library(ggpubr)
library(pracma)
library(princurve)
library(scDblFinder)
library(EBImage)
library(ks)
library(tiff)
library(imager)

source("spatial_functions.R")
source("segmentation_functions.R")
source("celltype_colours.R")
source("neighbourSegments.R")

path_subs = c("embryo5" = "TimEmbryos-021920",
              "embryo4" = "TimEmbryos-030320",
              "embryo6" = "TimEmbryos-020420",
              "embryo7" = "TimEmbryos-020420")
fov_subs = list("embryo5" = 0:6,
                "embryo4" = 0:5,
                "embryo6" = 0:7,
                "embryo7" = 8:15)

smFISH_filename = "../analysis_output/early/early_smFISH_quantification_3D.Rds"

if (!file.exists(smFISH_filename)) {
  ch_tab = matrix(1:24, nrow = 6, ncol = 4)
  
  quant_all = NULL
  
  # hybval = 17
  for (hybval in c(17:28)) {
    
    print(c("hybval", hybval))
    
    for (embryo in rev(names(path_subs))) {
      
      print(c("embryo", embryo))
      
      for (pos in fov_subs[[embryo]]) {
        
        print(c("pos", pos))
        
        seg = readRDS(paste0("../analysis_output/full_processed_data/seg/",embryo, "_Pos",
                             pos, "_seg.Rds"))
        
        sm = readTIFF(paste0("/Users/ghazan01/Nobackup/DataFromNico/",
                             gsub("TimEmbryos-", "", path_subs[embryo]),
                             "/Experiment/HybCycle_",
                             hybval,
                             "/MMStack_Pos",
                             pos,
                             ".ome.tif"),
                      all = TRUE, info = TRUE, indexed = TRUE, as.is = TRUE)
        
        
        for (ch in c(1:3)) {
          
          print(c("ch", ch))
          
          quant_sum_list = list()
          quant_area_list = list()
          quant_name_list = list()
          
          for (z in 2:6) {
            
            print(c("z", z))
            
            # view with auto brightness contrast
            # plot(as.cimg(sm[[2]]))
            
            ints = sm[[ch_tab[z,ch]]]
            ints_long = reshape::melt(ints)
            fit = lm(value ~ X1 + X2, data = ints_long)
            ints_norm <- matrix(fit$residuals + mean(ints_long[,3]), nrow(ints), ncol(ints))
            
            quant = computeFeatures(t(seg[,rev(seq_len(dim(seg)[2])),z]), ints,
                                    methods.ref = c("computeFeatures.basic"),
                                    methods.noref = c("computeFeatures.shape"))
            
            if (is.null(quant)) next
            
            quant_sum <- quant[,"x.Ba.b.q095"]
            # quant_sum <- quant[,"x.Ba.b.mean"]
            quant_name <- rownames(quant)
            
            quant_sum_list[[z]] <- quant_sum
            quant_name_list[[z]] <- quant_name
            
          }
          quant_sum_all = unlist(quant_sum_list)
          quant_name_all = unlist(quant_name_list)
          
          quant_sum = tapply(quant_sum_all, quant_name_all, FUN = max)
          
          quant_mean = quant_sum 
          
          quant = data.frame(b.mean = quant_mean)
          
          rownames(quant) <- paste0(embryo,"_","Pos",pos,"_cell",names(quant_mean))
          
          quant$hyb = hybval
          quant$embryo = embryo
          quant$pos = pos
          quant$channel = ch
          quant$uniqueID = rownames(quant)
          
          quant$x <- colData(sce)[rownames(quant), "dim1"]
          quant$y <- colData(sce)[rownames(quant), "dim2"]
          
          quant_all <- rbind(quant_all, quant)
        }
      }
      
    }
  }
  
  barcoding_info = read.csv("../github/SpatialMouseAtlas2020/barcoding_info/Sequential_readout-ID_v3_TL_ID_SG.csv",
                            header = TRUE)

  sequential_names = setNames(barcoding_info$X.NAME.,
                              paste0(as.character(barcoding_info$HybCycle), ".", as.character(barcoding_info$Channel)))
  
  quant_all$geneID = sequential_names[paste0(
    as.character(quant_all$hyb), ".", as.character(quant_all$channel)
  )]
  
  smFISH_feats = tapply(quant_all$b.mean, 
                        list(quant_all$geneID, quant_all$uniqueID),
                        mean)
  
  saveRDS(quant_all, file = "../analysis_output/early/early_smFISH_quantification_3D.Rds")
  saveRDS(smFISH_feats, file = "../analysis_output/early/early_smFISH_matrix_3D.Rds")
}


# some further visualisation (just for exploration)
exprs = smFISH_feats[,intersect(colnames(smFISH_feats), colnames(sce)), drop = FALSE]
cData = colData(sce[,intersect(colnames(smFISH_feats), colnames(sce))])

sme = SingleCellExperiment(assays = list(logcounts = (exprs)),
                           colData = cData,
                           reducedDims = list(spatial = cData[,c("dim1","dim2")]))
rownames(sme)

sme$total = colSums(logcounts(sme))

gene = "Nkx2-5"

sme <- sme[,order(assay(sme, "logcounts")[gene,])]

plotReducedDim(sme, "spatial", color_by = gene, other_fields = c("embryo", "pos"),
               point_size = 1) + 
  facet_wrap(~embryo, nrow = 1) +
  # facet_wrap(~embryo + pos, nrow = 1) +
  # scale_colour_gradient2(low = "grey", mid = "cornflowerblue", high = "black") +
  coord_fixed()

