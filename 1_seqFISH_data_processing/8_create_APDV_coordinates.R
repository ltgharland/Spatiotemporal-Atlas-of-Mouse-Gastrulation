library(princurve)
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(patchwork)

spe = readRDS("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/analysis_output/full_processed_data/seqFISH_spe.Rds")

# this section is run interactively as I need to use locator() function to draw out the regions
embryoval = "embryo3_2"

plot(reducedDim(spe, "spatial")[spe$embryo == gsub("_.*", "", embryoval),], asp = 1, pch = 16, cex = 0.5)
loc = locator() # the AP axis curve
emb = locator() # polygon of the cells to be considered inside embryo
cent = locator() # single point of the centre i.e. a point thats considered ventral

write.table(cbind(loc$x,loc$y), file = paste0("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/", embryoval, "_APDV.tsv"), sep = "\t", quote = FALSE)
write.table(cbind(emb$x,emb$y), file = paste0("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/", embryoval, "_emb.tsv"), sep = "\t", quote = FALSE)
write.table(cbind(cent$x,cent$y), file = paste0("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/", embryoval, "_cent.tsv"), sep = "\t", quote = FALSE)


####### below is code to create APDV coordinates for embryonic cells


APDV_all = matrix(NA, nrow = nrow(reducedDim(spe)), ncol = 2,
                  dimnames = list(colnames(spe), c("AP", "DV")))

e85s = c("embryo1_1", "embryo1_2", "embryo2", "embryo3_1", "embryo3_2")

for (embryoval in c(e85s, paste0("embryo", 4:7))) {
  
  
  loc = read.delim(paste0("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/", embryoval, "_APDV.tsv"), header = TRUE)
  emb = read.delim(paste0("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/", embryoval, "_emb.tsv"), header = TRUE)
  cent = read.delim(paste0("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/", embryoval, "_cent.tsv"), header = TRUE)
  
  if (embryoval %in% e85s) {
    loc <- loc[rev(seq_len(nrow(loc))),]
  }
  
  library(pracma)
  
  s = as.matrix(loc)
  x = reducedDim(spe, "spatial")[spe$embryo == gsub("_.*", "", embryoval),1:2]
  
  x_in = x[inpolygon(x[,1], x[,2], emb[,1], emb[,2]),]
  
  proj <- project_to_curve(x_in, s)
  
  # plot
  plot(x_in, asp = 1)
  lines(s, lwd = 5, col = "blue")
  # points(x_in[proj$ord < 20,])
  segments(x_in[, 1], x_in[, 2], proj$s[, 1], proj$s[, 2])
  
  # add a TRUE/FALSE if the point is closer to the centre than the fitted point
  addsign = function(proj, x_in, cent) {
    
    s0 = proj$s
    s0[,1] <- s0[,1] - cent[1,1]
    s0[,2] <- s0[,2] - cent[1,2]
    s_dists = sqrt(rowSums(s0^2))
    
    x0 = x_in
    x0[,1] <- x0[,1] - cent[1,1]
    x0[,2] <- x0[,2] - cent[1,2]
    x_dists = sqrt(rowSums(x0^2))
    
    ordered_s = proj$s[proj$ord,]
    proj$rank = setNames(seq_len(nrow(ordered_s)), rownames(ordered_s))[rownames(proj$s)]
    
    proj$nearCentre = s_dists <= x_dists
    proj$signedDist = proj$dist_ind * ifelse(proj$nearCentre, -1, 1)
    
    return(proj)
    
  }
  
  proj <- addsign(proj, x_in, cent)
  
  plot(x_in, asp = 1)
  lines(s, lwd = 5, col = "blue")
  segments(x_in[, 1], x_in[, 2], proj$s[, 1], proj$s[, 2], col = factor(proj$nearCentre))
  
  
  APDV = cbind(proj$rank, proj$signedDist)
  
  spe_sub = spe[, rownames(APDV)]
  
  reducedDim(spe_sub, "APDV") <- APDV
  
  source("/Users/ghazan01/Dropbox/Backup/Collaboration/Harvey_VCCRI_2022/CARTANA_heart/scripts/plotReducedDimGirafe.R")
  
  plotReducedDimGirafe(spe_sub, reducedDimNames = c("spatial", "APDV"))
  
  APDV_all[rownames(APDV), ] <- APDV
  
}

reducedDim(spe, "APDV") <- APDV_all

plotReducedDimGirafe(spe[, spe$embryo == "embryo1"], reducedDimNames = c("spatial", "APDV"))

saveRDS(APDV_all, file = "/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/Data/APDV/APDV.Rds")
