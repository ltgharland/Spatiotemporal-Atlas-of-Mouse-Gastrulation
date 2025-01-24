# merge MoleculeExperiment objects for every embryo and z-slice combination

library(MoleculeExperiment)
library(ggplot2)

# new merging function
mergeME = function(meList, molsName = "detected", bdsName = "cell") {
  
  # meList is a list (not named) of ME objects
  # molsName is the assay name for the molecules slot
  # bdsName is the assay name for the boundaries slot
  
  # currently can only handle merging MEs with just 1 assay for molecules and boundaries respectively
  
  molsList = lapply(meList, MoleculeExperiment::molecules, assayName = molsName, flatten = TRUE)
  
  bdsList = mapply(function(me, i) {
    bds = MoleculeExperiment::boundaries(me, assayName = bdsName, flatten = TRUE)
    bds$segment_id <- paste0(i,".",bds$segment_id)
    return(bds)
  },meList, seq_len(length(meList)),
  SIMPLIFY = FALSE)
  
  mols = do.call(rbind, molsList)
  bds = do.call(rbind, bdsList)
  
  molsList = dataframeToMEList(mols,
                               dfType = "molecules",
                               assayName = molsName,
                               factorCol = "feature_id")
  
  bdsList = dataframeToMEList(bds,
                              dfType = "boundaries",
                              assayName = bdsName,
                              factorCol = "segment_id")
  
  me = MoleculeExperiment(molecules = molsList, boundaries = bdsList)
  
  return(me)
}


# analysis

# now perform merging for all the embryo and z-slice combinations
library(dplyr)

dir = "../analysis_output/full_processed_data/me/"
outdir = "../analysis_output/full_processed_data/me_embryo_z/"

allMEnames = list.files(dir, pattern = "_me.Rds")
MEgroups = unlist(lapply(strsplit(allMEnames, "_"), function(x) paste0(x[1], "_", x[3])))

MEnames = split(allMEnames, MEgroups)

for (i in names(MEnames)) {
  
  nms = MEnames[[i]]
  
  MElist = sapply(paste0(dir,nms), readRDS, simplify = FALSE)
  
  me = mergeME(MElist)
  
  g = ggplot_me() + 
    geom_polygon_me(me, colour = "grey", fill = NA) + 
    geom_point_me(me, selectFeatures = c("T", "Eomes"), byColour = "feature_id", size = 0.01) + 
    scale_colour_manual(values = c("T" = "blue", "Eomes"= "red"))
  print(g)
  
  saveRDS(me, file = paste0(outdir, i, ".Rds"))
}



# generate graphs overlaid on the dapi

imgDir = "../stitching/stitched_rotated/"
meDir = "../analysis_output/full_processed_data/me_embryo_z/"

for (i in names(MEnames)) {
  
  me = readRDS(paste0(meDir, i, ".Rds"))
  imgPath = paste0(imgDir, i, ".tif")
  if (!file.exists(imgPath)) next
  img = readImage(imgPath)
  
  # make smaller by 4x factor
  fac = 4
  img_small <- EBImage::resize(img, w = dim(img)[1]/fac, h = dim(img)[2]/fac)
  
  g = ggplot_me() + 
    geom_raster_img(image = img_small,#/quantile(img_small,0.99),
                    pixelSize = 0.1108398*fac) + # manually copy-pasted
    geom_polygon_me(me, colour = "red", fill = NA, size = 0.1)
  ggsave(g, file = paste0("../Figures/early_segmentation_",i,".png"), height = 10, width = 10)
  
  
}

