# read in the mouse data structure csv
# last edited 19 December 2023

setwd("/Users/ghazan01/Dropbox/Backup/Projects/SpatialEmbryos/scripts/")
set.seed(2023)

####### functions ######

rotate_point <- function(x, y, angle_degrees) {
  # Convert angle from degrees to radians
  angle_radians <- angle_degrees * pi / 180
  
  # Calculate the new coordinates after rotation
  x_prime <- x * cos(angle_radians) - y * sin(angle_radians)
  y_prime <- x * sin(angle_radians) + y * cos(angle_radians)
  
  # Return the new coordinates as a vector
  return(c(x_prime, y_prime))
}

rotate_point_around_center <- function(x, y, cx, cy, theta) {
  # Translate the point and center to the origin
  translated_x <- x - cx
  translated_y <- y - cy
  
  # Rotate the translated point
  rotated_x <- translated_x * cos(theta) - translated_y * sin(theta)
  rotated_y <- translated_x * sin(theta) + translated_y * cos(theta)
  
  # Translate the rotated point back to the original position
  new_x <- rotated_x 
  new_y <- rotated_y 
  
  return(list(new_x, new_y))
}

bestMatch2Dto3D = function(data, jacthresh = 0.5) {
  
  require(igraph)
  
  dfList = list()
  
  for (zval in seq_len(dim(data)[3]-1)) { # because we have 6-1 contiguous pairs
    
    print(zval)
    i_vals = sort(unique(c(data[,,zval])))
    
    nexttab = table(c(data[,,zval+1]))
    
    out = list()
    
    for (i in i_vals) {
      
      if (i %% 20 == 0) print(i)
      
      inds = which(data[,,zval] == i)
      
      n_i = length(inds)
      
      nextvals = data[,,zval+1][inds]
      
      nextmajor = names(sort(table(nextvals), decreasing = TRUE))[1]
      
      n_j = nexttab[nextmajor]
      
      n_ij = sum(nextvals == nextmajor)
      
      jac_ij = n_ij/(n_i + n_j - n_ij)
      
      if (length(jac_ij) == 0) next
      if (jac_ij < jacthresh) next
      
      out[[i]] <- data.frame(z1 = zval, z2 = zval+1, cell_i = i, cell_j = nextmajor, jac = jac_ij)
      
    }
    
    dfList[[zval]] = do.call(rbind, out)
    
  }
  
  df = do.call(rbind, dfList)
  
  if (is.null(df)) {
    node_df = NULL
  } else {
    
    df$comp1 = paste0(df$z1,"_",df$cell_i)
    df$comp2 = paste0(df$z2,"_",df$cell_j)
    
    df_sub <- df
    
    df_sub <- df_sub[order(df_sub$jac, decreasing = TRUE),]
    df_sub <- df_sub[!duplicated(df_sub$comp2),]
    
    
    g = graph.edgelist(as.matrix(df_sub[,c("comp1","comp2")]), directed = FALSE)
    
    gc = components(g)
    
    node_df = data.frame(node = V(g)$name,
                         component = gc$membership)
  }
  
  cellstoadd = setdiff(
    unlist(sapply(seq_len(dim(data)[3]), function(i) paste0(i,"_",unique(data[,,i])), simplify = FALSE)),
    node_df$node)
  
  total_cells = length(cellstoadd) + # cells not merged with anything
    length(unique(node_df$component)) # merged components
  
  # total_cells
  maxval = max(0, max(node_df$component))
  
  full_df = rbind(node_df,
                  data.frame(node = cellstoadd,
                             component = seq_len(length(cellstoadd)) + maxval))
  
  # get the zslices x cells matrix
  mat = matrix(data = NA, ncol = 6, nrow = total_cells)
  colnames(mat) <- seq_len(ncol(mat))
  rownames(mat) <- paste0("cell",seq_len(total_cells))
  
  for (i in 1:nrow(full_df)) {
    
    # get the z-slice and cell name
    cellinfo = as.integer(strsplit(full_df[i,"node"], "_")[[1]])
    newcellname = full_df[i,"component"]
    
    mat[newcellname,cellinfo[1]] <- cellinfo[2]
    
  }
  
  return(mat)
}

rename2Dto3D = function(data, mat) {
  # data is a pixels x pixels x z-slices array with non-matching 2D labels
  # mat is a cells x z-slices matrix containing the new matching
  
  # want to get out a new pixels x pixels z-slices array
  newdata = data
  
  for (zval in seq_len(dim(data)[3])) {
    
    vec = mat[,zval]
    data_2D_z = data[,,zval]
    data_3D_z = match(c(data_2D_z), vec)
    newdata[,,zval] <- data_3D_z
    
  }
  
  return(newdata)
  
}

get3Dcentre = function(data_3D, mat) {
  
  centres = do.call(rbind,sapply(seq_len(nrow(mat)), function(i){
    ind = which(data_3D == i, arr = TRUE)
    colMeans(ind)
  }, simplify = FALSE))
  
  rownames(centres) <- rownames(mat)
  
  return(centres)
}

getJointSPE = function(speList, mat, centres) {
  
  allgenes = sort(unique(unlist(lapply(speList,rownames))))
  allcounts = matrix(0, nrow = length(allgenes), ncol = nrow(mat))
  rownames(allcounts) <- allgenes
  colnames(allcounts) <- rownames(mat)
  
  for (zval in seq_len(length(speList))) {
    
    if (is.null(speList[[zval]])) next
    
    spe = speList[[zval]]
    
    countsvals = as.matrix(counts(spe))
    colnames(countsvals) <- gsub("[0-9]\\.", "cell", colnames(countsvals))
    
    allcounts[rownames(countsvals), colnames(countsvals)] <- allcounts[rownames(countsvals), colnames(countsvals)] + countsvals
    
  }
  
  spe_3D = SpatialExperiment(assays = list(counts = allcounts),
                             colData = centres,
                             reducedDims = list(spatial = centres))
  
  return(spe_3D)
  
}


###### analysis #######


library(EBImage)
library(SpatialExperiment)
library(MoleculeExperiment)
library(abind)
library(ggplot2)
library(rhdf5)

structure = read.csv("../Data/data_structure.csv", row.names = 1)

for (ind in 1:nrow(structure)) {
  
  i = rownames(structure)[ind]
  
  print(i)
  
  if (file.exists(paste0("../analysis_output/full_processed_data/spe/", i, "_spe.Rds"))) next
  
  fovValue = structure[i,"fov"]
  
  fov = read.csv(structure[i,"fov_info_path"])
  fov$pos = fov$fov - 1
  image_width = structure[i,"width_pixels"]
  units_width = c(fov$bound_x_3 - fov$bound_x_1)[1]
  ratio = image_width / units_width
  start_x = min(fov$bound_x_1)*ratio
  start_y = min(fov$bound_y_1)*ratio
  pixel_pos_x = ceiling((fov$bound_x_1*ratio - start_x) + 1)
  pixel_pos_y = ceiling((fov$bound_y_1*ratio - start_y) + 1)
  total_range_x = ceiling(c(min(pixel_pos_x), max(pixel_pos_x) + image_width - 1))
  total_range_y = ceiling(c(min(pixel_pos_y), max(pixel_pos_y) + image_width - 1))
  centre_x = total_range_x[2]/2
  centre_y = total_range_y[2]/2
  #create a blank image
  stitched_image <- Image(0, colormode = "Grayscale", dim = c(total_range_x[2], total_range_y[2]))
  olddim = dim(stitched_image)
  ang = structure[i,"angle"]
  newcorners = rbind(rotate_point(olddim[1]/2, olddim[2]/2, ang),
                     rotate_point(-olddim[1]/2, olddim[2]/2, ang),
                     rotate_point(-olddim[1]/2, -olddim[2]/2, ang),
                     rotate_point(olddim[1]/2, -olddim[2]/2, ang))
  newspan = apply(newcorners,2,range)
  offset_x = min(newspan[,1])
  offset_y = min(newspan[,2])
  fov$pixel_x = pixel_pos_x
  fov$pixel_y = pixel_pos_y
  rownames(fov) <- paste0("fov.",fov$fov)
  
  pixel_to_micron = structure[i,"width_pixels"] / structure[i,"width_um"]
  
  if (structure[i,"segmentation_single"]) {
    
    if (structure[i,"segmentation_filetype"] == "tif") {
      
      img = readImage(structure[i,"segmentation_single_path"], as.is = TRUE)
      
      if (min(img@.Data) == 0) {
        data = img@.Data + 1
      }
      
    }
    
  }
  
  
  if (!structure[i,"segmentation_single"]) {
    
    if (structure[i,"segmentation_filetype"] == "h5") {
      
      # goal to get a 3D data array
      dataList = list()
      for (zval in 1:6) {
        segmentation_filename = structure[i,paste0("segmentation_multiple_path_z",zval)]
        if (is.na(segmentation_filename)) {
          data_z = matrix(NA,
                          nrow = structure[i,"segmentation_width"],
                          ncol = structure[i,"segmentation_width"])
          # next
        } else {
          data_z = h5read(file = segmentation_filename, h5ls(segmentation_filename)$name)
        }
        
        dataList[[zval]] <- data_z
        
      }
      data = abind(dataList, along = 3)
      
    }
    
  }
  
  mat = bestMatch2Dto3D(data, jacthresh = 0.5)
  data_3D = rename2Dto3D(data,mat)
  
  if (dim(data_3D)[1] != structure[i,"width_pixels"]) {
    data_3D <- EBImage::resize(data_3D,
                               w = structure[i,"width_pixels"],
                               h = structure[i,"width_pixels"],
                               filter = "none")
  }
  
  centres_local = as.data.frame(get3Dcentre(data_3D, mat))
  
  centres_local$x_global = centres_local$dim1 - 1 + fov[fovValue,"pixel_x"]
  centres_local$y_global = centres_local$dim2 - 1 + fov[fovValue,"pixel_y"]
  out = rotate_point_around_center(centres_local$x_global,
                                   centres_local$y_global,
                                   centre_x,
                                   centre_y,
                                   (ang/180) * pi)
  centres_local$x_global_affine <- out[[1]] - offset_x
  centres_local$y_global_affine <- out[[2]] - offset_y
  
  centres = centres_local[,c("x_global_affine", "y_global_affine", "dim3")]
  colnames(centres) <- c("dim1", "dim2", "dim3")
  
  if (structure[i,"flip_x"]) {
    centres$dim1 <- ceiling(-offset_x*2) - centres$dim1
  }
  
  if (structure[i,"flip_y"]) {
    centres$dim2 <- ceiling(-offset_y*2) - centres$dim2
  }
  
  centres$dim1 = centres$dim1 / pixel_to_micron
  centres$dim2 = centres$dim2 / pixel_to_micron
  
  # points
  if (structure[i,"molecules_single"]) {
    allpoints = read.csv(structure[i,"molecules_single_path"]) 
    allpointsList = split.data.frame(allpoints, allpoints$z)
  }
  
  speList = list()
  
  for (zval in 1:6) {
    
    print(zval)
    
    if (!structure[i,"molecules_single"]) {
      
      pointsfile = structure[i,paste0("molecules_path_z",zval)]
      if (is.na(pointsfile)) next
      
      pts = read.csv(pointsfile)
    } else {
      if (!as.character(zval) %in% names(allpointsList)) next
      pts = allpointsList[[as.character(zval)]]
    }
    
    # some filtering of the pts according to quality metrics
    # pts <- pts[pts$intensity > (7-pts$seeds)*100,]
    pts <- pts[pts$seeds >= structure[i,"min_seeds"],]
    
    pts$sample = 1
    pts$x_global = pts$x - 1 + fov[fovValue,"pixel_x"]
    pts$y_global = pts$y - 1 + fov[fovValue,"pixel_y"]
    out = rotate_point_around_center(pts$x_global,
                                     pts$y_global,
                                     centre_x,
                                     centre_y,
                                     (ang/180) * pi)
    pts$x_global_affine <- out[[1]] - offset_x
    pts$y_global_affine <- out[[2]] - offset_y
    
    if (structure[i,"flip_x"]) {
      pts$x_global_affine <- ceiling(-offset_x*2) - pts$x_global_affine
    }
    
    if (structure[i,"flip_y"]) {
      pts$y_global_affine <- ceiling(-offset_y*2) - pts$y_global_affine
    }
    
    pts$x_global_affine_um = pts$x_global_affine / pixel_to_micron
    pts$y_global_affine_um = pts$y_global_affine / pixel_to_micron
    
    moleculesMEList <- dataframeToMEList(pts,
                                         dfType = "molecules",
                                         assayName = "detected",
                                         sampleCol = "sample",
                                         factorCol = "geneID",
                                         xCol = "x_global_affine_um",
                                         yCol = "y_global_affine_um")
    
    if (length(unique(c(data_3D[,,zval]))) == 1) next
    
    boundariesMEList = readSegMask(c(0,ncol(data_3D),0,nrow(data_3D)),
                                   image = as.Image(data_3D[,rev(seq_len(dim(data_3D)[2])),zval]),
                                   background_value = 1,
                                   sample_id = 1)
    
    bdsValue = lapply(boundariesMEList[[1]][[1]], function(bds){
      df = bds
      df$x_global = df$x_location -1 + fov[fovValue,"pixel_x"]
      df$y_global = df$y_location - 1 + fov[fovValue,"pixel_y"]
      out = rotate_point_around_center(df$x_global,
                                       df$y_global,
                                       centre_x,
                                       centre_y,
                                       (ang/180) * pi)
      df$x_global_affine <- out[[1]] - offset_x
      df$y_global_affine <- out[[2]] - offset_y
      bds$x_location <- df$x_global_affine
      bds$y_location <- df$y_global_affine
      if (structure[i,"flip_x"]) {
        bds$x_location <- ceiling(-offset_x*2) - bds$x_location
      }
      if (structure[i,"flip_y"]) {
        bds$y_location <- ceiling(-offset_y*2) - bds$y_location
      }
      
      # convert from pixel to micron
      bds$x_location = bds$x_location / pixel_to_micron
      bds$y_location = bds$y_location / pixel_to_micron
      
      return(bds)
    })
    
    boundariesMEList[[1]][[1]] <- bdsValue
    
    me = MoleculeExperiment(molecules = moleculesMEList,
                            boundaries = boundariesMEList)
    
    g = ggplot_me() +
      geom_polygon_me(me, byFill = "segment_id", colour = "black") +
      geom_point_me(me, colour = "black", size = 0.1)
    ggsave(g, file = paste0("../analysis_output/full_processed_data/me/",i,"_z", zval,"_me_plot.png"),
           height = 12, width = 12)
    
    saveRDS(me, file = paste0("../analysis_output/full_processed_data/me/",i,"_z", zval,"_me.Rds"))
    
    spe = countMolecules(me)
    
    speList[[zval]] <- spe
    
  }
  
  spe_3D = getJointSPE(speList, mat, centres)
  
  colnames(spe_3D) <- paste0(i,"_",colnames(spe_3D))
  
  if (structure[i,"flip_x"]) {
    data_3D[,,] <- data_3D[,rev(seq_len(dim(data_3D)[2])),]
  }
  
  if (structure[i,"flip_y"]) {
    data_3D[,,] <- data_3D[rev(seq_len(dim(data_3D)[1])),,]
  }
  
  # calculate segmentation sizes
  tab = table(c(data_3D))
  names(tab) <- paste0(i,"_cell",names(tab))
  
  # add annotations
  spe_3D$embryo <- structure[i,"embryo"]
  spe_3D$pos <- structure[i,"pos"]
  spe_3D$fov <- structure[i,"fov"]
  spe_3D$Experiment <- structure[i,"Experiment"]
  spe_3D$DataStage <- structure[i,"DataStage"]
  spe_3D$Estage <- structure[i,"Estage"]
  spe_3D$pixel_width <- 1/pixel_to_micron
  spe_3D$area_um2 <- as.numeric(unclass(tab[colnames(spe_3D)])/(pixel_to_micron^2))
  
  saveRDS(spe_3D, file = paste0("../analysis_output/full_processed_data/spe/", i, "_spe.Rds"))
  
  saveRDS(data_3D, file = paste0("../analysis_output/full_processed_data/seg/", i, "_seg.Rds"))
  
}
