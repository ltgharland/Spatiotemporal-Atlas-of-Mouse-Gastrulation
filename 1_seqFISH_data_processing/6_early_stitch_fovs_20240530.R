# script to generate rotated and stitched images for every embryo and z-slice
# combination

library(EBImage)

## functions

rotate_point <- function(x, y, angle_degrees) {
  # Convert angle from degrees to radians
  angle_radians <- angle_degrees * pi / 180
  
  # Calculate the new coordinates after rotation
  x_prime <- x * cos(angle_radians) - y * sin(angle_radians)
  y_prime <- x * sin(angle_radians) + y * cos(angle_radians)
  
  # Return the new coordinates as a vector
  return(c(x_prime, y_prime))
}

## analysis

data_structure = read.csv("../Data/data_structure.csv", header = TRUE, colClasses = "character")

  for (emb in paste0("embryo",1:7)) {
  
  data_structure_sub = subset(data_structure, embryo %in% emb)
  fovs = as.integer(gsub("fov\\.","",data_structure_sub[,"fov"]))
  
  dir = as.character(data_structure_sub[1,"Experiment"])
  
  dirName = paste0("../Data/TimEmbryos-",dir)
  
  fov = read.csv(paste0(dirName, "/slide_explorer/fovinfo.csv"))
  
  fov <- fov[fovs,]
  
  fov$pos = fov$fov - 1
  
  image_width = 2048
  units_width = c(fov$bound_x_3 - fov$bound_x_1)[1]
  ratio = image_width / units_width
  start_x = min(fov$bound_x_1)*ratio
  start_y = min(fov$bound_y_1)*ratio
  
  pixel_pos_x = ceiling((fov$bound_x_1*ratio - start_x) + 1)
  pixel_pos_y = ceiling((fov$bound_y_1*ratio - start_y) + 1)
  
  total_range_x = ceiling(c(min(pixel_pos_x), max(pixel_pos_x) + image_width - 1))
  total_range_y = ceiling(c(min(pixel_pos_y), max(pixel_pos_y) + image_width - 1))
  
  for (zval in 1:6) {
    
    # create a blank image
    stitched_image <- Image(0, colormode = "Grayscale", dim = c(total_range_x[2], total_range_y[2]))
    
    # read in each of the tifs and populate the stitched image
    
    for (i in seq_len(length(fovs))) {
      
      print(i)
      
      pos = fov[as.character(fovs[i]),"pos"]
      
      tifname = paste0(dirName, "/HybCycle_0_fixed/MMStack_Pos",pos,".ome.tif")
      
      if ((i == 1 & emb %in% paste0("embryo", c(4,5,6))) | 
          emb %in% paste0("embryo", c(1,2,3))) {
        img = readImage(tifname)[,,18+zval]
      } else {
        img = readImage(tifname)[,,4*zval]
      }
      
      # fill the values
      stitched_image[seq(pixel_pos_x[i], by = 1, length.out = image_width),
                     seq(pixel_pos_y[i], by = 1, length.out = image_width)] <- img@.Data
      
    }
    
    writeImage(stitched_image, files = paste0("../stitching/stitched/",emb, "_", "z",zval,".tif"))
    writeImage(stitched_image/quantile(stitched_image, 0.99), files = paste0("../stitching/stitched/",emb, "_", "z",zval,".png"))
    
    olddim = dim(stitched_image)
    
    ang = as.numeric(data_structure_sub[1,"angle"])
    
    newcorners = rbind(rotate_point(olddim[1]/2, olddim[2]/2, ang),
                       rotate_point(-olddim[1]/2, olddim[2]/2, ang),
                       rotate_point(-olddim[1]/2, -olddim[2]/2, ang),
                       rotate_point(olddim[1]/2, -olddim[2]/2, ang))
    
    rotate_point(olddim[1]/2, olddim[2]/2, ang)
    
    newspan = apply(newcorners,2,range)
    
    newdim = ceiling(apply(newspan,2,diff))
    
    # Perform the rotation
    rotated_image <- rotate(stitched_image,
                            angle = ang,
                            output.dim = newdim)
    
    if (data_structure_sub[1,"flip_y"]) {
      rotated_image <- flip(rotated_image)
    }
    
    if (data_structure_sub[1,"flip_x"]) {
      rotated_image <- flop(rotated_image)
    }
    
    writeImage(rotated_image, files = paste0("../stitching/stitched_rotated/",emb, "_", "z",zval,".tif"))
    writeImage(rotated_image/quantile(rotated_image,0.99), files = paste0("../stitching/stitched_rotated/",emb, "_", "z",zval,".png"))
    
  }
  
}
