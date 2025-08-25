# install.packages('pacman')
rm(list = ls()); gc()
library(pacman) 

p_load(psych,curl,progress,
       foreign,gdal,psych,
       future,future.apply,
       geometry,outliers,
       lidR,
       raster,MASS,
       terra,
       devtools,
       tidyverse,
       exactextractr,TreeLS,BiocManager,
       sf,MASS,sp,Rdimtools,data.table,
       parallel, doParallel, foreach,corrplot,rlas,data.table
)
remotes::install_github('tiagodc/TreeLS')

# Set work directory
setwd(paste0(getwd(),'/Forest_Structural_Complexity/'))

################################################################################
#Get x and y coordinates from field plots
shapepath<-"data_Complete.shp"
fplots<-st_read(shapepath)
plot(fplots[1],  col = "blue")
chunk = "/Norm_tile/bc_092h001_2_4_4_xyes_8_utm10_20170713_norm.laz"
las <- readLAS(chunk, select = "xyzirc", filter = "-drop_class 7 -drop_point_redundant -drop_withheld")
las

crs(las)
#"EPSG:6653"="NAD83(CSRS) / UTM zone 10N + CGVD2013 height - CGG2013",which is the LiDAr coordinate system
fplots <- st_transform(fplots, crs = "EPSG:6653")
# get plot radio
fplots$radio<- sqrt(fplots$area*10000/3.14159265359)

hist(fplots$radio,30)


# Remove duplicates
filtplots <- fplots#[!duplicates, ]

# testplot<-fplots[c(1:15),]  #For testing the LiDAR Functions

radius<-15 #11, #25 # Tested a couple of plot radious based on their distribution, 15 m worked the best

# # Apply the function to the catalog
# # Set the number of cores for parallel processing
# Read the LAS catalog with a filter to keep only ground returns (classification == 2)
path<-paste0(getwd(),'/Norm_tile/')
# path<-'G:/Conservation Solution Lab/People/Luizmar/PhD_Luizmar/Norm_tile_23_09_24/'

ctg <- readLAScatalog(path, filter = " -drop_class 7 -drop_z_below 0.5 -drop_withheld -drop_z_above 100 -drop_point_redundant")
#
# # Set various options for the catalog processing
opt_output_files(ctg) <- ""
opt_wall_to_wall(ctg) <- TRUE
opt_progress(ctg) <- FALSE
opt_stop_early(ctg) <- FALSE


# # # Plot the catalog with chunks displayed
# plot(ctg, chunk = TRUE)
# # # Check for empty or corrupted files
# las_check(ctg)
#############################################################################
#set lidar processing resolution
res<-5
# Function to calculate forest complexity
# 1) -- 3D Canopy Entropy (Liu et al. 2022), adapted from De Conto et al. 2024
# https://github.com/LidarSu/Canopy_entropy
avg_point_dist = function(las, h=1){
  if (nrow(las@data) < 2) {
    warning("Not enough points in the LAS data to compute distances.")
    return(data.table(layer = integer(0), dist = numeric(0)))  # Return empty data.table
  }
  las@data$MinDist = nabor::knn(TreeLS:::las2xyz(las), k=2)$nn.dist[,2]
  z_pts = las$Z - min(las$Z)
  las@data$layer = as.integer(z_pts / h)
  avg_dst = las@data[order(layer),.(dist = mean(MinDist)), by=layer]
  return(avg_dst)
}

adaptative_resample = function(las, h=1, alpha = 0.05, step=1){
  avg_dst = avg_point_dist(las, h)

  noise_tol = mean(avg_dst$dist) + 2*sd(avg_dst$dist)
  avg_dst = avg_dst[dist < noise_tol,]
  if (nrow(avg_dst) < 3) return(NULL)
  max_dst = max(avg_dst$dist)
  
  for (ratio in seq(2,5,step)){
    vox = ratio * max_dst
    vlas = tlsSample(las, smp.voxelize(vox))
    avg_dst = avg_point_dist(vlas, h)
    mkt = trend::mk.test(avg_dst$dist)
    if(mkt$p.value > alpha) break
  }

  return(vlas)
}


# Function to calculate canopy entropy in the XY, XZ, and YZ planes
canopy_entropy = function(x,y,z, bw=0.2, grid_size=0.1) {
  las = suppressMessages(LAS(data.table(X=x,Y=y,Z=z), check=F))
  
  if (is.null(las) || nrow(las@data) <= 3) {
    return(list(CExy = NA_real_, CExz = NA_real_, CEyz = NA_real_, CExyz = NA_real_))
  }  
  las = adaptative_resample(las)
  
  if (is.null(las) || nrow(las@data) <= 3) {
    return(list(CExy = NA_real_, CExz = NA_real_, CEyz = NA_real_, CExyz = NA_real_))
  }  
  bounds = apply(las@data, 2, range)

  ce = c()
  
  planes = list(c('X', 'Y'), c('X', 'Z'), c('Y', 'Z'))
  
  for (plane in planes) {
    i = plane[1]
    j = plane[2]
    
    ni = 2 + (diff(bounds[, i]) + 8 * bw) %/% grid_size
    nj = 2 + (diff(bounds[, j]) + 8 * bw) %/% grid_size
    lims = c(bounds[, plane]) + (bw * 4 + (grid_size / 2)) * c(-1, 1, -1, 1)
    
    # KDE calculation
    den = ks::kde(las@data[, ..plane], h=bw, gridsize=c(ni, nj), xmin=lims[c(1, 3)], xmax=lims[c(2, 4)])
    if (is.null(den) || all(den$estimate <= 0)) {
      ce[paste(plane, collapse='')] = NA_real_
      next  # Skip to the next plane
    }
    
    den = den$estimate[den$estimate > 0]
    entropy = -1 * sum(den * log(den) * grid_size^2)
    
    ce[paste(plane, collapse='')] = entropy
  }
  
  ce['XYZ'] = sqrt(sum(ce^2, na.rm = TRUE))  # Use na.rm for safety
  return(list(
    CExy = ce[[1]],
    CExz = ce[[2]],
    CEyz = ce[[3]],
    CExyz = ce[[4]]
  ))
}

# 2) Coefficient of variation
# Custom function to calculate coefficient of variation (CV) and interquartile range (IQR)
metrics_custom <- function(Z) { 
  # Filter points by height
  z <- Z[Z>1]
  if (is.null(z) || length(z) == 0) 
    return(data.frame(CV = NA_real_, IQR = NA_real_))
  list(
    coef_var = sd(z) / mean(z) * 100,  # Coefficient of Variation (%)
    iqr = IQR(z)  # Interquartile Range
  )
}


# 3) Leaf area diversity
# FHD Calculation Function
lai_dvrt <- function(la_results) {
  # Convert la_results (list) to a numeric vector
  lai_values <- unlist(la_results)
  
  # Ensure there are no NA values and the sum of lai_values is not zero
  lai_values <- lai_values[!is.na(lai_values)]
  total_lai <- sum(lai_values)
  
  # Return NA if total LAI is zero to avoid division by zero
  if (total_lai == 0) return(NA_real_)
  
  # Calculate the proportion of returns in each layer (pi)
  pi <- lai_values / total_lai
  
  # Calculate FHD using the Shannon-Wiener index formula
  fhd <- -sum(pi * log(pi), na.rm = TRUE)
  
  return(fhd)
}

# Function to calculate metrics_lad and FHD
metrics_lad <- function(Z, zmin = NA, dz = 1, k = 0.5, z0 = 2) {
  # Filter points by height
  z <- Z
  zmax<-max(z)
  if (is.null(z) || length(z) == 0) 
    return(data.frame(laiherb = NA_real_,
                      laiShrub = NA_real_,
                      laiMid_Ht = NA_real_,
                      laiHigh_H = NA_real_,
                      laiVery_H = NA_real_,
                      fhd = NA_real_))

  if (!is.na(zmin)) z <- z[z > zmin]
  if (length(z) == 0) 
  return(data.frame(laiherb = NA_real_,
                                 laiShrub = NA_real_,
                                 laiMid_Ht = NA_real_,
                                 laiHigh_H = NA_real_,
                                 laiVery_H = NA_real_,
                                 fhd = NA_real_))  # Return NA for all 6 bands if no points
 hght1 = zmax*0.005
 hght2 = zmax*0.05
 hght3 = zmax*0.25
 hght4 = zmax*0.5
 hght5 = zmax*0.75
  # Define height ranges and strata names
  height_ranges <- list(
    Herb = c(hght1, hght2),
    Shrub = c(hght2, hght3),
    Mid_Height = c(hght3, hght4),
    High_Height = c(hght4, hght5),
    Very_High = c(hght5, Inf)
  ) 

  # Calculate LAI for each stratum
  lai_results <- sapply(names(height_ranges), function(stratum) {
    range <- height_ranges[[stratum]]
    filtered_z <- z[z >= range[1] & z < range[2]]
    
    if (length(filtered_z) > 2) {
      ladprofile <- lidR::LAD(filtered_z, dz = dz, k = k, z0 = z0)
      lai <- with(ladprofile, sum(lad, na.rm = TRUE))
    } else {
      lai <- NA_real_  # Return NA_real_ to keep the type consistent
    }
    
    return(lai)
  })
  
  # Calculate FHD using the lai_dvrt function
  fhd <- lai_dvrt(lai_results)
  
  # Return a list of LAI for each stratum and the FHD
  return(list(
    laiherb = lai_results[[1]],
    laiShrub = lai_results[[2]],
    laiMid_Ht = lai_results[[3]],
    laiHigh_H = lai_results[[4]],
    laiVery_H = lai_results[[5]],
    fhd = fhd
  ))
}


# 4) Median absolute deviation of height 
# Function to calculate HMAD
h_mad <- function(Z) {
  # Filter points by height
  z <- Z[Z>1]
  if (is.null(z) || length(z) <= 3) 
    return(data.frame(hmad = NA_real_))
  # Calculate mean height
  z_mean <- mean(z, na.rm = TRUE)
  
  # Calculate the absolute deviations from the mean
  abs_deviation <- abs(z - z_mean)
  
  # Calculate the median of the absolute deviations
  mad <- median(abs_deviation, na.rm = TRUE)
  
  # Calculate HMAD
  hmad <- 1.4826 * mad
  
  return(list(hmad=hmad))
}



# 5 Fractal Dimensions 
# Function to calculate HMAD
fd = function(X, Y, Z) {
  M = cbind(X, Y, Z)
  fdr = est.boxcount(M)$estdim  # Assuming est.boxcount() works correctly
  return(list(fdr=fdr))
}
# 

# 6) Canopy Rugosity = rumple_index(X, Y, Z)

# Addition LiDAR metrics
#7 Canopy cover
##################################################################################
canopy_cover <- function(x,y,z, threshold = 5) {
  las = suppressMessages(LAS(data.table(X=x,Y=y,Z=z), check=F))
  # Filter out points below the canopy height threshold (e.g., 2 meters)
  las_above_threshold <- filter_poi(las, Z > threshold)

  # Calculate the total number of points and the number of points above the threshold
  total_points <- nrow(las@data)
  canopy_points <- nrow(las_above_threshold@data)

  # Check if there are any points to avoid division by zero
  if (total_points == 0) {
    return(NA_real_)  # Return NA if there are no points
  }

  # Compute the canopy cover as the percentage of points above the threshold
  cover <- (canopy_points / total_points) * 100

  return(list(CCover = cover))
}


#Height metrics
HeightMetrics = function(Z){
  # Filter points by height
  z <- Z[Z<100]
  if (is.null(z) || length(z) <= 3) 
    return(data.frame(zmax = NA_real_, 
                      zmean = NA_real_,
                      z_sd = NA_real_,
                      zmedian = NA_real_,
                      z025quantile = NA_real_,
                      z075quantile = NA_real_,
                      z090quantile = NA_real_,
                      z095quantile = NA_real_,
                      z098quantile = NA_real_))
  
  heightmetrics = list(
    zmax = max(z), 
    zmean = mean(z),
    z_sd = sd(z),
    zmedian = median(z),
    z025quantile = quantile(z, 0.25),
    z075quantile = quantile(z, 0.75),
    z090quantile = quantile(z, 0.90),
    z095quantile = quantile(z, 0.95),
    z098quantile = quantile(z, 0.98)
  )
  return(heightmetrics)
}


#Plot metrics
set_lidr_threads(10)
plan(multisession, workers = 3)
res=5
radius = 15

# 1
# CExyz <- plot_metrics(ctg, ~canopy_entropy(X, Y, Z), filtplots[453,], radius = radius, fill = TRUE) #Ready
CExyz <- tryCatch(
  plot_metrics(ctg, ~canopy_entropy(X, Y, Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)

# plot(filtplots[['wsvha']],CExyz[[97]])
#2
# CV <- plot_metrics(ctg, ~metrics_custom(Height), filtplots, radius = radius, fill = TRUE) #Ready
CV <- tryCatch(
  plot_metrics(ctg,  ~metrics_custom(Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)

# nrow(CV[!is.na(CV$iqr),])
#3
# CR <- plot_metrics(ctg, ~rumple_index(X, Y, Height), filtplots, radius = radius, fill = TRUE) #Ready
CR <- tryCatch(
  plot_metrics(ctg,  ~rumple_index(X, Y, Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)

#4
# LADmet <- plot_metrics(ctg, ~metrics_lad(Height), filtplots, radius = radius, fill = TRUE)#Ready
LADmet <- tryCatch(
  plot_metrics(ctg,  ~metrics_lad(Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)

#5
# hmad <- plot_metrics(ctg, ~h_mad(Height),filtplots, radius = radius, fill = TRUE)#Ready
hmad <- tryCatch(
  plot_metrics(ctg,  ~h_mad(Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)


#6
# FrD <- plot_metrics(ctg, ~fd(X,Y,Height), filtplots, radius = radius, fill = TRUE) #ready
FrD <- tryCatch(
  plot_metrics(ctg,  ~fd(X,Y,Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)


#7
# CC <- plot_metrics(ctg, ~canopy_cover(X,Y,Height), filtplots, radius = radius, fill = TRUE) #ready
CC <- tryCatch(
  plot_metrics(ctg,  ~canopy_cover(X,Y,Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)


#8
# Hgt <- plot_metrics(ctg,~HeightMetrics(Height), filtplots, radius = radius, fill = TRUE) #Ready
Hgt <-  tryCatch(
  plot_metrics(ctg,  ~HeightMetrics(Z), filtplots, radius = radius, fill = TRUE),
  error = function(e) {
    message("Skipped empty plot: ", e$message)
    return(NULL)  # or return NA, if preferred
  }
)



### Merge all tables
#convert table into dataframes
FrDim_df <- as.data.frame(st_drop_geometry(FrD))
CExyz_df <- as.data.frame(st_drop_geometry(CExyz))
CV_df <- as.data.frame(st_drop_geometry(CV))
CR_df<- as.data.frame(st_drop_geometry(CR))
hmad_df<- as.data.frame(st_drop_geometry(hmad))
LAD_df<- as.data.frame(st_drop_geometry(LADmet))
CC_df<- as.data.frame(st_drop_geometry(CC))
Hgt_df<- as.data.frame(st_drop_geometry(Hgt))
nrow(CExyz_df)
#Find common names
comname<-names(CExyz_df)[names(CExyz_df)%in%names(CV_df)]

#merge tables
dat <- merge(CExyz_df, CV_df, by = comname, all = TRUE)

# dat<-cbind(CExyz_df,CV_df[,!(names(CV_df) %in% comname)])
dat<-merge(dat, FrDim_df, by = comname, all.x = TRUE)
dat <- merge(dat, CR_df, by = comname, all = TRUE)
# dat$CR<-CR_df[,!(names(CR_df)%in% comname)]
dat <- merge(dat, hmad_df, by = comname, all = TRUE)
# dat$hmad<-hmad_df[, !(names(hmad_df) %in% comname)]
dat <- merge(dat, LAD_df, by = comname, all = TRUE)
# dat<-cbind(dat,LAD_df[, !(names(LAD_df) %in% comname)])
dat <- merge(dat, CC_df, by = comname, all = TRUE)
# dat$CC<-CC_df[, !(names(CC_df) %in% comname)]
dat <- merge(dat, Hgt_df, by = comname, all = TRUE)
# dat<-cbind(dat,Hgt_df[, !(names(Hgt_df) %in% comname)])

#remove common columns
# dat <- dat[, !(names(dat) %in% comname)]

#Join with shape plots to create dataset
names(Hgt_df)
data<-dat
data<- merge(filtplots,dat, by = comname, all.y = TRUE)
st_write(data, "C2_plot_FSC_Full.shp",append=FALSE)
