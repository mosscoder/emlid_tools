library(terra)
library(sf)
library(tidyverse)
library(scales)

select <- dplyr::select

make_emlid_gcps <- function(dem, roi, gcp_num, crs=4326, buffer=30, buffer_crs=26911, wt_elevation=1, plt=TRUE){
  if(gcp_num < 5) return('You must specify 5 or more GCPs')
  
  roi_raw <- st_read(roi) %>% st_transform(crs)
  
  el_rast <- rast(dem)
  targ_ras_crs <- paste0('epsg:',crs)
  
  if(crs(el_rast) != targ_ras_crs){
    dem_crs_poly <- roi_raw %>%
      st_transform(buffer_crs) %>% 
      st_buffer(100) %>% 
      st_transform(crs = st_crs(el_rast))
    el_rast <- project(el_rast %>% crop(dem_crs_poly), targ_ras_crs)
  }
  
  roi_poly <- roi_raw %>% 
    st_transform(buffer_crs) %>% #probably need more elegant soltn. here
    st_buffer(-buffer) %>% 
    st_transform(crs)
  
  el_roi <- el_rast %>% crop(roi_poly) %>% mask(vect(roi_poly))
  
  el_roi_plt <- el_rast %>% crop(roi_raw)
  
  if(crs == 4326){
    x_name <- 'longitude'
    y_name <- 'latitude'
  } else {
    x_name <- 'easting'
    y_name <- 'northing'
  }
  
  corner_gcps <- roi_poly %>% 
    st_convex_hull() %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    select(X, Y) %>% 
    head(-1) %>% 
    arrange(-Y, X) %>% 
    mutate(name = paste0('corner_', seq_len(nrow(.)))) %>% 
    rename({{x_name}} := X, {{y_name}} := Y) %>% 
    select(name, everything()) %>% 
    mutate(elevation = terra::extract(el_rast, .[,2:3])[,2] %>% unlist())
  
  if(crs == 4326){
    corner_gcps <- corner_gcps %>% 
      rename(`ellipsoidal height` = elevation)
  }
  
  corner_buff <- corner_gcps %>% st_as_sf(coords = c(x_name, y_name), crs=crs) %>% 
    st_transform(buffer_crs) %>% #probably need more elegant soltn. here
    st_buffer(buffer) %>% 
    st_transform(crs)
  
  el_nested <- mask(el_roi, corner_buff %>% vect(), inverse = T)
  el_vals <- values(el_nested)
  area_check <- nrow(na.omit(el_vals))
  
  if(area_check > 0){
    xys <- xyFromCell(el_nested, seq_len(ncell(el_nested)))
    xyzs <- cbind(xys, el_vals)
    xyz_scaled <- scale(na.omit(xyzs))
    xyz_scaled[,3] <- xyz_scaled[,3]*wt_elevation 
    
    set.seed(123)
    el_grad <- kmeans(xyz_scaled, centers = gcp_num - 4 )$centers
    
    for(axis in 1:3){
      focal_axis_mean <- attributes(xyz_scaled)$`scaled:center`[axis]
      focal_axis_sd <- attributes(xyz_scaled)$`scaled:scale`[axis]
      el_grad[,axis] <- (el_grad[,axis] * focal_axis_sd) + focal_axis_mean
    }
    
    el_grad <- el_grad %>% 
      as.data.frame() %>% 
      mutate(name = paste0('strat_', round(.[,3]))) %>% 
      select(name, everything()) %>% 
      arrange(-y, x)
    
  } else {
    print('No viable area outside of corner buffers. Returning ROI center.')
    center_col <- round(ncol(el_roi)/2)
    center_row <- round(nrow(el_roi)/2)
    center_x <- xFromCol(el_roi, center_col)
    center_y <- yFromRow(el_roi, center_row)
    center_cell <- cellFromXY(el_roi, cbind(center_x, center_y))
    center_elevation <- terra::extract(el_roi, center_cell)[1] %>% unlist()
    el_grad <- data.frame(name = 'center',
                          x = center_x,
                          y = center_y,
                          elevation = center_elevation) 
  }
  
  colnames(el_grad) <- colnames(corner_gcps)
  out <- rbind(corner_gcps, el_grad)
  
  gcp_vect <- out %>% st_as_sf(coords=c(x_name,y_name),crs=crs) %>% vect()
  poly_vect <- roi_raw %>% vect()
  
  if(isTRUE(plt)){
    plot_cols_alpha <- alpha(viridis::viridis(100), 0.3)
    slope <- terrain(el_roi_plt, "slope", unit="radians")
    aspect <- terrain(el_roi_plt, "aspect", unit="radians")
    hill <- shade(slope, aspect, 10, 270)
    plot(hill, col=gray.colors(100), legend=FALSE, axes=FALSE)
    plot(el_roi_plt, col=plot_cols_alpha, legend=T, axes=FALSE, add = T)
    plot(gcp_vect, add = T, pch = 13)
    plot(poly_vect, border = 'red', add = T, legend=FALSE, axes=FALSE)
  }
  
  return(out)
}

# setwd('~/Downloads')
# roi_path <- 'topHouse_SW.kml'
# ras_path <- './46114f1/MISSOULA_2019_ClrkFrkBttrtRvr/HFDEM/46114f1_HFDEM.tif'
# 
# emlid_gcps <- make_emlid_gcps(dem = ras_path, #path to elevation model
#                               roi = roi_path, #path to polygon, could be any driver sf accepts
#                               gcp_num = 10, #target number of GCPs, must be >= 5
#                               wt_elevation = 1) #weight of elevation in relation to X and Y, 1 is equal weighting

#write.csv(emlid_gcps, 'gcp_for_emlid_flow.csv', row.names = F) #saves gcps in format for Emlid Flow