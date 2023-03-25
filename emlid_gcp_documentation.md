# Emlid Ground Control Points Generator

This R script generates ground control points (GCPs) for a region of interest (ROI) using a digital elevation model (DEM). The script incorporates the `terra`, `sf`, `tidyverse`, `scales`, and `viridis` packages to process spatial data and produce GCPs based on elevation and coordinate information.
## Dependencies

Install and load the required packages:

```R
install.packages(c("terra", "sf", "tidyverse", "scales", "viridis"))

library(terra)
library(sf)
library(tidyverse)
library(scales)
library(viridis)
```


## Function: make_emlid_gcps

This function takes several input parameters, including the DEM, ROI, number of GCPs, coordinate reference system (CRS), buffer size, elevation weighting, and an optional plotting flag. It returns a data frame containing the generated GCPs, their coordinates, and elevations.

```R
make_emlid_gcps <- function(dem, roi, gcp_num, crs=4326, buffer=30, buffer_crs=26911, wt_elevation=1, plt=TRUE){ ... }
```


### Input Parameters 
- `dem`: Path to the digital elevation model file (e.g., a GeoTIFF file). 
- `roi`: Path to the region of interest polygon file (e.g., a KML or shapefile). 
- `gcp_num`: Target number of GCPs (must be >= 5). 
- `crs`: Coordinate reference system for the input data (default: 4326). 
- `buffer`: Buffer distance around the ROI polygon corners (default: 30 meters). 
- `buffer_crs`: CRS for buffer calculations (default: 26911). 
- `wt_elevation`: Weight of elevation in relation to X and Y coordinates (default: 1 for equal weighting). 
- `plt`: A boolean flag to plot the DEM and GCPs (default: TRUE).
### Output

The function returns a data frame containing the GCPs, their coordinates, and elevations.
## Example Usage

To use the `make_emlid_gcps` function, provide the paths to the DEM and ROI files, and set the desired parameters:

```R
# Set paths to DEM and ROI files
roi_path <- 'path/to/your/roi.kml'
ras_path <- 'path/to/your/dem.tif'

# Generate GCPs
emlid_gcps <- make_emlid_gcps(dem = ras_path,
                              roi = roi_path,
                              gcp_num = 10,
                              wt_elevation = 1)
```


## Save GCPs to a CSV File

To save the generated GCPs to a CSV file, use the `write.csv` function:

```R
output_file <- "gcp_for_emlid_flow.csv"
write.csv(emlid_gcps, output_file, row.names = F)
```



This will save the GCPs to a file named `gcp_for_emlid_flow.csv`. You can change the `output_file` variable to save the CSV file with a different name or location.