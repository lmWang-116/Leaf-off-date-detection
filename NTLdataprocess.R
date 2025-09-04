library(raster)
library(zoo)
library(terra)  
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate)
library(forecast)
library(lubridate)
library(stringr)
library(geoTS)
library(mgcv)
# Read multi-band raster files
years <- c("your years")
input_path <- "your path"
output_path <- "your path"
for (year in years) {
  raster_file <- paste0(input_path,year,".tif") 
  r <- stack(raster_file) 
  cols <- r@ncols
  rows <- r@nrows
  band_names <- names(r)
  
  r_mat <- as.matrix(r)
  # get date
  dates <- sub("NTL", "", band_names)
  dates <- as.Date(dates, format="%Y_%m_%d")
  
  # groups by year
  years <- unique(format(dates, "%Y"))
  # Get all the dates of that year
  year_dates <- dates[format(dates, "%Y") == year]

  # Create the complete date sequence of the year 
  complete_dates <- seq(as.Date(paste0(year, "-01-01")),
                        as.Date(paste0(year, "-12-31")),
                        by = "day")

  # Check if there are any missing dates
  missing_dates <- complete_dates[!complete_dates %in% year_dates]

  # output
  cat("\n year:", year)
  cat("\n actual number of days:", length(year_dates))
  cat("\n totoal days:", length(complete_dates))

  if(length(missing_dates) > 0) {
    cat("\n missing dates:\n")
    print(format(missing_dates, "%Y_%m_%d"))
  } else {
    cat("\n complete\n")
  }
  
  # Complete the missing
  missing_all <- complete_dates[!complete_dates %in% dates]
  complete_names <- paste0("NTL", format(complete_dates, "%Y_%m_%d"))
  if(length(missing_all) > 0) {
    matrix_complete <- matrix(0, nrow = nrow(r_mat),
                              ncol = length(complete_dates))
    existing_indices <- match(names(r), complete_names)
    matrix_complete[, existing_indices] <- r_mat
  
    r_mat <- matrix_complete
    dim(r_mat)
  }
  # Output daily tif images
  for (i in 1:length(complete_names)) {
    band_name <- complete_names[i]
    j=1
    traits_m <- matrix(data = NA,ncol = cols,nrow = rows)
    for (k in seq(1,cols*rows,rows)) {
      traits_m[1:rows,j] <- r_mat[k:(k+rows-1),i]
      j = j + 1
    }
    r <- stack(raster_file)
    traits_output <- raster(traits_m,
                            xmn=r@extent@xmin,
                            xmx=r@extent@xmax,
                            ymn=r@extent@ymin,
                            ymx=r@extent@ymax,
                            crs = CRS("+proj=longlat +datum=WGS84"))
    values(traits_output) <- as.numeric(traits_m)
    writeRaster(traits_output,paste0(output_path,year,"/",band_name,".tif"),overwrite=T)
  }
}

