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
# 读取多波段栅格文件
for (year in 2023:2023) {
  raster_file <- paste0("F:/NTL/HDNTL/wlm/008-20464-Beijing-",year,".tif")  # 替换为你的文件路径
  r <- stack(raster_file)  # 使用raster包
  # 或者使用terra包
  # r <- rast(raster_file)
  cols <- r@ncols
  rows <- r@nrows
  # 获取波段名称
  band_names <- names(r)
  
  r_mat <- as.matrix(r)
  # 从波段名称中提取日期
  dates <- sub("NTL", "", band_names)
  dates <- as.Date(dates, format="%Y_%m_%d")
  
  # 按年份分组检查
  years <- unique(format(dates, "%Y"))
  # 获取该年的所有日期
  year_dates <- dates[format(dates, "%Y") == year]

  # 创建该年应有的完整日期序列
  complete_dates <- seq(as.Date(paste0(year, "-01-01")),
                        as.Date(paste0(year, "-12-31")),
                        by = "day")

  # 检查是否有缺失日期
  missing_dates <- complete_dates[!complete_dates %in% year_dates]

  # 输出结果
  cat("\n年份:", year)
  cat("\n实际天数:", length(year_dates))
  cat("\n应有天数:", length(complete_dates))

  if(length(missing_dates) > 0) {
    cat("\n缺失的日期:\n")
    print(format(missing_dates, "%Y_%m_%d"))
  } else {
    cat("\n数据完整\n")
  }
  
    # 如果需要补全数据，可以创建一个完整的日期序列
  missing_all <- complete_dates[!complete_dates %in% dates]
  complete_names <- paste0("NTL", format(complete_dates, "%Y_%m_%d"))
  if(length(missing_all) > 0) {
      # 创建新的完整波段名称列表
  
  
    # 创建新的完整矩阵
    matrix_complete <- matrix(0, nrow = nrow(r_mat),
                              ncol = length(complete_dates))
  
    # 将原始数据填入对应位置
    existing_indices <- match(names(r), complete_names)
    matrix_complete[, existing_indices] <- r_mat
  
    # 更新矩阵
    r_mat <- matrix_complete
  
    # 检查维度
    dim(r_mat)
  }
  # 循环遍历波段，输出每日tif图像
  for (i in 1:length(complete_names)) {
    # 获取当前波段名称
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
                            crs = CRS("+proj=longlat +datum=WGS84"))#定义空栅格
    values(traits_output) <- as.numeric(traits_m)
    writeRaster(traits_output,paste0("F:/NTL/HDNTL/Beijing","/",year,"/",band_name,".tif"),overwrite=T)
  }
}






path <- "F:/NTL/HDNTL/NewYork/2022"  # 替换为你的文件路径
setwd(path)
NTL <- list.files(pattern = '.tif$', full.names = TRUE)
NTL <- rast(NTL)
r <- stack(NTL)
dim_ntl <- dim(NTL)
row_ntl <- dim_ntl[1] * dim_ntl[2]
col_ntl <- dim_ntl[3]
NTL_matrix <- matrix(data = values(NTL), nrow = row_ntl, ncol = col_ntl)

NTL_matrix[NTL_matrix == 0] <- NA


# 初始化一个矩阵来存储年均值
annual_NTL <- matrix(data = NA, ncol = row_ntl, nrow = 1)

#计算每年的平均值
  # 计算每年的平均值
start_index <- 1
end_index <- 365  # 每年结束的索引
if (end_index > col_ntl) {
  end_index <- col_ntl  # 防止越界
}

#计算该年的平均值
for (j in 1:row_ntl) {
  if(all(is.na(NTL_matrix[j,])) | mean(NTL_matrix[j,], na.rm = TRUE) < 100){
    next
  }else{
    annual_NTL[1, j] <- mean(NTL_matrix[j,], na.rm = TRUE)
  }
}

# 循环遍历年均值矩阵，并输出为栅格文件
#获取每年的平均值
traits_m <- annual_NTL[1, ]

# 创建一个空栅格
traits_output <- rast(nrows = dim_ntl[1], ncols = dim_ntl[2],
                      xmin = ext(NTL)$xmin, xmax = ext(NTL)$xmax,
                      ymin = ext(NTL)$ymin, ymax = ext(NTL)$ymax,
                      crs = crs(NTL))  

# 将年均值填入栅格
values(traits_output) <- traits_m

# 写入TIF文件
writeRaster(traits_output,
            filename = paste0("F:/NTL/HDNTL/meanNTL_NewYork2022.tif"),
            overwrite = TRUE)






library(terra)

# 读取数据
building_height <- rast("F:/NTL/building_height/buildingheight1.tif")  # 替换为你的建筑高度数据路径
ntl_500m <- rast("F:/NTL/HDNTL/meanNTL_BJ2024_mask.tif")

# 裁剪建筑高度数据到NTL范围
building_cropped <- crop(building_height, ntl_500m)

building_resampled <- terra::resample(
  building_cropped,
  ntl_500m,
  method = "average"  # 计算均值
)

# 保存结果
writeRaster(building_resampled, "F:/NTL/building_height/building_height_500m_mean.tif", overwrite = TRUE)
