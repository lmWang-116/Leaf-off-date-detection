# Required packages
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
library(signal)  
library(Hmisc) 

# autumn phenology from daytime remote sensing data (e.g., EVI)
# set EVI data path
path <- "your path"  
setwd(path)
EVI <- list.files(pattern = '.tif$', full.names = TRUE)
EVI <- rast(EVI)
r <- stack(EVI)

dim_EVI <- dim(EVI)
row_EVI <- dim_EVI[1] * dim_EVI[2]
col_EVI <- dim_EVI[3]
EVI_matrix <- matrix(data = EVI*0.0001, nrow = row_EVI, ncol = col_EVI)


# 3-point max filter 
max_filter <- function(ts, window = 3) {
  rollapply(ts, width = window, 
            FUN = function(x) {
              if (all(is.na(x))) {
                return(NA)
              } else {
                return(max(x, na.rm = TRUE))
              }
            }, 
            fill = NA, align = "center", partial = TRUE)
}



# Function to process EVI time-series using adaptive Savitzky-Golay method
process_EVI_timeseries_adaptive_sg <- function(evi_matrix, time_intervals = 16, sg_m_range = 4:7, sg_d_range = 2:4, sg_m_fit = 4, sg_d_fit = 6, max_iter = 20) {
  # time_intervals: days between observations, default 8 for MOD09A1
  # sg_m_range: range for initial trend half-window
  # sg_d_range: range for initial polynomial degree
  # sg_m_fit: fixed m for iteration fitting (recommended 4)
  # sg_d_fit: fixed d for iteration fitting (recommended 6)
  # max_iter: maximum iterations to prevent infinite loop
  
  smoothed_matrix <- matrix(NA, nrow = nrow(evi_matrix), ncol = ncol(evi_matrix))
  
  n_time <- ncol(evi_matrix)
  
  for (i in 1:nrow(evi_matrix)) {
    print(i)
    pixel_ts <- evi_matrix[i, ]
    pixel_ts[pixel_ts < -0.05] <- NA
    pixel_ts <- max_filter(pixel_ts)
    pixel_ts <- c(pixel_ts[(length(pixel_ts)-1):length(pixel_ts)],pixel_ts,pixel_ts[1:2])
    
    if (!all(is.na(pixel_ts))) {
      # Step 1: Interpolate NA values by connecting flanking points with linear fit
      ts_original <- pixel_ts
      ts_filled <- ts_original  # Initialize ts_filled
      na_indices <- which(is.na(ts_original))
      
      if (length(na_indices) > 0) {
        for (na_idx in na_indices) {
          # Find the nearest non-NA points before and after
          left_idx <- max(which(!is.na(ts_original[1:(na_idx-1)])), 0)
          right_idx <- min(which(!is.na(ts_original[(na_idx+1):n_time])), na.rm = TRUE) + na_idx
          if (left_idx == 0 || right_idx > n_time) next  # Skip if no flanking points
          
          left_val <- ts_original[left_idx]
          right_val <- ts_original[right_idx]
          # Linear interpolation between left_idx and right_idx
          steps <- right_idx - left_idx
          slope <- (right_val - left_val) / steps
          for (j in 1:(steps-1)) {
            ts_filled[left_idx + j] <- left_val + slope * j
          }
        }
      }
      
      # Fill any remaining NAs at edges with LOCF
      ts_filled <- na.locf(ts_filled, na.rm = FALSE)
      ts_filled <- na.locf(ts_filled, fromLast = TRUE)
      
      N0 <- ts_filled  # N_i^0
      
      # Step 2: Fit initial long-term trend with optimal SG parameters (using periodic filter)
      best_sse <- Inf
      best_trend <- NULL
      for (m in sg_m_range) {
        n <- 2 * m + 1
        if (n > n_time) next  # Skip if window too large
        for (d in sg_d_range) {
          trend <- sgolayfilt(N0, p = d, n = n)
          sse <- sum((N0 - trend)^2, na.rm = TRUE)
          if (sse < best_sse) {
            best_sse <- sse
            best_trend <- trend
          }
        }
      }
      if (is.null(best_trend)) next  # If no valid fit, skip
      N_tr <- best_trend
      
      # Step 3: Determine weights W_i
      diffs <- abs(N0 - N_tr)
      d_max <- max(diffs, na.rm = TRUE)
      W <- ifelse(N0 >= N_tr, 1, ifelse(d_max > 0, 1 - (diffs / d_max), 1))
      
      # Initial envelope adjustment
      N_new <- ifelse(N0 >= N_tr, N0, N_tr)
      
      # Step 5: First fitting with SG (using periodic filter, m=4, d=6)
      n_fit <- 2 * sg_m_fit + 1
      if (n_fit > n_time) next  # Skip if window too large
      fitted <- sgolayfilt(N_new, p = sg_d_fit, n = n_fit)
      
      # Compute first F_k (k=1)
      F_values <- c()
      current_fitted <- fitted
      F_k <- sum(abs(current_fitted - N0) * W, na.rm = TRUE)
      F_values <- c(F_values, F_k)
      
      # Iteration
      converged <- FALSE
      iter <- 1
      while (!converged && iter < max_iter) {
        N_new <- ifelse(N0 >= current_fitted, N0, current_fitted)
        new_fitted <- sgolayfilt(N_new, p = sg_d_fit, n = n_fit)
        new_F <- sum(abs(new_fitted - N0) * W, na.rm = TRUE)
        F_values <- c(F_values, new_F)
        
        # Check stopping condition: F_{k-1} >= F_k <= F_{k+1}
        k <- length(F_values) - 1  # Current F_k index
        if (k >= 2 && F_values[k-1] >= F_values[k] && F_values[k] <= F_values[k+1]) {
          # Minimum at F_k, use the corresponding fitted
          smoothed_matrix[i, ] <- current_fitted[3:25]  # Fitted for F_k
          converged <- TRUE
        } else {
          current_fitted <- new_fitted
          iter <- iter + 1
        }
      }
      
      # If max iter reached, use last fitted
      if (!converged) {
        smoothed_matrix[i, ] <- current_fitted[3:25]
      }
    }
  }
  
  return(list(smoothed = smoothed_matrix))
}

# Example usage
results <- process_EVI_timeseries_adaptive_sg(EVI_matrix)
smoothed_rast <- results$smoothed

# logical function of autumn phenology
single_logistic <- function(t, a,b,c,d) {
  d + c / (1+exp(a+b*t))
}

#where a and b are fitting coefficients;
#c + d is the maximum value of a given autumn leaf phenology curve; d is the minimum value of the given autumn leaf phenology curve.
# perform logical fitting and derivative calculation
process_pixel <- function(evi,i) {
  initial_guess <- c(a = -15, b = 0.05, c = 0.2, d = 0.1)
  if(max(evi) > 0.1){
    doy <- 1 + (12:(length(evi)-1)) * 16

    # Define the parameter range
    lower_bounds <- c(-50,0.02,0.02,0.05)  # min
    upper_bounds <- c(-5,0.25,0.7,0.3) # max

    tryCatch({ 
      fit <- nls(evi[13:(length(evi))] ~ single_logistic(doy, a,b,c,d),
                 start = initial_guess,
                 lower = lower_bounds,
                 upper = upper_bounds,
                 algorithm = "port") 

      # Obtain the fitting parameters
      coefs <- coef(fit)
      a <- coefs["a"]
      b <- coefs["b"]
      c <- coefs["c"]
      d <- coefs["d"]

      # Generate a smooth fitting curve
      doy_fine <- seq(193,366,1)
      fitted_evi <- single_logistic(doy_fine, a,b,c,d)

      evi_max <- max(fitted_evi, na.rm = TRUE)
      evi_min <- min(fitted_evi, na.rm = TRUE)

      # Calculate the EVI ratio (SVR)
      evi_ratio <- (evi_max / abs(evi_min))

      if(evi_ratio > 1.5){
        first_derivative <- diff(fitted_evi) / diff(doy_fine)
        first_derivative <- c(first_derivative[1], first_derivative)

        second_derivative <- diff(first_derivative) / diff(doy_fine)
        second_derivative <- c(second_derivative[1], second_derivative)

        third_derivative <- diff(second_derivative) / diff(doy_fine)
        third_derivative <- c(third_derivative[1], third_derivative)

        # the extreme values of the third derivative
        minima_indices <- which(diff(sign(diff(third_derivative))) == 2) + 1  # min
        maxima_indices <- which(diff(sign(diff(third_derivative))) == -2) + 1 # max

        minima_doy <- doy_fine[minima_indices]
        maxima_doy <- doy_fine[maxima_indices]


        # middle of fall MOF
        mof_candidates <- maxima_doy[abs(maxima_doy - 300) < 60]
        mof_doy <- mof_candidates[which.max(mof_candidates)]

        # start of fall SOF
        if (!is.na(mof_doy)) {
          sof_candidates <- minima_doy[minima_doy < mof_doy]
          sof_doy <- ifelse(length(sof_candidates) > 0,
                            max(sof_candidates),
                            NA)
        } else {
          sof_doy <- NA
        }

        # end of fall EOF
        if (!is.na(mof_doy)) {
          eof_candidates <- minima_doy[minima_doy > mof_doy]
          if (length(eof_candidates) > 0 && min(eof_candidates) < mof_doy + 80 && eof_candidates <= 366) {
            eof_doy <- min(eof_candidates)
          } else {
            eof_doy <- NA
          }
        } else {
          eof_doy <- NA
        }
      }else{
        doy_fine <- NA
        fitted_evi <- NA
        third_derivative <- NA
        sof_doy <- NA
        mof_doy <- NA
        eof_doy <- NA
      }
      # return list
      list(
        doy_fine = doy_fine,
        fitted_evi = fitted_evi,
        third_derivative = third_derivative,
        sof_doy = sof_doy,
        mof_doy = mof_doy,
        eof_doy = eof_doy
      )

    }, error = function(e) {
      cat("Fitting failed:",i, conditionMessage(e), "\n")
      list(
        doy_fine = NA,
        fitted_evi = NA,
        third_derivative = NA,
        sof_doy = NA,
        mof_doy = NA,
        eof_doy = NA
      )
    })
  }
}




# autumn phenology derived from EVI
SOF <- matrix(NA, nrow = row_EVI, ncol = 1)
EOF <- matrix(NA, nrow = row_EVI, ncol = 1)
MOF <- matrix(NA, nrow = row_EVI, ncol = 1)
for (i in 1:row_EVI) {
  if(all(is.na(smoothed_rast[i,]))){
    next 
  }else{
    results <- process_pixel(smoothed_rast[i,],i)
    if(!is.null(results$eof_doy)){
      SOF[i,1] <- results$sof_doy
      EOF[i,1] <- results$eof_doy
      MOF[i,1] <- results$mof_doy
    }
  }
}

hist(SOF)
hist(EOF)
hist(MOF)







# SVR > 1.5

deciduous <- matrix(NA, nrow = row_EVI, ncol = 1)
for (i in 1:row_EVI) {
  if(all(is.na(smoothed_rast[i,]))){
    next 
  }else{
    if(smoothed_rast[i,13]/smoothed_rast[i,23]>1.5){
      deciduous[i,1] <- 1 
    }
  }
}


# determine whether the fitted curve conforms to the S-shaped feature, based on the second-order derivative shape.
is_s_curve_by_second_derivative <- function(x, t, min_peak_distance = 3) {

  if (length(x) != length(t)) {
    stop("wrong")
  }

  dx1 <- diff(x)
  dx2 <- diff(dx1)
  
  # Maximum and minimum points
  n <- min_peak_distance  # Window size to prevent adjacent fluctuations from being misjudged as extreme values
  is_maxima <- logical(length(dx2))
  is_minima <- logical(length(dx2))
  
  for (i in (n+1):(length(dx2)-n)) {
    current_val <- dx2[i]
    window_before <- dx2[(i-n):(i-1)]
    window_after <- dx2[(i+1):(i+n)]
    
    # If the current value is greater than all the values in the previous and subsequent Windows, it is a maximum value
    is_maxima[i] <- all(current_val > window_before) && all(current_val > window_after)
    
    # If the current value is less than all the values in the previous and subsequent Windows, it is a minimum value
    is_minima[i] <- all(current_val < window_before) && all(current_val < window_after)
  }
  
  # Obtain the index of the extreme points
  maxima_idx <- which(is_maxima)
  minima_idx <- which(is_minima)
  
  # Merge and sort the extreme points by time
  all_extrema <- data.frame(
    idx = c(maxima_idx, minima_idx),
    type = c(rep("max", length(maxima_idx)), rep("min", length(minima_idx))),
    value = c(dx2[maxima_idx], dx2[minima_idx])
  )
  
  if (nrow(all_extrema) > 0) {
    all_extrema <- all_extrema[order(all_extrema$idx), ]
    
    if (nrow(all_extrema) >= 2) {
      first_two_types <- all_extrema$type[1:2]
      
      return(all(first_two_types == c("max", "min")))
    }
  }
  
  return(FALSE) 
}



# detect leaf-off date from NTL data
NTL_path <- "your path"
path <- NTL_path
setwd(path)
NTL <- list.files(pattern = '.tif$', full.names = TRUE)
NTL <- stack(NTL)
dim_ntl <- dim(NTL)
row_ntl <- dim_ntl[1] * dim_ntl[2]
col_ntl <- dim_ntl[3]
NTL_matrix <- matrix(data = values(NTL)*0.01, nrow = row_ntl, ncol = col_ntl)

# Preprocess NTL_matrix and fill the NA values in rows that are not all NA to 0
preprocess_na <- function(matrix) {
  for (i in 1:nrow(matrix)) {
    row_data <- matrix[i, ]
    if (!all(is.na(row_data))) { 
      row_data[is.na(row_data)] <- 0 
      matrix[i, ] <- row_data
    }
  }
  return(matrix)
}

# remove outlier value
for (i in 1:nrow(NTL_matrix)) {
  threshold <- quantile(NTL_matrix[i, ], prob = 0.99, na.rm = TRUE)
  NTL_matrix[i, NTL_matrix[i, ] > threshold] <- 0
}


# Preprocess NTL_matrix
NTL_matrix <- preprocess_na(NTL_matrix)


# set the year you are interested
start_date <- c("2024-01-01")
dates <- seq(as.Date(start_date), by = "day", length.out = ncol(NTL_matrix))
row_ntl <- nrow(NTL_matrix)

# set the pixel index (here 70981 in Beijing as an example)
i= 70981 
ntls <- NTL_matrix[i,] 
year <- 2024
span = 0.8

# Step 1: Extract the data of the year
year_mask <- year(dates) == year
year_dates <- dates[year_mask]
year_ntl <- ntls[year_mask]

# Step 2: Extract season-end data (MOF-10 to EOF+30 was set)
if(is.na(EOF[i,]) & is.na(MOF[i,])){
  season_mask <- c((round(mean(MOF, na.rm = TRUE))-10):min(round(mean(EOF, na.rm = TRUE)) + 30, 366))
}else if(!is.na(EOF[i,]) & !is.na(MOF[i,])){
  season_mask <- c((MOF[i,]-10):min(EOF[i,] + 30, 366))
}else if(is.na(MOF[i,]) & !is.na(EOF[i,])){
  season_mask <- c((round(mean(MOF, na.rm = TRUE))-10):min(EOF[i,] + 30, 366))
}else if(!is.na(MOF[i,]) & is.na(EOF[i,])){
  season_mask <- c((MOF[i,])-10:min(round(mean(EOF, na.rm = TRUE)) + 30, 366))
}

# remove outlier value
s_dates <- year_dates[season_mask]
s_ntl <- year_ntl[season_mask]
s_ntl[s_ntl == 0] <- NA
q5 <- quantile(s_ntl[!is.na(s_ntl)], 0.05)
q95 <- quantile(s_ntl[!is.na(s_ntl)], 0.95)

s_ntl[s_ntl<q5 | s_ntl>q95] <- NA


# Step 3: three days maximum composition
window_size <- 3
max_ntl <- s_ntl  

max_ntl <- rollapply(s_ntl,
                     width = window_size,
                     FUN = function(x) {
                       if(all(is.na(x))) NA else max(x, na.rm = TRUE)
                     },
                     fill = NA,
                     align = "left")

max_ntl <- c(max_ntl, rep(NA, window_size - 1))[1:length(s_ntl)]
s_ntl <- max_ntl  


# Initialize the result
result <- c(NA, NA)

# Step 4: Effective observation point
valid_idx <- which(!is.na(s_ntl))
if (length(valid_idx) >= 0.4*length(season_mask)) {
  t <- yday(s_dates[valid_idx])
  y <- s_ntl[valid_idx]
  t_full <- yday(s_dates)
  
  # Step 5: LOESS regression
  loess_fit <- tryCatch({
    model <- loess(s_ntl ~ t_full, span = span,family = c("symmetric"))#
    predict(model, newdata = data.frame(t = t_full))
  }, error = function(e) NULL)
  
  # Step 6: logistic fitting
  if (!is.null(loess_fit)) {
    compressed_loess <- pmin(pmax(loess_fit, min(s_ntl, na.rm = TRUE)), 
                             max(s_ntl, na.rm = TRUE))
    
    # uniformization
    normalized_loess <- (compressed_loess - min(compressed_loess, na.rm = TRUE)) / 
      (max(compressed_loess, na.rm = TRUE) - min(compressed_loess, na.rm = TRUE) + 1e-5)
    
    # fitting
    logistic_model <- tryCatch({
      y_logit <- log((normalized_loess + 1e-5) / (1 - normalized_loess + 1e-5))
      glm(y_logit ~ t_full, family = gaussian())
    }, error = function(e) {
      message("Logistic fail:", e$message)
      NULL
    })
    
    # Predict and revert to the original data scale
    if (!is.null(logistic_model)) {
      logit_pred <- predict(logistic_model, newdata = data.frame(t = t_full))
      logistic_pred <- 1 / (1 + exp(-logit_pred))
      min_val <- min(compressed_loess, na.rm = TRUE)
      max_val <- max(compressed_loess, na.rm = TRUE)
      original_scale_pred <- logistic_pred * (max_val - min_val) + min_val
    }
  }
  
  
  if (!is.null(logistic_model) && is_s_curve_by_second_derivative(original_scale_pred, t_full) &&
      rcorr(cbind(loess_fit, original_scale_pred), type = "pearson")$r[1,2] > 0.7) {
    # Check whether there is sufficient data at the beginning and end of the fitting curve
    if (!(all(is.na(logistic_pred[1:10])) | all(is.na(logistic_pred[(length(season_mask)-10):length(season_mask)])))) {
      date <- which.max(diff(original_scale_pred))
      
      idx_first <- date
      
      leaf_off <- t_full[1] + idx_first
      result <- c(year, leaf_off)
    }
  }
}

result

if(!is.na(EOF[i,])){
  EOFs <- EOF[i,]
}else{
  EOFs <- round(mean(EOF,na.rm=T))
}


# plot
ggplot(data = data.frame(t_full = t_full, s_ntl = s_ntl), aes(x = t_full)) +
  # original data points
  geom_point(aes(y = s_ntl), alpha = 0.6, color = "gray30") +
  # LOESS regression 
  geom_line(aes(y = loess_fit), color = "black", linewidth = 0.8,linetype="dashed",alpha = 0.5) +
  # logistic fitting
  geom_line(aes(y = original_scale_pred), color = "#cc4400", linewidth = 1) +
  geom_vline(xintercept = leaf_off, color = "#cc4400", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = EOFs, color = "steelblue", linetype = "dashed", linewidth = 1) +
  annotate("text", 
           x = mean(t_full)+10, 
           y = max(s_ntl, na.rm = TRUE)*0.98,  
           label = paste("DOY=", result[2]),  
           color = "black", 
           size = 6) +  
  annotate("text", 
           x = mean(t_full),  
           y = max(s_ntl, na.rm = TRUE)*0.98,  
           label = paste("EOF", EOFs), 
           color = "black", 
           size = 6) +  
  # scale_y_continuous(breaks = seq(3,6,0.5))+
  scale_x_continuous(breaks = seq(300,360,10))+
  labs(x = "DOY", y = "NTL") +
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'right',
    text = element_text(family = "Arial")
  ) 

output_result_path <- "your path"
ggsave(paste0(output_result_path,".pdf"),device = cairo_pdf,
       width = 6,height = 4)
dev.off()


