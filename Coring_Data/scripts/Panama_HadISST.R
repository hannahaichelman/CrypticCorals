################################################################################
####### Panama Annual and Summer Mean SST CSV Construction #####################
####### Code written by Colleen B. Bove ########################################
################################################################################

# ------------------
#### Load libraries

### Used packages that need to be installed to run code:
needed_packages <- c("ncdf4", "raster", "tidyverse", "xts") # Specify necessary packages

not_installed <- needed_packages[!(needed_packages %in% installed.packages()[ , "Package"])] # Extract not installed packages
if(length(not_installed)) install.packages(not_installed) # Install not installed packages

library(ncdf4)
library(raster) 
library(tidyverse)
library(xts)


# ------------------
#### Other settings/defaults to set

## spatial subsetting for Caribbean region
Xmin <- -100
Xmax <- -55
Ymin <- 0
Ymax <- 40

## set the projection
proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

## Create base map for SST maps
map_base <- ggplot2::fortify(maps::map(fill = TRUE, plot = FALSE)) %>% 
  dplyr::rename(lon = long)


# ------------------
### setting file path

## HadISST_sst.nc: 
wd <- getwd()
HadISST_sst_path <- paste0(wd, "/HadISST_sst.nc")

## The current HadISST file was downloaded on 31 January 2024 from the Met Office website:https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html


# ------------------  
#### SST Extraction and Analysis (using HadISST data)

## read in the HadISST raster
had_data <- HadISST_sst_path # this is from the path defined at the start
data_brick <- brick(had_data) # read SST data as rasterbrick
data_brick <- crop(data_brick, extent(Xmin, Xmax, Ymin, Ymax)) # clip rasterbrick for Caribbean region only (see 'spatial subsetting' above)
plot(data_brick[[1:5]]) # plot a subset of the layers in the brick to check it worked

ncdf_file <- nc_open(had_data) # open up the netCDF

nc_lats <- ncdf_file$dim$latitude$vals #other way of looking into the data structure; extracting lats (HadISST)
nc_longs <- ncdf_file$dim$lon$vals #longs here (HadISST)
nc_longs[nc_longs>180] <- nc_longs[nc_longs>180]-360 # convert lat/lon

#time handling
raw_times <- seq(0, (length(ncdf_file$dim$time$vals)-1), by = 1) # total number of layers (these are monthly observations)
nc_times <- ymd("1870-01-01") + months(raw_times)
nc_close(ncdf_file) #close the file


## Add assumed dates or raster
data_brick <- setZ(data_brick, nc_times)
names(data_brick) <- nc_times

## create a dataframe for sampling locations
reef_gps <- data.frame("site" = c("CI",	"PD",	"SP",	"BS",	"CA",	"BN"),
                       "lon" =	c(-82.24351143,	-82.36403667,	-82.26453,	-82.092315,	-82.053951,	-82.176505),
                       "lat" =	c(9.265448571,	9.36924,	9.352434,	9.287438,	9.193858,	9.348495))

## convert sampling points to spatial points and applying projection to points
LongLat <- cbind(reef_gps$lon, reef_gps$lat) 
ReefPTS <- SpatialPoints(LongLat)
proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
projection(ReefPTS) <- proj

## extract data from the sst brick with the sampling points shifted
had_sst <- raster::extract(data_brick, ReefPTS)

## convert the data matrix to a final dataframe
had_sst_xts <- xts(t(had_sst), nc_times) # t is used to transpose the matrix because xts assumes dates are along the lines, not columns
had_sst_long <- data.frame(had_sst_xts)
had_sst_long$date <- rownames(had_sst_long)


# ------------------  
#### Final data manipulation for annual and summer mean SST

## Calculate annual mean SST for each year
year_sst_mean <- data.frame(apply.yearly(had_sst_xts, colMeans, na.rm = TRUE)) %>% # mean SST per year
  dplyr::select(X1) %>% # all sites are within the grid so just proceeding with one GPS coordinate
  rename(`Annual Mean SST` = X1) %>% 
  rownames_to_column("date") %>% 
  separate(date, c("Year", "month", "day")) %>% 
  dplyr::select(-month, -day)


## Function to identify summer months (here set from August to November) 
findSummer <- function(DATES) {
  Aug <- as.Date("2012-8-01",  format = "%Y-%m-%d") # August
  Nov <- as.Date("2012-11-30",  format = "%Y-%m-%d") # November
  
  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= Aug & d < Nov, "Yes", "No")
}


## Calculate mean summer SST for each year
had_sst_seasonal <- had_sst_long %>% 
  separate(date, c("Year", "month", "day")) %>% 
  rownames_to_column("date") %>% 
  mutate(isSummer = findSummer(had_sst_long$date)) %>% 
  group_by(Year, isSummer) %>% 
  mutate(`Annual Summer Mean` = mean(X1)) %>% 
  filter(month == "09") 


## Combine the annual and summer SST dataframes into single, full one
full_sst_df <- year_sst_mean %>% 
  left_join(had_sst_seasonal[c(8,12)], by = join_by(Year)) %>% 
  mutate(Year = as.numeric(Year)) %>% 
  pivot_longer(`Annual Mean SST`:`Annual Summer Mean`, names_to = "param", values_to = "SST")



# ------------------  
## Example plot
ggplot(data = full_sst_df, aes(x = Year, y = SST, colour = param, group = param)) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_manual(values = c("#00AFBB", "#C4961A")) +
  geom_point() +
  scale_x_continuous(breaks = seq(1870, 2023, 10)) +
  theme_classic() +
  theme(legend.title = element_blank())


#### Save the final dataframe
write.csv(full_sst_df, "Panama_SST_31Jan24.csv", row.names = FALSE)


################################################################################
###### Session information from the latest run on 31 January 2024:
################################################################################

# sessionInfo()
# > sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] splines   stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] transformr_0.1.4      fmsb_0.7.6            performance_0.10.8    modelr_0.1.11         broom_1.0.5          
# [6] sjPlot_2.8.15         ggpubr_0.6.0          gridExtra_2.3         kableExtra_1.3.4.9000 plotly_4.10.3        
# [11] foreach_1.5.2         mgcv_1.9-1            nlme_3.1-164          Hmisc_5.1-1           gganimate_1.0.8      
# [16] fields_15.2           spam_2.10-0           binr_1.1.1            gapminder_1.0.0       gifski_1.12.0-2      
# [21] RColorBrewer_1.1-3    xts_0.13.1            zoo_1.8-12            rnaturalearth_1.0.1   lubridate_1.9.3      
# [26] forcats_1.0.0         stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2           readr_2.1.4          
# [31] tidyr_1.3.0           tibble_3.2.1          tidyverse_2.0.0       openxlsx_4.2.5.2      viridis_0.6.5        
# [36] viridisLite_0.4.2     raster_3.6-26         sp_2.1-2              ggrepel_0.9.4         ggplot2_3.4.4        
# [41] ncdf4_1.22            sf_1.0-15            
# 
# loaded via a namespace (and not attached):
#   [1] rstudioapi_0.15.0  jsonlite_1.8.8     magrittr_2.0.3     estimability_1.4.1 nloptr_2.0.3       farver_2.1.1      
# [7] rmarkdown_2.25     vctrs_0.6.5        minqa_1.2.6        base64enc_0.1-3    terra_1.7-65       rstatix_0.7.2     
# [13] webshot_0.5.5      htmltools_0.5.7    progress_1.2.3     Formula_1.2-5      sjmisc_2.8.9       KernSmooth_2.23-22
# [19] htmlwidgets_1.6.4  emmeans_1.9.0      lifecycle_1.0.4    iterators_1.0.14   pkgconfig_2.0.3    sjlabelled_1.2.0  
# [25] Matrix_1.6-4       R6_2.5.1           fastmap_1.1.1      digest_0.6.33      colorspace_2.1-0   pkgload_1.3.3     
# [31] labeling_0.4.3     fansi_1.0.6        timechange_0.2.0   httr_1.4.7         abind_1.4-5        compiler_4.3.1    
# [37] proxy_0.4-27       withr_2.5.2        htmlTable_2.4.2    backports_1.4.1    carData_3.0-5      DBI_1.2.0         
# [43] maps_3.4.2         ggsignif_0.6.4     MASS_7.3-60        sjstats_0.18.2     classInt_0.4-10    tools_4.3.1       
# [49] units_0.8-5        foreign_0.8-86     zip_2.3.0          nnet_7.3-19        glue_1.6.2         grid_4.3.1        
# [55] checkmate_2.3.1    cluster_2.1.6      lpSolve_5.6.20     generics_0.1.3     gtable_0.3.4       tzdb_0.4.0        
# [61] class_7.3-22       data.table_1.14.10 hms_1.1.3          xml2_1.3.6         car_3.1-2          utf8_1.2.4        
# [67] pillar_1.9.0       tweenr_2.0.2       lattice_0.22-5     tidyselect_1.2.0   knitr_1.45         svglite_2.1.3     
# [73] xfun_0.41          stringi_1.8.3      boot_1.3-28.1      lazyeval_0.2.2     yaml_2.3.8         evaluate_0.23     
# [79] codetools_0.2-19   cli_3.6.2          rpart_4.1.23       xtable_1.8-4       systemfonts_1.0.5  munsell_0.5.0     
# [85] Rcpp_1.0.11        ggeffects_1.3.4    coda_0.19-4        prettyunits_1.2.0  bayestestR_0.13.1  dotCall64_1.1-1   
# [91] lme4_1.1-35.1      mvtnorm_1.2-4      scales_1.3.0       e1071_1.7-14       insight_0.19.7     crayon_1.5.2      
# [97] rlang_1.1.2        rvest_1.0.3 
