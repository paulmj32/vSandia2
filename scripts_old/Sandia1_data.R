#### LOAD PACKAGES and SET DIRECTORY 
library(corrplot)
library(viridis)
library(tidyverse)
library(tidycensus) 
library(sf)
library(terra)
library(stars)

setwd("~/Documents/01_VECTOR.nosync/Sandia")

#### CREATE COUNTY MAP
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 

# Get county map
options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
county_map = get_acs(geography = "county", state = state_list,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME, POPULATION) 
# st_crs(county_map)
# st_crs(county_map)$IsGeographic  
# st_crs(county_map)$units_gdal
# st_crs(county_map)$proj4string

# Project county map 
county_map_proj = county_map %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 

# Calculate area and population density of each county 
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) #population per sq-km 


#### LOAD STATIC ENVIRONMENTAL FACTORS
## NLCD
load(file = "./Data/county_map_nlcd.Rda")
county_map_nlcd = county_map_nlcd %>%
  select(GEOID, Barren:Total) %>%
  st_set_geometry(NULL)

## DEM
load(file = "./Data/county_map_dem.Rda")
county_map_dem = county_map_dem %>%
  select(GEOID, DEM_mean:DEM_max) %>%
  st_set_geometry(NULL)

## Root Zone 
load(file = "./Data/county_map_rz.Rda")
county_map_rz = county_map_rz %>%
  select(GEOID, RZ_mean:RZ_mode) %>%
  st_set_geometry(NULL)

# Join static data (Just GEOID and factors) 
county_map_static = county_map_area %>%
  dplyr::select(-NAME, -POPULATION, - AREA, -DENSITY) %>%
  st_set_geometry(NULL) %>% 
  inner_join(county_map_nlcd, by = "GEOID") %>%
  inner_join(county_map_dem, by = "GEOID") %>%
  inner_join(county_map_rz, by = "GEOID") %>%
  select(-Other, -Total)

# See if multicollinearity might be a problem
corX = county_map_static %>% dplyr::select(-GEOID) #only numeric
mycor = cor(corX)
col = viridis(100, direction = -1, option = "C")
corrplot(mycor, method = "circle", tl.col="black", tl.srt=45, tl.cex = 0.7, col = col, cl.cex = 0.7,
         order = "hclust", type = "upper", diag = T, mar = c(1, 1, 1, 1))
# DEM max, DEM mean, DEM med, and DEM min are all highly correlated ... just keep mean
county_map_static_CLEAN = county_map_static %>%
  dplyr::select(-DEM_max, -DEM_min, -DEM_med)


#### LOAD OTHER, DYNAMIC ENVIRONMENTAL FACTORS 
## SPI - monthly
load(file = "./Data/county_map_spi.Rda")
county_map_spi = county_map_spi %>%
  mutate(date_month = str_extract(Date, "^.{7}")) %>% # returns yyyy-mm
  dplyr::select(-NAME, -Date)

# See if multicollinearity might be a problem for SPI values
corX = county_map_spi %>% dplyr::select(-GEOID, -date_month) #only numeric
mycor = cor(corX)
col = viridis(100, direction = -1, option = "C")
corrplot(mycor, method = "circle", tl.col="black", tl.srt=45, tl.cex = 0.7, col = col, cl.cex = 0.7,
         order = "hclust", type = "upper", diag = T, mar = c(1, 1, 1, 1))
# spi_06 has risk of multicollinearity with spi_03 and spi_12 (r = 0.77) ... let's go ahead and drop
county_map_spi_CLEAN = county_map_spi %>%
  dplyr::select(-spi_06)

## Soil Moisture - hourly 
load(file = "./Data/county_map_soil.Rda")
county_map_soil = county_map_soil %>%
  dplyr::select(-NAME, -date_day) #drop NAME and date_day since not using either 
# See if multicollinearity might be a problem for soil moisture values
corX = county_map_soil %>% dplyr::select(-GEOID, -date_hour) #only numeric
mycor = cor(corX, use = "pairwise.complete.obs")
col = viridis(100, direction = -1, option = "C")
corrplot(mycor, method = "circle", tl.col="black", tl.srt=45, tl.cex = 0.7, col = col, cl.cex = 0.7,
         order = "hclust", type = "upper", diag = T, mar = c(1, 1, 1, 1))
#middle layer is highly correlated with top and bottom layer ... remove
county_map_soil_CLEAN = county_map_soil %>%
  dplyr::select(-soil10_40)


#### SOCIO-ECONOMIC FACTORS 
## County-level indicators
load(file = "./Data/final_data.Rda") # we already know multicollinearity is a problem 

# Use factor model
load(file = "./Data/factor_model.Rda") #load factor model
print(ff5, sort=T, cutoff=0.5)

# Select variables not loaded onto one of the five factors 
final_data_cut = final_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(c(GEOID ,
                  AMBULANCES , 
                  BROADBAND ,
                  BUILDING ,
                  BUILDINSP ,
                  BUSINORGS ,
                  BUSISIZE ,
                  CENSRESP ,
                  CHILDCARE ,
                  CIVICS ,
                  CIVILENG ,
                  COLLEGES ,
                  COMMFOOD ,
                  COMMHOUSE ,
                  CROPINS ,
                  DEATHS ,
                  EDEQUITY ,
                  EMPBLDG ,
                  EMPCIVENG ,
                  EMPENVCON ,
                  EMPENVORG ,
                  EMPFIRE ,
                  EMPHIWAY ,
                  EMPHLTH ,
                  EMPINSPECT ,
                  EMPINSUR ,
                  EMPLAND ,
                  EMPLANDSCAPE ,
                  EMPSCIENC ,
                  EMPSPTRN ,
                  EMPSSERVIC ,
                  EMPUNIV ,
                  ENVCONS ,
                  ENVORGS ,
                  FEDEMP ,
                  FIRESTATS ,
                  HAZMIT ,
                  HIGHWYENG ,
                  HOSPBEDS ,
                  HOUSEQUAL ,
                  INTERNET ,
                  LANDARCH ,
                  LANDSUB ,
                  LEGALSERV ,
                  LOCFOOD ,
                  MNTHLTH ,
                  NONPROFITS ,
                  NURSHOMES ,
                  PATTACHIM ,
                  PATTACHRES ,
                  PHYSICIANS ,
                  POPSTAB ,
                  PROFORGS ,
                  PROPINSUR ,
                  PROXCAP ,
                  QAGEWORK ,
                  QED12 ,
                  QEMPL ,
                  QEXTRCT ,
                  QFEMALE ,
                  QFEMLBR ,
                  QGROUPHSE ,
                  QMOHO ,
                  QNATIVE ,
                  QNRRES ,
                  QSERVIND ,
                  RADIO ,
                  RAILMILE ,
                  RECSPORTS ,
                  RELIGAFF ,
                  SCHOOLBUS ,
                  SCHOOLS ,
                  SCISERV ,
                  SPTRANSP ,
                  TELEPHONE ,
                  TELEVISION ,
                  TMPHOUSE ,
                  UTILITYCONS ,
                  WATEFF
                  #WETLAND #already have wetlands from NLCD data
  )) 

# See if multicollinearity might be a problem in remaining variables 
corX = final_data_cut %>% dplyr::select(-GEOID) #only numeric
mycor = cor(corX, use = "pairwise.complete.obs")
col = viridis(100, direction = -1, option = "C")
corrplot(mycor, method = "circle", tl.col="black", tl.srt=45, tl.cex = 0.7, col = col, cl.cex = 0.7,
         order = "hclust", type = "upper", diag = T, mar = c(1, 1, 1, 1))

which(mycor > 0.7 & mycor < 1, arr.ind = T) # FIRESTATS and EMPFIRE highly correlated (0.725) ... drop one

# Join county factor scores and cleaned additional variables 
load(file = "./Data/county_scores.Rda")
county_scores_CLEAN = county_scores %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-POPULATION, -NAME) %>%
  inner_join(final_data_cut, by = c("GEOID")) %>%
  dplyr::select(-FIRESTATS)


#### SAVE CLEANED DATA
save(county_map_area, county_map_static_CLEAN, county_map_spi_CLEAN, county_map_soil_CLEAN, county_scores_CLEAN,
     file = "./Data/SandiaVariables.Rda")


