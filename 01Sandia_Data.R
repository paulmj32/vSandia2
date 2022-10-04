## Load and join data sets
require(tidyverse)
require(tidycensus) 
require(tidymodels)
require(sf)
require(lubridate)
require(corrplot)

# ### Get census base map
# year = 2019
# mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
# options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
# state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
# county_map = get_acs(geography = "county", state = state_list,
#                      variables=c("B01003_001"), year = year, geometry = TRUE, 
#                      cache_table = TRUE)
# county_map = county_map %>%
#   mutate(POPULATION = estimate) %>%
#   dplyr::select(GEOID, NAME, POPULATION) 
# county_map_proj = county_map %>% 
#   st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070
# save(county_map_proj, file = "Data/processed/county_map_proj.Rda")

################################################################################
# LOAD DATA ####################################################################
################################################################################
# Population density (GEOID, POPULATION, DENSITY, and geometry)
load(file = "Data/processed/county_map_proj.Rda")
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) %>% #population per sq-km 
  dplyr::select(-c(AREA, NAME)) 

# NLCD (static)
load(file = "Data/processed/county_proj_nlcd.Rda")
county_map_nlcd = county_map_nlcd %>%
  dplyr::select(-c(NAME, POPULATION, ID, Other, Total)) %>%
  st_set_geometry(NULL)

# DEM (static)
load(file = "Data/processed/county_proj_dem.Rda")
county_map_dem = county_proj_dem %>%
  select(GEOID, DEM_mean, DEM_med, DEM_sd, DEM_min, DEM_max) %>%
  st_set_geometry(NULL)

# Root Zone (static)
load(file = "Data/processed/county_proj_rz.Rda")
county_map_rz = county_map_rz %>%
  select(GEOID, RZ_mean, RZ_med, RZ_mode) %>%
  st_set_geometry(NULL)

# SPI (monthly)
load(file = "Data/processed/county_proj_spi.Rda")
county_map_spi = county_proj_spi %>%
  mutate(date_ym = str_extract(Date, "^.{7}")) %>% # returns yyyy-mm
  dplyr::select(-NAME, -Date, -POPULATION) %>%
  st_set_geometry(NULL)

## Socio-Economic Data
# ## County-level indicators
load(file = "Data/processed/county_proj_social.Rda") 
county_map_social = county_proj_social %>%
  dplyr::select(-NAME, -POPULATION) %>%
  st_set_geometry(NULL)

## Outages
load(file = "Data/processed/county_proj_outages.Rda")
county_map_outages = county_proj_outages %>%
  dplyr::select(-NAME, -POPULATION) %>%
  st_set_geometry(NULL)

################################################################################
# JOIN AND PROCESS DATA ########################################################
################################################################################
## Join data frames 
county_map_JOIN = county_map_area %>%
  inner_join(county_map_outages, by = c("GEOID")) %>% #inner join b/c we don't care about counties w/o events
  left_join(county_map_spi, by = c("GEOID" = "GEOID", "start_ym" = "date_ym")) %>%
  left_join(county_map_nlcd, by = c("GEOID")) %>%
  left_join(county_map_dem, by = c("GEOID")) %>%
  left_join(county_map_rz, by = c("GEOID")) %>%
  left_join(county_map_social, by = c("GEOID"))

## Remove uneccessary and extraneous variables (e.g., WIND2M and WIND10M are basically the same)
# MERRA event categories that we use to derive ours
events_list = c(
  "drought",
  "excessive_heat",
  "extreme_cold_wind_chill",
  "flash_flood",
  "flood",
  "frost_freeze",
  "funnel_cloud",
  "hails",
  "heat",
  "heavy_rain",
  "heavy_snow",
  "high_wind",
  "hurricane",
  "hurricane_typhoon",
  "ice_storm",
  "lake_effect_snow",
  "lakeshore_flood",
  "sleet",
  "strong_wind",
  "tornado",
  "tropical_depression",
  "tropical_storm",
  "wildfire",
  "winter_storm",
  "winter_weather"
)

county_map_ALL = county_map_JOIN %>%
  dplyr::select(-c(outage_number, start_dt, start_ym, end_dt, end_ym)) %>% #these aren't used for ML
  #dplyr::select(-c(mean_cust_out, mean_frac_cust_out)) %>% #don't care about mean values during an outage (only Max)
  dplyr::select(-c(all_of(events_list))) %>% #strip out MERRA event categories 
  dplyr::select(-c(max_WIND2M, max_WIND50M, min_WIND2M, min_WIND50M, 
                   min_WIND10M, #linear combination of max_WIND10M and delta_WIND10M
                   RZ_mean, #highly correlated with RZ_med and truncated at upper values
                   min_humidity_2M, mean_humidity_2M, max_humidity_2M, #just use 10-meter for all of these
                   min_T2M, mean_T2M, max_T2M, delta_T2M, delta_WIND2M, delta_WIND50M)) %>%
  dplyr::select(-c(min_T10M_C, max_T10M_C, mean_T10M_C, min_T2M_C, max_T2M_C, mean_T2M_C, delta_T2M_C, delta_T10M_C)) %>% #already have these in Kelvin 
  dplyr::select(-c(WETLAND)) #already have Wetlands 
#save(county_map_ALL, file = "Data/processed/county_map_ALL.Rda")

## Checks to make sure MERRA definitions add up to our event indicator flags
# asd = county_map_JOIN %>%
#   filter(hurricanes == 1)
# asd2 = county_map_JOIN %>%
#   filter(hurricane > 0 | hurricane_typhoon > 0 | tropical_depression > 0 | tropical_storm > 0)

################################################################################
### CLEAN DATA #################################################################
################################################################################
sf_data = county_map_ALL

# Clean predictors
sf_pred = sf_data %>% 
  dplyr::select(c(duration_hr, DENSITY, droughts:WATEFF)) %>% #select only predictors (and 1 DV)
  st_set_geometry(NULL)
# sf_pred_weather = sf_pred %>% dplyr::select(droughts:delta_WIND10M) %>% dplyr::select(-c(hail, tornadoes))
# asd = cor(sf_pred_weather, use = "complete.obs")
# corrplot(asd, type = "full")

clean_recipe = recipe(duration_hr ~ . , data = sf_pred) %>%
  #step_rm(ln_hrs, contains("total")) %>% #remove total_ice and total_liquid 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  #step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors(), threshold = .875)  #removes highly correlated 
sf_pred_CLEAN = prep(clean_recipe) %>% juice() 
# clean_prep = prep(clean_recipe) #see changes
# sf_pred_CLEAN_weather = sf_pred_CLEAN %>% dplyr::select(droughts:delta_WIND10M)
# asd2 = cor(sf_pred_CLEAN_weather, use = "complete.obs")
# corrplot(asd2, type = "full")

# Data set with cleaned predictors and ALL events
sf_data_ALL = sf_data %>% 
  dplyr::select(c(GEOID, POPULATION, mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out)) %>%
  bind_cols(sf_pred_CLEAN)
#save(sf_data_ALL, file = "Data/processed/sf_data_ALL.Rda")

# Filter data to large events 
sf_data_CLEAN = sf_data_ALL %>%
  dplyr::filter(duration_hr >= 12) #>12 hrs (~95% quantile)
#save(sf_data_CLEAN, file = "Data/processed/sf_data_CLEAN.Rda")
