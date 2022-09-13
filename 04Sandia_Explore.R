#### Look at frequency of events based on socioeconomic status 
require(tidyverse)
require(tidycensus) 
require(tidymodels)
require(sf)
require(lubridate)
require(corrplot)

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
  left_join(county_map_nlcd, by = c("GEOID")) %>%
  left_join(county_map_dem, by = c("GEOID")) %>%
  left_join(county_map_rz, by = c("GEOID")) %>%
  left_join(county_map_social, by = c("GEOID"))

# MERRA event categories
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

our_events = c(
  "droughts",
  "extreme_cold",
  "extreme_heat",
  "floods",
  "hail",
  "high_winds",
  "hurricanes",
  "tornadoes",
  "wildfires",
  "winter_storms"
)

county_map_EXPLORE = county_map_JOIN %>%
  dplyr::select(-c(outage_number, start_dt, start_ym, end_dt, end_ym)) %>% #these aren't used for ML
  #dplyr::select(-c(mean_cust_out, mean_frac_cust_out)) %>% #don't care about mean values during an outage (only Max)
  dplyr::select(-c(all_of(our_events))) %>% #strip out event categories 
  dplyr::select(-c(max_WIND2M, max_WIND50M, min_WIND2M, min_WIND50M, 
                   min_humidity_2M, mean_humidity_2M, max_humidity_2M, #just use 10-meter for all of these
                   min_T2M, mean_T2M, max_T2M, delta_T2M, delta_WIND2M, delta_WIND50M)) %>%
  dplyr::select(-c(min_T10M_C, max_T10M_C, mean_T10M_C, min_T2M_C, max_T2M_C, mean_T2M_C, delta_T2M_C, delta_T10M_C)) %>% #already have these in Kelvin 
  dplyr::select(-c(WETLAND)) #already have Wetlands 

# Filter data to large events 
sf_data = county_map_EXPLORE %>%
  dplyr::filter(duration_hr >= 12) #>12 hrs (~95% quantile)

# Group and Summarize by county
df_data_group = sf_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(GEOID, all_of(events_list)) %>%
  group_by(GEOID) %>%
  summarise(across(all_of(events_list), ~ sum(.x, na.rm = TRUE)))

#Convert Merra groupings to ours 
df_events_group = df_data_group %>%
  mutate(droughts = drought) %>%
  dplyr::select(-c(all_of(events_list)))

#Group predictors by GEOID 
df_data_group_pred = sf_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(GEOID, DENSITY, Barren:WATEFF) %>%
  group_by(GEOID) %>%
  summarise_all(first)

# Clean predictors
clean_recipe = recipe(GEOID ~ . , data = df_data_group_pred) %>%
  #step_rm(ln_hrs, contains("total")) %>% #remove total_ice and total_liquid 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  #step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors(), threshold = .9)  #removes highly correlated 
df_pred_CLEAN = prep(clean_recipe) %>% juice() %>% dplyr::select(-GEOID)
clean_prep = prep(clean_recipe) #see changes

sf_data_EVENTS = sf_data %>%
  group_by(GEOID) %>%
  summarise(POPULATION = first(POPULATION)) %>%
  inner_join(df_events_group, by = c("GEOID")) %>%
  bind_cols(df_pred_CLEAN)
#save(sf_data_EVENTS, file = "Data/processed/sf_data_EVENTS.Rda")


  
  
