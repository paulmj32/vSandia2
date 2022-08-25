### Social / Resilience Factors 
library(tidyverse)
library(tidycensus)
library(sf)
library(terra)

### Get census base map
year = 2019
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
county_map = get_acs(geography = "county", state = state_list,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME, POPULATION) 
county_map_proj = county_map %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 

################################################################################
### Load indicators from Johnson et al., 2020
load("final_data.Rda")
final_data2 = final_data %>% 
  st_set_geometry(NULL) %>%
  dplyr::select(-POPULATION)

county_proj_social = county_map_proj %>%
  left_join(final_data2, by = c("GEOID","NAME"))
#save(county_proj_social, file = "county_proj_social.Rda")
