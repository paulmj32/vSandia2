### Power Outage Data 
library(tidyverse)
library(tidycensus)
library(sf)
library(terra)
library(DBI)
library(lubridate)

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
## Get SQL power outage data 
# https://stackoverflow.com/questions/9802680/importing-files-with-extension-sqlite-into-r
# https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html
# https://www.r-bloggers.com/2018/11/quick-guide-r-sqlite/
dbpath = "./Data/Outages/AGM_summary_CONUS_panel_Proc03Aug2022.sqlite"
mydb = RSQLite::dbConnect(drv = RSQLite::SQLite(), 
                          dbname = dbpath)
tables = dbListTables(mydb) #see tables
#myquery = dbGetQuery(mydb, 'SELECT * FROM summary LIMIT 10')
myquery = dbGetQuery(mydb, 'SELECT * FROM summary')
dbDisconnect(mydb)

# Format year and month for easier joins later
county_outages = myquery %>%
  mutate(start_dt = as.character(as_datetime(start_date))) %>%
  mutate(end_dt = as.character(as_datetime(end_date))) %>%
  mutate(start_ym = str_extract(start_dt, "^.{7}")) %>% 
  mutate(end_ym = str_extract(end_dt, "^.{7}")) %>% 
  dplyr::select(-c(start_date, end_date)) %>%
  relocate(c(start_dt, start_ym, end_dt, end_ym), .after = outage_number)

## Join tables
county_proj_outages = county_map_proj %>%
  inner_join(county_outages, by = c("GEOID" = "fips_code")) 
#save(county_proj_outages, file = "county_proj_outages.Rda")
