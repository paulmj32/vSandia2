## NLCD land use (30m resolution) 
# https://www.mrlc.gov/data/nlcd-2019-land-cover-conus
require(tidyverse)
require(tidycensus) 
require(terra)
require(sf)

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
# NLCD path 
nlcd_path = "./Data/nlcd_2019_land_cover_l48_20210604/nlcd_2019_land_cover_l48_20210604.img"

# Read in raster with Terra package 
nlcd_ras = terra::rast(nlcd_path)

# Downsample to 100m resolution 
#https://www.patrickbaylis.com/blog/2021-03-13-downsampling-magic/
nlcd_ras_100 = terra::aggregate(nlcd_ras, fact = 100/30, fun = "modal") 
nlcd_levels_100 = levels(nlcd_ras_100)

# Project raster to desired crs
nlcd_proj_100 = terra::project(nlcd_ras_100, paste("EPSG:", mycrs))

# Crop and mask raster 
nlcd_proj_100_crop = terra::crop(nlcd_proj_100, vect(county_map_proj))
nlcd_proj_100_mask = terra::mask(nlcd_proj_100_crop, vect(county_map_proj))

# Extract raster values to each county (values are categorical)
nlcd_proj_100_extract = terra::extract(x = nlcd_proj_100_mask, y = vect(county_map_proj))
#save(nlcd_proj_100_extract, file = "./Data/nlcd_proj_100_extract.Rda")
#load(file = "./Data/nlcd_proj_100_extract.Rda")

# Partition large dataset into percentiles of ID in order for R studio to process
# Assign land class based on NLCD land type
# https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend
# Then group based on ID and LANDCLASS, county number in each 
maxID = max(nlcd_proj_100_extract$ID) #3107
try01 = nlcd_proj_100_extract %>% filter(ID < 301) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try02 = nlcd_proj_100_extract %>% filter(ID < 601 & ID >= 301) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try03 = nlcd_proj_100_extract %>% filter(ID < 901 & ID >= 601) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try04 = nlcd_proj_100_extract %>% filter(ID < 1201 & ID >= 901) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try05 = nlcd_proj_100_extract %>% filter(ID < 1501 & ID >= 1201) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try06 = nlcd_proj_100_extract %>% filter(ID < 1801 & ID >= 1501) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try07 = nlcd_proj_100_extract %>% filter(ID < 2101 & ID >= 1801) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try08 = nlcd_proj_100_extract %>% filter(ID < 2401 & ID >= 2101) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try09 = nlcd_proj_100_extract %>% filter(ID < 2701 & ID >= 2401) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

try10 = nlcd_proj_100_extract %>% filter(ID <= (maxID) & ID >= 2701) %>%
  mutate(LANDCLASS = case_when(
    `NLCD Land Cover Class` %in% c("Open Water", "Perennial Snow/Ice") ~ "Water",
    `NLCD Land Cover Class` %in% c("Developed, Open Space", "Developed, Low Intensity", "Developed, Medium Intensity", "Developed, High Intensity") ~ "Developed",
    `NLCD Land Cover Class` %in% c("Barren Land") ~ "Barren",
    `NLCD Land Cover Class` %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest") ~ "Forest",
    `NLCD Land Cover Class` %in% c("Shrub/Scrub") ~ "Shrub",
    `NLCD Land Cover Class` %in% c("Herbaceous")~ "Herbaceous",
    `NLCD Land Cover Class` %in% c("Hay/Pasture", "Cultivated Crops") ~ "Cultivated",
    `NLCD Land Cover Class` %in% c("Woody Wetlands", "Emergent Herbaceous Wetlands") ~ "Wetlands",
    TRUE ~ "Other"
  )) %>%
  group_by(ID, LANDCLASS) %>%
  count()

# Stack data sets together 
nlcd_join_stack = try01 %>%
  bind_rows(try02) %>%
  bind_rows(try03) %>%
  bind_rows(try04) %>%
  bind_rows(try05) %>%
  bind_rows(try06) %>%
  bind_rows(try07) %>%
  bind_rows(try08) %>%
  bind_rows(try09) %>%
  bind_rows(try10) 
#save(nlcd_join_stack, file = "nlcd_join_stack.Rda")
#load(file = "nlcd_join_stack.Rda")

# Expand LANDCLASS into wide format 
nlcd_wide = nlcd_join_stack %>%
  pivot_wider(names_from = LANDCLASS, values_from = n)

# Column totals of land type 
nlcd_wide2 = nlcd_wide %>%
  replace(is.na(.), 0) %>%
  mutate(Total = rowSums(across(Barren:Other)))

# Percentages of land type
nlcd_perc = nlcd_wide2 %>%
  mutate(across(Barren:Total, ~ .x / Total))

# join to census shapefile
county_proj_nlcd = county_map_proj %>%
  bind_cols(nlcd_perc)
#save(county_proj_nlcd, file = "county_proj_nlcd.Rda")

# gg2 = ggplot(county_map_nlcd)+
#   geom_sf(aes(fill = Developed), color = NA) + 
#   scale_fill_viridis_c(option="plasma", na.value = "grey50") +
#   #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
#   #coord_sf(datum = NA) + #removes gridlines 
#   #guides(fill = "none") + #removes legend
#   #theme_minimal()  #removes background
#   theme_dark() +
#   labs(title = "NLCD Land Use", fill = "Developed\n(prop.)") + 
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#   )
# gg2
