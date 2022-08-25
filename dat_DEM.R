## DEM Data
# https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation?qt-science_center_objects=0#qt-science_center_objects
# go to EarthExplorer, Data Sets, Digital Elevation, GTOPO30, Results 
# download which sections you want 
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
# paths
dem_a_path = "./Data/GMTED2010/GMTED2010N10W090_300/10n090w_20101117_gmted_mea300.tif"
dem_b_path = "./Data/GMTED2010/GMTED2010N10W120_300/10n120w_20101117_gmted_mea300.tif"
dem_c_path = "./Data/GMTED2010/GMTED2010N30W090_300/30n090w_20101117_gmted_mea300.tif"
dem_d_path = "./Data/GMTED2010/GMTED2010N30W120_300/30n120w_20101117_gmted_mea300.tif"
dem_e_path = "./Data/GMTED2010/GMTED2010N30W150_300/30n150w_20101117_gmted_mea300.tif"

# load in rasters
dem_a = terra::rast(dem_a_path) 
dem_b = terra::rast(dem_b_path)
dem_c = terra::rast(dem_c_path)
dem_d = terra::rast(dem_d_path)
dem_e = terra::rast(dem_e_path)

# merge rasters
dem_merge1 = terra::merge(dem_a, dem_b) #merge rasters
dem_merge2 = terra::merge(dem_merge1, dem_c)
dem_merge3 = terra::merge(dem_merge2, dem_d)
dem_merge4 = terra::merge(dem_merge3, dem_e)

# project to coordinate system of county shapefile
dem_proj = terra::project(dem_merge4, paste("EPSG:", mycrs))

# crop extent to shapefile
dem_crop = terra::crop(dem_proj, vect(county_map_proj))

# mask values to shapefile (i.e., anything outside is NA) 
dem_final = terra::mask(dem_crop, vect(county_map_proj))

# extract values for each county
dem_extract = terra::extract(x = dem_final, y = vect(county_map_proj))
dem_extract_group = dem_extract %>% #group by ID
  group_by(ID) %>%
  summarize(DEM_mean = mean(`10n090w_20101117_gmted_mea300`, na.rm = TRUE),
            DEM_med = median(`10n090w_20101117_gmted_mea300`, na.rm = TRUE),
            DEM_sd = sd(`10n090w_20101117_gmted_mea300`, na.rm = TRUE),
            DEM_min = min(`10n090w_20101117_gmted_mea300`, na.rm = TRUE),
            DEM_max = max(`10n090w_20101117_gmted_mea300`, na.rm = TRUE)
  )

# join DEM to counties (bind columns because no GEOID in raster but order is preserved when extracting)
county_proj_dem = county_map_proj %>%
  bind_cols(dem_extract_group)
save(county_proj_dem, file = "county_proj_dem.Rda")


gg2 = ggplot(county_proj_dem)+
  geom_sf(aes(fill = DEM_mean), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
  #coord_sf(datum = NA) + #removes gridlines
  #guides(fill = "none") + #removes legend
  theme_minimal()  #removes background
# pdf("figure_dem.pdf", width = 7.48, height = 4.5)
# gg2
# dev.off()
