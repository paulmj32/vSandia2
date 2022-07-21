## DEM Data
# https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation?qt-science_center_objects=0#qt-science_center_objects
# go to EarthExplorer, Data Sets, Digital Elevation, GTOPO30, Results 
# download which sections you want 
dem_se_path = "./Data/gt30w100n40_dem/gt30w100n40.dem"
dem_ne_path = "./Data/gt30w100n90_dem/gt30w100n90.dem"
dem_sw_path = "./Data/gt30w140n40_dem/gt30w140n40.dem"
dem_nw_path = "./Data/gt30w140n90_dem/gt30w140n90.dem"

# load in rasters
dem_se = rast(dem_se_path)
dem_ne = rast(dem_ne_path)
dem_sw = rast(dem_sw_path)
dem_nw = rast(dem_nw_path)

# merge rasters
dem_merge1 = terra::merge(dem_se, dem_ne)
dem_merge2 = terra::merge(dem_merge1, dem_sw)
dem_merge3 = terra::merge(dem_merge2, dem_nw)

# project to coordinate system of county shapefile
dem_proj = terra::project(dem_merge3, paste("EPSG:", mycrs))

# crop extent to shapefile
dem_crop = crop(dem_proj, vect(county_map_proj))

# mask values to shapefile (i.e., anything outside is NA) 
dem_final = mask(dem_crop, vect(county_map_proj))

# extract values for each county
dem_extract = terra::extract(x = dem_final, y = vect(county_map_proj))
dem_extract_group = dem_extract %>%
  group_by(ID) %>%
  summarize(DEM_mean = mean(gt30w100n40, na.rm = TRUE),
            DEM_med = median(gt30w100n40, na.rm = TRUE),
            DEM_sd = sd(gt30w100n40, na.rm = TRUE),
            DEM_min = min(gt30w100n40, na.rm = TRUE),
            DEM_max = max(gt30w100n40, na.rm = TRUE)
            )

# join DEM to counties (bind columns because no GEOID in raster but order is preserved when extracting)
county_map_dem = county_map_proj %>%
  bind_cols(dem_extract_group)
#save(county_map_dem, file = "county_map_dem.Rda")

# gg2 = ggplot(county_map_dem)+
#   geom_sf(aes(fill = DEM_mean), color = NA) + 
#   scale_fill_viridis_c(option="plasma", na.value = "grey50") +
#   #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
#   #coord_sf(datum = NA) + #removes gridlines 
#   #guides(fill = "none") + #removes legend
#   theme_minimal()  #removes background
# pdf("figure_dem.pdf", width = 7.48, height = 4.5)
# gg2
# dev.off()
