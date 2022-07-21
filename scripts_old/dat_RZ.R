## Root zone data
# https://www.nrcs.usda.gov/wps/portal/nrcs/detail/soils/survey/geo/?cid=nrcs142p2_053628
# https://gdg.sc.egov.usda.gov/GDGHome_DirectDownLoad.aspx --> soil geographic database 


# NOTE: Too large to read in CONUS dataset (8GB Mac) ... need to run for each state 
# RZ_conus_gdb =  "/Users/paulmj/Downloads/gSSURGO_CONUS/gSSURGO_CONUS.gdb"
# RZ_conus = sf::st_read(dsn = RZ_conus_gdb, layer = "MUPOLYGON")
# 
# RZ_conus_sql = sf::st_read(dsn = RZ_conus_gdb,
#                            query = 'SELECT * FROM "MUPOLYGON" WHERE FID = 1',
#                            layer = "MUPOLYGON")

## Iterate raster processing instead  
# Read in state shapefile from state_list, group by key, join value table, rasterize, save raster 
try_list = state_list[1:48]
for (i in try_list){
  print(i)
  temp_wd = paste("./Data/July2020_gSSURGO_by_State/gSSURGO_", i, sep = "")   # change working directory to file 
  setwd(temp_wd) #set working directory 
  temp_gdb = paste("gSSURGO_", i, ".gdb", sep = "") # define geodatabase 
  temp_sf = sf::st_read(dsn = temp_gdb, layer = "MUPOLYGON") # read in MUPOLYGON layer 
  temp_group = temp_sf %>% group_by(MUKEY) %>% summarise(n = n()) #group by MUKEY 
  temp_val = sf::st_read(dsn = temp_gdb, layer = "Valu1") # read in Value1 table
  temp_join = temp_group %>% left_join(temp_val, by = c("MUKEY" = "mukey")) # join values
  temp_ras = st_rasterize(temp_join["rootznemc"], dx = 100, dy = 100) #100m x 100km resolution raster
  temp_write = paste(i, "_rootznemc100.tif", sep = "") #name of raster
  #temp_ras= st_rasterize(temp_join["rootznemc"], dx = 1000, dy = 1000) #1000m x 1000km resolution raster
  #temp_write = paste(i, "_rootznemc1000.tif", sep = "") #name of raster
  write_stars(temp_ras, dsn = temp_write) # write raster in working directory 
  rm(temp_gdb, temp_sf, temp_group, temp_val, temp_join, temp_ras) #clear working environment (for speed) 
}

## Merge rasters 
# https://stackoverflow.com/questions/50234139/using-mosaic-in-r-for-merge-multiple-geotiff
rast.names = list.files(path = "./Data/rootzone_rasters", full.names = T) # get file names 
rast.list = lapply(rast.names, rast) # read rasters into list 
rast.merge = do.call(merge, rast.list) #use do.call to merge them all at once 
#writeRaster(rast.merge, filename = "./rastmerge.tif")

# project to coordinate system of county shapefile
rz_proj = terra::project(rast.merge, paste("EPSG:", mycrs))

# crop extent to shapefile
rz_crop = crop(rz_proj, vect(county_map_proj))

# mask values to shapefile (i.e., anything outside is NA) 
rz_final = mask(rz_crop, vect(county_map_proj))

# extract values for each county
rz_extract = terra::extract(x = rz_final, y = vect(county_map_proj))
#save(rz_extract, file = "./Data/rz_extract.Rda")
load(file = "./Data/rz_extract.Rda")

# function for mode
mode = function(x) {
  ux = na.omit(unique(x) )
  tab = tabulate(match(x, ux)); ux[tab == max(tab)]
}

# group extracted data 
rz_extract_group = rz_extract %>%
  group_by(ID) %>%
  summarize(RZ_mean = mean(AL_rootznemc100, na.rm = TRUE), #AL_rootznemc100 is just the name of the variable for all states (bc it is the first raster layer)
            RZ_med = median(AL_rootznemc100, na.rm = TRUE), #AL_rootznemc100 is just the name of the variable for all states (bc it is the first raster layer)
            RZ_mode = mode(AL_rootznemc100) #AL_rootznemc100 is just the name of the variable for all states (bc it is the first raster layer)
            )

# join data to counties (bind columns because no GEOID in raster but order is preserved when extracting)
county_map_rz = county_map_proj %>%
  bind_cols(rz_extract_group)
#save(county_map_rz, file = "county_map_rz.Rda")

# gg2 = ggplot(county_map_rz)+
#   geom_sf(aes(fill = RZ_mean), color = NA) +
#   scale_fill_viridis_c(option="plasma", na.value = "grey50") +
#   #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
#   #coord_sf(datum = NA) + #removes gridlines
#   #guides(fill = "none") + #removes legend
#   theme_minimal()  #removes background
# pdf("figure_rz.pdf", width = 7.48, height = 4.5)
# gg2
# dev.off()
