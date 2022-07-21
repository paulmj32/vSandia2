## SPI (daily) 
# https://www.drought.gov/data-maps-tools/global-gridded-standardized-precipitation-index-spi-cmorph-daily
# historical data doesn't go back that far (and it's based on monthly SPI forecasts estimated at higher temporal resolution)

## SPI (monthly) and SPEI (monthly) 
# https://www.drought.gov/data-maps-tools/us-gridded-standardized-precipitation-index-spi-nclimgrid-monthly
# NWS has most recent observations of SPI (3mth, 6mth, 12mth, 24mth) and forecasts: https://www.cpc.ncep.noaa.gov/products/Drought/Monitoring/spi.shtml

library(tidyverse)
library(tidycensus)
library(sf)
library(raster)
library(terra)
setwd("~/Documents/01_VECTOR.nosync/Sandia")

### USER INPUTS
t_begin = "X2014.01.01" # starting value for data, in "XYYYY.MM.DD" format 
t_end = "X2019.12.01"
file.names = list.files(path = "./Data/SPI_ncdf", full.names = T, pattern = "\\.nc$") 

########################################################################
mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 
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

spi_list = ls()

for (i in 1:length(file.names)){
  file_nm = file.names[i]
  
  #read in ncdf4 file as raster brick
  spi_brick = raster::brick(file_nm)
  #spi_nc = nc_open("./Data/SPI_ncdf/nclimgrid-spi-pearson-03.nc")
  #spi_stars = read_stars("./Data/SPI_ncdf/nclimgrid-spi-pearson-03.nc")
  
  spi_names = names(spi_brick) #1522 names, each corresponding to time "XYYYY.MM.DD"
  raster::crs(spi_brick) #CRS arguments: +proj=longlat +datum=WGS84 +no_defs 
  
  # subset layers for what we want 
  spi_subset = raster::subset(spi_brick, which(spi_names == t_begin):which(spi_names == t_end))
  spi_names_sub = names(spi_subset)
  
  # project to my crs
  spi_proj = raster::projectRaster(spi_subset, crs = crs(county_map_proj)) #project raster into my crs
  
  # convert to SpatRaster (makes for easier projection and extraction) 
  spi_spat = terra::rast(spi_proj)
  
  # project, crop, and mask rasters
  spi_crop = terra::crop(spi_spat, vect(county_map_proj)) #crop to map
  spi_mask = terra::mask(spi_crop, vect(county_map_proj)) #mask to map (make anything outside NA)
  
  # extract SPI layer for each county 
  spi_ext = terra::extract(spi_mask, y = vect(county_map_proj)) 
  spi_ext_g = spi_ext %>% #group by ID 
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # rename columns based on date
  name1 = paste(substr(spi_names_sub, 2,5), "-", substr(spi_names_sub, 7,8), "-", substr(spi_names_sub, 10,11), sep = "")
  colnames(spi_ext_g) = c("ID", name1)
  
  # join to GEOID 
  county_map_spi_w = county_map_proj %>%
    bind_cols(spi_ext_g) %>%
    dplyr::select(-ID, -POPULATION) %>%
    st_set_geometry(NULL)
  
  # variable name 
  s = sub(".*pearson-", "", file_nm) #get spi variable name
  t = sub("*.nc", "", s)
  varname = paste("spi", t, sep = "")
  
  # wide to long format
  county_map_spi_l = county_map_spi_w %>%
    pivot_longer(paste(name1), names_to = "Date", values_to = paste(varname))
 
  # put SPI values into a list
  spi_list[[i]] = data.frame(county_map_spi_l[,4]) #just take SPI value for column in data.frame
}
  
spi_df = do.call(cbind.data.frame, spi_list[1:4])

county_map_spi = county_map_spi_l %>%
  dplyr::select(GEOID:Date) %>%
  bind_cols(spi_df)

save(county_map_spi, file = "./Data/county_map_spi.Rda")

# county_map_join = county_map_proj %>%
#   left_join(county_map_spi, by = c("GEOID","NAME")) %>%
#   filter(Date == name1[1])
# gg2 = ggplot(county_map_join)+
#   geom_sf(aes(fill = spi_24), color = NA) +
#   scale_fill_viridis_c(option="plasma", na.value = "grey50") +
#   #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
#   #coord_sf(datum = NA) + #removes gridlines
#   #guides(fill = "none") + #removes legend
#   theme_minimal()  #removes background
# pdf("figure_spi03.pdf", width = 7.48, height = 4.5)
# gg2
# dev.off()

