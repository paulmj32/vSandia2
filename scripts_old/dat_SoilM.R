## SOIL MOISTURE (hourly) 
# https://disc.gsfc.nasa.gov/datasets?keywords=NLDAS&page=1&subject=Soils
# We're using NOAH because it matches the soil layers we're interested in
# available in either grb or ncdf format ... grb is more current (hourly data with a 3-day lag ... perfect for forecasting)
# subset/Get data --> download file link list
# Follow instructions for downloading via wget or curl (I like latter): https://disc.gsfc.nasa.gov/data-access#mac_linux_wget
# example terminal code for MAC: cat '/Users/paulmj/Downloads/subset_NLDAS_NOAH0125_H_002_20211116_191020.txt' | tr -d '\r' | xargs -n 1 curl -LJO -n -c ~/.urs_cookies -b ~/.urs_cookies)
# make sure to change directory to whatever folder you want the files downloaded to before running

# soil_db = "/Users/paulmj/Downloads/NLDAS_NOAH_2018_hourly_grb/NLDAS_NOAH0125_H.A20180101.0000.002.grb"
# asd = read_stars(soil_db, along = "band")
# 
# # GRIB INFO (better for identifying names of layers but not reading in the file itself) 
# library(rNOMADS)
# GribInfo(soil_db, file.type = "grib1")
# soil0_10 = asd[ , , , 26] # units kg/m^2 
# soil10_40 = asd[ , , , 27]
# soil40_100 = asd[ , , , 28]

## Get list of .grb files 
grb.names = list.files(path = "/Users/paulmj/Downloads/NLDAS_NOAH_2018_hourly_grb", full.names = T,
                       pattern = "\\.grb$") # get file names and restrict pattern to ending in .grb 

## Iterate through .grb files (takes less memory) length(grb.names)
for (i in 7297:length(grb.names)){
  temp_i = grb.names[i]
  print(temp_i)
  
  # Read in data
  #temp_stars = read_stars(temp_i) #read in file as stars class
  temp_rast = rast(temp_i) #read in as raster layers
  
  # Project, crop, and mask rasters
  temp_proj = terra::project(temp_rast, paste("EPSG:", mycrs)) #transform raster into my CRS
  temp_crop = terra::crop(temp_proj, vect(county_map_proj)) #crop to map
  temp_mask = terra::mask(temp_crop, vect(county_map_proj)) #mask to map (make anything outside NA)
  
  # Get soil moisture layers (3 of them)
  temp_s26 = temp_mask[[26]] # soil moisture 0-10cm (kg/m^2)
  temp_s27 = temp_mask[[27]] # soil moisture 10-40cm (kg/m^2)
  temp_s28 = temp_mask[[28]] # soil moisture 40-100cm (kg/m^2)
  
  # Extract soil layer values by county 
  temp_extract26 = terra::extract(x = temp_s26, y = vect(county_map_proj)) %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  temp_extract27 = terra::extract(x = temp_s27, y = vect(county_map_proj)) %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  temp_extract28 = terra::extract(x = temp_s28, y = vect(county_map_proj)) %>%
    group_by(ID) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Join temp extractions by ID
  temp_extract = temp_extract26 %>%
    inner_join(temp_extract27, by = c("ID")) %>%
    inner_join(temp_extract28, by = c("ID"))
  colnames(temp_extract) = c("ID", "soil0_10", "soil10_40", "soil40_100")
  
  # Add date_hour (character) variable whose value is the file name, formatted to match outage data ("2018-01-01T00:00:00Z")
  temp_name1 = sub(".002.grb", "", sub(".*NLDAS_NOAH0125_H.A", "", temp_i)) #"YYYYMMDD.HHMM"
  temp_name2 = paste(substr(temp_name1, 1, 4), "-", substr(temp_name1, 5, 6), "-", substr(temp_name1, 7, 8), "T", substr(temp_name1, 10, 11), ":", substr(temp_name1, 12, 13), ":00Z", sep = "")
  temp_day = substr(temp_name2, 1,10)
  
  temp_out = temp_extract
  temp_out$date_hour = temp_name2
  temp_out$date_day = temp_day
  
  # Data-frame of GEOID, date-time (hourly and daily), and 3 soil moisture layer variables 
  temp_map_soil = county_map_proj %>%
    bind_cols(temp_out) %>%
    dplyr::select(-ID, -POPULATION) %>%
    st_set_geometry(NULL)
  
  # Use first dataset of each month as starting output and row bind subsequent days to it
  # if (i == 1) {
  #   county_map_soil = temp_map_soil
  # }
  if (i %in% c(1, 745, 1417, 2161, 2881, 3625, 4345, 5089, 5833, 6553, 7297, 8017)) {
    county_map_soil = temp_map_soil
  }
  else {
    county_map_soil = county_map_soil %>%
      bind_rows(temp_map_soil)
  }
  
  # Save county_map_soil data at end of every month 
  if (i == 744) {
    county_map_soil_jan = county_map_soil
    save(county_map_soil_jan, file = "./Data/county_map_soil_jan.Rda")
    rm(county_map_soil, county_map_soil_jan) #clear working environment for speed
  }
  if (i == 1416) {
    county_map_soil_feb = county_map_soil
    save(county_map_soil_feb, file = "./Data/county_map_soil_feb.Rda")
    rm(county_map_soil, county_map_soil_feb) #clear working environment for speed
  }
  if (i == 2160) {
    county_map_soil_mar = county_map_soil
    save(county_map_soil_mar, file = "./Data/county_map_soil_mar.Rda")
    rm(county_map_soil, county_map_soil_mar) #clear working environment for speed
  }
  if (i == 2880) {
    county_map_soil_apr = county_map_soil
    save(county_map_soil_apr, file = "./Data/county_map_soil_apr.Rda")
    rm(county_map_soil, county_map_soil_apr) #clear working environment for speed
  }
  if (i == 3624) {
    county_map_soil_may = county_map_soil
    save(county_map_soil_may, file = "./Data/county_map_soil_may.Rda")
    rm(county_map_soil, county_map_soil_may) #clear working environment for speed
  }
  if (i == 4344) {
    county_map_soil_jun = county_map_soil
    save(county_map_soil_jun, file = "./Data/county_map_soil_jun.Rda")
    rm(county_map_soil, county_map_soil_jun) #clear working environment for speed
  }
  if (i == 5088) {
    county_map_soil_jul = county_map_soil
    save(county_map_soil_jul, file = "./Data/county_map_soil_jul.Rda")
    rm(county_map_soil, county_map_soil_jul) #clear working environment for speed
  }
  if (i == 5832) {
    county_map_soil_aug = county_map_soil
    save(county_map_soil_aug, file = "./Data/county_map_soil_aug.Rda")
    rm(county_map_soil, county_map_soil_aug) #clear working environment for speed
  }
  if (i == 6552) {
    county_map_soil_sep = county_map_soil
    save(county_map_soil_sep, file = "./Data/county_map_soil_sep.Rda")
    rm(county_map_soil, county_map_soil_sep) #clear working environment for speed
  }
  if (i == 7296) {
    county_map_soil_oct = county_map_soil
    save(county_map_soil_oct, file = "./Data/county_map_soil_oct.Rda")
    rm(county_map_soil, county_map_soil_oct) #clear working environment for speed
  }
  if (i == 8016) {
    county_map_soil_nov = county_map_soil
    save(county_map_soil_nov, file = "./Data/county_map_soil_nov.Rda")
    rm(county_map_soil, county_map_soil_nov) #clear working environment for speed
  }
  if (i == 8760) {
    county_map_soil_dec = county_map_soil
    save(county_map_soil_dec, file = "./Data/county_map_soil_dec.Rda")
    rm(county_map_soil, county_map_soil_dec) #clear working environment for speed
  }
}

county_map_soil = county_map_soil_jan %>%
  bind_rows(county_map_soil_feb) %>%
  bind_rows(county_map_soil_mar) %>%
  bind_rows(county_map_soil_apr) %>%
  bind_rows(county_map_soil_may) %>%
  bind_rows(county_map_soil_jun) %>%
  bind_rows(county_map_soil_jul) %>%
  bind_rows(county_map_soil_aug) %>%
  bind_rows(county_map_soil_sep) %>%
  bind_rows(county_map_soil_oct) %>%
  bind_rows(county_map_soil_nov) %>%
  bind_rows(county_map_soil_dec) 
save(county_map_soil, file = "./Data/county_map_soil.Rda")

# asd = county_map_soil_jan %>%
#   filter(date_hour == temp_name2)
# gg2 = ggplot(asd)+
#   geom_sf(aes(fill = soil0_10), color = NA) +
#   scale_fill_viridis_c(option="plasma", na.value = "grey50") +
#   #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
#   #coord_sf(datum = NA) + #removes gridlines
#   #guides(fill = "none") + #removes legend
#   theme_minimal()  #removes background
# pdf("figure_rz.pdf", width = 7.48, height = 4.5)
# gg2
# # dev.off()
