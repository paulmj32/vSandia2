## https://cran.r-project.org/web/packages/hurricaneexposure/vignettes/hurricaneexposure.html

# library(drat)
# addRepo("geanders")
# install.packages("hurricaneexposuredata")

library(hurricaneexposuredata)
library(hurricaneexposure)

data("hurr_tracks")
head(hurr_tracks)

data("closest_dist")
head(closest_dist)

data("county_centers")
data("storm_events")

data("storm_winds")
data("rain")
data("ext_tracks_wind")

asd = map_counties("Katrina-2005", metric = "wind", wind_var = "sust_dur")

head(storm_winds)
head(ext_tracks_wind)
head(rain)
