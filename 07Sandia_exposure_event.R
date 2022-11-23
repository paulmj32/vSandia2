options(java.parameters = "-Xmx10g")
library(bartMachine)
library(tidyverse)
library(tidymodels)
library(tidycensus)
library(sf)
library(xgboost)
library(parallel)
library(doParallel)
library(vip)
library(spdep)
library(pdp)
# library(drat) # these are used to install hurricaneexposure
# addRepo("geanders")
# install.packages("hurricaneexposuredata")
library(hurricaneexposuredata)
library(hurricaneexposure)
library(spatialreg)
library(gstat)
library(ggpubr)
library(DBI)
library(lubridate)


################################################################################
## Get SQL power outage data 
# https://stackoverflow.com/questions/9802680/importing-files-with-extension-sqlite-into-r
# https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html
# https://www.r-bloggers.com/2018/11/quick-guide-r-sqlite/

dbpath = '/Users/paul/vSandia2/Data/Outages/AGM_noaa_CONUS_panel_Proc03Aug2022.sqlite'
mydb = RSQLite::dbConnect(drv = RSQLite::SQLite(), 
                          dbname = dbpath)
tables = dbListTables(mydb) #see tables
myquery = dbGetQuery(mydb, 'SELECT * FROM noaa LIMIT 10')
#myquery = dbGetQuery(mydb, 'SELECT * FROM summary')
dbDisconnect(mydb)

