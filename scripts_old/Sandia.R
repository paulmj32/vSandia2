options(java.parameters = "-Xmx6g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(lme4)

library(tidyverse)
library(tidycensus) 
library(sf)
library(terra)
library(stars)
library(sqldf)

setwd("~/Documents/01_VECTOR.nosync/Sandia")
mycrs = 5070

##### CONTIGUOUS US COUNTY MAP ################################################################
year=2019
options(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
county_map = get_acs(geography = "county", state = state_list,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME, POPULATION) 

# # county map info 
# st_crs(county_map)
# st_crs(county_map)$IsGeographic  
# st_crs(county_map)$units_gdal
# st_crs(county_map)$proj4string

## Project county map and calculate area and population density of each county 
county_map_proj = county_map %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 

county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% ## calculate area of each county (sq-meters); as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) #population per sq-km 

##### STATIC VARIABLES ######################################################################
## NLCD
load(file = "./Data/county_map_nlcd.Rda")
county_map_nlcd = county_map_nlcd %>%
  select(GEOID, Barren:Total) %>%
  st_set_geometry(NULL)

## DEM
load(file = "./Data/county_map_dem.Rda")
county_map_dem = county_map_dem %>%
  select(GEOID, DEM_mean:DEM_max) %>%
  st_set_geometry(NULL)

## Root Zone 
load(file = "./Data/county_map_rz.Rda")
county_map_rz = county_map_rz %>%
  select(GEOID, RZ_mean:RZ_mode) %>%
  st_set_geometry(NULL)

## Join data 
county_map_static = county_map_area %>%
  dplyr::select(-NAME, -POPULATION, - AREA, -DENSITY) %>%
  st_set_geometry(NULL) %>% 
  inner_join(county_map_nlcd, by = "GEOID") %>%
  inner_join(county_map_dem, by = "GEOID") %>%
  inner_join(county_map_rz, by = "GEOID") %>%
  select(-Other, -Total)
#save(county_map_static, file = "./Data/county_map_static.Rda")

##### DYNAMIC VARIABLES ######################################################################
## SPI - monthly
load(file = "./Data/county_map_spi.Rda")
county_map_spi = county_map_spi %>%
  mutate(date_month = str_extract(Date, "^.{7}")) %>% # returns yyyy-mm
  dplyr::select(-NAME, -Date)

## Soil Moisture - hourly 
load(file = "./Data/county_map_soil.Rda")
county_map_soil = county_map_soil %>%
  dplyr::select(-NAME, -date_day) #drop NAME and date_day since not using either 

##### SOCIO-ECONOMIC VARIABLES ###############################################################
# county factor scores
load(file = "./Data/county_scores.Rda")
county_scores = county_scores %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-POPULATION, -NAME)

# factor model 
load(file = "./Data/factor_model.Rda")
print(ff5, sort=T, cutoff=0.5)

# county raw data for indicators
load(file = "./Data/final_data.Rda")
final_data2 = final_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(c(GEOID ,
                AMBULANCES , # select all variables not loaded onto one of the five factors
                BROADBAND ,
                BUILDING ,
                BUILDINSP ,
                BUSINORGS ,
                BUSISIZE ,
                CENSRESP ,
                CHILDCARE ,
                CIVICS ,
                CIVILENG ,
                COLLEGES ,
                COMMFOOD ,
                COMMHOUSE ,
                CROPINS ,
                DEATHS ,
                EDEQUITY ,
                EMPBLDG ,
                EMPCIVENG ,
                EMPENVCON ,
                EMPENVORG ,
                EMPFIRE ,
                EMPHIWAY ,
                EMPHLTH ,
                EMPINSPECT ,
                EMPINSUR ,
                EMPLAND ,
                EMPLANDSCAPE ,
                EMPSCIENC ,
                EMPSPTRN ,
                EMPSSERVIC ,
                EMPUNIV ,
                ENVCONS ,
                ENVORGS ,
                FEDEMP ,
                FIRESTATS ,
                HAZMIT ,
                HIGHWYENG ,
                HOSPBEDS ,
                HOUSEQUAL ,
                INTERNET ,
                LANDARCH ,
                LANDSUB ,
                LEGALSERV ,
                LOCFOOD ,
                MNTHLTH ,
                NONPROFITS ,
                NURSHOMES ,
                PATTACHIM ,
                PATTACHRES ,
                PHYSICIANS ,
                POPSTAB ,
                PROFORGS ,
                PROPINSUR ,
                PROXCAP ,
                QAGEWORK ,
                QED12 ,
                QEMPL ,
                QEXTRCT ,
                QFEMALE ,
                QFEMLBR ,
                QGROUPHSE ,
                QMOHO ,
                QNATIVE ,
                QNRRES ,
                QSERVIND ,
                RADIO ,
                RAILMILE ,
                RECSPORTS ,
                RELIGAFF ,
                SCHOOLBUS ,
                SCHOOLS ,
                SCISERV ,
                SPTRANSP ,
                TELEPHONE ,
                TELEVISION ,
                TMPHOUSE ,
                UTILITYCONS ,
                WATEFF
                #WETLAND #already have wetlands
                )) 

all_soc_dat = county_scores %>%
  inner_join(final_data2, by = c("GEOID"))

# factor model itself - extend beyond 5 factors 
finalNA = which(is.na(final_data), arr.ind = T) #check for NAs
slice = c(finalNA[,1]) #county row index to remove from analysis 
slice = unique(slice); length(slice)  #remove 75 counties (~2%) 
scale2 = function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm) #z-score transform
scale3 = function(x, na.rm = FALSE) (x - min(x, na.rm = na.rm)) / (max(x, na.rm = na.rm) - min(x, na.rm = na.rm)) #min-max transform
yy = final_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-c(GEOID, NAME, POPULATION)) %>% 
  slice(-slice) %>% #remove NAs
  #mutate_all(~scale(.)) #this also works 
  mutate_all(~scale2(., na.rm = TRUE))
#mutate_all(~scale3(., na.rm = TRUE))
#mutate_all(list(~scale)) #standardize (mean 0, sd 1) #outdated code
summary(yy)

library(psych)
fa.parallel(yy, fm="ml") #parallel analysis and skree plots
#retain 10 factors based on Kaiser criterion
ff10 = factanal(yy, factors = 10, rotation = "varimax", scores = "regression") #fits using Maximum Likelihood
print(ff10, sort=T, cutoff=0.5)

#PCA - just remove all correlation from dataset 
pc129 = principal(yy, nfactors = 129, rotate = "varimax", scores = TRUE)
print(pc129, sort=T, cutoff=0.5)

county_pca = final_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(GEOID) %>%
  slice(-slice) %>%
  bind_cols(data.frame(pc129$scores))
county_pca_scores = county_map %>%
  st_set_geometry(NULL) %>%
  left_join(county_pca, by = "GEOID") %>% #put back into non-sliced indices
  dplyr::select(-NAME, -POPULATION)


##### OUTAGE DATA ############################################################################
## Read in outage data
#outages_csv = read.csv("./Data/SE_states_outage_merra_2018.csv", header = T)
#save(outages_csv, file = "./Data/outages.Rda")
load(file = "Data/outages.Rda")

## Format FIPS to character and include starting zero
outages_csv$fips_code = as.character(outages_csv$fips_code)
fips = outages_csv$fips_code
fips0 = str_pad(fips, 5, pad = "0")
outages_csv$fips_code = fips0

county_map_outages = outages_csv %>%
  mutate(rowID = row_number()) %>% # add row number to index events 
  mutate(date_day = str_extract(date_hour, "^.{10}")) %>% #returns first 10 characters of date_time (yyyy-mm-dd)
  mutate(date_month = str_extract(date_hour, "^.{7}")) # returns yyyy-mm

# Filter data 
county_map_outages_only = county_map_outages %>%
  filter(outage_status %in% c("pre", "start", "during", "end")) #filter events flagged as outages

## Calculate rolling 24-hr lagged min, max, avg, and sd values for all climate variables 

# # Add column of lagged outage variable as predictor (to capture temporal dependence between outputs)
# outages_lag1_index = county_map_outages_only$rowID - 1 
# outages_only$lag1 = outages_csv$hr_mean_frac_out[outages_lag1_index] # get t-1 customer outages
# 
# # Cleaned data-set for joining/statistical learning 
# county_outages = outages_only %>%
#   dplyr::relocate(date_day, date_month, .after = date_hour) %>%
#   dplyr::select(-c("state_fips", "state", "county", "year", "county_population", "hr_mean_customers_out", "outage_status", "outage_number", "rowID"))


##### JOIN ALL DATA #########################################################################
county_map_join = county_map_area %>%
  inner_join(county_map_outages_only, by = c("GEOID" = "fips_code")) %>% # join outage/hourly weather data 
  left_join(county_map_static, by = c("GEOID")) %>% #join static environmental variables
  #left_join(county_scores, by = c("GEOID")) %>% #join socio-economic variables (Factor Model scores) 
  #left_join(all_soc_dat, by = c("GEOID")) %>% #5 factors plus remaining variables
  left_join(county_pca_scores, by = c("GEOID")) %>% #join PCA scores (uncorrelated feature space)
  inner_join(county_map_soil, by = c("GEOID", "date_hour")) %>% #join soil moisture by county and hour time-stamp 
  inner_join(county_map_spi, by = c("GEOID", "date_month")) #join SPI by GEOID and month 
#save(county_map_join, file = "./Data/county_map_join.Rda")


# Group by event 
county_map_group = county_map_join %>%
  st_set_geometry(NULL) %>%
  group_by(outage_number, GEOID) %>%
  summarise(hours = n(), customers_out = max(hr_mean_customers_out), customers_out_perc = sum(hr_mean_customers_out) / sum(POPULATION),
            Density = mean(DENSITY), 
            PS_mean = mean(PS), PS_sd = sd(PS), #mean and sd of surface pressure
            SLP_mean = mean(SLP), PS_sd = sd(SLP), # same for sea level pressure
            QV_max = max(QV10M), # max specific humidity
            U_max = max(U10M),  
            V_max = max(V10M),
            WIND_max = max(WIND10M), 
            T_mean = mean(T10M), T_sd = sd(T10M), 
            TQI_mean = mean(TQI), TQL_mean = mean(TQL), TQV_mean = mean(TQV),
            soil10_mean = mean(soil0_10), soil40_mean = mean(soil10_40), soil100_mean = mean(soil40_100),
            spi03_mean = mean(spi_03), spi06_mean = mean(spi_06), spi12_mean = mean(spi_12),spi24_mean = mean(spi_24)
            )

county_map_group_join = county_map_group %>%
  left_join(county_map_static, by = c("GEOID")) %>% #join static environmental variables
  #left_join(county_scores, by = c("GEOID")) %>% #join socio-economic variables (Factor Model scores) 
  #left_join(all_soc_dat, by = c("GEOID")) %>% #5 factors plus remaining variables
  left_join(county_pca_scores, by = c("GEOID")) #join PCA scores (uncorrelated feature space)

#save(county_map_group_join, file = "./Data/county_map_group_join.Rda")
#load(file = "./Data/county_map_group_join.Rda")

##### MACHINE LEARNING ######################################################################
## BART model (with lag predictor) - BartMachine
rm(list=setdiff(ls(), "county_map_group_join")) # removes all variables from environment except 
gc() #garbage collect / free up memory 

# Modeling 3 different independent variables: ln(hours), Kilo customer_out_hours, and customer_out_perc (per event x county) 
df_bart = data.frame(county_map_group_join) %>%
  drop_na %>% 
  dplyr::filter(customers_out_perc < 1) %>% #remove bad data 
  mutate(customer_out_hrs = customers_out / 1000 * hours) #Kilo customer_out_hrs
y1 = log(df_bart$hours)
y2 = log(df_bart$customer_out_hrs)
y3 = log(df_bart$customers_out_perc)

X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -customers_out, -customers_out_perc, -customer_out_hrs, -hours)

set.seed(32) 
bartY1 = bartMachine(X, y1) 
bartY2 = bartMachine(X, y2) 
bartY3 = bartMachine(X, y3) 

#### PREDICTING DURATION - with customer outages ##########################
## Using actual outages 
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -customers_out_perc, -customer_out_hrs, -hours)
set.seed(32) 
bartY1 = bartMachine(X, y1) 

## restrict to only big outages (80 or 90 percentile)
df_bart = data.frame(county_map_group_join) %>%
  drop_na %>% 
  dplyr::filter(customers_out_perc < 1) %>% #remove bad data 
  dplyr::filter(hours > quantile(county_map_group_join$hours, .9)) %>%  
  mutate(customer_out_hrs = customers_out / 1000 * hours) #Kilo customer_out_hrs
y1 = log(df_bart$hours)
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -customers_out_perc, -customer_out_hrs, -hours)
set.seed(32) 
bartY1 = bartMachine(X, y1) 

## small events (<= 12 hours)
df_bart = data.frame(county_map_group_join) %>%
  drop_na %>% 
  dplyr::filter(customers_out_perc < 1) %>% #remove bad data 
  dplyr::filter(hours <= quantile(county_map_group_join$hours, .9)) %>%  
  mutate(customer_out_hrs = customers_out / 1000 * hours) #Kilo customer_out_hrs
y1 = log(df_bart$hours)
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -customers_out_perc, -customer_out_hrs, -hours)
set.seed(32) 
bartY1 = bartMachine(X, y1) 




bart = bartY1
Y = y1

# df_bart = county_map_join %>%
#   st_set_geometry(NULL) %>%
#   dplyr::select(-c(GEOID, NAME, POPULATION, AREA, date_hour, date_day, date_month)) %>%
#   dplyr::slice(1:10000)
# y = df_bart$hr_mean_frac_out
# X = df_bart; X$hr_mean_frac_out = NULL
# 
# set.seed(32) 
# bart = bartMachine(X, y, replace_missing_data_with_x_j_bar = TRUE) 

predictions = predict(bart, X)
CI = round(calc_credible_intervals(bart, X, ci_conf = 0.95), 2)
gg = dplyr::tibble(x_mean = predictions,
                   lower = CI[,1],
                   upper = CI[,2],
                   actual = Y
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))

plot_filtering_estimates <- function(df) {
  p <- ggplot(data = gg, aes(x = index)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "red") + #quantiles 
    geom_line(aes(y = x_mean), colour = "red", size = 0.4, alpha = 0.75) + #mean estimate
    geom_point(aes(y = actual), colour = "gray32", #actual observation points
               size = 0.9, shape = 16, alpha = 0.9) +
    ylab("log(outages_hrs)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Event (outage x county)") +
    ggtitle("Bayesian Additive Regression Tree") +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
plot_filtering_estimates(gg)

check_bart_error_assumptions(bart)
plot_convergence_diagnostics(bart)

varimp = investigate_var_importance(bart, num_replicates_for_avg = 5)

pd_plot(bartY1, j = "Factor2")

#bart_cv = bartMachineCV(X, y, replace_missing_data_with_x_j_bar = TRUE) 

## DBART with GEOID as  grouping variable

#### PLOTS ##########################################################################
gg = ggplot(county_map_static)+
  geom_sf(aes(fill = DEM_mean), color = NA) + 
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  #geom_sf(fill = NA, show.legend = F, color = "black", lwd = 0.005)+
  #coord_sf(datum = NA) + #removes gridlines 
  #guides(fill = "none") + #removes legend
  #theme_minimal() +  #removes background
  theme_dark() +
  #labs(title = "NLCD Land Use", fill = "Forest\n(prop.)") + 
  #labs(title = "Root-zone Depth", fill = "Mean (cm)") + 
  labs(title = "Elevation", fill = "Mean (m)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

pdf("US_dem_mean.pdf", width = 7.48, height = 4.5)
gg
dev.off()

## Exploratory analysis

#autocorrelation of one county through time (make sure to find example where outagedata) 

#boxplot of hr_mean_frac_out to see where outliers may be
summary(outages_csv$hr_mean_frac_out) #highly skewed right 
boxplot(outages_csv$hr_mean_frac_out)

#filter only outage events
outages_only = outages_csv %>%
  filter(hr_mean_frac_out >= 0.01) #could potentially use status instead? 
summary(outages_only$hr_mean_frac_out)
boxplot(outages_only$hr_mean_frac_out)

wind_10m = outages_csv %>%
  filter(WIND10M >= 4.69) #90th percentile

#order data with respect to mean frac_out
outages_only1 = outages_only %>% arrange(hr_mean_frac_out)
plot(outages_only1$hr_mean_frac_out)

outages_bad = outages_csv %>%
  filter(hr_mean_frac_out > .8)


## Try some statistical models
# BART
df_bart = outages_csv %>% 
  filter(hr_mean_frac_out >= 0.3) %>% #filter out non-outages
  dplyr::select(-c(state_fips, state, county, fips_code, year, date_hour, county_population, hr_mean_customers_out, outage_status, outage_number))
y = df_bart$hr_mean_frac_out
X = df_bart; X$hr_mean_frac_out = NULL

set.seed(32) #for some reason BART doesn't like set.seed() / the command doesn't afffect stochastic elements (just repeat)
bart_cv = bartMachineCV(X, y) #bartMachine CV win: k: 2 nu, q: 5, 0.99 m: 200 
print(bart_cv)

#predictions
predictions = predict(bart_cv, X)
rmse = sqrt(mean((predictions - df_bart$hr_mean_frac_out)^2)) 
rsq = 1 - sum((df_bart$hr_mean_frac_out - predictions)^2) / sum((df_bart$hr_mean_frac_out - mean(df_bart$hr_mean_frac_out))^2)

#check BART assumptions
check_bart_error_assumptions(bart_cv)
plot_convergence_diagnostics(bart_cv)

CI = round(calc_credible_intervals(bart_cv, X, ci_conf = 0.95), 2)
gg = dplyr::tibble(x_mean = predictions,
                   lower = CI[,1],
                   upper = CI[,2],
                   actual = df_bart$hr_mean_frac_out
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))

plot_filtering_estimates <- function(df) {
  p <- ggplot(data = gg, aes(x = index)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "red") + #quantiles 
    geom_point(aes(y = actual), colour = "gray32", #actual observation points
               size = 0.9, shape = 16, alpha = 0.9) +
    geom_line(aes(y = x_mean), colour = "red", size = 0.4) + #mean estimate
    ylab("hr_mean_frac_out") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Index") +
    ggtitle("Bayesian Additive Regression Tree") +
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
plot_filtering_estimates(gg)

summary = investigate_var_importance(bart_cv, num_replicates_for_avg = 5)
pd_plot(bart_cv, j = "PS")
pd_plot(bart_cv, j = "WIND2M")
pd_plot(bart_cv, j = "SLP")
pd_plot(bart_cv, j = "V2M")
pd_plot(bart_cv, j = "TQV")

library(corrplot)
library(viridis)
mycor = cor(X)
col = viridis(100, direction = -1, option = "C")
corrplot(mycor, method = "circle", tl.col="black", tl.srt=45, tl.cex = 0.7, col = col, cl.cex = 0.7,
         order = "hclust", type = "upper", diag = T, mar = c(1, 1, 1, 1))



## Get US county sf file
year=2018
options(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
county_map = get_acs(geography = "county", 
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)

# #join data with sf file
# county_outages = county_map %>%
#   left_join(outages_csv, by = c("GEOID" = "fips_code"))





#### DATA ###################################################################################


## Root zone data
TX_gdb =  "/Users/paulmj/Downloads/gSSURGO_TX/gSSURGO_TX.gdb"
# TX = sf::st_read(dsn = TX_gbd_ext, layer = "MUPOLYGON")
# TX_group = TX %>%
#   group_by(MUKEY) %>%
#   summarise(n = n())
# save(TX_group, file = "TX_group.Rda")
load("Data/TX_group.Rda")
TX_Valu1 = sf::st_read(dsn = TX_gdb, layer = "Valu1")
TX_group_val1 = TX_group %>% left_join(TX_Valu1, by = c("MUKEY" = "mukey"))

TX_rast_1k = st_rasterize(TX_group_val1["rootznemc"], dx = 1000, dy = 1000) #1km x 1km resolution 
TX_rast_100m = st_rasterize(TX_group_val1["rootznemc"], dx = 100, dy = 100) #100m x 100km resolution 

rr = ggplot() + 
  geom_stars(data = TX_rast_1k, aes(x = x, y = y, fill = rootznemc)) + 
  scale_fill_viridis_c(direction = -1, na.value = "gray") +
  theme_minimal()
rr

st_crs(TX_group_val1)
TX_map = st_transform(TX_group_val1, crs_map)

### Geometry base - US Counties 
year=2019
options(tigris_use_cache = TRUE) #to cache shapefiles for future sessions
state = "Texas"
county = c("Harris County")

county_map = get_acs(geography = "county", state = state,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>% select(GEOID, NAME) 

census_map = get_acs(geography = "tract", state = state, county = county,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
census_map = census_map %>% select(GEOID, NAME) 


## PROJECTION (USA Equal Area Conic: EPSG 5070) 
Houston_census = st_transform(census_map, crs = st_crs(TX_rast_100m))
TX_rast_100m_crop = st_crop(TX_rast_100m, Houston_census) #crop to Harris county 
TX_rast_100m_crop_mask = st_crop(TX_rast_100m_crop, Houston_census, crop = FALSE) #mask to Harris county

hh = ggplot() + 
  geom_stars(data = TX_rast_100m_crop_mask, aes(x = x, y = y, fill = rootznemc)) + 
  scale_fill_viridis_c(direction = -1, na.value = "transparent") +
  geom_sf(data = Houston_census, colour = alpha("white", 0.5), fill = NA, size = 0.45) +
  theme_dark() +
  labs(title = "Harris County: Root Zone Depth", fill = "Root Zone\n(cm)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        )
hh

## TREE DATA
tree_db = "/Users/paulmj/Downloads/L48_Totals/L48_Totals.gdb"
st_layers(tree_db)


crs_map = st_crs(census_map)
st_crs(census_map)$IsGeographic  
st_crs(census_map)$units_gdal
st_crs(census_map)$proj4string

## SOIL MOISTURE (monthly) 
soil_db = "/Users/paulmj/Downloads/NDLAS_NOAH_2018_monthly/NLDAS_NOAH0125_M.A201801.002.grb"
asd = read_stars(soil_db, along = "band")
st_dimensions(asd)

library(rNOMADS)
GribInfo(soil_db, file.type = "grib1")
soil0_10 = asd[ , , , 26]
soil10_40 = asd[ , , , 27]
soil40_100 = asd[ , , , 28]


