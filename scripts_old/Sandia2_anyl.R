#### LOAD PACKAGES and SET DIRECTORY 
options(java.parameters = "-Xmx7g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(lme4)
library(corrplot)
library(viridis)
library(tidyverse)
library(tidycensus) 
library(sf)
#library(sqldf)

setwd("~/Documents/01_VECTOR.nosync/Sandia")

#### LOAD VARIABLES
load("./Data/SandiaVariables.Rda") #from Sandia1_vars.R

#### POWER DATA
#outages_csv = read.csv("./Data/SE_states_outage_merra_2018.csv", header = T)
#save(outages_csv, file = "./Data/outages.Rda")
load(file = "Data/outages.Rda")

# Format FIPS to character and include starting zero
outages_csv$fips_code = as.character(outages_csv$fips_code)
fips = outages_csv$fips_code
fips0 = str_pad(fips, 5, pad = "0")
outages_csv$fips_code = fips0

# Add day and month to data
county_map_outages = outages_csv %>%
  mutate(rowID = row_number()) %>% # add row number to index events 
  mutate(date_day = str_extract(date_hour, "^.{10}")) %>% #returns first 10 characters of date_time (yyyy-mm-dd)
  mutate(date_month = str_extract(date_hour, "^.{7}")) # returns yyyy-mm

#### PROCESS POWER DATA
# Filter dates that are flagged as outages 
county_map_outages_FILTER = county_map_outages %>%
  filter(outage_status %in% c("pre", "start", "during", "end")) #filter events flagged as outages

# Join power outages by GEOID and time-stamp 
county_map_outages_JOIN = county_map_area %>%
  inner_join(county_map_outages_FILTER, by = c("GEOID" = "fips_code")) %>% # join outage/hourly weather data 
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) %>% #join socio-economic variables  
  left_join(county_map_soil_CLEAN, by = c("GEOID", "date_hour")) %>% #join soil moisture by county and hour time-stamp 
  left_join(county_map_spi_CLEAN, by = c("GEOID", "date_month")) #join SPI by GEOID and month 

####################################################################################################
#### ANALYSIS 1: SNAPSHOT BY EVENTS (i.e., independent variables correspond to one outage event)
####################################################################################################
# Group by event and aggregate dynamic variables 
outages_group = county_map_outages_JOIN %>%
  st_set_geometry(NULL) %>%
  group_by(outage_number, GEOID) %>%
  summarise(out_hrs= n(), out_maxcust = max(hr_mean_customers_out), out_percust = sum(hr_mean_customers_out) / sum(POPULATION),
            Density = mean(DENSITY), 
            PS_mean = mean(PS), PS_sd = sd(PS), #mean and sd of surface pressure
            SLP_mean = mean(SLP), PS_sd = sd(SLP), # same for sea level pressure
            QV_max = max(QV10M), # max specific humidity
            U_max = max(U10M),  
            V_max = max(V10M),
            WIND_max = max(WIND10M), 
            T_mean = mean(T10M), T_sd = sd(T10M), 
            TQI_mean = mean(TQI), TQL_mean = mean(TQL), TQV_mean = mean(TQV),
            soil10_mean = mean(soil0_10), soil100_mean = mean(soil40_100),
            spi03_mean = mean(spi_03), spi12_mean = mean(spi_12),spi24_mean = mean(spi_24)
  )

# Join back in the static environmental and socio-economic variables
county_outages_GROUP = outages_group %>%
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) #join socio-economic variables 


##### MODELING PART I: ln(hrs) #######################
### A) Big outage events (90th percentile in duration)
df_bart = data.frame(county_outages_GROUP) %>%
  drop_na %>% 
  dplyr::filter(out_hrs > quantile(county_outages_GROUP$out_hrs, .9)) #filter to big events 
y = log(df_bart$out_hrs) #take log to help deal with extreme values 
# 1. With just socio-economic, static environmental, SPI, and Soil Moisture data (i.e., publicly available)
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) %>%
  dplyr::select(-c(PS_mean:TQV_mean))
model_name = paste("Large Events - Static and Socio-economic Variables")
#2. With soci-economic, static, SPI, Soil Moisture, and dynamic weather variables
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) 
model_name = paste("Large Events - Static, Socio-economic, and Dynamic Variables")
#3. With soci-economic, static, SPI, Soil Moisture, dynamic weather variables, and max outages (Liu et al., 2007)
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_percust) 
model_name = paste("Large Events - Static, Socio-economic, and Dynamic Variables (with Outages)")
#4. Two-step model: model forecasted max custmer outages and then use this as an input (Liu et al., 2007)
# y1 = log(df_bart$out_maxcust)
# X1 = df_bart %>%
#   dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) 
# bart1 = bartMachine(X1, y1)
# model_name = paste("Large Events - Customer Outages")
# y1_predictions = predict(bart1, X1)
load(file = "y1_predictions.Rda")
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_percust, -out_maxcust) %>%
  bind_cols(data.frame(y1_predictions))
model_name = paste("Large Events - Static, Socio-economic, and Dynamic Variables (with Frcst Outages)")


### B) Small outage events (<=90th percentile in duration)
df_bart = data.frame(county_outages_GROUP) %>%
  drop_na %>% 
  dplyr::filter(out_hrs <= quantile(county_outages_GROUP$out_hrs, .9)) #filter to small events 
y = log(df_bart$out_hrs) #take log to help deal with extreme values 
# 1. With just socio-economic, static environmental, SPI, and Soil Moisture data (i.e., publicly available)
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) %>%
  dplyr::select(-c(PS_mean:TQV_mean))
model_name = paste("Small Events - Static and Socio-economic Variables")
#2. With soci-economic, static, SPI, Soil Moisture, and dynamic weather variables
#rm(list=setdiff(ls(), "county_outages_GROUP")) 
#gc() 
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_maxcust, -out_percust) 
model_name = paste("Small Events - Static, Socio-economic, and Dynamic Variables")
#3. With soci-economic, static, SPI, Soil Moisture, dynamic weather variables, and max outages (Liu et al. 2007)
#rm(list=setdiff(ls(), "county_outages_GROUP")) 
#gc() 
X = df_bart %>%
  dplyr::select(-outage_number, -GEOID, -out_hrs, -out_percust) 
model_name = paste("Small Events - Static, Socio-economic, and Dynamic Variables (with Outages)")


#### PLOTTING ##########################################################################################
bart = bartMachine(X, y) 
predictions = predict(bart, X)
CI = round(calc_credible_intervals(bart, X, ci_conf = 0.95), 2)
gg = dplyr::tibble(predictions = predictions,
                   lower = CI[,1],
                   upper = CI[,2],
                   actual = y
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
gg$Ymean = mean(gg$actual, na.rm = T)

rmse = sqrt(mean((gg$predictions - gg$actual)^2)) 
rsq = 1 - sum((gg$actual - gg$predictions)^2) / sum((gg$actual - mean(gg$actual))^2) 
lb1 = paste("R^2 == ", round(rsq, 3))


plot_filtering_estimates2 <- function(df) {
  p <- ggplot(data = gg, aes(x = index)) +
    theme_classic() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "indianred"), alpha = 0.5) + #credible intervals 
    geom_line(aes(y = predictions, colour = "red"), size = 0.25, alpha = .9) + #prediction point estimate
    geom_point(aes(y = actual, colour = "black"), size = 0.55, shape = 16, alpha = 0.9) + #actual observation points
    geom_line(aes(y = Ymean, colour = "blue"), size = 0.55, lty = "solid", alpha = 0.9) + #null model (mean only) 
    ylab("Outage Duration (log hours)") + 
    #ylab("Max Cust. Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle(model_name) +
    scale_fill_identity(name = "fill", guide = 'legend', labels = c('95% CI')) + 
    scale_colour_manual(name = 'colour',
                        values = c('black',"blue", "red"),
                        labels = c('Actual', 'Mean', 'BART'),
                        guide = guide_legend(
                          reverse = T,
                          override.aes = list(
                            linetype = c("solid", "solid","blank"),
                            shape = c(NA, NA, 16))
                        )) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.direction = "vertical",
          legend.box = "horizontal",
          legend.position = c(.25, .75) #x and y percentage
          ) +
    annotate("text", x = quantile(gg$index, 0.8), y = quantile(gg$actual, .05), label=lb1, parse=T, color="red", size = 3) 
  print(p)
}
plot_filtering_estimates2(gg)

check_bart_error_assumptions(bart)
plot_convergence_diagnostics(bart)

# Check spatial dependency of residuals
rr = data.frame(as.factor(df_bart$GEOID))
rr$resid = bart$residuals
colnames(rr) = c("GEOID", "resid")
plot(rr) # looks good but use varigoram below to confirm 
library(gstat)
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
county_centroid = st_centroid(county_map) # get center of counties 
county_lonlat = county_centroid %>% 
  mutate(X = unlist(map(county_centroid$geometry,1)),
       Y = unlist(map(county_centroid$geometry,2))) %>%
  dplyr::select(-NAME, -POPULATION) %>%
  inner_join(rr, by = c("GEOID")) %>%
  rename(Z = resid)
county_lonlat_sp = as_Spatial(county_lonlat)
vgram = variogram(Z~1, county_lonlat_sp)
plot(vgram)

#### FEATURE SELECTION ##################################################################################
# Kapelner and Bleich, 2016 
vs = var_selection_by_permute(bart, bottom_margin = 10, num_reps_for_avg = 20, num_permute_samples = 20, plot = F)
vs$important_vars_local_names
vs$important_vars_global_se_names

ii = interaction_investigator(bart, num_replicates_for_avg = 20, num_var_plot = 10, bottom_margin = 15)

## Model I.A.1
# > vs$important_vars_local_names
# [1] "soil10_mean"  "soil100_mean" "spi03_mean"   "CROPINS"      "spi24_mean"   "Wetlands"    
# [7] "MNTHLTH"      "spi12_mean"   "PROXCAP"      "WATEFF"      
# > vs$important_vars_global_se_names
# [1] "soil10_mean"  "soil100_mean" "spi03_mean"   "CROPINS"      "spi24_mean"  

## Model I.A.2
# > vs$important_vars_local_names
# [1] "SLP_mean"    "PS_sd"       "QV_max"      "T_sd"        "WIND_max"    "V_max"      
# [7] "TQI_mean"    "U_max"       "TQL_mean"    "TQV_mean"    "PS_mean"     "spi03_mean" 
# [13] "T_mean"      "soil10_mean"
# > vs$important_vars_global_se_names
# [1] "SLP_mean"   "PS_sd"      "QV_max"     "T_sd"       "WIND_max"   "V_max"      "TQI_mean"  
# [8] "U_max"      "spi03_mean"

## Model I.A.3
# > vs$important_vars_local_names
# [1] "SLP_mean"    "PS_sd"       "out_maxcust" "V_max"       "QV_max"      "T_sd"       
# [7] "WIND_max"    "U_max"       "TQI_mean"    "TQL_mean"    "PS_mean"     "spi03_mean" 
# [13] "TQV_mean"    "CROPINS"    
# > vs$important_vars_global_se_names
# [1] "SLP_mean"    "PS_sd"       "out_maxcust" "V_max"       "QV_max"      "T_sd"       
# [7] "WIND_max"    "U_max"       "TQI_mean"    "TQL_mean"    "spi03_mean" 

## Model I.A.4
# > vs$important_vars_local_names
# [1] "SLP_mean"       "PS_sd"          "y1_predictions" "V_max"          "QV_max"        
# [6] "T_sd"           "WIND_max"       "TQI_mean"       "U_max"          "PS_mean"       
# [11] "TQL_mean"       "TQV_mean"       "Density"       
# > vs$important_vars_global_se_names
# [1] "SLP_mean"       "PS_sd"          "y1_predictions" "V_max"          "QV_max"        
# [6] "T_sd"           "WIND_max"       "TQI_mean"       "U_max"      

## Model Customer Outages 
# > vs$important_vars_local_names
# [1] "PS_sd"     "SLP_mean"  "TQI_mean"  "WIND_max"  "T_sd"      "QV_max"    "BUSINORGS"
# [8] "T_mean"    "Developed" "U_max"     "EMPSCIENC" "Density"   "V_max"     "PROPINSUR"
# [15] "TQL_mean"  "ENVORGS"  
# > vs$important_vars_global_se_names
# [1] "PS_sd"     "SLP_mean"  "WIND_max"  "T_sd"      "BUSINORGS" "T_mean"   

