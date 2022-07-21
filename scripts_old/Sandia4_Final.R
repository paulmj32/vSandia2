#### LOAD PACKAGES and SET DIRECTORY 
library(tidyverse)
library(tidymodels)
library(tidycensus) 
library(sf)
library(corrplot)
library(viridis)
library(doParallel)
library(terra)
library(stars)
library(raster) #make sure ncdf4 package is installed 
library(lubridate)
library(spdep)
library(xgboost)
library(DBI)
library(vip)
library(pdp)
tidymodels_prefer()

setwd("~/Documents/01_VECTOR.nosync/Sandia")

#################################################################################################################
###### DATA #####################################################################################################
#################################################################################################################
## Get SQL power outage data 
# https://stackoverflow.com/questions/9802680/importing-files-with-extension-sqlite-into-r
# https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html
# https://www.r-bloggers.com/2018/11/quick-guide-r-sqlite/
mydb = RSQLite::dbConnect(drv = RSQLite::SQLite(), 
                          dbname = "./Data/Outages sqlite/AGM_full_CONUS_outages_merra_noaa_summary_Proc31Jan2022.sqlite")
tables = dbListTables(mydb) #see tables
myquery = dbGetQuery(mydb, 'SELECT * FROM summary')
dbDisconnect(mydb)
county_outages = myquery %>%
  mutate(start_dt = as.character(as_datetime(start_date))) %>%
  mutate(end_dt = as.character(as_datetime(end_date))) %>%
  mutate(start_ym = str_extract(start_dt, "^.{7}")) %>% 
  mutate(end_ym = str_extract(end_dt, "^.{7}")) %>% 
  dplyr::select(-c(start_date, end_date)) %>%
  relocate(c(start_dt, start_ym, end_dt, end_ym), .after = outage_number)

#asd = dbSendQuery(con, statement)
#asd2 = dbFetch(asd, 10) #to get first 10

## Get census map
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
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) %>% #population per sq-km 
  dplyr::select(-c(AREA, NAME))

## Get environmental data
# NLCD (static)
load(file = "./Data/county_map_nlcd.Rda")
county_map_nlcd = county_map_nlcd %>%
  dplyr::select(-c(NAME, POPULATION, ID, Other, Total)) %>%
  st_set_geometry(NULL)

# DEM (static)
load(file = "./Data/county_map_dem.Rda")
county_map_dem = county_map_dem %>%
  select(GEOID, DEM_mean:DEM_max) %>%
  st_set_geometry(NULL)

# Root Zone (static)
load(file = "./Data/county_map_rz.Rda")
county_map_rz = county_map_rz %>%
  select(GEOID, RZ_mean:RZ_mode) %>%
  st_set_geometry(NULL)

# SPI (monthly)
load(file = "./Data/county_map_spi.Rda")
county_map_spi = county_map_spi %>%
  mutate(date_ym = str_extract(Date, "^.{7}")) %>% # returns yyyy-mm
  dplyr::select(-NAME, -Date)

## Socio-Economic Data
# ## County-level indicators
load(file = "./Data/final_data.Rda") 
county_social = final_data %>%
  dplyr::select(-NAME, -POPULATION) %>%
  st_set_geometry(NULL)

## Join tables
county_map_JOIN = county_map_area %>%
  inner_join(county_outages, by = c("GEOID" = "fips_code")) %>%
  left_join(county_map_spi, by = c("GEOID" = "GEOID", "start_ym" = "date_ym")) %>%
  left_join(county_map_nlcd, by = c("GEOID")) %>%
  left_join(county_map_dem, by = c("GEOID")) %>%
  left_join(county_map_rz, by = c("GEOID")) %>%
  left_join(county_social, by = c("GEOID"))

county_map_ALL = county_map_JOIN %>%
  dplyr::select(-c(outage_number, start_dt, start_ym, end_dt, end_ym)) %>%
  dplyr::select(-c(max_WIND2M, max_WIND50M, min_humidity_2M, mean_humidity_2M, max_humidity_2M,
                   min_T2M, mean_T2M, max_T2M, delta_T2M, min_T10M_C:delta_T10M_C, WETLAND, NUCACC
                   )) #take extraneous variables
  
rm(list=setdiff(ls(), c("county_map_JOIN", "county_map_ALL"))) 
gc() 

##########################################################################################################
#### PRE-PROCESSING ######################################################################################
##########################################################################################################
## Clean data-frame
df_data = county_map_ALL %>%
  #dplyr::filter(duration_hr > quantile(county_map_ALL$duration_hr, .9)) %>% #filter to big events 
  dplyr::filter(duration_hr >= 12) %>% #filter to events >12 hrs ... 95% quantile for full dataset
  dplyr::filter(duration_hr < 3000) %>% #get rid of faulty data
  mutate(ln_cust = log(max_pct_affected * POPULATION), ln_hrs = log(duration_hr)) %>% # take log of DVs
  dplyr::select(-c(max_pct_affected, duration_hr, GEOID, POPULATION)) %>%
  relocate(c(ln_cust, ln_hrs)) %>%
  st_set_geometry(NULL)


num_cores = detectCores() - 1
#function to un-register parallel processing in doParallel
unregister_dopar = function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

## Split into training vs testing
set.seed(23)
df_split = initial_split(df_data, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

## Pre-processing
cust_recipe = recipe(ln_cust ~ . , data = df_data) %>%
  step_rm(ln_hrs, contains("total")) %>% 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  # step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors())  #removes highly correlated 
cust_prep = prep(cust_recipe) #see changes 
df_cust = prep(cust_recipe) %>% juice() #apply recipe to data frame 

hours_recipe = recipe(ln_hrs ~ . , data = df_data) %>%
  #step_rm(ln_cust) %>% 
  step_rm(ln_cust, contains("total")) %>% 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  # step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors())  #removes highly correlated 
hours_prep = prep(hours_recipe) #shows changes 
df_hours = prep(hours_recipe) %>% juice() #apply recipe to data frame

##########################################################################################################
#### MACHINE LEARNING ####################################################################################
##########################################################################################################
#X = df_hours %>% dplyr::select(-ln_hrs) %>% as.data.frame()
X = df_cust %>% dplyr::select(-ln_cust) %>% as.data.frame() 
#y = df_hours %>% dplyr::select(ln_hrs) %>% pull()
y = df_cust %>% dplyr::select(ln_cust) %>% pull()

#https://www.tidyverse.org/blog/2020/11/tune-parallel/
## Lasso/Ridge/ElasticNet 
show_model_info("linear_reg")
lre_model = linear_reg(penalty = tune(), mixture = tune()) %>% #lambda (penalty) and alpha/mixture (1 lasso, 0 ridge)
  set_engine("glmnet") %>%
  translate()
lre_work = workflow() %>% 
  add_recipe(cust_recipe) %>%
  #add_recipe(hours_recipe) %>%
  add_model(lre_model)
#lre_grid = dials::grid_regular(parameters(penalty(), mixture()), levels = c(5, 5))
set.seed(32); lre_grid = dials::grid_max_entropy(parameters(penalty(), mixture()), size = 40)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
lre_tune = lre_work %>%
  tune_grid(resamples = df_cv,
            grid = lre_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")
            ) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(lre_tune, metric = "rmse")
lre_tune_results = lre_tune %>% collect_metrics()
lre_best = lre_tune %>% select_best(metric = "rmse")
lre_fit = lre_work %>%
  finalize_workflow(lre_best) %>%
  last_fit(df_split)
lre_test = lre_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
lre_predictions = lre_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)

## Random Forest
show_model_info("rand_forest")
rf_model = rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression") %>%
  translate()
rf_work = workflow() %>%
  #add_recipe(cust_recipe) %>%
  add_recipe(hours_recipe) %>%
  add_model(rf_model)
#rf_grid = expand.grid(mtry = c(5, 10, 20), trees = c(100, 500, 1000), min_n = c(3, 5, 10))
set.seed(32); rf_grid = dials::grid_max_entropy(parameters(min_n(), trees(), finalize(mtry(), df_cust)), size = 60)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
rf_tune = rf_work %>%
  tune_grid(resamples = df_cv,
            grid = rf_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T)#, allow_par = T, parallel_over = "resamples") #parallel processing turns off verbose
  )
stopCluster(cl) 
unregister_dopar()
show_best(rf_tune, metric = "rmse")
rf_tune_results = rf_tune %>% collect_metrics()
rf_best = rf_tune %>% select_best(metric = "rmse")
rf_fit = rf_work %>%
  finalize_workflow(rf_best) %>%
  last_fit(df_split)
rf_test = rf_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
rf_predictions = rf_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)

## XGBoost
show_model_info("boost_tree")
xgb_model = boost_tree(mode = "regression", trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), mtry = tune()) %>%
  set_engine(engine = "xgboost") %>%
  translate()
xgb_work = workflow() %>%
  add_recipe(cust_recipe) %>%
  #add_recipe(hours_recipe) %>%
  add_model(xgb_model)
set.seed(32); xgb_grid = dials::grid_max_entropy(parameters(trees(), min_n(), tree_depth(), learn_rate(), loss_reduction(), finalize(mtry(), df_hours)), size = 100)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
xgb_tune = xgb_work %>%
  tune_grid(resamples = df_cv, 
            grid = xgb_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")
            ) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(xgb_tune, metric = "rmse")
xgb_tune_results = xgb_tune %>% collect_metrics()
xgb_best = xgb_tune %>% select_best(metric = "rmse")
xgb_fit = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  last_fit(df_split)
xgb_test = xgb_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
xgb_predictions = xgb_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)

## Testing results
#y_test = df_hours %>% slice(-df_split$in_id) %>% dplyr::select(ln_hrs) %>% pull()
y_test = df_cust %>% slice(-df_split$in_id) %>% dplyr::select(ln_cust) %>% pull()

rsq_lre = paste(lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_xgb = paste(xgb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_rf = paste(rf_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))

cverror_lre = paste(show_best(lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_xgb = paste(show_best(xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_rf = paste(show_best(rf_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

### PLOTTING 
gg = dplyr::tibble(actual = y_test,
                   eNet = as.vector(lre_predictions$.pred),
                   #rf = as.vector(rf_predictions$.pred),
                   xgb = as.vector(xgb_predictions$.pred)
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
#gg$mean = mean(gg$actual, na.rm = T)
gg_actual = gg %>% dplyr::select(index, actual)
gg_pred = gg %>% dplyr::select(-actual) %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
gg_all = gg %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")


color_vec= c("black", "#FDE725FF", "#35B779FF")
color_vec = c("black", "#FDE725FF", "#35B779FF", "#440154FF")
lty_vec = c(1, 1, 1)
lty_vec = c(1, 1, 1, 1)
alpha_vec = c(1, 0.7, 0.7)
alpha_vec = c(1, 0.7, 0.7, 0.7)

# ggplot (for all)
plot_filtering_estimates2 <- function(df) {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
    geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
    scale_color_manual(
      values = color_vec, 
      labels = c("Actual",
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 #bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
                 ),
    name = element_blank()) +
    scale_linetype_manual(
      values = lty_vec,
      labels = c("Actual",
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 #bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")         
                 ),
      name = element_blank()) + 
    scale_alpha_manual(
      values = alpha_vec,
      labels = c("Actual",
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 #bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
                 ),
      name = element_blank()) + 
    ylab("Max Cust Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle("County-Level Predictions: Test Sample\nNo 'Total' Variables") + 
    # guides(
    #   color = guide_legend(order = 2),
    #   shape = guide_legend(order = 1),
    #   linetype = guide_legend(order = 2)
    # ) + 
    theme(legend.spacing.y = unit(-0.25, "cm"),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = c(.225, .8),
          plot.title = element_text(hjust = 0.5)
    )
  
  print(p)
}
plot_filtering_estimates2(gg)


##########################################################################################################
##### VARIABLE IMPORTANCE - FINAL MODEL ##################################################################
##########################################################################################################
#https://www.rebeccabarter.com/blog/2020-03-25_machine_learning/#variable-importance
#https://xgboost.readthedocs.io/en/stable/R-package/discoverYourData.html
#https://bgreenwell.github.io/pdp/articles/pdp-example-xgboost.html
#https://bgreenwell.github.io/pdp/articles/pdp-extending.html
#https://christophm.github.io/interpretable-ml-book/ice.html
# https://datascience.stackexchange.com/questions/15305/how-does-xgboost-learn-what-are-the-inputs-for-missing-values
final_model = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_data)
final_obj = extract_fit_parsnip(final_model)$fit
importance = xgb.importance(model = final_obj)
head(importance, n = 15); vip(final_obj, n = 20)

pdp_NUCACC = pdp::partial(final_obj, pred.var = "NUCACC", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X, type = "regression")
pdp_delpress = pdp::partial(final_obj, pred.var = "delta_pressure", ice = F, center = F,
                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                     train = X, type = "regression")
pdp_minhumid = pdp::partial(final_obj, pred.var = "min_humidity_10M", ice = F, center = F,
                   plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                   train = X, type = "regression")
pdp_totliq = pdp::partial(final_obj, pred.var = "total_liquid", ice = F, center = F,
                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                     train = X, type = "regression")
pdp_deltemp = pdp::partial(final_obj, pred.var = "delta_T10M", ice = F, center = F,
                        plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                        train = X, type = "regression")
grid.arrange(pdp_totvapor, pdp_delpress, pdp_minhumid, pdp_totliq, pdp_deltemp, ncol = 3)


# p1 = pdp::partial(final_obj, pred.var = c("PS_sd", "SLP_mean"), 
#                   plot = T, chull= T, rug = T, plot.engine = "ggplot2",
#                   train = x, type = "regression")
# p2 = pdp::partial(final_obj, pred.var = c("PS_sd", "V_max"), 
#                   plot = T, chull= T, plot.engine = "ggplot2",
#                   train = x, type = "regression")
# p3 = pdp::partial(final_obj, pred.var = c("SLP_mean", "V_max"), 
#                   plot = T, chull= T, plot.engine = "ggplot2",
#                   train = x, type = "regression")

### View tree diagram 
# xgb.plot.tree(model = final_obj, trees = 1:2) 

## MAPS
df_plot = county_map_ALL %>%
  dplyr::filter(duration_hr >= 12) %>% #filter to events >12 hrs ... 95% quantile for full dataset
  dplyr::filter(duration_hr < 3000) #get rid of faulty data
X = df_plot %>%
  st_set_geometry(NULL) %>%
  dplyr::select(final_obj$feature_names) %>%
  as.matrix()
predictions = predict(final_obj, X)
df_plot2 = df_plot %>%
  dplyr::select(GEOID) %>%
  st_set_geometry(NULL) %>%
  mutate(fit = predictions) %>%
  group_by(GEOID) %>%
  summarise(fit_mean = mean(fit))


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
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) %>% #population per sq-km 
  dplyr::select(-c(AREA, NAME))

county_map_risk = county_map_area %>%
  left_join(df_plot2, by = c("GEOID"))
gg2 = ggplot(county_map_risk)+
  geom_sf(aes(fill = fit_mean), color = NA) + 
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  theme_dark() +
  labs(title = "Power Outage Prediction - All Events", fill = "Mean Duration\nln(Hours)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg2)
