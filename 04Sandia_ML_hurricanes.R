#### Machine learning 
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

################################################################################
#### PRE-PROCESSING ############################################################
################################################################################
# Parallel processing setup
num_cores = detectCores() - 1
unregister_dopar = function() { #function to un-register parallel processing in doParallel
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

our_events = c(
  "droughts",
  "extreme_cold",
  "extreme_heat",
  "floods",
  #"hail", #no events
  "high_winds",
  "hurricanes",
  #"tornadoes", #no events
  "wildfires",
  "winter_storms"
)

# Load data and select final DVs and GROUP by GEOID
load(file = "Data/processed/sf_data_ALL_nototal.Rda")

# Summarize events by counts (or total hours)
df_data_EVENTS = sf_data_ALL %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-c(mean_cust_out:max_frac_cust_out, max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  group_by(GEOID) %>% 
  #summarise(across(all_of(our_events), ~ sum(.x, na.rm = TRUE))) #total counts of each event
  summarise(across(all_of(our_events), ~ sum(.x * duration_hr, na.rm = TRUE))) #total hrs of each event

# Data-frame of static predictors 
load(file = "Data/processed/county_map_proj.Rda")
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) %>% #population per sq-km 
  dplyr::select(-c(AREA, NAME)) 
load(file = "Data/processed/county_proj_nlcd.Rda")
county_map_nlcd = county_map_nlcd %>%
  dplyr::select(-c(NAME, POPULATION, ID, Other, Total)) %>%
  st_set_geometry(NULL)
load(file = "Data/processed/county_proj_dem.Rda")
county_map_dem = county_proj_dem %>%
  select(GEOID, DEM_mean, DEM_med, DEM_sd, DEM_min, DEM_max) %>%
  st_set_geometry(NULL)
load(file = "Data/processed/county_proj_rz.Rda")
county_map_rz = county_map_rz %>%
  select(GEOID, RZ_mean, RZ_med, RZ_mode) %>%
  st_set_geometry(NULL)
load(file = "Data/processed/county_proj_social.Rda") 
county_map_social = county_proj_social %>%
  dplyr::select(-NAME, -POPULATION) %>%
  st_set_geometry(NULL)
county_map_JOIN = county_map_area %>%
  left_join(county_map_nlcd, by = c("GEOID")) %>%
  left_join(county_map_dem, by = c("GEOID")) %>%
  left_join(county_map_rz, by = c("GEOID")) %>%
  left_join(county_map_social, by = c("GEOID"))

clean_recipe = recipe(duration_hr ~ . , data = sf_pred) %>%
  step_rm(contains("total")) %>% #remove total_ice and total_liquid 
  step_impute_knn(all_predictors()) 










# Summarize other variables by just taking first value (all the same)
df_data_OTHER = sf_data_ALL %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-duration_hr) %>%
  dplyr::select(-c(mean_cust_out:max_frac_cust_out, max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  dplyr::select(-c(all_of(our_events), "avalanche")) %>%
  group_by(GEOID) %>% 
  summarise(across(everything(), first))

# Combine summarized dataframes 
df_data_CLEAN = df_data_EVENTS %>%
  mutate(across(all_of(our_events), ~ .x / 5 )) %>% #annualize
  #mutate(across(all_of(our_events), ~ jitter(log(.x / 5 + 0.001)))) %>% #annualize, take log, and add noise
  inner_join(df_data_OTHER, by = "GEOID")

# Split into training vs testing
set.seed(23)
df_split = initial_split(df_data_CLEAN, prop = 0.80, strata = "hurricanes")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)
df_preds = df_data_CLEAN %>% dplyr::select(-c(GEOID, POPULATION, all_of(our_events)))

# Recipes for different events 
recipe_hurricanes = recipe(hurricanes ~ . , data = df_data_CLEAN) %>% 
  step_rm(GEOID, POPULATION, droughts, extreme_cold, extreme_heat, floods, high_winds, wildfires, winter_storms) %>% 
  step_naomit(hurricanes)

recipe_winter = recipe(winter_storms ~ . , data = df_data_CLEAN) %>% 
  step_rm(GEOID, POPULATION, droughts, extreme_cold, extreme_heat, floods, high_winds, wildfires, hurricanes) %>% 
  step_naomit(winter_storms)

recipe_highwind = recipe(high_winds ~ . , data = df_data_CLEAN) %>% 
  step_rm(GEOID, POPULATION, droughts, extreme_cold, extreme_heat, floods, winter_storms, wildfires, hurricanes) %>% 
  step_naomit(high_winds)

recipe_droughts = recipe(droughts ~ . , data = df_data_CLEAN) %>% 
  step_rm(GEOID, POPULATION, hurricanes, extreme_cold, extreme_heat, floods, high_winds, wildfires, winter_storms) %>% 
  step_naomit(droughts)

recipe_floods = recipe(floods ~ . , data = df_data_CLEAN) %>% 
  step_rm(GEOID, POPULATION, droughts, extreme_cold, extreme_heat, hurricanes, high_winds, wildfires, winter_storms) %>% 
  step_naomit(floods)

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
### Define which recipe you want to use 
recipe_mine = recipe_hurricanes

y_train = df_train %>% pull(hurricanes) %>% na.omit() 
y_test = df_test %>% pull(hurricanes) %>% na.omit() 
X_train = df_train %>% dplyr::select(-c(GEOID, POPULATION, all_of(our_events))) %>% na.omit()
X_test = df_test %>% dplyr::select(-c(GEOID, POPULATION, all_of(our_events))) %>% na.omit()

### Lasso, Ridge Regression, and Elastic Net ###################################
#https://www.tidyverse.org/blog/2020/11/tune-parallel/
show_model_info("linear_reg")
lre_model = linear_reg(penalty = tune(), mixture = tune()) %>% #lambda (penalty) and alpha/mixture (1 lasso, 0 ridge)
  set_engine("glmnet") %>%
  translate()
lre_work = workflow() %>% 
  add_recipe(recipe_mine) %>%
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

rsq_lre = paste(lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
cverror_lre = paste(show_best(lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

### XGBoost ####################################################################
show_model_info("boost_tree")
xgb_model = boost_tree(mode = "regression", trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), mtry = tune()) %>%
  set_engine(engine = "xgboost") %>%
  translate()
xgb_work = workflow() %>%
  add_recipe(recipe_mine) %>%
  add_model(xgb_model)
set.seed(32); xgb_grid = dials::grid_max_entropy(parameters(trees(), min_n(), tree_depth(), learn_rate(), loss_reduction(), finalize(mtry(), df_preds)), size = 100)
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

rsq_xgb = paste(xgb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
cverror_xgb = paste(show_best(xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
  
### PLOTTING - Duration Hours
gg = dplyr::tibble(actual = y_test, 
                   eNet = as.vector(lre_predictions$.pred),
                   #rf = as.vector(rf_predictions$.pred),
                   xgb = as.vector(xgb_predictions$.pred)
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
#gg = gg %>% filter(eNet < 20)
#gg$mean = mean(gg$actual, na.rm = T)
gg_actual = gg %>% dplyr::select(index, actual)
gg_pred = gg %>% dplyr::select(-actual) %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
gg_all = gg %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
color_vec= c("black", "#FDE725FF", "#35B779FF")
#color_vec = c("black", "#FDE725FF", "#35B779FF", "#440154FF")
lty_vec = c(1, 1, 1)
#lty_vec = c(1, 1, 1, 1)
alpha_vec = c(1, 0.7, 0.7)
#alpha_vec = c(1, 0.7, 0.7, 0.7)

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
    ylab("Expected Outage Hours per Year (Hurricanes)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("County Index") +
    #ggtitle("County-Level Predictions: Test Sample\nNo 'Total' Variables") + 
    ggtitle("County-Level Predictions: Test Sample") + 
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
#https://towardsdatascience.com/be-careful-when-interpreting-your-features-importance-in-xgboost-6e16132588e7
#https://stackoverflow.com/questions/57360703/feature-importance-gain-in-xgboost
#https://xgboost.readthedocs.io/en/stable/R-package/discoverYourData.html
#https://bgreenwell.github.io/pdp/articles/pdp-example-xgboost.html
#https://bgreenwell.github.io/pdp/articles/pdp-extending.html
#https://christophm.github.io/interpretable-ml-book/pdp.html
#https://christophm.github.io/interpretable-ml-book/ice.html
#https://medium.com/dataman-in-ai/how-is-the-partial-dependent-plot-computed-8d2001a0e556
#https://datascience.stackexchange.com/questions/15305/how-does-xgboost-learn-what-are-the-inputs-for-missing-values
final_lre = lre_work %>%
  finalize_workflow(lre_best) %>%
  fit(df_data_CLEAN)
print(tidy(final_lre), n = 150)

final_model = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_data_CLEAN)
final_obj = extract_fit_parsnip(final_model)$fit
importance = xgb.importance(model = final_obj)
head(importance, n = 20) 
vip(final_obj, n = 20)

##########################################################################################################
##### PARTIAL DEPENDENCE PLOTS (PDPs) - FINAL MODEL ######################################################
##########################################################################################################
X_train = df_train %>% dplyr::select(-c(GEOID, POPULATION, all_of(our_events)))
pdp1 = pdp::partial(final_obj, pred.var = "Wetlands", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp2 = pdp::partial(final_obj, pred.var = "DEM_sd", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp3= pdp::partial(final_obj, pred.var = "DEM_min", ice = F, center = F,
                   plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                   train = X_train, type = "regression")
pdp4 = pdp::partial(final_obj, pred.var = "QMOHO", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp5 = pdp::partial(final_obj, pred.var = "GENDINC", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp6 = pdp::partial(final_obj, pred.var = "QBLACK", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp7 = pdp::partial(final_obj, pred.var = "PATTACHRES", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp8 = pdp::partial(final_obj, pred.var = "QHLTH65", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")

grid.arrange(pdp1, pdp2, pdp3, pdp4, ncol = 2)
grid.arrange(pdp5, pdp6, pdp7, pdp8, ncol = 2)


##########################################################################################################
##### CREATE MAPS  #######################################################################################
##########################################################################################################
# Base map
load(file = "Data/processed/county_map_proj.Rda")
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) %>% #population per sq-km 
  dplyr::select(-c(AREA, NAME)) 

# NLCD (static)
load(file = "Data/processed/county_proj_nlcd.Rda")
county_map_nlcd = county_map_nlcd %>%
  dplyr::select(-c(NAME, POPULATION, ID, Other, Total)) %>%
  st_set_geometry(NULL)

# DEM (static)
load(file = "Data/processed/county_proj_dem.Rda")
county_map_dem = county_proj_dem %>%
  select(GEOID, DEM_mean, DEM_med, DEM_sd, DEM_min, DEM_max) %>%
  st_set_geometry(NULL)

# Root Zone (static)
load(file = "Data/processed/county_proj_rz.Rda")
county_map_rz = county_map_rz %>%
  select(GEOID, RZ_mean, RZ_med, RZ_mode) %>%
  st_set_geometry(NULL)

# Socio-Economic Data
load(file = "Data/processed/county_proj_social.Rda") 
county_map_social = county_proj_social %>%
  dplyr::select(-NAME, -POPULATION) %>%
  st_set_geometry(NULL)

# Join data frames 
county_map_JOIN = county_map_area %>%
  left_join(county_map_nlcd, by = c("GEOID")) %>%
  left_join(county_map_dem, by = c("GEOID")) %>%
  left_join(county_map_rz, by = c("GEOID")) %>%
  left_join(county_map_social, by = c("GEOID"))

# Select GEOID and predictors in final mobel (final_obj)
df_join = county_map_JOIN %>%
  dplyr::select(GEOID, final_obj$feature_names) %>%
  st_set_geometry(NULL)

# Impute missing predictors 
impute_recipe = recipe(GEOID ~ ., data = df_join) %>%
  step_impute_knn(all_predictors())

# Matrix of predictors
X = prep(impute_recipe) %>%
  juice() %>%
  dplyr::select(-GEOID) %>%
  as.matrix()

# Predict 
predictions = predict(final_obj, X)

# Join predictions to SF
sf_map = county_map_JOIN %>%
  dplyr::select(GEOID, POPULATION) %>%
  mutate(pred = predictions)

gg3 = ggplot()+
  geom_sf(data = sf_map, aes(fill = pred), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  geom_sf(data = county_map_proj, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  labs(title = "Expected Annual Power Outages \n-- Hurricanes --", fill = "Hours") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)


##########################################################################################################
##### HURRICANE EXPOSURE  ################################################################################
##########################################################################################################
## https://cran.r-project.org/web/packages/hurricaneexposure/vignettes/hurricaneexposure.html
county_list = sf_map$GEOID
hurricane_exposure = county_wind(counties = c(county_list), start_year = 2015, end_year = 2019, wind_limit = 0.01)
hurricane_summ = hurricane_exposure %>%
  group_by(fips) %>%
  summarize(AVG_WIND = mean(vmax_sust))
hurricane_map = sf_map %>%
  left_join(hurricane_summ, by = c("GEOID" = "fips"))

gg3 = ggplot()+
  geom_sf(data = hurricane_map, aes(fill = AVG_WIND), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  geom_sf(data = county_map_proj, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  labs(title = "Hurricane Exposure \n2015 - 2019", fill = "Mean \nWind (m/s)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)

cor(hurricane_map$AVG_WIND, hurricane_map$pred, use = "complete.obs")
# 0.635

## Hurricane deep dive
county_zoom_list = c("FL", "GA", "AL", "SC", "NC")
county_zoom = get_acs(geography = "county", state = county_zoom_list,
                      variables=c("B01003_001"), year = 2019, geometry = TRUE, 
                      cache_table = TRUE) %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME) %>%
  st_transform(5070) 

county_zoom_map = county_zoom %>%
  inner_join((st_set_geometry(hurricane_map, NULL)), by = "GEOID") %>%
  inner_join(df_join, by = "GEOID")

gg4 = ggplot()+
  geom_sf(data = county_zoom_map, aes(fill = QBLACK), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  geom_sf(data = county_zoom_map, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  labs(title = "Hurricane Exposure (2015 - 2019)", fill = "Mean \nWind (m/s)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg4)
