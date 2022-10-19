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
#sort(sapply(ls(),function(x){object.size(get(x))}), decreasing = T)
# Parallel processing setup
num_cores = detectCores() - 1
unregister_dopar = function() { #function to un-register parallel processing in doParallel
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
set_bart_machine_num_cores(num_cores = num_cores)

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

# Load data and select final DVs 
load(file = "Data/processed/sf_data_ALL_nototal.Rda")
sf_data = sf_data_ALL %>%
  dplyr::filter(hurricanes >= 1) %>% #filter to event of interest
  dplyr::select(-c(POPULATION, mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out, avalanche)) 
rm(list=c("sf_data_ALL")) 
gc() 

# Summarize events by counts (or total hours)
df_data_EVENTS = sf_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-c(max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  group_by(GEOID) %>% 
  #summarise(across(all_of(our_events), ~ sum(.x, na.rm = TRUE))) #total counts of each event
  summarise(across(all_of(our_events), ~ sum(.x * duration_hr, na.rm = TRUE))) #total hrs of each event

# Dataframe of predictors 
df_data_PRED = sf_data %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-duration_hr) %>%
  dplyr::select(-c(max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  dplyr::select(-c(all_of(our_events))) %>%
  group_by(GEOID) %>% 
  summarise(across(everything(), first))
df_preds = df_data_PRED %>% select(-GEOID)

# Combine summarized dataframes 
df_data = df_data_EVENTS %>%
  #mutate(across(all_of(our_events), ~ .x / 5 )) %>% #annualize
  mutate(across(all_of(our_events), ~ log(.x / 5 + 0.001))) %>% #annualize and take log
  inner_join(df_data_PRED, by = "GEOID")

# Split into training vs testing
set.seed(23)
df_split = initial_split(df_data, prop = 0.80, strata = "hurricanes")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

# Recipes for tidymodels 
recipe_hurricanes = recipe(hurricanes ~ . , data = df_data) %>% 
  step_rm(GEOID, droughts, extreme_cold, extreme_heat, floods, high_winds, wildfires, winter_storms) %>% 
  step_naomit(hurricanes)

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
### Define which recipe you want to use 
recipe_mine = recipe_hurricanes

y_train = df_train %>% pull(hurricanes) %>% na.omit() 
y_test = df_test %>% pull(hurricanes) %>% na.omit() 
X_train = df_train %>% dplyr::select(-c(GEOID, all_of(our_events))) %>% na.omit()
X_test = df_test %>% dplyr::select(-c(GEOID, all_of(our_events))) %>% na.omit()

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

### BART #######################################################################
bart_fit = bartMachineCV(data.frame(X_train), y_train, k_folds = 10, serialize = T) 
bart_predictions = predict(bart_fit, data.frame(X_test))

rsq_bart = format(round(1 - sum((y_test - bart_predictions)^2) / sum((y_test - mean(y_test))^2), 3), nsmall = 3)
bart_cv = k_fold_cv(data.frame(X_train), y_train, k_folds = 10)
cverror_bart = format(round(bart_cv$rmse, 3), nsmall = 3)


##########################################################################################################
##### PLOTTING ###########################################################################################
##########################################################################################################
gg = dplyr::tibble(actual = y_test, 
                   eNet = as.vector(lre_predictions$.pred),
                   bart = as.vector(bart_predictions),
                   #rf = as.vector(rf_predictions$.pred),
                   xgb = as.vector(xgb_predictions$.pred)
)
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
gg_actual = gg %>% dplyr::select(index, actual)
gg_pred = gg %>% dplyr::select(-actual) %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
gg_all = gg %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")

color_vec = c("black", "#FDE725FF", "#35B779FF", "#440154FF")
lty_vec = c(1, 1, 1, 1)
alpha_vec = c(1, 0.6, 0.6, 0.6)

plot_filtering_estimates2 <- function(df) {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
    geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
    scale_color_manual(
      values = color_vec, 
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"),
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
      ),
      name = element_blank()) +
    scale_linetype_manual(
      values = lty_vec,
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"),
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")         
      ),
      name = element_blank()) + 
    scale_alpha_manual(
      values = alpha_vec,
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"),
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
      ),
      name = element_blank()) + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("County") +
    ylab("Annual Outages (ln - hrs)") + 
    ggtitle("Power Outages Due to Tropical Storms: Test Sample") + 
    # guides(
    #   color = guide_legend(order = 2),
    #   shape = guide_legend(order = 1),
    #   linetype = guide_legend(order = 2)
    # ) + 
    theme(legend.spacing.y = unit(-0.25, "cm"),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = c(.225, .85),
          plot.title = element_text(hjust = 0.5)
    )
  
  print(p)
}
plot_filtering_estimates2(gg)


final_model = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_data)
final_obj = extract_fit_parsnip(final_model)$fit
importance = xgb.importance(model = final_obj)
head(importance, n = 20) 
vip(final_obj, n = 20)

##########################################################################################################
##### EXPOSURE MAPS  #####################################################################################
##########################################################################################################
# Hurricane exposure
load(file = "Data/processed/county_map_proj.Rda")
hurricane_exposure = county_wind(counties = as.vector(county_map_proj$GEOID), 
                                 start_year = 2015, end_year = 2019, wind_limit = 17.4)
hurricane_group = hurricane_exposure %>%
  group_by(fips) %>%
  summarize(avg_sust_spd = mean(vmax_sust), avg_gust_spd = mean(vmax_gust), sum_sust_hrs = sum(sust_dur) / 60, sum_gust_hrs = sum(gust_dur) / 60)

# Annual outage hours - all US
hrs = df_data %>%
  dplyr::select(GEOID, hurricanes) %>%
  mutate(hrs = exp(hurricanes)) %>%
  group_by(GEOID) %>%
  summarise(Hours = sum(hrs) / 5) #annual hours

# Annual outage hours - zoom 
county_zoom_list = c("TX", "LA", "MS", "AL", "FL", "GA", "SC", "NC", "VA")
county_zoom = get_acs(geography = "county", state = county_zoom_list,
                      variables=c("B01003_001"), year = 2019, geometry = TRUE, 
                      cache_table = TRUE) %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME) %>%
  st_transform(5070) 
hr_zoom = county_zoom %>%
  left_join(hrs, by = c("GEOID")) %>%
  left_join(hurricane_group, by = c("GEOID" = "fips"))

gg3 = ggplot()+
  #geom_sf(data = hr_zoom, aes(fill = Hours), color = NA) +
  geom_sf(data = hr_zoom, aes(fill = sum_sust_hrs / 5), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey30") +
  geom_sf(data = county_zoom, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  #labs(title = "Annual Power Outages\nTropical Storms (2015 - 2019)", fill = "Hours") +
  #labs(title = "Annual Exposure to Sustained Winds \u2265 34 knots\nTropical Storms (2015 - 2019)", fill = "Hours") +
  #labs(title = expression(atop("Annual Power Outages", paste("Tropical Storms (2015 - 2019)"))), fill = "Hours") +
  labs(title = expression(atop("Annual Exposure to Sustained"~Winds>=34~"knots", paste("Tropical Storms (2015 - 2019)"))), fill = "Hours") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)


################################################################################
#### ML - EXPOSURE WINDS #######################################################
################################################################################
load(file = "Data/processed/sf_data_ALL_nototal.Rda")
sf_data_wind = sf_data_ALL %>%
  dplyr::filter(GEOID %in% as.vector(hurricane_group$fips)) %>% #filter to event of interest
  dplyr::select(-c(POPULATION, mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out, avalanche)) 
rm(list=c("sf_data_ALL")) 
gc() 

df_wind_PRED = sf_data_wind %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-duration_hr) %>%
  dplyr::select(-c(max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  dplyr::select(-c(all_of(our_events))) %>%
  group_by(GEOID) %>% 
  summarise(across(everything(), first))
df_wind_preds = df_wind_PRED %>% select(-GEOID)

df_wind_EXP = hurricane_group %>%
  filter(sum_sust_hrs > 0) %>%
  mutate(ln_sust = log(sum_sust_hrs / 5 + 0.001), ln_gust = log(sum_gust_hrs / 5 + 0.001)) %>%
  dplyr::select(c(fips, ln_sust, ln_gust))

df_data_WIND = df_wind_PRED %>%
  inner_join(df_wind_EXP, by = c("GEOID" = "fips")) %>%
  relocate(ln_sust, ln_gust)

set.seed(23)
df_split_wind = initial_split(df_data_WIND, prop = 0.80, strata = "ln_sust")
df_train_wind = training(df_split_wind)
df_test_wind = testing(df_split_wind)  
df_cv_wind = vfold_cv(df_train_wind, v = 10, repeats = 1)

recipe_sust = recipe(ln_sust ~ . , data = df_data_WIND) %>%
  step_rm(GEOID, ln_gust)
 
recipe_wind = recipe_sust
y_train_wind = df_train_wind %>% pull(ln_sust) %>% na.omit() 
y_test_wind = df_test_wind %>% pull(ln_sust) %>% na.omit() 
X_train_wind = df_train_wind %>% dplyr::select(-c(GEOID, ln_sust, ln_gust)) %>% na.omit()
X_test_wind = df_test_wind %>% dplyr::select(-c(GEOID, ln_sust, ln_gust)) %>% na.omit()

show_model_info("linear_reg")
wind_lre_model = linear_reg(penalty = tune(), mixture = tune()) %>% #lambda (penalty) and alpha/mixture (1 lasso, 0 ridge)
  set_engine("glmnet") %>%
  translate()
wind_lre_work = workflow() %>% 
  add_recipe(recipe_wind) %>%
  add_model(wind_lre_model)
set.seed(32); wind_lre_grid = dials::grid_max_entropy(parameters(penalty(), mixture()), size = 40)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
wind_lre_tune = wind_lre_work %>%
  tune_grid(resamples = df_cv_wind,
            grid = lre_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")
  ) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(wind_lre_tune, metric = "rmse")
wind_lre_tune_results = wind_lre_tune %>% collect_metrics()
wind_lre_best = wind_lre_tune %>% select_best(metric = "rmse")
wind_lre_fit = wind_lre_work %>%
  finalize_workflow(wind_lre_best) %>%
  last_fit(df_split_wind)
wind_lre_test = wind_lre_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
wind_lre_predictions = wind_lre_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)
rsq_lre_wind = paste(wind_lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
cverror_lre_wind = paste(show_best(wind_lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

show_model_info("boost_tree")
wind_xgb_model = boost_tree(mode = "regression", trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune(), mtry = tune()) %>%
  set_engine(engine = "xgboost") %>%
  translate()
wind_xgb_work = workflow() %>%
  add_recipe(recipe_wind) %>%
  add_model(wind_xgb_model)
set.seed(32); wind_xgb_grid = dials::grid_max_entropy(parameters(trees(), min_n(), tree_depth(), learn_rate(), loss_reduction(), finalize(mtry(), df_wind_preds)), size = 100)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
wind_xgb_tune = wind_xgb_work %>%
  tune_grid(resamples = df_cv_wind, 
            grid = wind_xgb_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")
  ) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(wind_xgb_tune, metric = "rmse")
wind_xgb_tune_results = wind_xgb_tune %>% collect_metrics()
wind_xgb_best = wind_xgb_tune %>% select_best(metric = "rmse")
wind_xgb_fit = wind_xgb_work %>%
  finalize_workflow(wind_xgb_best) %>%
  last_fit(df_split_wind)
wind_xgb_test = wind_xgb_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
wind_xgb_predictions = wind_xgb_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)
rsq_xgb_wind = paste(wind_xgb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
cverror_xgb_wind = paste(show_best(wind_xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

wind_bart_fit = bartMachineCV(data.frame(X_train_wind), y_train_wind, k_folds = 10, serialize = T) 
wind_bart_predictions = predict(wind_bart_fit, data.frame(X_test_wind))
rsq_bart_wind = format(round(1 - sum((y_test_wind - wind_bart_predictions)^2) / sum((y_test_wind - mean(y_test_wind))^2), 3), nsmall = 3)
bart_cv_wind = k_fold_cv(data.frame(X_train_wind), y_train_wind, k_folds = 10)
cverror_bart_wind = format(round(bart_cv_wind$rmse, 3), nsmall = 3)

gg_wind = dplyr::tibble(actual = y_test_wind, 
                   eNet = as.vector(wind_lre_predictions$.pred),
                   bart = as.vector(wind_bart_predictions),
                   #rf = as.vector(rf_predictions$.pred),
                   xgb = as.vector(wind_xgb_predictions$.pred)
)
gg_wind = arrange(gg_wind, actual)
gg_wind$index = seq.int(nrow(gg_wind))
gg_wind_actual = gg_wind %>% dplyr::select(index, actual)
gg_wind_pred = gg_wind %>% dplyr::select(-actual) %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
gg_wind_all = gg_wind %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")

color_vec = c("black", "#FDE725FF", "#35B779FF", "#440154FF")
lty_vec = c(1, 1, 1, 1)
alpha_vec = c(1, 0.6, 0.6, 0.6)

plot_filtering_estimates2 <- function(df) {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg_wind$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
    geom_line(data = gg_wind_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
    scale_color_manual(
      values = color_vec, 
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart_wind) * ")"),
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre_wind) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb_wind) * ")")
      ),
      name = element_blank()) +
    scale_linetype_manual(
      values = lty_vec,
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart_wind) * ")"),
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre_wind) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb_wind) * ")")         
      ),
      name = element_blank()) + 
    scale_alpha_manual(
      values = alpha_vec,
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart_wind) * ")"),
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre_wind) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb_wind) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb_wind) * ")")
      ),
      name = element_blank()) + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("County") +
    ylab(expression("Annual Sustained"~Winds>=34~"knots (ln - hrs)")) +
    #ylab("Annual Sustained Winds >= 34 knots (ln - hrs)") + 
    ggtitle("Exposure to Tropical Storms: Test Sample") + 
    # guides(
    #   color = guide_legend(order = 2),
    #   shape = guide_legend(order = 1),
    #   linetype = guide_legend(order = 2)
    # ) + 
    theme(legend.spacing.y = unit(-0.25, "cm"),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = c(.225, .85),
          plot.title = element_text(hjust = 0.5)
    )
  
  print(p)
}
plot_filtering_estimates2(gg)

wind_final_model = wind_xgb_work %>%
  finalize_workflow(wind_xgb_best) %>%
  fit(df_data_WIND)
wind_final_obj = extract_fit_parsnip(wind_final_model)$fit
wind_importance = xgb.importance(model = wind_final_obj)
head(wind_importance, n = 20) 
vip(wind_final_obj, n = 20)


