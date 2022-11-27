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

my_event = "hurricanes"

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
  dplyr::filter(eval(parse(text = my_event)) >= 1) %>% #filter to event of interest
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
set.seed(32)
df_split = initial_split(df_data, prop = 0.80, strata = eval(my_event))
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

# Recipes for tidymodels 
if (my_event == "hurricanes") {
  recipe_hurricanes = recipe(hurricanes ~ . , data = df_data) %>% 
    step_rm(GEOID, droughts, extreme_cold, extreme_heat, floods, high_winds, wildfires, winter_storms) %>% 
    step_naomit()
  recipe_mine = recipe_hurricanes
} else if (my_event == "winter_storms") {
  recipe_winter = recipe(winter_storms ~ . , data = df_data) %>% 
    step_rm(GEOID, droughts, extreme_cold, extreme_heat, floods, high_winds, wildfires, hurricanes) %>% 
    step_naomit()
  recipe_mine = recipe_winter
} else if (my_event == "high_winds") {
  recipe_highwind = recipe(high_winds ~ . , data = df_data) %>% 
    step_rm(GEOID, droughts, extreme_cold, extreme_heat, floods, winter_storms, wildfires, hurricanes) %>% 
    step_naomit()
  recipe_mine = recipe_highwind
} else if (my_event == "droughts") {
  recipe_droughts = recipe(droughts ~ . , data = df_data) %>% 
    step_rm(GEOID, hurricanes, extreme_cold, extreme_heat, floods, high_winds, wildfires, winter_storms) %>% 
    step_naomit()
  recipe_mine = recipe_droughts
} else if (my_event == "floods") {
  recipe_floods = recipe(floods ~ . , data = df_data) %>% 
    step_rm(GEOID, droughts, extreme_cold, extreme_heat, hurricanes, high_winds, wildfires, winter_storms) %>% 
    step_naomit()
  recipe_mine = recipe_floods
}

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
y_train = df_train %>% pull(eval(parse(text = my_event))) %>% na.omit() 
y_test = df_test %>% pull(eval(parse(text = my_event))) %>% na.omit() 
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
    ylab("Annual Outage Hours (ln)") + 
    ggtitle(paste("Power Outages Due to ", str_to_title(gsub("_", " ", my_event)), ": Test Sample", sep = "")) +
    #ggtitle("Power Outages Due to Tropical Storms: Test Sample") + 
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
# Base map
load(file = "Data/processed/county_map_proj.Rda")

# Annual outage hours - all US
hrs = df_data %>%
  dplyr::select(GEOID, eval(my_event)) %>%
  mutate(hrs = exp(eval(parse(text = my_event)))) %>%
  group_by(GEOID) %>%
  summarise(Hours = sum(hrs) / 5) #annual hours

county_hrs = county_map_proj %>%
  left_join(hrs, by = c("GEOID"))

# # Annual outage hours - zoom 
county_zoom_list = c("TX", "LA", "MS", "AL", "FL", "GA", "SC", "NC", "VA")
county_zoom = get_acs(geography = "county", state = county_zoom_list,
                      variables=c("B01003_001"), year = 2019, geometry = TRUE,
                      cache_table = TRUE) %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME) %>%
  st_transform(5070)
# hr_zoom = county_zoom %>%
#   left_join(hrs, by = c("GEOID")) %>%
#   left_join(hurricane_group, by = c("GEOID" = "fips"))
geoid_zoom = county_zoom$GEOID
county_hrs_zoom = county_hrs %>%
  dplyr::filter(GEOID %in% geoid_zoom)

county_plot = county_hrs
gg3 = ggplot()+
  geom_sf(data = county_plot, aes(fill = Hours), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey30") +
  geom_sf(data = county_plot, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  labs(title = paste("Annual Power Outages Due to ", str_to_title(gsub("_", " ", my_event)), " (2015 - 2019)", sep = ""), fill = "Hours") +
  #labs(title = expression(atop("Annual Power Outages", paste("Tropical Storms (2015 - 2019)"))), fill = "Hours") +
  #labs(title = expression(atop("Annual Exposure to Sustained"~Winds>=34~"knots", paste("Tropical Storms (2015 - 2019)"))), fill = "Hours") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)




