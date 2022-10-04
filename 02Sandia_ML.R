#### Machine learning 
library(tidyverse)
library(tidymodels)
library(sf)
library(parallel)
library(doParallel)

################################################################################
#### PRE-PROCESSING ############################################################
################################################################################
# Parallel processing setup
num_cores = detectCores() - 1
unregister_dopar = function() { #function to un-register parallel processing in doParallel
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

# Load data and select final DVs
load(file = "Data/processed/sf_data_CLEAN.Rda")
df_data_CLEAN = sf_data_CLEAN %>%
  st_set_geometry(NULL) %>%
  mutate(ln_hrs = log(duration_hr)) %>% 
  mutate(ln_cust = log(max_cust_out)) %>% # take log of DVs# take log of DVs
  mutate(pct_cust = max_frac_cust_out) %>%
  dplyr::select(-c(POPULATION, mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out, duration_hr)) %>%
  relocate(c(ln_hrs, ln_cust, pct_cust))

# Split into training vs testing
set.seed(23)
df_split = initial_split(df_data_CLEAN, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)
df_preds = df_data_CLEAN %>% dplyr::select(-c(ln_hrs, ln_cust, pct_cust, GEOID))

# Recipes for tidymodels 
recipe_hrs = recipe(ln_hrs ~ . , data = df_data_CLEAN) %>% step_rm(ln_cust, pct_cust, GEOID)
recipe_cust = recipe(ln_cust ~ . , data = df_data_CLEAN) %>% step_rm(ln_hrs, pct_cust, GEOID)
recipe_pct = recipe(pct_cust ~ . , data = df_data_CLEAN) %>% step_rm(ln_hrs, ln_cust, GEOID) %>% step_naomit(pct_cust)

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
### Define which recipe you want to use 
recipe_mine = recipe_hrs

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

### Random Forest ##############################################################
show_model_info("rand_forest")
rf_model = rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression") %>%
  translate()
rf_work = workflow() %>%
  add_recipe(recipe_mine) %>%
  add_model(rf_model)
#rf_grid = expand.grid(mtry = c(5, 10, 20), trees = c(100, 500, 1000), min_n = c(3, 5, 10))
set.seed(32); rf_grid = dials::grid_max_entropy(parameters(min_n(), trees(), finalize(mtry(), df_preds)), size = 60)
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

