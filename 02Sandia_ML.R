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
load(file = "sf_data_CLEAN.Rda")
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
recipe_mine = recipe_pct

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




# ## Testing results
# #y_test = df_hours %>% slice(-df_split$in_id) %>% dplyr::select(ln_hrs) %>% pull()
# y_test = df_test$ln_hrs
# 
# rsq_lre = paste(lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
# rsq_xgb = paste(xgb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
# rsq_rf = paste(rf_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
# 
# cverror_lre = paste(show_best(lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
# cverror_xgb = paste(show_best(xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
# cverror_rf = paste(show_best(rf_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
# 
# ### PLOTTING 
# gg = dplyr::tibble(actual = df_test$ln_hrs,
#                    eNet = as.vector(lre_predictions$.pred),
#                    rf = as.vector(rf_predictions$.pred),
#                    xgb = as.vector(xgb_predictions$.pred)
# )
# gg = arrange(gg, actual)
# gg$index = seq.int(nrow(gg))
# #gg = gg %>% filter(eNet < 20)
# #gg$mean = mean(gg$actual, na.rm = T)
# gg_actual = gg %>% dplyr::select(index, actual)
# gg_pred = gg %>% dplyr::select(-actual) %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
# gg_all = gg %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
# color_vec= c("black", "#FDE725FF", "#35B779FF")
# color_vec = c("black", "#FDE725FF", "#35B779FF", "#440154FF")
# lty_vec = c(1, 1, 1)
# lty_vec = c(1, 1, 1, 1)
# alpha_vec = c(1, 0.7, 0.7)
# alpha_vec = c(1, 0.7, 0.7, 0.7)
# 
# plot_filtering_estimates2 <- function(df) {
#   p = ggplot() + 
#     theme_classic() + 
#     geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
#     geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
#     scale_color_manual(
#       values = color_vec, 
#       labels = c("Actual",
#                  bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
#                  bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
#                  bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
#       ),
#       name = element_blank()) +
#     scale_linetype_manual(
#       values = lty_vec,
#       labels = c("Actual",
#                  bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
#                  bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
#                  bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")         
#       ),
#       name = element_blank()) + 
#     scale_alpha_manual(
#       values = alpha_vec,
#       labels = c("Actual",
#                  bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
#                  bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
#                  bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
#       ),
#       name = element_blank()) + 
#     ylab("Outage Duration (log-hrs)") + 
#     scale_y_continuous(labels = function(x) paste0(x), limits = c(2, 10)) +
#     xlab("Outage Index (event x county)") +
#     #ggtitle("County-Level Predictions: Test Sample\nNo 'Total' Variables") + 
#     ggtitle("County-Level Predictions: Test Sample") + 
#     # guides(
#     #   color = guide_legend(order = 2),
#     #   shape = guide_legend(order = 1),
#     #   linetype = guide_legend(order = 2)
#     # ) + 
#     theme(legend.spacing.y = unit(-0.25, "cm"),
#           legend.direction = "vertical",
#           legend.box = "vertical",
#           legend.position = c(.225, .8),
#           plot.title = element_text(hjust = 0.5)
#     )
#   
#   print(p)
# }
# plot_filtering_estimates2(gg)

