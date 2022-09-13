#### Machine learning 
library(tidyverse)
library(tidymodels)
library(sf)
library(parallel)
library(doParallel)
library(xgboost)
library(vip)
library(spdep)
library(pdp)

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
load(file = "Data/processed/sf_data_CLEAN.Rda")
df_data_EVENTS = sf_data_CLEAN %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-duration_hr) %>%
  dplyr::select(-c(mean_cust_out:max_frac_cust_out, max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  group_by(GEOID) %>% 
  summarise(across(all_of(our_events), ~ sum(.x, na.rm = TRUE)))

df_data_OTHER = sf_data_CLEAN %>%
  st_set_geometry(NULL) %>%
  dplyr::select(-duration_hr) %>%
  dplyr::select(-c(mean_cust_out:max_frac_cust_out, max_WIND10M:spi24)) %>% #strip out all dynamic factors 
  dplyr::select(-c(all_of(our_events), "avalanche")) %>%
  group_by(GEOID) %>% 
  summarise(across(everything(), first))

df_data_CLEAN = df_data_EVENTS %>%
  inner_join(df_data_OTHER, by = "GEOID")

# Split into training vs testing
set.seed(23)
df_split = initial_split(df_data_CLEAN, prop = 0.80, strata = "winter_storms")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)
df_preds = df_data_CLEAN %>% dplyr::select(-c(GEOID, POPULATION, all_of(our_events)))

# Recipes for tidymodels 
recipe_winter = recipe(winter_storms ~ . , data = df_data_CLEAN) %>% 
  step_rm(GEOID, POPULATION, droughts, extreme_cold, extreme_heat, floods, hurricanes, wildfires, high_winds) %>% 
  step_naomit(winter_storms)

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
### Define which recipe you want to use 
recipe_mine = recipe_winter
y_test = df_test$winter_storms

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
final_lre = lre_work %>%
  finalize_workflow(lre_best) %>%
  fit(df_data_CLEAN)
print(tidy(final_lre), n = 150)

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
gg = dplyr::tibble(actual = y_test / 5, #divide by 5 to annualize
                   eNet = as.vector(lre_predictions$.pred) / 5,
                   #rf = as.vector(rf_predictions$.pred),
                   xgb = as.vector(xgb_predictions$.pred) / 5
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
    ylab("Severe Outages due to WinterStorms (per year)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
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
#https://christophm.github.io/interpretable-ml-book/ice.html
#https://datascience.stackexchange.com/questions/15305/how-does-xgboost-learn-what-are-the-inputs-for-missing-values
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
pdp1 = pdp::partial(final_obj, pred.var = "Forest", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp2 = pdp::partial(final_obj, pred.var = "DEM_max", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp3= pdp::partial(final_obj, pred.var = "DEM_sd", ice = F, center = F,
                   plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                   train = X_train, type = "regression")
pdp4 = pdp::partial(final_obj, pred.var = "QAGE17", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp5 = pdp::partial(final_obj, pred.var = "Wetlands", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp6 = pdp::partial(final_obj, pred.var = "DEM_min", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp7 = pdp::partial(final_obj, pred.var = "HDENS", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp8 = pdp::partial(final_obj, pred.var = "RZ_med", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")

grid.arrange(pdp1, pdp2, pdp3, pdp4, ncol = 2)
grid.arrange(pdp5, pdp6, pdp7, pdp8, ncol = 2)

