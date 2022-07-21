#### LOAD PACKAGES and SET DIRECTORY 
options(java.parameters = "-Xmx7g")
library(bartMachine)
set_bart_machine_num_cores(4)
library(tidyverse)
library(tidymodels)
library(tidycensus) 
library(sf)
library(lme4)
library(corrplot)
library(viridis)
library(doParallel)
#library(sqldf)

tidymodels_prefer()
setwd("~/Documents/01_VECTOR.nosync/Sandia")
num_cores = detectCores() - 1

#function to un-register parallel processing in doParallel
unregister_dopar = function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#### LOAD VARIABLES ##############################################################################
load("./Data/SandiaVariables.Rda") #from Sandia1_data.R
load(file = "Data/outages.Rda")
outages_csv$fips_code = as.character(outages_csv$fips_code)
fips = outages_csv$fips_code
fips0 = str_pad(fips, 5, pad = "0")
outages_csv$fips_code = fips0
county_map_outages = outages_csv %>%
  mutate(rowID = row_number()) %>% # add row number to index events 
  mutate(date_day = str_extract(date_hour, "^.{10}")) %>% #returns first 10 characters of date_time (yyyy-mm-dd)
  mutate(date_month = str_extract(date_hour, "^.{7}")) # returns yyyy-mm
county_map_outages_FILTER = county_map_outages %>%
  filter(outage_status %in% c("pre", "start", "during", "end")) #filter events flagged as outages
county_map_outages_JOIN = county_map_area %>%
  inner_join(county_map_outages_FILTER, by = c("GEOID" = "fips_code")) %>% # join outage/hourly weather data 
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) %>% #join socio-economic variables  
  left_join(county_map_soil_CLEAN, by = c("GEOID", "date_hour")) %>% #join soil moisture by county and hour time-stamp 
  left_join(county_map_spi_CLEAN, by = c("GEOID", "date_month")) #join SPI by GEOID and month 
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
county_outages_GROUP = outages_group %>%
  left_join(county_map_static_CLEAN, by = c("GEOID")) %>% #join static environmental variables
  left_join(county_scores_CLEAN, by = c("GEOID")) #join socio-economic variables 

##########################################################################################################
#### MACHINE LEARNING ####################################################################################
##########################################################################################################
df_data = data.frame(county_outages_GROUP) %>%
  dplyr::filter(out_hrs > quantile(county_outages_GROUP$out_hrs, .9)) %>% #filter to big events 
  dplyr::select(-outage_number, -GEOID, -out_percust) %>% #get rid of variables not used 
  mutate(out_hrs = log(out_hrs), out_maxcust = log(out_maxcust)) %>% #take log of DVs
  rename(ln_hrs = out_hrs, ln_cust = out_maxcust) #rename to emphasize ln 
rm(list=setdiff(ls(), "df_data")) 
gc() 

## Split into training vs testing
set.seed(32)
df_split = initial_split(df_data, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

## Pre-processing (recipe)
hours_recipe = recipe(ln_hrs ~ . , data = df_data) %>%
  step_rm(ln_cust) %>% 
  step_impute_knn(all_predictors()) %>% #knn impute missing predictors (if any)
  # step_normalize(all_predictors()) %>% #z-score standardize all predictors (important for PLS or NN)
  step_zv(all_predictors()) %>% #removes predictors of single value 
  step_corr(all_predictors())  #removes highly correlated 
prep(hours_recipe) #shows changes 
hours_juice = prep(hours_recipe) %>% juice() #view prepared dataset 

## BART (not part of tidymodels yet) 
df_bart = prep(hours_recipe) %>% juice()
df_bart_train = df_bart %>% slice(df_split$in_id) %>% dplyr::select(-ln_hrs)
X = data.frame(df_bart_train)
y = df_bart %>% dplyr::select(ln_hrs) %>% slice(df_split$in_id) %>% pull()
bart_fit = bartMachineCV(X, y, k_folds = 10) #bartMachine CV win: k: 2, nu: 3, q: 0.99, num_trees: 50 

df_bart_test = df_bart %>% slice(-df_split$in_id) %>% dplyr::select(-ln_hrs)
X_test = data.frame(df_bart_test)
y_test = df_bart %>% slice(-df_split$in_id) %>% dplyr::select(ln_hrs) %>% pull()
bart_pred = predict(bart_fit, X_test)
bart_rmse = sqrt(mean((bart_pred - y_test)^2)) 
bart_rsq = 1 - sum((y_test - bart_pred)^2) / sum((y_test - mean(y_test))^2) 
bart_CI = round(calc_credible_intervals(bart_fit, X_test, ci_conf = 0.95), 2)

## Random Forest
#https://www.rebeccabarter.com/blog/2020-03-25_machine_learning/
show_model_info("rand_forest")
rf_model = rand_forest(mtry = tune(), trees = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "permutation") %>%
  set_mode("regression") %>%
  translate()
rf_work = workflow() %>%
  add_recipe(hours_recipe) %>%
  add_model(rf_model)
rf_grid = expand.grid(mtry = c(1, 3, 5), trees = c(100, 500, 1000), min_n = c(3, 5, 10))
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
rf_tune = rf_work %>%
  tune_grid(resamples = df_cv,
            grid = rf_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")) #parallel processing turns off verbose
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

## Lasso/Ridge/ElasticNet 
#https://dnield.com/posts/tidymodels-intro/ 
show_model_info("linear_reg")
lre_model = linear_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  translate()
lre_work = workflow() %>% 
  add_recipe(hours_recipe) %>%
  add_model(lre_model)
lre_grid = grid_regular(parameters(penalty(), mixture()), levels = c(5, 5))
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
lre_tune = lre_work %>%
  tune_grid(resamples = df_cv,
            grid = lre_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")) #parallel processing turns off verbose
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

## GBM 
#https://www.r-bloggers.com/2020/05/using-xgboost-with-tidymodels/
# https://datascience.stackexchange.com/questions/15305/how-does-xgboost-learn-what-are-the-inputs-for-missing-values
# https://www.capitalone.com/tech/machine-learning/how-to-control-your-xgboost-model/
show_model_info("boost_tree")
gb_model = boost_tree(mode = "regression", trees = 1000, 
                      min_n = tune(), tree_depth = tune(), learn_rate = tune(), loss_reduction = tune()) %>%
  set_engine(engine = "xgboost") %>%
  translate()
gb_work = workflow() %>%
  add_recipe(hours_recipe) %>%
  add_model(gb_model)
gb_grid = dials::grid_max_entropy(parameters(min_n(), tree_depth(), learn_rate(), loss_reduction()), size = 100)
cl = makeCluster(num_cores, type = "FORK")
registerDoParallel(cl, cores = num_cores)
gb_tune = gb_work %>%
  tune_grid(resamples = df_cv,
            grid = gb_grid,
            metrics = metric_set(yardstick::rmse, yardstick::rsq),
            control = tune::control_grid(verbose = T, allow_par = T, parallel_over = "resamples")) #parallel processing turns off verbose
stopCluster(cl) 
unregister_dopar()
show_best(gb_tune, metric = "rmse")
gb_tune_results = gb_tune %>% collect_metrics()
gb_best = gb_tune %>% select_best(metric = "rmse")
gb_fit = gb_work %>%
  finalize_workflow(gb_best) %>%
  last_fit(df_split)
gb_test = gb_fit %>% collect_metrics() #metrics evaluated on test sample (b/c last_fit() function) 
gb_predictions = gb_fit %>% collect_predictions() #predictions for test sample (b/c last_fit() function)


##########################################################################################################
###### PLOTTING ##########################################################################################
##########################################################################################################
gg = dplyr::tibble(actual = y_test,
                   lasso = as.vector(lre_predictions$.pred),
                   rf = as.vector(rf_predictions$.pred),
                   xgb = as.vector(gb_predictions$.pred),
                   bart = bart_pred
  )
gg = arrange(gg, actual)
gg$index = seq.int(nrow(gg))
#gg$mean = mean(gg$actual, na.rm = T)
gg_actual = gg %>% dplyr::select(index, actual)
gg_pred = gg %>% dplyr::select(-actual) %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")
gg_all = gg %>% pivot_longer(!index, names_to = "Model", values_to = "ypred")

#Test sample r-squared
rsq_lre = paste(lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_rf = paste(rf_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_gbm = paste(gb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_bart = paste(round(bart_rsq, 3) %>% format(nsmall = 3))

#Cross-validation errors
cverror_lre = paste(show_best(lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_rf = paste(show_best(rf_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_gbm = paste(show_best(gb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_bart = paste(data.frame(bart_fit$cv_stats) %>% dplyr::slice(1) %>% pull(oos_error) %>% round(3) %>% format(nsmall = 3))

#color_vec = c("black", "#bef7ff", "#86b7ff", "#4d76ff", "#1536ff")
#color_vec = c("black", "#FC4E07", "#1E88E5", "forestgreen", "#E7B800")
color_vec = c("black", "#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")

# ggplot (for all)
plot_filtering_estimates2 <- function(df) {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50") +
    # geom_point(data = gg_all, aes(x = index, y = ypred, color = Model, shape = Model), alpha = 1) +
    # scale_shape_manual(
    #   name = element_blank(),
    #   values = c(16, 2, 3, 4, 5),
    #   labels = c("Actual",
    #              bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"), 
    #              bquote("Lasso (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
    #              bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
    #              bquote("XGB (" * R^2 ~ "=" ~ .(rsq_gbm) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_gbm) * ")")
    #   )
    # ) +
    scale_color_manual(
      values = color_vec, 
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"), 
                 bquote("Lasso (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_gbm) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_gbm) * ")")
      ),
      name = element_blank()) +
    geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model), alpha = 0.8) +
    #scale_color_viridis(discrete = T, option = "D", direction = 1) + 
    # scale_color_manual(values = c("#FC4E07", "#1E88E5", "forestgreen", "#E7B800"), 
    #                    labels = c(bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"), 
    #                               bquote("Lasso (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
    #                               bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
    #                               bquote("XGB (" * R^2 ~ "=" ~ .(rsq_gbm) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_gbm) * ")")
    #                    ),
    #                    name = element_blank()) +
    scale_linetype_manual(
      values = c(1, 1, 1, 1, 1),
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"), 
                 bquote("Lasso (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_gbm) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_gbm) * ")")         
                 ),
      name = element_blank()) + 
    #scale_color_locuszoom() + 
    ylab("Outage Duration (log hours)") + 
    #ylab("Max Cust. Outages (log)") + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Outage Index (event x county)") +
    ggtitle("County-Level Predictions: Test Sample") + 
    # guides(
    #   color = guide_legend(order = 2),
    #   shape = guide_legend(order = 1),
    #   linetype = guide_legend(order = 2)
    # ) + 
    theme(legend.spacing.y = unit(-0.25, "cm"),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.position = c(.215, .8),
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
library(xgboost)
library(vip)
library(pdp)
final_model = gb_work %>%
  finalize_workflow(gb_best) %>%
  fit(df_data)
final_obj = extract_fit_parsnip(final_model)$fit
importance = xgb.importance(model = final_obj)
head(importance, n = 15); vip(final_obj, n = 20)
#xgb.plot.importance(importance_matrix = importance)
x = hours_juice %>% dplyr::select(-ln_hrs)
p_ps = pdp::partial(final_obj, pred.var = "PS_sd", ice = T, center = F,
                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                     train = x, type = "regression")
p_slp = pdp::partial(final_obj, pred.var = "SLP_mean", ice = T, center = F,
                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                     train = x, type = "regression")
p_v = pdp::partial(final_obj, pred.var = "V_max", ice = T, center = F,
                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                     train = x, type = "regression")
p_tqv = pdp::partial(final_obj, pred.var = "TQV_mean", ice = T, center = F,
                   plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                   train = x, type = "regression")
grid.arrange(p_ps, p_slp, p_v, p_tqv, ncol = 2)

p_wind =  pdp::partial(final_obj, pred.var = "WIND_max", ice = T, center = F,
                       plot = T, rug= T, alpha = 0.1, plot.engine = "ggplot2",
                       train = x, type = "regression") +
  geom_vline(xintercept = 8.9, color = "blue") + 
  ggtitle("PDP and ICE - Max Wind Speed (m/s)") + 
  annotate("text", x =12.5, y = 6, label="Tropical Storm Wind Forecast", color="blue", size = 3) 

p1 = pdp::partial(final_obj, pred.var = c("PS_sd", "SLP_mean"), 
                  plot = T, chull= T, rug = T, plot.engine = "ggplot2",
                  train = x, type = "regression")
p2 = pdp::partial(final_obj, pred.var = c("PS_sd", "V_max"), 
                  plot = T, chull= T, plot.engine = "ggplot2",
                  train = x, type = "regression")
p3 = pdp::partial(final_obj, pred.var = c("SLP_mean", "V_max"), 
                  plot = T, chull= T, plot.engine = "ggplot2",
                  train = x, type = "regression")

### View tree diagram 
xgb.plot.tree(model = final_obj, trees = 1:3) 
  