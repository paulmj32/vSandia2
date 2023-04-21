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
library(grid)
library(gridExtra)
library(cowplot)

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
setwd("/Users/paul/vSandia2")
load(file = "Data/processed/sf_data_ALL_nototal.Rda")

sf_data = sf_data_ALL %>%
  dplyr::filter(hurricanes >= 1) %>% #filter to event of interest
  mutate(ln_hrs = log(duration_hr)) %>% 
  mutate(ln_cust = log(max_cust_out)) %>% 
  #mutate(pct_cust = max_frac_cust_out) %>%
  mutate(pct_cust = log(max_frac_cust_out)) %>%
  dplyr::select(-c(POPULATION, mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out, 
                   duration_hr, all_of(our_events), avalanche)) %>%
  relocate(c(ln_hrs, ln_cust, pct_cust))

rm(list=c("sf_data_ALL")) 
gc() 

df_data = sf_data %>%
  st_set_geometry(NULL)

# #Demo for hiba's class
# set.seed(32)
# df_hiba = initial_split(df_data, prop = 0.4, strata = "ln_hrs")
# df_hiba2 = training(df_hiba)
# df_hurricanes = df_hiba2 %>%
#   dplyr::select(ln_hrs, max_WIND10M:max_T10M, delta_pressure, spi03, spi12, Barren:RZ_mode
#                 ,QPOVTY, QMINORITY, QFEMALE, QED12, QSSBEN
#                 )
# save(df_hurricanes, file = "Data/df_hurricanes.Rda")  
  
# Split into training vs testing
set.seed(23)
df_split = initial_split(df_data, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)
df_preds = df_data %>% dplyr::select(-c(ln_hrs, ln_cust, pct_cust, GEOID))

# Recipes for tidymodels 
recipe_hrs = recipe(ln_hrs ~ . , data = df_data) %>% step_rm(ln_cust, pct_cust, GEOID) %>% step_naomit(ln_hrs)
recipe_pct = recipe(pct_cust ~ . , data = df_data) %>% step_rm(ln_hrs, ln_cust, GEOID) %>% step_naomit(pct_cust) 

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
### Define which recipe and responses you want to use 
recipe_mine = recipe_hrs

## Recipe - ln_hrs
# y_train = df_train %>% pull(ln_hrs) %>% na.omit()
# y_test = df_test %>% pull(ln_hrs) %>% na.omit()
# X_train = df_train %>% dplyr::select(-c(GEOID, ln_hrs, ln_cust, pct_cust)) %>% na.omit()
# X_test = df_test %>% dplyr::select(-c(GEOID, ln_hrs, ln_cust, pct_cust)) %>% na.omit()

## Recipe - pct_cust
y_train = df_train %>% na.omit() %>% pull(pct_cust)
y_test = df_test %>% na.omit() %>% pull(pct_cust)
X_train = df_train %>% na.omit() %>% dplyr::select(-c(GEOID, ln_hrs, ln_cust, pct_cust))
X_test = df_test %>% na.omit() %>% dplyr::select(-c(GEOID, ln_hrs, ln_cust, pct_cust))
na_test = which(is.na(df_test$pct_cust))
  
### Lasso, Ridge Regression, and Elastic Net ###################################
#https://www.tidyverse.org/blog/2020/11/tune-parallel/
show_model_info("linear_reg")
lre_model = linear_reg(penalty = tune(), mixture = tune()) %>% #lambda (penalty) and alpha/mixture (1 lasso, 0 ridge)
  set_engine("glmnet") %>%
  translate()
lre_work = workflow() %>% 
  add_recipe(recipe_hrs_lre) %>%
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

#regular MLR
df_mlr =  X_train %>%
  mutate(y = y_train)
mlr = lm(y ~ ., data = df_mlr)
summary(mlr)
mlr_predictions = predict(mlr, newdata = X_test)
rsq_mlr = format(round(cor(mlr_predictions, y_test)^2, 3), nsmall = 3)
asd = vip(mlr, n = 20)
#imp.new = asd$data$Importance/sum(asd$data$Importance)
#asd$data$Importance = imp.new
plot(asd)
asdf = data.frame(mlr$coefficients) 


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
#save(bart_fit, file = "bart_fit.Rda")

rsq_bart = format(round(1 - sum((y_test - bart_predictions)^2) / sum((y_test - mean(y_test))^2), 3), nsmall = 3)
bart_cv = k_fold_cv(data.frame(X_train), y_train, k_folds = 10)
cverror_bart = format(round(bart_cv$rmse, 3), nsmall = 3)

### SPATIAL REGRESSION #########################################################
#https://rpubs.com/quarcs-lab/tutorial-spatial-regression
# sf_train_sreg = county_map_proj %>%
#   inner_join(df_data, by = "GEOID") %>%
#   dplyr::select(-c(GEOID, NAME, POPULATION, ln_cust, pct_cust))
# neighbors = poly2nb(sf_train_sreg, queen = T) #find county neighbors 
# weights = nb2listw(include.self(neighbors)) #include self for weights
# sdm.eq = as.formula(paste("ln_hrs", 
#                           paste(colnames(sf_train_sreg)[2:length(colnames(sf_train_sreg))-1], collapse = " + "), 
#                           sep = " ~ "))
# asd = spatialreg::lagsarlm(formula = sdm.eq, data = sf_train_sreg, listw = weights, type = "mixed")

##########################################################################################################
##### PLOTTING ###########################################################################################
##########################################################################################################
if (length(na_test) > 0) {
  gg = dplyr::tibble(actual = y_test, 
                     #eNet = as.vector(lre_predictions$.pred)[-c(na_test)],
                     bart = as.vector(bart_predictions),
                     mlr = as.vector(mlr_predictions),
                     #rf = as.vector(rf_predictions$.pred),
                     xgb = as.vector(xgb_predictions$.pred)[-c(na_test)]
  )
} else {
  gg = dplyr::tibble(actual = y_test, 
                     #eNet = as.vector(lre_predictions$.pred),
                     bart = as.vector(bart_predictions),
                     mlr = as.vector(mlr_predictions),
                     #rf = as.vector(rf_predictions$.pred),
                     xgb = as.vector(xgb_predictions$.pred)
  )
}
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
                 bquote("MLR (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
      ),
      name = element_blank()) +
    scale_linetype_manual(
      values = lty_vec,
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"),
                 bquote("MLR (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")         
      ),
      name = element_blank()) + 
    scale_alpha_manual(
      values = alpha_vec,
      labels = c("Actual",
                 bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_bart) * ")"),
                 bquote("MLR (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
      ),
      name = element_blank()) + 
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("Index (County x Event)") +
    ylab("Hours (ln)") + 
    ggtitle("WinterStorm Outage Duration: Test Sample") + 
    #ylab("Pct. Max Customers Out (ln)") + 
    #ggtitle("Hurricane Outage Impact: Test Sample") + 
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
  fit(df_data)
print(tidy(final_lre), n = 150)
final_lre_model = extract_fit_parsnip(final_lre)$fit
lre_coef = predict(final_lre_model, type = "coefficients", s = lre_best$penalty)
df_lre_coef = data.frame(var = lre_coef@Dimnames[[1]]) %>%
  mutate(s1 = 0)
df_lre_coef[lre_coef@i + 1, 2] = lre_coef@x 
df_lre_coef2 = df_lre_coef %>%
  mutate(rel.s1 = abs(s1) / sum(abs(df_lre_coef$s1)))
vip(final_lre_model, n = 20)

final_model = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_data)
final_obj = extract_fit_parsnip(final_model)$fit
importance = xgb.importance(model = final_obj)
head(importance, n = 20) 
vip(final_obj, n = 20)

bart_vimp = investigate_var_importance(bart_fit)


df_imp = importance %>% 
  mutate(cat = case_when(
    Feature %in% c("delta_WIND10M", "delta_vapor", "delta_T10M", "delta_pressure", 
                   "max_humidity_10M", "max_T10M", "max_WIND10M",
                   "spi03", "spi06", "spi12", "spi24") ~ "Weather",
    Feature %in% c("Barren", "Cultivated", "Forest", "Herbaceous", "Shrub", "Water", "Wetlands",
                   "DEM_sd", "DEM_min", "RZ_med", "RZ_mode") ~ "Environment",
    TRUE ~ "Socio-Economic"
  )) %>%
  mutate(Feature = case_when(Feature == "max_humidity_10M" ~ "max_humidity10M", TRUE ~ as.character(Feature))) %>% #name correction
  dplyr::arrange(desc(Gain)) %>%
  dplyr::slice(1:20) %>%
  mutate(Feature1 = str_replace(Feature, "_", " ")) %>%
  mutate(Feature2 = tolower(Feature1)) %>%
  mutate(Feature_clean = str_to_title(Feature2)) %>%
  dplyr::select(-c(Feature1, Feature2))

df_imp_mlr = asd$data %>%
  mutate(Feature = Variable) %>%
  mutate(cat = case_when(
    Feature %in% c("delta_WIND10M", "delta_vapor", "delta_T10M", "delta_pressure", 
                   "max_humidity_10M", "max_T10M", "max_WIND10M",
                   "spi03", "spi06", "spi12", "spi24") ~ "Weather",
    Feature %in% c("Barren", "Cultivated", "Forest", "Herbaceous", "Shrub", "Water", "Wetlands",
                   "DEM_sd", "DEM_min", "RZ_med", "RZ_mode") ~ "Environment",
    TRUE ~ "Socio-Economic"
  )) %>%
  mutate(Feature = case_when(Feature == "max_humidity_10M" ~ "max_humidity10M", TRUE ~ as.character(Feature))) %>% #name correction
  dplyr::arrange(desc(Importance)) %>%
  dplyr::slice(1:20) %>%
  mutate(Feature1 = str_replace(Feature, "_", " ")) %>%
  mutate(Feature2 = tolower(Feature1)) %>%
  mutate(Feature_clean = str_to_title(Feature2)) %>%
  dplyr::select(-c(Feature1, Feature2))

fill_vec = c("#CEB966", "#A379BB", "#6BB1C9")  
pp2 = ggplot() + 
  theme_classic() + 
  geom_col(data = df_imp, aes(x = Gain, y = fct_reorder(Feature_clean, Gain), fill = cat, color = cat), 
          alpha = .75) +
  scale_fill_manual(values = fill_vec) +
  scale_color_manual(values = fill_vec) +
  xlab("Importance") +
  ylab(element_blank()) + 
  ggtitle("Variable Importance: Hurricane Outage Duration") + 
  #ggtitle("BART Variable Importance: Hurricane Max% Impacted") + 
  labs(fill = "Variable Type", color = "Variable Type") +
  guides(fill = guide_legend(byrow = T), color = guide_legend(byrow = T)) +
  theme(
        plot.title = element_text(hjust = 0.5),
        legend.spacing.y = unit(0.33, "cm"),
        legend.position = c(.8, .4)
  )
pp2


df_imp_bart = data.frame(as.list(bart_vimp[[1]])) %>% 
  gather(key = "Feature", value = "Inclusion") %>%
  arrange(desc(Inclusion)) %>%
  mutate(cat = case_when(
    Feature %in% c("delta_WIND10M", "delta_vapor", "delta_T10M", "delta_pressure", 
                   "max_humidity_10M", "max_T10M", "max_WIND10M",
                   "spi03", "spi06", "spi12", "spi24") ~ "Weather",
    Feature %in% c("Barren", "Cultivated", "Forest", "Herbaceous", "Shrub", "Water", "Wetlands",
                   "DEM_sd", "DEM_min", "RZ_med", "RZ_mode") ~ "Environment",
    TRUE ~ "Socio-Economic"
  )) %>%
  mutate(Feature = case_when(Feature == "max_humidity_10M" ~ "max_humidity10M", TRUE ~ as.character(Feature))) %>% #name correction
  dplyr::slice(1:20) %>%
  mutate(Feature1 = str_replace(Feature, "_", " ")) %>%
  mutate(Feature2 = tolower(Feature1)) %>%
  mutate(Feature_clean = str_to_title(Feature2)) %>%
  dplyr::select(-c(Feature1, Feature2))
bb2 = ggplot() + 
  theme_classic() + 
  geom_col(data = df_imp_bart, aes(x = Inclusion, y = fct_reorder(Feature_clean, Inclusion), fill = cat, color = cat), 
           alpha = .75) +
  scale_fill_manual(values = fill_vec) +
  scale_color_manual(values = fill_vec) +
  xlab("Inclusion") +
  ylab(element_blank()) + 
  ggtitle("BART Variable Importance: Hurricane Outage Duration") + 
  #ggtitle("BART Variable Importance: Hurricane Max% Impacted") + 
  labs(fill = "Variable Type", color = "Variable Type") +
  guides(fill = guide_legend(byrow = T), color = guide_legend(byrow = T)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0.33, "cm"),
    legend.position = c(.8, .4)
  )
bb2

##########################################################################################################
##### PARTIAL DEPENDENCE PLOTS (PDPs) - FINAL MODEL ######################################################
##########################################################################################################
pdp1 = pdp::partial(final_obj, pred.var = "delta_WIND10M", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp1 = ggplot() +
  theme_gray() +
  geom_line(data = pdp1, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = delta_WIND10M), sides = "b", alpha = 0.33, color = "black") +
  xlab("Delta Wind10m")
ggpdp1

pdp2 = pdp::partial(final_obj, pred.var = "delta_vapor", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp2 = ggplot() +
  theme_gray() +
  geom_line(data = pdp2, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = delta_vapor), sides = "b", alpha = 0.33, color = "black") +
  xlab("Delta Vapor")
ggpdp2

pdp3 = pdp::partial(final_obj, pred.var = "delta_T10M", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp3 = ggplot() +
  theme_gray() +
  geom_line(data = pdp3, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = delta_T10M), sides = "b", alpha = 0.33, color = "black") +
  xlab("Delta T10m")
ggpdp3

pdp4 = pdp::partial(final_obj, pred.var = "delta_pressure", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp4 = ggplot() +
  theme_gray() +
  geom_line(data = pdp4, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = delta_pressure), sides = "b", alpha = 0.33, color = "black") +
  xlab("Delta Pressure")
ggpdp4

grid.arrange(ggpdp1, ggpdp2, ggpdp3, ggpdp4, ncol = 2,
             top = textGrob("Partial Dependence Plots: XGBoost Outage Duration"))

##########################################################################################################
##### CREATE MAPS  #######################################################################################
##########################################################################################################
# # Base county map 
# load(file = "Data/processed/county_map_proj.Rda")
# 
# # Hurricane exposure
# hurricane_exposure = county_wind(counties = as.vector(county_map_proj$GEOID), 
#                                  start_year = 2015, end_year = 2019, wind_limit = 17.4)
# hurricane_group = hurricane_exposure %>%
#   group_by(fips) %>%
#   summarize(avg_sust_spd = mean(vmax_sust), avg_gust_spd = mean(vmax_gust), sum_sust_hrs = sum(sust_dur) / 60, sum_gust_hrs = sum(gust_dur) / 60)
# 
# # Annual outage hours - all US
# hrs = df_data %>%
#   dplyr::select(GEOID, ln_hrs) %>%
#   mutate(hrs = exp(ln_hrs)) %>%
#   group_by(GEOID) %>%
#   summarise(Hours = sum(hrs))
# hr_map = county_map_proj %>%
#   left_join(hrs, by = c("GEOID")) %>%
#   left_join(hurricane_group, by = c("GEOID" = "fips"))
# 
# gg3 = ggplot()+
#   geom_sf(data = hr_map, aes(fill = Hours), color = NA) +
#   #geom_sf(data = hr_map, aes(fill = sum_sust_hrs), color = NA) +
#   scale_fill_viridis_c(option="plasma", na.value = "grey30") +
#   geom_sf(data = county_map_proj, fill = NA, color = "black", lwd = 0.1) + 
#   theme_dark() +
#   labs(title = "Power Outages\nHurricanes (2015 - 2019)", fill = "Hours") +
#   #labs(title = "Duration Sustained Winds > 20 m/s\nHurricanes (2015-2019)", fill = "Hours") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#   )
# print(gg3)
# 
# # Annual outage hours - zoom 
# county_zoom_list = c("TX", "LA", "MS", "AL", "FL", "GA", "SC", "NC", "VA")
# county_zoom = get_acs(geography = "county", state = county_zoom_list,
#                       variables=c("B01003_001"), year = 2019, geometry = TRUE, 
#                       cache_table = TRUE) %>%
#   mutate(POPULATION = estimate) %>%
#   dplyr::select(GEOID, NAME) %>%
#   st_transform(5070) 
# hr_zoom = county_zoom %>%
#   left_join(hrs, by = c("GEOID")) %>%
#   left_join(hurricane_group, by = c("GEOID" = "fips"))
# 
# gg3 = ggplot()+
#   geom_sf(data = hr_zoom, aes(fill = Hours), color = NA) +
#   #geom_sf(data = hr_zoom, aes(fill = sum_sust_hrs), color = NA) +
#   scale_fill_viridis_c(option="plasma", na.value = "grey30") +
#   geom_sf(data = county_zoom, fill = NA, color = "black", lwd = 0.1) + 
#   theme_dark() +
#   labs(title = "Power Outages\nHurricanes (2015 - 2019)", fill = "Hours") +
#   #labs(title = "Duration Sustained Winds > 20 m/s\nHurricanes (2015-2019)", fill = "Hours") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#   )
# print(gg3)


##########################################################################################################
##### CHECK MODEL ASSUMPTIONS  ###########################################################################
##########################################################################################################
# Residuals
res_xgb = as.vector(xgb_predictions$.pred) - y_test
res_mlr = as.vector(mlr_predictions) - y_test
res = res_mlr

# Normality of residuals 
# http://www.sthda.com/english/wiki/normality-test-in-r
qqnorm(res)
qqline(res)
res_shap= shapiro.test(res)
gg_norm = ggpubr::ggqqplot(res, color = '#00AFBB') +
  labs(title = "Normal Q-Q Plot") +
    theme(plot.title = element_text(hjust = 0.5)
          ) +
  annotate("text", x = -1.5, y = 0.6, label = paste("Shapiro-Wilk: ", format(round(res_shap$p.value, 3), nsmall = 3), sep = ""), color = '#00AFBB')
gg_norm

# Heteroskedacity of residuals
fitted = as.vector(xgb_predictions$.pred) 
plot(fitted, res_xgb)
gg_hetero_data = data.frame(fitted, res)
gg_hetero = ggplot(gg_hetero_data) +
  theme_classic() +
  geom_point(aes(x = fitted, y = res), color = "#E7B800") +
  geom_hline(yintercept = 0, linetype="dashed", color = "#E7B800", alpha = 0.85) +
  xlab("Fitted") +
  ylab("Residuals") +
  labs(title = "Residuals vs. Fitted") +
  theme(plot.title = element_text(hjust = 0.5)
  ) 
gg_hetero  

bart_assmp = check_bart_error_assumptions(bart_fit)

# Spatial dependency of residuals
# Base county map 
load(file = "Data/processed/county_map_proj.Rda")
rr = data.frame(as.factor(df_test$GEOID))
rr$res = res
colnames(rr) = c("GEOID", "resid")
plot(rr) # looks good but use varigoram below to confirm 

county_centroid = st_centroid(county_map_proj) # get center of counties 
county_lonlat = county_centroid %>% 
  mutate(X = unlist(map(county_centroid$geometry,1)),
         Y = unlist(map(county_centroid$geometry,2))) %>%
  dplyr::select(-NAME, -POPULATION) %>%
  inner_join(rr, by = c("GEOID")) %>%
  rename(Z = resid)
county_lonlat_sp = as_Spatial(county_lonlat)
vgram = variogram(Z~1, county_lonlat_sp)
plot(vgram)

gg_vgram_data = data.frame(Distance = vgram$dist / 1000, Gamma = vgram$gamma) #distance from m to km
gg_vgram = ggplot(gg_vgram_data, aes(x = Distance, y = Gamma)) +
  geom_line(linetype = "dashed", color = "#FC4E07") + 
  geom_point(color = "#FC4E07") +
  theme_classic() +
  xlab("Distance (km)") + 
  ylim(0, round(1.1* max(vgram$gamma), 2))+
  labs(title = "Variogram of Residuals") +
  theme(plot.title = element_text(hjust = 0.5)
  ) 
gg_vgram  


########################################################################################################
### FIGURES FOR PUBLICATION
########################################################################################################
pdf_x = 8.5 #sra 8.5 x 3.65 in
pdf_y = 3.65

#Figure model accuracy
gg_acc = ggplot() + 
  theme_classic() + 
  geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
  geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
  scale_color_manual(
    values = color_vec, 
    labels = c("Actual",
               bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * ")"),
               bquote("MLR (" * R^2 ~ "=" ~ .(rsq_mlr) * ")"), 
               bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * ")")
    ),
    name = element_blank()) +
  scale_linetype_manual(
    values = lty_vec,
    labels = c("Actual",
               bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * ")"),
               bquote("MLR (" * R^2 ~ "=" ~ .(rsq_mlr) * ")"), 
               bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * ")")         
    ),
    name = element_blank()) + 
  scale_alpha_manual(
    values = alpha_vec,
    labels = c("Actual",
               bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * ")"),
               bquote("MLR (" * R^2 ~ "=" ~ .(rsq_mlr) * ")"), 
               bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * ")")
    ),
    name = element_blank()) + 
  scale_y_continuous(labels = function(x) paste0(x)) +
  xlab("Index (County x Tropical Cyclone)") +
  ylab("Hours (ln)") + 
  ggtitle("Outage Duration: Test Sample") + 
  #ylab("Max %Customers Out (ln)") + 
  #ggtitle("Customers w/o Power: Test Sample") + 
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
gg_acc
#gg_acc_pct = gg_acc
#save(gg_acc_pct, file = "Figures/Paper/gg_acc_pct.Rda")
load("Figures/Paper/gg_acc_pct.Rda")
cow_acc = cowplot::plot_grid(gg_acc, gg_acc_pct, ncol = 2, scale = 0.95)
cow_acc
pdf("Figures/Paper/test.pdf", width = pdf_x, height = pdf_y) #SRA ppt
cow_acc
dev.off()

#Figure model assumptions
cow_asmp1 = cowplot::plot_grid(gg_norm, gg_hetero, nrow = 2)
cow_asmp = cowplot::plot_grid(cow_asmp1, gg_vgram, ncol = 2, rel_widths = c(3,2), scale = 0.95)
cow_asmp
dev.off()
pdf("Figures/Publish/assumptions.pdf", width = pdf_x, height = pdf_y) #SRA ppt 8.5 x 3.65
cow_asmp
dev.off()

#Figure variable importance 
df_vv = df_imp %>%
  dplyr::slice(1:20)
fill_vec = c("#CEB966", "#A379BB", "#6BB1C9")  
vv2 = ggplot() + 
  theme_classic() + 
  geom_col(data = df_vv, aes(x = Gain, y = fct_reorder(Feature_clean, Gain), fill = cat, color = cat), 
           alpha = .75) +
  scale_fill_manual(values = fill_vec) +
  scale_color_manual(values = fill_vec) +
  xlab("Relative Gain") +
  ylab(element_blank()) + 
  ggtitle("XGBoost: Outage Duration") + 
  labs(fill = "Variable Type", color = "Variable Type") +
  guides(fill = guide_legend(byrow = T), color = guide_legend(byrow = T)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0.33, "cm"),
    legend.position = c(.8, .4)
  )
vv2

vv2mlr = ggplot() + 
  theme_classic() + 
  geom_col(data = df_imp_mlr, aes(x = Importance, y = fct_reorder(Feature_clean, Importance), fill = cat, color = cat), 
           alpha = .75) +
  scale_fill_manual(values = fill_vec) +
  scale_color_manual(values = fill_vec) +
  xlab("Absolute t-value") +
  ylab(element_blank()) + 
  ggtitle("MLR: Outage Duration") + 
  labs(fill = "Variable Type", color = "Variable Type") +
  guides(fill = guide_legend(byrow = T), color = guide_legend(byrow = T)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0.33, "cm"),
    legend.position = c(.8, .4)
  )
vv2mlr

#vv2_pct = vv2
#save(vv2_pct, file = "Figures/Publish/vv2_pct.Rda")
#load("Figures/Publish/vv2_pct.Rda")
cow_vimp = cowplot::plot_grid(vv2, vv2mlr, ncol = 2, scale = 0.95)
cow_vimp
pdf("Figures/Paper/vimp.pdf", width = pdf_x, height = pdf_y) #SRA ppt 8.5 x 3.65
cow_vimp
dev.off()

## Figure Five factors
ff1 = gg_acc + ggtitle("Hurricane Outage Duration: Test Sample\n--Five Factors--")
ff2 = vv2 + ggtitle("Hurricane Outage Duration\n--Five Factors--")
cow_ff = cowplot::plot_grid(ff1, ff2, ncol = 2, scale = 0.95)
cow_ff
pdf("Figures/Publish/five.pdf", width = pdf_x, height = pdf_y) #SRA ppt 8.5 x 3.65
cow_ff
dev.off()




