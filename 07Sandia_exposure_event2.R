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
library(DBI)
library(lubridate)
library(grid)
library(gridExtra)
library(cowplot)

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
df_data1 = sf_data_ALL %>%
  dplyr::select(-c(POPULATION,
                   mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out, 
                   duration_hr,
                   avalanche, #erroneous
                   all_of(our_events),
                   max_WIND10M:spi24
                   )) %>%
  st_set_geometry(NULL) %>%
  group_by(GEOID) %>%
  summarise(across(everything(), mean))
rm(list=c("sf_data_ALL")) 
gc() 

## Get SQL data for NOAA event hours
# https://stackoverflow.com/questions/9802680/importing-files-with-extension-sqlite-into-r
# https://cran.r-project.org/web/packages/RSQLite/vignettes/RSQLite.html
# https://www.r-bloggers.com/2018/11/quick-guide-r-sqlite/
# my_statement = 'SELECT fips_code, SUM(droughts), SUM(extreme_cold), SUM(extreme_heat), SUM(floods), SUM(hails), SUM(high_winds), SUM(hurricanes), SUM(tornadoes), SUM(wildfires), SUM(winter_storms), SUM(avalanche), SUM(blizzard), SUM(coastal_flood), SUM(cold_wind_chill), SUM(dense_smoke), SUM(drought), SUM(excessive_heat), SUM(extreme_cold_wind_chill), SUM(flash_flood), SUM(flood), SUM(frost_freeze), SUM(funnel_cloud), SUM(hail), SUM(heat), SUM(heavy_rain), SUM(heavy_snow), SUM(high_wind), SUM(hurricane), SUM(hurricane_typhoon), SUM(ice_storm), SUM(lake_effect_snow), SUM(lakeshore_flood), SUM(sleet), SUM(strong_wind), SUM(tornado), SUM(tropical_depression), SUM(tropical_storm), SUM(wildfire), SUM(winter_storm), SUM(winter_weather) FROM noaa GROUP BY fips_code'
# dbpath = '/Users/paul/vSandia2/Data/Outages/AGM_noaa_CONUS_panel_Proc03Aug2022.sqlite'
# mydb = RSQLite::dbConnect(drv = RSQLite::SQLite(), dbname = dbpath)
# tables = dbListTables(mydb) #see tables
# myquery_10 = dbGetQuery(mydb, 'SELECT * FROM noaa LIMIT 10')
# #myquery = dbGetQuery(mydb, statement = my_statement)
# dbDisconnect(mydb)
# 
# #save(myquery, file = "Data/processed/myquery_hazards.Rda")
# load(file = "Data/processed/myquery_hazards.Rda")
# columns = colnames(myquery_10[,3:ncol(myquery_10)])
# df_hazards = myquery
# colnames(df_hazards) = c("GEOID", columns)
# df_hazards_long = df_hazards %>%
#   gather(key = "hazard", value = "minutes", -GEOID) %>%
#   mutate(hours = minutes / 60)
# df_hazards_clean = df_hazards_long %>%
#   mutate(hazard2 = case_when(
#     hazard %in% c('drought') ~ "droughts",
#     hazard %in% c('cold_wind_chill', 'extreme_cold_wind_chill', 'frost_freeze') ~ "extreme_cold",
#     hazard %in% c('excessive_heat', 'heat') ~ "extreme_heat",
#     hazard %in% c('coastal_flood', 'flash_flood', 'flood', 'heavy_rain', 'lakeshore_flood') ~ "floods",
#     hazard %in% c('hail') ~ "hails",
#     hazard %in% c('high_wind', 'strong_wind') ~ "high_winds",
#     hazard %in% c('hurricane', 'hurricane_typhoon', 'tropical_depression', 'tropical_storm') ~ "hurricanes",
#     hazard %in% c('funnel_cloud', 'tornado') ~ "tornadoes",
#     hazard %in% c('dense_smoke', 'wildfire') ~ "wild_fires",
#     hazard %in% c('blizzard', 'avalanche', 'heavy_snow', 'ice_storm', 'lake_effect_snow', 'sleet', 'winter_storm', 'winter_weather') ~ "winter_storms"
#   )) %>%
#   dplyr::select(-minutes) %>%
#   drop_na() #remove binary indicators 

#save(df_hazards_clean, file = "Data/processed/df_hazards_clean.Rda")
load(file = "Data/processed/df_hazards_clean.Rda")
hazard_count = df_hazards_clean %>%
  group_by(hazard2) %>%
  summarise(Hours = sum(hours)) 

# Data-frame of hazard of interest
df_hazard = df_hazards_clean %>%
  dplyr::select(-hazard) %>%
  group_by(GEOID, hazard2) %>%
  summarise(Hours = sum(hours)) %>%
  dplyr::filter(hazard2 == my_event)
rm(list=c("df_hazards_clean")) 
gc() 

# Combine data for final data-frame 
df_data = df_data1 %>%
  inner_join(df_hazard, by = c("GEOID")) %>%
  dplyr::filter(Hours > 0) %>%
  mutate(lnHrs = log(Hours / 5 + 0.001)) %>% #take log of annualized hours
  dplyr::select(-c(hazard2, Hours)) %>%
  relocate(lnHrs)

df_preds = df_data %>% 
  dplyr::select(-c(lnHrs, GEOID))

# Split into training vs testing
set.seed(32)
df_split = initial_split(df_data, prop = 0.80, strata = lnHrs)
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)

# Recipe for ML
recipe_mine = recipe(lnHrs ~ . , data = df_data) %>% 
  step_rm(GEOID) %>% 
  step_naomit()

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
y_train = df_train %>% pull(lnHrs) %>% na.omit() 
y_test = df_test %>% pull(lnHrs) %>% na.omit() 
X_train = df_train %>% dplyr::select(-c(GEOID, lnHrs)) %>% na.omit()
X_test = df_test %>% dplyr::select(-c(GEOID, lnHrs)) %>% na.omit()

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


#regular MLR
df_mlr =  X_train %>%
  mutate(y = y_train)
mlr = lm(y ~ ., data = df_mlr)
summary(mlr)
mlr_predictions = predict(mlr, newdata = X_test)
rsq_mlr = format(round(cor(mlr_predictions, y_test)^2, 3), nsmall = 3)
asd = vip(mlr, n = 20)
imp.new = asd$data$Importance/sum(asd$data$Importance)
asd$data$Importance = imp.new
plot(asd)

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
                   #eNet = as.vector(lre_predictions$.pred),
                   bart = as.vector(bart_predictions),
                   mlr = as.vector(mlr_predictions),
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
    ylab("Annual Hours (ln)") + 
    ggtitle(paste("Exposure to ", str_to_title(gsub("_", " ", my_event)), ": Test Sample", sep = "")) +
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

fill_vec = c("#CEB966", "#A379BB", "#6BB1C9")  
pp2 = ggplot() + 
  theme_classic() + 
  geom_col(data = df_imp, aes(x = Gain, y = fct_reorder(Feature_clean, Gain), fill = cat, color = cat), 
           alpha = .75) +
  scale_fill_manual(values = fill_vec) +
  scale_color_manual(values = fill_vec) +
  xlab("Importance") +
  ylab(element_blank()) + 
  ggtitle(paste("Variable Importance (XGB): Exposure to ", str_to_title(gsub("_", " ", my_event)), sep = "")) +
  labs(fill = "Variable Type", color = "Variable Type") +
  guides(fill = guide_legend(byrow = T), color = guide_legend(byrow = T)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0.33, "cm"),
    legend.position = c(.8, .4)
  )
pp2

##########################################################################################################
##### EXPOSURE MAPS  #####################################################################################
##########################################################################################################
# Base map
load(file = "Data/processed/county_map_proj.Rda")

# Annual outage hours - all US log(Hours / 5 + 0.001)
hrs = df_data %>%
  dplyr::select(GEOID, lnHrs) %>%
  mutate(Hours = exp(lnHrs) - 0.001) %>%
  group_by(GEOID) 

county_hrs = county_map_proj %>%
  left_join(hrs, by = c("GEOID"))

county_plot = county_hrs
gg3 = ggplot()+
  geom_sf(data = county_plot, aes(fill = Hours), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey30") +
  geom_sf(data = county_plot, fill = NA, color = "black", lwd = 0.01) + 
  theme_dark() +
  #labs(title = paste("Annual ", str_remove(str_to_title(gsub("_", " ", my_event)),"s"), " Exposure (2015 - 2019)", sep = ""), fill = "Hours") +
  #labs(title = expression(atop("Annual Power Outages", paste("Tropical Storms (2015 - 2019)"))), fill = "Hours") +
  labs(title = expression(atop("Annual Exposure to Sust."~Winds>=34~"kt", paste("Tropical Cyclones (2015 - 2019)"))), fill = "Hours") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)

##########################################################################################################
##### PARTIAL DEPENDENCE PLOTS (PDPs) - FINAL MODEL ######################################################
##########################################################################################################
pdp = pdp::partial(final_obj, pred.var = "QHLTH", ice = F, center = F,
                   plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                   train = X_train, type = "regression")
pdp

pdp1 = pdp::partial(final_obj, pred.var = "Forest", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp1 = ggplot() +
  theme_gray() +
  geom_line(data = pdp1, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = Forest), sides = "b", alpha = 0.33, color = "black") +
  xlab("Forest")
ggpdp1

pdp2 = pdp::partial(final_obj, pred.var = "DEM_min", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp2 = ggplot() +
  theme_gray() +
  geom_line(data = pdp2, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = DEM_min), sides = "b", alpha = 0.33, color = "black") +
  xlab("Dem Min")
ggpdp2

pdp3 = pdp::partial(final_obj, pred.var = "PROXMET", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1) 
ggpdp3 = ggplot() +
  theme_gray() +
  geom_line(data = pdp3, aes(x = IV / 1000, y = yhat)) +
  geom_rug(data = X_train, aes(x = PROXMET / 1000), sides = "b", alpha = 0.33, color = "black") +
  xlab("Proxmet")
ggpdp3

pdp4 = pdp::partial(final_obj, pred.var = "QSPANISH", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp4 = ggplot() +
  theme_gray() +
  geom_line(data = pdp4, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = QSPANISH), sides = "b", alpha = 0.33, color = "black") +
  xlab("Qspanish")
ggpdp4

pdp5 = pdp::partial(final_obj, pred.var = "QHLTH65", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp5 = ggplot() +
  theme_gray() +
  geom_line(data = pdp5, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = QHLTH65), sides = "b", alpha = 0.33, color = "black") +
  xlab("Qhlth65")
ggpdp5

pdp6 = pdp::partial(final_obj, pred.var = "CROPINS", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp6 = ggplot() +
  theme_gray() +
  geom_line(data = pdp6, aes(x = IV * 100, y = yhat)) +
  geom_rug(data = X_train, aes(x = CROPINS * 100), sides = "b", alpha = 0.33, color = "black") +
  xlab("Cropins") 
ggpdp6

pdp7 = pdp::partial(final_obj, pred.var = "GENDINC", ice = F, center = F,
                    plot = F, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression") %>%
  rename(IV = 1)
ggpdp7 = ggplot() +
  theme_gray() +
  geom_line(data = pdp7, aes(x = IV, y = yhat)) +
  geom_rug(data = X_train, aes(x = GENDINC), sides = "b", alpha = 0.33, color = "black") +
  xlab("Gendinc")
ggpdp7

ggpdps=grid.arrange(ggpdp3, ggpdp4, ggpdp5, ggpdp6, ncol = 2,
             top = textGrob("Partial Dependence Plots: Socio-Economic"))



########################################################################################################
### FIGURES FOR PUBLICATION
########################################################################################################
pdf_x = 8.5 #sra 8.5 x 3.65 in
pdf_y = 3.65

#Figure model accuracy
gg_acc  = ggplot() + 
  theme_classic() + 
  geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
  geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
  scale_color_manual(
    values = color_vec, 
    labels = c("Actual",
               bquote("BART (" * R^2 ~ "=" ~ .(rsq_bart) * ")"),
               bquote("MLR (" * R^2 ~ "=" ~ .(rsq_mlr)  * ")"), 
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
  xlab("County") +
  ylab("Hours (ln)") + 
  ggtitle(paste("Exposure: Test Sample", sep = "")) +
  # guides(
  #   color = guide_legend(order = 2),
  #   shape = guide_legend(order = 1),
  #   linetype = guide_legend(order = 2)
  # ) + 
  theme(legend.spacing.y = unit(-0.25, "cm"),
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.position = c(.8, .2),
        plot.title = element_text(hjust = 0.5)
  )
gg_acc
cow_haz = cowplot::plot_grid(gg3, gg_acc, ncol = 2, scale = 0.95, rel_widths = c(1.33,1))
pdf("Figures/Paper/hazard.pdf", width = pdf_x, height = pdf_y) #SRA ppt
cow_haz
dev.off()

# Figure Variable importance and PDPs
df_vv = df_imp %>%
  dplyr::slice(1:20) 
fill_vec = c("#CEB966", "#A379BB", "#6BB1C9")  
vv2 = ggplot() + 
  theme_classic() + 
  geom_col(data = df_vv, aes(x = Gain, y = fct_reorder(Feature_clean, Gain), fill = cat, color = cat), 
           alpha = .75) +
  scale_fill_manual(values = fill_vec) +
  scale_color_manual(values = fill_vec) +
  xlab("Variable Importance") +
  ylab(element_blank()) + 
  ggtitle("Tropical Cyclone Exposure") + 
  labs(fill = "Variable Type", color = "Variable Type") +
  guides(fill = guide_legend(byrow = T), color = guide_legend(byrow = T)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0.33, "cm"),
    legend.position = c(.8, .4)
  )
vv2

cow_pdps = cowplot::plot_grid(vv2, ggpdps, ncol = 2, scale = 0.95)
cow_pdps
pdf("Figures/Publish/pdps.pdf", width = pdf_x, height = pdf_y) #SRA ppt
cow_pdps
dev.off()

