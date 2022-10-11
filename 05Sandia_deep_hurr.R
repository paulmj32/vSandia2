#### Machine learning 
library(tidyverse)
library(tidymodels)
library(sf)
library(parallel)
library(doParallel)
library(xgboost)
library(vip)
library(spdep)
library(spatialreg)
library(pdp)
# library(drat)
# addRepo("geanders")
# install.packages("hurricaneexposuredata")
library(hurricaneexposuredata)
library(hurricaneexposure)
library(tidycensus)
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

# Load data and select final DVs 
load(file = "Data/processed/sf_data_ALL_nototal.Rda")

sf_data = sf_data_ALL %>%
  dplyr::filter(hurricanes >= 1) %>% #filter to event of interest
  mutate(ln_hrs = log(duration_hr)) %>% 
  mutate(ln_cust = log(max_cust_out)) %>% 
  mutate(pct_cust = max_frac_cust_out) %>%
  dplyr::select(-c(POPULATION, mean_cust_out, mean_frac_cust_out, max_cust_out, max_frac_cust_out, 
                   duration_hr, all_of(our_events), avalanche)) %>%
  relocate(c(ln_hrs, ln_cust, pct_cust))

rm(list=c("sf_data_ALL")) 
gc() 

df_data = sf_data %>%
  st_set_geometry(NULL) 

# Split into training vs testing
set.seed(23)
df_split = initial_split(df_data, prop = 0.80, strata = "ln_hrs")
df_train = training(df_split)
df_test = testing(df_split)  
df_cv = vfold_cv(df_train, v = 10, repeats = 1)
df_preds = df_data %>% dplyr::select(-c(ln_hrs, ln_cust, pct_cust, GEOID))

# Recipes for tidymodels 
recipe_hrs = recipe(ln_hrs ~ . , data = df_data) %>% step_rm(ln_cust, pct_cust, GEOID)
recipe_pct = recipe(pct_cust ~ . , data = df_data) %>% step_rm(ln_hrs, ln_cust, GEOID) %>% step_naomit(pct_cust)

################################################################################
#### MACHINE LEARNING ##########################################################
################################################################################
### Define which recipe you want to use 
recipe_mine = recipe_hrs
y_test = df_test$ln_hrs

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
    scale_y_continuous(labels = function(x) paste0(x)) +
    xlab("County x Event Index") +
    ylab("Hours (ln)") + 
    ggtitle("Hurricane Outage Duration: Test Sample") + 
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

final_model = xgb_work %>%
  finalize_workflow(xgb_best) %>%
  fit(df_data)
final_obj = extract_fit_parsnip(final_model)$fit
importance = xgb.importance(model = final_obj)
head(importance, n = 20) 
vip(final_obj, n = 20)

final_lre = lre_work %>%
  finalize_workflow(lre_best) %>%
  fit(df_data)
print(tidy(final_lre), n = 150)

##########################################################################################################
##### PARTIAL DEPENDENCE PLOTS (PDPs) - FINAL MODEL ######################################################
##########################################################################################################
X_train = df_train %>% dplyr::select(-c(GEOID, ln_hrs, ln_cust, pct_cust))
# pdp1 = pdp::partial(final_obj, pred.var = "Wetlands", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# pdp2 = pdp::partial(final_obj, pred.var = "DEM_sd", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# pdp3= pdp::partial(final_obj, pred.var = "DEM_min", ice = F, center = F,
#                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                    train = X_train, type = "regression")
# pdp4 = pdp::partial(final_obj, pred.var = "QMOHO", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# pdp5 = pdp::partial(final_obj, pred.var = "GENDINC", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# pdp6 = pdp::partial(final_obj, pred.var = "QBLACK", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# pdp7 = pdp::partial(final_obj, pred.var = "PATTACHRES", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# pdp8 = pdp::partial(final_obj, pred.var = "QHLTH65", ice = F, center = F,
#                     plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
#                     train = X_train, type = "regression")
# 
# grid.arrange(pdp1, pdp2, pdp3, pdp4, ncol = 2)
# grid.arrange(pdp5, pdp6, pdp7, pdp8, ncol = 2)


##########################################################################################################
##### CREATE MAPS  #######################################################################################
##########################################################################################################
# Base county map 
load(file = "Data/processed/county_map_proj.Rda")

# Hurricane exposure
hurricane_exposure = county_wind(counties = as.vector(county_map_proj$GEOID), 
                                 start_year = 2015, end_year = 2019, wind_limit = 17.4)
hurricane_group = hurricane_exposure %>%
  group_by(fips) %>%
  summarize(avg_sust_spd = mean(vmax_sust), avg_gust_spd = mean(vmax_gust), sum_sust_hrs = sum(sust_dur) / 60, sum_gust_hrs = sum(gust_dur) / 60)

# Annual outage hours - all US
hrs = df_data %>%
  dplyr::select(GEOID, ln_hrs) %>%
  mutate(hrs = exp(ln_hrs)) %>%
  group_by(GEOID) %>%
  summarise(Hours = sum(hrs))
hr_map = county_map_proj %>%
  left_join(hrs, by = c("GEOID")) %>%
  left_join(hurricane_group, by = c("GEOID" = "fips"))

gg3 = ggplot()+
  geom_sf(data = hr_map, aes(fill = Hours), color = NA) +
  #geom_sf(data = hr_map, aes(fill = sum_sust_hrs), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey30") +
  geom_sf(data = county_map_proj, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  labs(title = "Power Outages\nHurricanes (2015 - 2019)", fill = "Hours") +
  #labs(title = "Duration Sustained Winds > 20 m/s\nHurricanes (2015-2019)", fill = "Hours") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)

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
  geom_sf(data = hr_zoom, aes(fill = Hours), color = NA) +
  #geom_sf(data = hr_zoom, aes(fill = sum_sust_hrs), color = NA) +
  scale_fill_viridis_c(option="plasma", na.value = "grey30") +
  geom_sf(data = county_zoom, fill = NA, color = "black", lwd = 0.1) + 
  theme_dark() +
  labs(title = "Power Outages\nHurricanes (2015 - 2019)", fill = "Hours") +
  #labs(title = "Duration Sustained Winds > 20 m/s\nHurricanes (2015-2019)", fill = "Hours") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg3)


##########################################################################################################
##### CHECK MODEL ASSUMPTIONS  ###########################################################################
##########################################################################################################
# Residuals
res_xgb = as.vector(xgb_predictions$.pred) - y_test
res_lre = as.vector(lre_predictions$.pred) - y_test
res = res_xgb

# Normality of residuals 
# http://www.sthda.com/english/wiki/normality-test-in-r
qqnorm(res)
qqline(res)
shapiro.test(res)

# Heteroskedacity of residuals
fitted = as.vector(xgb_predictions$.pred) 
plot(fitted, res_xgb)

# Spatial dependency of residuals
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




