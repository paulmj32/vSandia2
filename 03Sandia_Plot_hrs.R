#### PLOT AND ANALYZE RESULTS OF ML 
library(tidyverse)
library(tidymodels)
library(sf)
library(parallel)
library(doParallel)
library(xgboost)
library(vip)
library(spdep)
library(pdp)


## Testing results
#y_test = df_hours %>% slice(-df_split$in_id) %>% dplyr::select(ln_hrs) %>% pull()
y_test = df_test$ln_hrs

rsq_lre = paste(lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_xgb = paste(xgb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_rf = paste(rf_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
#NOTE: the eNet Rs-squared is 0.736 (hrs) without without the large event

cverror_lre = paste(show_best(lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_xgb = paste(show_best(xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_rf = paste(show_best(rf_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

### PLOTTING - Duration Hours
gg = dplyr::tibble(actual = y_test,
                   eNet = as.vector(lre_predictions$.pred),
                   rf = as.vector(rf_predictions$.pred),
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
color_vec = c("black", "#FDE725FF", "#35B779FF", "#440154FF")
lty_vec = c(1, 1, 1)
lty_vec = c(1, 1, 1, 1)
alpha_vec = c(1, 0.7, 0.7)
alpha_vec = c(1, 0.7, 0.7, 0.7)

plot_filtering_estimates2 <- function(df) {
  p = ggplot() + 
    theme_classic() + 
    geom_hline(yintercept = mean(gg$actual, na.rm = T), linetype="dashed", color = "gray50", alpha = 0.85) +
    geom_line(data = gg_all, aes(x = index, y = ypred, color = Model, lty = Model, alpha = Model)) +
    scale_color_manual(
      values = color_vec, 
      labels = c("Actual",
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
      ),
      name = element_blank()) +
    scale_linetype_manual(
      values = lty_vec,
      labels = c("Actual",
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")         
      ),
      name = element_blank()) + 
    scale_alpha_manual(
      values = alpha_vec,
      labels = c("Actual",
                 bquote("eNET (" * R^2 ~ "=" ~ .(rsq_lre) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_lre) * ")"), 
                 bquote("RF (" * R^2 ~ "=" ~ .(rsq_rf) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_rf) * ")"), 
                 bquote("XGB (" * R^2 ~ "=" ~ .(rsq_xgb) * "," ~ RMSE[cv] ~ "=" ~ .(cverror_xgb) * ")")
      ),
      name = element_blank()) + 
    ylab("Outage Duration (ln-hrs)") + 
    scale_y_continuous(labels = function(x) paste0(x), limits = c(2, 10)) +
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
X_train = df_train %>% dplyr::select(-c(ln_hrs, ln_cust, pct_cust, GEOID))
pdp1 = pdp::partial(final_obj, pred.var = "total_liquid", ice = F, center = F,
                          plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                          train = X_train, type = "regression")
pdp2 = pdp::partial(final_obj, pred.var = "delta_pressure", ice = F, center = F,
                            plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                            train = X_train, type = "regression")
pdp3= pdp::partial(final_obj, pred.var = "delta_T10M", ice = F, center = F,
                            plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                            train = X_train, type = "regression")
pdp4 = pdp::partial(final_obj, pred.var = "min_WIND10M", ice = F, center = F,
                          plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                          train = X_train, type = "regression")
pdp5 = pdp::partial(final_obj, pred.var = "total_ice", ice = F, center = F,
                           plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                           train = X_train, type = "regression")
pdp6 = pdp::partial(final_obj, pred.var = "delta_WIND10M", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp7 = pdp::partial(final_obj, pred.var = "min_vapor", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")
pdp8 = pdp::partial(final_obj, pred.var = "delta_vapor", ice = F, center = F,
                    plot = T, rug= T, alpha = 0.1, #plot.engine = "ggplot2",
                    train = X_train, type = "regression")

grid.arrange(pdp1, pdp2, pdp3, pdp4, ncol = 2)
grid.arrange(pdp5, pdp6, pdp7, pdp8, ncol = 2)

pdp1_4 = pdp::partial(final_obj, pred.var = c("total_liquid", "min_WIND10M"),
                  plot = T, chull= T, rug = T, plot.engine = "ggplot2",
                  train = X_train, type = "regression")

# p1 = pdp::partial(final_obj, pred.var = c("PS_sd", "SLP_mean"), 
#                   plot = T, chull= T, rug = T, plot.engine = "ggplot2",
#                   train = x, type = "regression")
# p2 = pdp::partial(final_obj, pred.var = c("PS_sd", "V_max"), 
#                   plot = T, chull= T, plot.engine = "ggplot2",
#                   train = x, type = "regression")
# p3 = pdp::partial(final_obj, pred.var = c("SLP_mean", "V_max"), 
#                   plot = T, chull= T, plot.engine = "ggplot2",
#                   train = x, type = "regression")

### View tree diagram 
# xgb.plot.tree(model = final_obj, trees = 1:2) 





## MAPS
df_plot = county_map_ALL %>%
  dplyr::filter(duration_hr >= 12) %>% #filter to events >12 hrs ... 95% quantile for full dataset
  dplyr::filter(duration_hr < 3000) #get rid of faulty data
X = df_plot %>%
  st_set_geometry(NULL) %>%
  dplyr::select(final_obj$feature_names) %>%
  as.matrix()
predictions = predict(final_obj, X)
df_plot2 = df_plot %>%
  dplyr::select(GEOID) %>%
  st_set_geometry(NULL) %>%
  mutate(fit = predictions) %>%
  group_by(GEOID) %>%
  summarise(fit_mean = mean(fit))


mycrs = 5070 #chose projected coordinate system: EPSG 5070 NAD83 Conus Albers
year=2019 # year for county boundaries 
options(tigris_use_cache = TRUE) #cache shapefiles for future sessions
state_list = c("AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
county_map = get_acs(geography = "county", state = state_list,
                     variables=c("B01003_001"), year = year, geometry = TRUE, 
                     cache_table = TRUE)
county_map = county_map %>%
  mutate(POPULATION = estimate) %>%
  dplyr::select(GEOID, NAME, POPULATION) 
county_map_proj = county_map %>% 
  st_transform(mycrs) # project to Conic Equal Area Albers, EPSG:5070 
county_map_area = county_map_proj %>%  
  mutate(AREA = as.vector(st_area(county_map_proj))) %>% #sq-meters; as.vector removes units suffix 
  mutate(DENSITY = POPULATION / AREA * 1000^2) %>% #population per sq-km 
  dplyr::select(-c(AREA, NAME))

county_map_risk = county_map_area %>%
  left_join(df_plot2, by = c("GEOID"))
gg2 = ggplot(county_map_risk)+
  geom_sf(aes(fill = fit_mean), color = NA) + 
  scale_fill_viridis_c(option="plasma", na.value = "grey50") +
  theme_dark() +
  labs(title = "Power Outage Prediction - All Events", fill = "Mean Duration\nln(Hours)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )
print(gg2)
