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
y_test = df_test$ln_cust

rsq_lre = paste(lre_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_xgb = paste(xgb_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
rsq_rf = paste(rf_test %>% dplyr::filter(.metric == "rsq") %>% pull(.estimate) %>% round(3) %>% format(nsmall = 3))
#NOTE: the eNet Rs-squared is 0.736 (hrs) without without the large event

cverror_lre = paste(show_best(lre_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_xgb = paste(show_best(xgb_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))
cverror_rf = paste(show_best(rf_tune, metric = "rmse") %>% dplyr::slice(1) %>% pull(mean) %>% round(3) %>% format(nsmall = 3))

### PLOTTING - Max Customer Outages
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
    ylab("Max Customers Out (ln)") + 
    scale_y_continuous(labels = function(x) paste0(x), limits = c(0,15)) +
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
