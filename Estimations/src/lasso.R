## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     March 7th, 2023
##
## This version:      March 7th, 2023
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##  Lasso Estimation                                                                                 ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    Set - UP                                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(pacman)
p_load(tidyverse, glmnet)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Download                                                                                            ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") %>%
  filter(year < 2013 & year > 2003) %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         dpto = as.factor(dpto),
         month = as.factor(month))

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Lasso                                                                                           ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define response variable
y <- merge_data.df %>%
  dplyr::select(ruv_desplazamiento_forzado) %>%
  pull()

# Violencia

violencia <- data.matrix(merge_data.df[, c('ruv_amenaza', "ruv_desaparicion_forzada", "ruv_combates", "ruv_secuestro", "ruv_perdida_bienes",
               "ruv_minas", "ruv_tortura", "ruv_homicidio", "ruv_reclutamiento_menores", "cnmh_violencia_sexual",
               "ruv_abandono_despojo", "cnmh_masacres", "ruv_violencia_sexual", "cnmh_atentado", "cnmh_d_bienes", 
               "cnmh_desaparicion", "cnmh_enfrentamientos", "cnmh_minas", "cnmh_secuestro", "cnmh_reclutamiento",
               "cnmh_ataque_poblacion",'cnmh_asesinatos_selectivos', "codmpio", "year")])

cv_model <- cv.glmnet(x = violencia, y = y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(violencia, y, alpha = 1, lambda = best_lambda)
coef(best_model)

violencia_lasso <- data.matrix(merge_data.df[, c('ruv_amenaza', 'cnmh_desaparicion', 'ruv_combates',
                                                 'cnmh_d_bienes', 'cnmh_minas', 'ruv_homicidio', 
                                                 'ruv_violencia_sexual', 'ruv_abandono_despojo',
                                                 'cnmh_reclutamiento', 'cnmh_atentado', 'cnmh_ataque_poblacion',
                                                 'codmpio', 'year')])

cv_model <- cv.glmnet(x = violencia_lasso, y = y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(violencia_lasso, y, alpha = 1, lambda = best_lambda)
final_violencia <- coef(best_model)

# Economicas

economicas <- data.matrix(merge_data.df[, c('night_lights', "vegetation", "rainFall", 'codmpio', 'year')])

cv_model <- cv.glmnet(x = economicas, y = y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(economicas, y, alpha = 1, lambda = best_lambda)
coef(best_model)

# Final

final_lasso <- data.matrix(merge_data.df[, c('ruv_amenaza', 'cnmh_desaparicion', 'ruv_combates',
                                                 'cnmh_d_bienes', 'cnmh_minas', 'ruv_homicidio', 
                                                 'ruv_violencia_sexual', 'ruv_abandono_despojo',
                                                 'cnmh_reclutamiento', 'cnmh_atentado', 'cnmh_ataque_poblacion',
                                                 'night_lights', "vegetation", "rainFall", 'elecciones','codmpio','year')])

cv_model <- cv.glmnet(x = final_lasso, y = y, alpha = 1)
best_lambda <- cv_model$lambda.min
best_model <- glmnet(final_lasso, y, alpha = 1, lambda = best_lambda)
a <- coef(best_model)

variables <- list(a)

saveRDS(variables, file = "Data/Merge/output/lasso_variables.rds")

