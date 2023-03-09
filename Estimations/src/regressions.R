## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     February 2nd, 2023
##
## This version:      March 7th, 2023
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##  Model Estimation                                                                                 ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    Set - UP                                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(pacman)
p_load(tidyverse, sandwich, lmtest, ivreg, corrr, modelsummary, kableExtra, gt, tibble, stargazer)

cluster_errors.fn <- function(reg) {
  m1coeffs_std <- data.frame(summary(reg)$coefficients)
  coi_indices <- which(!startsWith(row.names(m1coeffs_std), 'codmpio'))
  m1coeffs_std[coi_indices,]
  
  m1coeffs_cl <- coeftest(reg, vcov = vcovCL, cluster = ~codmpio)
  a <- m1coeffs_cl[coi_indices,]
  estimation <- data.frame(a)
  
  estimation <- rownames_to_column(estimation, "term")
  
  estimation <- estimation %>%
    rename(estimate = Estimate, std.error = Std..Error, p.value = Pr...t..)
  
  return(estimation)
}

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


antimerge_data.df <- readRDS("Data/Merge/output/antimerge_data.rds") %>%
  filter(year < 2013 & year > 2003) %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         month = as.factor(month)) %>%
  drop_na()

lassoVariables <- readRDS("Data/Merge/output/lasso_variables.rds")[[1]]
lassoVariables  # Variables to select the controls



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Correlations                                                                                        ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

controles <- c('ruv_amenaza', 'cnmh_desaparicion', 'ruv_combates',
               'cnmh_d_bienes', 'cnmh_minas', 'ruv_homicidio', 
               'ruv_violencia_sexual', 'ruv_abandono_despojo',
               'cnmh_reclutamiento', 'cnmh_atentado', 'cnmh_ataque_poblacion',
               'night_lights', "vegetation", 'elecciones')
controles_fe <- c('ruv_amenaza', 'cnmh_desaparicion', 'ruv_combates',
                  'cnmh_d_bienes', 'cnmh_minas', 'ruv_homicidio', 
                  'ruv_violencia_sexual', 'ruv_abandono_despojo', "rainFall",
                  'cnmh_reclutamiento', 'cnmh_atentado', 'cnmh_ataque_poblacion',
                  'night_lights', 'elecciones','codmpio','year')


correlaciones.df <- merge_data.df %>%
  dplyr::select(starts_with(controles), spraying, night_lights, windIV20MERRA, windSpeedFLDAS, ruv_desplazamiento_forzado)

correlaciones <- correlate(x = correlaciones.df[,-1]) 

# Homicidio (RUV), Reclutamiento (RUV), Amenazas (RUV), Desaparición (RUV) and Despojo (RUV) are the variables with the highest correlation

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Fixed Effects Regression                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Regresion desplazamiento ~ aspersion
lm.reg <- lm(data = merge_data.df, 
             formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying"))
summary(lm.reg)
cluster_errors.fn(lm.reg)

# Regresion desplazamiento ~ aspersion con controles

fe.reg <- lm(data = merge_data.df, 
             formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying +", paste(controles, collapse = "+")))
summary(fe.reg)
cluster_errors.fn(fe.reg)


# Regresion desplazamiento ~ aspersion con controles y efectos fijos

final.reg <- lm(data = merge_data.df, 
          formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying +", paste(controles_fe, collapse = "+")))
summary(final.reg)
cluster_errors.fn(final.reg)

# m1coeffs_std <- data.frame(summary(final.reg)$coefficients)
# coi_indices <- which(!startsWith(row.names(m1coeffs_std), 'codmpio'))
# m1coeffs_std[coi_indices,]
# 
# m1coeffs_cl <- coeftest(final.reg, vcov = vcovCL, cluster = ~codmpio)
# m1coeffs_cl[coi_indices,]
# (m1cis <- coefci(final.reg , parm = coi_indices, vcov = vcovCL,
#                  cluster = ~codmpio, level = 0.90))
# Clusterizar errores mostrarlo

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   4. Instrumental Variables Regression                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windIV.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ spraying | windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(windIV.reg)
cluster_errors.fn(windIV.reg)

# Efectos heterogeneos 

windIV.dpto.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ spraying*dpto | windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(windIV.dpto.reg)
cluster_errors.fn(windIV.dpto.reg)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5.  Robustness Checks                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# First Stage

modelsFS <- list(
  "(1)"     = lm (data = merge_data.df,
                  formula = paste0("spraying ~ windSpeedFLDAS")),
  "(2)" = lm (data = merge_data.df,
              formula = paste0("spraying ~ windSpeedFLDAS +", paste(controles, collapse = "+"))),
  "(3)"     = lm (data = merge_data.df,
                  formula = paste0("spraying ~ windSpeedFLDAS +", paste(controles_fe, collapse = "+")))
)

coef_rename_ws <- "Velocidad viento"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nobs",      "N",             0,
  "r.squared", "R<sup>2</sup>", 2,)
new_rows <- tibble::tribble(
  ~term,~"(1)",~"(2)",~"(3)",
  "Controles","No","Sí","Sí",
  "Efectos Fijos", "No","No","Sí",)

modelsummary(modelsFS, vcov = ~codmpio, estimate = "{estimate}{stars}", ,
             coef_omit = c(-2), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 1, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Aspersiones aéreas", 
             notes = "Los errores estándar son robustos y están corregidos por clusters de hogar.  *** p<0.01, ** p<0.05, * p<0.1.")

# Restricción de exclusión

modelsRE <- list(
  "Municipios con Aspérsion Aéreas" = lm(merge_data.df,
                                         formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  windSpeedFLDAS +", paste(controles_fe, collapse = "+"))),
  "Municipios sin Aspérsion Aéreas" = lm(antimerge_data.df,
                                         formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  windSpeedFLDAS +", paste(controles_fe, collapse = "+"))))

coef_rename_ws <- "Velocidad viento"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nobs",      "N",             0,
  "r.squared", "R<sup>2</sup>", 2,)
new_rows <- tibble::tribble(
  ~term,~"Municipios con Aspérsion Aéreas",~"Municipios sin Aspérsion Aéreas",
  "Controles","Sí","Sí",
  "Efectos Fijos", "Sí","Sí")

modelsummary(modelsRE, vcov = ~codmpio, estimate = "{estimate}{stars}", ,
             coef_omit = c(-2), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 3, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Variable dependiente: *Desplazamiento forzado*", 
             notes = "*** p<0.01, ** p<0.05, * p<0.1.
             Los errores estándar son robustos y están corregidos por clusters de hogar.")

# linear model
firstStage.reg <- lm (data = merge_data.df,
                            formula = paste0("spraying ~ windSpeedFLDAS"))

statisticsFS <- summary(firstStage.reg)
glFS <- data.frame(N.obs = c(nrow(merge_data.df)), 
                 "Adjusted R-squared" = c(round(statisticsFS$adj.r.squared,2)))

firstStage <- cluster_errors.fn(firstStage.reg) %>%
  filter(term %in% "windSpeedFLDAS") %>%
  mutate(term = if_else(term %in% "windSpeedFLDAS", "Velocidad del viento", term)) %>%
  cbind(glFS)

# model with controls

firstStageControl.reg <- lm (data = merge_data.df,
                      formula = paste0("spraying ~ windSpeedFLDAS +", paste(controles, collapse = "+")))

statisticsFSC <- summary(firstStageControl.reg)
glFSC <- data.frame(N.obs = c(nrow(merge_data.df)), 
                   "R.squared-adjusted" = c(round(statisticsFSC$adj.r.squared,2)))

firstStageControl <- cluster_errors.fn(firstStageControl.reg) %>%
  filter(term %in% "windSpeedFLDAS") %>%
  mutate(term = if_else(term %in% "windSpeedFLDAS", "Velocidad del viento", term)) %>%
  cbind(glFSC)

# model with controls and FE

firstStage.reg.final <- lm (data = merge_data.df,
                      formula = paste0("spraying ~ windSpeedFLDAS +", paste(controles_fe, collapse = "+")))

statisticsFSCF <- summary(firstStage.reg.final)
glFSCF <- data.frame(N.obs = c(nrow(merge_data.df)), 
                   "R.squared-adjusted" = c(round(statisticsFSCF$adj.r.squared,2)))

firstStageFinal <- cluster_errors.fn(firstStage.reg.final) %>%
  filter(term %in% "windSpeedFLDAS") %>%
  mutate(term = if_else(term %in% "windSpeedFLDAS", "Velocidad del viento", term))

# Restricción de exclusión

secondCheck <- lm(merge_data.df,
                  formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  windSpeedFLDAS +", paste(controles_fe, collapse = "+")))
summary(secondCheck)
secondCheck <- cluster_errors.fn(secondCheck)

secondCheck_anti <- lm(antimerge_data.df, 
                  formula = paste0("log(ruv_desplazamiento_forzado + min(ruv_desplazamiento_forzado[ruv_desplazamiento_forzado>0])/2) ~  windSpeedFLDAS +", paste(controles_fe, collapse = "+")))
summary(secondCheck_anti)
secondCheck_anti <- cluster_errors.fn(secondCheck_anti)

# Añadir acumulado

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   6.  Tables                                                                  ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# First Stage

class(firstStage) <- "lm"
class(firstStageControl) <- "lm"

mod <- list(firstStage,firstStageControl)

a <- modelsummary(models = mod)

mod1 <- list(
  tidy = firstStageControl,
  glance = glFSC
)

class(mod1) <- "modelsummary_list"

a <- modelsummary(models = c(mod,mod1), estimate = c("{estimate}{stars}"), fmt = 1, output = "markdown")



### Valores atipicos, por encima de dos desviaciones estandar
