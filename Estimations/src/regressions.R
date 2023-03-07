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
## This version:      February 2nd, 2023
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
p_load(tidyverse, sandwich, lmtest, ivreg, corrr)

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
##   1. Correlations                                                                                        ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
controles <- c("ruv_amenaza", "ruv_desaparicion_forzada", "ruv_combates", "ruv_secuestro",
               "ruv_minas", "ruv_tortura", "ruv_homicidio", "ruv_reclutamiento_menores", "cnmh_violencia_sexual",
               "ruv_abandono_despojo", "cnmh_masacres", "elecciones", "vegetation",
               "night_lights")
controles_fe <- c("ruv_amenaza", "ruv_desaparicion_forzada", "ruv_combates", "ruv_secuestro",
               "ruv_minas", "ruv_tortura", "ruv_reclutamiento_menores", "cnmh_violencia_sexual",
               "ruv_abandono_despojo", "cnmh_masacres", "elecciones", "vegetation",
               "night_lights", "codmpio", "year")

correlaciones.df <- merge_data.df %>%
  select(starts_with(controles), spraying, night_lights, windIV20MERRA, windSpeedFLDAS)

correlaciones <- correlate(x = correlaciones.df[,-1]) 

# Homicidio (RUV), Reclutamiento (RUV), Amenazas (RUV), Desaparición (RUV) and Despojo (RUV) are the variables with the highest correlation

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Fixed Effects Regression                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
lm.reg <- lm(data = merge_data.df, 
             formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying"))
summary(lm.reg)

fe.reg <- lm(data = merge_data.df, 
             formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying +", paste(controles, collapse = "+")))
summary(fe.reg)

final.reg <- lm(data = merge_data.df, 
          formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying +", paste(controles_fe, collapse = "+")))
summary(final.reg)

m1coeffs_std <- data.frame(summary(fe.reg)$coefficients)
coi_indices <- which(!startsWith(row.names(m1coeffs_std), 'codmpio'))
m1coeffs_std[coi_indices,]

m1coeffs_cl <- coeftest(fe.reg, vcov = vcovCL, cluster = ~codmpio)
m1coeffs_cl[coi_indices,]
(m1cis <- coefci(fe.reg, parm = coi_indices, vcov = vcovCL,
                 cluster = ~codmpio, level = 0.90))
# Clusterizar errores mostrarlo

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Instrumental Variables Regression                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

windIV.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ spraying | windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(windIV.reg)

m1coeffs_std <- data.frame(summary(windIV.reg)$coefficients)
coi_indices <- which(!startsWith(row.names(m1coeffs_std), 'codmpio'))
m1coeffs_std[coi_indices,]

m1coeffs_cl <- coeftest(windIV.reg, vcov = vcovCL, cluster = ~codmpio)
m1coeffs_cl[coi_indices,]
(m1cis <- coefci(windIV.reg, parm = coi_indices, vcov = vcovCL,
                 cluster = ~codmpio, level = 0.90))

# Robustness Checks

# First Stage

firstStage.reg <- lm (data = merge_data.df,
                      formula = paste0("spraying ~ windSpeedFLDAS +", paste(controles_fe, collapse = "+")))
summary(firstStage.reg)

# Añadir acumulado
