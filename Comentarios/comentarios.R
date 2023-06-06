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
## This version:      June 3rd, 2023
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
p_load(tidyverse, sandwich, lmtest, ivreg, corrr, modelsummary, kableExtra, gt, 
       tibble, stargazer, plm,  ggpubr, showtext, patchwork, ggh4x, knitr, flextable)

controles_fe_pop <- c('night_lights', 
                      "rainFall","vegetation_norm", 'windIV10RMBOS',
                      'ruv_abandono_despojo_pop','ruv_combates_pop', 'ruv_homicidio_pop',
                      'cnmh_minas_pop', 'cnmh_reclutamiento_pop','cnmh_desaparicion_pop')

controles_fe_3month <- c('night_lights', 
                         "rainFall","vegetation",'windIV10RMBOS',
                         'sum_combates_pop', 'sum_despojo_pop', 'sum_minas_pop', 'sum_reclutamiento_pop', 
                         'sum_homicidio_pop', 'sum_desaparicion_pop')


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Download                                                                                            ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") 


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Aggregate effects                                                                                ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

agg.df <- merge_data.df %>%
  mutate(het = if_else(spraying_norm > 0 & lag1_spraying > 0, 1, 0))  %>%
  mutate(des = sum_desplazamiento_pop) %>%
  mutate(sum_spraying_new = if_else(het == 1, lag1_spraying + spraying_norm, spraying_norm))

agg.reg <- plm(data = agg.df,
               formula =   
                 paste0("des ~ sum_spraying_new +",
                        paste(controles_fe_3month, collapse = "+"),"|",
                        paste("windSpeedRMBOS +"),
                        paste(controles_fe_3month, collapse = "+")) ,
               effect = "twoways", 
               model = "within", 
               index=c("date", "codmpio"))

AGG <- coeftest(agg.reg, vcov=vcovHC(agg.reg, type="sss", cluster="group")) 
AGG

stargazer(
  agg.reg,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "sum_spraying",
  se = list(AGG[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si", "Si", "Si"),
                   c("Controles", "Si", "No", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = list(AGG[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Heterogeneous effects                                                                                ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

het.df <- merge_data.df %>%
  mutate(het = if_else(spraying_norm > 0 & (lag3_spraying > 0 | lag2_spraying > 0 | lag1_spraying > 0), 1, 0))  %>%
  mutate(des = sum_desplazamiento_pop) 

het.reg <- plm(data = het.df,
               formula =   
                 paste0("des ~ sum_spraying + het + sum_spraying*het +",
                        paste(controles_fe_3month, collapse = "+"),"|",
                        paste("windSpeedRMBOS + het + windSpeedRMBOS*het +"),
                        paste(controles_fe_3month, collapse = "+")) ,
               effect = "individual", 
               model = "within", 
               index=c("date", "codmpio", "query"))

HET <- coeftest(het.reg, vcov=vcovHC(het.reg, type="sss", cluster="group")) 
HET

stargazer(
  het.reg,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "sum_spraying",
  se = list(HET[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si", "Si", "Si"),
                   c("Controles", "Si", "No", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = list(HET[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)
