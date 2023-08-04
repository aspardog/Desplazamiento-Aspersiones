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
## This version:      March 26th, 2023
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

controles_fe_pop <- c('night_lights', 'vegetation',
                      "rainFall",
                      'ruv_abandono_despojo_pop','ruv_combates_pop',
                      'cnmh_minas_pop', 'cnmh_reclutamiento_pop', 'cnmh_desaparicion_pop')

controles_fe_3month <- c('night_lights', 
                         "rainFall","vegetation", 
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
##   2. Efectos heterogeneos                                                                ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2.1. Aspersión                                                                         ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IVASP <- merge_data.df %>%
  filter(spraying > 0) %>%
  group_by(codmpio) %>%
  summarise(asp_total = mean(spraying_norm, na.rm = T)) %>%
  ungroup() %>%
  mutate(codmpio = as.character(codmpio))

cuartiles <- quantile(IVASP$asp_total, probs = seq(0, 1, 0.25))
IVASP$cuartil <- cut(IVASP$asp_total, cuartiles, labels = FALSE)
IVASP$cuartil <- as.factor(IVASP$cuartil)

IVASP.df <- merge_data.df %>%
  mutate(codmpio = as.character(codmpio)) %>%
  left_join(IVASP, by = c("codmpio")) %>%
  mutate(des = lag1_ruv_desplazamiento_forzado_pop,
         desp2 = lag2_ruv_desplazamiento_forzado_pop,
         desp3 = lag3_ruv_desplazamiento_forzado_pop) 

IVASP.reg1 <- plm(data= IVASP.df,
                 formula =   
                   paste0("des ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                          paste(controles_fe_pop, collapse = "+"),"|",
                          paste("windSpeedRMBOS + cuartil + windSpeedRMBOS*cuartil +"),
                          paste(controles_fe_pop, collapse = "+")) ,
                 effect = "individual", 
                 model = "within", 
                 index=c("date", "codmpio", "query"))
summary(IVASP.reg1)
IVASPHE1 <- coeftest(IVASP.reg1, vcov=vcovHC(IVASP.reg1, type="sss", cluster="group")) 
IVASPHEIC1 <- confint(IVASPHE1, level = 0.8)


IVASP.reg2 <- plm(data= IVASP.df,
                 formula =   
                   paste0("desp2 ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                          paste(controles_fe_pop, collapse = "+"),"|",
                          paste("windSpeedRMBOS + cuartil + windSpeedRMBOS*cuartil +"),
                          paste(controles_fe_pop, collapse = "+")) ,
                 effect = "individual", 
                 model = "within", 
                 index=c("date", "codmpio", "query"))
summary(IVASP.reg2)
IVASPHE2 <- coeftest(IVASP.reg2, vcov=vcovHC(IVASP.reg2, type="sss", cluster="group")) 
IVASPHEIC2 <- confint(IVASPHE2, level = 0.8)

IVASP.reg3 <- plm(data= IVASP.df,
                 formula =   
                   paste0("desp3 ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                          paste(controles_fe_pop, collapse = "+"),"|",
                          paste("windSpeedRMBOS + cuartil + windSpeedRMBOS*cuartil +"),
                          paste(controles_fe_pop, collapse = "+")) ,
                 effect = "individual", 
                 model = "within", 
                 index=c("date", "codmpio", "query"))
summary(IVASP.reg3)
IVASPHE3 <- coeftest(IVASP.reg3, vcov=vcovHC(IVASP.reg3, type="sss", cluster="group")) 
IVASPHEIC3 <- confint(IVASPHE3, level = 0.8)

stargazer(
  IVASP.reg1, IVASP.reg2, IVASP.reg3, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = c("spraying_norm","cuartil"),
  se = list(IVASPHE1[, 2], IVASPHE2[, 2], IVASPHE3[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si"),
                   c("Controles", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "5pt",
  keep.stat = c("rsq", "n"), 
  p = list(IVASPHE1[, 4], IVASPHE2[, 4], IVASPHE3[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2.2. Coca                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVASP <- merge_data.df %>%
  filter(cultivos > 0) %>%
  group_by(codmpio) %>%
  mutate(cultivos_norm = cultivos) %>%
  summarise(asp_total = mean(cultivos_norm, na.rm = T)) %>%
  ungroup()

cuartiles <- quantile(IVASP$asp_total, probs = seq(0, 1, 0.25))

IVASP$cuartil <- cut(IVASP$asp_total, cuartiles, labels = FALSE)
IVASP$cuartil <- as.factor(IVASP$cuartil)

IVASP.df <- merge_data.df %>%
  left_join(IVASP, by = "codmpio") %>%
  mutate(des = lag1_ruv_desplazamiento_forzado_pop,
         desp2 = lag2_ruv_desplazamiento_forzado_pop,
         desp3 = lag3_ruv_desplazamiento_forzado_pop) 

IVASP.reg1 <- plm(data= IVASP.df,
                  formula =   
                    paste0("des ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS + cuartil + spraying_norm*cuartil +"),
                           paste(controles_fe_pop, collapse = "+")) ,
                  effect = "individual", 
                  model = "within", 
                  index=c("date", "codmpio", "query"))
summary(IVASP.reg1)
IVASPHE1 <- coeftest(IVASP.reg1, vcov=vcovHC(IVASP.reg1, type="sss", cluster="group")) 
IVASPHEIC1 <- confint(IVASPHE1, level = 0.8)

IVASP.reg2 <- plm(data= IVASP.df,
                  formula =   
                    paste0("desp2 ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS + cuartil + spraying_norm*cuartil+"),
                           paste(controles_fe_pop, collapse = "+")) ,
                  effect = "individual", 
                  model = "within", 
                  index=c("date", "codmpio", "query"))
summary(IVASP.reg2)
IVASPHE2 <- coeftest(IVASP.reg2, vcov=vcovHC(IVASP.reg2, type="sss", cluster="group")) 
IVASPHEIC2 <- confint(IVASPHE2, level = 0.8)

IVASP.reg3 <- plm(data= IVASP.df,
                  formula =   
                    paste0("desp3 ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS + cuartil + spraying_norm*cuartil +"),
                           paste(controles_fe_pop, collapse = "+")) ,
                  effect = "individual", 
                  model = "within", 
                  index=c("date", "codmpio", "query"))
summary(IVASP.reg3)
IVASPHE3 <- coeftest(IVASP.reg3, vcov=vcovHC(IVASP.reg3, type="sss", cluster="group")) 
IVASPHEIC3 <- confint(IVASPHE3, level = 0.8)

stargazer(
  IVASP.reg1, IVASP.reg2, IVASP.reg3, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = c("spraying_norm","cuartil"),
  se = list(IVASPHE1[, 2], IVASPHE2[, 2], IVASPHE3[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si"),
                   c("Controles", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "5pt",
  keep.stat = c("rsq", "n"), 
  p = list(IVASPHE1[, 4], IVASPHE2[, 4], IVASPHE3[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2.3. Coca - Aspersión                                                                ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVEF <- merge_data.df %>%
  filter(spraying > 0) %>%
  mutate(cultivos_norm = (cultivos/mpio_area)*100) %>%
  group_by(year) %>%
  mutate(meanYearAspersion = mean(spraying_norm, na.rm = T),
         meanYearCultivos  = mean(cultivos_norm, na.rm = T)) %>%
  group_by(codmpio, year) %>%
  mutate(meanAspersion = mean(spraying_norm, na.rm = T),
         meanCultivos  = mean(cultivos_norm, na.rm = T)) %>%
  ungroup()  %>%
  group_by(codmpio, year) %>%
  summarise(hetEffects = if_else(meanAspersion > meanYearAspersion & meanCultivos > meanYearCultivos, "High-High",
                              if_else(meanAspersion > meanYearAspersion & meanCultivos < meanYearCultivos, "High-Low",
                                      if_else(meanAspersion < meanYearAspersion & meanCultivos > meanYearCultivos, "Low-High",
                                              if_else(meanAspersion < meanYearAspersion & meanCultivos < meanYearCultivos, "Low-Low", NA_character_))))) %>%
  distinct()
  
IVHF <- merge_data.df %>%
  left_join(IVEF, by = c("codmpio", "year")) %>%
  mutate(desp1 = lag1_ruv_desplazamiento_forzado_pop,
         desp2 = lag2_ruv_desplazamiento_forzado_pop,
         desp3 = lag3_ruv_desplazamiento_forzado_pop) %>%
  drop_na(hetEffects)

IVHReg1 <- plm(data= IVHF,
              formula =   
                paste0("desp1 ~ spraying_norm + hetEffects + spraying_norm*hetEffects +",
                       paste(controles_fe_pop, collapse = "+"),"|",
                       paste("windSpeedRMBOS + hetEffects + windSpeedRMBOS*hetEffects +"),
                       paste(controles_fe_pop, collapse = "+")),
              effect = "individual", 
              model = "within", 
              index=c("codmpio", "date", "query"))
summary(IVHReg1)
IVHRegHE1 <- coeftest(IVHReg1, vcov=vcovHC(IVHReg1, type="sss", cluster="group")) 
IVHRegHEIC1 <- confint(IVHRegHE1, level = 0.8)

IVHReg2 <- plm(data= IVHF,
              formula =   
                paste0("desp2 ~ spraying_norm + hetEffects + spraying_norm*hetEffects +",
                       paste(controles_fe_pop, collapse = "+"),"|",
                       paste("windSpeedRMBOS + hetEffects + windSpeedRMBOS*hetEffects +"),
                       paste(controles_fe_pop, collapse = "+")),
              effect = "individual", 
              model = "within", 
              index=c("codmpio", "date", "query"))
summary(IVHReg2)
IVHRegHE2 <- coeftest(IVHReg2, vcov=vcovHC(IVHReg2, type="sss", cluster="group")) 
IVHRegHEIC2 <- confint(IVHRegHE2, level = 0.8)


IVHReg3 <- plm(data= IVHF,
               formula =   
                 paste0("desp3 ~ spraying_norm + hetEffects + spraying_norm*hetEffects +",
                        paste(controles_fe_pop, collapse = "+"),"|",
                        paste("windSpeedRMBOS + hetEffects + windSpeedRMBOS*hetEffects +"),
                        paste(controles_fe_pop, collapse = "+")),
               effect = "individual", 
               model = "within", 
               index=c("codmpio", "date", "query"))
summary(IVHReg3)
IVHRegHE3 <- coeftest(IVHReg3, vcov=vcovHC(IVHReg3, type="sss", cluster="group")) 
IVHRegHEIC3 <- confint(IVHRegHE3, level = 0.8)

