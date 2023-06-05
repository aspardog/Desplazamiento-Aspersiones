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

controles_fe_pop <- c('night_lights', 
                      "rainFall","vegetation", 
                      'ruv_abandono_despojo_pop','ruv_combates_pop', 'ruv_homicidio_pop',
                      'cnmh_minas_pop', 'cnmh_reclutamiento_pop','cnmh_desaparicion_pop')

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

antimerge_data.df <- readRDS("Data/Merge/output/antimerge_data.rds") 

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. First Stage                                                               ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

firstStage <- plm(data = merge_data.df,
                  formula = paste0("spraying_norm ~  windSpeedRMBOS +", paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", model = "within", index=c("codmpio", "date", "query"))
summary(firstStage)
stdFirstStage <- coeftest(firstStage, vcov=vcovHC(firstStage, type="sss", cluster="group"))  

stargazer(
  firstStage, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Aspersiones aéreas",
  keep = "windSpeedRMBOS",
  se = list(stdFirstStage[, 2]), 
  title = "Resultados primera etapa", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = TRUE, 
  covariate.labels = c("Velocidad del Viento"), 
  omit = "Constant",
  add.lines = list(c("F estadístico", "21,208"),
                   c("F estadístico efectivo", "17,952"),
                   c("Efectos Fijos", "Sí"),
                   c("Controles", "Sí")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = stdFirstStage[, 4],
  decimal.mark = ",",
  notes.align = "l"
)
# Nota: Los errores fueron clusterizado a nivel municipal. Los valores dentro de los parentesis representan la desviación estandar. Paralelamente se aplicaron efectos fijos por municipio, año-mes y núcleo. Las variables de control asociadas a la violencia se tomaron en tasas por 100 habitantes, estas son: desaparición forzada, reclutamiento de menores, minas, combates y despojo. Las variables de control geograficas son: choques de viento, indice de vegetación y niveles de lluvia. La variable de control asociada al desarrollo economico es la intensidad de luminosidad del municipio. Además, los niveles de signifancia se ven representados de la siguiente manera: $^{*}$p$<$0,1; $^{**}$p$<$0,05; $^{***}$p$<$0,01  

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Restricción de exclusión                                                                ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

restExcl1 <- plm(data = antimerge_data.df,
                 formula = paste0("lag1_ruv_desplazamiento_forzado_pop ~  windSpeedRMBOS +", paste(controles_fe_pop, collapse = "+")),
                 effect = "twoways", model = "within", index=c("codmpio", "date", "query"))

summary(restExcl1)

restExcl1std <- coeftest(restExcl1, vcov=vcovHC(restExcl1, type="sss", cluster="group"))  


restExcl2 <- plm(data = merge_data.df,
                 formula = paste0("lag1_ruv_desplazamiento_forzado_pop ~  windSpeedRMBOS +", paste(controles_fe_pop, collapse = "+")),
                 effect = "twoways", model = "within", index=c("codmpio", "date", "query"))
summary(restExcl2)
restExcl2std <- coeftest(restExcl2, vcov=vcovHC(restExcl2, type="sss", cluster="group"))  

stargazer(
  restExcl1, restExcl2,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "windSpeedRMBOS",
  se = list(restExcl1std[, 2], restExcl2std[, 2]), 
  title = "Resultados primera etapa", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = TRUE, 
  covariate.labels = c("Velocidad del Viento"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Sí"),
                   c("Controles", "Sí")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = list(restExcl1std[, 4], restExcl2std[,4]),
  decimal.mark = ",",
  notes.align = "l"
)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   4. Efectos dinámicos                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# t+1

IV1MonthFE <- plm(data= merge_data.df,
                  formula =  
                    paste0("lag1_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date", "query"))
summary(IV1MonthFE)
IVFER <- coeftest(IV1MonthFE, vcov=vcovHC(IV1MonthFE, type="sss", cluster="group"))  
IVFECI <- confint(IVFER, level = 0.9)

IVFET1 <- data.frame("Estimate" = IVFER[[1]],
                     "std" = IVFER[[1,2]],
                     "p.value" = IVFER[[1,4]],
                     "label" = "t+1",
                     "lower" = IVFECI[[1]],
                     "upper" = IVFECI[[1,2]])

# t+2

dataLag2.df <- merge_data.df %>%
  drop_na(lag2_ruv_desplazamiento_forzado)

IV2MonthFE <- plm(data= dataLag2.df,
                  formula =  
                    paste0("lag2_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date", "query"))
summary(IV2MonthFE)
IVFE2R <- coeftest(IV2MonthFE, vcov=vcovHC(IV2MonthFE, type="sss", cluster="group"))  
IVFE2CI <- confint(IVFE2R, level = 0.9)

IVFET2 <- data.frame("Estimate" = IVFE2R[[1]],
                     "std" = IVFE2R[[1,2]],
                     "p.value" = IVFE2R[[1,4]],
                     "label" = "t+2",
                     "lower" = IVFE2CI[[1]],
                     "upper" = IVFE2CI[[1,2]]) 
# t+3

dataLag3.df <- merge_data.df %>%
  drop_na(lag3_ruv_desplazamiento_forzado)

IV3MonthFE <- plm(data= dataLag3.df,
                  formula =  
                    paste0("lag3_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date", "query"))
summary(IV3MonthFE)
IVFE3R <- coeftest(IV3MonthFE, vcov=vcovHC(IV3MonthFE, type="sss", cluster="group"), parm = ci) 
IVFE3CI <- confint(IVFE3R, level = 0.9)

IVFET3 <- data.frame("Estimate" = IVFE3R[[1]],
                     "std" = IVFE3R[[1,2]],
                     "p.value" = IVFE3R[[1,4]],
                     "label" = "t+3",
                     "lower" = IVFE3CI[[1]],
                     "upper" = IVFE3CI[[1,2]])

stargazer(
  IV1MonthFE, IV2MonthFE, IV3MonthFE,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list( IVFER[, 2], IVFE2R[, 2], IVFE3R[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si", "Si", "Si"),
                   c("Controles", "Si", "Si", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "5pt",
  keep.stat = c("rsq", "n"), 
  p = list(IVFER[, 4], IVFE2R[, 4], IVFE3R[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

# t+4

dataLag4.df <- merge_data.df %>%
  drop_na(lag4_ruv_desplazamiento_forzado)

IV4MonthFE <- plm(data= dataLag4.df,
                  formula =  
                    paste0("lag4_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date", "query"))
summary(IV4MonthFE)
IVFE4R <- coeftest(IV4MonthFE, vcov=vcovHC(IV4MonthFE, type="sss", cluster="group"))  
IVFE4CI <- confint(IVFE4R, level = 0.9)

IVFET4 <- data.frame("Estimate" = IVFE4R[[1]],
                     "std" = IVFE4R[[1,2]],
                     "p.value" = IVFE4R[[1,4]],
                     "label" = "t+4",
                     "lower" = IVFE4CI[[1]],
                     "upper" = IVFE4CI[[1,2]])
#t

dataT.df <- merge_data.df 

IVTMonthFE <- plm(data= dataT.df,
                  formula =  
                    paste0("ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date", "query"))
summary(IVTMonthFE)
IVFETR <- coeftest(IVTMonthFE, vcov=vcovHC(IVTMonthFE, type="sss", cluster="group"))  
IVFETCI <- confint(IVFETR, level = 0.9)

IVFET  <- data.frame("Estimate" = IVFETR[[1]],
                     "std" = IVFETR[[1,2]],
                     "p.value" = IVFETR[[1,4]],
                     "label" = "t",
                     "lower" = IVFETCI[[1]],
                     "upper" = IVFETCI[[1,2]])
#t-1

dataLagM1.df <- merge_data.df %>%
  arrange(order_value) %>%
  group_by(codmpio) %>%
  mutate(lagM1_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado_pop, n=1)) %>%
  ungroup() %>%
  drop_na(lagM1_ruv_desplazamiento_forzado)

IVM1MonthFE <- plm(data= dataLagM1.df,
                   formula =  
                     paste0("lagM1_ruv_desplazamiento_forzado ~ spraying_norm +",
                            paste(controles_fe_pop, collapse = "+"),"|",
                            paste("windSpeedRMBOS +"),
                            paste(controles_fe_pop, collapse = "+")),
                   effect = "twoways", 
                   model = "within", 
                   index=c("codmpio", "date", "query"))
summary(IVM1MonthFE)
IVFEM1R <- coeftest(IVM1MonthFE, vcov=vcovHC(IVM1MonthFE, type="sss", cluster="group"))  
IVFEM1CI <- confint(IVFEM1R, level = 0.9)

IVFEM1  <- data.frame("Estimate" = IVFEM1R[[1]],
                      "std" = IVFEM1R[[1,2]],
                      "p.value" = IVFEM1R[[1,4]],
                      "label" = "t-1",
                      "lower" = IVFEM1CI[[1]],
                      "upper" = IVFEM1CI[[1,2]])
#t-2

dataLagM2.df <- merge_data.df %>%
  arrange(order_value) %>%
  group_by(codmpio) %>%
  mutate(lagM2_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado_pop, n=2)) %>%
  ungroup() %>%
  drop_na(lagM2_ruv_desplazamiento_forzado)

IVM2MonthFE <- plm(data= dataLagM2.df,
                   formula =  
                     paste0("lagM2_ruv_desplazamiento_forzado ~ spraying +",
                            paste(controles_fe_pop, collapse = "+"),"|",
                            paste("windSpeedRMBOS +"),
                            paste(controles_fe_pop, collapse = "+")),
                   effect = "twoways", 
                   model = "within", 
                   index=c("codmpio", "date", "query"))
summary(IVM2MonthFE)
IVFEM2R <- coeftest(IVM2MonthFE, vcov=vcovHC(IVM2MonthFE, type="sss", cluster="group"))  
IVFEM2CI <- confint(IVFEM2R, level = 0.9)

IVFEM2  <- data.frame("Estimate" = IVFEM2R[[1]],
                      "std" = IVFEM2R[[1,2]],
                      "p.value" = IVFEM2R[[1,4]],
                      "label" = "t-2",
                      "lower" = IVFEM2CI[[1]],
                      "upper" = IVFEM2CI[[1,2]])

# Create the data frame

dinamicsEffects <- rbind(IVFET1, IVFET2, IVFET3, IVFET4, IVFET, IVFEM1, IVFEM2) %>%
  mutate(estimate = round(Estimate,3),
         lower = round(lower,3),
         upper = round(upper,3),
         variable = "Aspersiones Aéreas") %>%
  mutate(order_value = case_when(
    label == "t-2" ~ 1,
    label == "t-1" ~ 2,
    label == "t" ~ 3,
    label == "t+1" ~ 4,
    label == "t+2" ~ 5,
    label == "t+3" ~ 6,
    label == "t+4" ~ 7,
    
  ))

# Chart

dinamicPlot <- ggplot(dinamicsEffects, aes(x = reorder(label, order_value), y = estimate)) +
  geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
  geom_errorbar(aes(x = reorder(label, order_value),  ymin = lower, ymax = upper),
                lwd = 0.5, position = position_dodge(width = .7), 
                stat = "identity", color = "black")+
  geom_point(aes(x = reorder(label, order_value), y = estimate), 
             size = 3.5, position = position_dodge(width = 1), color = "black") +
  geom_point(aes(x = reorder(label, order_value), y = estimate), 
             size = 2.5, position = position_dodge(width = 1), color = "white") +
  scale_y_continuous(limits = c(-3, 3),
                     breaks = seq(-3, 3, by = 1),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-3","-2","-1", "0", "1", "2", "3")) +
  labs(x = "Mes",
       y = "Efecto aspersiones sobre desplazamiento") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(legend.position = "none",
        panel.background   = element_blank(),
        panel.grid.major.x = element_line(colour = "#d1cfd1", 
                                          size = 0.5, linetype = "dashed"),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 12, margin   = margin(0, 20, 0, 10), vjust = 1),
        axis.title.x = element_text(size = 12, margin   = margin(20, 0, 10, 0), vjust = 0),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.minor.y = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1),
        axis.line.x.bottom = element_line(linetype = "solid", size = 1));dinamicPlot

ggsave(dinamicPlot, filename = "Visualizations/output/DinamicEffectsNorm.png", dpi = 320, width = 10, height = 7.5)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Efectos heterogeneos                                                                ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVEF <- merge_data.df %>%
  mutate(cultivos_norm = (cultivos/mpio_area)) %>%
  mutate(meanYearAspersion = mean(spraying_norm, na.rm = T),
         meanYearCultivos  = mean(cultivos_norm, na.rm = T)) %>%
  group_by(codmpio) %>%
  mutate(meanAspersion = mean(spraying_norm, na.rm = T),
         meanCultivos  = mean(cultivos_norm, na.rm = T)) %>%
  ungroup() 

IVHF <- IVEF %>%
  mutate(hetEffects = if_else(meanAspersion > meanYearAspersion & meanCultivos > meanYearCultivos, "1High-High",
                              if_else(meanAspersion > meanYearAspersion & meanCultivos < meanYearCultivos, "High-Low",
                                      if_else(meanAspersion < meanYearAspersion & meanCultivos > meanYearCultivos, "Low-High",
                                              if_else(meanAspersion < meanYearAspersion & meanCultivos < meanYearCultivos, "Low-Low", NA_character_))))) %>%
  mutate(quintilesAspersion = as.factor(ntile(spraying, 5)),
         hetEffects = as.factor(hetEffects),
         desp = lag1_ruv_desplazamiento_forzado_pop/100)

IVHReg <- plm(data= IVHF,
              formula =   
                paste0("desp ~ spraying_norm + hetEffects + spraying_norm*hetEffects +",
                       paste(controles_fe_pop, collapse = "+"),"|",
                       paste("windSpeedRMBOS + hetEffects + windSpeedRMBOS*hetEffects +"),
                       paste(controles_fe_pop, collapse = "+")),
              effect = "individual", 
              model = "within", 
              index=c("codmpio", "date", "query"))
summary(IVHReg)
IVHRegHE <- coeftest(IVHReg, vcov=vcovHC(IVHReg, type="sss", cluster="group")) 
IVHRegHEIC <- confint(IVHRegHE, level = 0.8)

stargazer(
  IVHReg, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(IVHRegHE[, 2]), 
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
  p = list(IVHRegHE[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

estimate <- as.data.frame(IVHRegHE[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:hetEffectsLow-Low", "spraying_norm:hetEffectsHigh-Low", "spraying_norm:hetEffectsLow-High"))
CI <- as.data.frame(IVHRegHEIC[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:hetEffectsLow-Low", "spraying_norm:hetEffectsHigh-Low", "spraying_norm:hetEffectsLow-High")) %>%
  rename("lower" = "10 %", "upper" = "90 %")

data2plot <- estimate %>%
  inner_join(CI, by = "variable") %>%
  mutate(
    order_value =
      case_when(
        variable == "spraying_norm:hetEffectsLow-High" ~ 1,
        variable == "spraying_norm:hetEffectsHigh-Low" ~ 2,
        variable == "spraying_norm:hetEffectsLow-Low" ~ 3,
      ),
    variable = 
      case_when(
        variable == "spraying_norm:hetEffectsLow-High" ~ "Baja Aspersión - Altos Cultivos",
        variable == "spraying_norm:hetEffectsHigh-Low" ~ "Alta Aspersión - Bajos Cultivos",
        variable == "spraying_norm:hetEffectsLow-Low" ~ "Baja Aspersión - Bajos Cultivos",
      )
  )

HECocaPlot <- ggplot(data2plot, aes(x = reorder(variable, order_value), y = Estimate)) +
  geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
  geom_linerange(aes(x = reorder(variable, order_value),  ymin = lower, ymax = upper),
                 lwd = 0.5, position = position_dodge(width = .7), 
                 stat = "identity", color = "#003b8a")+
  geom_point(aes(x =reorder(variable, order_value), y = Estimate), 
             size = 3.5, position = position_dodge(width = 1), color = "#003b8a") +
  geom_point(aes(x = reorder(variable, order_value), y = Estimate), 
             size = 2.5, position = position_dodge(width = 1), color = "white") +
  scale_y_continuous(limits = c(-0.25, 0.25),
                     breaks = seq(-0.25, 0.25, by = 0.125),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-0.25", "-0.125", "0", "0.125","0.25")) +
  coord_flip() +
  labs(x = "Nivel de aspersión y cultivos de coca",
       y = "Efecto aspersiones sobre desplazamiento") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(legend.position = "none",
        panel.background   = element_blank(),
        panel.grid.major.x = element_line(colour = "#d1cfd1", 
                                          size = 0.5, linetype = "dashed"),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 12, margin   = margin(0, 20, 0, 10), vjust = 1),
        axis.title.x = element_text(size = 12, margin   = margin(20, 0, 10, 0), vjust = 0),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.minor.y = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1),
        axis.line.x.bottom = element_line(linetype = "solid", size = 1));HECocaPlot 
ggsave(HECocaPlot, filename = "Visualizations/output/CultivosAspersionHENorm.png", dpi = 320, width = 7.5, height = 7.5)

# Quintiles intensidad


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
  mutate(des = lag1_ruv_desplazamiento_forzado_pop/100) 

IVASP.reg <- plm(data= IVASP.df,
                 formula =   
                   paste0("des ~ spraying_norm + cuartil + spraying_norm*cuartil +",
                          paste(controles_fe_pop, collapse = "+"),"|",
                          paste("windSpeedRMBOS + cuartil + windSpeedRMBOS*cuartil +"),
                          paste(controles_fe_pop, collapse = "+")) ,
                 effect = "individual", 
                 model = "within", 
                 index=c("date", "codmpio", "query"))
summary(IVASP.reg)
IVASPHE <- coeftest(IVASP.reg, vcov=vcovHC(IVASP.reg, type="sss", cluster="group")) 
IVASPHEIC <- confint(IVASPHE, level = 0.8)

stargazer(
  IVASP.reg, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(IVASPHE[, 2]), 
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
  p = list(IVASPHE[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

estimate <- as.data.frame(IVASPHE[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:cuartil2", "spraying_norm:cuartil3", "spraying_norm:cuartil4"))
CI <- as.data.frame(IVASPHEIC[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:cuartil2", "spraying_norm:cuartil3", "spraying_norm:cuartil4")) %>%
  rename("lower" = "10 %", "upper" = "90 %")

data2plot <- estimate %>%
  inner_join(CI, by = "variable") %>%
  mutate(
    order_value =
      case_when(
        variable == "spraying_norm:cuartil2" ~ 1,
        variable == "spraying_norm:cuartil3" ~ 2,
        variable == "spraying_norm:cuartil4" ~ 3
      ),
    variable = 
      case_when(
        variable == "spraying_norm:cuartil2" ~ "Segundo cuartil",
        variable == "spraying_norm:cuartil3" ~ "Tercer cuartil",
        variable == "spraying_norm:cuartil4" ~ "Cuarto cuartil"      
      )
  )

HEPlot <- ggplot(data2plot, aes(x = reorder(variable, order_value), y = Estimate)) +
  geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
  geom_linerange(aes(x = reorder(variable, order_value),  ymin = lower, ymax = upper),
                 lwd = 0.5, position = position_dodge(width = .7), 
                 stat = "identity", color = "#003b8a")+
  geom_point(aes(x =reorder(variable, order_value), y = Estimate), 
             size = 3.5, position = position_dodge(width = 1), color = "#003b8a") +
  geom_point(aes(x = reorder(variable, order_value), y = Estimate), 
             size = 2.5, position = position_dodge(width = 1), color = "white") +
  scale_y_continuous(limits = c(-2, 2),
                     breaks = seq(-2, 2, by = 1),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-2", "-1", "0", "1","2")) +
  coord_flip() +
  labs(x = "Quintil de intensidad de aspersión",
       y = "Efecto aspersiones sobre desplazamiento") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(legend.position = "none",
        panel.background   = element_blank(),
        panel.grid.major.x = element_line(colour = "#d1cfd1", 
                                          size = 0.5, linetype = "dashed"),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 12, margin   = margin(0, 20, 0, 10), vjust = 1),
        axis.title.x = element_text(size = 12, margin   = margin(20, 0, 10, 0), vjust = 0),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.minor.y = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1),
        axis.line.x.bottom = element_line(linetype = "solid", size = 1));HEPlot 
ggsave(HEPlot, filename = "Visualizations/output/AspersionHENorm.png", dpi = 320, width = 7.5, height = 7.5)

# Quintiles intensidad COCA


IVASP <- merge_data.df %>%
  filter(cultivos > 0) %>%
  group_by(codmpio) %>%
  mutate(cultivos_norm = (cultivos/mpio_area)*100) %>%
  summarise(asp_total = mean(cultivos_norm, na.rm = T)) %>%
  ungroup()

quintiles <- quantile(IVASP$asp_total, probs = seq(0, 1, 0.25))
IVASP$quintil <- cut(IVASP$asp_total, quintiles, labels = FALSE)
IVASP$quintil <- as.factor(IVASP$quintil)

IVASP.df <- merge_data.df %>%
  left_join(IVASP, by = "codmpio") %>%
  mutate(des = lag1_ruv_desplazamiento_forzado_pop/100) 

IVASP.reg <- plm(data= IVASP.df,
                 formula =   
                   paste0("des ~ spraying_norm + quintil + spraying_norm*quintil +",
                          paste(controles_fe_pop, collapse = "+"),"|",
                          paste("windSpeedRMBOS + quintil + windSpeedRMBOS*quintil +"),
                          paste(controles_fe_pop, collapse = "+")) ,
                 effect = "individual", 
                 model = "within", 
                 index=c("date", "codmpio", "query"))
summary(IVASP.reg)
IVASPHE <- coeftest(IVASP.reg, vcov=vcovHC(IVASP.reg, type="sss", cluster="group")) 
IVASPHEIC <- confint(IVASPHE, level = 0.8)

stargazer(
  IVASP.reg, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(IVASPHE[, 2]), 
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
  p = list(IVASPHE[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

estimate <- as.data.frame(IVASPHE[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:quintil2", "spraying_norm:quintil3", "spraying_norm:quintil4"))
CI <- as.data.frame(IVASPHEIC[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:quintil2", "spraying_norm:quintil3", "spraying_norm:quintil4")) %>%
  rename("lower" = "10 %", "upper" = "90 %")

data2plot <- estimate %>%
  inner_join(CI, by = "variable") %>%
  mutate(
    order_value =
      case_when(
        variable == "spraying_norm:quintil2" ~ 1,
        variable == "spraying_norm:quintil3" ~ 2,
        variable == "spraying_norm:quintil4" ~ 3,
      ),
    variable = 
      case_when(
        variable == "spraying_norm:quintil2" ~ "Segundo quintil",
        variable == "spraying_norm:quintil3" ~ "Tercer quintil",
        variable == "spraying_norm:quintil4" ~ "Cuarto quintil"
      )
  )

HEPlot <- ggplot(data2plot, aes(x = reorder(variable, order_value), y = Estimate)) +
  geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
  geom_linerange(aes(x = reorder(variable, order_value),  ymin = lower, ymax = upper),
                 lwd = 0.5, position = position_dodge(width = .7), 
                 stat = "identity", color = "#003b8a")+
  geom_point(aes(x =reorder(variable, order_value), y = Estimate), 
             size = 3.5, position = position_dodge(width = 1), color = "#003b8a") +
  geom_point(aes(x = reorder(variable, order_value), y = Estimate), 
             size = 2.5, position = position_dodge(width = 1), color = "white") +
  scale_y_continuous(limits = c(-2, 2),
                     breaks = seq(-2, 2, by = 1),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-2", "-1", "0", "1","2")) +
  coord_flip() +
  labs(x = "Quintil de niveles cultivos de coca",
       y = "Efecto aspersiones sobre desplazamiento") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(legend.position = "none",
        panel.background   = element_blank(),
        panel.grid.major.x = element_line(colour = "#d1cfd1", 
                                          size = 0.5, linetype = "dashed"),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 12, margin   = margin(0, 20, 0, 10), vjust = 1),
        axis.title.x = element_text(size = 12, margin   = margin(20, 0, 10, 0), vjust = 0),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.minor.y = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1),
        axis.line.x.bottom = element_line(linetype = "solid", size = 1));HEPlot 
ggsave(HEPlot, filename = "Visualizations/output/CultivosHENorm.png", dpi = 320, width = 7.5, height = 7.5)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   6. Trimestral                                                              ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Fixed effects 3 Months

Month3FE <- plm(data= merge_data.df,
                formula =  
                  paste0("sum_desplazamiento_pop ~ spraying_norm +",
                         paste(controles_fe_3month, collapse = "+")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date", "query"))
summary(Month3FE)
Month3MCOstd <- coeftest(Month3FE, vcov=vcovHC(Month3FE, type="sss", cluster="group"))  

# IV 3 Months without controls

IV3Month <- plm(data= merge_data.df,
                formula = 
                  paste0("sum_desplazamiento_pop ~ spraying_norm | windSpeedRMBOS"),
                effect = "twoways",
                model = "within",
                index=c("codmpio", "date", "query"))
summary(IV3Month)
IV3Mstd <- coeftest(IV3Month, vcov=vcovHC(IV3Month, type="sss", cluster="group")) 

# IV 3 months with controls

IV3MonthFE <- plm(data= merge_data.df,
                  formula = 
                    paste0("sum_desplazamiento_pop ~ spraying_norm +",
                           paste(controles_fe_3month, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_3month, collapse = "+")),
                  effect = "twoways",
                  model = "within",
                  index=c("codmpio", "date", "query"))
summary(IV3MonthFE)
IVFE3Mstd <- coeftest(IV3MonthFE, vcov=vcovHC(IV3MonthFE, type="sss", cluster="group"))  

stargazer(
  Month3FE, IV3Month, IV3MonthFE,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(Month3MCOstd[, 2], IV3Mstd[, 2], IVFE3Mstd[, 2]), 
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
  p = list(Month3MCOstd[, 4], IV3Mstd[,4], IVFE3Mstd[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

stargazer(
  IV3MonthFE,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = c("spraying_norm", "sum_combates_pop", 'sum_despojo_pop', 'sum_minas_pop', 'sum_reclutamiento_pop', 
           'sum_homicidio_pop', 'sum_desaparicion_pop'),
  se = list(IVFE3Mstd[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas", "Tasa de combates", "Tasa de despojo", "Tasa de minas", "Tasa de reclutamiento",
                       "Tasa de homicidio", "Tasa de desaparición forzada"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si", "Si", "Si"),
                   c("Controles", "Si", "No", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = list(IVFE3Mstd[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   7. Robustez                                                            ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# IV 3 Months without controls

IV3Month <- plm(data= merge_data.df,
                formula = 
                  paste0("sum_desplazamiento_pop ~ spraying_norm | windSpeedRMBOS"),
                effect = "twoways",
                model = "within",
                index=c("codmpio", "date", "query"))
summary(IV3Month)
IV3Mstd <- coeftest(IV3Month, vcov=vcovHC(IV3Month, type="sss", cluster="group")) 

controles_geograficos <- c("rainFall","vegetation", 'windIV10RMBOS')
controles_violencia <- c(
  'sum_combates_pop', 'sum_despojo_pop', 
  'sum_minas_pop', 'sum_reclutamiento_pop', 'sum_homicidio_pop',
  'sum_desaparicion_pop')
controles_economicos <- c('night_lights')

# IV 3 months geograficos

IV3MonthG <- plm(data= merge_data.df,
                 formula = 
                   paste0("sum_desplazamiento_pop ~ spraying_norm +",
                          paste(controles_geograficos, collapse = "+"),"|",
                          paste("windSpeedRMBOS +"),
                          paste(controles_geograficos, collapse = "+")),
                 effect = "twoways",
                 model = "within",
                 index=c("codmpio", "date"))
summary(IV3MonthG)
IV3MGstd <- coeftest(IV3MonthG, vcov=vcovHC(IV3MonthG, type="sss", cluster="group"))  

# IV 3 months violencia

IV3MonthV <- plm(data= merge_data.df,
                 formula = 
                   paste0("sum_desplazamiento_pop ~ spraying_norm +",
                          paste(controles_violencia, collapse = "+"),"|",
                          paste("windSpeedRMBOS +"),
                          paste(controles_violencia, collapse = "+")),
                 effect = "twoways",
                 model = "within",
                 index=c("codmpio", "date"))
summary(IV3MonthV)
IV3MVstd <- coeftest(IV3MonthV, vcov=vcovHC(IV3MonthV, type="sss", cluster="group"))  

# IV 3 months economicos

IV3MonthE <- plm(data= merge_data.df,
                 formula = 
                   paste0("sum_desplazamiento_pop ~ spraying_norm +",
                          paste(controles_economicos, collapse = "+"),"|",
                          paste("windSpeedRMBOS +"),
                          paste(controles_economicos, collapse = "+")),
                 effect = "twoways",
                 model = "within",
                 index=c("codmpio", "date"))
summary(IV3MonthE)
IV3MEstd <- coeftest(IV3MonthE, vcov=vcovHC(IV3MonthE, type="sss", cluster="group"))  

# IV 3 months with controls

IV3MonthFE <- plm(data= merge_data.df,
                  formula = 
                    paste0("sum_desplazamiento_pop ~ spraying_norm +",
                           paste(controles_fe_3month, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_3month, collapse = "+")),
                  effect = "twoways",
                  model = "within",
                  index=c("codmpio", "date", "query"))
summary(IV3MonthFE)
IVFE3Mstd <- coeftest(IV3MonthFE, vcov=vcovHC(IV3MonthFE, type="sss", cluster="group"))  

stargazer(
  IV3Month, IV3MonthG, IV3MonthV, IV3MonthE, IV3MonthFE,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(IV3Mstd[, 2], IV3MGstd[, 2], IV3MVstd[, 2], IV3MEstd[, 2], IVFE3Mstd[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = TRUE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si", "Si", "Si", "Si", "Si"),
                   c("Controles geográficos", "No", "Si", "No", "No", "Si"),
                   c("Controles violencia", "No", "No", "Si", "No", "Si"),
                   c("Controles economicos", "No", "No", "No", "Si", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = list(IV3Mstd[, 4], IV3MGstd[,4], IV3MVstd[, 4], IV3MEstd[, 4], IVFE3Mstd[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

