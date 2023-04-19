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
EH_panel <- function(mainData = data2plot,
                     line_color = "#003b8a",
                     line_size  = 0.5,
                     point_color = "#003b8a",
                     point_size   = 2.5) {
  
  plot <- ggplot(mainData, aes(x = reorder(label, order_value), y = estimate)) +
    geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
    geom_errorbar(aes(x = reorder(label, order_value),  ymin = lower, ymax = upper),
                  lwd = line_size, position = position_dodge(width = .7), 
                  stat = "identity", color = line_color)+
    geom_point(aes(x = reorder(label, order_value), y = estimate), 
               size = point_size, position = position_dodge(width = .7), color = point_color) +
    geom_point(aes(x = reorder(label, order_value), y = estimate), 
               size = 2, position = position_dodge(width = .7), color = "white") +
    scale_y_continuous(limits = c(-0.2, 0.2),
                       breaks = seq(-0.2, 0.2, by = 0.1),
                       expand = expansion(mult = 0.025), position = "left",
                       labels = c("-0.2", "-0.1", "0", "0.1","0.2"))+
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
          axis.title       = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x  = element_blank(),
          panel.grid.minor.y = element_blank())
  
  return(plot)
}

controles_fe_pop <- c('night_lights', "rainFall","vegetation", 'ruv_abandono_despojo_pop',
                      'ruv_combates_pop', 'ruv_homicidio_pop','cnmh_minas_pop', 
                      'cnmh_desaparicion_pop', 'windIV10RMBOS')

controles_fe_3month <- c('night_lights', "rainFall","vegetation", 'sum_combates_pop', 'sum_despojo_pop', 
                         'sum_minas_pop', 'sum_reclutamiento_pop', 'windIV10RMBOS')


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
                  formula = paste0("spraying_norm ~  windSpeedFLDAS +", paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(firstStage)
stdFirstStage <- coeftest(firstStage, vcov=vcovHC(firstStage, type="sss", cluster="group"))  

stargazer(
  firstStage, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Aspersiones aéreas",
  keep = "windSpeedFLDAS",
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
                 formula = paste0("lag1_ruv_desplazamiento_forzado_pop ~  windSpeedFLDAS +", paste(controles_fe_pop, collapse = "+")),
                 effect = "twoways", model = "within", index=c("codmpio", "date", "query"))

summary(restExcl1)

restExcl1std <- coeftest(restExcl1, vcov=vcovHC(restExcl1, type="sss", cluster="group"))  


restExcl2 <- plm(data = merge_data.df,
                 formula = paste0("lag1_ruv_desplazamiento_forzado_pop ~  windSpeedFLDAS +", paste(controles_fe_pop, collapse = "+")),
                 effect = "twoways", model = "within", index=c("codmpio", "date", "query"))
summary(restExcl2)
restExcl2std <- coeftest(restExcl2, vcov=vcovHC(restExcl2, type="sss", cluster="group"))  

stargazer(
  restExcl1, restExcl2,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "windSpeedFLDAS",
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
##   5. Instrumental Variables Regression                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Fixed effects

Month1FE <- plm(data= merge_data.df,
                formula =  
                  paste0("lag1_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                         paste(controles_fe_pop, collapse = "+")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date", "query"))
summary(Month1FE)
Month1MCOstd <- coeftest(Month1FE, vcov=vcovHC(Month1FE, type="sss", cluster="group"))  

# One month

IV1Month <- plm(data= merge_data.df,
                formula =  
                  paste0("lag1_ruv_desplazamiento_forzado_pop ~ spraying_norm |",
                         paste(" windSpeedFLDAS")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date", "query"))
summary(IV1Month)
IVstd <- coeftest(IV1Month, vcov=vcovHC(IV1Month, type="sss", cluster="group"))  


# One month controles

IV1MonthFE <- plm(data= merge_data.df,
                  formula =  
                    paste0("lag1_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste(" windSpeedFLDAS  +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date", "query"), diagnostics = TRUE)
summary(IV1MonthFE)
IVFEstd <- coeftest(IV1MonthFE, vcov=vcovHC(IV1MonthFE, type="sss", cluster="group"))  

stargazer(
  Month1FE, IV1Month, IV1MonthFE,
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(Month1MCOstd[, 2], IVstd[, 2], IVFEstd[, 2]), 
  title = "Resultados primera etapa", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = TRUE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si", "Si", "Si"),
                   c("Controles", "Si", "No", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "7pt",
  keep.stat = c("rsq", "n"), 
  p = list(Month1MCOstd[, 4], IVstd[,4], IVFEstd[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
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
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(IV1MonthFE)
IVFE <- coeftest(IV1MonthFE, vcov=vcovHC(IV1MonthFE, type="sss", cluster="group"))  
IVFECI <- confint(IVFE, level = 0.9)

IVFET1 <- data.frame("Estimate" = IVFE[[1]],
                     "std" = IVFE[[1,2]],
                     "p.value" = IVFE[[1,4]],
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
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(IV2MonthFE)
IVFE2 <- coeftest(IV2MonthFE, vcov=vcovHC(IV2MonthFE, type="sss", cluster="group"))  
IVFE2CI <- confint(IVFE2, level = 0.9)

IVFET2 <- data.frame("Estimate" = IVFE2[[1]],
                     "std" = IVFE2[[1,2]],
                     "p.value" = IVFE2[[1,4]],
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
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(IV3MonthFE)
IVFE3 <- coeftest(IV3MonthFE, vcov=vcovHC(IV3MonthFE, type="sss", cluster="group"), parm = ci) 
IVFE3CI <- confint(IVFE3, level = 0.9)

IVFET3 <- data.frame("Estimate" = IVFE3[[1]],
                     "std" = IVFE3[[1,2]],
                     "p.value" = IVFE3[[1,4]],
                     "label" = "t+3",
                     "lower" = IVFE3CI[[1]],
                     "upper" = IVFE3CI[[1,2]])
# t+4

dataLag4.df <- merge_data.df %>%
  drop_na(lag4_ruv_desplazamiento_forzado)

IV4MonthFE <- plm(data= dataLag4.df,
                  formula =  
                    paste0("lag4_ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(IV4MonthFE)
IVFE4 <- coeftest(IV4MonthFE, vcov=vcovHC(IV4MonthFE, type="sss", cluster="group"))  
IVFE4CI <- confint(IVFE4, level = 0.9)

IVFET4 <- data.frame("Estimate" = IVFE4[[1]],
                     "std" = IVFE4[[1,2]],
                     "p.value" = IVFE4[[1,4]],
                     "label" = "t+4",
                     "lower" = IVFE4CI[[1]],
                     "upper" = IVFE4CI[[1,2]])
#t

dataT.df <- merge_data.df 

IVTMonthFE <- plm(data= dataT.df,
                  formula =  
                    paste0("ruv_desplazamiento_forzado_pop ~ spraying_norm +",
                           paste(controles_fe_pop, collapse = "+"),"|",
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_pop, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(IVTMonthFE)
IVFET <- coeftest(IVTMonthFE, vcov=vcovHC(IVTMonthFE, type="sss", cluster="group"))  
IVFETCI <- confint(IVFET, level = 0.9)

IVFET  <- data.frame("Estimate" = IVFET[[1]],
                     "std" = IVFET[[1,2]],
                     "p.value" = IVFET[[1,4]],
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
                            paste("windSpeedFLDAS +"),
                            paste(controles_fe_pop, collapse = "+")),
                   effect = "twoways", 
                   model = "within", 
                   index=c("codmpio", "date"))
summary(IVM1MonthFE)
IVFEM1 <- coeftest(IVM1MonthFE, vcov=vcovHC(IVM1MonthFE, type="sss", cluster="group"))  
IVFEM1CI <- confint(IVFEM1, level = 0.9)

IVFEM1  <- data.frame("Estimate" = IVFEM1[[1]],
                      "std" = IVFEM1[[1,2]],
                      "p.value" = IVFEM1[[1,4]],
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
                            paste("windSpeedFLDAS +"),
                            paste(controles_fe_pop, collapse = "+")),
                   effect = "twoways", 
                   model = "within", 
                   index=c("codmpio", "date"))
summary(IVM2MonthFE)
IVFEM2 <- coeftest(IVM2MonthFE, vcov=vcovHC(IVM2MonthFE, type="sss", cluster="group"))  
IVFEM2CI <- confint(IVFEM2, level = 0.9)

IVFEM2  <- data.frame("Estimate" = IVFEM2[[1]],
                      "std" = IVFEM2[[1,2]],
                      "p.value" = IVFEM2[[1,4]],
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
##   5. Trimestral                                                              ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
controles_fe_3month <- c('night_lights', "rainFall","vegetation", 'sum_combates_pop', 'sum_despojo_pop', 
                         'sum_minas_pop', 'sum_reclutamiento_pop', 'sum_homicidio')
# Fixed effects 3 Months

Month3FE <- plm(data= merge_data.df,
                formula =  
                  paste0("sum_desplazamiento_pop ~ spraying_norm +",
                         paste(controles_fe_3month, collapse = "+")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date"))
summary(Month3FE)
Month3MCOstd <- coeftest(Month3FE, vcov=vcovHC(Month3FE, type="sss", cluster="group"))  

# IV 3 Months without controls

IV3Month <- plm(data= merge_data.df,
                formula = 
                  paste0("sum_desplazamiento_pop ~ spraying_norm | windSpeedFLDAS"),
                effect = "twoways",
                model = "within",
                index=c("codmpio", "date"))
summary(IV3Month)
IV3Mstd <- coeftest(IV3Month, vcov=vcovHC(IV3Month, type="sss", cluster="group")) 

# IV 3 months with controls

IV3MonthFE <- plm(data= merge_data.df,
                  formula = 
                    paste0("sum_desplazamiento_pop ~ spraying_norm +",
                           paste(controles_fe_3month, collapse = "+"),"|",
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_3month, collapse = "+")),
                  effect = "twoways",
                  model = "within",
                  index=c("codmpio", "date"))
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