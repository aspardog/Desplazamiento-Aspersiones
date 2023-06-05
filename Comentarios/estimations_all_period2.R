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

merge_data.df <- readRDS("Data/Merge/output/merge_data_comments.rds") %>%
  mutate(filtro = if_else(str_detect(pattern = "2013", date), 1, 
                          if_else(str_detect(pattern = "2014", date), 1, 0))) %>%
  filter(filtro == 0)

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

