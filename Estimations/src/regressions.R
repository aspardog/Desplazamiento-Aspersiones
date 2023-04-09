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

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Download                                                                                            ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") 
antimerge_data.df <- readRDS("Data/Merge/output/antimerge_data.rds") 

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Correlations                                                                                        ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

controles_fe <- c('vegetation', 'night_lights', "rainFall", 'ruv_desaparicion_forzada', 'ruv_homicidio',
                  'ruv_combates', 'ruv_abandono_despojo', 'cnmh_minas', 'cnmh_reclutamiento', 'windIV10RMBOS')

controles_fe_3month <- c('vegetation', 'night_lights', "rainFall", 'sum_homicidio',
                         'sum_combates', 'sum_despojo', 'sum_minas', 'sum_reclutamiento')

ef <- c('year', 'codmpio', 'month')

#'night_lights'
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Fixed Effects Regression with controls                                                                           ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Regresion desplazamiento ~ aspersion con controles y efectos fijos

#spraying_rate

# finalMonth.reg <- plm(formula =  paste0("desplazamiento_log ~ spraying +", 
#                                         paste(controles_fe, collapse = "+")),
#                       data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
# summary(finalMonth.reg)
# coeftest(finalMonth.reg, vcov=vcovHC(finalMonth.reg, type="sss", cluster="group"))  
# 
# merge_data.df$fittedMonth = predict(finalMonth.reg)
# 
# scatter_month <- merge_data.df %>%
#   filter(spraying > 0) %>%
#   mutate(aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
#   ggplot(data = ., aes(x = aspersiones_log, y = fittedMonth)) +
#   geom_point(mapping=aes(x = aspersiones_log, y = desplazamiento_log)) +
#   geom_smooth(method=lm, aes(y = fittedMonth), color="#C42126", se= T, size = 1) +
#   stat_cor(method = "pearson") +
#   labs(subtitle = "Correlación mensual entre el desplazamiento forzado \ny aspersión aérea",
#        x ="Logaritmo de aspersiones con glifosato", 
#        y = "Logaritmo de desplazamientos forzados") +
#   theme(panel.background   = element_blank(),
#         panel.grid.major   = element_blank(),
#         axis.ticks  = element_blank(),
#         plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
#         axis.text = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
#         axis.title = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
#         plot.background = element_rect(fill = "white", colour = "white"))
# 
# #haven::write_dta(finalRegData.df, path = "Data/Merge/output/merge_data.dta")
# #saveRDS(finalRegData.df, file = "Data/Merge/output/merge_data.rds")
# 
# # 3 month
# 
# final3Month.reg <- plm(formula =  paste0("sum_desplazamiento_log ~ spraying + cultivos +",
#                          paste(controles_fe_3month, collapse = "+")),
#                        data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
# summary(final3Month.reg)
# coeftest(final3Month.reg, vcov=vcovHC(final3Month.reg, type="sss", cluster="group")) 
# 
# merge_data.df$fitted3Month = predict(final3Month.reg)
# 
# scatter_3month <- merge_data.df %>%
#   filter(spraying > 0) %>%
#   mutate(aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
#   ggplot(data = ., aes(x = aspersiones_log, y = fitted3Month)) +
#   geom_point(mapping=aes(x = aspersiones_log, y = sum_desplazamiento_log)) +
#   stat_smooth(method=lm, aes(y = fitted3Month, x = aspersiones_log), color="#C42126", se= T, size = 1) +
#   stat_cor(method = "pearson") +
#   labs(subtitle = "Correlación entre el desplazamiento trimestral \nacumulado y aspersión aérea con glifosato*",
#        x ="Logaritmo de aspersiones con glifosato", 
#        y = "Logaritmo de desplazamientos forzados") +
#   theme(panel.background   = element_blank(),
#         panel.grid.major   = element_blank(),
#         axis.ticks  = element_blank(),
#         plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
#         axis.text = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
#         axis.title = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
#         plot.background = element_rect(fill = "white", colour = "white"));scatter_3month
# 
# figures <- list()
# figures[["Panel A"]] <- scatter_month
# figures[["Panel B"]] <- scatter_3month
# 
# figureScatter <- figures[["Panel A"]] + plot_spacer() + figures[["Panel B"]]  + 
#   plot_layout(ncol = 3, nrow = 1,
#               widths = unit(c(10,2,10), "cm"),
#               heights = unit(15, "cm"))+ 
#   plot_annotation(caption = "*El primer mes de desplazamiento forzado acumulado se cuenta a partir del mes en el que asperjó",
#                   theme = theme(plot.caption = element_text(hjust = 0, size = 8, margin = margin(20, 0, 0, 0))));figureScatter
# 
# ggsave(figureScatter, filename = "Visualizations/output/Scatter.png", dpi = 320, width = 10, height = 7.5)


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
##   4. First Stage                                                               ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

firstStage <- plm(data = merge_data.df,
                  formula = paste0("spraying_norm ~ windIVDaysRMBOS +  windSpeedFLDAS +", paste(controles_fe, collapse = "+")),
                  effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(firstStage)
coeftest(firstStage, vcov=vcovHC(firstStage, type="sss", cluster="group"))  

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Restricción de exclusión                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

restExcl1 <- plm(data = antimerge_data.df,
                formula = paste0("desplazamiento_log ~  windSpeedRMBOS + windIVDaysRMBOS +", paste(controles_fe, collapse = "+")),
                effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(restExcl1)
coeftest(restExcl1, vcov=vcovHC(restExcl1, type="sss", cluster="group"))  


restExcl2 <- plm(data = merge_data.df,
                formula = paste0("desplazamiento_log ~  windSpeedRMBOS + windIVDaysRMBOS+", paste(controles_fe, collapse = "+")),
                effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(restExcl2)
coeftest(restExcl2, vcov=vcovHC(restExcl2, type="sss", cluster="group"))  

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Instrumental Variables Regression                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# Fixed effects

Month1FE <- plm(data= merge_data.df,
                  formula =  
                    paste0("desplazamiento_log ~ spraying +",
                           paste(controles_fe, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(Month1FE)
Month1MCO <- coeftest(Month1FE, vcov=vcovHC(Month1FE, type="sss", cluster="group"))  

# One month

IV1Month <- plm(data= merge_data.df,
                formula =  
                  paste0("desplazamiento_log ~ spraying |",
                         paste("windSpeedRMBOS")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date"))
summary(IV1Month)
IV <- coeftest(IV1Month, vcov=vcovHC(IV1Month, type="sss", cluster="group"))  


# One month controles

IV1MonthFE <- plm(data= merge_data.df,
                formula =  
                  paste0("desplazamiento_log ~ spraying +",
                  paste(controles_fe, collapse = "+"),"|",
                  paste("windSpeedRMBOS  +"),
                  paste(controles_fe, collapse = "+")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date"), diagnostics = TRUE, cluster = "codmpio")
summary(IV1MonthFE)
IVFE <- coeftest(IV1MonthFE, vcov=vcovHC(IV1MonthFE, type="sss", cluster="group"))  

models <- list(
  "(1)"    = Month1MCO,
  "(2)"    = IV,
  "(3)"   = IVFE
)

Nota <- "Los errores estándar son robustos y están corregidos por clusters a nivel municipal.\n *** p<0.01, ** p<0.05, * p<0.1."

s <- stargazer(models, 
          title="Relación Aspersiones Aéreas y Desplazamiento Forzado",
          align=TRUE, 
          dep.var.labels=c("Desplazamiento Forzado Interno"),
          covariate.labels=c("Aspersiones Aéreas con Glifosato"), 
          no.space = TRUE, keep = "spraying", package = T,
          keep.stat = c("rsq", "n"), notes.align = "l")


# coef_rename_ws <- "Aspersiones Aéreas \ncon Glifosato"
# 
# gm <- tibble::tribble(
#   ~raw,        ~clean,          ~fmt,
#   "nobs",      "N",             0,
#   "r.squared", "R<sup>2</sup>", 2,)
# new_rows <- tibble::tribble(
#   ~term,~"(1)",~"(2)",~"(3)",
#   "Controles","Sí","No","Sí",
#   "Efectos Fijos", "Sí","Sí","Sí")
# 
# nota <- "Los errores estándar son robustos y están corregidos por clusters a nivel municipal. *** p<0.01, ** p<0.05, * p<0.1."
# tab <- modelsummary(models, 
#                     estimate = "{estimate}{stars}",
#                     coef_omit = c(-1), 
#                     stars = c('*' = .1, '**' = .05, '***' = 0.01),
#                     output = "latex", fmt = 6, coef_rename = coef_rename_ws, 
#                     gof_map = gm, add_rows = new_rows, 
#                     title = "Aspersiones aéreas",
#                     add_table_notes = nota %>% 
#                       as_chunk() %>% 
#                       add_row() %>% 
#                       row_spec(c(1), font_size = 8))
#   
# table_3 <- tab %>%
#   tab_footnote(footnote = md("Los errores estándar son robustos y están corregidos por clusters a nivel municipal. *** p<0.01, ** p<0.05, * p<0.1."))



# Fixed effects 3 Months

Month3FE <- plm(data= merge_data.df,
                formula =  
                  paste0("sum_desplazamiento_log ~ spraying +",
                         paste(controles_fe_3month, collapse = "+")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date"))
summary(Month3FE)
Month3MCO <- coeftest(Month3FE, vcov=vcovHC(Month3FE, type="sss", cluster="group"))  

# IV 3 Months without controls

IV3Month <- plm(data= merge_data.df,
                formula = 
                  paste0("sum_desplazamiento_log ~ spraying | windSpeedRMBOS"),
                effect = "twoways",
                model = "within",
                index=c("codmpio", "date"))
summary(IV3Month)
IV3M <- coeftest(IV3Month, vcov=vcovHC(IV3Month, type="sss", cluster="group")) 

# IV 3 months with controls

IV3MonthFE <- plm(data= merge_data.df,
                  formula = 
                    paste0("sum_desplazamiento_log ~ spraying +",
                           paste(controles_fe_3month, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe_3month, collapse = "+")),
                  effect = "twoways",
                  model = "within",
                  index=c("codmpio", "date"))
summary(IV3MonthFE)
IVFE3M <- coeftest(IV3MonthFE, vcov=vcovHC(IV3MonthFE, type="sss", cluster="group"))  

models3Month <- list(
  "(MCO)"    = Month3MCO,
  "(IV WC)"  = IV3M,
  "(IV)"     = IVFE3M
)

coef_rename_ws <- "Aspersiones Aéreas \ncon Glifosato"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nobs",      "N",             0,
  "r.squared", "R<sup>2</sup>", 2,)
new_rows <- tibble::tribble(
  ~term,~"(MCO)",~"(IV WC)",~"(IV)",
  "Controles","Sí","No","Sí",
  "Efectos Fijos", "Sí","Sí","Sí")

modelsummary(models3Month, estimate = "{estimate}{stars}", ,
             coef_omit = c(-1), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 5, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Aspersiones aéreas", 
             notes = "*** p<0.01, ** p<0.05, * p<0.1.Los errores estándar son robustos y están corregidos por clusters a nivel municipal.")


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   6. Efectos dinámicos                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# t+1

IV1MonthFE <- plm(data= merge_data.df,
                  formula =  
                    paste0("desplazamiento_log ~ spraying +",
                           paste(controles_fe, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe, collapse = "+")),
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
  drop_na(lag2_ruv_desplazamiento_forzado) %>%
  mutate(desplazamiento_log2 = log(lag2_ruv_desplazamiento_forzado + quantile(lag2_ruv_desplazamiento_forzado, .25)^2/quantile(lag2_ruv_desplazamiento_forzado, .75)))

IV2MonthFE <- plm(data= dataLag2.df,
                  formula =  
                    paste0("desplazamiento_log2 ~ spraying +",
                           paste(controles_fe, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe, collapse = "+")),
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
  drop_na(lag3_ruv_desplazamiento_forzado) %>%
  mutate(desplazamiento_log3 = log(lag3_ruv_desplazamiento_forzado + quantile(lag3_ruv_desplazamiento_forzado, .25)^2/quantile(lag3_ruv_desplazamiento_forzado, .75)))

IV3MonthFE <- plm(data= dataLag3.df,
                  formula =  
                    paste0("desplazamiento_log3 ~ spraying +",
                           paste(controles_fe, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe, collapse = "+")),
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
  drop_na(lag4_ruv_desplazamiento_forzado) %>%
  mutate(desplazamiento_log4 = log(lag4_ruv_desplazamiento_forzado + quantile(lag4_ruv_desplazamiento_forzado, .25)^2/quantile(lag4_ruv_desplazamiento_forzado, .75)))

IV4MonthFE <- plm(data= dataLag4.df,
                  formula =  
                    paste0("desplazamiento_log4 ~ spraying +",
                           paste(controles_fe, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe, collapse = "+")),
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

dataT.df <- merge_data.df %>%
  mutate(desplazamiento_logT = log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)))

IVTMonthFE <- plm(data= dataT.df,
                  formula =  
                    paste0("desplazamiento_logT ~ spraying +",
                           paste(controles_fe, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe, collapse = "+")),
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
  mutate(lagM1_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado, n=1)) %>%
  ungroup() %>%
  drop_na(lagM1_ruv_desplazamiento_forzado) %>%
  mutate(desplazamiento_logM1 = log(lagM1_ruv_desplazamiento_forzado + quantile(lagM1_ruv_desplazamiento_forzado, .25)^2/quantile(lagM1_ruv_desplazamiento_forzado, .75)))

IVM1MonthFE <- plm(data= dataLagM1.df,
                  formula =  
                    paste0("desplazamiento_logM1 ~ spraying +",
                           paste(controles_fe, collapse = "+"),"|",
                           paste("windSpeedRMBOS +"),
                           paste(controles_fe, collapse = "+")),
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
  mutate(lagM2_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado, n=2)) %>%
  ungroup() %>%
  drop_na(lagM2_ruv_desplazamiento_forzado) %>%
  mutate(desplazamiento_logM2 = log(lagM2_ruv_desplazamiento_forzado + quantile(lagM2_ruv_desplazamiento_forzado, .25)^2/quantile(lagM2_ruv_desplazamiento_forzado, .75)))

IVM2MonthFE <- plm(data= dataLagM2.df,
                   formula =  
                     paste0("desplazamiento_logM2 ~ spraying +",
                            paste(controles_fe, collapse = "+"),"|",
                            paste("windSpeedRMBOS +"),
                            paste(controles_fe, collapse = "+")),
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
  mutate(estimate = round(Estimate*100,3),
         lower = round(lower*100,3),
         upper = round(upper*100,2),
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
  scale_y_continuous(limits = c(-0.1, 0.1),
                     breaks = seq(-0.1, 0.1, by = 0.05),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-0.1", "-0.05", "0", "0.05","0.1")) +
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

ggsave(dinamicPlot, filename = "Visualizations/output/DinamicEffects.png", dpi = 320, width = 10, height = 7.5)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   7. Efectos heterogeneos                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVEF <- merge_data.df %>%
  group_by(year) %>%
  mutate(meanYearAspersion = mean(spraying, na.rm = T),
         meanYearCultivos  = mean(cultivos, na.rm = T)) %>%
  ungroup() %>%
  group_by(codmpio, year) %>%
  mutate(meanAspersion = mean(spraying, na.rm = T),
         meanCultivos  = mean(cultivos, na.rm = T)) %>%
  ungroup() 

IVHF <- IVEF %>%
  mutate(hetEffects = if_else(meanAspersion > meanYearAspersion & meanCultivos > meanYearCultivos, "High-High",
                              if_else(meanAspersion > meanYearAspersion & meanCultivos < meanYearCultivos, "High-Low",
                                      if_else(meanAspersion < meanYearAspersion & meanCultivos > meanYearCultivos, "Low-High",
                                              if_else(meanAspersion < meanYearAspersion & meanCultivos < meanYearCultivos, "1Low-Low", NA_character_))))) %>%
  mutate(quintilesAspersion = as.factor(ntile(spraying, 5)),
         hetEffects = as.factor(hetEffects))

IVHReg <- plm(data= IVHF,
               formula = paste0("desplazamiento_log ~ spraying + hetEffects + spraying*hetEffects +",
                                paste(controles_fe, collapse = "+"),"|",
                                paste("windSpeedRMBOS + hetEffects + windSpeedRMBOS*hetEffects +"),
                                paste(controles_fe, collapse = "+")),
               effect = "twoways", 
               model = "within", 
               index=c("codmpio", "date"))
summary(IVHReg)
IVHRegHE <- coeftest(IVHReg, vcov=vcovHC(IVHReg, type="sss", cluster="group")) 
IVHRegHEIC <- confint(IVHReg, vcov=vcovHC(IVHReg, type="sss", cluster="group"), level = 0.9)

estimate <- as.data.frame(IVHRegHE[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying:hetEffectsHigh-High", "spraying:hetEffectsHigh-Low", "spraying:hetEffectsLow-High"))
CI <- as.data.frame(IVHRegHEIC[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying:hetEffectsHigh-High", "spraying:hetEffectsHigh-Low", "spraying:hetEffectsLow-High")) %>%
  rename("lower" = "5 %", "upper" = "95 %")

data2plot <- estimate %>%
  inner_join(CI, by = "variable") %>%
  mutate(
    order_value =
      case_when(
        variable == "spraying:hetEffectsLow-High" ~ 1,
        variable == "spraying:hetEffectsHigh-Low" ~ 2,
        variable == "spraying:hetEffectsHigh-High" ~ 3,
      ),
    variable = 
      case_when(
        variable == "spraying:hetEffectsLow-High" ~ "Baja Aspersión - Altos Cultivos",
        variable == "spraying:hetEffectsHigh-Low" ~ "Alta Aspersión - Bajos Cultivos",
        variable == "spraying:hetEffectsHigh-High" ~ "Alta Aspersión - Altos Cultivos",
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
  scale_y_continuous(limits = c(-0.1, 0.1),
                     breaks = seq(-0.1, 0.1, by = 0.05),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-0.1", "-0.05", "0", "0.05","0.1")) +
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
ggsave(HECocaPlot, filename = "Visualizations/output/CultivosHE.png", dpi = 320, width = 10, height = 7.5)


# Quintiles intensidad


IVASP <- merge_data.df %>%
  group_by(codmpio) %>%
  summarise(asp_total = sum(spraying_norm, na.rm = T)) %>%
  ungroup()

quintiles <- quantile(IVASP$asp_total, probs = seq(0, 1, 0.2))
IVASP$quintil <- cut(IVASP$asp_total, quintiles, labels = FALSE)
IVASP$quintil <- as.factor(IVASP$quintil)

IVASP <- merge_data.df %>%
  left_join(IVASP, by = "codmpio")

IVASP.reg <- plm(data= IVASP,
             formula = paste0("desplazamiento_log ~ quintil + spraying + quintil*spraying +",
                              paste(controles_fe, collapse = "+"),"|",
                              paste("windSpeedRMBOS*quintil + quintil +"),
                              paste(controles_fe, collapse = "+")),
             effect = "individual", 
             model = "random", 
             index=c("codmpio", "date"))
summary(IVASP.reg)
IVASPHE <- coeftest(IVASP.reg, vcov=vcovHC(IVASP.reg, type="sss", cluster="group")) 
IVASPHEIC <- confint(IVASP.reg, vcov=vcovHC(IVASP.reg, type="sss", cluster="group"), level = 0.9)

estimate <- as.data.frame(IVASPHE[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("quintil2:spraying", "quintil3:spraying", "quintil4:spraying", "quintil5:spraying"))
CI <- as.data.frame(IVASPHEIC[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("quintil2:spraying", "quintil3:spraying", "quintil4:spraying", "quintil5:spraying")) %>%
  rename("lower" = "5 %", "upper" = "95 %")

data2plot <- estimate %>%
  inner_join(CI, by = "variable") %>%
  mutate(
    order_value =
      case_when(
        variable == "quintil2:spraying" ~ 1,
        variable == "quintil3:spraying" ~ 2,
        variable == "quintil4:spraying" ~ 3,
        variable == "quintil5:spraying" ~ 4
      ),
    variable = 
      case_when(
        variable == "quintil2:spraying" ~ "Segundo quintil",
        variable == "quintil3:spraying" ~ "Tercer quintil",
        variable == "quintil4:spraying" ~ "Cuarto quintil",
        variable == "quintil5:spraying" ~ "Quinto quintil"
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
  scale_y_continuous(limits = c(-0.1, 0.1),
                     breaks = seq(-0.1, 0.1, by = 0.05),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("-0.1", "-0.05", "0", "0.05","0.1")) +
  coord_flip() +
  labs(x = "Quintil",
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
ggsave(HEPlot, filename = "Visualizations/output/AspersionHE.png", dpi = 320, width = 10, height = 7.5)

