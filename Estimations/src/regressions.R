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
       tibble, stargazer, plm,  ggpubr, showtext, patchwork)

EH_panel <- function(mainData = data2plot,
                     line_color = "#003b8a",
                     line_size  = 2,
                     point_color = "#003b8a",
                     point_size   = 4) {
  
  plot <- ggplot(mainData, aes(x = term, y = estimate)) +
    geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
    geom_linerange(aes(x = term,  ymin = lower, ymax = upper),
                   lwd = line_size, position = position_dodge(width = .7), 
                   stat = "identity", color = line_color)+
    geom_point(aes(x = term, y = estimate), 
               size = point_size, position = position_dodge(width = .7), color = point_color) +
    geom_point(aes(x = term, y = estimate), 
               size = 2, position = position_dodge(width = .7), color = "white") +
    scale_y_continuous(limits = c(-0.02, 0.02),
                       breaks = seq(-0.02, 0.02, by = 0.01),
                       expand = expansion(mult = 0.025), position = "right",
                       labels = c("-0.02", "- 0.01", "0", "0.01","0.02"))+
    theme(panel.background   = element_blank(),
          plot.background    = element_blank(),
          panel.grid.major   = element_line(size     = 0.25,
                                            colour   = "#5e5c5a",
                                            linetype = "dashed"),
          panel.grid.minor   = element_blank(),
          axis.ticks  = element_blank(),
          plot.margin  = unit(c(0, 0, 0, 0), "points")) +
    coord_flip() +
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

controles_fe <- c('vegetation', 'cultivos', 'night_lights', "rainFall",
                  'lag1_ruv_combates', 'lag1_ruv_abandono_despojo', 'lag1_cnmh_minas', 'lag1_cnmh_reclutamiento', 'lag1_ruv_homicidio')

controles_fe_3month <- c('vegetation', 'cultivos', 'night_lights', "rainFall",
                         'sum_combates', 'sum_despojo', 'sum_minas', 'sum_reclutamiento', 'sum_homicidio')

ef <- c('year', 'codmpio', 'month')

#'night_lights'
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Fixed Effects Regression with controls                                                                           ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Regresion desplazamiento ~ aspersion con controles y efectos fijos

# One month

finalMonth.reg <- plm(formula =  desplazamiento_log ~
                        spraying + vegetation + cultivos + rainFall + night_lights + lag1_ruv_combates + lag1_ruv_abandono_despojo + lag1_cnmh_minas + lag1_cnmh_reclutamiento + lag1_ruv_homicidio, 
                      data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(finalMonth.reg)
coeftest(finalMonth.reg, vcov=vcovHC(finalMonth.reg, type="sss", cluster="group"))  

finalRegData <- as.data.frame(finalMonth.reg[["model"]]) %>%
  mutate(desplazamiento_log = as.double(desplazamiento_log),
         vegetation = as.double(vegetation),
         night_lights = as.double(night_lights))

codmpioData <- merge_data.df %>%
  select(codmpio, date, year, month, desplazamiento_log, vegetation, night_lights, windSpeedRMBOS, windSpeedFLDAS,
         windIV05RMBOS, windIV10RMBOS, windIV15RMBOS, windIV20RMBOS)

finalRegData.df <- finalRegData %>%
  left_join(codmpioData, by = c("desplazamiento_log", "vegetation", "night_lights"))

finalRegData$fittedMonth = predict(finalMonth.reg)

scatter_month <- finalRegData %>%
  filter(spraying > 0) %>%
  mutate(aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = fittedMonth)) +
  geom_point(mapping=aes(x = aspersiones_log, y = desplazamiento_log)) +
  geom_smooth(method=lm, aes(y = fittedMonth), color="#C42126", se= T, size = 1) +
  stat_cor(method = "pearson") +
  labs(subtitle = "Correlación mensual entre el desplazamiento forzado \ny aspersión aérea",
       x ="Logaritmo de aspersiones con glifosato", 
       y = "Logaritmo de desplazamientos forzados") +
  theme(panel.background   = element_blank(),
        panel.grid.major   = element_blank(),
        axis.ticks  = element_blank(),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
        axis.title = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
        plot.background = element_rect(fill = "white", colour = "white"))

haven::write_dta(finalRegData.df, path = "Data/Merge/output/merge_data.dta")

# 3 month

final3Month.reg <- plm(formula =  sum_desplazamiento_log ~ 
                         spraying + vegetation + cultivos + night_lights + sum_combates + sum_despojo + sum_minas + sum_reclutamiento + sum_homicidio,
                       data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(final3Month.reg)
coeftest(final3Month.reg, vcov=vcovHC(final3Month.reg, type="sss", cluster="group"))  

finalRegData <- as.data.frame(final3Month.reg[["model"]])

finalRegData$fitted3Month = predict(final3Month.reg)

scatter_3month <- finalRegData %>%
  filter(spraying > 0) %>%
  mutate(aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = fitted3Month)) +
  geom_point(mapping=aes(x = aspersiones_log, y = sum_desplazamiento_log)) +
  stat_smooth(method=lm, aes(y = fitted3Month, x = aspersiones_log), color="#C42126", se= T, size = 1) +
  stat_cor(method = "pearson") +
  labs(subtitle = "Correlación entre el desplazamiento trimestral \nacumulado y aspersión aérea con glifosato*",
       x ="Logaritmo de aspersiones con glifosato", 
       y = "Logaritmo de desplazamientos forzados") +
  theme(panel.background   = element_blank(),
        panel.grid.major   = element_blank(),
        axis.ticks  = element_blank(),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
        axis.title = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
        plot.background = element_rect(fill = "white", colour = "white"));scatter_3month

figures <- list()
figures[["Panel A"]] <- scatter_month
figures[["Panel B"]] <- scatter_3month

figureScatter <- figures[["Panel A"]] + plot_spacer() + figures[["Panel B"]]  + 
  plot_layout(ncol = 3, nrow = 1,
              widths = unit(c(10,2,10), "cm"),
              heights = unit(15, "cm"))+ 
  plot_annotation(caption = "*El primer mes de desplazamiento forzado acumulado se cuenta a partir del mes en el que asperjó",
                  theme = theme(plot.caption = element_text(hjust = 0, size = 8, margin = margin(20, 0, 0, 0))));figureScatter

ggsave(figureScatter, filename = "Visualizations/output/Scatter.png", dpi = 320, width = 10, height = 7.5)


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

firstStage <- plm(data = finalRegData.df,
                  formula = paste0("spraying ~  windSpeedRMBOS + ", paste(controles_fe, collapse = "+")),
                  effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(firstStage)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Restricción de exclusión                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

restExcl1 <- plm(data = antimerge_data.df,
                formula = paste0("desplazamiento_log ~  windSpeedRMBOS + ", paste(controles_fe, collapse = "+")),
                effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(restExcl1)

restExcl2 <- plm(data = finalRegData.df,
                formula = paste0("desplazamiento_log ~  windSpeedRMBOS +", paste(controles_fe, collapse = "+")),
                effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(restExcl2)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Instrumental Variables Regression                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# One month

IV1Month <- plm(data= merge_data.df,
                formula =  
                  paste0("desplazamiento_log ~ spraying +",
                  paste(controles_fe, collapse = "+"),"|",
                  paste("windSpeedFLDAS + "),
                  paste(controles_fe, collapse = "+")),
                effect = "twoways", 
                model = "within", 
                index=c("codmpio", "date"))
summary(IV1Month)
coeftest(IV1Month, vcov=vcovHC(IV1Month, type="sss", cluster="group"))  

IV3Month <- plm(data= merge_data.df,
                  formula =  
                    paste0("desplazamiento_log ~ spraying +",
                           paste(controles_fe_3month, collapse = "+"),"|",
                           paste("windSpeedFLDAS +"),
                           paste(controles_fe_3month, collapse = "+")),
                  effect = "twoways", 
                  model = "within", 
                  index=c("codmpio", "date"))
summary(IV3Month)
coeftest(IV3Month, vcov=vcovHC(IV3Month, type="sss", cluster="group"))  


windIV.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~",  paste(controles_fe, collapse = "+"), paste(" | spraying | windSpeedFLDAS + rainFall")))
summary(windIV.reg)
cluster_errors.fn(windIV.reg)


models <- list(
  "(MCO)"    = lm(data = merge_data.df, 
                             formula = paste0("log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)) ~ 0 + spraying +", paste(controles_fe, collapse = "+"))),
  "(IV)"     = ivreg(data = merge_data.df,
                               formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ 0 + spraying | windSpeedFLDAS + ", paste(controles_fe, collapse = "+"))),
  "(IV)"     =  ivreg(data = merge_data.df,
                                   formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ 0 + spraying | windSpeedFLDAS + ", paste(ef, collapse = "+")))
)

coef_rename_ws <- "Aspersiones Aéreas \ncon Glifosato"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nobs",      "N",             0,
  "r.squared", "R<sup>2</sup>", 2,)
new_rows <- tibble::tribble(
  ~term,~"(MCO)",~"(IV)",~"(IV)",
  "Controles","Sí","Sí","No",
  "Efectos Fijos", "Sí","Sí","Sí")

modelsummary(models, vcov = ~codmpio, estimate = "{estimate}{stars}", ,
             coef_omit = c(-2), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 6, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Aspersiones aéreas", 
             notes = "*** p<0.01, ** p<0.05, * p<0.1.Los errores estándar son robustos y están corregidos por clusters de hogar.")

# 3 month

windIV3.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)) ~",  paste(controles_fe, collapse = "+"), paste("| spraying | windSpeedFLDAS")))
summary(windIV3.reg)
cluster_errors.fn(windIV3.reg)


models <- list(
  "(MCO)"    = lm(data = merge_data.df, 
                  formula = paste0("log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)) ~  0 + spraying +", paste(controles_fe_3month, collapse = "+"))),
  "(IV)"     = ivreg(data = merge_data.df,
                     formula = paste0("log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)) ~ 0 + spraying | windSpeedFLDAS + ", paste(controles_fe_3month, collapse = "+"))),
  "(IV)"     =  ivreg(data = merge_data.df,
                      formula = paste0("log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)) ~ 0 + spraying | windSpeedFLDAS + ", paste(ef, collapse = "+")))
)

coef_rename_ws <- "Aspersiones Aéreas \ncon Glifosato"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nobs",      "N",             0,
  "r.squared", "R<sup>2</sup>", 2,)
new_rows <- tibble::tribble(
  ~term,~"(MCO)",~"(IV)",~"(IV)",
  "Controles","Sí","Sí","No",
  "Efectos Fijos", "Sí","Sí","Sí")

modelsummary(models, vcov = ~codmpio, estimate = "{estimate}{stars}", ,
             coef_omit = c(-2), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 6, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Aspersiones aéreas", 
             notes = "*** p<0.01, ** p<0.05, * p<0.1.Los errores estándar son robustos y están corregidos por clusters de hogar.")


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Efectos heterogeneos                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

IVEF <- merge_data.df %>%
  group_by(codmpio) %>%
  summarise(aspersion = sum(spraying, na.rm = T)) %>%
  mutate(quintilesAspersion = ntile(aspersion, 5))

IVRegData.df <- merge_data.df %>%
  left_join(IVEF, by = "codmpio") %>%
  mutate(quintilesAspersion = as.factor(quintilesAspersion)) 
  

# Quintiles
IV1Month <- plm(data= IVRegData.df,
                formula =  
                  paste0("desplazamiento_log ~ spraying*quintilesAspersion + quintilesAspersion +",
                         paste(controles_fe, collapse = "+"),"|",
                         paste("windSpeedFLDAS*quintilesAspersion + quintilesAspersion +"),
                         paste(controles_fe, collapse = "+")),
                effect = "individual", 
                model = "within", 
                index=c("codmpio", "date"))
summary(IV1Month)
coeftest(IV1Month, vcov=vcovHC(IV1Month, type="sss", cluster="group"))  

IV3Month <- plm(data= IVRegData.df,
                formula =  
                  paste0("sum_desplazamiento_log ~ spraying*quintilesAspersion + quintilesAspersion +",
                         paste(controles_fe_3month, collapse = "+"),"|",
                         paste("windSpeedFLDAS*quintilesAspersion + quintilesAspersion +"),
                         paste(controles_fe_3month, collapse = "+")),
                effect = "individual", 
                model = "within", 
                index=c("codmpio", "date"))
summary(IV3Month)


alpha <- 0.05

data2plot <- IVEH %>%
  #filter(p.value < 0.05) %>%
  mutate(filtro = if_else(str_detect(pattern = "spraying:", term), 1, 0)) %>%
  filter(filtro == 1) %>%
  mutate(term = str_replace(pattern = "spraying:as.factor\\(regiones\\)", replacement = "", term)) %>%
  mutate( lower = estimate - qt(1- alpha/2, (n() - 1))*std.error/sqrt(n()),
          upper = estimate + qt(1- alpha/2, (n() - 1))*std.error/sqrt(n()))

regiones_estimation <- EH_panel()

ggsave(regiones_estimation, filename = "Visualizations/output/EH_model.png", dpi = 320, width = 3.5, height = 3.5)

# Option 2

dpto_estimations <- ivreg(data = merge_data.df,
                          formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying*dpto |  windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(dpto_estimations)
DPTO <- cluster_errors.fn(dpto_estimations) 

data2plot <- DPTO %>%
  filter(p.value < 0.05) %>%
  mutate(filtro = if_else(str_detect(pattern = "spraying:", term), 1, 0)) %>%
  filter(filtro == 1) %>%
  mutate(term = str_replace(pattern = "spraying:dpto", replacement = "", term)) %>%
  mutate( lower = estimate - qt(1- alpha/2, (n() - 1))*std.error/sqrt(n()),
          upper = estimate + qt(1- alpha/2, (n() - 1))*std.error/sqrt(n()))

dpto_estimation <- EH_panel() +
  scale_y_continuous(limits = c(-0.2, 0.2),
                     breaks = seq(-0.2, 0.2, by = 0.1),
                     expand = expansion(mult = 0.025), position = "right",
                     labels = c("-0.2", "- 0.1", "0", "0.1","0.2"))

ggsave(dpto_estimation, filename = "Visualizations/output/EH_DPTO_model.png", dpi = 320, width = 3.5, height = 3.5)


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   6.  Robustness Checks                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# First Stage

modelsFS <- list(
  
  "(1)" = lm (data = merge_data.df,
              formula = paste0("spraying ~  windSpeedFLDAS +", paste(controles_fe, collapse = "+")))
  )

coef_rename_ws <- "Velocidad viento"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt, ~omit,
  "nobs",      "N",             0, F,
  "r.squared", "R<sup>2</sup>", 2, F,
  "fstatistic", "F",  3, F)
new_rows <- tibble::tribble(
  ~term,~"(1)",
  "Controles","Sí",
  "Efectos Fijos", "Sí")

modelsummary(modelsFS, vcov = ~codmpio, estimate = "{estimate}{stars}", ,
             coef_omit = c(-2), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 1, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Aspersiones aéreas", 
             notes = "Los errores estándar son robustos y están corregidos por clusters de hogar.  *** p<0.01, ** p<0.05, * p<0.1.")


# Restricción de exclusión

controles_fe <- c('ruv_amenaza', 'ruv_combates', 'vegetation',
                  'ruv_abandono_despojo', "rainFall", "ruv_homicidio",
                  'cnmh_minas', 'cnmh_reclutamiento', "night_lights",
                  'cnmh_ataque_poblacion','codmpio','year', 'month')

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


### Valores atipicos, por encima de dos desviaciones estandar
