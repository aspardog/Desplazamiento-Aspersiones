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
## This version:      March 9th, 2023
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
p_load(tidyverse, sandwich, lmtest, ivreg, corrr, modelsummary, kableExtra, gt, tibble, stargazer, plm)

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

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") %>%
  dplyr::filter(year < 2013 & year > 2003) %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         dpto = as.factor(dpto),
         month = as.factor(month),
         counter = 1) %>%
  group_by(counter) %>%
  mutate(meanDF = mean(ruv_desplazamiento_forzado, na.rm = T), 
         stdDF = sd(ruv_desplazamiento_forzado, na.rm = T)) %>%
  ungroup() %>%
  mutate(outliers = if_else(ruv_desplazamiento_forzado > meanDF+2*stdDF, 1, 0)) %>%
  filter(outliers == 0)

antimerge_data.df <- readRDS("Data/Merge/output/antimerge_data.rds") %>%
  filter(year < 2013 & year > 2003) %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         month = as.factor(month),
         counter = 1) %>%
  group_by(counter) %>%
  mutate(meanDF = mean(ruv_desplazamiento_forzado, na.rm = T), 
         stdDF = sd(ruv_desplazamiento_forzado, na.rm = T)) %>%
  ungroup() %>%
  mutate(outliers = if_else(ruv_desplazamiento_forzado > meanDF+2*stdDF, 1, 0)) %>%
  filter(outliers == 0)

lassoVariables <- readRDS("Data/Merge/output/lasso_variables.rds")[[1]]
lassoVariables  # Variables to select the controls



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Correlations                                                                                        ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

controles <- c('ruv_amenaza', 'ruv_combates', 
               'ruv_abandono_despojo', "rainFall", "ruv_homicidio",
               'cnmh_minas', 'cnmh_reclutamiento',
               'cnmh_ataque_poblacion')

controles_fe <- c('ruv_amenaza', 'ruv_combates', 'vegetation',
                  'ruv_abandono_despojo', "rainFall", "ruv_homicidio",
                  'cnmh_minas', 'cnmh_reclutamiento',
                  'cnmh_ataque_poblacion','codmpio','year', 'month')

ef <- c('year', 'codmpio', 'month')

controles_fe_3month <- c('sum_amenaza', 'sum_combates', 
                         'sum_despojo', "rainFall", "sum_homicidio",
                         'sum_minas', 'sum_reclutamiento',
                         'sum_ataque','codmpio','year', 'month')


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Fixed Effects Regression with controls                                                                           ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Regresion desplazamiento ~ aspersion con controles y efectos fijos

# One month

finalMonth.reg <- lm(data = merge_data.df, 
          formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying +", paste(controles_fe, collapse = "+")))
summary(finalMonth.reg)
a <- cluster_errors.fn(finalMonth.reg)
merge_data.df$fittedMonth = predict(finalMonth.reg)

prueba <- plm(data = merge_data.df,
              formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ spraying | windSpeedFLDAS | + ruv_amenaza"), 
              effect = "twoways", model = "within", random.method = "ht", inst.method = "baltagi")
summary(prueba)
punbalancedness(merge_data.df)

prueba2 <- plm(data = merge_data.df,
               formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying +", paste(controles, collapse = "+")), 
               effect = "twoways", model = "within")


scatter_month <- merge_data.df %>%
  filter(spraying > 0) %>%
  mutate(desplazamiento_forzado_log = log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)),
         aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = fitted)) +
  geom_point(mapping=aes(x = aspersiones_log, y = desplazamiento_forzado_log)) +
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



# 3 month

final3Month.reg <- lm(data = merge_data.df, 
                     formula = paste0("log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)) ~   spraying +", paste(controles_fe_3month, collapse = "+")))
summary(final3Month.reg)
a3 <- cluster_errors.fn(final3Month.reg)

merge_data.df$fitted3Month = predict(final3Month.reg)

scatter_3month <- merge_data.df %>%
  filter(spraying > 0) %>%
  mutate(desplazamiento_forzado_log = log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)),
         aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = fitted)) +
  geom_point(mapping=aes(x = aspersiones_log, y = desplazamiento_forzado_log)) +
  geom_smooth(method=lm, aes(y = fitted3Month), color="#C42126", se= T, size = 1) +
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
        plot.background = element_rect(fill = "white", colour = "white"))

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
##   4. Instrumental Variables Regression                                                                   ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# One month

windIV.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying | windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(windIV.reg)
cluster_errors.fn(windIV.reg)
library(fixest)

controles_fe <- c('ruv_amenaza', 'ruv_combates', 'vegetation',
                  'ruv_abandono_despojo', "rainFall", "ruv_homicidio",
                  'cnmh_minas', 'cnmh_reclutamiento',
                  'cnmh_ataque_poblacion','codmpio','year', 'month')

windIV.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~",  paste(controles_fe, collapse = "+"), paste(" | spraying | windSpeedFLDAS")))
summary(windIV.reg)
cluster_errors.fn(windIV.reg)

feols(data = merge_data.df,
      fml = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~",  paste(controles_fe, collapse = "+"), paste("+ 0 | spraying | windSpeedFLDAS")), panel.id = ~ codmpio + year)

windIV.reg <- ivreg(data = merge_data.df,
                    formula = paste0("ruv_desplazamiento_forzado ~",  paste(controles_fe, collapse = "+"), paste("| spraying | windSpeedFLDAS")))
summary(windIV.reg)

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

# Option 1

merge_data.df <- merge_data.df %>%
  mutate(regiones = case_when(
    dpto == "AMAZONAS" ~"Amazonia",
    dpto == "CAQUETÁ" ~ "Amazonia",
    dpto == "GUAVIARE" ~ "Amazonia",
    dpto == "PUTUMAYO" ~ "Amazonia",
    dpto == "VAUPÉS" ~ "Amazonia",
    dpto == "ARAUCA" ~ "Orinoquía", 
    dpto == "CASANARE" ~ "Orinoquía", 
    dpto == "META" ~ "Orinoquía", 
    dpto == "NORTE DE SANTANDER SANTANDER" ~ "Orinoquía", 
    dpto == "VICHADA" ~ "Amazonia",
    dpto == "CHOCÓ" ~ "Pacífico", 
    dpto == "NARIÑO" ~ "Pacífico",
    dpto == "VALLE DEL CAUCA" ~ "Pacífico", 
    dpto == "CAUCA" ~ "Pacífico",
    dpto == "ATLÁNTICO" ~ "Caribe",
    dpto == "BOLÍVAR" ~ "Caribe",
    dpto == "CESAR" ~ "Caribe", 
    dpto == "CÓRDOBA" ~ "Caribe", 
    dpto == "LA GUAJIRA" ~ "Caribe", 
    dpto == "MAGDALENA" ~ "Caribe",
    dpto == "SUCRE" ~ "Caribe",
    dpto == "ANTIOQUIA" ~ "Pacífico", 
    dpto == "CALDAS" ~ "Andina", 
    dpto == "CUNDINAMARCA" ~ "Andina", 
    dpto == "HUILA" ~ "Andina",
    dpto == "QUINDIO" ~ "Andina", 
    dpto == "RISARALDA" ~ "Andina", 
    dpto == "SANTANDER" ~ "Andina",
    dpto == "TOLIMA" ~ "Andina",
  )) %>%
  mutate(regiones = if_else(regiones %in% "Pacífico", "AAPacífico", regiones))

# Quintiles
windIV.regional.reg <- ivreg(data = merge_data.df,
                             formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ spraying*intensidad |  windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(windIV.regional.reg)
cluster_errors.fn(windIV.regional.reg) 

# Regiones
windIV.regional.reg <- ivreg(data = merge_data.df,
                    formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~  spraying*as.factor(regiones) |  windSpeedFLDAS + ", paste(controles_fe, collapse = "+")))
summary(windIV.regional.reg)

IVEH <- cluster_errors.fn(windIV.regional.reg) 
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
