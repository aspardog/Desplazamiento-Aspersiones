## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     June 3rd, 2023
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
p_load(tidyverse, sandwich, lmtest, ivreg, corrr, modelsummary, kableExtra, gt, sf,
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


mainPath <- ("Data/Download/output/")

# Main DataBases

ruv.df         <- readRDS(paste0(mainPath, "Victimization/ruv_victimization.rds"))
aspersiones.df <- readRDS(paste0(mainPath, "Spraying/aspersiones.rds")) %>%
  dplyr::select(!c(year, month, dpto, mpio)) 

# IV DataBase

FLDAS.df <- readRDS(paste0(mainPath, "Wind/windFLDAS.rds"))
MERRA.df <- readRDS(paste0(mainPath, "Wind/windMERRA2.rds"))
RMOBS.df <- readRDS(paste0(mainPath, "Wind/windRMBOS.rds"))

# Controls DataBase
colombia.sf <- st_read('Data/Download/input/ShapeFiles/Municipio_ctm12.shp') %>%
  st_transform(crs = 4326) %>%
  select(mpio_ccnct, mpio_cnmbr, dpto_ccdgo, dpto_cnmbr, mpio_narea) %>%
  rename(codmpio = mpio_ccnct, coddpto = dpto_ccdgo, mpio = mpio_cnmbr, dpto = dpto_cnmbr, mpio_area = mpio_narea) %>%
  mutate(codmpio = as.numeric(codmpio)) %>%
  st_drop_geometry()

cnmh.df          <- readRDS(paste0(mainPath, "Victimization/cnmh_victimization.rds"))
nigthLights.df   <- readRDS(paste0(mainPath, "Lights/nightLights.rds")) %>%
  dplyr::select(!c(year, month))
vegetation.df    <- readRDS(paste0(mainPath, "Vegetation/vegetation.rds")) %>%
  dplyr::select(!c(year, month))
rainfall.df      <- readRDS(paste0(mainPath, "Rainfall/rainFall.rds"))%>%
  dplyr::select(!c(year, month))
distance.df      <- readRDS(paste0(mainPath, "Airports/airports_distance.rds"))
coca.df          <- readRDS(paste0(mainPath, "Coca/coca.rds")) %>%
  mutate(codmpio = as.numeric(codmpio))
population.df    <- readRDS(paste0(mainPath, "Population/population.rds")) %>%
  mutate(codmpio = as.numeric(codmpio))

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Build the IV                                                                                       ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# This function adapts the date according to the merge database

windDate <- function(data) {
  
  wind.df <- data %>%
    mutate(year = as.numeric(format(date, format = "%Y")),
           month = as.numeric(format(date, format = "%m"))) %>%
    mutate(month   = if_else(month < 10, paste0("0", month), as.character(month)),
           date    = paste0(year,month)) %>%
    filter(year > 2002 & year < 2015)
  
}

# This function creates the wind shocks IV

iv_windShocks <- function(data) {
  
  wind.df <- data %>%
    mutate(year = as.numeric(format(date, format = "%Y")),
           month = as.numeric(format(date, format = "%m"))) %>%
    mutate(month   = if_else(month < 10, paste0("0", month), as.character(month)),
           date    = paste0(year,month)) %>%
    filter(year > 2003 & year < 2013) %>%
    mutate(wind_speed = if_else(wind_speed == "NaN", NA_real_, wind_speed))
  
  wind_iv.df <- wind.df  %>%
    group_by(month, codmpio) %>%
    mutate(meanWind = mean(wind_speed, na.rm = T), 
           stdWind = sd(wind_speed, na.rm = T)) %>%
    mutate(wind_shocks05 = if_else(wind_speed > meanWind+0.5*stdWind, 1, 0),
           wind_shocks10 = if_else(wind_speed > meanWind+stdWind, 1, 0),
           wind_shocks15 = if_else(wind_speed > meanWind+1.5*stdWind, 1, 0),
           wind_shocks20 = if_else(wind_speed > meanWind+2*stdWind, 1, 0),
           wind_shocks25 = if_else(wind_speed > meanWind+2.5*stdWind, 1, 0),
           wind_shocks30 = if_else(wind_speed > meanWind+3*stdWind, 1, 0)) %>%
    mutate(wind_days = if_else(wind_speed > 1.93, 1, 0)) %>%
    ungroup() %>%
    distinct() %>%
    group_by(date, codmpio) %>%
    summarise(windIV05 = sum(wind_shocks05, na.rm = T),
              windIV10 = sum(wind_shocks10, na.rm = T),
              windIV15 = sum(wind_shocks15, na.rm = T),
              windIV20 = sum(wind_shocks20, na.rm = T),
              windIV25 = sum(wind_shocks25, na.rm = T),
              windIV30 = sum(wind_shocks30, na.rm = T),
              windIVDays = sum(wind_days, na.rm = T)) %>%
    arrange(codmpio, date)
  
}

# Arrange the dates

windRMBOS.df <- windDate(RMOBS.df) %>%
  group_by(date, codmpio) %>%
  summarise(windSpeedRMBOS = mean(wind_speed))
windMERRA.df <- windDate(MERRA.df) %>%
  group_by(date, codmpio) %>%
  summarise(windSpeedMERRA = mean(wind_speed))
windFLDAS.df <- windDate(FLDAS.df) %>%
  rename(windSpeedFLDAS = wind_speed) 

# Create IV shocks
windShocksRMOBS <- iv_windShocks(data = RMOBS.df) %>%
  rename(windIV05RMBOS = windIV05,
         windIV10RMBOS = windIV10,
         windIV15RMBOS = windIV15,
         windIV20RMBOS = windIV20,
         windIV25RMBOS = windIV25,
         windIV30RMBOS = windIV30,
         windIVDaysRMBOS = windIVDays
  )
windShocksMERRA <- iv_windShocks(data = MERRA.df) %>%
  rename(windIV05MERRA = windIV05,
         windIV10MERRA = windIV10,
         windIV15MERRA = windIV15,
         windIV20MERRA = windIV20,
         windIV25MERRA = windIV25,
         windIV30MERRA = windIV30,
         windIVDaysMERRA = windIVDays)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   3. Merge                                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Merge

setdiff(ruv.df$date,aspersiones.df$date)

merge_ruv_aspersion.df <- ruv.df %>%
  inner_join(aspersiones.df, by = c("codmpio", "date")) 
#mutate(across(!c(mpio, dpto, date, longitud, latitud),
# ~if_else(is.na(.x), 0, .x)))

merge_data.df <- merge_ruv_aspersion.df %>%
  left_join(y = cnmh.df, by = c("codmpio", "date")) %>%
  mutate(across(starts_with("cnmh"),
                ~if_else(is.na(.x), 0, .x))) %>%
  left_join(y = nigthLights.df, by = c("codmpio", "date")) %>%
  left_join(y = windRMBOS.df, by = c("codmpio", "date")) %>%
  left_join(y = windMERRA.df, by = c("codmpio", "date")) %>%
  left_join(y = windFLDAS.df, by = c("codmpio", "date")) %>%
  left_join(y = windShocksRMOBS, by = c("codmpio", "date")) %>%
  left_join(y = windShocksMERRA, by = c("codmpio", "date")) %>%
  left_join(y = vegetation.df, by = c("codmpio", "date")) %>%
  left_join(y = rainfall.df, by = c("codmpio", "date")) %>%
  left_join(y = coca.df, by = c("codmpio", "year")) %>%
  dplyr::select(!year) %>%
  left_join(y = distance.df, by = c("codmpio")) %>%
  left_join(y = colombia.sf, by = c("codmpio")) %>%
  left_join(y = population.df, by = c("codmpio"))

# Cleaning year

merge_data.df <- merge_data.df %>%
  mutate(year = as.numeric(year(ym(date)))) %>%
  filter(year < 2015 & year > 2002)

merge_data_lags.df <- merge_data.df %>%
  mutate(order_value = (as.numeric(date))) %>% 
  mutate(ruv_desplazamiento_forzado_pop = (ruv_desplazamiento_forzado/population)*100,
         ruv_homicidio_pop          = (ruv_homicidio/population)*100,
         ruv_amenaza_pop            = (ruv_amenaza/population)*100,
         ruv_combates_pop           = (ruv_combates/population)*100,
         ruv_abandono_despojo_pop   = (ruv_abandono_despojo/population)*100,
         cnmh_minas_pop             = (cnmh_minas/population)*100,
         cnmh_reclutamiento_pop     = (cnmh_reclutamiento/population)*100,
         cnmh_ataque_poblacion_pop  = (cnmh_ataque_poblacion/population)*100,
         cnmh_asesinatos_selectivos_pop = (cnmh_asesinatos_selectivos/population)*100,
         cnmh_atentado_pop          = (cnmh_atentado/population)*100,
         cnmh_d_bienes_pop          = (cnmh_d_bienes/population)*100,
         cnmh_desaparicion_pop      = (cnmh_desaparicion/population)*100,
         cnmh_enfrentamientos_pop   = (cnmh_enfrentamientos/population)*100,
         cnmh_masacres_pop          = (cnmh_masacres/population)*100,
         cnmh_secuestro_pop         = (cnmh_secuestro/population)*100,
         cnmh_violencia_sexual_pop  = (cnmh_violencia_sexual/population)*100,
         ruv_tortura_pop            = (ruv_tortura/population)*100,
         ruv_desaparicion_forzada_pop = (ruv_desaparicion_forzada/population)*100,
         ruv_secuestro_pop          = (ruv_secuestro/population)*100,
         ruv_violencia_sexual_pop   = (ruv_violencia_sexual/population)*100,
         ruv_perdida_bienes_pop     = (ruv_perdida_bienes/population)*100,
         mpio_area = mpio_area*100, # Trasnformación a hectareas
         spraying_norm = round((spraying/mpio_area)*100,2),
         vegetation_norm = round((vegetation/mpio_area)*100,2),
  ) %>%
  drop_na(cultivos, windSpeedRMBOS) %>%
  arrange(-order_value) %>%
  group_by(codmpio) %>%
  mutate(lag1_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado, n=1),
         lag2_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado, n=2),
         lag3_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado, n=3),
         lag4_ruv_desplazamiento_forzado = dplyr::lag(ruv_desplazamiento_forzado, n=4),
         lag1_ruv_homicidio = dplyr::lag(ruv_homicidio, n=1),
         lag2_ruv_homicidio = dplyr::lag(ruv_homicidio, n=2),
         lag3_ruv_homicidio = dplyr::lag(ruv_homicidio, n=3),
         lag1_ruv_amenaza = dplyr::lag(ruv_amenaza, n=1),
         lag2_ruv_amenaza = dplyr::lag(ruv_amenaza, n=2),
         lag3_ruv_amenaza = dplyr::lag(ruv_amenaza, n=3),
         lag1_ruv_combates = dplyr::lag(ruv_combates, n=1),
         lag2_ruv_combates = dplyr::lag(ruv_combates, n=2),
         lag3_ruv_combates = dplyr::lag(ruv_combates, n=3),
         lag1_ruv_abandono_despojo = dplyr::lag(ruv_abandono_despojo, n=1),
         lag2_ruv_abandono_despojo = dplyr::lag(ruv_abandono_despojo, n=2),
         lag3_ruv_abandono_despojo = dplyr::lag(ruv_abandono_despojo, n=3),
         lag1_cnmh_minas = dplyr::lag(cnmh_minas, n=1),
         lag2_cnmh_minas = dplyr::lag(cnmh_minas, n=2),
         lag3_cnmh_minas = dplyr::lag(cnmh_minas, n=3),
         lag1_cnmh_reclutamiento = dplyr::lag(cnmh_reclutamiento, n=1),
         lag2_cnmh_reclutamiento = dplyr::lag(cnmh_reclutamiento, n=2),
         lag3_cnmh_reclutamiento = dplyr::lag(cnmh_reclutamiento, n=3),
         lag1_cnmh_ataque_poblacion = dplyr::lag(cnmh_ataque_poblacion, n=1),
         lag2_cnmh_ataque_poblacion = dplyr::lag(cnmh_ataque_poblacion, n=2),
         lag3_cnmh_ataque_poblacion = dplyr::lag(cnmh_ataque_poblacion, n=3),
         lag1_cnmh_asesinatos_selectivos = dplyr::lag(cnmh_asesinatos_selectivos, n=1),
         lag2_cnmh_asesinatos_selectivos = dplyr::lag(cnmh_asesinatos_selectivos, n=2),
         lag3_cnmh_asesinatos_selectivos = dplyr::lag(cnmh_asesinatos_selectivos, n=3),
         lag1_cnmh_atentado = dplyr::lag(cnmh_atentado, n=1),
         lag2_cnmh_atentado = dplyr::lag(cnmh_atentado, n=2),
         lag3_cnmh_atentado = dplyr::lag(cnmh_atentado, n=3),
         lag1_cnmh_d_bienes = dplyr::lag(cnmh_d_bienes, n=1),
         lag2_cnmh_d_bienes = dplyr::lag(cnmh_d_bienes, n=2),
         lag3_cnmh_d_bienes = dplyr::lag(cnmh_d_bienes, n=3),
         lag1_cnmh_desaparicion = dplyr::lag(cnmh_desaparicion, n=1),
         lag2_cnmh_desaparicion = dplyr::lag(cnmh_desaparicion, n=2),
         lag3_cnmh_desaparicion = dplyr::lag(cnmh_desaparicion, n=3),
         lag1_cnmh_enfrentamientos = dplyr::lag(cnmh_enfrentamientos, n=1),
         lag2_cnmh_enfrentamientos = dplyr::lag(cnmh_enfrentamientos, n=2),
         lag3_cnmh_enfrentamientos = dplyr::lag(cnmh_enfrentamientos, n=3),
         lag1_cnmh_masacres = dplyr::lag(cnmh_masacres, n=1),
         lag2_cnmh_masacres = dplyr::lag(cnmh_masacres, n=2),
         lag3_cnmh_masacres = dplyr::lag(cnmh_masacres, n=3),
         lag1_cnmh_secuestro = dplyr::lag(cnmh_secuestro, n=1),
         lag2_cnmh_secuestro = dplyr::lag(cnmh_secuestro, n=2),
         lag3_cnmh_secuestro = dplyr::lag(cnmh_secuestro, n=3),
         lag1_cnmh_violencia_sexual = dplyr::lag(cnmh_violencia_sexual, n=1),
         lag2_cnmh_violencia_sexual = dplyr::lag(cnmh_violencia_sexual, n=2),
         lag3_cnmh_violencia_sexual = dplyr::lag(cnmh_violencia_sexual, n=3),
         lag1_ruv_tortura = dplyr::lag(ruv_tortura, n=1),
         lag2_ruv_tortura = dplyr::lag(ruv_tortura, n=2),
         lag3_ruv_tortura = dplyr::lag(ruv_tortura, n=3),
         lag1_ruv_desaparicion_forzada = dplyr::lag(ruv_desaparicion_forzada, n=1),
         lag2_ruv_desaparicion_forzada = dplyr::lag(ruv_desaparicion_forzada, n=2),
         lag3_ruv_desaparicion_forzada = dplyr::lag(ruv_desaparicion_forzada, n=3),
         lag1_ruv_secuestro = dplyr::lag(ruv_secuestro, n=1),
         lag2_ruv_secuestro = dplyr::lag(ruv_secuestro, n=2),
         lag3_ruv_secuestro = dplyr::lag(ruv_secuestro, n=3),
         lag1_ruv_violencia_sexual = dplyr::lag(ruv_violencia_sexual, n=1),
         lag2_ruv_violencia_sexual = dplyr::lag(ruv_violencia_sexual, n=2),
         lag3_ruv_violencia_sexual = dplyr::lag(ruv_violencia_sexual, n=3),
         lag1_ruv_perdida_bienes = dplyr::lag(ruv_perdida_bienes, n=1),
         lag2_ruv_perdida_bienes = dplyr::lag(ruv_perdida_bienes, n=2),
         lag3_ruv_perdida_bienes = dplyr::lag(ruv_perdida_bienes, n=3),
         lag1_vegetation = dplyr::lag(vegetation, n = 1),
         lag2_vegetation = dplyr::lag(vegetation, n = 2),
         lag3_vegetation = dplyr::lag(vegetation, n = 3),
         lag1_vegetation_norm = dplyr::lag(vegetation_norm, n = 1),
         lag2_vegetation_norm = dplyr::lag(vegetation_norm, n = 2),
         lag3_vegetation_norm = dplyr::lag(vegetation_norm, n = 3),
         lag1_ruv_desplazamiento_forzado_pop = dplyr::lag(ruv_desplazamiento_forzado_pop, n=1),
         lag2_ruv_desplazamiento_forzado_pop = dplyr::lag(ruv_desplazamiento_forzado_pop, n=2),
         lag3_ruv_desplazamiento_forzado_pop = dplyr::lag(ruv_desplazamiento_forzado_pop, n=3),
         lag4_ruv_desplazamiento_forzado_pop = dplyr::lag(ruv_desplazamiento_forzado_pop, n=4),
         lag1_ruv_homicidio_pop = dplyr::lag(ruv_homicidio_pop, n=1),
         lag2_ruv_homicidio_pop = dplyr::lag(ruv_homicidio_pop, n=2),
         lag3_ruv_homicidio_pop = dplyr::lag(ruv_homicidio_pop, n=3),
         lag1_ruv_amenaza_pop = dplyr::lag(ruv_amenaza_pop, n=1),
         lag2_ruv_amenaza_pop = dplyr::lag(ruv_amenaza_pop, n=2),
         lag3_ruv_amenaza_pop = dplyr::lag(ruv_amenaza_pop, n=3),
         lag1_ruv_combates_pop = dplyr::lag(ruv_combates_pop, n=1),
         lag2_ruv_combates_pop = dplyr::lag(ruv_combates_pop, n=2),
         lag3_ruv_combates_pop = dplyr::lag(ruv_combates_pop, n=3),
         lag1_ruv_abandono_despojo_pop = dplyr::lag(ruv_abandono_despojo_pop, n=1),
         lag2_ruv_abandono_despojo_pop = dplyr::lag(ruv_abandono_despojo_pop, n=2),
         lag3_ruv_abandono_despojo_pop = dplyr::lag(ruv_abandono_despojo_pop, n=3),
         lag1_cnmh_minas_pop = dplyr::lag(cnmh_minas_pop, n=1),
         lag2_cnmh_minas_pop = dplyr::lag(cnmh_minas_pop, n=2),
         lag3_cnmh_minas_pop = dplyr::lag(cnmh_minas_pop, n=3),
         lag1_cnmh_reclutamiento_pop = dplyr::lag(cnmh_reclutamiento_pop, n=1),
         lag2_cnmh_reclutamiento_pop = dplyr::lag(cnmh_reclutamiento_pop, n=2),
         lag3_cnmh_reclutamiento_pop = dplyr::lag(cnmh_reclutamiento_pop, n=3),
         lag1_cnmh_ataque_poblacion_pop = dplyr::lag(cnmh_ataque_poblacion_pop, n=1),
         lag2_cnmh_ataque_poblacion_pop = dplyr::lag(cnmh_ataque_poblacion_pop, n=2),
         lag3_cnmh_ataque_poblacion_pop = dplyr::lag(cnmh_ataque_poblacion_pop, n=3),
         lag1_cnmh_asesinatos_selectivos_pop = dplyr::lag(cnmh_asesinatos_selectivos_pop, n=1),
         lag2_cnmh_asesinatos_selectivos_pop = dplyr::lag(cnmh_asesinatos_selectivos_pop, n=2),
         lag3_cnmh_asesinatos_selectivos_pop = dplyr::lag(cnmh_asesinatos_selectivos_pop, n=3),
         lag1_cnmh_atentado_pop = dplyr::lag(cnmh_atentado_pop, n=1),
         lag2_cnmh_atentado_pop = dplyr::lag(cnmh_atentado_pop, n=2),
         lag3_cnmh_atentado_pop = dplyr::lag(cnmh_atentado_pop, n=3),
         lag1_cnmh_d_bienes_pop = dplyr::lag(cnmh_d_bienes_pop, n=1),
         lag2_cnmh_d_bienes_pop = dplyr::lag(cnmh_d_bienes_pop, n=2),
         lag3_cnmh_d_bienes_pop = dplyr::lag(cnmh_d_bienes_pop, n=3),
         lag1_cnmh_desaparicion_pop = dplyr::lag(cnmh_desaparicion_pop, n=1),
         lag2_cnmh_desaparicion_pop = dplyr::lag(cnmh_desaparicion_pop, n=2),
         lag3_cnmh_desaparicion_pop = dplyr::lag(cnmh_desaparicion_pop, n=3),
         lag1_cnmh_enfrentamientos_pop = dplyr::lag(cnmh_enfrentamientos_pop, n=1),
         lag2_cnmh_enfrentamientos_pop = dplyr::lag(cnmh_enfrentamientos_pop, n=2),
         lag3_cnmh_enfrentamientos_pop = dplyr::lag(cnmh_enfrentamientos_pop, n=3),
         lag1_cnmh_masacres_pop = dplyr::lag(cnmh_masacres_pop, n=1),
         lag2_cnmh_masacres_pop = dplyr::lag(cnmh_masacres_pop, n=2),
         lag3_cnmh_masacres_pop = dplyr::lag(cnmh_masacres_pop, n=3),
         lag1_cnmh_secuestro_pop = dplyr::lag(cnmh_secuestro_pop, n=1),
         lag2_cnmh_secuestro_pop = dplyr::lag(cnmh_secuestro_pop, n=2),
         lag3_cnmh_secuestro_pop = dplyr::lag(cnmh_secuestro_pop, n=3),
         lag1_cnmh_violencia_sexual_pop = dplyr::lag(cnmh_violencia_sexual_pop, n=1),
         lag2_cnmh_violencia_sexual_pop = dplyr::lag(cnmh_violencia_sexual_pop, n=2),
         lag3_cnmh_violencia_sexual_pop = dplyr::lag(cnmh_violencia_sexual_pop, n=3),
         lag1_ruv_tortura_pop = dplyr::lag(ruv_tortura_pop, n=1),
         lag2_ruv_tortura_pop = dplyr::lag(ruv_tortura_pop, n=2),
         lag3_ruv_tortura_pop = dplyr::lag(ruv_tortura_pop, n=3),
         lag1_ruv_desaparicion_forzada_pop = dplyr::lag(ruv_desaparicion_forzada_pop, n=1),
         lag2_ruv_desaparicion_forzada_pop = dplyr::lag(ruv_desaparicion_forzada_pop, n=2),
         lag3_ruv_desaparicion_forzada_pop = dplyr::lag(ruv_desaparicion_forzada_pop, n=3),
         lag1_ruv_secuestro_pop = dplyr::lag(ruv_secuestro_pop, n=1),
         lag2_ruv_secuestro_pop = dplyr::lag(ruv_secuestro_pop, n=2),
         lag3_ruv_secuestro_pop = dplyr::lag(ruv_secuestro_pop, n=3),
         lag1_ruv_violencia_sexual_pop = dplyr::lag(ruv_violencia_sexual_pop, n=1),
         lag2_ruv_violencia_sexual_pop = dplyr::lag(ruv_violencia_sexual_pop, n=2),
         lag3_ruv_violencia_sexual_pop = dplyr::lag(ruv_violencia_sexual_pop, n=3),
         lag1_ruv_perdida_bienes_pop = dplyr::lag(ruv_perdida_bienes_pop, n=1),
         lag2_ruv_perdida_bienes_pop = dplyr::lag(ruv_perdida_bienes_pop, n=2),
         lag3_ruv_perdida_bienes_pop = dplyr::lag(ruv_perdida_bienes_pop, n=3)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(sum_desplazamiento = sum(lag1_ruv_desplazamiento_forzado, lag2_ruv_desplazamiento_forzado, lag3_ruv_desplazamiento_forzado, na.rm = T),
         sum_homicidio      = sum(lag1_ruv_homicidio, lag2_ruv_homicidio, lag3_ruv_homicidio, na.rm = T),
         sum_amenaza        = sum(lag1_ruv_amenaza, lag2_ruv_amenaza, lag3_ruv_amenaza, na.rm = T),
         sum_combates       = sum(lag1_ruv_combates, lag2_ruv_combates, lag3_ruv_combates, na.rm = T),
         sum_despojo        = sum(lag1_ruv_abandono_despojo, lag2_ruv_abandono_despojo, lag3_ruv_abandono_despojo, na.rm = T),
         sum_minas          = sum(lag1_cnmh_minas, lag2_cnmh_minas, lag3_cnmh_minas, na.rm = T),
         sum_reclutamiento  = sum(lag1_cnmh_reclutamiento, lag2_cnmh_reclutamiento, lag3_cnmh_reclutamiento, na.rm = T),
         sum_ataque         = sum(lag1_cnmh_ataque_poblacion, lag2_cnmh_ataque_poblacion, lag3_cnmh_ataque_poblacion, na.rm = T),
         sum_desaparicion   = sum(lag1_cnmh_desaparicion, lag2_cnmh_desaparicion, lag3_cnmh_desaparicion, na.rm = T),
         sum_desplazamiento_pop = (sum_desplazamiento/population)*100,         
         sum_homicidio_pop      = (sum_homicidio/population)*100,
         sum_amenaza_pop        = (sum_amenaza/population)*100,
         sum_combates_pop       = (sum_combates/population)*100,
         sum_despojo_pop        = (sum_despojo/population)*100,
         sum_minas_pop          = (sum_minas/population)*100,
         sum_reclutamiento_pop  = (sum_reclutamiento/population)*100,
         sum_ataque_pop         = (sum_ataque/population)*100,
         sum_desaparicion_pop   = (sum_desaparicion/population)*100) %>%
  ungroup() %>%
  filter(!is.na(lag1_ruv_desplazamiento_forzado)) %>%
  mutate(sum_desplazamiento_log = log(sum_desplazamiento + quantile(sum_desplazamiento, .25)^2/quantile(sum_desplazamiento, .75)),
         desplazamiento_log = log(lag1_ruv_desplazamiento_forzado + quantile(lag1_ruv_desplazamiento_forzado, .25)^2/quantile(lag1_ruv_desplazamiento_forzado, .75)),
         desplazamiento_log_pop = log(lag1_ruv_desplazamiento_forzado_pop + quantile(lag1_ruv_desplazamiento_forzado_pop, .25)^2/quantile(lag1_ruv_desplazamiento_forzado_pop, .75)))
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   5. Outliers                                                                                           ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_data_final.df <- merge_data_lags.df %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         dpto = as.factor(dpto),
         month = as.factor(month),
         counter = 1,
         date = as.factor(date)) %>%
  filter(desplazamiento_log_pop != "Inf") %>%
  group_by(counter) %>%
  mutate(meanDF = mean(ruv_desplazamiento_forzado_pop, na.rm = T), 
         stdDF = sd(ruv_desplazamiento_forzado_pop, na.rm = T),
         sprayingMean = mean(spraying, na.rm = T), 
         sprayingStd = sd(spraying, na.rm = T)) %>%
  ungroup() %>%
  mutate(outliers = if_else(ruv_desplazamiento_forzado_pop > meanDF+2*stdDF, 1, 0),
         spraying_outliers = if_else(spraying > sprayingMean+2*sprayingStd, 1, 0)) %>%
  filter(outliers == 0) %>%
  #filter(spraying_outliers == 0) %>%
  distinct(codmpio, date, .keep_all = T)

saveRDS(merge_data_final.df, file = "Data/Merge/output/merge_data_comments.rds")

