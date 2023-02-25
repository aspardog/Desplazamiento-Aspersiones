library(tidyverse)
library(plm)
library(sandwich)
library(lmtest)
library(ivreg)


correlaciones.df <- merge_data.df %>%
  select(starts_with("ruv"), spraying, starts_with("cnmh"), night_lights)

a <- cor(correlaciones.df)

merge_data <- merge_data.df %>%
  filter(year > 2003 & year < 2013) %>%
  mutate(ydep = log(ruv_desplazamiento_forzado)) %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         dpto = as.factor(dpto),
         month = as.factor(month))
controles <- c("ruv_amenaza", "ruv_desaparicion_forzada",
               "ruv_minas", "ruv_tortura",  "cnmh_violencia_sexual",
               "ruv_lesiones_fis", "ruv_lesiones_psico", "ruv_abandono_despojo", 
               "cnmh_atentado","cnmh_masacres", "cnmh_asesinatos_selectivos",
               "cnmh_enfrentamientos", "night_lights", "year", "codmpio")

reg <- lm(data = merge_data, 
          formula = paste0("log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75)) ~  windIV20MERRA +", paste(controles, collapse = "+")))

summary(reg)
m1coeffs_std <- data.frame(summary(reg)$coefficients)
coi_indices <- which(!startsWith(row.names(m1coeffs_std), 'dpto'))
m1coeffs_std[coi_indices,]

m1coeffs_cl <- coeftest(reg, vcov = vcovCL, cluster = ~codmpio)
m1coeffs_cl[coi_indices,]
(m1cis <- coefci(reg, parm = coi_indices, vcov = vcovCL,
                cluster = ~codmpio, level = 0.9))


# Second prueba
merge_data <- merge_data %>%
  mutate(fe = paste0(dpto, year))


fixed_effects <- plm::plm(paste0("log(ruv_desplazamiento_forzado + min(ruv_desplazamiento_forzado[ruv_desplazamiento_forzado > 0])/2) ~ windSpeedFLDAS + ", paste(controles, collapse = "+")),
                     data = merge_data, 
                     index = c("year"),
                     model = "within",
                     effect = "twoways")
summary(fixed_effects)

vcov <- fixed_effects$vcov

m1coeffs_cl <- coeftest(fixed_effects, vcov = vcov, cluster = ~codmpio)
m1coeffs_cl[coi_indices,]

summary(m1coeffs_cl)

# Clusterizar errores estandar a nivel de municipio
# Libreria sandwich 
# Efectos fijos, ajustar errores estandar para clusterizar a nivel municipal
prueba <- lm(formula = ruv_desplazamiento_forzado ~ windIV20MERRA, data = merge_data)
instrumento <- lm(formula = spraying ~  windIV20MERRA, data = merge_data)
summary(instrumento)

a <- ivreg::ivreg(formula = paste0("log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~ spraying | windIV20MERRA + windSpeedFLDAS +", paste(controles, collapse = "+")), 
                  data = merge_data)
summary(a)
