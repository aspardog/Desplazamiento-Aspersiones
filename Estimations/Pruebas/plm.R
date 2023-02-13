library(tidyverse)
library(plm)
library(sandwich)
library(lmtest)

correlaciones.df <- merge_data %>%
  select(starts_with("ruv"), spraying)

a <- cor(correlaciones.df)

merge_data <- merge_data %>%
  mutate(codmpio = as.factor(codmpio),
         year = as.factor(year),
         dpto = as.factor(dpto),
         month = as.factor(month))
controles <- c("ruv_amenaza", "ruv_combates", "ruv_desaparicion_forzada", "ruv_homicidio", 
"ruv_lesiones_fis", "ruv_lesiones_psico","ruv_minas", "ruv_tortura", "ruv_violencia_sexual", 
"ruv_abandono_despojo", "year", "dpto")

reg <- lm(data = merge_data, 
          formula = paste0("ruv_desplazamiento_forzado ~ spraying + ", paste(controles, collapse = "+")))
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


fixed_effects <- plm::plm(paste0("ruv_desplazamiento_forzado ~ spraying + ", paste(controles, collapse = "+")),
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



