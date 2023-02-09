library(tidyverse, plm)
correlaciones.df <- merge_data %>%
  select(starts_with("ruv"), spraying)

a <- cor(correlaciones.df)

controles <- c("ruv_amenaza", "ruv_combates", "ruv_desaparicion_forzada", "ruv_homicidio", 
"ruv_lesiones_fis", "ruv_lesiones_psico","ruv_minas", "ruv_tortura", "ruv_violencia_sexual", "ruv_abandono_despojo")

reg <- lm(data = merge_data, formula = paste0("ruv_desplazamiento_forzado ~ spraying + ", paste(controles, collapse = "+")))
summary(reg)

fixed_effects <- plm::plm(paste0("ruv_desplazamiento_forzado ~ spraying + ", paste(controles, collapse = "+")),
                     data = merge_data, 
                     index = c("codmpio","year"),
                     model = "within",
                     effect = "twoways")
summary(fixed_effects)

# Clusterizar errores estandar a nivel de municipio
# Libreria sandwich 
# Efectos fijos, ajustar errores estandar para clusterizar a nivel municipal



