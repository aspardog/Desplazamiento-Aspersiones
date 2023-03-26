model_plm <- plm(formula =  desplazamiento_log ~
                   spraying + vegetation + cultivos + night_lights + ruv_combates + ruv_abandono_despojo + cnmh_minas + cnmh_reclutamiento + ruv_homicidio| 
                   windSpeedRMBOS + rainFall + vegetation + cultivos + night_lights  + ruv_combates + ruv_abandono_despojo + cnmh_minas + cnmh_reclutamiento + ruv_homicidio, 
                 data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(model_plm)
coeftest(model_plm, vcov=vcovHC(model_plm, type="sss", cluster="group"))  


model_plm3 <- plm(formula =  sum_desplazamiento_log ~
                    spraying + vegetation + cultivos + night_lights + sum_combates + sum_despojo + sum_minas + sum_reclutamiento + sum_homicidio| 
                    windSpeedRMBOS + rainFall + vegetation + cultivos + night_lights +  sum_combates + sum_despojo + sum_minas + sum_reclutamiento + sum_homicidio,
                 data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(model_plm3)

coeftest(model_plm3, vcov=vcovHC(model_plm3, type="sss", cluster="group"))  


model_plm <- plm(formula =  desplazamiento_log ~
                   spraying + vegetation + cultivos + night_lights + ruv_amenaza + ruv_combates + ruv_abandono_despojo + cnmh_minas + cnmh_reclutamiento + ruv_homicidio| 
                   windSpeedFLDAS + windSpeedFLDAS*dpto + rainFall + rainFall*dpto + vegetation + cultivos + night_lights + ruv_amenaza + ruv_combates + ruv_abandono_despojo + cnmh_minas + cnmh_reclutamiento + ruv_homicidio, 
                 data= merge_data.df, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(model_plm)
a <- coeftest(model_plm, vcov=vcovHC(model_plm, type="sss", cluster="group"))  


models <- list(
  "(MCO)"    = coeftest(model_plm, vcov=vcovHC(model_plm, type="sss", cluster="group"))  ,
  "(IV)"     = coeftest(model_plm3, vcov=vcovHC(model_plm, type="sss", cluster="group"))  
)

coef_rename_ws <- "Aspersiones Aéreas \ncon Glifosato"

gm <- tibble::tribble(
  ~raw,        ~clean,          ~fmt,
  "nobs",      "N",             0,
  "r.squared", "R<sup>2</sup>", 2,)
new_rows <- tibble::tribble(
  ~term,~"(MCO)",~"(IV)",
  "Controles","Sí","Sí",
  "Efectos Fijos", "Sí","Sí")

modelsummary(models, estimate = "{estimate}{stars}", ,
             coef_omit = c(-1), stars = c('*' = .1, '**' = .05, '***' = 0.01), 
             output = "markdown", fmt = 5, coef_rename = coef_rename_ws, 
             gof_map = gm, add_rows = new_rows, title = "Aspersiones aéreas", 
             notes = "*** p<0.01, ** p<0.05, * p<0.1.Los errores estándar son robustos y están corregidos por clusters de hogar.")



model_plm <- plm(formula =  log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)) ~
                   spraying*dpto + vegetation + night_lights + rainFall + ruv_amenaza + ruv_combates + ruv_abandono_despojo + ruv_homicidio + cnmh_minas + cnmh_reclutamiento | windSpeedFLDAS*dpto + rainFall + vegetation + night_lights + ruv_amenaza + ruv_combates + ruv_abandono_despojo + ruv_homicidio + cnmh_minas + cnmh_reclutamiento , 
                 data= prueba, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(model_plm)

prueba <- merge_data.df %>%
  distinct(codmpio, date, .keep_all = T)
    
model_plm <- plm(formula =  ruv_desplazamiento_forzado ~
                   spraying + cultivos + night_lights + ruv_amenaza + ruv_combates + ruv_abandono_despojo + ruv_homicidio + cnmh_minas + cnmh_reclutamiento | windSpeedFLDAS + cultivos + + rainFall + night_lights + ruv_amenaza + ruv_combates + ruv_abandono_despojo + ruv_homicidio + cnmh_minas + cnmh_reclutamiento , 
                 data= prueba, effect = "twoways", model = "within", index=c("codmpio", "date"))
summary(model_plm)
