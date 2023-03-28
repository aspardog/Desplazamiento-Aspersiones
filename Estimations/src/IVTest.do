 tab year, gen(dr)
 tab month, gen(dm)
 tab codmpio, gen(dmpio)
 tab date, gen(date)
 
 xtset codmpio date

  xtivreg2 desplazamiento_log (spraying = windSpeedRMBOS) vegetation cultivos night_lights rainFall lag1_ruv_combates lag1_ruv_abandono_despojo lag1_cnmh_minas lag1_cnmh_reclutamiento lag1_ruv_homicidio, fe robust cluster(codmpio)
  weakivtest, level(0.1)

 
 ivreg2 desplazamiento_log vegetation cultivos night_lights rainFall ruv_combates ruv_abandono_despojo ruv_desaparicion_forzada cnmh_minas cnmh_reclutamiento ruv_homicidio  (spraying = windSpeedRMBOS windIV30RMBOS) i.codmpio i.date, robust cluster(codmpio)
  weakivtest, level(0.1)

  ivreg2 sum_desplazamiento_log cultivos vegetation sum_combates sum_despojo sum_minas sum_reclutamiento sum_homicidio night_lights  dm* dr* dmpio* (spraying = windSpeedRMBOS), robust cluster(codmpio)
  weakivtest, level(0.1)
 

 
**#
ivreg2 desplazamiento_log ruv_amenaza ruv_combates ruv_abandono_despojo ruv_homicidio  cnmh_minas cnmh_reclutamiento cnmh_ataque_poblacion  ivreg2 desplazamiento_log ruv_amenaza ruv_combates ruv_abandono_despojo ruv_homicidio  cnmh_minas cnmh_reclutamiento cnmh_ataque_poblacion night_lights  i.month i.codmpio i.year (spraying = windSpeedR rainFall), robust cluster(codmpio)
  weakivtest, level(0.1)
  
 night_lights  i.month i.codmpio i.year (spraying = windSpeedFLDAS rainFall dist2airport), robust cluster(codmpio)
  weakivtest, level(0.1)

  
  
xtdata desplazamiento_log ruv_amenaza ruv_combates ruv_abandono_despojo ruv_homicidio  cnmh_minas cnmh_reclutamiento cnmh_ataque_poblacion night_lights spraying windSpeedRMBOS rainFall, i(codmpio month year) fe clear 
 ivreg2 desplazamiento_log ruv_amenaza ruv_combates ruv_abandono_despojo ruv_homicidio  cnmh_minas cnmh_reclutamiento cnmh_ataque_poblacion night_lights (spraying = windSpeedRMBOS rainFall), robust cluster(codmpio)
  weakivtest, level(0.1)
