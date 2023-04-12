## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     February 25th, 2023
##
## This version:      February 25th, 2023
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Visualizations                                                                                       ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    Set - UP                                                                                             ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(pacman)
p_load(sf, sp, ggplot2, tidyverse, lubridate, RColorBrewer)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Dowloading and Processing data                                                                      ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ruv.df <- readRDS("Data/Download/output/Victimization/ruv_victimization.rds") %>%
  mutate(year = as.numeric(year(ym(date)))) %>%
  filter(year > 2003 & year < 2013)

colombia.sf <- st_read("Data/Download/input/ShapeFiles/Municipio_ctm12.shp") %>%
  mutate(codmpio = as.numeric(mpio_ccnct)) %>%
  filter(codmpio != 88001) %>%
  filter(codmpio != 88564) 

population.df    <- readRDS("Data/Download/output/Population/population.rds") %>%
  mutate(codmpio = as.numeric(codmpio))

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    2. Calculation                                                                         ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

desplazamiento.df <- ruv.df %>%
  group_by(codmpio) %>%
  summarise(desplazamientos = sum(ruv_desplazamiento_forzado, na.rm = T))

data2plot <- colombia.sf %>%
  left_join(desplazamiento.df, by = "codmpio") %>%
  select(codmpio, desplazamientos)

desplazamientos_pop <- ruv.df %>%
  left_join(population.df, by = "codmpio") %>%
  filter(population != 0) %>%
  mutate(desplazamiento_pop = (ruv_desplazamiento_forzado/population)*100) %>%
  group_by(codmpio) %>%
  summarise(desplazamientos = mean(desplazamiento_pop, na.rm = T)) 

desplazamientos_pop  <- colombia.sf %>%
  left_join(desplazamientos_pop, by = "codmpio")


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    2. Graph                                                                          ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mapas <- function(data, variable, grupo){
  
  mapa <- data %>%
    group_by({{grupo}}) %>%
    dplyr::summarize(counter = sum({{variable}}*100, na.rm = T))
  
  mapa <- mapa %>%
    mutate(breaks = cut(counter, b = unique(round((quantile(counter, probs= seq(0,1,1/5), na.rm = T)),0)), dig.lab = 7), include.lowest = TRUE) %>%
    arrange(counter)
  
  col <- 5
  fichas <- colorRampPalette(brewer.pal(9, "OrRd"))(col)
  
  breaks2 <- unique(pull(mapa,breaks))
  
  breaks2 <- na.omit(breaks2)
  
  a <- ggplot(mapa) +
    geom_sf(aes(fill = breaks),color ="black",size = 0.1)+
    scale_fill_manual(values = fichas,
                      limits = breaks2, na.value = "white") +
    theme_map()+
    labs(fill = "Desplazamientos Forzados",
         caption = "Fuente: RUV") +
    theme(legend.key.size = unit(1, "cm"),
          legend.position = "right",
          legend.key = element_blank(),
          legend.background = element_blank(),
          axis.line = element_blank(),
          legend.title=element_text(size=12,face="bold"), 
          legend.text=element_text(size=10),
          plot.title=element_text(size=14, hjust=0.8, vjust=1, face="bold"),
          plot.subtitle=element_text(size=10, hjust=0.6, vjust=1),
          plot.caption=element_text(size=8, hjust=0.01))
  
  return(a)
}

a <- mapas(data = data2plot, variable = desplazamientos, grupo = codmpio);a 

ggsave(a, filename = "Visualizations/output/desplazamientos.png", dpi = 320, width = 10, height = 10)

b <- mapas(data = desplazamientos_pop, variable = desplazamientos, grupo = codmpio);b
ggsave(b, filename = "Visualizations/output/desplazamientosNorm.png", dpi = 320, width = 10, height = 10)
