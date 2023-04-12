## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     February 12th, 2023
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
p_load(sf, tmap, RColorBrewer, cartogram, geodaData, tidyverse, ggthemes)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Dowloading and Processing data                                                                      ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load data
merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") %>%
  mutate(codmpio = as.character(codmpio),
         mpio_area = mpio_area*100, # Trasnformación a hectareas
         coca_norm = round((cultivos/mpio_area)*10000,2)
  ) %>%
  group_by(codmpio) %>%
  summarise(aspersiones = mean(spraying_norm, na.rm = T),
            coca = mean(coca_norm, na.rm = T)) 

colombia.sf <- st_read("Data/Download/input/ShapeFiles/Municipio_ctm12.shp") %>%
  mutate(codmpio = as.character(mpio_ccnct)) %>%
  filter(codmpio != 88001) %>%
  filter(codmpio != 88564) 

merge_data.sf <- colombia.sf %>%
  left_join(merge_data.df, by = "codmpio") %>%
  select(codmpio, aspersiones, coca) %>%
  mutate(aspersiones = if_else(aspersiones == 0 , NA_real_, aspersiones),
         aspersiones = round(aspersiones*100,2),
         coca = round(coca,2)
         )

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    2. Graph                                                                          ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
mapas <- function(data,variable, grupo){
  
  mapa <- data %>%
    group_by({{grupo}}) %>%
    dplyr::summarize(counter = sum({{variable}}, na.rm = T))
  
  mapa <- mapa %>%
    mutate(breaks = cut(counter, b = unique(quantile(counter, probs= seq(0,1,1/35), na.rm = T)), dig.lab = 7), include.lowest = TRUE) %>%
    arrange(counter)
  
  col <- 5
  fichas <- colorRampPalette(brewer.pal(9, "BuGn"))(col)
  
  breaks2 <- unique(pull(mapa,breaks))
  
  breaks2 <- na.omit(breaks2)
  
  a <- ggplot(mapa) +
    geom_sf(aes(fill = breaks),color ="black",size = 0.1)+
    scale_fill_manual(values = fichas,
                      limits = breaks2, na.value = "white") +
    theme_map()+
    labs(fill = "Intensidad de hectáreas asperjadas",
         caption = "Fuente: DIRAN") +
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
a <- mapas(data = merge_data.sf, variable = aspersiones, grupo = codmpio);a
b <- mapas(data = merge_data.sf, variable = coca, grupo = codmpio);b
### PLOT

ggsave(a, filename = "Visualizations/output/aspersionesNorm.png", dpi = 320, width = 10, height = 10)


