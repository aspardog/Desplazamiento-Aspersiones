## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     March 4th, 2023
##
## This version:      March 4th, 2023
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
p_load(sf, tmap, RColorBrewer, cartogram, geodaData, tidyverse, ggthemes, lubridate, patchwork)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Dowloading and Processing data                                                                      ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load data
merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") %>%
  filter(year < 2013) %>%
  filter(year > 2003)


colombia.sf <- st_read("Data/Download/input/ShapeFiles/Municipio_ctm12.shp") %>%
  mutate(codmpio = as.numeric(mpio_ccnct)) %>%
  filter(codmpio != 88001) %>%
  filter(codmpio != 88564) 

merge_data.sf <- colombia.sf %>%
  left_join(merge_data.df, by = "codmpio") %>%
  select(codmpio, spraying, date, windSpeedFLDAS) %>%
  mutate(aspersiones = if_else(spraying == 0 , NA_real_, spraying),
         aspersiones = round(aspersiones, 0)) %>%
  filter(aspersiones > 0) %>%
  arrange(date) %>%
  group_by(date) %>%
  summarise(intensidad = mean(aspersiones, na.rm = T), windSpeed = mean(windSpeedFLDAS, na.rm = T)) %>%
  st_drop_geometry()

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   2. Graph                                                                     ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data2plot <- merge_data.sf %>%
  mutate(year= as.numeric(year(ym(date))),
         date = parse_date_time(date, "ym")) %>%
  group_by(year) %>%
  summarise(intensidad = mean(intensidad), windSpeed = mean(windSpeed)) %>%
  mutate(logIntensidad = log(intensidad + quantile(intensidad, .25)^2/quantile(intensidad, .75)),
         logWindSpeed  = log(windSpeed + quantile(windSpeed, .25)^2/quantile(windSpeed, .75))) %>%
  rename(Año = year)

# Option 1

coeff <- 5

a <- ggplot(data2plot, aes(x=Año)) +
  
  geom_line( aes(y=logWindSpeed, linetype = "Velocidad del Viento"), color = "#427AA1", show.legend = T) + 
  geom_point(aes(y=logWindSpeed), color = "#427AA1") +
  geom_line( aes(y=logIntensidad/5, linetype = "Aspersiones Aéreas"), color = "#679436", show.legend = T) +
  geom_point(aes(y=logIntensidad/5), color = "#679436") +
  

  scale_y_continuous(
    
    # Features of the first axis
    name = "Log(Velocidad del Viento)", breaks = seq(0, 2, by = 0.1),
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Log(Aspersiones aereas)", breaks = seq(6, 8, by = 0.5)), 
    expand = expansion(mult = c(0.25, 0.05))
  ) +
  scale_x_continuous(breaks = seq(2004, 2012, by = 1),limits = c(2004,2012)) +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.title.y       = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C",
                                          margin   = margin(0, 0, 0, 10)),
        axis.title.x       = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C",
                                          margin   = margin(0, 0, 0, 10)),
        axis.text.y        = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C"),
        axis.text.x = element_text(face   = "plain",
                                   size   = 3.514598*.pt,
                                   color  = "#524F4C"),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#d1cfd1"),
        axis.line.x        = element_line(color    = "#d1cfd1"), 
        legend.position = "top", 
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent"), legend.title = element_blank()
        );a

ggsave(a, filename = "Visualizations/output/viento_aspersion.png", dpi = 320, width = 7.5, height = 5)

# Option 2

plt1 <- ggplot(data2plot, aes(x = year)) +
  geom_line( aes(y=windSpeed), color = "blue") + 
  geom_point(aes(y=windSpeed), color = "blue") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.title.y       = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C",
                                          margin   = margin(0, 0, 0, 10)),
        axis.title.x       = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C",
                                          margin   = margin(0, 0, 0, 10)),
        axis.text.y        = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C"),
        axis.text.x = element_text(face   = "plain",
                                   size   = 3.514598*.pt,
                                   color  = "#524F4C"),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#d1cfd1"),
        axis.line.x        = element_line(color    = "#d1cfd1"), legend.position = "bottom");plt1

plt2 <- ggplot(data2plot, aes(x = year)) +
  geom_line( aes(y=intensidad), color = "blue") + 
  geom_point(aes(y=intensidad), color = "blue") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.title.y       = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C",
                                          margin   = margin(0, 0, 0, 10)),
        axis.title.x       = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C",
                                          margin   = margin(0, 0, 0, 10)),
        axis.text.y        = element_text(face     = "plain",
                                          size     = 3.514598*.pt,
                                          color    = "#524F4C"),
        axis.text.x = element_text(face   = "plain",
                                   size   = 3.514598*.pt,
                                   color  = "#524F4C"),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#d1cfd1"),
        axis.line.x        = element_line(color    = "#d1cfd1"), legend.position = "bottom");plt2

plt1 + plt2  
  

