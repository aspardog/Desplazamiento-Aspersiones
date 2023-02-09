## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
## Script:            Tesis Maestria
##
## Author(s):        Santiago Pardo        (santiagopardo03@gmail.com)
##
## Dependencies:      Universidad de los Andes
##
## Creation date:     February 3rd, 2023
##
## This version:      February 6th, 2023
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
p_load(ggplot2, tidyverse, ggthemes, ggpubr)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    1. Processing data and visualize                                                                     ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds")

# ================================== Scatter month to month ======================================================

scatter_month <- merge_data.df %>%
  dplyr::select(codmpio, dpto, date, ruv_desplazamiento_forzado, spraying) %>%
  filter(spraying > 0) %>%
  mutate(desplazamiento_forzado_log = log(ruv_desplazamiento_forzado),
         aspersiones_log = log(spraying)) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = desplazamiento_forzado_log)) +
  geom_point(mapping=aes(x = aspersiones_log, y = desplazamiento_forzado_log)) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = TRUE,
              size = 1, 
              fill="#69b3a2") +
  stat_cor(method = "pearson") +
  theme(panel.background   = element_blank(),
        panel.grid.major   = element_blank(),
        axis.ticks  = element_blank(),
        plot.background = element_rect(fill = "white", colour = "white"));scatter_month

ggsave(scatter_month, filename = "Visualizations/output/scatter_month.png", dpi = 320, width = 10, height = 10)

# ================================== Scatter aggregation 3 months to month =========================================

merge_data.df$ruv_desplazamiento_forzado_1 <- sapply(1:nrow(merge_data.df), function(x) merge_data.df$ruv_desplazamiento_forzado[x-1]) # Adding the first lag
merge_data.df$ruv_desplazamiento_forzado_2 <- sapply(1:nrow(merge_data.df), function(x) merge_data.df$ruv_desplazamiento_forzado[x-2]) # Adding the second lag

scatter_3month <- merge_data.df[-c(1,2,3),] %>% # Removing the observations with lag
  filter(spraying > 0) %>%
  mutate(sum_des = ruv_desplazamiento_forzado + as.numeric(ruv_desplazamiento_forzado_1) + as.numeric(ruv_desplazamiento_forzado_2),
         sum_des_log = log(sum_des),
         aspersiones_log = log(spraying)) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = sum_des_log)) +
  geom_point(mapping=aes(x = aspersiones_log, y = sum_des_log)) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = TRUE,
              size = 1,
              fill="#69b3a2") +
  stat_cor(method = "pearson") +
  theme(panel.background   = element_blank(),
        panel.grid.major   = element_blank(),
        axis.ticks  = element_blank(),
        plot.background = element_rect(fill = "white", colour = "white"));scatter_3month

ggsave(scatter_3month, filename = "Visualizations/output/scatter_3month.png", dpi = 320, width = 10, height = 10)


