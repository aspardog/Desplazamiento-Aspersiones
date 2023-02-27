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
p_load(ggplot2, tidyverse, ggthemes, ggpubr, showtext, patchwork)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    1. Processing data and visualize                                                                     ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds")

# ================================== Scatter month to month ======================================================

scatter_month <- merge_data.df %>%
  filter(year != 2003) %>%
  dplyr::select(codmpio, dpto, date, ruv_desplazamiento_forzado, spraying) %>%
  filter(spraying > 0) %>%
  mutate(desplazamiento_forzado_log = log(ruv_desplazamiento_forzado + quantile(ruv_desplazamiento_forzado, .25)^2/quantile(ruv_desplazamiento_forzado, .75)),
         aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = desplazamiento_forzado_log)) +
  geom_point(mapping=aes(x = aspersiones_log, y = desplazamiento_forzado_log)) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = TRUE,
              size = 1, 
              fill="#69b3a2") +
  stat_cor(method = "pearson") +
  labs(subtitle = "Correlación mensual entre el desplazamiento forzado \ny aspersión aérea",
       x ="Logaritmo de aspersiones con glifosato", 
       y = "Logaritmo de desplazamientos forzados") +
  theme(panel.background   = element_blank(),
        panel.grid.major   = element_blank(),
        axis.ticks  = element_blank(),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
        axis.title = element_text(size = 8, margin   = margin(10, 20, 20, 0)),
        plot.background = element_rect(fill = "white", colour = "white"))

# ================================== Scatter aggregation 3 months to month =========================================

merge_data.df$ruv_desplazamiento_forzado_1 <- sapply(1:nrow(merge_data.df), function(x) merge_data.df$ruv_desplazamiento_forzado[x-1]) # Adding the first lag
merge_data.df$ruv_desplazamiento_forzado_2 <- sapply(1:nrow(merge_data.df), function(x) merge_data.df$ruv_desplazamiento_forzado[x-2]) # Adding the second lag

scatter_3month <- merge_data.df[-c(1,2,3),] %>% # Removing the observations with lag
  filter(spraying > 0) %>%
  mutate(sum_des = ruv_desplazamiento_forzado + as.numeric(ruv_desplazamiento_forzado_1) + as.numeric(ruv_desplazamiento_forzado_2),
         sum_des_log = log(sum_des + quantile(sum_des, .25)^2/quantile(sum_des, .75)),
         aspersiones_log = log(spraying + quantile(spraying, .25)^2/quantile(spraying, .75))) %>%
  ggplot(data = ., aes(x = aspersiones_log, y = sum_des_log)) +
  geom_point(mapping=aes(x = aspersiones_log, y = sum_des_log)) +
  stat_smooth(method = "lm",
              col = "#C42126",
              se = TRUE,
              size = 1,
              fill="#69b3a2") +
  stat_cor(method = "pearson") +
  labs(subtitle = "Correlación entre el desplazamiento trimestral \nacumulado y aspersión aérea con glifosato*", 
       x ="Logaritmo de aspersiones con glifosato", 
       y = "Logaritmo de desplazamientos forzados acumulados") +
  theme(panel.background   = element_blank(),
        panel.grid.major   = element_blank(),
        axis.ticks  = element_blank(),
        plot.subtitle = element_text(size = 12, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 8, margin   = margin(10, 20, 0, 0)),
        axis.title = element_text(size = 8, margin   = margin(10, 20, 0, 0)),
        plot.background = element_rect(fill = "white", colour = "white"))

figures <- list()
figures[["Panel A"]] <- scatter_month
figures[["Panel B"]] <- scatter_3month

figureScatter <- figures[["Panel A"]] + plot_spacer() + figures[["Panel B"]]  + 
  plot_layout(ncol = 3, nrow = 1,
              widths = unit(c(10,2,10), "cm"),
              heights = unit(15, "cm"))+ 
  plot_annotation(caption = "*El primer mes de desplazamiento forzado acumulado se cuenta a partir del mes en el que asperjó",
                  theme = theme(plot.caption = element_text(hjust = 0, size = 8, margin = margin(20, 0, 0, 0))));figureScatter

ggsave(figureScatter, filename = "Visualizations/output/Scatter.png", dpi = 320, width = 10, height = 7.5)
