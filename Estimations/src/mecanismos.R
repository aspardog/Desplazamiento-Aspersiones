## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Download                                                                                            ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
controles_fe_3month <- c('night_lights', 
                         "rainFall","vegetation", 'windIV10RMBOS',
                         'sum_combates_pop', 'sum_despojo_pop', 'sum_minas_pop', 'sum_reclutamiento_pop', 
                         'sum_homicidio_pop', 'sum_desaparicion_pop')

merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") 

counter <- merge_data.df %>%
  filter(spraying > 0) %>%
  mutate(counter = 1) %>%
  group_by(codmpio) %>%
  summarise(counter_spraying = sum(counter, na.rm = T))

mergeCounter <- merge_data.df %>%
  left_join(y = counter, by = "codmpio") %>%
  mutate(intensidad = ntile(counter_spraying, n = 4),
         intensidad = as.factor(intensidad),
         sum_desplazamiento_pop = sum_desplazamiento_pop/10)

IV3MonthFE <- plm(data= mergeCounter,
                  formula = 
                    paste0("sum_desplazamiento_pop ~ spraying_norm*intensidad + spraying_norm + intensidad +",
                           paste(controles_fe_3month, collapse = "+"),"|",
                           paste("windSpeedRMBOS*intensidad + spraying_norm + intensidad +"),
                           paste(controles_fe_3month, collapse = "+")),
                  effect = "twoways",
                  model = "within",
                  index=c("codmpio", "date", "query"))
summary(IV3MonthFE)
IVFE3Mstd <- coeftest(IV3MonthFE, vcov=vcovHC(IV3MonthFE, type="sss", cluster="group"))  
IVASPHEIC <- confint(IV3MonthFE, level = 0.8)

estimate <- as.data.frame(IVFE3Mstd[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:intensidad2", "spraying_norm:intensidad3", "spraying_norm:intensidad4"))
CI <- as.data.frame(IVASPHEIC[,]) %>%
  rownames_to_column(var = "variable") %>%
  filter(variable %in% c("spraying_norm:intensidad2", "spraying_norm:intensidad3", "spraying_norm:intensidad4")) %>%
  rename("lower" = "10 %", "upper" = "90 %")

data2plot <- estimate %>%
  inner_join(CI, by = "variable") %>%
  mutate(
    order_value =
      case_when(
        variable == "spraying_norm:intensidad2" ~ 1,
        variable == "spraying_norm:intensidad3" ~ 2,
        variable == "spraying_norm:intensidad4" ~ 3,
      ),
    variable = 
      case_when(
        variable == "spraying_norm:intensidad2" ~ "Segundo quintil",
        variable == "spraying_norm:intensidad3" ~ "Tercer quintil",
        variable == "spraying_norm:intensidad4" ~ "Cuarto quintil"
      )
  )

MPlot <- ggplot(data2plot, aes(x = reorder(variable, order_value), y = Estimate)) +
  geom_hline(yintercept = 0, lty = 1, color = "#fa4d57", lwd = 1)  +
  geom_linerange(aes(x = reorder(variable, order_value),  ymin = lower, ymax = upper),
                 lwd = 0.5, position = position_dodge(width = .7), 
                 stat = "identity", color = "#003b8a")+
  geom_point(aes(x =reorder(variable, order_value), y = Estimate), 
             size = 3.5, position = position_dodge(width = 1), color = "#003b8a") +
  geom_point(aes(x = reorder(variable, order_value), y = Estimate), 
             size = 2.5, position = position_dodge(width = 1), color = "white") +
  scale_y_continuous(limits = c(0, 5),
                     breaks = seq(0, 5, by = 1),
                     expand = expansion(mult = 0.025), position = "left",
                     labels = c("0", "1","2", "3", "4", "5")) +
  coord_flip() +
  labs(x = "Cuartiles de aplicación del programa",
       y = "Efecto aspersiones sobre desplazamiento") +
  theme(panel.background   = element_blank(),
        plot.background    = element_blank(),
        panel.grid.major   = element_line(size     = 0.25,
                                          colour   = "#5e5c5a",
                                          linetype = "dashed"),
        panel.grid.minor   = element_blank(),
        axis.ticks  = element_blank(),
        plot.margin  = unit(c(0, 0, 0, 0), "points")) +
  theme(legend.position = "none",
        panel.background   = element_blank(),
        panel.grid.major.x = element_line(colour = "#d1cfd1", 
                                          size = 0.5, linetype = "dashed"),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 12, margin   = margin(0, 20, 0, 10), vjust = 1),
        axis.title.x = element_text(size = 12, margin   = margin(20, 0, 10, 0), vjust = 0),
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 10),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.minor.y = element_blank(),
        ggh4x.axis.ticks.length.minor = rel(1),
        axis.line.x.bottom = element_line(linetype = "solid", size = 1));MPlot 
ggsave(MPlot, filename = "Visualizations/output/Mecanismos.png", dpi = 320, width = 7.5, height = 7.5)

stargazer(
  IV3MonthFE, 
  type = "latex",
  dep.var.labels= "Aspersiones aéreas",
  dep.var.caption = "Variable dependiente: Desplazamiento Forzado",
  keep = "spraying_norm",
  se = list(IVFE3Mstd[, 2]), 
  title = "Efecto agregado de las aspersiones aéreas sobre el desplazamiento forzado", 
  align = TRUE, 
  dep.var.labels.include = FALSE, 
  no.space = FALSE, 
  covariate.labels = c("Aspersiones aéreas"), 
  omit = "Constant",
  add.lines = list(c("Efectos Fijos", "Si"),
                   c("Controles", "Si")),
  #column.labels = c("MCO"),
  column.sep.width = "5pt",
  keep.stat = c("rsq", "n"), 
  p = list(IVFE3Mstd[, 4]),
  decimal.mark = ",",
  notes.align = "l",
  notes.append = F
)

