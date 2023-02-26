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
p_load(stringr, spdep, sfdep, rgdal, magrittr, ggplot2, tidyverse, sf, circlize, Matrix, spatialreg, ggthemes)

# Programming some functions

# Bivariate Moran's I
moran_I <- function(x, y = NULL, W){
  if(is.null(y)) y = x
  
  xp <- scale(x)[, 1]
  yp <- scale(y)[, 1]
  W[which(is.na(W))] <- 0
  n <- nrow(W)
  
  global <- (xp%*%W%*%yp)/(n - 1)
  local  <- (xp*W%*%yp)
  
  list(global = global, local  = as.numeric(local))
}

# Permutations for the Bivariate Moran's I
simula_moran <- function(x, y = NULL, W, nsims = 5000){
  
  if(is.null(y)) y = x
  
  n   = nrow(W)
  IDs = 1:n
  
  xp <- scale(x)[, 1]
  W[which(is.na(W))] <- 0
  
  global_sims = NULL
  local_sims  = matrix(NA, nrow = n, ncol=nsims)
  
  ID_sample = sample(IDs, size = n*nsims, replace = T)
  
  y_s = y[ID_sample]
  y_s = matrix(y_s, nrow = n, ncol = nsims)
  y_s <- (y_s - apply(y_s, 1, mean))/apply(y_s, 1, sd)
  
  global_sims  <- as.numeric( (xp%*%W%*%y_s)/(n - 1) )
  local_sims  <- (xp*W%*%y_s)
  
  list(global_sims = global_sims,
       local_sims  = local_sims)
}


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##   1. Dowloading and Processing data                                                                      ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load data
merge_data.df <- readRDS("Data/Merge/output/merge_data.rds") %>% 
  filter(year != 2003)

colombia.sf <- st_read("Data/Download/input/ShapeFiles/Municipio_ctm12.shp") %>%
  mutate(codmpio = as.numeric(mpio_ccnct)) %>%
  filter(codmpio != 88001) %>%
  filter(codmpio != 88564) 

merge_data.sf <- colombia.sf %>%
  left_join(merge_data.df, by = "codmpio") %>%
  select(codmpio, ruv_desplazamiento_forzado, spraying, dpto_ccdgo) %>%
  mutate(ruv_desplazamiento_forzado = if_else(is.na(ruv_desplazamiento_forzado), 0, ruv_desplazamiento_forzado),
         spraying = if_else(is.na(spraying), 0, spraying)) %>%
  group_by(codmpio, dpto_ccdgo) %>%
  summarise(desplazamiento_forzado = sum(ruv_desplazamiento_forzado, na.rm = T),
            aspersiones = sum(spraying, na.rm = T))

# Variables to use in the correlation: white and black population in each census track

x <- merge_data.sf$desplazamiento_forzado
y <- merge_data.sf$aspersiones

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    2. Moran Index Calculation                                                                          ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# adjancency matrix ============================================================================================

# Here I create the adjancency matrix that connect all the towns and also weight the value of the variables

# Adjacency Matrix (Queen)

nb <- poly2nb(merge_data.sf) # Create the neighbour links
spdep::is.symmetric.nb(nb)
coords <- st_coordinates(st_centroid(st_geometry(merge_data.sf)))
plot(nb, coords, col="grey") # Check that all the towns are connected

lw <- nb2listw(nb, style = "B", zero.policy = T)
lw$style

W  <- as(lw, "symmetricMatrix") # This convert the list object to a simmetry matrix for this step is necessary the spatialreg package
all(W == t(W)) # Verify the simmetry

# Moran calculation ============================================================================================

m <- moran_I(x, y, W)

# Global Moral
global_moran <- m[[1]][1]
#> 1.983673

# Local values
m_i <- m[[2]] 
# local simulations
local_sims <- simula_moran(x, y, W)$local_sims

# global pseudo p-value  
# get all simulated global moran
global_sims <- simula_moran(x, y, W)$global_sims

# Identifying the significant values 
alpha <- .2  # for a 80% confidence interval
probs <- c(alpha/2, 1-alpha/2)
intervals <- t( apply(local_sims, 1, function(x) quantile(x, probs=probs)))
sig       <- ( m_i < intervals[,1] )  | ( m_i > intervals[,2] )

#Create the map  ============================================================================================
# Preparing for plotting

# Convert shape file into sf object
map_sf     <- st_as_sf(merge_data.sf)
map_sf$sig <- sig


# Identifying the LISA clusters
xp <- scale(x)[,1]
yp <- scale(y)[,1]

val <-  as.vector(W%*%yp)

patterns <- as.character(interaction(xp > 0, val > 0))

patterns <- patterns %>% 
  str_replace_all("TRUE","High") %>% 
  str_replace_all("FALSE","Low")

patterns[map_sf$sig==0] <- "Not significant"
map_sf$patterns <- patterns


# Rename LISA clusters
map_sf$patterns2 <- factor(map_sf$patterns, levels=c("High.High", "High.Low", "Low.High", "Low.Low", "Not significant"),
                           labels=c("Mayor Desplazamiento - Mayores Aspersiones", "Mayor Desplazamiento - Menores Aspersiones", 
                                    "Menor Desplazamiento - Mayores Aspersiones", "Menor Desplazamiento - Menores Aspersiones", 
                                    "No Significativo"))

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##
##    3. Graph                                                                          ----
##
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### PLOT

a <- ggplot(map_sf) +
  geom_sf(aes(fill=patterns2), color = "grey") +
  scale_fill_manual(values = c("#C1121F", "#FDF0D5", "#9CBED3", "#005D8F", "#E2E4E6")) + 
  guides(fill = guide_legend(title="LISA clusters")) +
  coord_sf(expand = FALSE) +
  theme_map() +
  theme(legend.key.size = unit(1, "cm"),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.line = element_blank(),
        legend.text=element_text(size=10),
        legend.title = element_text(size = 12, face = "bold"),
        plot.background = element_rect(fill = "white", colour = "white"));a

ggsave(a, filename = "Visualizations/output/moran_bivariado.png", dpi = 320, width = 10, height = 10)


