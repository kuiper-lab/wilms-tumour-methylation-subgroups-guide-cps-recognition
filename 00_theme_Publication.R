# ==============================================================================
# Script: Publication Theme
# Author: Franziska Oberhammer
# Date: 2026-09-15
#
# Description:
# Custom ggplot2 theme and colour scales for publication-quality figures.
# The theme uses Helvetica font, minimal backgrounds and borders, and a
# consistent colour palette for fill and colour aesthetics.
# Adapted from Koundinya Desiraju, (https://rpubs.com/koundy/71792)
#
# Functions:
#   theme_Publication().       - Custom publication-ready ggplot2 theme
#   scale_fill_Publication()   - Custom discrete fill colour scale
#   scale_colour_Publication() - Custom discrete colour scale
#
# Packages:
#   ggplot2
#   ggthemes
#   scales
#   grid
#
# Usage:
#   ggplot(data, aes(x, y, fill = group)) +
#     geom_bar(stat = "identity") +
#     cotheme_Publication() +
#     scale_fill_Publication()
# ==============================================================================


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
(theme_foundation(base_size = base_size) +
  theme(
    text = element_text(family = "Helvetica"),
    
    plot.title = element_text(
      family = "Helvetica",
      face = "bold",
      size = 7,
      hjust = 0.5
    ),
    
    # no background, no border
    panel.background = element_blank(),
    plot.background  = element_blank(),
    legend.key       = element_blank(),
    panel.border     = element_blank(),
    
    axis.title = element_text(
      family = "Helvetica",
      face = "bold",
      size = 7
    ),
    axis.title.y = element_text(
      family = "Helvetica",
      face = "bold",
      size = 7,
      angle = 90,
      vjust = 2
    ),
    axis.title.x = element_text(
      family = "Helvetica",
      face = "bold",
      size = 7,
      vjust = -0.2
    ),
    
    axis.text = element_text(
      family = "Helvetica",
      face = "plain",
      size = 6
    ),
    
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(),
    panel.grid.major = element_line(colour = "#f0f0f0"),
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.size = unit(0.2, "cm"),
    legend.margin = unit(0, "cm"),
    
    legend.title = element_text(
      family = "Helvetica",
      face = "plain",
      size = 7
    ),
    legend.text = element_text(
      family = "Helvetica",
      face = "plain",
      size = 6
    ),
    
    plot.margin = unit(c(10, 5, 5, 5), "mm"),
    strip.background = element_rect(colour = NA, fill = NA),
    strip.text = element_text(
      family = "Helvetica",
      face = "bold"
    )
  ))
}


scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}