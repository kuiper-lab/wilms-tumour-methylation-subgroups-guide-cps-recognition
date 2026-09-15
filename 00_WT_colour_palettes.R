
# Colour palettes for WT project
# 15.09.2026
# Franziska Oberhammer

colours_WT <- list()

# Colour scale for histology
colours_WT$histology <- c(Blastemal = "#E69F00",
                          Completely_necrotic = "#56B4E9",
                          DA = "#009E73", 
                          Epithelial = "#D3D3D3", 
                          Mixed = "#0072B2",
                          `Non-anaplastic` = "#D55E00",
                          Non_anaplastic = "#D55E00",
                          Regressive = "#F0E442",
                          Stromal = "#000000", 
                          WT_with_rhabdomyoblastic_differentiation = "#CA16FD")


# Colour scale for histology
colours_WT$histology_colors <- c(
  "Completely_necrotic" = "#cc4400",
  "Regressive" = "#d66915",
  "DA" = "#e08e29",
  "Non_anaplastic" = "#f0c761",
  "Mixed" = "#ffff99",
  "Epithelial" = "#c2fcff",
  "Blastemal" = "#3890bc",
  "Stromal" = "#1c489a"
  # NA will be grey by default if not assigned
)

# Colour scale for beta values
colours_WT$meth <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101)

# Colour palette for hyper/hypomethylation in IC1 and IC2
colours_WT$IC1_IC2 = c( "hyper" = "#A31621",
                        "pUPD" = "#29335C",
                        "intermediate" = "#FFD131",
                        "normal" = "#999999")

# Colour palette for methylation clusters
colours_WT$DNAm_clustering2 = c(`1a` = "#8C1834",
                               `1b` = "#D1244C",
                               `2a` = "#87BBD8",
                               `2b` = "#34658A",
                               `2c` = "#2F4858")

# Color palette for BWSp categories
colours_WT$bwsp_colors <- c(
  "No" = "#E7EFEB",
  "Mosaic" = "#A5CDBB",
  "Clinical only" = "#7DA192",
  "Yes" = "#517164"
)


