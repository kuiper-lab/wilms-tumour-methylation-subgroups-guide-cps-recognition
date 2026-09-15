# Colour palettes for somatic profiling WT project

colours_WT <- list()

# meta data colour palettes ---------------------------------------------------------------------
# colour palette Sample_Group
colours_WT$sampleGroup <- c(No_genetic_predisposition = "#999999",
                            BWSp = "#E69F00",
                            BWSp_clinical_only = "#E66C00",
                            WT1_mutation = "#56B4E9",
                            DIS3L2_mutation = "#009E73", 
                            Other = "black")

# colour palette Sample_Group (with names as strings, as some plots need the vector to be like this)
colours_WT$sampleGroup_str <- c("No_genetic_predisposition" = "#999999",
                                "BWSp" = "#E69F00",
                                "BWSp_clinical_only" = "#E66C00",
                                "WT1_mutation" = "#56B4E9",
                                "DIS3L2_mutation" = "#009E73", 
                                "Other" = "black")

# colour scale histology
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


colours_WT$Sex <- c("Female" = "firebrick",
                    "Male" = "royalblue3")


# general colour scales ------------------------------------------------------------------------------------------
colours_WT$yesNo <- c("Yes" = "darkgrey", "No" = "lightgrey")

# colour blind palette that can be used for everything
colours_WT$general <- c("#E69F00",
                        "#56B4E9",
                        "#009E73", 
                        "#D3D3D3", 
                        "#0072B2",
                        "#D55E00",
                        "#F0E442",
                        "#000000")

colours_WT$tumourContent <- colorRamp2(c(0,50, 100), c("red", "yellow", "blue"))


# colour scale data availability plot------------------------------------------------------------------------------------------
colours_WT$cbPalette_DataAvailabilityPlot <- c(No_genetic_predisposition = "#999999",
                                               BWSp = "#E69F00",
                                               BWSp_clinical_only = "#E66C00",
                                               WT1_mutation = "#56B4E9",
                                               DIS3L2_mutation = "#CC79A7", 
                                               Other = "black",
                                               Prospective = "#FFCA33",
                                               Retrospective = "#FF5733",
                                               `FALSE` = "#D6D7D6", 
                                               `TRUE` = "#298F49",
                                               No = "#D6D7D6", 
                                               Yes = "#298F49",
                                               Sophie1 = "#581845", 
                                               Sophie2  = "#C70039",
                                               Diagnostics = "#FF7233",
                                               unknown = "pink") 



# colour palette for CNV information----------------------------------------------------------------------------------------------------
colours_WT$CNV <- c("high gain" = "#196F3D",
                    "gain" = "#4169E1",
                    "mild gains" = "#7DCEA0",
                    "gain?" = "#D4EFDF",
                    "gain_artefact?" = "#E9F7EF",
                    "gain_loss" = "#EFDE50",
                    "gain_loss?" = "#FAF4C3",
                    "loss_artefact?" = "#F9EBEA",
                    "high loss_artefact?" = "#FADBD8",
                    "loss?" = "#F1948A",
                    "loss_start" = "#B88A85",
                    "loss" = "#E74C3C",
                    "high loss" = "#B03A2E", 
                    "post mortem!" = "grey",
                    "didnt_check" = "grey",
                    
                    "neutral" = "lightgrey", 
                    "LOH" = "#F0DA4D", #technically this should be CN-LOH!!!!!
                    "fragmented" = "white",
                    `NA` = "#E2E2E0"
)



# CNV colour list for CNV heatmap annotations 
colours_WT$heatmap_colour_list <-
  list(CNV_Chr1p = colours_WT$CNV,
       CNV_Chr1p_combined = colours_WT$CNV,
       CNV_Chr1q = colours_WT$CNV,
       CNV_Chr1q_partial = colours_WT$CNV,
       
       CNV_Chr2p = colours_WT$CNV,
       CNV_Chr2p_partial = colours_WT$CNV,
       CNV_Chr2q = colours_WT$CNV,
       CNV_Chr2q_partial = colours_WT$CNV,
       
       CNV_Chr3p = colours_WT$CNV,
       CNV_Chr3p_partial = colours_WT$CNV,
       CNV_Chr3q = colours_WT$CNV,
       CNV_Chr3q_partial = colours_WT$CNV,
       
       CNV_Chr4p = colours_WT$CNV,
       CNV_Chr4p_partial = colours_WT$CNV,
       CNV_Chr4q = colours_WT$CNV,
       CNV_Chr4q_partial = colours_WT$CNV,
       
       CNV_Chr5p = colours_WT$CNV,
       CNV_Chr5p_partial = colours_WT$CNV,
       CNV_Chr5q = colours_WT$CNV,
       CNV_Chr5q_partial = colours_WT$CNV,
       
       CNV_Chr6p = colours_WT$CNV,
       CNV_Chr6p_partial = colours_WT$CNV,
       CNV_Chr6q = colours_WT$CNV,
       CNV_Chr6q_partial = colours_WT$CNV,
       
       CNV_Chr7p = colours_WT$CNV,
       CNV_Chr7p_partial = colours_WT$CNV,
       CNV_Chr7q = colours_WT$CNV,
       CNV_Chr7q_partial = colours_WT$CNV,
       
       CNV_Chr8p = colours_WT$CNV,
       CNV_Chr8p_partial = colours_WT$CNV,
       CNV_Chr8q = colours_WT$CNV,
       CNV_Chr8q_partial = colours_WT$CNV,
       
       CNV_Chr9p = colours_WT$CNV,
       CNV_Chr9p_partial = colours_WT$CNV,
       CNV_Chr9q = colours_WT$CNV,
       CNV_Chr9q_partial = colours_WT$CNV,
       
       CNV_Chr10p = colours_WT$CNV,
       CNV_Chr10p_partial = colours_WT$CNV,
       CNV_Chr10q = colours_WT$CNV,
       CNV_Chr10q_partial = colours_WT$CNV,
       
       CNV_Chr11p = colours_WT$CNV,
       CNV_Chr11p_partial = colours_WT$CNV,
       CNV_Chr11q = colours_WT$CNV,
       CNV_Chr11q_partial = colours_WT$CNV,
       
       CNV_Chr12p = colours_WT$CNV,
       CNV_Chr12p_partial = colours_WT$CNV,
       CNV_Chr12q = colours_WT$CNV,
       CNV_Chr12q_partial = colours_WT$CNV,
       
       CNV_Chr13p = colours_WT$CNV,
       CNV_Chr13p_partial = colours_WT$CNV,
       CNV_Chr13q = colours_WT$CNV,
       CNV_Chr13q_partial = colours_WT$CNV,
       
       CNV_Chr14p = colours_WT$CNV,
       CNV_Chr14p_partial = colours_WT$CNV,
       CNV_Chr14q = colours_WT$CNV,
       CNV_Chr14q_partial = colours_WT$CNV,
       
       CNV_Chr15p = colours_WT$CNV,
       CNV_Chr15p_partial = colours_WT$CNV,
       CNV_Chr15q = colours_WT$CNV,
       CNV_Chr15q_partial = colours_WT$CNV,
       
       CNV_Chr16p = colours_WT$CNV,
       CNV_Chr16p_partial = colours_WT$CNV,
       CNV_Chr16q = colours_WT$CNV,
       CNV_Chr16q_partial = colours_WT$CNV,
       
       CNV_Chr17p = colours_WT$CNV,
       CNV_Chr17p_partial = colours_WT$CNV,
       CNV_Chr17q = colours_WT$CNV,
       CNV_Chr17q_partial = colours_WT$CNV,
       
       CNV_Chr18p = colours_WT$CNV,
       CNV_Chr18p_partial = colours_WT$CNV,
       CNV_Chr18q = colours_WT$CNV,
       CNV_Chr18q_partial = colours_WT$CNV,
       
       CNV_Chr19p = colours_WT$CNV,
       CNV_Chr19p_partial = colours_WT$CNV,
       CNV_Chr19q = colours_WT$CNV,
       CNV_Chr19q_partial = colours_WT$CNV,
       
       CNV_Chr20p = colours_WT$CNV,
       CNV_Chr20p_partial = colours_WT$CNV,
       CNV_Chr20q = colours_WT$CNV,
       CNV_Chr20q_partial = colours_WT$CNV,
       
       CNV_Chr21p = colours_WT$CNV,
       CNV_Chr21p_partial = colours_WT$CNV,
       CNV_Chr21q = colours_WT$CNV,
       CNV_Chr21q_partial = colours_WT$CNV,
       
       CNV_Chr22p = colours_WT$CNV,
       CNV_Chr22p_partial = colours_WT$CNV,
       CNV_Chr22q = colours_WT$CNV,
       CNV_Chr22q_partial = colours_WT$CNV
  )






# colour scales and palettes for  methylation--------------------------------------------------------------------------------
# colour scale for beta values
colours_WT$meth <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101)

colours_WT$meth_MSMLPA_hyper <- c("Hypermethylated" = "red", "Normal"= "grey")
colours_WT$meth_MSMLPA_hypo <- c("Hypomethylated" = "blue", "Normal"= "grey")


# colour palette for hyper/hypomethylation in IC1 and IC2
colours_WT$IC1_IC2 = c( "hyper" = "#A31621",
                        "pUPD" = "#29335C",
                        "intermediate" = "#FFD131",
                        "normal" = "#999999")

# colours for DNAm based clustering
# colours_WT$DNAm_clustering = c(`1.a` = "#E69F00",
#                                `1.b` = "#F0E442",
#                                `2.a` = "#56B4E9",
#                                `2.b` = "#8B65F6", 
#                                `2.c` = "#0333FF")

colours_WT$DNAm_clustering = c(`1.a` = "#8C1834",
                               `1.b` = "#D1244C",
                               `2.a` = "#87BBD8",
                               `2.b` = "#34658A",
                               `2.c` = "#2F4858")

colours_WT$DNAm_clustering2 = c(`1a` = "#8C1834",
                               `1b` = "#D1244C",
                               `2a` = "#87BBD8",
                               `2b` = "#34658A",
                               `2c` = "#2F4858")

# colours for DNAm based clustering using only 2 clusters
colours_WT$DNAm_clustering_2cluster =c(`1` = "#E69F00",
                                       `2` = "#0333FF")

# colours for DNAm based clustering from very first clustering based on 96 DNAm samples
colours_WT$old_clustering <- c( "cluster1" = "#E69F00",
                                "cluster1(otherDNAmsample)" = "#E69F00",
                                "cluster2" = "#56B4E9",
                                "cluster3" = "#009E73",
                                "cluster4" = "#F0E442",
                                "postmortem" = "#0072B2" # there shouldnt be any sample with postmortem anymore, i removed them all
)

# Define color palette for BWSp categories
colours_WT$bwsp_colors <- c(
  "No" = "#E7EFEB",
  "Mosaic" = "#A5CDBB",
  "Clinical only" = "#7DA192",
  "Yes" = "#517164"
)

# colour palettes for specific genes ------------------------------------------------------------------------------------------
# WT1 germline PV
colours_WT$WT1_germline_variant = c("Yes" = "#785EF0",
                                    "no_WGS_WES_data" = "white", #this shouldnt be anywhere the case, because germline WT information is not based on our WGS/WES data!
                                    "No" = "grey")

# WT1 somatic and germline PV
colours_WT$WT1_heatmap <- c("somatic" = "#785EF0", #the names of this vector are not ideal, but works for now
                            "Yes" = "#1D0A70",
                            "germline" = "#1D0A70", 
                            "No" = "grey",
                            "no_WGS_WES_data" = "white"
)

# WT1 somatic PV
colours_WT$WT1_somatic_variant = c("Yes" = "#FE6100",
                                   "no_WGS_WES_data" = "white",
                                   "No" = "grey")

# WT1 deletion 
colours_WT$WT1_deletion_CNV = c("Yes" = "darkblue",
                                "No" = "grey")

# CTNNB1 somatic PV
colours_WT$CTNNB1_somatic_variant = c("Yes" = "#FBC449",
                                      "no_WGS_WES_data" = "white",
                                      "No" = "grey",
                                      "somatic" = "#FBC449")

# AMER1 somatic PV
colours_WT$AMER1_somatic_variant = c("Yes" = "darkgreen",
                                     "no_WGS_WES_data" = "white",
                                     "No" = "grey")

# TP53 somatic PV
colours_WT$TP53_somatic_variant = c("Yes" = "darkred",
                                    "no_WGS_WES_data" = "white",
                                    "No" = "grey")

colours_WT$bubble_WT1_CTNNB1 <- c(`WT1 somatic PV` = "#785EF0", 
                                  `CTNNB1 somatic PV` = "#FBC449",
                                  `WT1 germline PV` = "#1D0A70")




# GSEA colour scale ---------------------------------------------------------------
colours_WT$ggplot_GSEA <- scale_fill_gradient2(low="blue", mid = "white", high = "red")

