# title: "CNV analysis"
# author: "Franziska Oberhammer"
# date: "2026-09-10"
# use: script to generate CNV plots based on EPIC data 

# load libraries ----------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library("minfi") 
library("conumee2")  # devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
library("minfiData") 
library('minfiDataEPIC')
library(circlize)
library(dplyr)
library(rlang)
library(qs)
library(biomaRt)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38) # remotes::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#remotes::install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest")



# configuration ----------------------------------------------------------------------------------------------------------------------------------
parameters <- list()

# CNV analysis folder
parameters$CNV_FOLDER <- "/hpc/pmc_kuiper/Wilms/METHYLATION/ANALYSIS/CNV"

## folder with the idat files of the samples that need to be analyzed ##
parameters$dataFolder <- "/hpc/pmc_kuiper/Wilms/METHYLATION/DATA/PMCPR001ACC/"

## folder with a set of healthy control samples to use for normalisation ##
## control samples were downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246337
parameters$controlFolder <- "/hpc/pmc_kuiper/Wilms/METHYLATION/DATA/BIN/controlSamples_GSE246337"

# sample sheet
parameters$sampleSheetPath <- "/hpc/pmc_kuiper/Wilms/METHYLATION/ANALYSIS/PROCESSED_DATA/sample_sheet_updated_20250923.rds"

# colour palette
source(file.path("/hpc/pmc_kuiper/Wilms/METHYLATION/CODE", "00_WT_colour_palettes.R"))


# Set up functions -----------------------------------------------------------------------------------------------------------------------------  

### prepare data and sample sheets --- 
listEPICsamples <- function(folder){
  files <- list.files(folder,recursive = F,full.names=F,pattern = ".idat")
  files <- sub("_Grn.idat","",files)
  files <- sub("_Red.idat","",files)
  return(unique(files))
}

# generate control sample sheet
generateSampleSheet <- function(dataFolder,sampleNames){
  sampleSheet <- as.data.frame(matrix(nrow=length(sampleNames),ncol=3))
  colnames(sampleSheet) <- c("Sample_Name","Sentrix_ID","Sentrix_Position")
  sampleSheet$Sample_Name <- sampleNames
  sampleSheet$Sentrix_ID <- sub(".*_(\\d+)_.*", "\\1", sampleNames)
  sampleSheet$Sentrix_Position <- sub(".*_(R\\d+C\\d+)", "\\1", sampleNames)
  return(sampleSheet)
}

### load and preprocessing EPIC data ----
loadAndPreprocessEPICsamplesToCNVdata <- function(dataFolder, sampleSheet, EPICVersion){
  RGset <- minfi::read.metharray.exp(file.path(dataFolder),  sampleSheet, force = TRUE) 
  
  if (EPICVersion == "EPICv1"){
    annotation(RGset) <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
    RGset_normalised <- minfi::normalize.illumina.control(rgSet = RGset)
    RGset_normalised_backgroundCorrected <- minfi::bgcorrect.illumina(RGset_normalised)
    MSet <- minfi::preprocessRaw(RGset_normalised_backgroundCorrected)
    
  } else if (EPICVersion == "EPICv2"){
    annotation(RGset) <- c(array = "IlluminaHumanMethylationEPIC", annotation = "20a1.hg38")  
    RGset_normalised <- minfi::normalize.illumina.control(rgSet = RGset)
    RGset_normalised_backgroundCorrected <- minfi::bgcorrect.illumina(RGset_normalised)
    annotation(RGset_normalised_backgroundCorrected) <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")  
    MSet <- minfi::preprocessRaw(RGset_normalised_backgroundCorrected)
    
    # some probes in EPIC array v2 have replicates. For these ones we just take the first CpG in the dataframe. 
    all_CpGs <- rownames(MSet@assays@data$Meth)  
    all_CpGs_shortenend <- gsub("_.*$", "", all_CpGs)
    all_CpGs_unique <- all_CpGs[!duplicated(all_CpGs_shortenend)]
    
    MSet <- MSet[all_CpGs_unique,]
    
    #rename probes so that they are named like the EPICv1 probes (without the suffix)
    rownames(MSet) <-  gsub("_.*$", "", rownames(MSet))
  }
  
  CNVdata <- conumee2::CNV.load(MSet)
  return(CNVdata = CNVdata)
}



# process CNV data 
segmentData <- function(patient, control, annotation){
  print("CNV fit -------------------------")
  x <- conumee2::CNV.fit(patient, control, annotation) 
  
  print("CNV bin -------------------------")
  x <- conumee2::CNV.bin(x) 
  
  print("CNV bin - remove NAs ------------")
  # Remove NA values for all samples, otherwise segmentation doesnt work
  x@bin$variance <- lapply(x@bin$variance, function(sample) sample[!is.na(sample)])
  x@bin$ratio <- lapply(x@bin$ratio, function(sample) sample[!is.na(sample)])
  x@anno@bins     <- x@anno@bins[names(x@anno@bins) %in% names(x@bin$variance[[1]])]
  x@bin$variance  <- lapply(x@bin$variance, function(sample) sample[names(sample) %in% names(x@anno@bins)])
  x@bin$ratio     <- lapply(x@bin$ratio, function(sample) sample[names(sample) %in% names(x@anno@bins)])
  
  # continue with CNV pipeline  
  print("CNV detail ----------------------")
  x <- conumee2::CNV.detail(x) 
  
  print("CNV segment ---------------------")
  x <- conumee2::CNV.segment(x,alpha=0.005)
  return(x)
}



# load patient data ----------------------------------------------------------------------------------------------------------------------------------------
## sample sheets --------
sample_sheet_general <- readRDS(parameters$sampleSheetPath)$SAMPLESHEET %>% as.data.frame()

sample_sheets <- list()

sample_sheets$sample_sheet <- sample_sheet_general %>% 
  mutate(Sample_Name = Internal_Patient_ID, 
         Basename = barcode) %>% 
  filter(!is.na(Basename)) %>% 
  filter(Internal_Patient_ID != "EB")

# split in EPIC v1 and v2 cohort
sample_sheets$sample_sheet_EPICV1 <- sample_sheets$sample_sheet %>% 
  filter(EPIC_Version =="EPIC_V1") 


sample_sheets$sample_sheet_EPICV2 <- sample_sheets$sample_sheet %>% 
  filter(EPIC_Version =="EPIC_V2")




## load & preprocess samples to CNV data  ------------------------------------------------------------------------------------------------------------------------
patient_data <- list()

### Epic v1 ----
# Uncomment if you want to rerun it, takes quite long though
# patient_data$patient_V1_data <- loadAndPreprocessEPICsamplesToCNVdata(dataFolder = parameters$dataFolder, sampleSheet = sample_sheets$sample_sheet_EPICV1, EPICVersion = "EPICv1")
# qsave(patient_data$patient_V1_data, file.path(parameters$CNV_FOLDER, "CNVdata_EPICV1_InternalPatientIDs_hg38.qs"))
# patient_data$patient_V1_data <- qread(file.path(parameters$CNV_FOLDER, "CNVdata_EPICV1_InternalPatientIDs.qs")) (with hg19 manifest)
patient_data$patient_V1_data <- qread(file.path(parameters$CNV_FOLDER, "CNVdata_EPICV1_InternalPatientIDs_hg38.qs")) #(with hg38 manifest)

### Epic v2 ----
#uncomment if you want to rerun it
#patient_data$patient_V2_data <- loadAndPreprocessEPICsamplesToCNVdata(dataFolder = parameters$dataFolder,sampleSheet = sample_sheets$sample_sheet_EPICV2, EPICVersion = "EPICv2")
#qsave(patient_data$patient_V2_data, file.path(parameters$CNV_FOLDER, "CNVdata_EPICV2_InternalPatientIDs.qs"))
patient_data$patient_V2_data <- qread(file.path(parameters$CNV_FOLDER, "CNVdata_EPICV2_InternalPatientIDs.qs"))


### Combine V1 & V2 ----
# combine CNV data from V1 and V2 samples into one object
probes_on_both_arrays <- intersect(patient_data$patient_V1_data@intensity %>% rownames(), patient_data$patient_V2_data@intensity %>% rownames())

patient_data$patient_V1_data_overlapping_probes <- patient_data$patient_V1_data
patient_data$patient_V1_data_overlapping_probes@intensity <- patient_data$patient_V1_data@intensity[probes_on_both_arrays,]

patient_data$patient_V2_data_overlapping_probes <- patient_data$patient_V2_data
patient_data$patient_V2_data_overlapping_probes@intensity <- patient_data$patient_V2_data@intensity[probes_on_both_arrays,]

# before combining the two dataframes, check whether order of cpgs is the same
stopifnot(rownames(patient_data$patient_V1_data_overlapping_probes@intensity) == rownames(patient_data$patient_V2_data_overlapping_probes@intensity))


#combine dataframes
patient_data$patient_combined <- patient_data$patient_V1_data_overlapping_probes # generate dummy CNV.data
patient_data$patient_combined@intensity <- cbind(patient_data$patient_V1_data_overlapping_probes@intensity, 
                                                 patient_data$patient_V2_data_overlapping_probes@intensity)

#order like in sample sheet
patient_data$patient_combined@intensity <- patient_data$patient_combined@intensity[,sample_sheets$sample_sheet$Internal_Patient_ID]


# prepare normal data ------------------------------------------------------------------
## load & preprocess controls ----
control_data <- list()

# list all normal samples that we have
control_data$control_sampleNames <- listEPICsamples(parameters$controlFolder)

# in the control data we have 500 samples. we dont use all of them, just a random selection
set.seed(42)
control_data$control_sampleNames_50RandomSamples<- control_data$control_sampleNames[sample(1:500, 50)]

# generate sample sheet
control_data$sample_sheet <- generateSampleSheet(dataFolder = parameters$controlFolder, 
                                                 sampleNames = control_data$control_sampleNames_50RandomSamples) %>% 
  mutate(Basename = Sample_Name)


#uncomment if you want to rerun the preprocessing.
# control_data$CNVdata <- loadAndPreprocessEPICsamplesToCNVdata(dataFolder = parameters$controlFolder, sampleSheet = control_data$sample_sheet, EPICVersion = "EPICv2")
#qsave(x = control_data$CNVdata, file = file.path(parameters$CNV_FOLDER,  "CNV_control_samples.qs"))
control_data$CNVdata <- qread(file = file.path(parameters$CNV_FOLDER,  "CNV_control_samples.qs"))

# subset control CNV data for probes that are on both arrays
control_data$CNVdata_overlappingProbes <- control_data$CNVdata
control_data$CNVdata_overlappingProbes@intensity <- control_data$CNVdata@intensity[probes_on_both_arrays,]


# data processing ------------------------------------------------------------------------------------
## segmentation of the genome --- 
# uncomment if you want to rerun segmentation, otherwise load precalculated data
# cnv_result_woDetail   <- segmentData(patient = patient_data$patient_combined, control = control_data$CNVdata)
# cnv_result_withDetail <- segmentData(patient = patient_data$patient_combined, control = control_data$CNVdata)

## save cnv results
# qsave(cnv_result_woDetail, file.path(parameters$CNV_FOLDER, "cnv_result_woDetail_allSamples_epicV1AndV2hg38.qs"))
# qsave(cnv_result_withDetail, file.path(parameters$CNV_FOLDER, "cnv_result_withDetail_allSamples_epicV1AndV2hg38.qs"))

# load cnv results
cnv_result_woDetail <- qread(file.path(parameters$CNV_FOLDER, "cnv_result_woDetail_allSamples_epicV1AndV2hg38.qs"))
cnv_result_withDetail <- qread(file.path(parameters$CNV_FOLDER, "cnv_result_withDetail_allSamples_epicV1AndV2hg38.qs"))

# make plots -----
## plots per sample ----
# uncomment if you want to remake them (aka resave to the directory)
# CNV.genomeplot(cnv_result_withDetail, output = "pdf", directory = file.path(parameters$CNV_FOLDER, "plots", "Version_4.0", "withDetails", "chr11"), chr = "chr11", height = 12)
# CNV.genomeplot(cnv_result_withDetail, output = "pdf", directory = file.path(parameters$CNV_FOLDER, "plots", "Version_4.0", "withDetails", "genomewide"), height = 12)
# CNV.genomeplot(cnv_result_woDetail, output = "pdf", directory = file.path(parameters$CNV_FOLDER, "plots", "Version_4.0", "woDetails", "genomewide"), height = 12)



## summary plots per cluster -----------------------------------------------------------------------------------------------------------------------------
# parameters <- list()
parameters$PROJECT_FOLDER               <- "/hpc/pmc_kuiper/Wilms/METHYLATION" 
parameters$ANALYSIS_FOLDER              <- file.path(parameters$PROJECT_FOLDER, "ANALYSIS") 
parameters$PROCESSED_DATA_FOLDER        <- file.path(parameters$ANALYSIS_FOLDER, "PROCESSED_DATA")
parameters$CURRENT_SAMPLESHEET          <- file.path(parameters$PROCESSED_DATA_FOLDER, 
                                                     "masterSampleSheet_after_figure1.rds")
sample_sheet                            <- readRDS(parameters$CURRENT_SAMPLESHEET) %>% as.data.frame()

# define functions----------
generate_CNV_summaryplot <- function(cluster_label, sample_sheet, grouping_column, cnv_result) {
  # Get sample IDs for the given cluster
  samples <- sample_sheet %>%
    filter({{ grouping_column }} == cluster_label) %>%
    pull(Internal_Patient_ID)
  
  nr_of_samples <- length(samples)
  
  # Plot
  conumee2::CNV.summaryplot(cnv_result[samples], main = paste0(cluster_label, " (n = ", nr_of_samples, ")"))
}


# Helper function to write PDF plots
save_CNV_plot <- function(cluster_label, grouping_column, filename, sample_sheet, cnv_result, height = 4, width = 10) {
  pdf(file = filename, height = height, width = width)
  generate_CNV_summaryplot(cluster_label, sample_sheet, {{ grouping_column }}, cnv_result)
  dev.off()
}

# Output folder
output_dir <- file.path(parameters$CNV_FOLDER, "plots", "manuscript", "summary")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define clusters to plot
clusters <- c("1a", "1b", "2a", "2b", "2c")

# Loop over clusters
for (cl in clusters) {
  filename <- file.path(output_dir, paste0("20260605_CNV_summary_", gsub("\\.", "", cl), ".pdf"))
  save_CNV_plot(cl, DNAm_5_cluster_letters, filename, sample_sheet, cnv_result_woDetail)
}


#' Generate CNV summary table for a sample cluster
#' This function extracts CNV segments for all samples in a specified cluster,
#' computes disjoint genomic segments, counts the number of gains, losses,
#' and balanced regions per segment, and returns a table with both counts
#' and percentages of samples exhibiting each alteration.
#' Unlike `generate_CNV_summaryplot()`, this version does not plot but
#' outputs the underlying data used for plotting.
#'
#' @param cluster_label Name of the cluster of interest (e.g. "ClusterA")
#' @param sample_sheet Data frame with sample metadata (must contain `Internal_Patient_ID`)
#' @param grouping_column Column in `sample_sheet` defining clusters
#' @param cnv_result CNV analysis object
#' @param threshold Numeric, CNV calling threshold (default: 0.1)
#'
#' @return A data frame of disjoint genomic segments with counts and percentages
#'         of samples showing gains, losses, or balanced CNVs
## Generate CNV summary table function ----
generate_CNV_summarytable <- function(cluster_label, sample_sheet, grouping_column, cnv_result, threshold = 0.1) {
  # --- Get sample IDs for the given cluster ---
  samples <- sample_sheet %>%
    dplyr::filter({{ grouping_column }} == cluster_label) %>%
    dplyr::pull(Internal_Patient_ID)
  
  nr_of_samples <- length(samples)
  
  # --- Export CNV thresholded calls ---
  y <- CNV.write(cnv_result[samples], what = "threshold", threshold = threshold)
  
  # --- Build genomic segments ---
  segments.i <- GenomicRanges::GRanges(seqnames = y$Chromosome,
                                       ranges   = IRanges::IRanges(y$Start_Position, y$End_Position))
  segments_chromosomes <- GenomicRanges::GRanges(seqnames = cnv_result@anno@genome$chr,
                                                 ranges   = IRanges::IRanges(start = 1, end = cnv_result@anno@genome$size))
  segments <- c(segments.i, segments_chromosomes)
  d_segments <- as.data.frame(GenomicRanges::disjoin(segments))
  
  # --- Count alterations per segment ---
  overview <- as.data.frame(matrix(nrow = 0, ncol = 4))
  for (i in 1:nrow(d_segments)) {
    x <- d_segments[i,]
    involved_segments <- y[y$Chromosome == x$seqnames &
                             y$Start_Position <= x$start &
                             y$End_Position   >= x$end, ]
    balanced <- sum(involved_segments$Alteration == "balanced")
    gain     <- sum(involved_segments$Alteration == "gain")
    loss     <- sum(involved_segments$Alteration == "loss")
    overview <- rbind(overview, c(as.character(x$seqnames), balanced, gain, loss))
  }
  
  # --- Format overview ---
  colnames(overview) <- c("disjoined_segment", "count_balanced", "count_gains", "count_losses")
  overview$count_balanced <- as.numeric(overview$count_balanced)
  overview$count_gains    <- as.numeric(overview$count_gains)
  overview$count_losses   <- as.numeric(overview$count_losses)
  
  # --- Add percentages to d_segments ---
  d_segments$gains    <- overview$count_gains    / nr_of_samples * 100
  d_segments$losses   <- overview$count_losses   / nr_of_samples * 100
  d_segments$balanced <- overview$count_balanced / nr_of_samples * 100
  
  # --- Return final table ---
  return(d_segments)
}

# get summary table per cluster --
CNV_summaryTableList <- list()
CNV_summaryTableList$CNV_summaryTable_1a <- generate_CNV_summarytable("1a", sample_sheet, DNAm_5_cluster_letters, cnv_result_woDetail)
CNV_summaryTableList$CNV_summaryTable_1b <- generate_CNV_summarytable("1b", sample_sheet, DNAm_5_cluster_letters, cnv_result_woDetail)
CNV_summaryTableList$CNV_summaryTable_2a <- generate_CNV_summarytable("2a", sample_sheet, DNAm_5_cluster_letters, cnv_result_woDetail)
CNV_summaryTableList$CNV_summaryTable_2b <- generate_CNV_summarytable("2b", sample_sheet, DNAm_5_cluster_letters, cnv_result_woDetail)
CNV_summaryTableList$CNV_summaryTable_2c <- generate_CNV_summarytable("2c", sample_sheet, DNAm_5_cluster_letters, cnv_result_woDetail)


# important events per cluster (at least 20% of samples have a aberations)
important_events_per_cluster <- lapply(CNV_summaryTableList, function(df) {
  df %>%
    filter(gains >= 20 | losses >= 20) %>%
    arrange(desc(pmax(gains, losses)))
})
important_events_per_cluster$CNV_summaryTable_1a
important_events_per_cluster$CNV_summaryTable_1b
important_events_per_cluster$CNV_summaryTable_2a
important_events_per_cluster$CNV_summaryTable_2b
important_events_per_cluster$CNV_summaryTable_2c



important_events_per_cluster_combined <- bind_rows(important_events_per_cluster, .id = 'id')
all_events_per_cluster_combined <- bind_rows(CNV_summaryTableList, .id = 'id')


# write.table(x = important_events_per_cluster_combined, file = file.path(parameters$CNV_FOLDER, "importantEventsPerCluster_patientRemoved.csv"), row.names = F, quote = F, col.names = T)
# write.table(x = all_events_per_cluster_combined, file = file.path(parameters$CNV_FOLDER, "allEventsPerCluster_patientRemoved.tsv"), row.names = F, quote = F, col.names = T, sep = "\t")
