# title: "CNV analysis of validation cohort"
# author: "Franziska Oberhammer"
# date: "2026-09-10"


# load libraries
library(tidyverse)
library(minfi) 
library(circlize)
library(conumee2)  # devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
library(minfiData) 
library(minfiDataEPIC)
library(circlize)
library(dplyr)
library(stringr)
library(rlang)
library(qs)
library(RColorBrewer)
library(biomaRt)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b5.hg38) # remotes::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#remotes::install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest")

# setup ---------------------------------------------------------------------------------------------

## script to generate CNV plots based on EPIC data ##
parameters <- list()

# CNV analysis folder
parameters$CNV_FOLDER <- "/hpc/pmc_kuiper/Wilms/METHYLATION/ANALYSIS/CNV/validationCohort"

## folder with the idat files of the samples that need to be analyzed 
parameters$dataFolderInternal <- "/hpc/pmc_kuiper/Wilms/METHYLATION/DATA/PMCPR001ACC"
parameters$dataFolderExternal <- "/hpc/pmc_kuiper/Wilms/EXTERNAL_DATA/2025_Treger_CancerDiscovery_SomaticFootprintPredisposition/DATA/METHYLATION_EPIC/GSE292896_idat"

## folder with a set of healthy control samples to use for normalisation ##
## control sampels were downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246337
parameters$controlFolder <- "/hpc/pmc_kuiper/Wilms/METHYLATION/DATA/BIN/controlSamples_GSE246337"

# sample sheet
parameters$sampleSheetPathInternal <- "/hpc/pmc_kuiper/Wilms/METHYLATION/DATA/SAMPLE_SHEETS/20260910_Wilms_master_samplelist.xlsx"
parameters$sampleSheetPathExternal <- "/hpc/pmc_kuiper/Wilms/EXTERNAL_DATA/2025_Treger_CancerDiscovery_SomaticFootprintPredisposition/DATA/METHYLATION_EPIC/GSE292896_series_matrix.txt"


source(file.path("/hpc/pmc_kuiper/Wilms/METHYLATION/CODE", "00_WT_colour_palettes.R"))


# Set up functions -------------------------------------------------------------------------------------

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
    
    # remove EPIC array v2 probe replicates
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

# Prepare annotation ---------------------------------------------------------------------
# create annotation objects
anno <- conumee2::CNV.create_anno(array_type = c("EPICv2"), 
                                  genome = "hg38")  
annov1 <- conumee2::CNV.create_anno(array_type = c("EPIC"))  

# Load External Patient information, make sample sheet  ----------------------------------------------------------------------------------------------------------------------------------------

# Read the file as text
lines_External <- readLines(parameters$sampleSheetPathExternal)

# Keep only the lines starting with !Sample_
sample_lines_External <- grep("^!Sample_", lines_External, value = TRUE)

# Split each line into fields by tab
sample_split_External <- strsplit(sample_lines_External, "\t")

# Extract the row names (first element, without the "!" prefix)
rownames_External <- sapply(sample_split_External, `[`, 1)
rownames_External <- sub("^!", "", rownames_External)

# Extract the sample values (skip first column, which is the label)
values_External <- lapply(sample_split_External, function(x) x[-1])

# Turn into a data.frame: rows = sample attributes, columns = samples
samplesheet_External <- as.data.frame(do.call(rbind, values_External), stringsAsFactors_External = FALSE)

# Transpose so samples are rows and attributes are columns
samplesheet_External <- as.data.frame(t(samplesheet_External), stringsAsFactors = FALSE)

samplesheet_External[,c(1, 2, 8, 28)] %>% head()

samplesheet_External$Sample_Name <- sub("DNA from ", "", gsub("\"", "", samplesheet_External$V1))
samplesheet_External$Basename <- sub("_Grn.idat.gz", "", sub(".*/", "", gsub("\"", "", samplesheet_External$V28)))



# get colnumber of samples we want to investigate (where we know WT1 status)
sampleNamesWithWT1Info_External <- c("PD40713f", "PD40713h", "PD40713i", "PD40713j", "PD40713k", "PD48686a", "PD48687a", 
                                     "PD48688a", "PD48690a", "PD48691a", "PD48692a", "PD48696a", "PD48696c", "PD48696d", 
                                     "PD48699a", "PD48699c", "PD48700a", "PD48700c", "PD48700e", "PD48701a", "PD48701c", 
                                     "PD48705a", "PD48706a", "PD48708a", "PD48709a", "PD48710a", "PD48712a", "PD48713a", 
                                     "PD48714a", "PD48715a", "PD48717a", "PD48718a", "PD48718c", "PD48719a", "PD48720a", 
                                     "PD48721a", "PD48722a", "PD48723a", "PD48724a", "PD48725a", "PD48726c", "PD48726d", 
                                     "PD48726f", "PD48726h", "PD49171a", "PD49173a", "PD49174a", "PD49176a", "PD49177a",
                                     "PD49183a", "PD49184a", "PD49184c", "PD49185a", "PD49186a", "PD49186c", "PD49187a", 
                                     "PD49189a", "PD49189c", "PD49189e", "PD49191a", "PD49191c", "PD49193a", "PD49193c", 
                                     "PD49195a", "PD49197a", "PD49197c", "PD49198a", "PD49199a", "PD49199c", "PD49200a", 
                                     "PD49203a", "PD49204a", "PD49209a", "PD49210a", "PD49211a", "PD49213a", "PD49214a", 
                                     "PD49215a", "PD49216a", "PD49216c", "PD49217a", "PD49217c", "PD49219a", "PD49223a", 
                                     "PD49348m", "PD49348n", "PD49348o", "PD49348p", "PD49348r", "PD49348s", "PD49348t", 
                                     "PD49348u", "PD50589a", "PD50589c", "PD50589d", "PD50589f", "PD50589k", "PD50589l", 
                                     "PD50589m", "PD50589n", "PD50590a", "PD50590p", "PD50591a", "PD50591c", "PD50591d", 
                                     "PD50593a", "PD50594a", "PD50594d", "PD50596a", "PD50596c", "PD50596d", "PD50596e", 
                                     "PD50599e", "PD50599f", "PD50600g", "PD50600h", "PD50600i", "PD50602a", "PD50602c", 
                                     "PD50602d", "PD50643d", "PD50643e", "PD50643f", "PD50643h", "PD50643j", "PD50643k", 
                                     "PD50643l", "PD50662a", "PD50663a", "PD50665a", "PD50667a", "PD50669a", "PD50675a", 
                                     "PD50678a", "PD50682a", "PD50683a", "PD50686a", "PD50695a", "PD50696a", "PD50697a", 
                                     "PD50698a", "PD50699a", "PD50707a", "PD50708a", "PD50709a", "PD50711a", "PD50713a", 
                                     "PD50715a", "PD50716a", "PD50717a", "PD50719a", "PD50720a", "PD50724a", "PD50726a", 
                                     "PD50726c", "PD50728a", "PD50730a", "PD50733a", "PD50734a", "PD51622a", "PD51622c", 
                                     "PD51622d", "PD52201a", "PD52201c", "PD52201d", "PD52201f", "PD52202a", "PD52202c", 
                                     "PD52205a", "PD52205c", "PD52205d", "PD52209a", "PD52209c", "PD52210a", "PD52211c",
                                     "PD52211e", "PD52215a", "PD52218d", "PD52218f", "PD52219a", "PD52219c", "PD52220a",
                                     "PD52220c", "PD52221a", "PD52221c", "PD52222a", "PD52227a", "PD52227c", "PD52227d", 
                                     "PD52228a", "PD52230d", "PD52230e", "PD52235a", "PD52235c", "PD52236a", "PD52236c", 
                                     "PD52239a", "PD52239c", "PD52239d", "PD52239e", "PD52241a", "PD52241d", "PD52241f", 
                                     "PD52241h", "PD52241j", "PD52241k", "PD52241l", "PD53619a", "PD53619c", "PD53619e", 
                                     "PD53640a", "PD53641a", "PD53641c", "PD53643a", "PD53643c", "PD53647e", "PD53647f", 
                                     "PD53647g", "PD53651a", "PD53651c", "PD53651d", "PD53651e", "PD53654a", "PD53654c", 
                                     "PD53658a", "PD53658c", "PD53964a", "PD53964c", "PD53968a", "PD53969a", "PD53970a", 
                                     "PD54838a", "PD54840a", "PD54846f", "PD54846g", "PD54846h", "PD54846i")

samplesheet_filtered_External <- samplesheet_External %>% 
  filter(Sample_Name %in% sampleNamesWithWT1Info_External)



# Internal validation cohort sample sheet ------------------------------------------------------------
samplesheet_Internal <- readxl::read_xlsx(parameters$sampleSheetPathInternal) %>%
  filter(Cohort == "Test") %>% 
  filter(Methylation_Status == "DNAm_usedinClassification_Status20260902")  

# Add columns so CNV processing works
samplesheet_Internal <- samplesheet_Internal %>% 
  mutate(Basename = Methylation_barcode, 
         Sample_Name = Internal_Patient_ID)


## load & preprocess samples to CNV data  ----------------------------------------------------------------------------------------------------------------------------------------------------------------
patient_data <- list()

### Epic v1 ----
# uncomment if you want to rerun it
# patient_data$patient_V1_data_External <- loadAndPreprocessEPICsamplesToCNVdata(dataFolder = parameters$dataFolderExternal, sampleSheet = samplesheet_filtered_External, EPICVersion = "EPICv1")
# qsave(patient_data$patient_V1_data_External, file.path(parameters$CNV_FOLDER, "CNVdata_EPICV1_ExternalValidationCohort_hg38.qs"))
patient_data$patient_V1_data_External <- qread(file.path(parameters$CNV_FOLDER, "CNVdata_EPICV1_ExternalValidationCohort_hg38.qs"))

# uncomment if you want to rerun it
# patient_data$patient_V2_data_Internal <- loadAndPreprocessEPICsamplesToCNVdata(dataFolder = parameters$dataFolderInternal , sampleSheet = samplesheet_Internal, EPICVersion = "EPICv2")
# names(patient_data$patient_V2_data_Internal) <- samplesheet_Internal$M_Label[match(names(patient_data$patient_V2_data_Internal), samplesheet_Internal$Methylation_barcode)]
# qsave(patient_data$patient_V2_data_Internal, file.path(parameters$CNV_FOLDER, "CNVdata_EPICV2_InternalValidationCohort_hg38.qs"))
patient_data$patient_V2_data_Internal <- qread(file.path(parameters$CNV_FOLDER, "CNVdata_EPICV2_InternalValidationCohort_hg38.qs"))


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
#qsave(x = control_data$CNVdata, file = file.path(/hpc/pmc_kuiper/Wilms/METHYLATION/ANALYSIS/CNV/,  "CNV_control_samples.qs"))
control_data$CNVdata <- qread(file = file.path("/hpc/pmc_kuiper/Wilms/METHYLATION/ANALYSIS/CNV/",  "CNV_control_samples.qs"))





# data processing ------------------------------------------------------------------------------------
## segmentation of the genome ---
RERUN_SEGMENTATION <- FALSE
if(RERUN_SEGMENTATION == TRUE){
  # Internal cohort 
  cnv_result_woDetail_Internal   <- segmentData(patient = patient_data$patient_V2_data_Internal, control = control_data$CNVdata, annotation = anno)
  
  
  # # External cohort
  # External cohort - intersect EPICv1 probes across annotations and datasets
  # there was no hg38 annotation available for EPICv1 at the time of writing this code. Therefore we had to work around that. 
  probes_anno <- annov1@probes@ranges %>% names()
  probes_control <- control_data$CNVdata@intensity %>% rownames()
  probes_external <- patient_data$patient_V1_data_External@intensity %>% rownames()

  probes_controlAndAnno <- intersect(probes_anno,probes_control)
  probes_controlAndAnnoAndExternal <- intersect(probes_controlAndAnno, probes_external)
  
  
  # Subset datasets to overlapping probes
  patient_data$patient_V1_data_External_overlapping_probes <- patient_data$patient_V1_data_External
  patient_data$patient_V1_data_External_overlapping_probes@intensity <- patient_data$patient_V1_data_External_overlapping_probes@intensity[probes_controlAndAnnoAndExternal,]
  dim(patient_data$patient_V1_data_External_overlapping_probes@intensity)

  control_data$CNVdata_subset <-  control_data$CNVdata
  control_data$CNVdata_subset@intensity <- control_data$CNVdata_subset@intensity[probes_controlAndAnnoAndExternal,]
  dim(control_data$CNVdata_subset@intensity)

  annov1Subset <- annov1
  annov1Subset@probes <- annov1Subset@probes[probes_controlAndAnnoAndExternal,]
  length(annov1Subset@probes@ranges@NAMES)

  # Run CNV segmentation
  cnv_result_woDetail_External   <- segmentData(patient = patient_data$patient_V1_data_External_overlapping_probes, control = control_data$CNVdata_subset, annotation = annov1Subset)
  
  # ## save cnv results
  qsave(cnv_result_woDetail_External, file.path(parameters$CNV_FOLDER, "cnv_result_woDetail_ExternalCohort.qs"))
  qsave(cnv_result_woDetail_Internal, file.path(parameters$CNV_FOLDER, "cnv_result_woDetail_InternalCohort.qs"))
  # 
  
} else {
  # load cnv results
  cnv_result_woDetail_External <- qread(file.path(parameters$CNV_FOLDER, "cnv_result_woDetail_ExternalCohort.qs"))
  cnv_result_woDetail_Internal <- qread(file.path(parameters$CNV_FOLDER, "cnv_result_woDetail_InternalCohort.qs"))
  
}


# make plots -----
## plots per sample ----
# uncomment if you want to remake them (aka resave to the directory)
# CNV.genomeplot(cnv_result_woDetail_Internal, output = "pdf", directory = file.path(parameters$CNV_FOLDER, "plots", "woDetails", "genomewide"), height = 12)
# CNV.genomeplot(cnv_result_woDetail_External, output = "pdf", directory = file.path(parameters$CNV_FOLDER, "plots", "woDetails", "genomewide"), height = 12)


