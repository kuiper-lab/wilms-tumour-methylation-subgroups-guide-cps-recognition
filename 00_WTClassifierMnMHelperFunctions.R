# Classifier helper functions from MnM (https://github.com/princessmaximacenter/MnM/tree/main) that are needed for WT manuscript
# 15.09.2026
# Most functions are exactly the same as in the repository
# The major changes were in createFeatureDF & runMinorityClassifier


checkFormatInputData <- function(sampleColumn,
                                 classColumn,
                                 higherClassColumn,
                                 domainColumn,
                                 metaDataRef,
                                 countDataRef,
                                 outputDir = "NA",
                                 saveModel) {
  
  if (sampleColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the sample IDs is not present within metaDataRef. Please check the sampleColumn.")
  } else if (classColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor subtype labels is not present within metaDataRef. Please check the classColumn")
  } else if (higherClassColumn %notin% base::colnames(metaDataRef)){
    base::stop("The column you specified for the tumor type labels is not present within metaDataRef. Please check the higherClassColumn")
  } else if (domainColumn %notin% base::colnames(metaDataRef)) {
    base::stop("The column you specified for the tumor domain labels is not present within metaDataRef. Please check the domainColumn")
  }
  
  # Make sure the metadata and count data are in the right format and same order
  if (base::nrow(metaDataRef) != base::ncol(countDataRef)) {
    base::stop("The number of samples do not match between the metadata and the count data. Please make sure you include all same samples in both objects.")
  } else if (base::all(metaDataRef[, sampleColumn] %notin% base::colnames(countDataRef))) {
    base::stop("Your input data is not as required. Please make sure your sample IDs are stored in the sampleColumn, and in the column names of the count data")
  }
  
  
  if (base::is.numeric(countDataRef) != T) {
    base::stop("Your input data is not as required. Please make sure your countDataRef object only contains numerical count data and is a matrix.")
    
  }
  
  # Include a statement to store the classColumn, higherClassColumn and domainColumn
  base::cat(base::paste0("The column used for tumor subtypes labels within the metadata, used for model training purposes, is: ",
                         classColumn, '\nThis column contains values such as: \n'))
  base::cat(paste0(base::unique(metaDataRef[,classColumn])[1:3]), "\n")
  
  base::cat(base::paste0("\nThe column used for tumor type labels within the metadata, is: ",
                         higherClassColumn,'\nThis column contains values such as: \t'))
  base::cat(paste0(base::unique(metaDataRef[,higherClassColumn])[1:3]), "\n\n")
  
  base::cat(base::paste0("\nThe column used for tumor domain labels within the metadata, is: ",
                         domainColumn, '\nThis column contains values such as: \t'))
  base::cat(base::unique(metaDataRef[,domainColumn])[1:3])
  base::cat(paste0("\n\nIf any of these are incorrect, specify a different 'classColumn' (subtype),",
                   "\n'higherClassColumn' (tumor type) or 'domainColumn' (domain) to function as labels.\n\n"))
  
  if (saveModel == T & !base::dir.exists(outputDir)) {
    checkDirectory <- base::tryCatch(base::dir.create(outputDir))
    if (checkDirectory == F) {
      base::stop(base::paste0("The directory you want the classification to be saved in cannot be created due to an error in the directory path.",
                              " Please check the spelling of your specified outputDir - it is probable the parent-directory does not exist."))
    }
  }
  
}


getAnovaResults <- function(allMaterial,
                            allGenes,
                            classColumn) {
  anovaValues <- base::data.frame(allGenes = allGenes, p_val = NA, F_val = NA)
  for (i in base::seq(1:base::length(allGenes))) {
    ANOVARes <- stats::aov(allMaterial[,allGenes[i]] ~ allMaterial[, classColumn])
    anovaValues[i,"F_val"] <- base::summary(ANOVARes)[[1]][["F value"]][1]
    anovaValues[i,"p_val"] <- base::summary(ANOVARes)[[1]][["Pr(>F)"]][1]
  }
  return(anovaValues)
}

selectAnovaGenes <- function(metaDataRef,
                             countDataRef,
                             nANOVAgenes,
                             classColumn
) {
  
  countDataRef <- base::t(countDataRef) %>% base::as.data.frame()
  allGenes <- base::colnames(countDataRef)
  
  classesWith2 <- base::table(metaDataRef[,classColumn])[base::table(metaDataRef[,classColumn]) == 2] %>%
    base::names(.)
  
  if (base::length(classesWith2) > 0 ) {
    countDataRef$class <- base::as.factor(metaDataRef[rownames(countDataRef),classColumn])
    countDataRef <- createExtraData(countDataRef, classesWith2)
    metaDataRef$class <- base::as.character(metaDataRef[,classColumn])
    metaDataRef <- createExtraMetaData(metaDataRef = metaDataRef,
                                       classesWith2 = classesWith2)
  }
  
  countDataRef$Sample <- base::rownames(countDataRef)
  
  metaDataRef$Sample <- base::rownames(metaDataRef)
  allMaterial <- dplyr::left_join(metaDataRef, countDataRef, by = "Sample")
  
  # Perform an ANOVA test on all RNA-transcripts
  results <- getAnovaResults(allMaterial = allMaterial,
                             allGenes = allGenes,
                             classColumn = classColumn)
  
  # Arrange the RNA-transcripts so that the RNA-transcripts with the highest F-scores are at the top
  filterResults <- results %>% dplyr::arrange(dplyr::desc(F_val))
  
  # Select the top n F-score RNA-transcripts
  interestingAnovaGenes <- utils::head(filterResults$allGenes,
                                       n = nANOVAgenes)
  return(interestingAnovaGenes)
}

obtainTrainData <- function(metaDataRef, classColumn, maxSamplesPerType = 50, nModels = 100) {
  
  samplesTrainDefList <- base::list()
  metaDataRef[, classColumn] <- base::as.factor(metaDataRef[, classColumn])
  
  typesInFold <- base::table(metaDataRef[, classColumn])
  typeProbs <- 1/base::sqrt(typesInFold)
  
  for ( j in c(1:nModels)){
    #set.seed(j)
    samplesTrainDef <- c()
    
    samplesTrain <- base::unique(base::sample(base::rownames(metaDataRef),
                                              prob=typeProbs[metaDataRef[ , classColumn]],
                                              replace = T))
    
    includedTypes <- base::unique(metaDataRef[samplesTrain, classColumn])
    
    # Look which tumor types are not present in your selected dataset
    missingType <- metaDataRef[!(metaDataRef[ , classColumn] %in% includedTypes), ]
    
    # If there are tumor types not present, we loop through them all to add
    # one entry for each tumor type to the training samples.
    if (base::length(rownames(missingType)) > 0) {
      missingTumors <- base::unique(missingType[,classColumn])
      
      for(i in base::seq(1:base::length(missingTumors))) {
        currentTumorType <- missingType[missingType[ , classColumn] == missingTumors[i],]
        samplesTrain <- c(samplesTrain, base::sample(rownames(currentTumorType), size = 1))
      }
    }
    
    #
    nSamplesPerType <- base::table(metaDataRef[samplesTrain,classColumn])
    
    
    for (n in c(1:base::length(nSamplesPerType))){
      curSamples <- samplesTrain[metaDataRef[samplesTrain,classColumn] ==
                                   base::names(nSamplesPerType)[n]]
      
      # Entries that are present more often than the specified maxSamplesPerType are downsampled to maxSamplesPerType
      if (nSamplesPerType[n] > maxSamplesPerType){
        
        samplesTrainDef <- c(samplesTrainDef,base::sample(curSamples,maxSamplesPerType,replace = F))
        
      }else{
        samplesTrainDef <- c(samplesTrainDef,curSamples)
      }
    }
    samplesTrainDefList[[j]] <- samplesTrainDef
  }
  return(samplesTrainDefList)
}

reduceFeatures <- function(dataTrain,
                           samplesTrainDefList,
                           ntree = 500,
                           nModels =10,
                           nFeatures = 300,
                           nANOVAgenes) {
  
  modelList <- list()
  nModels <- base::min(nModels, 100)
  for (i in base::seq(1:nModels)) {
    samplesTrainDef <- samplesTrainDefList[[i]]
    
    train.data <- dataTrain[base::rownames(dataTrain) %in% samplesTrainDef,]
    
    train.category <- base::as.character(train.data$class)
    
    train.data %<>% dplyr::select(-c("class"))
    
    classesVal <- base::table(train.category)
    probabilityClasses <- 1/classesVal
    
    classwt <- base::as.numeric(probabilityClasses)
    
    model <- randomForest::randomForest(x = train.data, y = as.factor(train.category),
                                        importance = T, ntree = ntree,
                                        proximity = F, classwt = classwt)
    modelList[[i]] <- model
  }
  
  accuracyValuesDF <- createFeatureDF(modelList=modelList,
                                      whichAccuracyMeasure = "MeanDecreaseAccuracy",
                                      nANOVAgenes=nANOVAgenes)
  
  meanAccuracyValuesDF <- base::apply(accuracyValuesDF, 1, base::mean)
  
  topFeatures <- meanAccuracyValuesDF %>%
    base::sort(decreasing = T) %>%
    utils::head(n = nFeatures)
  
  topFeaturesNamesAccuracy <- base::names(topFeatures)
  
  return(topFeaturesNamesAccuracy)
}

createFeatureDF <- function(modelList,
                            whichAccuracyMeasure,
                            nANOVAgenes) {
  ##Create an empty dataframe
  empty_df <- base::data.frame()
  
  #Recursively extract column values from ModelList
  for (i in base::seq_along(modelList)) {
    
    columnValues <- base::as.data.frame(modelList[[i]][["importance"]]) %>%
      dplyr::select(tidyselect::all_of(whichAccuracyMeasure)) %>%
      dplyr::arrange(dplyr::desc(.data[[whichAccuracyMeasure]])) %>% 
      dplyr::slice(1:nANOVAgenes) %>%
      tibble::rownames_to_column(var = "Gene") %>%
      dplyr::mutate(model = i) %>%
      dplyr::rename(value = whichAccuracyMeasure) 
    
    #Append to Dataframe
    empty_df <- base::rbind(empty_df, columnValues)
  }
  
  #Create create a wide DF for comparing features
  accuracyValuesDF <- empty_df %>% tidyr::pivot_wider(names_from = "model")
  
  #Store columnName as rownames
  accuracyValuesDF <- accuracyValuesDF %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var="Gene")
  
  return(accuracyValuesDF)
}

obtainModelsMinorityClassifier <- function(dataTrain,
                                           samplesTrainDefList,
                                           nModels = 100,
                                           ntree = 500
) {
  modelList <- base::list()
  
  for (i in base::seq(1:nModels)) {
    base::print(base::paste("Working on model", i))
    trainSamples <- samplesTrainDefList[[i]]
    
    train.data <- dataTrain[base::rownames(dataTrain) %in% trainSamples,]
    
    train.category <- base::as.character(train.data$class)
    train.data %<>% dplyr::select(-c("class"))
    
    # Determine class weights for the weighted RF
    classesVal <- base::table(train.category)
    probabilityClasses <- 1/classesVal
    classwt <- base::as.numeric(probabilityClasses)
    
    # Generate RF model on training subset
    model <- randomForest::randomForest(x = train.data,
                                        y = base::as.factor(train.category),
                                        importance = T,
                                        ntree = ntree,
                                        proximity = F,
                                        classwt = classwt)
    
    modelList[[i]] <- model
  }
  
  return(modelList)
}

predictTest <- function(modelList, testData) {
  for (i in base::seq(1:base::length(modelList))) {
    model <- modelList[[i]]
    prediction <- predict(model, newdata=testData)
    
    if (i == 1) {
      result <- base::data.frame(fold1 = prediction)
    } else {
      result[, base::paste0("fold", i)] <- prediction
    }
  }
  
  return(result)
}

convertResultToClassification <- function(result,
                                          metaDataRef,
                                          addOriginalCall,
                                          classColumn = NA) {
  
  if (base::nrow(result) == 1 || typeof(base::apply(result, 1, base::table)) == 'integer') {
    randomVector <- base::paste0("fake", 1:base::ncol(result)) %>%
      base::as.data.frame() %>%
      base::t() %>% base::as.data.frame()
    base::colnames(randomVector) <- base::colnames(result)
    result1 <- base::rbind(result, randomVector)
    probability <-  base::apply(result1, 1, base::table)
    probability <- probability[-base::length(probability)]
  } else {
    # Find out how often a certain tumor type prediction is made for a specific sample
    probability <- base::apply(result, 1, base::table)
  }
  # Locate the position of the highest probability
  positions <- base::lapply(probability, base::which.max)
  positions <- base::unlist(positions)
  
  bestFit <- base::data.frame(predict = base::rep(NA, times = base::length(result$fold1)))
  probabilityScores <- base::vector()
  
  # Extract the different calls being made for each sample
  mostAppearingNames <- base::lapply(probability, base::names)
  
  # Store the one with the highest probability score into the bestFit dataframe
  for (j in base::seq(1:base::length(mostAppearingNames))) {
    numberPositions <- base::as.numeric(positions[j])
    probabilityScores[j] <- probability[[j]][numberPositions]
    bestFit[j,] <- mostAppearingNames[[j]][numberPositions]
  }
  
  
  # Store the bestFit, the originalCall and the accompanying probability score within the final dataframe.
  
  if (addOriginalCall == T) {
    # Look at the original calls for each test sample
    originalCall <- metaDataRef[base::rownames(result),classColumn]
    
    classifications <- base::cbind(predict = bestFit,
                                   originalCall = originalCall,
                                   probability = probabilityScores)
  } else {
    classifications <- base::cbind(predict = bestFit,
                                   probability = probabilityScores)
  }
  
  
  if (base::nrow(result) == 1) {
    classifications <- classifications[1, , drop = F]
    
  }
  # Make sure that the classifications still have their accompanying biomaterial_id
  base::rownames(classifications) <- base::rownames(result)
  
  classificationList <- base::list(classifications = classifications,
                                   probabilityList = probability,
                                   metaDataRef = metaDataRef
  )
  
  return(classificationList)
  
}







checkFormatTestData <- function(countDataNew,
                                countDataRef,
                                outputDir,
                                saveModel) {
  
  if (is.null(countDataRef)) {
    base::stop("You probably are using an old version of the model that is no longer compatible with MnM.\n\nPlease download the latest version from https://zenodo.org/records/14167359. ")
  }
  # Check whether the genes are within the rows of countDataNew
  geneOverlap <- sum(rownames(countDataNew) %in% rownames(countDataRef))
  
  if(geneOverlap == 0) {
    base::stop("Your input data is not as required. Please make sure your genes are in the rownames of countDataNew and are in the format of HGNC gene names.")
  } else if(geneOverlap < 0.6 * nrow(countDataRef)) {
    base::cat("Please note that there is less than 60% overlap between the genes within the reference cohort and countDataNew.\n")
    
  }
  
  # Check whether counts are supplied within the RNA-seq counts
  if (base::is.numeric(countDataNew) != T) {
    base::stop("Your input data is not as required. Please make sure your countDataNew object only contains numerical count data and is a matrix.")
  }
  
  # Generate the directory, if not possible abort
  if (saveModel == T & !base::dir.exists(outputDir)) {
    checkDirectory <- base::tryCatch(base::dir.create(outputDir))
    if (checkDirectory == F) {
      base::stop("The directory you want the classification to be saved in cannot be created due to an error in the directory path. Please check the spelling of your specified outputDir.")
    }
  }
  
}

calculateMissingGenes <- function(countDataRef,
                                  countDataNew,
                                  neededGenes,
                                  whichK = 3) {
  
  
  
  countDataNew <- as.matrix(countDataNew)
  countDataNew[is.nan(countDataNew)] <- NA # This needed to be added, otherwise script failed
  countDataNew[is.infinite(countDataNew)] <- NA # This needed to be added, otherwise script failed
  
  
  missingGenes <- neededGenes[neededGenes %notin% base::rownames(countDataNew)]
  halfmissingGenes <- rownames(countDataNew)[!complete.cases(countDataNew)]
  
  
  missingDF <- matrix(NA, nrow = length(missingGenes), ncol = base::ncol(countDataNew))
  base::rownames(missingDF) <- missingGenes
  base::colnames(missingDF) <- base::colnames(countDataNew)
  newDataFrame <- base::rbind(countDataNew,
                              missingDF)
  
  countDataNewSubset <- newDataFrame[neededGenes, , drop = F]
  countDataRefSubset <- countDataRef[neededGenes, , drop = F]
  
  hush=function(code){
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }
  
  # Subdivide dataset in chunks so that there are always more reference samples than samples to impute
  
  howManySplits <- base::ceiling(base::ncol(newDataFrame) / 500)
  
  
  splits <- base::rep(1:howManySplits, each = 500)
  splits <- splits[1:base::ncol(newDataFrame)]
  if (base::requireNamespace("impute") == F & base::requireNamespace("BiocManager") == F) {
    base::install.packages("BiocManager")
  } else if (base::requireNamespace("impute") == F) {
    BiocManager::install("impute")
  }
  
  for (i in unique(splits)) {
    combiSubset <- base::cbind(sesame::BetaValueToMValue(countDataRefSubset), countDataNewSubset[,splits == i, drop = F]) 
    #imputedDataComplete <- hush(impute::impute.knn(base::as.matrix(combiSubset), k = whichK)) # removed this, it opens a connections which leads to the results not being printed to console anymore
    combiSubset <- base::as.matrix(combiSubset)
    
    imputedDataComplete <- impute::impute.knn(combiSubset, k = whichK)
    
    imputedData <- imputedDataComplete$data[c(halfmissingGenes, missingGenes),base::colnames(countDataNewSubset[,splits == i, drop = F]), drop = F]
    
    if (i == 1) {
      imputedDataTotal <- imputedData
    } else {
      imputedDataTotal <- base::cbind(imputedDataTotal, imputedData)
    }
  }
  
  countDataNewTotal <-  base::rbind(countDataNew[setdiff(rownames(countDataNew), rownames(imputedDataTotal)), ], imputedDataTotal)
  
  return(countDataNewTotal)
}





runMinorityClassifier <- function(createdModelsMinorityList,
                                  MethDataNew,
                                  outputDir,
                                  sample_sheet,
                                  saveModel = TRUE,
                                  whichKimputation = 3) {
  
  
  source("/hpc/pmc_kuiper/Wilms/METHYLATION/CODE/00_WTClassifierMnMHelperFunctions.R")
  
  `%notin%` <<- Negate(`%in%`)
  
  # Prepare output dir
  outputDirPredictions <- file.path(outputDir)
  
  # Track missing values per sample (optional plot)
  missingValuesPerSample <- colSums(is.na(MethDataNew))
  hist(missingValuesPerSample)
  
  # Convert beta to M-values # important for methylation data
  countDataNew <- sesame::BetaValueToMValue(MethDataNew)
  
  # Missing feature imputation
  neededGenes <- createdModelsMinorityList$reducedFeatures
  missingGenes <- neededGenes[neededGenes %notin% rownames(countDataNew)]

  if (length(missingGenes) > 0 | sum(is.na(countDataNew)) > 0) {
    cat(paste0(
      "There are ", length(missingGenes),
      " genes missing from the dataset for classification.\nImputing their values.\n"
    ))

    countDataRef <- createdModelsMinorityList$countDataRef

    countDataNew <- calculateMissingGenes(
      countDataNew = countDataNew,
      neededGenes = neededGenes,
      countDataRef = countDataRef,
      whichK = whichKimputation
    )
  }
  
  # Transpose to sample rows × features
  dataNew <- t(countDataNew) %>% as.data.frame()
  dataNew <- dataNew[, createdModelsMinorityList$reducedFeatures, drop = FALSE]
  
  # Predict per-model
  result <- predictTest(
    modelList = createdModelsMinorityList$modelList,
    testData = dataNew
  )
  
  cat("\nFinished with classifying results.")
  
  # Majority vote
  classificationList <- convertResultToClassification(
    result = result,
    metaDataRef = createdModelsMinorityList$metaDataRef,
    addOriginalCall = FALSE
  )
  
  # Add metadata
  classificationList$metaDataRun <- sample_sheet
  
  # Save RDS
  if (saveModel) {
    filename <- file.path(outputDir, "minorityClassifierResult.rds")
    saveRDS(classificationList, file = filename)
    cat(paste0("\nSaved classification results to ", filename, "\n"))
  }
  
  return(classificationList)
}

