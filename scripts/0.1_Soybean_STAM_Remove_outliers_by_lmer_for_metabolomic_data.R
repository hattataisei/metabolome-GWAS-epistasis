##########################################################################################
######  Title: 0.1_Soybean_STAM_Remove_outliers_by_lmer_for_metabolomic_data        ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/09 (Created), 2025/01/07 (Last Updated)                       ######
##########################################################################################





###### 1. Settings ######
##### 1.0. Reset workspace ######
# rm(list=ls())



##### 1.1. Setting working directory to the "projectName" directory #####
# cropName <- "soybean"
# project <- "STAM"
# os <- osVersion
#
# isRproject <- function(path = getwd()) {
#   files <- list.files(path)
#
#   if (length(grep(".Rproj", files)) >= 1) {
#     out <- TRUE
#   } else {
#     out <-  FALSE
#   }
#   return(out)
# }

scriptID <- "0.1"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
# nPC <- 6

dirMidSTAMViewDataRemoveOutlier <- paste0(dirMidSTAMBase, scriptID,
                                                 "_View_data_and_remove_outlier/")
dir.create(dirMidSTAMViewDataRemoveOutlier)



thresSd <- 5


# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
install.packages("openxlsx")

require(openxlsx)
require(data.table)
require(lme4)
# require(heritability)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### Create list of 188 metabolites for year of 2017 and 2018 and so on ######
metabListAll <- read.xlsx("raw_data/extra/CREST_WT_265compounds_for_KEGG.xlsx")
See(metabListAll)
metabListAll <- metabListAll[, 1:3]
rownames(metabListAll) <- metabListAll[, 2]

metab188Metab <- read.csv("raw_data/phenotype/2017_Tottori_May_Metabolome.csv")
See(metab188Metab, coln = 15)
metabNames188Metab <- colnames(metab188Metab[, 10:ncol(metab188Metab)])
See(metabNames188Metab)

metabList188Metab <- metabListAll[metabNames188Metab, ]
See(metabList188Metab)

fileName <- paste0("raw_data/extra/0.1_List_of_188_metablites_not_ordered_in_wheter_flavonoid_or_not.csv")
write.csv(x = metabList188Metab, file = fileName, row.names = FALSE)

##### A. Check distributions of all matabolites #####
#### A.1 Without removing outliers ####
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/"))


metab2017Raw <- read.csv("raw_data/phenotype/2017_Tottori_May_Metabolome.csv")
See(metab2017Raw, rown = 6, coln = 12)

metabStart <- 10
metabEnd <- ncol(metab2017Raw)
metabNames <- colnames(metab2017Raw)[metabStart:metabEnd]

nMetab <- metabEnd - metabStart + 1

metab2017Control <- metab2017Raw[metab2017Raw$block == "C", ]
See(metab2017Control, coln = 10)
metab2017ControlOnlyMetab <- metab2017Control[, metabStart:metabEnd]
See(metab2017ControlOnlyMetab)

metab2017Drought <- metab2017Raw[metab2017Raw$block == "D", ]
metab2017DroughtOnlyMetab <- metab2017Drought[, metabStart:metabEnd]


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])


### Control
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/"))
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid/"))
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/", scriptID, "_Non_Flavonoid/"))
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/"))


metabNo <- 1
for(metabNo in 1:nMetab){
  metabNow <- metab2017ControlOnlyMetab[, metabNo]
  metabName <- metabNames[metabNo]

  if (metabName %in% metabNamesFlavonoid){
    pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
    hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
    # hist(metabNow, freq = FALSE, probability = TRUE)
    dev.off()

    } else {
      pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
    hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
    # hist(metabNow, freq = FALSE, probability = TRUE)
    dev.off()

  }
}


for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
  metabNow <- metab2017ControlOnlyMetab[, metabName]

  pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
  hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
  dev.off()

}



### Drought
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/"))
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid/"))
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/", scriptID, "_Non_Flavonoid/"))
dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/"))


metabNo <- 1
for(metabNo in 1:nMetab){
  metabNow <- metab2017DroughtOnlyMetab[, metabNo]
  metabName <- metabNames[metabNo]

  if (metabName %in% metabNamesFlavonoid){
    pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
    hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
    # hist(metabNow, freq = FALSE, probability = TRUE)
    dev.off()

  } else {
    pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
    hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
    # hist(metabNow, freq = FALSE, probability = TRUE)
    dev.off()

  }
}


for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
  metabNow <- metab2017DroughtOnlyMetab[, metabName]

  pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_With_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
  hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
  dev.off()

}



#### 2017
###### 2. Estimate broad-sense heritability for Metabolomic data in 2017 ######
##### 2.1. Read Metabolomic data in 2017 into R #####
metab2017Raw <- read.csv("raw_data/phenotype/2017_Tottori_May_Metabolome.csv")
See(metab2017Raw, rown = 6, coln = 12)


metabStart <- 10
metabEnd <- ncol(metab2017Raw)
metabNames <- colnames(metab2017Raw)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2017Raw$variety)
blockNames <- unique(metab2017Raw$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2017Raw$ind))
nBlock <- length(blockNames)



##### 2.3. Perform lmer for all Metabolomic data in 2017 #####
#### 2.3.1 Perform lmer for all Metabolomic data in 2017 ####
lmerVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)
lmerResidMat2017 <- matrix(NA, nrow = nrow(metab2017Raw), ncol = nMetab)

rownames(lmerVarHerits2017) <- metabNames
colnames(lmerVarHerits2017) <- varNames
rownames(lmerRanefMat2017) <- varietyNames
colnames(lmerRanefMat2017) <- metabNames
rownames(lmerResidMat2017) <- rownames(metab2017Raw)
colnames(lmerResidMat2017) <- metabNames

for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2017Raw[, metabNo]

  lmerRes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                      data = metab2017Raw)

  lmerVars2017 <- data.frame(VarCorr(lmerRes2017))$vcov
  lmerVu2017 <- lmerVars2017[1]
  lmerVe2017 <- lmerVars2017[2]

  lmerHeritInd2017 <- lmerVu2017 / (lmerVu2017 + lmerVe2017)
  lmerHeritLine2017 <- lmerVu2017 / (lmerVu2017 + lmerVe2017 / (nRep * nBlock))

  lmerRanef2017 <- ranef(lmerRes2017)$variety[varietyNames, ]
  lmerResid2017 <- residuals(lmerRes2017)


  lmerVarHerits2017[metabNo - metabStart + 1, ] <-
    c(lmerVars2017, lmerHeritInd2017, lmerHeritLine2017)
  lmerRanefMat2017[, metabNo - metabStart + 1] <- lmerRanef2017
  lmerResidMat2017[, metabNo - metabStart + 1] <- lmerResid2017
}

upperThres2017 <- apply(lmerResidMat2017, 2, mean) + thresSd * apply(lmerResidMat2017, 2, sd)
lowerThres2017 <- apply(lmerResidMat2017, 2, mean) - thresSd * apply(lmerResidMat2017, 2, sd)


outlier2017 <- sapply(1:nMetab, function(metabNo) {
  return((lmerResidMat2017[, metabNo] > upperThres2017[metabNo]) |
           (lmerResidMat2017[, metabNo] < lowerThres2017[metabNo]))
})

table(outlier2017)


#### 2.3.2 Perform lmer for all Metabolomic data for control in 2017 ####
metabControl2017 <- metab2017Raw[metab2017Raw$block == "C", ]

lmerVarHeritsControl2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatControl2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)
lmerResidMatControl2017 <- matrix(NA, nrow = nrow(metabControl2017), ncol = nMetab)

rownames(lmerVarHeritsControl2017) <- metabNames
colnames(lmerVarHeritsControl2017) <- varNames
rownames(lmerRanefMatControl2017) <- varietyNames
colnames(lmerRanefMatControl2017) <- metabNames
rownames(lmerResidMatControl2017) <- rownames(metabControl2017)
colnames(lmerResidMatControl2017) <- metabNames


for (metabNo in metabStart:metabEnd){
  metabControlNow <- metabControl2017[, metabNo]

  lmerControlRes2017 <- lmer(formula = metabControlNow ? (1 | variety),
                             data = metabControl2017)

  lmerControlVars2017 <- data.frame(VarCorr(lmerControlRes2017))$vcov
  lmerControlVu2017 <- lmerControlVars2017[1]
  lmerControlVe2017 <- lmerControlVars2017[2]

  lmerControlHeritInd2017 <- lmerControlVu2017 / (lmerControlVu2017 + lmerControlVe2017)
  lmerControlHeritLine2017 <- lmerControlVu2017 / (lmerControlVu2017 + lmerControlVe2017 / nRep)

  lmerControlRanef2017 <- ranef(lmerControlRes2017)$variety[varietyNames, ]
  lmerResidControl2017 <- residuals(lmerControlRes2017)

  lmerVarHeritsControl2017[metabNo - metabStart + 1, ] <-
    c(lmerControlVars2017, lmerControlHeritInd2017, lmerControlHeritLine2017)
  lmerRanefMatControl2017[, metabNo - metabStart + 1] <- lmerControlRanef2017
  lmerResidMatControl2017[, metabNo - metabStart + 1] <- lmerResidControl2017

}


upperThresControl2017 <- apply(lmerResidMatControl2017, 2, mean) + thresSd * apply(lmerResidMatControl2017, 2, sd)
lowerThresControl2017 <- apply(lmerResidMatControl2017, 2, mean) - thresSd * apply(lmerResidMatControl2017, 2, sd)


outlierControl2017 <- sapply(1:nMetab, function(metabNo) {
  return((lmerResidMatControl2017[, metabNo] > upperThresControl2017[metabNo]) |
           (lmerResidMatControl2017[, metabNo] < lowerThresControl2017[metabNo]))
})

table(outlierControl2017)

#### 2.3.3  Perform lmer for all Metabolomic data for drought in 2017 ####
metabDrought2017 <- metab2017Raw[metab2017Raw$block == "D", ]

lmerVarHeritsDrought2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatDrought2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)
lmerResidMatDrought2017 <- matrix(NA, nrow = nrow(metabDrought2017), ncol = nMetab)

rownames(lmerVarHeritsDrought2017) <- metabNames
colnames(lmerVarHeritsDrought2017) <- varNames
rownames(lmerRanefMatDrought2017) <- varietyNames
colnames(lmerRanefMatDrought2017) <- metabNames
rownames(lmerResidMatDrought2017) <- rownames(metabDrought2017)
colnames(lmerResidMatDrought2017) <- metabNames


for (metabNo in metabStart:metabEnd){
  metabDroughtNow <- metabDrought2017[, metabNo]

  lmerDroughtRes2017 <- lmer(formula = metabDroughtNow ? (1 | variety),
                             data = metabDrought2017)

  lmerDroughtVars2017 <- data.frame(VarCorr(lmerDroughtRes2017))$vcov
  lmerDroughtVu2017 <- lmerDroughtVars2017[1]
  lmerDroughtVe2017 <- lmerDroughtVars2017[2]

  lmerDroughtHeritInd2017 <- lmerDroughtVu2017 / (lmerDroughtVu2017 + lmerDroughtVe2017)
  lmerDroughtHeritLine2017 <- lmerDroughtVu2017 / (lmerDroughtVu2017 + lmerDroughtVe2017 / nRep)

  lmerDroughtRanef2017 <- ranef(lmerDroughtRes2017)$variety[varietyNames, ]
  lmerResidDrought2017 <- residuals(lmerDroughtRes2017)

  lmerVarHeritsDrought2017[metabNo - metabStart + 1, ] <-
    c(lmerDroughtVars2017, lmerDroughtHeritInd2017, lmerDroughtHeritLine2017)
  lmerRanefMatDrought2017[, metabNo - metabStart + 1] <- lmerDroughtRanef2017
  lmerResidMatDrought2017[, metabNo - metabStart + 1] <- lmerResidDrought2017

}


upperThresDrought2017 <- apply(lmerResidMatDrought2017, 2, mean) + thresSd * apply(lmerResidMatDrought2017, 2, sd)
lowerThresDrought2017 <- apply(lmerResidMatDrought2017, 2, mean) - thresSd * apply(lmerResidMatDrought2017, 2, sd)


outlierDrought2017 <- sapply(1:nMetab, function(metabNo) {
  return((lmerResidMatDrought2017[, metabNo] > upperThresDrought2017[metabNo]) |
           (lmerResidMatDrought2017[, metabNo] < lowerThresDrought2017[metabNo]))
})

table(outlierDrought2017)



#### 2.3.4 Gather outlier information for 2017 data ####
colnames(outlier2017) <- colnames(outlierControl2017) <-
  colnames(outlierDrought2017) <- metabNames

outlierCD2017 <- rbind(outlierControl2017,
                       outlierDrought2017)

outlierAll2017 <- outlier2017 | outlierCD2017

table(apply(outlierAll2017, 2, sum))

metabOnly2017 <- metab2017Raw[, metabStart:metabEnd]
metabOnly2017[outlierAll2017] <- NA


metab2017NoOL <- metab2017Raw
metab2017NoOL[, metabStart:metabEnd] <- metabOnly2017


fileNameNoOL2017 <- paste0("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_thres=",
                           thresSd, ".csv")
write.csv(x = metab2017NoOL, file = fileNameNoOL2017, quote = FALSE, row.names = FALSE)






# ##### B. Check distributions of all matabolites #####
# #### B.1. No outliers ####
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/"))
#
#
# metab2017RawNoOutlier <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
# See(metab2017RawNoOutlier, rown = 6, coln = 12)
#
# metabStart <- 10
# metabEnd <- ncol(metab2017RawNoOutlier)
# metabNames <- colnames(metab2017Raw)[metabStart:metabEnd]
#
# table(is.na(metab2017RawNoOutlier[, metabStart:metabEnd]))
#
# nMetab <- metabEnd - metabStart + 1
#
# metab2017Control <- metab2017RawNoOutlier[metab2017RawNoOutlier$block == "C", ]
# See(metab2017Control, coln = 10)
# metab2017ControlOnlyMetab <- metab2017Control[, metabStart:metabEnd]
# See(metab2017ControlOnlyMetab)
#
# metab2017Drought <- metab2017RawNoOutlier[metab2017RawNoOutlier$block == "D", ]
# metab2017DroughtOnlyMetab <- metab2017Drought[, metabStart:metabEnd]
#
#
# metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
# metabNamesFlavonoid <- metabFlavonoid[, "Name"]
# metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
# See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
# metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])
#
#
# ### Control
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/"))
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017ControlOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   }
# }
#
#
# for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
#   metabNow <- metab2017ControlOnlyMetab[, metabName]
#
#   pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
#
#
#
# ### Drought
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/"))
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017DroughtOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   }
# }
#
#
# for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
#   metabNow <- metab2017DroughtOnlyMetab[, metabName]
#
#   pdf(paste0(dirMidSTAMViewDataRemoveOutlier, scriptID, "_Withoout_outliers/", scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }




##### 3. Extracting data of metabolites only related to Flavonoid pathway in 2017 into R #####
#### 3.1. Read Metabolome data without Outlier in 2017 into R ####
metab2017Raw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
See(metab2017Raw, rown = 6, coln = 12)

OnlyMetab2017 <- metab2017Raw[, 10:ncol(metab2017Raw)]
See(OnlyMetab2017)


#### 3.2. Extracting data of metabolites only related to Flavonoid pathway in 2017 ####
CRESTMetabAnnotation <- read.xlsx("raw_data/extra/CREST_WT_265compounds_for_KEGG.xlsx")
metabOnlyRelatedToFlavonoidPathway <- CRESTMetabAnnotation[CRESTMetabAnnotation$XX == "a", ]
head(metabOnlyRelatedToFlavonoidPathway)
metabSampleIDRelatedToFlavonoidPathway <- metabOnlyRelatedToFlavonoidPathway[, 2]
metabAnnotationRelatedToFlavonoidPathway <- metabOnlyRelatedToFlavonoidPathway[, 3]

names(metabSampleIDRelatedToFlavonoidPathway) <- metabAnnotationRelatedToFlavonoidPathway


CommonMetabNames <- colnames(OnlyMetab2017[, colnames(OnlyMetab2017)%in% metabSampleIDRelatedToFlavonoidPathway])
table(colnames(OnlyMetab2017)%in% metabSampleIDRelatedToFlavonoidPathway)
FlavonoidMetab <- OnlyMetab2017[, colnames(OnlyMetab2017)%in% metabSampleIDRelatedToFlavonoidPathway]
See(FlavonoidMetab)

FlavonoidMetab <- cbind(metab2017Raw[, 1:9], FlavonoidMetab)
head(FlavonoidMetab)
See(FlavonoidMetab, coln = 15)


CommonAnnotation <- metabSampleIDRelatedToFlavonoidPathway[metabSampleIDRelatedToFlavonoidPathway %in% CommonMetabNames]
length(CommonAnnotation)
Annotation <- names(CommonAnnotation)
CommonMetabNamesAndAnnotationFlavo <- cbind(CommonMetabNames, Annotation)
colnames(CommonMetabNamesAndAnnotationFlavo)[1] <- "Name"

CommonMetabNamesAndAnnotationFlavo
See(CommonMetabNamesAndAnnotationFlavo)


fileNameMetabNoOutlierFlavonoid <- paste0("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway", ".csv")

fileNameFlavonoidMetabNamesAndAnnotation <- paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation", ".csv")


write.csv(x = FlavonoidMetab, file = fileNameMetabNoOutlierFlavonoid, quote = FALSE, row.names = FALSE)

write.csv(x = CommonMetabNamesAndAnnotationFlavo, file = fileNameFlavonoidMetabNamesAndAnnotation, row.names = FALSE)
# write.csv(x = CommonMetabNamesAndAnnotationFlavo, file = fileNameFlavonoidMetabNamesAndAnnotation, quote = FALSE, row.names = FALSE)






##### 2018
###### 4. Estimate broad-sense heritability for Metabolomic data in 2018 ######
##### 4.1. Read Metabolomic data in 2018 into R #####
metab2018Raw <- read.csv("raw_data/phenotype/2018_Tottori_May_Metabolome.csv")
See(metab2018Raw, rown = 6, coln = 12)


metabStart <- 10
metabEnd <- ncol(metab2018Raw)
metabNames <- colnames(metab2018Raw)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2018Raw$variety)
blockNames <- unique(metab2018Raw$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2018Raw$ind))
nBlock <- length(blockNames)



##### 4.3. Perform lmer for all Metabolomic data in 2018 #####
#### 4.3.1 Perform lmer for all Metabolomic data in 2018 ####
lmerVarHerits2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMat2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)
lmerResidMat2018 <- matrix(NA, nrow = nrow(metab2018Raw), ncol = nMetab)

rownames(lmerVarHerits2018) <- metabNames
colnames(lmerVarHerits2018) <- varNames
rownames(lmerRanefMat2018) <- varietyNames
colnames(lmerRanefMat2018) <- metabNames
rownames(lmerResidMat2018) <- rownames(metab2018Raw)
colnames(lmerResidMat2018) <- metabNames

for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2018Raw[, metabNo]

  lmerRes2018 <- lmer(formula = metabNow ? (1 | variety) + block,
                      data = metab2018Raw)

  lmerVars2018 <- data.frame(VarCorr(lmerRes2018))$vcov
  lmerVu2018 <- lmerVars2018[1]
  lmerVe2018 <- lmerVars2018[2]

  lmerHeritInd2018 <- lmerVu2018 / (lmerVu2018 + lmerVe2018)
  lmerHeritLine2018 <- lmerVu2018 / (lmerVu2018 + lmerVe2018 / (nRep * nBlock))

  lmerRanef2018 <- ranef(lmerRes2018)$variety[varietyNames, ]
  lmerResid2018 <- residuals(lmerRes2018)


  lmerVarHerits2018[metabNo - metabStart + 1, ] <-
    c(lmerVars2018, lmerHeritInd2018, lmerHeritLine2018)
  lmerRanefMat2018[, metabNo - metabStart + 1] <- lmerRanef2018
  lmerResidMat2018[, metabNo - metabStart + 1] <- lmerResid2018
}

upperThres2018 <- apply(lmerResidMat2018, 2, mean) + thresSd * apply(lmerResidMat2018, 2, sd)
lowerThres2018 <- apply(lmerResidMat2018, 2, mean) - thresSd * apply(lmerResidMat2018, 2, sd)


outlier2018 <- sapply(1:nMetab, function(metabNo) {
  return((lmerResidMat2018[, metabNo] > upperThres2018[metabNo]) |
           (lmerResidMat2018[, metabNo] < lowerThres2018[metabNo]))
})

table(outlier2018)


#### 4.3.2 Perform lmer for all Metabolomic data for control in 2018 ####
metabControl2018 <- metab2018Raw[metab2018Raw$block == "C", ]

lmerVarHeritsControl2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatControl2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)
lmerResidMatControl2018 <- matrix(NA, nrow = nrow(metabControl2018), ncol = nMetab)

rownames(lmerVarHeritsControl2018) <- metabNames
colnames(lmerVarHeritsControl2018) <- varNames
rownames(lmerRanefMatControl2018) <- varietyNames
colnames(lmerRanefMatControl2018) <- metabNames
rownames(lmerResidMatControl2018) <- rownames(metabControl2018)
colnames(lmerResidMatControl2018) <- metabNames


for (metabNo in metabStart:metabEnd){
  metabControlNow <- metabControl2018[, metabNo]

  lmerControlRes2018 <- lmer(formula = metabControlNow ? (1 | variety),
                             data = metabControl2018)

  lmerControlVars2018 <- data.frame(VarCorr(lmerControlRes2018))$vcov
  lmerControlVu2018 <- lmerControlVars2018[1]
  lmerControlVe2018 <- lmerControlVars2018[2]

  lmerControlHeritInd2018 <- lmerControlVu2018 / (lmerControlVu2018 + lmerControlVe2018)
  lmerControlHeritLine2018 <- lmerControlVu2018 / (lmerControlVu2018 + lmerControlVe2018 / nRep)

  lmerControlRanef2018 <- ranef(lmerControlRes2018)$variety[varietyNames, ]
  lmerResidControl2018 <- residuals(lmerControlRes2018)

  lmerVarHeritsControl2018[metabNo - metabStart + 1, ] <-
    c(lmerControlVars2018, lmerControlHeritInd2018, lmerControlHeritLine2018)
  lmerRanefMatControl2018[, metabNo - metabStart + 1] <- lmerControlRanef2018
  lmerResidMatControl2018[, metabNo - metabStart + 1] <- lmerResidControl2018

}


upperThresControl2018 <- apply(lmerResidMatControl2018, 2, mean) + thresSd * apply(lmerResidMatControl2018, 2, sd)
lowerThresControl2018 <- apply(lmerResidMatControl2018, 2, mean) - thresSd * apply(lmerResidMatControl2018, 2, sd)


outlierControl2018 <- sapply(1:nMetab, function(metabNo) {
  return((lmerResidMatControl2018[, metabNo] > upperThresControl2018[metabNo]) |
           (lmerResidMatControl2018[, metabNo] < lowerThresControl2018[metabNo]))
})

table(outlierControl2018)

#### 4.3.3  Perform lmer for all Metabolomic data for drought in 2018 ####
metabDrought2018 <- metab2018Raw[metab2018Raw$block == "D", ]

lmerVarHeritsDrought2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatDrought2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)
lmerResidMatDrought2018 <- matrix(NA, nrow = nrow(metabDrought2018), ncol = nMetab)

rownames(lmerVarHeritsDrought2018) <- metabNames
colnames(lmerVarHeritsDrought2018) <- varNames
rownames(lmerRanefMatDrought2018) <- varietyNames
colnames(lmerRanefMatDrought2018) <- metabNames
rownames(lmerResidMatDrought2018) <- rownames(metabDrought2018)
colnames(lmerResidMatDrought2018) <- metabNames


for (metabNo in metabStart:metabEnd){
  metabDroughtNow <- metabDrought2018[, metabNo]

  lmerDroughtRes2018 <- lmer(formula = metabDroughtNow ? (1 | variety),
                             data = metabDrought2018)

  lmerDroughtVars2018 <- data.frame(VarCorr(lmerDroughtRes2018))$vcov
  lmerDroughtVu2018 <- lmerDroughtVars2018[1]
  lmerDroughtVe2018 <- lmerDroughtVars2018[2]

  lmerDroughtHeritInd2018 <- lmerDroughtVu2018 / (lmerDroughtVu2018 + lmerDroughtVe2018)
  lmerDroughtHeritLine2018 <- lmerDroughtVu2018 / (lmerDroughtVu2018 + lmerDroughtVe2018 / nRep)

  lmerDroughtRanef2018 <- ranef(lmerDroughtRes2018)$variety[varietyNames, ]
  lmerResidDrought2018 <- residuals(lmerDroughtRes2018)

  lmerVarHeritsDrought2018[metabNo - metabStart + 1, ] <-
    c(lmerDroughtVars2018, lmerDroughtHeritInd2018, lmerDroughtHeritLine2018)
  lmerRanefMatDrought2018[, metabNo - metabStart + 1] <- lmerDroughtRanef2018
  lmerResidMatDrought2018[, metabNo - metabStart + 1] <- lmerResidDrought2018

}


upperThresDrought2018 <- apply(lmerResidMatDrought2018, 2, mean) + thresSd * apply(lmerResidMatDrought2018, 2, sd)
lowerThresDrought2018 <- apply(lmerResidMatDrought2018, 2, mean) - thresSd * apply(lmerResidMatDrought2018, 2, sd)


outlierDrought2018 <- sapply(1:nMetab, function(metabNo) {
  return((lmerResidMatDrought2018[, metabNo] > upperThresDrought2018[metabNo]) |
           (lmerResidMatDrought2018[, metabNo] < lowerThresDrought2018[metabNo]))
})

table(outlierDrought2018)



#### 4.3.4 Gather outlier information for 2018 data ####
colnames(outlier2018) <- colnames(outlierControl2018) <-
  colnames(outlierDrought2018) <- metabNames

outlierCD2018 <- rbind(outlierControl2018,
                       outlierDrought2018)

outlierAll2018 <- outlier2018 | outlierCD2018

table(apply(outlierAll2018, 2, sum))

metabOnly2018 <- metab2018Raw[, metabStart:metabEnd]
metabOnly2018[outlierAll2018] <- NA


metab2018NoOL <- metab2018Raw
metab2018NoOL[, metabStart:metabEnd] <- metabOnly2018


# fileNameNoOL2018 <- paste0("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier_thres=",
fileNameNoOL2018 <- paste0("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier", ".csv")
write.csv(x = metab2018NoOL, file = fileNameNoOL2018, quote = FALSE, row.names = FALSE)




##### 5. Extracting data of metabolites only related to Flavonoid pathway in 2018 into R #####
#### 5.1. Read Metabolome data without Outlier in 2018 into R ####
metab2018Raw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier.csv")
See(metab2018Raw, rown = 6, coln = 12)

OnlyMetab2018 <- metab2018Raw[, 10:ncol(metab2018Raw)]
See(OnlyMetab2018)


#### 5.2. Extracting data of metabolites only related to Flavonoid pathway in 2018 ####
CRESTMetabAnnotation <- read.xlsx("raw_data/extra/CREST_WT_265compounds_for_KEGG.xlsx")
metabOnlyRelatedToFlavonoidPathway <- CRESTMetabAnnotation[CRESTMetabAnnotation$XX == "a", ]
head(metabOnlyRelatedToFlavonoidPathway)
metabSampleIDRelatedToFlavonoidPathway <- metabOnlyRelatedToFlavonoidPathway[, 2]
metabAnnotationRelatedToFlavonoidPathway <- metabOnlyRelatedToFlavonoidPathway[, 3]

names(metabSampleIDRelatedToFlavonoidPathway) <- metabAnnotationRelatedToFlavonoidPathway


CommonMetabNames <- colnames(OnlyMetab2018[, colnames(OnlyMetab2018)%in% metabSampleIDRelatedToFlavonoidPathway])
table(colnames(OnlyMetab2018)%in% metabSampleIDRelatedToFlavonoidPathway)
FlavonoidMetab <- OnlyMetab2018[, colnames(OnlyMetab2018)%in% metabSampleIDRelatedToFlavonoidPathway]
See(FlavonoidMetab)

FlavonoidMetab <- cbind(metab2018Raw[, 1:9], FlavonoidMetab)
head(FlavonoidMetab)
See(FlavonoidMetab, coln = 15)


CommonAnnotation <- metabSampleIDRelatedToFlavonoidPathway[metabSampleIDRelatedToFlavonoidPathway %in% CommonMetabNames]
length(CommonAnnotation)
Annotation <- names(CommonAnnotation)
CommonMetabNamesAndAnnotationFlavo <- cbind(CommonMetabNames, Annotation)
colnames(CommonMetabNamesAndAnnotationFlavo)[1] <- "Name"

CommonMetabNamesAndAnnotationFlavo
See(CommonMetabNamesAndAnnotationFlavo)


fileNameMetabNoOutlierFlavonoid <- paste0("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway", ".csv")

fileNameFlavonoidMetabNamesAndAnnotation <- paste0("data/extra/2018_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation", ".csv")


write.csv(x = FlavonoidMetab, file = fileNameMetabNoOutlierFlavonoid, quote = FALSE, row.names = FALSE)

write.csv(x = CommonMetabNamesAndAnnotationFlavo, file = fileNameFlavonoidMetabNamesAndAnnotation, row.names = FALSE)
# write.csv(x = CommonMetabNamesAndAnnotationFlavo, file = fileNameFlavonoidMetabNamesAndAnnotation, quote = FALSE, row.names = FALSE)












