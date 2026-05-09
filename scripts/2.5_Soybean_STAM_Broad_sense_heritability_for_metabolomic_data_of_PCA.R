##########################################################################################
######  Title: 2.5_Soybean_STAM_Broad_sense_heritability_for_metabolomic_PCA_data   ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/09 (Created), 2026/02/19 (Last Updated)                       ######
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

scriptID <- "2.5"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMBSHPCA <- paste0(dirMidSTAMBase, scriptID,
                        "_BSH_for_PCA/")
dir.create(dirMidSTAMBSHPCA)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
require(data.table)
require(lme4)
# require(heritability)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### 2. Estimate broad-sense heritability for PC Score based on Metabolomic data in 2017 ######
##### 2.1. Read Metabolomic data in 2017 into R #####
metabPCA2017 <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_2017.csv")
See(metabPCA2017, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCA2017)
metabNames <- colnames(metabPCA2017)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCA2017$variety)
blockNames <- unique(metabPCA2017$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCA2017$ind))
nBlock <- length(blockNames)



##### 2.2. Perform lmer for X00006 in Metabolomic data in 2017 #####
#### 2.2.1 Perform lmer for X00006 in all Metabolomic data in 2017 ####
#lmerPCARes2017 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCA2017)
#summary(lmerPCARes2017)

#lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
#lmerPCAVu2017 <- lmerPCAVars2017[1]
#lmerPCAVe2017 <- lmerPCAVars2017[2]

#lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
#lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

#lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames,]





##### 2.3. Perform lmer for all of PC Score in 2017 #####
#### 2.3.1 Perform lmer for PC Score based on all Metabolomic data in 2017 ####
lmerPCAVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2017) <- metabNames
colnames(lmerPCAVarHerits2017) <- varNames
rownames(lmerPCARanefMat2017) <- varietyNames
colnames(lmerPCARanefMat2017) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2017[, metabNo]

  lmerPCARes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                      data = metabPCA2017)

  lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
  lmerPCAVu2017 <- lmerPCAVars2017[1]
  lmerPCAVe2017 <- lmerPCAVars2017[2]

  lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
  lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

  lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames, ]


  lmerPCAVarHerits2017[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2017, lmerPCAHeritInd2017, lmerPCAHeritLine2017)
  lmerPCARanefMat2017[, metabNo - metabStart + 1] <- lmerPCARanef2017
}





### write csv files for each treatment x herit, ranef
fileNamePCA2017 <- paste0(dirMidSTAMBSHPCA, scriptID, "_lmer_genotypic_values_for_PC_Score_in_2017.csv")
#fileNameControl2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2017.csv")
#fileNameDrought2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2017.csv")
#fileNameCPlusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2017.csv")
#fileNameCMinusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2017.csv")


write.csv(x = lmerPCARanefMat2017,
          file = fileNamePCA2017)

#write.csv(x = lmerRanefMatControl2017,
#          file = fileNameControl2017)

#write.csv(x = lmerRanefMatDrought2017,
#          file = fileNameDrought2017)

#write.csv(x = lmerRanefMatControl2017 + lmerRanefMatDrought2017,
#          file = fileNameCPlusD2017)

#write.csv(x = lmerRanefMatControl2017 - lmerRanefMatDrought2017,
#          file = fileNameCMinusD2017)



#plot(round(apply(lmerRanefMatControl2017, 2, mean), 5),
#     round(apply(lmerRanefMatDrought2017, 2, mean), 5))
#plot(round(apply(lmerRanefMatControl2017, 2, var), 5),
#     round(apply(lmerRanefMatDrought2017, 2, var), 5))
#which((apply(lmerRanefMatControl2017, 2, var) < 1000) &
#        (apply(lmerRanefMatDrought2017, 2, var) > 10000))
#lmerRanefMatCPlusD2017 <- lmerRanefMatControl2017 + lmerRanefMatDrought2017
#lmerRanefMatCMinusD2017 <- lmerRanefMatControl2017 - lmerRanefMatDrought2017


#lmerVarHeritsTotalControlDrought2017<- cbind(lmerVarHerits2017[, 3],
#lmerVarHeritsControl2017[, 3],
#lmerVarHeritsDrought2017[, 3])




#### 2.3.2. Heritability ####
### individual heritability
pdf(paste0(dirMidSTAMBSHPCA, scriptID, "_Individual_heritability_PCA_for_188_metab_in_2017.pdf"))
hist(lmerPCAVarHerits2017[, 3], xlim = range(0:1), main = "Individual heritability", xlab = "heritability" )
dev.off()
#hist(lmerVarHeritsControl2017[, 3], xlim = range(0:1))
#hist(lmerVarHeritsDrought2017[, 3], xlim = range(0:1))


### variety heritability
pdf(paste0(dirMidSTAMBSHPCA, scriptID, "_Variety_heritability_PCA_for_188_metab_in_2017.pdf"))
hist(lmerPCAVarHerits2017[, 4], xlim = range(0:1), main = "Variety heritability", xlab = "heritability")
dev.off()
#hist(lmerVarHeritsControl2017[, 4], xlim = range(0:1))
#hist(lmerVarHeritsDrought2017[, 4], xlim = range(0:1))




#### 2.3.3. Histogram of PCs ####
gvPCMetab2017 <- read.csv("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_for_PC_Score_in_2017.csv", row.names = 1)
See(gvPCMetab2017)

metabNames <- colnames(gvPCMetab2017)
nMetab <- ncol(gvPCMetab2017)


### Total
dir.create(paste0(dirMidSTAMBSHPCA, scriptID, "_Histogram_of_BLUP_for_PCs_in_2017/"))
dirMidSTAMBSHPCABLUPHist2017 <- paste0(dirMidSTAMBSHPCA, scriptID, "_Histogram_of_BLUP_for_PCs_in_2017/")

dir.create(paste0(dirMidSTAMBSHPCABLUPHist2017, scriptID, "_Total/"))

for(metabNo in 1:nMetab){
  metabNow <- gvPCMetab2017[, metabNo]
  metabName <- metabNames[metabNo]

  pdf(paste0(dirMidSTAMBSHPCABLUPHist2017, scriptID, "_Total/", scriptID, "_", metabName, ".pdf"))
  hist(metabNow, xlim = c(min(metabNow), max(metabNow)))
  # hist(metabNow, freq = FALSE, probability = TRUE)
  dev.off()

}






##### 2.4. nPC=100 ####
# Read data
nPC <- 100

metabPCA2017 <- read.csv(paste0("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_nPC=", nPC, "_2017.csv"))
See(metabPCA2017, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCA2017)
metabNames <- colnames(metabPCA2017)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCA2017$variety)
blockNames <- unique(metabPCA2017$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCA2017$ind))
nBlock <- length(blockNames)


# Perform lmer for all of PC Score in 2017
lmerPCAVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2017) <- metabNames
colnames(lmerPCAVarHerits2017) <- varNames
rownames(lmerPCARanefMat2017) <- varietyNames
colnames(lmerPCARanefMat2017) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2017[, metabNo]

  lmerPCARes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCA2017)

  lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
  lmerPCAVu2017 <- lmerPCAVars2017[1]
  lmerPCAVe2017 <- lmerPCAVars2017[2]

  lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
  lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

  lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames, ]


  lmerPCAVarHerits2017[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2017, lmerPCAHeritInd2017, lmerPCAHeritLine2017)
  lmerPCARanefMat2017[, metabNo - metabStart + 1] <- lmerPCARanef2017
}


# write csv files for each treatment x herit, ranef
fileNamePCA2017 <- paste0(dirMidSTAMBSHPCA, scriptID, "_lmer_genotypic_values_for_PC_Score_nPC=", nPC, "_in_2017.csv")

write.csv(x = lmerPCARanefMat2017,
          file = fileNamePCA2017)







##### 2018
###### 3. Estimate broad-sense heritability for PC Score based on Metabolomic data in 2018 ######
##### 3.1. Read Metabolomic data in 2018 into R #####
metabPCA2018 <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_2018.csv")
See(metabPCA2018, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCA2018)
metabNames <- colnames(metabPCA2018)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCA2018$variety)
blockNames <- unique(metabPCA2018$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCA2018$ind))
nBlock <- length(blockNames)



##### 3.2. Perform lmer for X00006 in Metabolomic data in 2018 #####
#### 3.2.1 Perform lmer for X00006 in all Metabolomic data in 2018 ####
#lmerPCARes2018 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCA2018)
#summary(lmerPCARes2018)

#lmerPCAVars2018 <- data.frame(VarCorr(lmerPCARes2018))$vcov
#lmerPCAVu2018 <- lmerPCAVars2018[1]
#lmerPCAVe2018 <- lmerPCAVars2018[2]

#lmerPCAHeritInd2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018)
#lmerPCAHeritLine2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018 / (nRep * nBlock))

#lmerPCARanef2018 <- ranef(lmerPCARes2018)$variety[varietyNames,]





##### 3.3. Perform lmer for all of PC Score in 2018 #####
#### 3.3.1 Perform lmer for PC Score base on all Metabolomic data in 2018 ####
lmerPCAVarHerits2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2018) <- metabNames
colnames(lmerPCAVarHerits2018) <- varNames
rownames(lmerPCARanefMat2018) <- varietyNames
colnames(lmerPCARanefMat2018) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2018[, metabNo]

  lmerPCARes2018 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCA2018)

  lmerPCAVars2018 <- data.frame(VarCorr(lmerPCARes2018))$vcov
  lmerPCAVu2018 <- lmerPCAVars2018[1]
  lmerPCAVe2018 <- lmerPCAVars2018[2]

  lmerPCAHeritInd2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018)
  lmerPCAHeritLine2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018 / (nRep * nBlock))

  lmerPCARanef2018 <- ranef(lmerPCARes2018)$variety[varietyNames, ]


  lmerPCAVarHerits2018[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2018, lmerPCAHeritInd2018, lmerPCAHeritLine2018)
  lmerPCARanefMat2018[, metabNo - metabStart + 1] <- lmerPCARanef2018
}





### write csv files for each treatment x herit, ranef
fileNamePCA2018 <- paste0(dirMidSTAMBSHPCA, scriptID, "_lmer_genotypic_values_for_PC_Score_in_2018.csv")
#fileNameControl2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2018.csv")
#fileNameDrought2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2018.csv")
#fileNameCPlusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2018.csv")
#fileNameCMinusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2018.csv")


write.csv(x = lmerPCARanefMat2018,
          file = fileNamePCA2018)

#write.csv(x = lmerRanefMatControl2018,
#          file = fileNameControl2018)

#write.csv(x = lmerRanefMatDrought2018,
#          file = fileNameDrought2018)

#write.csv(x = lmerRanefMatControl2018 + lmerRanefMatDrought2018,
#          file = fileNameCPlusD2018)

#write.csv(x = lmerRanefMatControl2018 - lmerRanefMatDrought2018,
#          file = fileNameCMinusD2018)



#plot(round(apply(lmerRanefMatControl2018, 2, mean), 5),
#     round(apply(lmerRanefMatDrought2018, 2, mean), 5))
#plot(round(apply(lmerRanefMatControl2018, 2, var), 5),
#     round(apply(lmerRanefMatDrought2018, 2, var), 5))
#which((apply(lmerRanefMatControl2018, 2, var) < 1000) &
#        (apply(lmerRanefMatDrought2018, 2, var) > 10000))
#lmerRanefMatCPlusD2018 <- lmerRanefMatControl2018 + lmerRanefMatDrought2018
#lmerRanefMatCMinusD2018 <- lmerRanefMatControl2018 - lmerRanefMatDrought2018


#lmerVarHeritsTotalControlDrought2018<- cbind(lmerVarHerits2018[, 3],
#lmerVarHeritsControl2018[, 3],
#lmerVarHeritsDrought2018[, 3])




### individual heritability
pdf(paste0(dirMidSTAMBSHPCA, scriptID, "_Individual_heritability_PCA_for_188_metab_in_2018.pdf"))
hist(lmerPCAVarHerits2018[, 3], xlim = range(0:1), main = "Individual heritability", xlab = "heritability" )
dev.off()
#hist(lmerVarHeritsControl2018[, 3], xlim = range(0:1))
#hist(lmerVarHeritsDrought2018[, 3], xlim = range(0:1))


### variety heritability
pdf(paste0(dirMidSTAMBSHPCA, scriptID, "_Variety_heritability_PCA_for_188_metab_in_2018.pdf"))
hist(lmerPCAVarHerits2018[, 4], xlim = range(0:1), main = "Variety heritability", xlab = "heritability")
dev.off()
#hist(lmerVarHeritsControl2018[, 4], xlim = range(0:1))
#hist(lmerVarHeritsDrought2018[, 4], xlim = range(0:1))






##### 2017 + 2018


#########################





###### 5. Arithmetic mean of Metabolomic data in 2017 ######
##### 5.1. Arithmetic mean of all Metabolomic data in 2017 #####
metabAM2017 <- apply(X = metab2017Raw[, metabStart:metabEnd],
                     MARGIN = 2,
                     FUN = function (metabNow) {
                       tapply(X = metabNow,
                              INDEX = metab2017Raw$variety,
                              FUN = mean, na.rm = TRUE)[varietyNames]
                     })
metabAMScaled2017 <- scale(metabAM2017, center = T, scale = F)

##### 5.2. Arithmetic mean of Metabolomic data for control in 2017 #####
metabAMControl2017 <- apply(X = metabControl2017[, metabStart:metabEnd],
                            MARGIN = 2,
                            FUN = function (metabNow) {
                              tapply(X = metabNow,
                                     INDEX = metabControl2017$variety,
                                     FUN = mean, na.rm = TRUE)[varietyNames]
                            })
metabAMControlScaled2017 <- scale(metabAMControl2017, center = T, scale = F)

##### 5.3. Arithmetic mean of Metabolomic data for drought in 2017 #####
metabAMDrought2017 <- apply(X = metabDrought2017[, metabStart:metabEnd],
                            MARGIN = 2,
                            FUN = function (metabNow) {
                              tapply(X = metabNow,
                                     INDEX = metabDrought2017$variety,
                                     FUN = mean, na.rm = TRUE)[varietyNames]
                            })
metabAMDroughtScaled2017 <- scale(metabAMDrought2017, center = T, scale = F)






###### 6. Estimate broad-sense heritability for PC Score based on Root Metabolomic data in 2020 ######
##### 6.1. Read Metabolomic data in 2020 into R #####
nPC <- 20

metabPCA2020 <- read.csv(paste0("midstream/2.3_PCA/2.3_pcaMethods_PCA_nPC_", nPC, "_root_metab_2020.csv"))
See(metabPCA2020, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCA2020)
metabNames <- colnames(metabPCA2020)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCA2020$variety)
blockNames <- unique(metabPCA2020$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- 1
nBlock <- length(blockNames)



##### 6.2. Perform lmer for X00006 in Metabolomic data in 2020 #####
#### 6.2.1 Perform lmer for X00006 in all Metabolomic data in 2020 ####
#lmerPCARes2020 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCA2020)
#summary(lmerPCARes2020)

#lmerPCAVars2020 <- data.frame(VarCorr(lmerPCARes2020))$vcov
#lmerPCAVu2020 <- lmerPCAVars2020[1]
#lmerPCAVe2020 <- lmerPCAVars2020[2]

#lmerPCAHeritInd2020 <- lmerPCAVu2020 / (lmerPCAVu2020 + lmerPCAVe2020)
#lmerPCAHeritLine2020 <- lmerPCAVu2020 / (lmerPCAVu2020 + lmerPCAVe2020 / (nRep * nBlock))

#lmerPCARanef2020 <- ranef(lmerPCARes2020)$variety[varietyNames,]





##### 6.3. Perform lmer for all of PC Score for Root Metabolomic data in 2020 #####
#### 6.3.1 Perform lmer for PC Score based on all Root Metabolomic data in 2020 ####
lmerPCAVarHerits2020 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2020 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2020) <- metabNames
colnames(lmerPCAVarHerits2020) <- varNames
rownames(lmerPCARanefMat2020) <- varietyNames
colnames(lmerPCARanefMat2020) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2020[, metabNo]

  lmerPCARes2020 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCA2020)

  lmerPCAVars2020 <- data.frame(VarCorr(lmerPCARes2020))$vcov
  lmerPCAVu2020 <- lmerPCAVars2020[1]
  lmerPCAVe2020 <- lmerPCAVars2020[2]

  lmerPCAHeritInd2020 <- lmerPCAVu2020 / (lmerPCAVu2020 + lmerPCAVe2020)
  lmerPCAHeritLine2020 <- lmerPCAVu2020 / (lmerPCAVu2020 + lmerPCAVe2020 / (nRep * nBlock))

  lmerPCARanef2020 <- ranef(lmerPCARes2020)$variety[varietyNames, ]


  lmerPCAVarHerits2020[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2020, lmerPCAHeritInd2020, lmerPCAHeritLine2020)
  lmerPCARanefMat2020[, metabNo - metabStart + 1] <- lmerPCARanef2020
}

See(lmerPCARanefMat2020)



#### 6.3.2 Perform lmer for PC Score based on Control Metabolomic data in 2020 ####
metabPCA2020 <- metabPCA2020[metabPCA2020$block == "W1", ]

lmerPCAVarHerits2020 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2020 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2020) <- metabNames
colnames(lmerPCAVarHerits2020) <- varNames
rownames(lmerPCARanefMat2020) <- varietyNames
colnames(lmerPCARanefMat2020) <- metabNames

### NG, for no biological reprecations
for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2020[, metabNo]

  lmerPCARes2020 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCA2020)

  lmerPCAVars2020 <- data.frame(VarCorr(lmerPCARes2020))$vcov
  lmerPCAVu2020 <- lmerPCAVars2020[1]
  lmerPCAVe2020 <- lmerPCAVars2020[2]

  lmerPCAHeritInd2020 <- lmerPCAVu2020 / (lmerPCAVu2020 + lmerPCAVe2020)
  lmerPCAHeritLine2020 <- lmerPCAVu2020 / (lmerPCAVu2020 + lmerPCAVe2020 / (nRep * nBlock))

  lmerPCARanef2020 <- ranef(lmerPCARes2020)$variety[varietyNames, ]


  lmerPCAVarHerits2020[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2020, lmerPCAHeritInd2020, lmerPCAHeritLine2020)
  lmerPCARanefMat2020[, metabNo - metabStart + 1] <- lmerPCARanef2020
}

See(lmerPCARanefMat2020)






### write csv files for each treatment x herit, ranef
fileNamePCA2020 <- paste0(dirMidSTAMBSHPCA, scriptID, "_lmer_genotypic_values_for_PC_Score_based_on_root_metabolome_data_in_2020_nPC_", nPC, ".csv")
#fileNameControl2020 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2020.csv")
#fileNameDrought2020 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2020.csv")
#fileNameCPlusD2020 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2020.csv")
#fileNameCMinusD2020 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2020.csv")


write.csv(x = lmerPCARanefMat2020,
          file = fileNamePCA2020)

#write.csv(x = lmerRanefMatControl2020,
#          file = fileNameControl2020)

#write.csv(x = lmerRanefMatDrought2020,
#          file = fileNameDrought2020)

#write.csv(x = lmerRanefMatControl2020 + lmerRanefMatDrought2020,
#          file = fileNameCPlusD2020)

#write.csv(x = lmerRanefMatControl2020 - lmerRanefMatDrought2020,
#          file = fileNameCMinusD2020)



#plot(round(apply(lmerRanefMatControl2020, 2, mean), 5),
#     round(apply(lmerRanefMatDrought2020, 2, mean), 5))
#plot(round(apply(lmerRanefMatControl2020, 2, var), 5),
#     round(apply(lmerRanefMatDrought2020, 2, var), 5))
#which((apply(lmerRanefMatControl2020, 2, var) < 1000) &
#        (apply(lmerRanefMatDrought2020, 2, var) > 10000))
#lmerRanefMatCPlusD2020 <- lmerRanefMatControl2020 + lmerRanefMatDrought2020
#lmerRanefMatCMinusD2020 <- lmerRanefMatControl2020 - lmerRanefMatDrought2020


#lmerVarHeritsTotalControlDrought2020<- cbind(lmerVarHerits2020[, 3],
#lmerVarHeritsControl2020[, 3],
#lmerVarHeritsDrought2020[, 3])


# individual heritability
hist(lmerPCAVarHerits2020[, 3], xlim = range(0:1), breaks = 10)
#hist(lmerVarHeritsControl2020[, 3], xlim = range(0:1))
#hist(lmerVarHeritsDrought2020[, 3], xlim = range(0:1))

# variety heritability
hist(lmerPCAVarHerits2020[, 4], xlim = range(0:1), breaks = 10)
#hist(lmerVarHeritsControl2020[, 4], xlim = range(0:1))
#hist(lmerVarHeritsDrought2020[, 4], xlim = range(0:1))




###### 10. Perform lmer for all of PC Score based on prcomp() and ppca() in 2017 ######
##### 10.1. Using PC Scores by prcomp() #####
#### 10.1.1 Read Metabolomic data in 2017 into R ####
metabPCA2017Prcomp <- read.csv("midstream/2.3_PCA/2.3_prcomp_PCA_metab_without_NA_2017.csv")
See(metabPCA2017Prcomp, rown = 3, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCA2017Prcomp)
metabNames <- colnames(metabPCA2017Prcomp)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCA2017Prcomp$variety)
blockNames <- unique(metabPCA2017Prcomp$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCA2017Prcomp$ind))
nBlock <- length(blockNames)



#### 10.1.2 Perform lmer for PC Score based on all Metabolomic data in 2017 ####
lmerPCAVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2017) <- metabNames
colnames(lmerPCAVarHerits2017) <- varNames
rownames(lmerPCARanefMat2017) <- varietyNames
colnames(lmerPCARanefMat2017) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2017Prcomp[, metabNo]

  lmerPCARes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCA2017Prcomp)

  lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
  lmerPCAVu2017 <- lmerPCAVars2017[1]
  lmerPCAVe2017 <- lmerPCAVars2017[2]

  lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
  lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

  lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames, ]


  lmerPCAVarHerits2017[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2017, lmerPCAHeritInd2017, lmerPCAHeritLine2017)
  lmerPCARanefMat2017[, metabNo - metabStart + 1] <- lmerPCARanef2017
}


### write csv files for each treatment x herit, ranef
fileNamePCA2017 <- paste0(dirMidSTAMBSHPCA, scriptID, "_lmer_genotypic_values_for_PC_Score_by_prcomp_based_on_data_without_NA_in_2017.csv")

write.csv(x = lmerPCARanefMat2017,
          file = fileNamePCA2017)



##### 10.2. Using PC Scores by ppca() #####
metabPCA2017Ppca <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_without_NA_metab_2017.csv")
See(metabPCA2017Ppca, rown = 3, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCA2017Ppca)
metabNames <- colnames(metabPCA2017Ppca)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCA2017Ppca$variety)
blockNames <- unique(metabPCA2017Ppca$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCA2017Ppca$ind))
nBlock <- length(blockNames)



#### 10.1.2 Perform lmer for PC Score based on all Metabolomic data in 2017 ####
lmerPCAVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2017) <- metabNames
colnames(lmerPCAVarHerits2017) <- varNames
rownames(lmerPCARanefMat2017) <- varietyNames
colnames(lmerPCARanefMat2017) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2017Ppca[, metabNo]

  lmerPCARes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCA2017Ppca)

  lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
  lmerPCAVu2017 <- lmerPCAVars2017[1]
  lmerPCAVe2017 <- lmerPCAVars2017[2]

  lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
  lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

  lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames, ]


  lmerPCAVarHerits2017[metabNo - metabStart + 1, ] <-
    c(lmerPCAVars2017, lmerPCAHeritInd2017, lmerPCAHeritLine2017)
  lmerPCARanefMat2017[, metabNo - metabStart + 1] <- lmerPCARanef2017
}



### write csv files for each treatment x herit, ranef
fileNamePCA2017 <- paste0(dirMidSTAMBSHPCA, scriptID, "_lmer_genotypic_values_for_PC_Score_by_ppca_based_on_data_without_NA_in_2017.csv")

write.csv(x = lmerPCARanefMat2017,
          file = fileNamePCA2017)










########################
##### 2019


##### 2019 + 2020


#########################




