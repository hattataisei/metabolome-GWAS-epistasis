##########################################################################################
######  Title: 2.18_Soybean_STAM_Broad_sense_heritability_for_metabolomic_PCA_data_with_or_without_Flavonoid   ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/04 (Created), 2025/01/11 (Last Updated)                       ######
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

scriptID <- "2.18"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMBSHPCAFlavonoid <- paste0(dirMidSTAMBase, scriptID,
                           "_BSH_for_PCA_with_or_without_flavonoid/")
dir.create(dirMidSTAMBSHPCAFlavonoid)
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





#### 2017
###### 2. Estimate broad-sense heritability for PC Score based on Flavonoid-related Metabolomic data in 2017 ######
##### 2.1. Read Flavonoid-related Metabolomic data in 2017 into R #####
metabPCAFlavonoid2017 <- read.csv("midstream/2.16_PCA_with_or_without_flavonoid/2.16_pcaMethods_PCA_metab_related_to_flavonoid_pathway_2017.csv")
See(metabPCAFlavonoid2017, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCAFlavonoid2017)
metabNames <- colnames(metabPCAFlavonoid2017)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCAFlavonoid2017$variety)
blockNames <- unique(metabPCAFlavonoid2017$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCAFlavonoid2017$ind))
nBlock <- length(blockNames)



##### 2.2. Perform lmer for X00006 in Metabolomic data in 2017 #####
#### 2.2.1 Perform lmer for X00006 in all Metabolomic data in 2017 ####
#lmerPCARes2017 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCAFlavonoid2017)
#summary(lmerPCARes2017)

#lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
#lmerPCAVu2017 <- lmerPCAVars2017[1]
#lmerPCAVe2017 <- lmerPCAVars2017[2]

#lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
#lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

#lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames,]





##### 2.3. Perform lmer for all of PC Score in 2017 #####
#### 2.3.1 Perform lmer for PC Score base on Flavonoid-related Metabolomic data in 2017 ####
lmerPCAVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2017) <- metabNames
colnames(lmerPCAVarHerits2017) <- varNames
rownames(lmerPCARanefMat2017) <- varietyNames
colnames(lmerPCARanefMat2017) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCAFlavonoid2017[, metabNo]

  lmerPCARes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCAFlavonoid2017)

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

See(lmerPCARanefMat2017)


See(lmerPCAVarHerits2017)



pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Individual_heritability_83_metab_for_flavonoid-related_2017.pdf"))
hist(lmerPCAVarHerits2017[, "heritInd"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Flavonoid-related")
dev.off()


pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Variety_heritability_83_metab_for_flavonoid-related_2017.pdf"))
hist(lmerPCAVarHerits2017[, "heritLine"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Flavonoid-related")
dev.off()






### write csv files for each treatment x herit, ranef
fileNamePCAFlavonoid2017 <- paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_lmer_genotypic_values_for_PC_Score_for_flavonoid_related_metab_in_2017.csv")
#fileNameControl2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2017.csv")
#fileNameDrought2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2017.csv")
#fileNameCPlusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2017.csv")
#fileNameCMinusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2017.csv")


write.csv(x = lmerPCARanefMat2017,
          file = fileNamePCAFlavonoid2017)

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




#### 2.3.2. Histogram of PCs ####
gvPCMetab2017 <- read.csv("midstream/2.18_BSH_for_PCA_with_or_without_flavonoid/2.18_lmer_genotypic_values_for_PC_Score_for_flavonoid_related_metab_in_2017.csv", row.names = 1)
See(gvPCMetab2017)

metabNames <- colnames(gvPCMetab2017)
nMetab <- ncol(gvPCMetab2017)


### Total
dir.create(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Histogram_of_BLUP_for_PCs_Flavonoid_in_2017/"))
dirMidSTAMBSHPCABLUPFlavonoidHist2017 <- paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Histogram_of_BLUP_for_PCs_Flavonoid_in_2017/")

dir.create(paste0(dirMidSTAMBSHPCABLUPFlavonoidHist2017, scriptID, "_Total/"))

for(metabNo in 1:nMetab){
  metabNow <- gvPCMetab2017[, metabNo]
  metabName <- metabNames[metabNo]

  pdf(paste0(dirMidSTAMBSHPCABLUPFlavonoidHist2017, scriptID, "_Total/", scriptID, "_", metabName, ".pdf"))
  hist(metabNow, xlim = c(min(metabNow), max(metabNow)))
  # hist(metabNow, freq = FALSE, probability = TRUE)
  dev.off()

}







###### 3. Estimate broad-sense heritability for PC Score based on Flavonoid Non-related Metabolomic data in 2017 ######
##### 3.1. Read Flavonoid Non-related Metabolomic data in 2017 into R #####
metabPCAWithoutFlavonoid2017 <- read.csv("midstream/2.16_PCA_with_or_without_flavonoid/2.16_pcaMethods_PCA_metab_Not_related_to_flavonoid_pathway_2017.csv")
See(metabPCAWithoutFlavonoid2017, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCAWithoutFlavonoid2017)
metabNames <- colnames(metabPCAWithoutFlavonoid2017)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCAWithoutFlavonoid2017$variety)
blockNames <- unique(metabPCAWithoutFlavonoid2017$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCAWithoutFlavonoid2017$ind))
nBlock <- length(blockNames)



##### 3.2. Perform lmer for X00006 in Metabolomic data in 2017 #####
#### 3.2.1 Perform lmer for X00006 in all Metabolomic data in 2017 ####
#lmerPCARes2017 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCAWithoutFlavonoid2017)
#summary(lmerPCARes2017)

#lmerPCAVars2017 <- data.frame(VarCorr(lmerPCARes2017))$vcov
#lmerPCAVu2017 <- lmerPCAVars2017[1]
#lmerPCAVe2017 <- lmerPCAVars2017[2]

#lmerPCAHeritInd2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017)
#lmerPCAHeritLine2017 <- lmerPCAVu2017 / (lmerPCAVu2017 + lmerPCAVe2017 / (nRep * nBlock))

#lmerPCARanef2017 <- ranef(lmerPCARes2017)$variety[varietyNames,]





##### 3.3. Perform lmer for all of PC Score in 2017 #####
#### 3.3.1 Perform lmer for PC Score base on Flavonoid Non-related Metabolomic data in 2017 ####
lmerPCAVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2017) <- metabNames
colnames(lmerPCAVarHerits2017) <- varNames
rownames(lmerPCARanefMat2017) <- varietyNames
colnames(lmerPCARanefMat2017) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCAWithoutFlavonoid2017[, metabNo]

  lmerPCARes2017 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCAWithoutFlavonoid2017)

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

See(lmerPCARanefMat2017)
See(lmerPCAVarHerits2017)



pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Individual_heritability_83_metab_for_Non_flavonoid-related_2017.pdf"))
hist(lmerPCAVarHerits2017[, "heritInd"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Non Flavonoid-related")
dev.off()

pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Variety_heritability_83_metab_for_Non_flavonoid-related_2017.pdf"))
hist(lmerPCAVarHerits2017[, "heritLine"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Non Flavonoid-related")
dev.off()





### write csv files for each treatment x herit, ranef
fileNamePCAWithoutFlavonoid2017 <- paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_lmer_genotypic_values_for_PC_Score_for_flavonoid_non_related_metab_in_2017.csv")
#fileNameControl2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2017.csv")
#fileNameDrought2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2017.csv")
#fileNameCPlusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2017.csv")
#fileNameCMinusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2017.csv")


write.csv(x = lmerPCARanefMat2017,
          file = fileNamePCAWithoutFlavonoid2017)

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


# individual heritability
hist(lmerPCAVarHerits2017[, 3], xlim = range(0:1), breaks = 10)
#hist(lmerVarHeritsControl2017[, 3], xlim = range(0:1))
#hist(lmerVarHeritsDrought2017[, 3], xlim = range(0:1))

# variety heritability
hist(lmerPCAVarHerits2017[, 4], xlim = range(0:1), breaks = 10)
#hist(lmerVarHeritsControl2017[, 4], xlim = range(0:1))
#hist(lmerVarHeritsDrought2017[, 4], xlim = range(0:1))






##### 2018
###### 4. Estimate broad-sense heritability for PC Score based on Flavonoid-related Metabolomic data in 2018 ######
##### 4.1. Read Flavonoid-related Metabolomic data in 2018 into R #####
metabPCAFlavonoid2018 <- read.csv("midstream/2.16_PCA_with_or_without_flavonoid/2.16_pcaMethods_PCA_metab_related_to_flavonoid_pathway_2018.csv")
See(metabPCAFlavonoid2018, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCAFlavonoid2018)
metabNames <- colnames(metabPCAFlavonoid2018)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCAFlavonoid2018$variety)
blockNames <- unique(metabPCAFlavonoid2018$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCAFlavonoid2018$ind))
nBlock <- length(blockNames)



##### 4.2. Perform lmer for X00006 in Metabolomic data in 2018 #####
#### 4.2.1 Perform lmer for X00006 in all Metabolomic data in 2018 ####
#lmerPCARes2018 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCAFlavonoid2018)
#summary(lmerPCARes2018)

#lmerPCAVars2018 <- data.frame(VarCorr(lmerPCARes2018))$vcov
#lmerPCAVu2018 <- lmerPCAVars2018[1]
#lmerPCAVe2018 <- lmerPCAVars2018[2]

#lmerPCAHeritInd2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018)
#lmerPCAHeritLine2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018 / (nRep * nBlock))

#lmerPCARanef2018 <- ranef(lmerPCARes2018)$variety[varietyNames,]





##### 4.3. Perform lmer for all of PC Score in 2018 #####
#### 4.3.1 Perform lmer for PC Score base on Flavonoid-related Metabolomic data in 2018 ####
lmerPCAVarHerits2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2018) <- metabNames
colnames(lmerPCAVarHerits2018) <- varNames
rownames(lmerPCARanefMat2018) <- varietyNames
colnames(lmerPCARanefMat2018) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCAFlavonoid2018[, metabNo]

  lmerPCARes2018 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCAFlavonoid2018)

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

See(lmerPCARanefMat2018)


See(lmerPCAVarHerits2018)



pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Individual_heritability_83_metab_for_flavonoid-related_2018.pdf"))
hist(lmerPCAVarHerits2018[, "heritInd"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Flavonoid-related")
dev.off()


pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Variety_heritability_83_metab_for_flavonoid-related_2018.pdf"))
hist(lmerPCAVarHerits2018[, "heritLine"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Flavonoid-related")
dev.off()






### write csv files for each treatment x herit, ranef
fileNamePCAFlavonoid2018 <- paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_lmer_genotypic_values_for_PC_Score_for_flavonoid_related_metab_in_2018.csv")
#fileNameControl2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2018.csv")
#fileNameDrought2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2018.csv")
#fileNameCPlusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2018.csv")
#fileNameCMinusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2018.csv")


write.csv(x = lmerPCARanefMat2018,
          file = fileNamePCAFlavonoid2018)

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








###### 5. Estimate broad-sense heritability for PC Score based on Flavonoid Non-related Metabolomic data in 2018 ######
##### 5.1. Read Flavonoid Non-related Metabolomic data in 2018 into R #####
metabPCAWithoutFlavonoid2018 <- read.csv("midstream/2.16_PCA_with_or_without_flavonoid/2.16_pcaMethods_PCA_metab_Not_related_to_flavonoid_pathway_2018.csv")
See(metabPCAWithoutFlavonoid2018, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metabPCAWithoutFlavonoid2018)
metabNames <- colnames(metabPCAWithoutFlavonoid2018)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metabPCAWithoutFlavonoid2018$variety)
blockNames <- unique(metabPCAWithoutFlavonoid2018$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metabPCAWithoutFlavonoid2018$ind))
nBlock <- length(blockNames)



##### 5.2. Perform lmer for X00006 in Metabolomic data in 2018 #####
#### 5.2.1 Perform lmer for X00006 in all Metabolomic data in 2018 ####
#lmerPCARes2018 <- lmer(formula = X00006 ? (1 | variety) + block, data = metabPCAWithoutFlavonoid2018)
#summary(lmerPCARes2018)

#lmerPCAVars2018 <- data.frame(VarCorr(lmerPCARes2018))$vcov
#lmerPCAVu2018 <- lmerPCAVars2018[1]
#lmerPCAVe2018 <- lmerPCAVars2018[2]

#lmerPCAHeritInd2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018)
#lmerPCAHeritLine2018 <- lmerPCAVu2018 / (lmerPCAVu2018 + lmerPCAVe2018 / (nRep * nBlock))

#lmerPCARanef2018 <- ranef(lmerPCARes2018)$variety[varietyNames,]





##### 5.3. Perform lmer for all of PC Score in 2018 #####
#### 5.3.1 Perform lmer for PC Score base on Flavonoid Non-related Metabolomic data in 2018 ####
lmerPCAVarHerits2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerPCARanefMat2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerPCAVarHerits2018) <- metabNames
colnames(lmerPCAVarHerits2018) <- varNames
rownames(lmerPCARanefMat2018) <- varietyNames
colnames(lmerPCARanefMat2018) <- metabNames


for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCAWithoutFlavonoid2018[, metabNo]

  lmerPCARes2018 <- lmer(formula = metabNow ? (1 | variety) + block,
                         data = metabPCAWithoutFlavonoid2018)

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

See(lmerPCARanefMat2018)
See(lmerPCAVarHerits2018)



pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Individual_heritability_83_metab_for_Non_flavonoid-related_2018.pdf"))
hist(lmerPCAVarHerits2018[, "heritInd"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Non Flavonoid-related")
dev.off()

pdf(paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_Variety_heritability_83_metab_for_Non_flavonoid-related_2018.pdf"))
hist(lmerPCAVarHerits2018[, "heritLine"], xlim = c(0, 1),breaks = seq(0, 1, 0.1), xlab = "Heritability", main = "Non Flavonoid-related")
dev.off()





### write csv files for each treatment x herit, ranef
fileNamePCAWithoutFlavonoid2018 <- paste0(dirMidSTAMBSHPCAFlavonoid, scriptID, "_lmer_genotypic_values_for_PC_Score_for_flavonoid_non_related_metab_in_2018.csv")
#fileNameControl2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2018.csv")
#fileNameDrought2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2018.csv")
#fileNameCPlusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2018.csv")
#fileNameCMinusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2018.csv")


write.csv(x = lmerPCARanefMat2018,
          file = fileNamePCAWithoutFlavonoid2018)

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


# individual heritability
hist(lmerPCAVarHerits2018[, 3], xlim = range(0:1), breaks = 10)
#hist(lmerVarHeritsControl2018[, 3], xlim = range(0:1))
#hist(lmerVarHeritsDrought2018[, 3], xlim = range(0:1))

# variety heritability
hist(lmerPCAVarHerits2018[, 4], xlim = range(0:1), breaks = 10)
#hist(lmerVarHeritsControl2018[, 4], xlim = range(0:1))
#hist(lmerVarHeritsDrought2018[, 4], xlim = range(0:1))






##### 2017 + 2018


#########################





###### 6. Arithmetic mean of Metabolomic data in 2017 ######
##### 6.1. Arithmetic mean of all Metabolomic data in 2017 #####
metabAM2017 <- apply(X = metab2017Raw[, metabStart:metabEnd],
                     MARGIN = 2,
                     FUN = function (metabNow) {
                       tapply(X = metabNow,
                              INDEX = metab2017Raw$variety,
                              FUN = mean, na.rm = TRUE)[varietyNames]
                     })
metabAMScaled2017 <- scale(metabAM2017, center = T, scale = F)

##### 6.2. Arithmetic mean of Metabolomic data for control in 2017 #####
metabAMControl2017 <- apply(X = metabControl2017[, metabStart:metabEnd],
                            MARGIN = 2,
                            FUN = function (metabNow) {
                              tapply(X = metabNow,
                                     INDEX = metabControl2017$variety,
                                     FUN = mean, na.rm = TRUE)[varietyNames]
                            })
metabAMControlScaled2017 <- scale(metabAMControl2017, center = T, scale = F)

##### 6.3. Arithmetic mean of Metabolomic data for drought in 2017 #####
metabAMDrought2017 <- apply(X = metabDrought2017[, metabStart:metabEnd],
                            MARGIN = 2,
                            FUN = function (metabNow) {
                              tapply(X = metabNow,
                                     INDEX = metabDrought2017$variety,
                                     FUN = mean, na.rm = TRUE)[varietyNames]
                            })
metabAMDroughtScaled2017 <- scale(metabAMDrought2017, center = T, scale = F)











