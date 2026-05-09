##########################################################################################
######  Title: 2.2_Soybean_STAM_Broad_sense_heritability_for_metabolomic_data       ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/09 (Created), 2026/01/13 (Last Updated)                       ######
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

scriptID <- "2.2"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMBSH <- paste0(dirMidSTAMBase, scriptID, "_BSH/")
dir.create(dirMidSTAMBSH)

# dirMidSTAMBSHlmerGv <- paste0(dirMidSTAMBase, scriptID, "_BSH/", scriptID, "_lmer_genotypic_values_2017/")
# dir.create(dirMidSTAMBSHlmerGv)

dirMidSTAMBSHHeritHist2017 <- paste0(dirMidSTAMBase, scriptID, "_BSH/", scriptID, "_Heritability_of_individual_and_line_2017/")
dir.create(dirMidSTAMBSHHeritHist2017)
dirMidSTAMBSHHeritHist2018 <- paste0(dirMidSTAMBase, scriptID, "_BSH/", scriptID, "_Heritability_of_individual_and_line_2018/")
dir.create(dirMidSTAMBSHHeritHist2018)

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




###### 2. Estimate broad-sense heritability for Metabolomic data in 2017, using No outlier raw data ######
##### 2.1. Read Metabolomic data in 2017 into R #####
metab2017Raw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
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



##### 2.2. Perform lmer for X00006 in Metabolomic data in 2017 #####
#### 2.2.1 Perform lmer for X00006 in all Metabolomic data in 2017 ####
# lmerRes2017 <- lmer(formula = X00006 ? (1 | variety) + block,
#                     data = metab2017Raw)
# summary(lmerRes2017)
#
# lmerVars2017 <- data.frame(VarCorr(lmerRes2017))$vcov
# lmerVu2017 <- lmerVars2017[1]
# lmerVe2017 <- lmerVars2017[2]
#
# lmerHeritInd2017 <- lmerVu2017 / (lmerVu2017 + lmerVe2017)
# lmerHeritLine2017 <- lmerVu2017 / (lmerVu2017 + lmerVe2017 / (nRep * nBlock))
#
# lmerRanef2017 <- ranef(lmerRes2017)$variety[varietyNames,]
#
#
#
# #### 2.2.2 Perform lmer for X00006 in Metabolomic data for control in 2017 ####
# metabControl2017 <- metab2017Raw[metab2017Raw$block == "C", ]
# lmerControlRes2017 <- lmer(formula = X00006 ? (1 | variety),
#                            data = metabControl2017)
#
# lmerControlVars2017 <- data.frame(VarCorr(lmerControlRes2017))$vcov
# lmerControlVu2017 <- lmerControlVars2017[1]
# lmerControlVe2017 <- lmerControlVars2017[2]
#
# lmerControlHeritInd2017 <- lmerControlVu2017 / (lmerControlVu2017 + lmerControlVe2017)
# lmerControlHeritLine2017 <- lmerControlVu2017 / (lmerControlVu2017 + lmerControlVe2017 / nRep)
#
# lmerControlRanef2017 <- ranef(lmerControlRes2017)$variety[varietyNames, ]
#
#
# #### 2.2.3  Perform lmer for X00006 in Metabolomic data for drought in 2017 ####
# metabDrought2017 <- metab2017Raw[metab2017Raw$block == "D", ]
# lmerDroughtRes2017 <- lmer(formula = X00006 ? (1 | variety),
#                            data = metabDrought2017)
#
# lmerDroughtVars2017 <- data.frame(VarCorr(lmerDroughtRes2017))$vcov
# lmerDroughtVu2017 <- lmerDroughtVars2017[1]
# lmerDroughtVe2017 <- lmerDroughtVars2017[2]
#
# lmerDroughtHeritInd2017 <- lmerDroughtVu2017 / (lmerDroughtVu2017 + lmerDroughtVe2017)
# lmerDroughtHeritLine2017 <- lmerDroughtVu2017 / (lmerDroughtVu2017 + lmerDroughtVe2017 / nRep)
#
# lmerDroughtRanef2017 <- ranef(lmerDroughtRes2017)$variety[varietyNames, ]



##### 2.3. Perform lmer for all Metabolomic data in 2017 #####
#### 2.3.1 Perform lmer for all Metabolomic data in 2017 ####
lmerVarHerits2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

### extra
lmerInterceptPlusRanefMat2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)



rownames(lmerVarHerits2017) <- metabNames
colnames(lmerVarHerits2017) <- varNames
rownames(lmerRanefMat2017) <- varietyNames
colnames(lmerRanefMat2017) <- metabNames

rownames(lmerInterceptPlusRanefMat2017) <- varietyNames
colnames(lmerInterceptPlusRanefMat2017) <- metabNames


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


  lmerVarHerits2017[metabNo - metabStart + 1, ] <-
    c(lmerVars2017, lmerHeritInd2017, lmerHeritLine2017)
  lmerRanefMat2017[, metabNo - metabStart + 1] <- lmerRanef2017


  lmerInterceptPlusRanefMat2017[, metabNo - metabStart + 1] <- coef(lmerRes2017)$variety[varietyNames, 1]

}


See(lmerInterceptPlusRanefMat2017, coln = 188)
table(lmerInterceptPlusRanefMat2017 > 0)
table(apply(lmerInterceptPlusRanefMat2017, 2, function(x)any(x < 0)))
apply(lmerInterceptPlusRanefMat2017, 2, min)

metabNamesMinMinus <- metabNames[apply(lmerInterceptPlusRanefMat2017, 2, function(x)any(x < 0))]
metabNamesMinMinus <- na.omit(metabNamesMinMinus)
lmerInterceptPlusRanefMat2017[, metabNamesMinMinus]



#### 2.3.2 Perform lmer for all Metabolomic data for control in 2017 ####
metabControl2017 <- metab2017Raw[metab2017Raw$block == "C", ]

lmerVarHeritsControl2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatControl2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsControl2017) <- metabNames
colnames(lmerVarHeritsControl2017) <- varNames
rownames(lmerRanefMatControl2017) <- varietyNames
colnames(lmerRanefMatControl2017) <- metabNames



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


  lmerVarHeritsControl2017[metabNo - metabStart + 1, ] <-
    c(lmerControlVars2017, lmerControlHeritInd2017, lmerControlHeritLine2017)
  lmerRanefMatControl2017[, metabNo - metabStart + 1] <- lmerControlRanef2017
}



#### 2.3.3  Perform lmer for all Metabolomic data for drought in 2017 ####
metabDrought2017 <- metab2017Raw[metab2017Raw$block == "D", ]

lmerVarHeritsDrought2017 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatDrought2017 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsDrought2017) <- metabNames
colnames(lmerVarHeritsDrought2017) <- varNames
rownames(lmerRanefMatDrought2017) <- varietyNames
colnames(lmerRanefMatDrought2017) <- metabNames



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


  lmerVarHeritsDrought2017[metabNo - metabStart + 1, ] <-
    c(lmerDroughtVars2017, lmerDroughtHeritInd2017, lmerDroughtHeritLine2017)
  lmerRanefMatDrought2017[, metabNo - metabStart + 1] <- lmerDroughtRanef2017
}


lmerRanefMatCPlusD2017 <- lmerRanefMatControl2017 + lmerRanefMatDrought2017
lmerRanefMatCMinusD2017 <- lmerRanefMatControl2017 - lmerRanefMatDrought2017



### write csv files for each treatment x herit, ranef
fileNameTotal2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_2017.csv")
fileNameControl2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2017.csv")
fileNameDrought2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2017.csv")
fileNameCPlusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2017.csv")
fileNameCMinusD2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2017.csv")

fileNameTotalInterceptPlusRanef2017 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_intercept_plus_ranef_2017.csv")


write.csv(x = lmerRanefMat2017,
          file = fileNameTotal2017)
write.csv(x = lmerRanefMatControl2017,
          file = fileNameControl2017)
write.csv(x = lmerRanefMatDrought2017,
          file = fileNameDrought2017)
write.csv(x = lmerRanefMatCPlusD2017,
          file = fileNameCPlusD2017)
write.csv(x = lmerRanefMatCMinusD2017,
          file = fileNameCMinusD2017)

write.csv(x = lmerInterceptPlusRanefMat2017,
          file = fileNameTotalInterceptPlusRanef2017)




plot(round(apply(lmerRanefMatControl2017, 2, var), 5),
     round(apply(lmerRanefMatDrought2017, 2, var), 5))
plot(round(apply(lmerRanefMatControl2017, 2, mean), 5),
     round(apply(lmerRanefMatDrought2017, 2, mean), 5))
which((apply(lmerRanefMatControl2017, 2, var) < 1000) &
        (apply(lmerRanefMatDrought2017, 2, var) > 10000))


### ratio of metabolites ( heritability > 0.5 )
nrow(lmerVarHerits2017[lmerVarHerits2017[, 4] > 0.5, ] ) / nrow(lmerVarHerits2017)
nrow(lmerVarHeritsControl2017[lmerVarHeritsControl2017[, 4] > 0.5, ] ) / nrow(lmerVarHeritsControl2017)
nrow(lmerVarHeritsDrought2017[lmerVarHeritsDrought2017[, 4] > 0.5, ] ) / nrow(lmerVarHeritsDrought2017)



#### 2.3.4.  Extract flavonoid > 0.9 line heritability and data, and write csv files ####
metab2017FlavonoidRaw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
See(metab2017FlavonoidRaw, rown = 6, coln = 12)
metab2017FlavonoidNames <- colnames(metab2017FlavonoidRaw[, 10:ncol(metab2017FlavonoidRaw)])
lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017[metab2017FlavonoidNames, ]
See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017FlavonoidMoreThan0.9Heritability[lmerVarHerits2017FlavonoidMoreThan0.9Heritability$heritLine > 0.9, ]
See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)

varietyNames <- rownames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017FlavonoidMoreThan0.9Heritability$heritLine
names(lmerVarHerits2017FlavonoidMoreThan0.9Heritability) <- varietyNames
lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
colnames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability) <- "heritLine"
See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)


onlyMetab2017Flavonoid <- metab2017FlavonoidRaw[, 10:ncol(metab2017FlavonoidRaw)]
onlyMetabFlavonoidMoreThan0.9Heritability <- onlyMetab2017Flavonoid[, rownames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)]
See(onlyMetabFlavonoidMoreThan0.9Heritability)
metab2017FlavonoidMoreThan0.9HeritabilityRaw <- cbind(metab2017FlavonoidRaw[, 1:9], onlyMetabFlavonoidMoreThan0.9Heritability)
See(metab2017FlavonoidMoreThan0.9HeritabilityRaw, coln = 12)



fileNameFlavonoidAnnotationMoreThan0.9Heritability2017 <- paste0(dirMidSTAMBSH, scriptID, "_Flavonoid_>0.9_heritability_2017.csv")
fileNameFlavonoidMoreThan0.9Heritability2017Raw <- paste0("data/phenotype/", "2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")


write.csv(x = lmerVarHerits2017FlavonoidMoreThan0.9Heritability, file = fileNameFlavonoidAnnotationMoreThan0.9Heritability2017)
write.csv(x = metab2017FlavonoidMoreThan0.9HeritabilityRaw, file = fileNameFlavonoidMoreThan0.9Heritability2017Raw)



##### 2.4.  Histogram of Broad sense heritability in 2017 #####
#### 2.4.1 Total ####
See(lmerVarHerits2017)
class(lmerVarHerits2017)

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_individual_for_Total.pdf"))
hist(lmerVarHerits2017[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

class(lmerVarHerits2017)
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_line_for_Total.pdf"))
hist(lmerVarHerits2017[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()



#### 2.4.2 Control ####
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_individual_for_Control.pdf"))
hist(lmerVarHeritsControl2017[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_line_for_Control.pdf"))
hist(lmerVarHeritsControl2017[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()




#### 2.4.3 Drought ####
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_individual_for_Drought.pdf"))
hist(lmerVarHeritsControl2017[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_line_for_Drought.pdf"))
hist(lmerVarHeritsDrought2017[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()




#### 2.4.4 Flavonoid ####
flavonoidMetab <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
See(flavonoidMetab)
flavonoidMetabNames <- flavonoidMetab[, 1]

lmerVarHerits2017DF <- as.data.frame((lmerVarHerits2017))
lmerVarHerits2017FlavonoidDF <- lmerVarHerits2017DF[flavonoidMetabNames, ]
See(lmerVarHerits2017FlavonoidDF)
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_for_Total.pdf"))
hist(lmerVarHerits2017FlavonoidDF$heritLine, ylim = c(0,50), xlab =  "_Heritability of Line", main = "Flavonoid")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_for_Total.pdf"))
hist(lmerVarHerits2017FlavonoidDF$heritInd, ylim = c(0,70), xlab = "Heritability of Individual", main = "Flavonoid")
dev.off()

# metabHeritsLineHigher2017 <- lmerVarHerits2017DF[lmerVarHerits2017DF$heritLine > 0.9, ]



#### 2.4.5 Not Flavonoid ####
flavonoidMetab <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
See(flavonoidMetab)
flavonoidMetabNames <- flavonoidMetab[, 1]

lmerVarHerits2017DF <- as.data.frame((lmerVarHerits2017))
lmerVarHerits2017NonFlavonoidDF <- lmerVarHerits2017DF[(!rownames(lmerVarHerits2017DF) %in% flavonoidMetabNames), ]
See(lmerVarHerits2017NonFlavonoidDF)
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_non-related_for_Total.pdf"))
hist(lmerVarHerits2017NonFlavonoidDF$heritLine, xlim = c(0,1), ylim = c(0,50), xlab = "Heritability of Line", main = "Flavonoid non-related")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_non-related_for_Total.pdf"))
hist(lmerVarHerits2017NonFlavonoidDF$heritInd, xlim = c(0,1), ylim = c(0,70), breaks = seq(0,1,0.1), xlab = "Heritability of Individual", main = "Flavonoid non-related")
dev.off()




# col.has.na <- apply(lmerRanefMatDrought2017, 2, function(x){any(is.na(x))})
# which(col.has.na)
# table(col.has.na)
#
# lmerRanefMat2017NACol <- lmerRanefMat2017[, col.has.na]
# lmerRanefMat2017NACol
#
# metab2017Raw["UA-4805", ]
#
#
# row.has.na <- apply(lmerRanefMat2017, 1, function(x){any(is.na(x))})
# which(row.has.na)




##### 2.5. Histogram of BLUP for each variety of metabolites in 2017 #####
gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
See(gvMetab2017Total)
table(is.na(gvMetab2017Total))
gvMetab2017Total <- na.omit(gvMetab2017Total)
See(gvMetab2017Total)

metabNames <- colnames(gvMetab2017Total)
nMetab <- ncol(gvMetab2017Total)

### to edit
# metab2017Control <- metab2017Raw[metab2017Raw$block == "C", ]
# See(metab2017Control, coln = 10)
# metab2017ControlOnlyMetab <- metab2017Control[, metabStart:metabEnd]
# See(metab2017ControlOnlyMetab)
#
# metab2017Drought <- metab2017Raw[metab2017Raw$block == "D", ]
# metab2017DroughtOnlyMetab <- metab2017Drought[, metabStart:metabEnd]
###

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])


### Total
dir.create(paste0(dirMidSTAMBSH, scriptID, "_Histogram_of_BLUP_for_each_metabolite_in_2017/"))
dirMidSTAMBSHBLUPHist2017 <- paste0(dirMidSTAMBSH, scriptID, "_Histogram_of_BLUP_for_each_metabolite_in_2017/")

dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/"))
dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid/"))
dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Non_Flavonoid/"))
dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid_heritability_>0.9/"))


metabNo <- 1
for(metabNo in 1:nMetab){
  metabNow <- gvMetab2017Total[, metabNo]
  metabName <- metabNames[metabNo]

  if (metabName %in% metabNamesFlavonoid){
    pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
    hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
    # hist(metabNow, freq = FALSE, probability = TRUE)
    dev.off()

  } else {
    pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
    hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
    # hist(metabNow, freq = FALSE, probability = TRUE)
    dev.off()

  }
}


for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
  metabNow <- gvMetab2017Total[, metabName]

  pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
  hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
  dev.off()

}



### to edit
# ### Control
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017ControlOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
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
#   pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
#
#
#
# ### Drought
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017DroughtOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
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
#   pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
###












##### 2018
###### 3. Estimate broad-sense heritability for Metabolomic data in 2018, using No outlier raw data ######
##### 3.1. Read Metabolomic data in 2018 into R #####
metab2018Raw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier.csv")
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



##### 3.2. Perform lmer for X00006 in Metabolomic data in 2018 #####
#### 3.2.1 Perform lmer for X00006 in all Metabolomic data in 2018 ####



##### 3.3. Perform lmer for all Metabolomic data in 2018 #####
#### 3.3.1 Perform lmer for all Metabolomic data in 2018 ####
lmerVarHerits2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMat2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHerits2018) <- metabNames
colnames(lmerVarHerits2018) <- varNames
rownames(lmerRanefMat2018) <- varietyNames
colnames(lmerRanefMat2018) <- metabNames


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


  lmerVarHerits2018[metabNo - metabStart + 1, ] <-
    c(lmerVars2018, lmerHeritInd2018, lmerHeritLine2018)
  lmerRanefMat2018[, metabNo - metabStart + 1] <- lmerRanef2018
}




#### 3.3.2 Perform lmer for all Metabolomic data for control in 2018 ####
metabControl2018 <- metab2018Raw[metab2018Raw$block == "C", ]

lmerVarHeritsControl2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatControl2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsControl2018) <- metabNames
colnames(lmerVarHeritsControl2018) <- varNames
rownames(lmerRanefMatControl2018) <- varietyNames
colnames(lmerRanefMatControl2018) <- metabNames



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


  lmerVarHeritsControl2018[metabNo - metabStart + 1, ] <-
    c(lmerControlVars2018, lmerControlHeritInd2018, lmerControlHeritLine2018)
  lmerRanefMatControl2018[, metabNo - metabStart + 1] <- lmerControlRanef2018
}




#### 3.3.3  Perform lmer for all Metabolomic data for drought in 2018 ####
metabDrought2018 <- metab2018Raw[metab2018Raw$block == "D", ]

lmerVarHeritsDrought2018 <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatDrought2018 <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsDrought2018) <- metabNames
colnames(lmerVarHeritsDrought2018) <- varNames
rownames(lmerRanefMatDrought2018) <- varietyNames
colnames(lmerRanefMatDrought2018) <- metabNames



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


  lmerVarHeritsDrought2018[metabNo - metabStart + 1, ] <-
    c(lmerDroughtVars2018, lmerDroughtHeritInd2018, lmerDroughtHeritLine2018)
  lmerRanefMatDrought2018[, metabNo - metabStart + 1] <- lmerDroughtRanef2018
}


lmerRanefMatCPlusD2018 <- lmerRanefMatControl2018 + lmerRanefMatDrought2018
lmerRanefMatCMinusD2018 <- lmerRanefMatControl2018 - lmerRanefMatDrought2018



### write csv files for each treatment x herit, ranef
fileNameTotal2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_2018.csv")
fileNameControl2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2018.csv")
fileNameDrought2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2018.csv")
fileNameCPlusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2018.csv")
fileNameCMinusD2018 <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2018.csv")


write.csv(x = lmerRanefMat2018,
          file = fileNameTotal2018)

write.csv(x = lmerRanefMatControl2018,
          file = fileNameControl2018)

write.csv(x = lmerRanefMatDrought2018,
          file = fileNameDrought2018)

write.csv(x = lmerRanefMatCPlusD2018,
          file = fileNameCPlusD2018)

write.csv(x = lmerRanefMatCMinusD2018,
          file = fileNameCMinusD2018)



plot(round(apply(lmerRanefMatControl2018, 2, var), 5),
     round(apply(lmerRanefMatDrought2018, 2, var), 5))
plot(round(apply(lmerRanefMatControl2018, 2, mean), 5),
     round(apply(lmerRanefMatDrought2018, 2, mean), 5))
which((apply(lmerRanefMatControl2018, 2, var) < 1000) &
        (apply(lmerRanefMatDrought2018, 2, var) > 10000))


### ratio of metabolites ( heritability > 0.5 )
nrow(lmerVarHerits2018[lmerVarHerits2018[, 4] > 0.5, ] ) / nrow(lmerVarHerits2018)
nrow(lmerVarHeritsControl2018[lmerVarHeritsControl2018[, 4] > 0.5, ] ) / nrow(lmerVarHeritsControl2018)
nrow(lmerVarHeritsDrought2018[lmerVarHeritsDrought2018[, 4] > 0.5, ] ) / nrow(lmerVarHeritsDrought2018)
See(lmerVarHerits2018)
See(lmerVarHeritsControl2018)
See(lmerVarHeritsDrought2018)


#### 3.3.4.  Extract flavonoid > 0.9 line heritability and data, and write csv files ####
metab2018FlavonoidRaw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
See(metab2018FlavonoidRaw, rown = 6, coln = 12)
metab2018FlavonoidNames <- colnames(metab2018FlavonoidRaw[, 10:ncol(metab2018FlavonoidRaw)])
lmerVarHerits2018FlavonoidMoreThan0.9Heritability <- lmerVarHerits2018[metab2018FlavonoidNames, ]
See(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)
lmerVarHerits2018FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)
lmerVarHerits2018FlavonoidMoreThan0.9Heritability <- lmerVarHerits2018FlavonoidMoreThan0.9Heritability[lmerVarHerits2018FlavonoidMoreThan0.9Heritability$heritLine > 0.9, ]
See(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)

varietyNames <- rownames(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)
lmerVarHerits2018FlavonoidMoreThan0.9Heritability <- lmerVarHerits2018FlavonoidMoreThan0.9Heritability$heritLine
names(lmerVarHerits2018FlavonoidMoreThan0.9Heritability) <- varietyNames
lmerVarHerits2018FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)
colnames(lmerVarHerits2018FlavonoidMoreThan0.9Heritability) <- "heritLine"
See(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)


onlyMetab2018Flavonoid <- metab2018FlavonoidRaw[, 10:ncol(metab2018FlavonoidRaw)]
onlyMetabFlavonoidMoreThan0.9Heritability <- onlyMetab2018Flavonoid[, rownames(lmerVarHerits2018FlavonoidMoreThan0.9Heritability)]
See(onlyMetabFlavonoidMoreThan0.9Heritability)
metab2018FlavonoidMoreThan0.9HeritabilityRaw <- cbind(metab2018FlavonoidRaw[, 1:9], onlyMetabFlavonoidMoreThan0.9Heritability)
See(metab2018FlavonoidMoreThan0.9HeritabilityRaw, coln = 12)



fileNameFlavonoidAnnotationMoreThan0.9Heritability2018 <- paste0(dirMidSTAMBSH, scriptID, "_Flavonoid_>0.9_heritability_2018.csv")
fileNameFlavonoidMoreThan0.9Heritability2018Raw <- paste0("data/phenotype/", "2018_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")


write.csv(x = lmerVarHerits2018FlavonoidMoreThan0.9Heritability, file = fileNameFlavonoidAnnotationMoreThan0.9Heritability2018)
write.csv(x = metab2018FlavonoidMoreThan0.9HeritabilityRaw, file = fileNameFlavonoidMoreThan0.9Heritability2018Raw)



##### 3.4.  Histogram of Broad sense heritability in 2018 #####
#### 3.4.1 Total ####
See(lmerVarHerits2018)
class(lmerVarHerits2018)

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_individual_for_Total.pdf"))
hist(lmerVarHerits2018[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

class(lmerVarHerits2018)
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_line_for_Total.pdf"))
hist(lmerVarHerits2018[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()



#### 3.4.2 Control ####
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_individual_for_Control.pdf"))
hist(lmerVarHeritsControl2018[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_line_for_Control.pdf"))
hist(lmerVarHeritsControl2018[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()



#### 3.4.3 Drought ####
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_individual_for_Drought.pdf"))
hist(lmerVarHeritsControl2018[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_line_for_Drought.pdf"))
hist(lmerVarHeritsDrought2018[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()



#### 3.4.4 Flavonoid ####
flavonoidMetab <- read.csv("data/extra/2018_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
See(flavonoidMetab)
flavonoidMetabNames <- flavonoidMetab[, 1]

lmerVarHerits2018DF <- as.data.frame((lmerVarHerits2018))
lmerVarHerits2018FlavonoidDF <- lmerVarHerits2018DF[flavonoidMetabNames, ]
See(lmerVarHerits2018FlavonoidDF)
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid.pdf"))
hist(lmerVarHerits2018FlavonoidDF$heritLine, ylim = c(0,50), xlab =  "_Heritability of Line", main = "Flavonoid")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid.pdf"))
hist(lmerVarHerits2018FlavonoidDF$heritInd, ylim = c(0,70), xlab = "Heritability of Individual", main = "Flavonoid")
dev.off()

# metabHeritsLineHigher2018 <- lmerVarHerits2018DF[lmerVarHerits2018DF$heritLine > 0.9, ]
?hist


#### 3.4.5 Not Flavonoid ####
flavonoidMetab <- read.csv("data/extra/2018_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
See(flavonoidMetab)
flavonoidMetabNames <- flavonoidMetab[, 1]

lmerVarHerits2018DF <- as.data.frame((lmerVarHerits2018))
lmerVarHerits2018NonFlavonoidDF <- lmerVarHerits2018DF[(!rownames(lmerVarHerits2018DF) %in% flavonoidMetabNames), ]
See(lmerVarHerits2018NonFlavonoidDF)
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_non-related.pdf"))
hist(lmerVarHerits2018NonFlavonoidDF$heritLine, xlim = c(0,1), ylim = c(0,50), xlab = "Heritability of Line", main = "Flavonoid non-related")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_non-related.pdf"))
hist(lmerVarHerits2018NonFlavonoidDF$heritInd, xlim = c(0,1), ylim = c(0,70), breaks = seq(0,1,0.1), xlab = "Heritability of Individual", main = "Flavonoid non-related")
dev.off()











###### 4. Estimate broad-sense heritability for Metabolomic data in 2017, using Box-Cox data ######
##### 4.1. Read Metabolomic Box-Cox data in 2017 into R #####
metab2017BoxCox <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_BoxCox.csv", row.names = 1)
See(metab2017BoxCox, rown = 6, coln = 12)


metabStart <- 10
metabEnd <- ncol(metab2017BoxCox)
metabNames <- colnames(metab2017BoxCox)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2017BoxCox$variety)
blockNames <- unique(metab2017BoxCox$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2017BoxCox$ind))
nBlock <- length(blockNames)



##### 4.2. Perform lmer for X00006 in Metabolomic data in 2017 #####
#### 4.2.1 Perform lmer for X00006 in all Metabolomic data in 2017 ####
# lmerRes2017 <- lmer(formula = X00006 ? (1 | variety) + block,
#                     data = metab2017Raw)
# summary(lmerRes2017)
#
# lmerVars2017 <- data.frame(VarCorr(lmerRes2017))$vcov
# lmerVu2017 <- lmerVars2017[1]
# lmerVe2017 <- lmerVars2017[2]
#
# lmerHeritInd2017 <- lmerVu2017 / (lmerVu2017 + lmerVe2017)
# lmerHeritLine2017 <- lmerVu2017 / (lmerVu2017 + lmerVe2017 / (nRep * nBlock))
#
# lmerRanef2017 <- ranef(lmerRes2017)$variety[varietyNames,]
#
#
#
# #### 4.2.2 Perform lmer for X00006 in Metabolomic data for control in 2017 ####
# metabControl2017 <- metab2017Raw[metab2017Raw$block == "C", ]
# lmerControlRes2017 <- lmer(formula = X00006 ? (1 | variety),
#                            data = metabControl2017)
#
# lmerControlVars2017 <- data.frame(VarCorr(lmerControlRes2017))$vcov
# lmerControlVu2017 <- lmerControlVars2017[1]
# lmerControlVe2017 <- lmerControlVars2017[2]
#
# lmerControlHeritInd2017 <- lmerControlVu2017 / (lmerControlVu2017 + lmerControlVe2017)
# lmerControlHeritLine2017 <- lmerControlVu2017 / (lmerControlVu2017 + lmerControlVe2017 / nRep)
#
# lmerControlRanef2017 <- ranef(lmerControlRes2017)$variety[varietyNames, ]
#
#
# #### 4.2.3  Perform lmer for X00006 in Metabolomic data for drought in 2017 ####
# metabDrought2017 <- metab2017Raw[metab2017Raw$block == "D", ]
# lmerDroughtRes2017 <- lmer(formula = X00006 ? (1 | variety),
#                            data = metabDrought2017)
#
# lmerDroughtVars2017 <- data.frame(VarCorr(lmerDroughtRes2017))$vcov
# lmerDroughtVu2017 <- lmerDroughtVars2017[1]
# lmerDroughtVe2017 <- lmerDroughtVars2017[2]
#
# lmerDroughtHeritInd2017 <- lmerDroughtVu2017 / (lmerDroughtVu2017 + lmerDroughtVe2017)
# lmerDroughtHeritLine2017 <- lmerDroughtVu2017 / (lmerDroughtVu2017 + lmerDroughtVe2017 / nRep)
#
# lmerDroughtRanef2017 <- ranef(lmerDroughtRes2017)$variety[varietyNames, ]





##### 4.3. Perform lmer for all Metabolomic data in 2017 #####
#### 4.3.1 Perform lmer for all Metabolomic data in 2017 ####
lmerVarHerits2017BoxCox <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMat2017BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)

### extra
lmerInterceptPlusRanefMat2017BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)
###



rownames(lmerVarHerits2017BoxCox) <- metabNames
colnames(lmerVarHerits2017BoxCox) <- varNames
rownames(lmerRanefMat2017BoxCox) <- varietyNames
colnames(lmerRanefMat2017BoxCox) <- metabNames


### extra
rownames(lmerInterceptPlusRanefMat2017BoxCox) <- varietyNames
colnames(lmerInterceptPlusRanefMat2017BoxCox) <- metabNames
###


for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2017BoxCox[, metabNo]

  lmerRes2017BoxCox <- lmer(formula = metabNow ? (1 | variety) + block,
                      data = metab2017BoxCox)

  lmerVars2017BoxCox <- data.frame(VarCorr(lmerRes2017BoxCox))$vcov
  lmerVu2017BoxCox <- lmerVars2017BoxCox[1]
  lmerVe2017BoxCox <- lmerVars2017BoxCox[2]


  lmerHeritInd2017BoxCox <- lmerVu2017BoxCox / (lmerVu2017BoxCox + lmerVe2017BoxCox)
  lmerHeritLine2017BoxCox <- lmerVu2017BoxCox / (lmerVu2017BoxCox + lmerVe2017BoxCox / (nRep * nBlock))

  lmerRanef2017BoxCox <- ranef(lmerRes2017BoxCox)$variety[varietyNames, ]


  lmerVarHerits2017BoxCox[metabNo - metabStart + 1, ] <-
    c(lmerVars2017BoxCox, lmerHeritInd2017BoxCox, lmerHeritLine2017BoxCox)
  lmerRanefMat2017BoxCox[, metabNo - metabStart + 1] <- lmerRanef2017BoxCox

  ###
  lmerInterceptPlusRanefMat2017BoxCox[, metabNo - metabStart + 1] <- coef(lmerRes2017BoxCox)$variety[varietyNames, 1]
  ###

}


###
See(lmerInterceptPlusRanefMat2017BoxCox, coln = 188)
table(lmerInterceptPlusRanefMat2017BoxCox > 0)
table(apply(lmerInterceptPlusRanefMat2017BoxCox, 2, function(x)any(x < 0)))
apply(lmerInterceptPlusRanefMat2017BoxCox, 2, min)

metabNamesMinMinus <- metabNames[apply(lmerInterceptPlusRanefMat2017BoxCox, 2, function(x)any(x < 0))]
metabNamesMinMinus <- na.omit(metabNamesMinMinus)
lmerInterceptPlusRanefMat2017BoxCox[, metabNamesMinMinus]
###




#### 4.3.2 Perform lmer for all Metabolomic data for control in 2017 ####
metabControl2017BoxCox <- metab2017BoxCox[metab2017BoxCox$block == "C", ]

lmerVarHeritsControl2017BoxCox <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatControl2017BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsControl2017BoxCox) <- metabNames
colnames(lmerVarHeritsControl2017BoxCox) <- varNames
rownames(lmerRanefMatControl2017BoxCox) <- varietyNames
colnames(lmerRanefMatControl2017BoxCox) <- metabNames



for (metabNo in metabStart:metabEnd){
  metabControlNow <- metabControl2017BoxCox[, metabNo]

  lmerControlRes2017BoxCox <- lmer(formula = metabControlNow ? (1 | variety),
                             data = metabControl2017BoxCox)

  lmerControlVars2017BoxCox <- data.frame(VarCorr(lmerControlRes2017BoxCox))$vcov
  lmerControlVu2017BoxCox <- lmerControlVars2017BoxCox[1]
  lmerControlVe2017BoxCox <- lmerControlVars2017BoxCox[2]

  lmerControlHeritInd2017BoxCox <- lmerControlVu2017BoxCox / (lmerControlVu2017BoxCox + lmerControlVe2017BoxCox)
  lmerControlHeritLine2017BoxCox <- lmerControlVu2017BoxCox / (lmerControlVu2017BoxCox + lmerControlVe2017BoxCox / nRep)

  lmerControlRanef2017BoxCox <- ranef(lmerControlRes2017BoxCox)$variety[varietyNames, ]


  lmerVarHeritsControl2017BoxCox[metabNo - metabStart + 1, ] <-
    c(lmerControlVars2017BoxCox, lmerControlHeritInd2017BoxCox, lmerControlHeritLine2017BoxCox)
  lmerRanefMatControl2017BoxCox[, metabNo - metabStart + 1] <- lmerControlRanef2017BoxCox
}



#### 4.3.3  Perform lmer for all Metabolomic data for drought in 2017 ####
metabDrought2017BoxCox <- metab2017BoxCox[metab2017BoxCox$block == "D", ]

lmerVarHeritsDrought2017BoxCox <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatDrought2017BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsDrought2017BoxCox) <- metabNames
colnames(lmerVarHeritsDrought2017BoxCox) <- varNames
rownames(lmerRanefMatDrought2017BoxCox) <- varietyNames
colnames(lmerRanefMatDrought2017BoxCox) <- metabNames



for (metabNo in metabStart:metabEnd){
  metabDroughtNow <- metabDrought2017BoxCox[, metabNo]

  lmerDroughtRes2017BoxCox <- lmer(formula = metabDroughtNow ? (1 | variety),
                             data = metabDrought2017BoxCox)

  lmerDroughtVars2017BoxCox <- data.frame(VarCorr(lmerDroughtRes2017BoxCox))$vcov
  lmerDroughtVu2017BoxCox <- lmerDroughtVars2017BoxCox[1]
  lmerDroughtVe2017BoxCox <- lmerDroughtVars2017BoxCox[2]

  lmerDroughtHeritInd2017BoxCox <- lmerDroughtVu2017BoxCox / (lmerDroughtVu2017BoxCox + lmerDroughtVe2017BoxCox)
  lmerDroughtHeritLine2017BoxCox <- lmerDroughtVu2017BoxCox / (lmerDroughtVu2017BoxCox + lmerDroughtVe2017BoxCox / nRep)

  lmerDroughtRanef2017BoxCox <- ranef(lmerDroughtRes2017BoxCox)$variety[varietyNames, ]


  lmerVarHeritsDrought2017BoxCox[metabNo - metabStart + 1, ] <-
    c(lmerDroughtVars2017BoxCox, lmerDroughtHeritInd2017BoxCox, lmerDroughtHeritLine2017BoxCox)
  lmerRanefMatDrought2017BoxCox[, metabNo - metabStart + 1] <- lmerDroughtRanef2017BoxCox
}


###
lmerRanefMatCPlusD2017BoxCox <- lmerRanefMatControl2017 + lmerRanefMatDrought2017
lmerRanefMatCMinusD2017BoxCox <- lmerRanefMatControl2017 - lmerRanefMatDrought2017



### write csv files for each treatment x herit, ranef
fileNameTotal2017BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_2017_BoxCox.csv")
fileNameControl2017BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2017_BoxCox.csv")
fileNameDrought2017BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2017_BoxCox.csv")
fileNameCPlusD2017BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2017_BoxCox.csv")
fileNameCMinusD2017BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2017_BoxCox.csv")

fileNameTotalInterceptPlusRanef2017BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_intercept_plus_ranef_2017_BoxCox.csv")


write.csv(x = lmerRanefMat2017BoxCox,
          file = fileNameTotal2017BoxCox)
write.csv(x = lmerRanefMatControl2017BoxCox,
          file = fileNameControl2017BoxCox)
write.csv(x = lmerRanefMatDrought2017BoxCox,
          file = fileNameDrought2017BoxCox)
write.csv(x = lmerRanefMatCPlusD2017BoxCox,
          file = fileNameCPlusD2017BoxCox)
write.csv(x = lmerRanefMatCMinusD2017BoxCox,
          file = fileNameCMinusD2017BoxCox)

write.csv(x = lmerInterceptPlusRanefMat2017BoxCox,
          file = fileNameTotalInterceptPlusRanef2017BoxCox)




plot(round(apply(lmerRanefMatControl2017BoxCox, 2, var), 5),
     round(apply(lmerRanefMatDrought2017BoxCox, 2, var), 5))
plot(round(apply(lmerRanefMatControl2017BoxCox, 2, mean), 5),
     round(apply(lmerRanefMatDrought2017BoxCox, 2, mean), 5))
which((apply(lmerRanefMatControl2017BoxCox, 2, var) < 1e+10) &
        (apply(lmerRanefMatDrought2017BoxCox, 2, var) < 1e+10))
(apply(lmerRanefMatControl2017BoxCox, 2, var) < 1e+10) &
table((apply(lmerRanefMatControl2017BoxCox, 2, var) < 1e+10) &
        (apply(lmerRanefMatDrought2017BoxCox, 2, var) < 1e+10))
table(which((apply(lmerRanefMatControl2017BoxCox, 2, var) < 1e+10) &
        (apply(lmerRanefMatDrought2017BoxCox, 2, var) < 1e+10)))


### ratio of metabolites ( heritability > 0.5 )
hist(lmerVarHerits2017BoxCox[, 4])
hist(lmerVarHeritsControl2017BoxCox[, 4])
hist(lmerVarHeritsDrought2017BoxCox[, 4])

nrow(lmerVarHerits2017BoxCox[lmerVarHerits2017BoxCox[, 4] > 0.5, ] ) / nrow(lmerVarHerits2017BoxCox)
nrow(lmerVarHeritsControl2017BoxCox[lmerVarHeritsControl2017BoxCox[, 4] > 0.5, ] ) / nrow(lmerVarHeritsControl2017BoxCox)
nrow(lmerVarHeritsDrought2017BoxCox[lmerVarHeritsDrought2017BoxCox[, 4] > 0.5, ] ) / nrow(lmerVarHeritsDrought2017BoxCox)


#### to edit ####
# #### 2.3.4.  Extract flavonoid > 0.9 line heritability and data, and write csv files ####
# metab2017FlavonoidRaw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
# See(metab2017FlavonoidRaw, rown = 6, coln = 12)
# metab2017FlavonoidNames <- colnames(metab2017FlavonoidRaw[, 10:ncol(metab2017FlavonoidRaw)])
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017[metab2017FlavonoidNames, ]
# See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017FlavonoidMoreThan0.9Heritability[lmerVarHerits2017FlavonoidMoreThan0.9Heritability$heritLine > 0.9, ]
# See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
#
# varietyNames <- rownames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017FlavonoidMoreThan0.9Heritability$heritLine
# names(lmerVarHerits2017FlavonoidMoreThan0.9Heritability) <- varietyNames
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# colnames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability) <- "heritLine"
# See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
#
#
# onlyMetab2017Flavonoid <- metab2017FlavonoidRaw[, 10:ncol(metab2017FlavonoidRaw)]
# onlyMetabFlavonoidMoreThan0.9Heritability <- onlyMetab2017Flavonoid[, rownames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)]
# See(onlyMetabFlavonoidMoreThan0.9Heritability)
# metab2017FlavonoidMoreThan0.9HeritabilityRaw <- cbind(metab2017FlavonoidRaw[, 1:9], onlyMetabFlavonoidMoreThan0.9Heritability)
# See(metab2017FlavonoidMoreThan0.9HeritabilityRaw, coln = 12)
#
#
#
# fileNameFlavonoidAnnotationMoreThan0.9Heritability2017 <- paste0(dirMidSTAMBSH, scriptID, "_Flavonoid_>0.9_heritability_2017.csv")
# fileNameFlavonoidMoreThan0.9Heritability2017Raw <- paste0("data/phenotype/", "2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
#
#
# write.csv(x = lmerVarHerits2017FlavonoidMoreThan0.9Heritability, file = fileNameFlavonoidAnnotationMoreThan0.9Heritability2017)
# write.csv(x = metab2017FlavonoidMoreThan0.9HeritabilityRaw, file = fileNameFlavonoidMoreThan0.9Heritability2017Raw)
#
#
#
##### 4.4.  Histogram of Broad sense heritability in 2017 #####
#### 4.4.1 Total ####
See(lmerVarHerits2017BoxCox)
class(lmerVarHerits2017BoxCox)

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_individual_for_Total_BoxCox.pdf"))
hist(lmerVarHerits2017BoxCox[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

class(lmerVarHerits2017BoxCox)
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_line_for_Total_BoxCox.pdf"))
hist(lmerVarHerits2017BoxCox[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()



#### 4.4.2 Control ####
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_individual_for_Control_BoxCox.pdf"))
hist(lmerVarHeritsControl2017BoxCox[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_line_for_Control_BoxCox.pdf"))
hist(lmerVarHeritsControl2017BoxCox[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()




#### 4.4.3 Drought ####
pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_individual_for_Drought_BoxCox.pdf"))
hist(lmerVarHeritsDrought2017BoxCox[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_Heritability_of_line_for_Drought_BoxCox.pdf"))
hist(lmerVarHeritsDrought2017BoxCox[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()




# #### 2.4.4 Flavonoid ####
# flavonoidMetab <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
# See(flavonoidMetab)
# flavonoidMetabNames <- flavonoidMetab[, 1]
#
# lmerVarHerits2017DF <- as.data.frame((lmerVarHerits2017))
# lmerVarHerits2017FlavonoidDF <- lmerVarHerits2017DF[flavonoidMetabNames, ]
# See(lmerVarHerits2017FlavonoidDF)
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_for_Total.pdf"))
# hist(lmerVarHerits2017FlavonoidDF$heritLine, ylim = c(0,50), xlab =  "_Heritability of Line", main = "Flavonoid")
# dev.off()
#
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_for_Total.pdf"))
# hist(lmerVarHerits2017FlavonoidDF$heritInd, ylim = c(0,70), xlab = "Heritability of Individual", main = "Flavonoid")
# dev.off()
#
# # metabHeritsLineHigher2017 <- lmerVarHerits2017DF[lmerVarHerits2017DF$heritLine > 0.9, ]
#
#
#
# #### 2.4.5 Not Flavonoid ####
# flavonoidMetab <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
# See(flavonoidMetab)
# flavonoidMetabNames <- flavonoidMetab[, 1]
#
# lmerVarHerits2017DF <- as.data.frame((lmerVarHerits2017))
# lmerVarHerits2017NonFlavonoidDF <- lmerVarHerits2017DF[(!rownames(lmerVarHerits2017DF) %in% flavonoidMetabNames), ]
# See(lmerVarHerits2017NonFlavonoidDF)
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_non-related_for_Total.pdf"))
# hist(lmerVarHerits2017NonFlavonoidDF$heritLine, xlim = c(0,1), ylim = c(0,50), xlab = "Heritability of Line", main = "Flavonoid non-related")
# dev.off()
#
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_non-related_for_Total.pdf"))
# hist(lmerVarHerits2017NonFlavonoidDF$heritInd, xlim = c(0,1), ylim = c(0,70), breaks = seq(0,1,0.1), xlab = "Heritability of Individual", main = "Flavonoid non-related")
# dev.off()
#
#
#
#
# # col.has.na <- apply(lmerRanefMatDrought2017, 2, function(x){any(is.na(x))})
# # which(col.has.na)
# # table(col.has.na)
# #
# # lmerRanefMat2017NACol <- lmerRanefMat2017[, col.has.na]
# # lmerRanefMat2017NACol
# #
# # metab2017Raw["UA-4805", ]
# #
# #
# # row.has.na <- apply(lmerRanefMat2017, 1, function(x){any(is.na(x))})
# # which(row.has.na)
#
#
#
#
# ##### 2.5. Histogram of BLUP for each variety of metabolites in 2017 #####
# gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
# See(gvMetab2017Total)
# table(is.na(gvMetab2017Total))
# gvMetab2017Total <- na.omit(gvMetab2017Total)
# See(gvMetab2017Total)
#
# metabNames <- colnames(gvMetab2017Total)
# nMetab <- ncol(gvMetab2017Total)
#
# ### to edit
# # metab2017Control <- metab2017Raw[metab2017Raw$block == "C", ]
# # See(metab2017Control, coln = 10)
# # metab2017ControlOnlyMetab <- metab2017Control[, metabStart:metabEnd]
# # See(metab2017ControlOnlyMetab)
# #
# # metab2017Drought <- metab2017Raw[metab2017Raw$block == "D", ]
# # metab2017DroughtOnlyMetab <- metab2017Drought[, metabStart:metabEnd]
# ###
#
# metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
# metabNamesFlavonoid <- metabFlavonoid[, "Name"]
# metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
# See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
# metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])
#
#
# ### Total
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Histogram_of_BLUP_for_each_metabolite_in_2017/"))
# dirMidSTAMBSHBLUPHist2017 <- paste0(dirMidSTAMBSH, scriptID, "_Histogram_of_BLUP_for_each_metabolite_in_2017/")
#
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/"))
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- gvMetab2017Total[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   }
# }
#
#
# for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
#   metabNow <- gvMetab2017Total[, metabName]
#
#   pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
#   dev.off()
#
# }



### to edit
# ### Control
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017ControlOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
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
#   pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
#
#
#
# ### Drought
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017DroughtOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
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
#   pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
###










###### 5. Estimate broad-sense heritability for Metabolomic data in 2018, using Box-Cox data ######
##### 5.1. Read Metabolomic Box-Cox data in 2018 into R #####
metab2018BoxCox <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_BoxCox.csv", row.names = 1)
See(metab2018BoxCox, rown = 6, coln = 12)


metabStart <- 10
metabEnd <- ncol(metab2018BoxCox)
metabNames <- colnames(metab2018BoxCox)[metabStart:metabEnd]
varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2018BoxCox$variety)
blockNames <- unique(metab2018BoxCox$block)

nMetab <- metabEnd - metabStart + 1
nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2018BoxCox$ind))
nBlock <- length(blockNames)



##### 5.2. Perform lmer for X00006 in Metabolomic data in 2018 #####
#### 5.2.1 Perform lmer for X00006 in all Metabolomic data in 2018 ####
# lmerRes2018 <- lmer(formula = X00006 ? (1 | variety) + block,
#                     data = metab2018Raw)
# summary(lmerRes2018)
#
# lmerVars2018 <- data.frame(VarCorr(lmerRes2018))$vcov
# lmerVu2018 <- lmerVars2018[1]
# lmerVe2018 <- lmerVars2018[2]
#
# lmerHeritInd2018 <- lmerVu2018 / (lmerVu2018 + lmerVe2018)
# lmerHeritLine2018 <- lmerVu2018 / (lmerVu2018 + lmerVe2018 / (nRep * nBlock))
#
# lmerRanef2018 <- ranef(lmerRes2018)$variety[varietyNames,]
#
#
#
# #### 5.2.2 Perform lmer for X00006 in Metabolomic data for control in 2018 ####
# metabControl2018 <- metab2018Raw[metab2018Raw$block == "C", ]
# lmerControlRes2018 <- lmer(formula = X00006 ? (1 | variety),
#                            data = metabControl2018)
#
# lmerControlVars2018 <- data.frame(VarCorr(lmerControlRes2018))$vcov
# lmerControlVu2018 <- lmerControlVars2018[1]
# lmerControlVe2018 <- lmerControlVars2018[2]
#
# lmerControlHeritInd2018 <- lmerControlVu2018 / (lmerControlVu2018 + lmerControlVe2018)
# lmerControlHeritLine2018 <- lmerControlVu2018 / (lmerControlVu2018 + lmerControlVe2018 / nRep)
#
# lmerControlRanef2018 <- ranef(lmerControlRes2018)$variety[varietyNames, ]
#
#
# #### 5.2.3  Perform lmer for X00006 in Metabolomic data for drought in 2018 ####
# metabDrought2018 <- metab2018Raw[metab2018Raw$block == "D", ]
# lmerDroughtRes2018 <- lmer(formula = X00006 ? (1 | variety),
#                            data = metabDrought2018)
#
# lmerDroughtVars2018 <- data.frame(VarCorr(lmerDroughtRes2018))$vcov
# lmerDroughtVu2018 <- lmerDroughtVars2018[1]
# lmerDroughtVe2018 <- lmerDroughtVars2018[2]
#
# lmerDroughtHeritInd2018 <- lmerDroughtVu2018 / (lmerDroughtVu2018 + lmerDroughtVe2018)
# lmerDroughtHeritLine2018 <- lmerDroughtVu2018 / (lmerDroughtVu2018 + lmerDroughtVe2018 / nRep)
#
# lmerDroughtRanef2018 <- ranef(lmerDroughtRes2018)$variety[varietyNames, ]





##### 5.3. Perform lmer for all Metabolomic data in 2018 #####
#### 5.3.1 Perform lmer for all Metabolomic data in 2018 ####
lmerVarHerits2018BoxCox <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMat2018BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)

### extra
lmerInterceptPlusRanefMat2018BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)
###



rownames(lmerVarHerits2018BoxCox) <- metabNames
colnames(lmerVarHerits2018BoxCox) <- varNames
rownames(lmerRanefMat2018BoxCox) <- varietyNames
colnames(lmerRanefMat2018BoxCox) <- metabNames


### extra
rownames(lmerInterceptPlusRanefMat2018BoxCox) <- varietyNames
colnames(lmerInterceptPlusRanefMat2018BoxCox) <- metabNames
###


for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2018BoxCox[, metabNo]

  lmerRes2018BoxCox <- lmer(formula = metabNow ? (1 | variety) + block,
                            data = metab2018BoxCox)

  lmerVars2018BoxCox <- data.frame(VarCorr(lmerRes2018BoxCox))$vcov
  lmerVu2018BoxCox <- lmerVars2018BoxCox[1]
  lmerVe2018BoxCox <- lmerVars2018BoxCox[2]


  lmerHeritInd2018BoxCox <- lmerVu2018BoxCox / (lmerVu2018BoxCox + lmerVe2018BoxCox)
  lmerHeritLine2018BoxCox <- lmerVu2018BoxCox / (lmerVu2018BoxCox + lmerVe2018BoxCox / (nRep * nBlock))

  lmerRanef2018BoxCox <- ranef(lmerRes2018BoxCox)$variety[varietyNames, ]


  lmerVarHerits2018BoxCox[metabNo - metabStart + 1, ] <-
    c(lmerVars2018BoxCox, lmerHeritInd2018BoxCox, lmerHeritLine2018BoxCox)
  lmerRanefMat2018BoxCox[, metabNo - metabStart + 1] <- lmerRanef2018BoxCox

  ###
  lmerInterceptPlusRanefMat2018BoxCox[, metabNo - metabStart + 1] <- coef(lmerRes2018BoxCox)$variety[varietyNames, 1]
  ###

}


###
See(lmerInterceptPlusRanefMat2018BoxCox, coln = 188)
table(lmerInterceptPlusRanefMat2018BoxCox > 0)
table(apply(lmerInterceptPlusRanefMat2018BoxCox, 2, function(x)any(x < 0)))
apply(lmerInterceptPlusRanefMat2018BoxCox, 2, min)

metabNamesMinMinus <- metabNames[apply(lmerInterceptPlusRanefMat2018BoxCox, 2, function(x)any(x < 0))]
metabNamesMinMinus <- na.omit(metabNamesMinMinus)
lmerInterceptPlusRanefMat2018BoxCox[, metabNamesMinMinus]
###




#### 5.3.2 Perform lmer for all Metabolomic data for control in 2018 ####
metabControl2018BoxCox <- metab2018BoxCox[metab2018BoxCox$block == "C", ]

lmerVarHeritsControl2018BoxCox <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatControl2018BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsControl2018BoxCox) <- metabNames
colnames(lmerVarHeritsControl2018BoxCox) <- varNames
rownames(lmerRanefMatControl2018BoxCox) <- varietyNames
colnames(lmerRanefMatControl2018BoxCox) <- metabNames



for (metabNo in metabStart:metabEnd){
  metabControlNow <- metabControl2018BoxCox[, metabNo]

  lmerControlRes2018BoxCox <- lmer(formula = metabControlNow ? (1 | variety),
                                   data = metabControl2018BoxCox)

  lmerControlVars2018BoxCox <- data.frame(VarCorr(lmerControlRes2018BoxCox))$vcov
  lmerControlVu2018BoxCox <- lmerControlVars2018BoxCox[1]
  lmerControlVe2018BoxCox <- lmerControlVars2018BoxCox[2]

  lmerControlHeritInd2018BoxCox <- lmerControlVu2018BoxCox / (lmerControlVu2018BoxCox + lmerControlVe2018BoxCox)
  lmerControlHeritLine2018BoxCox <- lmerControlVu2018BoxCox / (lmerControlVu2018BoxCox + lmerControlVe2018BoxCox / nRep)

  lmerControlRanef2018BoxCox <- ranef(lmerControlRes2018BoxCox)$variety[varietyNames, ]


  lmerVarHeritsControl2018BoxCox[metabNo - metabStart + 1, ] <-
    c(lmerControlVars2018BoxCox, lmerControlHeritInd2018BoxCox, lmerControlHeritLine2018BoxCox)
  lmerRanefMatControl2018BoxCox[, metabNo - metabStart + 1] <- lmerControlRanef2018BoxCox
}



#### 5.3.3  Perform lmer for all Metabolomic data for drought in 2018 ####
metabDrought2018BoxCox <- metab2018BoxCox[metab2018BoxCox$block == "D", ]

lmerVarHeritsDrought2018BoxCox <- matrix(NA, nrow = nMetab, ncol = nVar)
lmerRanefMatDrought2018BoxCox <- matrix(NA, nrow = nVariety, ncol = nMetab)

rownames(lmerVarHeritsDrought2018BoxCox) <- metabNames
colnames(lmerVarHeritsDrought2018BoxCox) <- varNames
rownames(lmerRanefMatDrought2018BoxCox) <- varietyNames
colnames(lmerRanefMatDrought2018BoxCox) <- metabNames



for (metabNo in metabStart:metabEnd){
  metabDroughtNow <- metabDrought2018BoxCox[, metabNo]

  lmerDroughtRes2018BoxCox <- lmer(formula = metabDroughtNow ? (1 | variety),
                                   data = metabDrought2018BoxCox)

  lmerDroughtVars2018BoxCox <- data.frame(VarCorr(lmerDroughtRes2018BoxCox))$vcov
  lmerDroughtVu2018BoxCox <- lmerDroughtVars2018BoxCox[1]
  lmerDroughtVe2018BoxCox <- lmerDroughtVars2018BoxCox[2]

  lmerDroughtHeritInd2018BoxCox <- lmerDroughtVu2018BoxCox / (lmerDroughtVu2018BoxCox + lmerDroughtVe2018BoxCox)
  lmerDroughtHeritLine2018BoxCox <- lmerDroughtVu2018BoxCox / (lmerDroughtVu2018BoxCox + lmerDroughtVe2018BoxCox / nRep)

  lmerDroughtRanef2018BoxCox <- ranef(lmerDroughtRes2018BoxCox)$variety[varietyNames, ]


  lmerVarHeritsDrought2018BoxCox[metabNo - metabStart + 1, ] <-
    c(lmerDroughtVars2018BoxCox, lmerDroughtHeritInd2018BoxCox, lmerDroughtHeritLine2018BoxCox)
  lmerRanefMatDrought2018BoxCox[, metabNo - metabStart + 1] <- lmerDroughtRanef2018BoxCox
}


###
lmerRanefMatCPlusD2018BoxCox <- lmerRanefMatControl2018 + lmerRanefMatDrought2018
lmerRanefMatCMinusD2018BoxCox <- lmerRanefMatControl2018 - lmerRanefMatDrought2018



### write csv files for each treatment x herit, ranef
fileNameTotal2018BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_2018_BoxCox.csv")
fileNameControl2018BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Control_2018_BoxCox.csv")
fileNameDrought2018BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Drought_2018_BoxCox.csv")
fileNameCPlusD2018BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CPlusD_2018_BoxCox.csv")
fileNameCMinusD2018BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_CMinusD_2018_BoxCox.csv")

fileNameTotalInterceptPlusRanef2018BoxCox <- paste0(dirMidSTAMBSH, scriptID, "_lmer_genotypic_values_Total_intercept_plus_ranef_2018_BoxCox.csv")


write.csv(x = lmerRanefMat2018BoxCox,
          file = fileNameTotal2018BoxCox)
write.csv(x = lmerRanefMatControl2018BoxCox,
          file = fileNameControl2018BoxCox)
write.csv(x = lmerRanefMatDrought2018BoxCox,
          file = fileNameDrought2018BoxCox)
write.csv(x = lmerRanefMatCPlusD2018BoxCox,
          file = fileNameCPlusD2018BoxCox)
write.csv(x = lmerRanefMatCMinusD2018BoxCox,
          file = fileNameCMinusD2018BoxCox)

write.csv(x = lmerInterceptPlusRanefMat2018BoxCox,
          file = fileNameTotalInterceptPlusRanef2018BoxCox)




plot(round(apply(lmerRanefMatControl2018BoxCox, 2, var), 5),
     round(apply(lmerRanefMatDrought2018BoxCox, 2, var), 5))
plot(round(apply(lmerRanefMatControl2018BoxCox, 2, mean), 5),
     round(apply(lmerRanefMatDrought2018BoxCox, 2, mean), 5))
which((apply(lmerRanefMatControl2018BoxCox, 2, var) < 1e+10) &
        (apply(lmerRanefMatDrought2018BoxCox, 2, var) < 1e+10))
(apply(lmerRanefMatControl2018BoxCox, 2, var) < 1e+10) &
  table((apply(lmerRanefMatControl2018BoxCox, 2, var) < 1e+10) &
          (apply(lmerRanefMatDrought2018BoxCox, 2, var) < 1e+10))
table(which((apply(lmerRanefMatControl2018BoxCox, 2, var) < 1e+10) &
              (apply(lmerRanefMatDrought2018BoxCox, 2, var) < 1e+10)))


### ratio of metabolites ( heritability > 0.5 )
hist(lmerVarHerits2018BoxCox[, 4])
hist(lmerVarHeritsControl2018BoxCox[, 4])
hist(lmerVarHeritsDrought2018BoxCox[, 4])

nrow(lmerVarHerits2018BoxCox[lmerVarHerits2018BoxCox[, 4] > 0.5, ] ) / nrow(lmerVarHerits2018BoxCox)
nrow(lmerVarHeritsControl2018BoxCox[lmerVarHeritsControl2018BoxCox[, 4] > 0.5, ] ) / nrow(lmerVarHeritsControl2018BoxCox)
nrow(lmerVarHeritsDrought2018BoxCox[lmerVarHeritsDrought2018BoxCox[, 4] > 0.5, ] ) / nrow(lmerVarHeritsDrought2018BoxCox)


#### to edit ####
# #### 2.3.4.  Extract flavonoid > 0.9 line heritability and data, and write csv files ####
# metab2017FlavonoidRaw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
# See(metab2017FlavonoidRaw, rown = 6, coln = 12)
# metab2017FlavonoidNames <- colnames(metab2017FlavonoidRaw[, 10:ncol(metab2017FlavonoidRaw)])
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017[metab2017FlavonoidNames, ]
# See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017FlavonoidMoreThan0.9Heritability[lmerVarHerits2017FlavonoidMoreThan0.9Heritability$heritLine > 0.9, ]
# See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
#
# varietyNames <- rownames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- lmerVarHerits2017FlavonoidMoreThan0.9Heritability$heritLine
# names(lmerVarHerits2017FlavonoidMoreThan0.9Heritability) <- varietyNames
# lmerVarHerits2017FlavonoidMoreThan0.9Heritability <- as.data.frame(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
# colnames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability) <- "heritLine"
# See(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)
#
#
# onlyMetab2017Flavonoid <- metab2017FlavonoidRaw[, 10:ncol(metab2017FlavonoidRaw)]
# onlyMetabFlavonoidMoreThan0.9Heritability <- onlyMetab2017Flavonoid[, rownames(lmerVarHerits2017FlavonoidMoreThan0.9Heritability)]
# See(onlyMetabFlavonoidMoreThan0.9Heritability)
# metab2017FlavonoidMoreThan0.9HeritabilityRaw <- cbind(metab2017FlavonoidRaw[, 1:9], onlyMetabFlavonoidMoreThan0.9Heritability)
# See(metab2017FlavonoidMoreThan0.9HeritabilityRaw, coln = 12)
#
#
#
# fileNameFlavonoidAnnotationMoreThan0.9Heritability2017 <- paste0(dirMidSTAMBSH, scriptID, "_Flavonoid_>0.9_heritability_2017.csv")
# fileNameFlavonoidMoreThan0.9Heritability2017Raw <- paste0("data/phenotype/", "2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
#
#
# write.csv(x = lmerVarHerits2017FlavonoidMoreThan0.9Heritability, file = fileNameFlavonoidAnnotationMoreThan0.9Heritability2017)
# write.csv(x = metab2017FlavonoidMoreThan0.9HeritabilityRaw, file = fileNameFlavonoidMoreThan0.9Heritability2017Raw)
#
#
#
##### 5.4.  Histogram of Broad sense heritability in 2018 #####
#### 5.4.1 Total ####
See(lmerVarHerits2018BoxCox)
class(lmerVarHerits2018BoxCox)

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_individual_for_Total_BoxCox.pdf"))
hist(lmerVarHerits2018BoxCox[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

class(lmerVarHerits2018BoxCox)
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_line_for_Total_BoxCox.pdf"))
hist(lmerVarHerits2018BoxCox[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()



#### 5.4.2 Control ####
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_individual_for_Control_BoxCox.pdf"))
hist(lmerVarHeritsControl2018BoxCox[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_line_for_Control_BoxCox.pdf"))
hist(lmerVarHeritsControl2018BoxCox[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()




#### 5.4.3 Drought ####
pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_individual_for_Drought_BoxCox.pdf"))
hist(lmerVarHeritsDrought2018BoxCox[, "heritInd"], ylim = c(0,70),xlab = "Heritability of Individual")
dev.off()

pdf(paste0(dirMidSTAMBSHHeritHist2018, scriptID, "_Histogram_of_Heritability_of_line_for_Drought_BoxCox.pdf"))
hist(lmerVarHeritsDrought2018BoxCox[, "heritLine"], ylim = c(0,50), xlab = "Heritability of Line")
dev.off()




# #### 2.4.4 Flavonoid ####
# flavonoidMetab <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
# See(flavonoidMetab)
# flavonoidMetabNames <- flavonoidMetab[, 1]
#
# lmerVarHerits2017DF <- as.data.frame((lmerVarHerits2017))
# lmerVarHerits2017FlavonoidDF <- lmerVarHerits2017DF[flavonoidMetabNames, ]
# See(lmerVarHerits2017FlavonoidDF)
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_for_Total.pdf"))
# hist(lmerVarHerits2017FlavonoidDF$heritLine, ylim = c(0,50), xlab =  "_Heritability of Line", main = "Flavonoid")
# dev.off()
#
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_for_Total.pdf"))
# hist(lmerVarHerits2017FlavonoidDF$heritInd, ylim = c(0,70), xlab = "Heritability of Individual", main = "Flavonoid")
# dev.off()
#
# # metabHeritsLineHigher2017 <- lmerVarHerits2017DF[lmerVarHerits2017DF$heritLine > 0.9, ]
#
#
#
# #### 2.4.5 Not Flavonoid ####
# flavonoidMetab <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
# See(flavonoidMetab)
# flavonoidMetabNames <- flavonoidMetab[, 1]
#
# lmerVarHerits2017DF <- as.data.frame((lmerVarHerits2017))
# lmerVarHerits2017NonFlavonoidDF <- lmerVarHerits2017DF[(!rownames(lmerVarHerits2017DF) %in% flavonoidMetabNames), ]
# See(lmerVarHerits2017NonFlavonoidDF)
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_line_for_flavonoid_non-related_for_Total.pdf"))
# hist(lmerVarHerits2017NonFlavonoidDF$heritLine, xlim = c(0,1), ylim = c(0,50), xlab = "Heritability of Line", main = "Flavonoid non-related")
# dev.off()
#
# pdf(paste0(dirMidSTAMBSHHeritHist2017, scriptID, "_Histogram_of_heritability_of_individual_for_flavonoid_non-related_for_Total.pdf"))
# hist(lmerVarHerits2017NonFlavonoidDF$heritInd, xlim = c(0,1), ylim = c(0,70), breaks = seq(0,1,0.1), xlab = "Heritability of Individual", main = "Flavonoid non-related")
# dev.off()
#
#
#
#
# # col.has.na <- apply(lmerRanefMatDrought2017, 2, function(x){any(is.na(x))})
# # which(col.has.na)
# # table(col.has.na)
# #
# # lmerRanefMat2017NACol <- lmerRanefMat2017[, col.has.na]
# # lmerRanefMat2017NACol
# #
# # metab2017Raw["UA-4805", ]
# #
# #
# # row.has.na <- apply(lmerRanefMat2017, 1, function(x){any(is.na(x))})
# # which(row.has.na)
#
#
#
#
# ##### 2.5. Histogram of BLUP for each variety of metabolites in 2017 #####
# gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
# See(gvMetab2017Total)
# table(is.na(gvMetab2017Total))
# gvMetab2017Total <- na.omit(gvMetab2017Total)
# See(gvMetab2017Total)
#
# metabNames <- colnames(gvMetab2017Total)
# nMetab <- ncol(gvMetab2017Total)
#
# ### to edit
# # metab2017Control <- metab2017Raw[metab2017Raw$block == "C", ]
# # See(metab2017Control, coln = 10)
# # metab2017ControlOnlyMetab <- metab2017Control[, metabStart:metabEnd]
# # See(metab2017ControlOnlyMetab)
# #
# # metab2017Drought <- metab2017Raw[metab2017Raw$block == "D", ]
# # metab2017DroughtOnlyMetab <- metab2017Drought[, metabStart:metabEnd]
# ###
#
# metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
# metabNamesFlavonoid <- metabFlavonoid[, "Name"]
# metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
# See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
# metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])
#
#
# ### Total
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Histogram_of_BLUP_for_each_metabolite_in_2017/"))
# dirMidSTAMBSHBLUPHist2017 <- paste0(dirMidSTAMBSH, scriptID, "_Histogram_of_BLUP_for_each_metabolite_in_2017/")
#
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/"))
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- gvMetab2017Total[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   }
# }
#
#
# for(metabName in metabNamesFlavonoidHeritabilityMoreThan0.9){
#   metabNow <- gvMetab2017Total[, metabName]
#
#   pdf(paste0(dirMidSTAMBSHBLUPHist2017, scriptID, "_Total/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 50, col = "blue")
#   dev.off()
#
# }



### to edit
# ### Control
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017ControlOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
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
#   pdf(paste0(dirMidSTAMBSH, scriptID, "_Control/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
#
#
#
# ### Drought
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Non_Flavonoid/"))
# dir.create(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/"))
#
#
# metabNo <- 1
# for(metabNo in 1:nMetab){
#   metabNow <- metab2017DroughtOnlyMetab[, metabNo]
#   metabName <- metabNames[metabNo]
#
#   if (metabName %in% metabNamesFlavonoid){
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid/", scriptID, "_", metabName, ".pdf"))
#     hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#     # hist(metabNow, freq = FALSE, probability = TRUE)
#     dev.off()
#
#   } else {
#     pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Non_Flavonoid/", scriptID, "_", metabName, ".pdf"))
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
#   pdf(paste0(dirMidSTAMBSH, scriptID, "_Drought/", scriptID, "_Flavonoid_heritability_>0.9/", scriptID, "_", metabName, ".pdf"))
#   hist(metabNow, xlim = c(min(metabNow), max(metabNow)), breaks = 200, col = "blue")
#   dev.off()
#
# }
###













#########################
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

# table(is.na(metab2017Raw))
# table(is.na(metab2017Raw[, 10:ncol(metab2017Raw)]))
# See(metab2017Raw, rown = 1, coln = 210)
# See(metab2017Raw, coln = 13)
# apply(metab2017Raw, 1, function(x)any(is.na(x)))
# table(apply(metab2017Raw[, 10:ncol(metab2017Raw)], 1, function(x)any(is.na(x))))
# apply(metab2017Raw, 2, function(x)any(is.na(x)))
# table(apply(metab2017Raw, 2, function(x)any(is.na(x))))




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







