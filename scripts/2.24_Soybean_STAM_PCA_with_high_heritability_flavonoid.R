##########################################################################################
######  Title: 2.24_Soybean_STAM_PCA_with_high_heritability_flavonoid ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/23 (Created), 2024/08/17 (Last Updated)                       ######
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

scriptID <- "2.24"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
# nPC <- 6


dirMidSTAMPCAHighHeritabilityFlavonoid <- paste0(dirMidSTAMBase, scriptID,
                        "_PCA_with_high_heritability_flavonoid/")
dir.create(dirMidSTAMPCAHighHeritabilityFlavonoid)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(pcaMethods)
require(plotly)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)







###### 2. Perform PCA for flavonoid > 0.9 heritability in 2017 ######
##### 2.1. Read data of Flavonoid >0.9 heritability in 2017 into R #####
metab2017FlavonoidMoreThan0.9HeritabilityRaw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metab2017FlavonoidMoreThan0.9HeritabilityRaw, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metab2017FlavonoidMoreThan0.9HeritabilityRaw)
metabNames <- colnames(metab2017FlavonoidMoreThan0.9HeritabilityRaw)[metabStart:metabEnd]
# varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2017FlavonoidMoreThan0.9HeritabilityRaw$variety)
blockNames <- unique(metab2017FlavonoidMoreThan0.9HeritabilityRaw$block)

nMetab <- metabEnd - metabStart + 1
# nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2017FlavonoidMoreThan0.9HeritabilityRaw$ind))
nBlock <- length(blockNames)


##### 2.2. Perform PCA for Flavonoid >0.9 heritability data in 2017 #####
#### 2.2.1. nPC = 6 ####
nPC <- 6
PCAMethod <- "pcaMethods"

metab2017FlavonoidMoreThan0.9HeritabilityScaled <- scale(metab2017FlavonoidMoreThan0.9HeritabilityRaw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotalFlavonoidMoreThan0.9Heritability2017 <- pcaMethods::pca(metab2017FlavonoidMoreThan0.9HeritabilityScaled,
                                            method = "ppca",
                                            nPcs = nPC,
                                            completeObs = FALSE)
pcaR2TotalFlavonoidMoreThan0.9Heritability2017 <- pcaResTotalFlavonoidMoreThan0.9Heritability2017@R2
barplot(pcaR2TotalFlavonoidMoreThan0.9Heritability2017)
pcaR2CumTotalFlavonoidMoreThan0.9Heritability2017 <- pcaResTotalFlavonoidMoreThan0.9Heritability2017@R2cum
pcaScoresTotalFlavonoidMoreThan0.9Heritability2017 <- pcaResTotalFlavonoidMoreThan0.9Heritability2017@scores

head(pcaScoresTotalFlavonoidMoreThan0.9Heritability2017)
table(is.na(pcaScoresTotalFlavonoidMoreThan0.9Heritability2017))

metabPCAFlavonoidMoreThan0.9Heritability2017 <- data.frame(metab2017FlavonoidMoreThan0.9HeritabilityRaw[, 1:(metabStart - 1)],
                                    pcaScoresTotalFlavonoidMoreThan0.9Heritability2017)
head(metabPCAFlavonoidMoreThan0.9Heritability2017)
See(metabPCAFlavonoidMoreThan0.9Heritability2017)

fileName <- paste0(dirMidSTAMPCAHighHeritabilityFlavonoid, scriptID, "_", PCAMethod, "_PCA_", "nPC=", nPC , "_", "flavonoid_metab_>0.9_Heritability_2017.csv")
write.csv(x = metabPCAFlavonoidMoreThan0.9Heritability2017, file = fileName)




#### 2.2.2. nPC = 20 ####
nPC <- 20
PCAMethod <- "pcaMethods"

metab2017FlavonoidMoreThan0.9HeritabilityScaled <- scale(metab2017FlavonoidMoreThan0.9HeritabilityRaw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotalFlavonoidMoreThan0.9Heritability2017 <- pcaMethods::pca(metab2017FlavonoidMoreThan0.9HeritabilityScaled,
                                                                   method = "ppca",
                                                                   nPcs = nPC,
                                                                   completeObs = FALSE)
pcaR2TotalFlavonoidMoreThan0.9Heritability2017 <- pcaResTotalFlavonoidMoreThan0.9Heritability2017@R2
barplot(pcaR2TotalFlavonoidMoreThan0.9Heritability2017)
pcaR2CumTotalFlavonoidMoreThan0.9Heritability2017 <- pcaResTotalFlavonoidMoreThan0.9Heritability2017@R2cum
pcaScoresTotalFlavonoidMoreThan0.9Heritability2017 <- pcaResTotalFlavonoidMoreThan0.9Heritability2017@scores

head(pcaScoresTotalFlavonoidMoreThan0.9Heritability2017)
table(is.na(pcaScoresTotalFlavonoidMoreThan0.9Heritability2017))

metabPCAFlavonoidMoreThan0.9Heritability2017 <- data.frame(metab2017FlavonoidMoreThan0.9HeritabilityRaw[, 1:(metabStart - 1)],
                                                           pcaScoresTotalFlavonoidMoreThan0.9Heritability2017)
head(metabPCAFlavonoidMoreThan0.9Heritability2017)
See(metabPCAFlavonoidMoreThan0.9Heritability2017)

fileName <- paste0(dirMidSTAMPCAHighHeritabilityFlavonoid, scriptID, "_", PCAMethod, "_PCA_", "nPC=", nPC , "_", "flavonoid_metab_>0.9_Heritability_2017.csv")
write.csv(x = metabPCAFlavonoidMoreThan0.9Heritability2017, file = fileName)






###### 3. Perform PCA for flavonoid > 0.9 heritability in 2018 ######
##### 3.1. Read data of Flavonoid >0.9 heritability in 2018 into R #####
metab2018FlavonoidMoreThan0.9HeritabilityRaw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metab2018FlavonoidMoreThan0.9HeritabilityRaw, rown = 6, coln = 12)


metabStart <- 11
metabEnd <- ncol(metab2018FlavonoidMoreThan0.9HeritabilityRaw)
metabNames <- colnames(metab2018FlavonoidMoreThan0.9HeritabilityRaw)[metabStart:metabEnd]
# varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2018FlavonoidMoreThan0.9HeritabilityRaw$variety)
blockNames <- unique(metab2018FlavonoidMoreThan0.9HeritabilityRaw$block)

nMetab <- metabEnd - metabStart + 1
# nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2018FlavonoidMoreThan0.9HeritabilityRaw$ind))
nBlock <- length(blockNames)


##### 3.2. Perform PCA for Flavonoid >0.9 heritability data in 2018 #####
#### 3.2.1. nPC = 6 ####
nPC <- 6
PCAMethod <- "pcaMethods"

metab2018FlavonoidMoreThan0.9HeritabilityScaled <- scale(metab2018FlavonoidMoreThan0.9HeritabilityRaw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotalFlavonoidMoreThan0.9Heritability2018 <- pcaMethods::pca(metab2018FlavonoidMoreThan0.9HeritabilityScaled,
                                                                   method = "ppca",
                                                                   nPcs = nPC,
                                                                   completeObs = FALSE)
pcaR2TotalFlavonoidMoreThan0.9Heritability2018 <- pcaResTotalFlavonoidMoreThan0.9Heritability2018@R2
barplot(pcaR2TotalFlavonoidMoreThan0.9Heritability2018)
pcaR2CumTotalFlavonoidMoreThan0.9Heritability2018 <- pcaResTotalFlavonoidMoreThan0.9Heritability2018@R2cum
pcaScoresTotalFlavonoidMoreThan0.9Heritability2018 <- pcaResTotalFlavonoidMoreThan0.9Heritability2018@scores

head(pcaScoresTotalFlavonoidMoreThan0.9Heritability2018)
table(is.na(pcaScoresTotalFlavonoidMoreThan0.9Heritability2018))

metabPCAFlavonoidMoreThan0.9Heritability2018 <- data.frame(metab2018FlavonoidMoreThan0.9HeritabilityRaw[, 1:(metabStart - 1)],
                                                           pcaScoresTotalFlavonoidMoreThan0.9Heritability2018)
head(metabPCAFlavonoidMoreThan0.9Heritability2018)
See(metabPCAFlavonoidMoreThan0.9Heritability2018)

fileName <- paste0(dirMidSTAMPCAHighHeritabilityFlavonoid, scriptID, "_", PCAMethod, "_PCA_", "nPC=", nPC , "_", "flavonoid_metab_>0.9_Heritability_2018.csv")
write.csv(x = metabPCAFlavonoidMoreThan0.9Heritability2018, file = fileName)




#### 3.2.2. nPC = 20 ####
nPC <- 20
PCAMethod <- "pcaMethods"

metab2018FlavonoidMoreThan0.9HeritabilityScaled <- scale(metab2018FlavonoidMoreThan0.9HeritabilityRaw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotalFlavonoidMoreThan0.9Heritability2018 <- pcaMethods::pca(metab2018FlavonoidMoreThan0.9HeritabilityScaled,
                                                                   method = "ppca",
                                                                   nPcs = nPC,
                                                                   completeObs = FALSE)
pcaR2TotalFlavonoidMoreThan0.9Heritability2018 <- pcaResTotalFlavonoidMoreThan0.9Heritability2018@R2
barplot(pcaR2TotalFlavonoidMoreThan0.9Heritability2018)
pcaR2CumTotalFlavonoidMoreThan0.9Heritability2018 <- pcaResTotalFlavonoidMoreThan0.9Heritability2018@R2cum
pcaScoresTotalFlavonoidMoreThan0.9Heritability2018 <- pcaResTotalFlavonoidMoreThan0.9Heritability2018@scores

head(pcaScoresTotalFlavonoidMoreThan0.9Heritability2018)
table(is.na(pcaScoresTotalFlavonoidMoreThan0.9Heritability2018))

metabPCAFlavonoidMoreThan0.9Heritability2018 <- data.frame(metab2018FlavonoidMoreThan0.9HeritabilityRaw[, 1:(metabStart - 1)],
                                                           pcaScoresTotalFlavonoidMoreThan0.9Heritability2018)
head(metabPCAFlavonoidMoreThan0.9Heritability2018)
See(metabPCAFlavonoidMoreThan0.9Heritability2018)

fileName <- paste0(dirMidSTAMPCAHighHeritabilityFlavonoid, scriptID, "_", PCAMethod, "_PCA_", "nPC=", nPC , "_", "flavonoid_metab_>0.9_Heritability_2018.csv")
write.csv(x = metabPCAFlavonoidMoreThan0.9Heritability2018, file = fileName)


