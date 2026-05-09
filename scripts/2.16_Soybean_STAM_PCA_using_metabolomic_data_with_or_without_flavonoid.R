##########################################################################################
######  Title: 2.16_Soybean_STAM_PCA_using_metabolomic_data_with_or_without_flavonoid######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/04 (Created), 2024/04/04 (Last Updated)                       ######
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

scriptID <- "2.16"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
nPC <- 6


dirMidSTAMPCA <- paste0(dirMidSTAMBase, scriptID,
                        "_PCA_with_or_without_flavonoid/")
dir.create(dirMidSTAMPCA)
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







###### 2. Perform PCA for data of metabolites only related to Flavonoid pathway in 2017 ######
##### 2.1. Read data of metabolites only related to Flavonoid pathway in 2017 into R #####
metab2017FlavonoidRaw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
See(metab2017FlavonoidRaw, rown = 6, coln = 12)

metabStart <- 10
metabEnd <- ncol(metab2017FlavonoidRaw)
metabNames <- colnames(metab2017FlavonoidRaw)[metabStart:metabEnd]
# varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2017FlavonoidRaw$variety)
blockNames <- unique(metab2017FlavonoidRaw$block)

nMetab <- metabEnd - metabStart + 1
# nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2017FlavonoidRaw$ind))
nBlock <- length(blockNames)


##### 2.2. Perform PCA in Flavonoid Metabolomic data in 2017 #####
metab2017FlavonoidScaled <- scale(metab2017FlavonoidRaw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotalFlavonoid2017 <- pcaMethods::pca(metab2017FlavonoidScaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2TotalFlavonoid2017 <- pcaResTotalFlavonoid2017@R2
barplot(pcaR2TotalFlavonoid2017)
pcaR2CumTotalFlavonoid2017 <- pcaResTotalFlavonoid2017@R2cum
pcaScoresTotalFlavonoid2017 <- pcaResTotalFlavonoid2017@scores

head(pcaScoresTotalFlavonoid2017)
table(is.na(pcaScoresTotalFlavonoid2017))

metabPCAFlavonoid2017 <- data.frame(metab2017FlavonoidRaw[, 1:(metabStart - 1)],
                           pcaScoresTotalFlavonoid2017)
head(metabPCAFlavonoid2017)
See(metabPCAFlavonoid2017)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_related_to_flavonoid_pathway_2017.csv")
write.csv(x = metabPCAFlavonoid2017, file = fileName)




###### 3. Perform PCA for data in 2017 Without metabolites related to Flavonoid pathway  ######
##### 3.1. Read data in 2017 and Extract metabolites Not related to Flavonoid pathway  into R #####
metab2017Raw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
metab2017FlavonoidRaw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
See(metab2017Raw, rown = 6, coln = 12)
See(metab2017FlavonoidRaw, rown = 6, coln = 12)

OnlyMetab2017 <- metab2017Raw[, 10:ncol(metab2017Raw)]
OnlyMetab2017Flavonoid <- metab2017FlavonoidRaw[, metabStart:metabEnd]
metab2017WithoutFlavonoid <- OnlyMetab2017[, !(colnames(OnlyMetab2017) %in% colnames(OnlyMetab2017Flavonoid))]
head(metab2017WithoutFlavonoid)
See(metab2017WithoutFlavonoid)

# metabStart <- 10
# metabEnd <- ncol(metab2017FlavonoidRaw)
# metabNames <- colnames(metab2017FlavonoidRaw)[metabStart:metabEnd]
# # varNames <- c("Vu", "Ve", "heritInd", "heritLine")
# varietyNames <- unique(metab2017FlavonoidRaw$variety)
# blockNames <- unique(metab2017FlavonoidRaw$block)
#
# nMetab <- metabEnd - metabStart + 1
# # nVar <- length(varNames)
# nVariety <- length(varietyNames)
# nRep <- length(unique(metab2017FlavonoidRaw$ind))
# nBlock <- length(blockNames)
#
#


##### 3.2. Perform PCA Without Flavonoid Metabolomic data in 2017 #####
metab2017WithoutFlavonoidScaled <- scale(metab2017WithoutFlavonoid, center = TRUE, scale = TRUE)
pcaResTotalWithoutFlavonoid2017 <- pcaMethods::pca(metab2017WithoutFlavonoidScaled,
                                            method = "ppca",
                                            nPcs = nPC,
                                            completeObs = FALSE)
pcaR2TotalWithoutFlavonoid2017 <- pcaResTotalWithoutFlavonoid2017@R2
barplot(pcaR2TotalWithoutFlavonoid2017)
pcaR2CumTotalWithoutFlavonoid2017 <- pcaResTotalWithoutFlavonoid2017@R2cum
pcaScoresTotalWithoutFlavonoid2017 <- pcaResTotalWithoutFlavonoid2017@scores

class(pcaResTotalWithoutFlavonoid2017)
head(pcaScoresTotalWithoutFlavonoid2017)
table(is.na(pcaScoresTotalWithoutFlavonoid2017))

metabPCAWithoutFlavonoid2017 <- data.frame(metab2017FlavonoidRaw[, 1:(metabStart - 1)],
                                    pcaScoresTotalWithoutFlavonoid2017)
head(metabPCAFlavonoid2017)
See(metabPCAFlavonoid2017)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_Not_related_to_flavonoid_pathway_2017.csv")

write.csv(x = metabPCAWithoutFlavonoid2017, file = fileName)





###### 4. Perform PCA for data of metabolites only related to Flavonoid pathway in 2018 ######
##### 4.1. Read data of metabolites only related to Flavonoid pathway in 2018 into R #####
metab2018FlavonoidRaw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
See(metab2018FlavonoidRaw, rown = 6, coln = 12)

metabStart <- 10
metabEnd <- ncol(metab2018FlavonoidRaw)
metabNames <- colnames(metab2018FlavonoidRaw)[metabStart:metabEnd]
# varNames <- c("Vu", "Ve", "heritInd", "heritLine")
varietyNames <- unique(metab2018FlavonoidRaw$variety)
blockNames <- unique(metab2018FlavonoidRaw$block)

nMetab <- metabEnd - metabStart + 1
# nVar <- length(varNames)
nVariety <- length(varietyNames)
nRep <- length(unique(metab2018FlavonoidRaw$ind))
nBlock <- length(blockNames)


##### 4.2. Perform PCA in Flavonoid Metabolomic data in 2018 #####
metab2018FlavonoidScaled <- scale(metab2018FlavonoidRaw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotalFlavonoid2018 <- pcaMethods::pca(metab2018FlavonoidScaled,
                                            method = "ppca",
                                            nPcs = nPC,
                                            completeObs = FALSE)
pcaR2TotalFlavonoid2018 <- pcaResTotalFlavonoid2018@R2
barplot(pcaR2TotalFlavonoid2018)
pcaR2CumTotalFlavonoid2018 <- pcaResTotalFlavonoid2018@R2cum
pcaScoresTotalFlavonoid2018 <- pcaResTotalFlavonoid2018@scores

head(pcaScoresTotalFlavonoid2018)
table(is.na(pcaScoresTotalFlavonoid2018))

metabPCAFlavonoid2018 <- data.frame(metab2018FlavonoidRaw[, 1:(metabStart - 1)],
                                    pcaScoresTotalFlavonoid2018)
head(metabPCAFlavonoid2018)
See(metabPCAFlavonoid2018)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_related_to_flavonoid_pathway_2018.csv")
write.csv(x = metabPCAFlavonoid2018, file = fileName)




###### 5. Perform PCA for data in 2018 Without metabolites related to Flavonoid pathway  ######
##### 5.1. Read data in 2018 and Extract metabolites Not related to Flavonoid pathway  into R #####
metab2018Raw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier.csv")
metab2018FlavonoidRaw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
See(metab2018Raw, rown = 6, coln = 12)
See(metab2018FlavonoidRaw, rown = 6, coln = 12)

OnlyMetab2018 <- metab2018Raw[, 10:ncol(metab2018Raw)]
OnlyMetab2018Flavonoid <- metab2018FlavonoidRaw[, metabStart:metabEnd]
metab2018WithoutFlavonoid <- OnlyMetab2018[, !(colnames(OnlyMetab2018) %in% colnames(OnlyMetab2018Flavonoid))]
head(metab2018WithoutFlavonoid)
See(metab2018WithoutFlavonoid)

# metabStart <- 10
# metabEnd <- ncol(metab2018FlavonoidRaw)
# metabNames <- colnames(metab2018FlavonoidRaw)[metabStart:metabEnd]
# # varNames <- c("Vu", "Ve", "heritInd", "heritLine")
# varietyNames <- unique(metab2018FlavonoidRaw$variety)
# blockNames <- unique(metab2018FlavonoidRaw$block)
#
# nMetab <- metabEnd - metabStart + 1
# # nVar <- length(varNames)
# nVariety <- length(varietyNames)
# nRep <- length(unique(metab2018FlavonoidRaw$ind))
# nBlock <- length(blockNames)
#
#


##### 5.2. Perform PCA Without Flavonoid Metabolomic data in 2018 #####
metab2018WithoutFlavonoidScaled <- scale(metab2018WithoutFlavonoid, center = TRUE, scale = TRUE)
pcaResTotalWithoutFlavonoid2018 <- pcaMethods::pca(metab2018WithoutFlavonoidScaled,
                                                   method = "ppca",
                                                   nPcs = nPC,
                                                   completeObs = FALSE)
pcaR2TotalWithoutFlavonoid2018 <- pcaResTotalWithoutFlavonoid2018@R2
barplot(pcaR2TotalWithoutFlavonoid2018)
pcaR2CumTotalWithoutFlavonoid2018 <- pcaResTotalWithoutFlavonoid2018@R2cum
pcaScoresTotalWithoutFlavonoid2018 <- pcaResTotalWithoutFlavonoid2018@scores

class(pcaResTotalWithoutFlavonoid2018)
head(pcaScoresTotalWithoutFlavonoid2018)
table(is.na(pcaScoresTotalWithoutFlavonoid2018))

metabPCAWithoutFlavonoid2018 <- data.frame(metab2018FlavonoidRaw[, 1:(metabStart - 1)],
                                           pcaScoresTotalWithoutFlavonoid2018)
head(metabPCAFlavonoid2018)
See(metabPCAFlavonoid2018)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_Not_related_to_flavonoid_pathway_2018.csv")

write.csv(x = metabPCAWithoutFlavonoid2018, file = fileName)







