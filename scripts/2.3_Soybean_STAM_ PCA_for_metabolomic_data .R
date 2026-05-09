##########################################################################################
######  Title: 2.3_Soybean_STAM_PCA_for_metabolomic_data                            ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/13 (Created), 2026/02/19 (Last Updated)                       ######
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

scriptID <- "2.3"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
# nPC <- 6


dirMidSTAMPCA <- paste0(dirMidSTAMBase, scriptID,
                          "_PCA/")
dir.create(dirMidSTAMPCA)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaMethods")

install.packages("pcaMethods")
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(pcaMethods)
require(plotly)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### 2. Perform PCA and adjust p-values for Metabolomic data in 2017 ######
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


##### 2.2. Perform PCA in Metabolomic data in 2017 #####
nPC <- 6

metab2017Scaled <- scale(metab2017Raw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotal2017 <- pcaMethods::pca(metab2017Scaled,
                               method = "ppca",
                               nPcs = nPC,
                               completeObs = FALSE)
pcaR2Total2017 <- pcaResTotal2017@R2
barplot(pcaR2Total2017)
pcaR2CumTotal2017 <- pcaResTotal2017@R2cum
pcaScoresTotal2017 <- pcaResTotal2017@scores

head(pcaScoresTotal2017)
table(is.na(pcaScoresTotal2017))

metabPCA2017 <- data.frame(metab2017Raw[, 1:(metabStart - 1)],
                           pcaScoresTotal2017)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_2017.csv")
write.csv(x = metabPCA2017, file = fileName)




plotly::plot_ly(data = metabPCA2017,
                x = ? PC1,
                y = ? PC2,
                z = ? PC3,
                color = ? block,
                type = "scatter3d",
                mode = "markers",
                hoverinfo = "text",
                text = apply(metabPCA2017, 1, function(l) {
                  paste(names(l), ":", l, collapse = "\n")
                }))

plotly::plot_ly(data = metabPCA2017,
                x = ? PC1,
                y = ? PC2,
                color = ? block,
                type = "scatter",
                mode = "markers",
                hoverinfo = "text",
                text = apply(metabPCA2017, 1, function(l) {
                  paste(names(l), ":", l, collapse = "\n")
                }))

plotly::plot_ly(data = metabPCA2017,
                x = ? PC3,
                y = ? PC4,
                color = ? block,
                type = "scatter",
                mode = "markers",
                hoverinfo = "text",
                text = apply(metabPCA2017, 1, function(l) {
                  paste(names(l), ":", l, collapse = "\n")
                }))



# Figure #
nPC <- 30
pcaResTotal2017 <- pcaMethods::pca(metab2017Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2017 <- pcaResTotal2017@R2
pdf(paste0(dirMidSTAMPCA, scriptID, "_PCA_Explained_variance_of_Metab_2017.pdf"))
barplot(pcaR2Total2017, xlab = "PC", ylab = "R^2")
dev.off()



#### nPC = 100 ####
nPC <- 100

metab2017Scaled <- scale(metab2017Raw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotal2017 <- pcaMethods::pca(metab2017Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2017 <- pcaResTotal2017@R2
barplot(pcaR2Total2017)
pcaR2CumTotal2017 <- pcaResTotal2017@R2cum
pcaScoresTotal2017 <- pcaResTotal2017@scores

head(pcaScoresTotal2017)
table(is.na(pcaScoresTotal2017))

metabPCA2017 <- data.frame(metab2017Raw[, 1:(metabStart - 1)],
                           pcaScoresTotal2017)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_nPC=", nPC, "_2017.csv")
write.csv(x = metabPCA2017, file = fileName)



###### 3. Perform PCA and adjust p-values for Metabolomic data in 2018 ######
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


##### 3.2. Perform PCA in Metabolomic data in 2018 #####
nPC <- 6

metab2018Scaled <- scale(metab2018Raw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotal2018 <- pcaMethods::pca(metab2018Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2018 <- pcaResTotal2018@R2
barplot(pcaR2Total2018)
pcaR2CumTotal2018 <- pcaResTotal2018@R2cum
pcaScoresTotal2018 <- pcaResTotal2018@scores

head(pcaScoresTotal2018)
table(is.na(pcaScoresTotal2018))

metabPCA2018 <- data.frame(metab2018Raw[, 1:(metabStart - 1)],
                           pcaScoresTotal2018)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_metab_2018.csv")
write.csv(x = metabPCA2018, file = fileName)




plotly::plot_ly(data = metabPCA2018,
                x = ? PC1,
                y = ? PC2,
                z = ? PC3,
                color = ? block,
                type = "scatter3d",
                mode = "markers",
                hoverinfo = "text",
                text = apply(metabPCA2018, 1, function(l) {
                  paste(names(l), ":", l, collapse = "\n")
                }))

plotly::plot_ly(data = metabPCA2018,
                x = ? PC1,
                y = ? PC2,
                color = ? block,
                type = "scatter",
                mode = "markers",
                hoverinfo = "text",
                text = apply(metabPCA2018, 1, function(l) {
                  paste(names(l), ":", l, collapse = "\n")
                }))

plotly::plot_ly(data = metabPCA2018,
                x = ? PC3,
                y = ? PC4,
                color = ? block,
                type = "scatter",
                mode = "markers",
                hoverinfo = "text",
                text = apply(metabPCA2018, 1, function(l) {
                  paste(names(l), ":", l, collapse = "\n")
                }))



# Figure #
nPC <- 30
pcaResTotal2018 <- pcaMethods::pca(metab2018Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2018 <- pcaResTotal2018@R2
pdf(paste0(dirMidSTAMPCA, scriptID, "_PCA_Explained_variance_of_Metab_2018.pdf"))
barplot(pcaR2Total2018, xlab = "PC", ylab = "R^2")
dev.off()






###### 4. Perform PCA and adjust p-values for Root Metabolomic data in 2020 ######
##### 4.1. Read all Root Metabolomic data in 2020 into R #####
metab2020Raw <- read.csv("data/phenotype/2020_Tottori_Main_Metabolome_No_Outlier_thres=5.csv")
See(metab2020Raw, rown = 6, coln = 12)

metab2020Raw <- metab2020Raw[metab2020Raw$block == "W1" | metab2020Raw$block == "W4", ]


metabStart <- 10
metabEnd <- ncol(metab2020Raw)
# metabNames <- colnames(metab2020Raw)[metabStart:metabEnd]
# varNames <- c("Vu", "Ve", "heritInd", "heritLine")
# varietyNames <- unique(metab2020Raw$variety)
# blockNames <- unique(metab2020Raw$block)
#
# nMetab <- metabEnd - metabStart + 1
# nVar <- length(varNames)
# nVariety <- length(varietyNames)
# nRep <- length(unique(metab2020Raw$ind))
# nBlock <- length(blockNames)


##### 4.2. Perform PCA in all Root Metabolomic data in 2020 #####
#### 4.2.1. nPC = 6 ####
nPC <- 6
metab2020Scaled <- scale(metab2020Raw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotal2020 <- pcaMethods::pca(metab2020Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2020 <- pcaResTotal2020@R2
barplot(pcaR2Total2020)
pcaR2CumTotal2020 <- pcaResTotal2020@R2cum
pcaScoresTotal2020 <- pcaResTotal2020@scores

head(pcaScoresTotal2020)
table(is.na(pcaScoresTotal2020))

metabPCA2020 <- data.frame(metab2020Raw[, 1:(metabStart - 1)],
                           pcaScoresTotal2020)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_nPC_6_root_metab_2020.csv")
write.csv(x = metabPCA2020, file = fileName)



#### 4.2.2. nPC = 20 ####
nPC <- 20
metab2020Scaled <- scale(metab2020Raw[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResTotal2020 <- pcaMethods::pca(metab2020Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2020 <- pcaResTotal2020@R2
barplot(pcaR2Total2020)
pcaR2CumTotal2020 <- pcaResTotal2020@R2cum
pcaScoresTotal2020 <- pcaResTotal2020@scores

head(pcaScoresTotal2020)
table(is.na(pcaScoresTotal2020))

metabPCA2020 <- data.frame(metab2020Raw[, 1:(metabStart - 1)],
                           pcaScoresTotal2020)


fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_nPC_20_root_metab_2020.csv")
write.csv(x = metabPCA2020, file = fileName)



##### 4.3. Perform PCA in Root Metabolomic data for Control and Drought in 2020 #####

#### 4.3.1. Control, nPC = 6 ####
nPC <- 6

metab2020Control <- metab2020Raw[metab2020Raw$block == "W1", ]
See(metab2020Control, coln = 12)
metab2020Scaled <- scale(metab2020Control[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResControl2020 <- pcaMethods::pca(metab2020Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Control2020 <- pcaResControl2020@R2
barplot(pcaR2Control2020)
pcaR2CumControl2020 <- pcaResControl2020@R2cum
pcaScoresControl2020 <- pcaResControl2020@scores

head(pcaScoresControl2020)
table(is.na(pcaScoresControl2020))

metabPCA2020 <- data.frame(metab2020Control[, 1:(metabStart - 1)],
                           pcaScoresControl2020)
See(metabPCA2020)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_nPC_6_root_metab_Control_2020.csv")
write.csv(x = metabPCA2020, file = fileName)



#### 4.3.2. Control, nPC = 20 ####
nPC <- 20

metab2020Control <- metab2020Raw[metab2020Raw$block == "W1", ]
metab2020Scaled <- scale(metab2020Control[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResControl2020 <- pcaMethods::pca(metab2020Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Control2020 <- pcaResControl2020@R2
barplot(pcaR2Control2020)
pcaR2CumControl2020 <- pcaResControl2020@R2cum
pcaScoresControl2020 <- pcaResControl2020@scores

head(pcaScoresControl2020)
table(is.na(pcaScoresControl2020))

metabPCA2020 <- data.frame(metab2020Control[, 1:(metabStart - 1)],
                           pcaScoresControl2020)
See(metabPCA2020)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_nPC_20_root_metab_Control_2020.csv")
write.csv(x = metabPCA2020, file = fileName)



#### 4.4.1. Drought, nPC = 6 ####
nPC <- 6

metab2020Drought <- metab2020Raw[metab2020Raw$block == "W4", ]
See(metab2020Drought, coln = 12)
metab2020Scaled <- scale(metab2020Drought[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResDrought2020 <- pcaMethods::pca(metab2020Scaled,
                                     method = "ppca",
                                     nPcs = nPC,
                                     completeObs = FALSE)
pcaR2Drought2020 <- pcaResDrought2020@R2
barplot(pcaR2Drought2020)
pcaR2CumDrought2020 <- pcaResDrought2020@R2cum
pcaScoresDrought2020 <- pcaResDrought2020@scores

head(pcaScoresDrought2020)
table(is.na(pcaScoresDrought2020))

metabPCA2020 <- data.frame(metab2020Drought[, 1:(metabStart - 1)],
                           pcaScoresDrought2020)
See(metabPCA2020)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_nPC_6_root_metab_Drought_2020.csv")
write.csv(x = metabPCA2020, file = fileName)





#### 4.4.2. Drought, nPC = 20 ####
nPC <- 20

metab2020Drought <- metab2020Raw[metab2020Raw$block == "W4", ]
See(metab2020Drought, coln = 12)
metab2020Scaled <- scale(metab2020Drought[, metabStart:metabEnd], center = TRUE, scale = TRUE)
pcaResDrought2020 <- pcaMethods::pca(metab2020Scaled,
                                     method = "ppca",
                                     nPcs = nPC,
                                     completeObs = FALSE)
pcaR2Drought2020 <- pcaResDrought2020@R2
barplot(pcaR2Drought2020)
pcaR2CumDrought2020 <- pcaResDrought2020@R2cum
pcaScoresDrought2020 <- pcaResDrought2020@scores

head(pcaScoresDrought2020)
table(is.na(pcaScoresDrought2020))

metabPCA2020 <- data.frame(metab2020Drought[, 1:(metabStart - 1)],
                           pcaScoresDrought2020)
See(metabPCA2020)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_nPC_20_root_metab_Drought_2020.csv")
write.csv(x = metabPCA2020, file = fileName)





###### 10. Perform PCA by prcomp() and ppca() and adjust p-values for Metabolomic data in 2017 ######
##### 10.1. Read Metabolomic data in 2017 into R #####
metab2017Raw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
See(metab2017Raw, rown = 6, coln = 12)


metabStart <- 10
metabEnd <- ncol(metab2017Raw)
metabNames <- colnames(metab2017Raw)[metabStart:metabEnd]
# varNames <- c("Vu", "Ve", "heritInd", "heritLine")
# varietyNames <- unique(metab2017Raw$variety)
# blockNames <- unique(metab2017Raw$block)

nMetab <- metabEnd - metabStart + 1
# nVar <- length(varNames)
# nVariety <- length(varietyNames)
# nRep <- length(unique(metab2017Raw$ind))
# nBlock <- length(blockNames)

table(is.na(metab2017Raw[, metabStart:metabEnd]))
See(na.omit(metab2017Raw[, metabStart:metabEnd]))
metab2017RawOnlyMetabNoNA <- na.omit(metab2017Raw[, metabStart:metabEnd])
See(metab2017RawOnlyMetabNoNA)

apply(metab2017Raw[, metabStart:metabEnd], 1, function(x){
  table(is.na(x))[2] >= 1
  })
table(apply(metab2017Raw[, metabStart:metabEnd], 1, function(x){
  table(is.na(x))[2] >= 1
  }))

# naNumMetab2017Raw <- which(is.na(apply(metab2017Raw[, metabStart:metabEnd], 1, function(x){
#   table(is.na(x))[2] >= 1
#   })))
# metab2017RawNoNA <- metab2017Raw[naNumMetab2017Raw, ]
# See(metab2017RawNoNA)

##### 10.2. Perform PCA in Metabolomic data in 2017 #####
#### 10.2.1 by prcomp() ####
# nPC <- 6
naNumMetab2017Raw <- which(is.na(apply(metab2017Raw[, metabStart:metabEnd], 1, function(x){
  table(is.na(x))[2] >= 1
})))
metab2017RawNoNA <- metab2017Raw[naNumMetab2017Raw, ]
See(metab2017RawNoNA)


metab2017Scaled <- scale(metab2017RawOnlyMetabNoNA , center = TRUE, scale = TRUE)
pcaResTotal2017Prcomp <- prcomp(metab2017Scaled)
summary(pcaResTotal2017Prcomp)

explainedVariance <- (pcaResTotal2017Prcomp$sdev)^2 / sum((pcaResTotal2017Prcomp$sdev)^2)
cumulativeVariance <- cumsum(explained_variance)
variance <- data.frame(
  Principal_Component = paste0("PC", 1:length(explained_variance)),
  Explained_Variance = explained_variance,
  Cumulative_Variance = cumulative_variance
)

barplot(variance$Explained_Variance[1:30])

pcaScoresTotal2017Prcomp <- pcaResTotal2017Prcomp$x
See(pcaScoresTotal2017Prcomp)
table(is.na(pcaScoresTotal2017))

metabPCA2017Prcomp <- data.frame(metab2017RawNoNA[, 1:(metabStart - 1)],
                                 pcaScoresTotal2017)
See(metabPCA2017Prcomp, coln = 12)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_prcomp_PCA_metab_without_NA_2017.csv")
write.csv(x = metabPCA2017Prcomp, file = fileName)



#### 10.2.2. by ppca()  ####
nPC <- 100

naNumMetab2017Raw <- which(is.na(apply(metab2017Raw[, metabStart:metabEnd], 1, function(x){
  table(is.na(x))[2] >= 1
})))
metab2017RawNoNA <- metab2017Raw[naNumMetab2017Raw, ]
See(metab2017RawNoNA)


metab2017Scaled <- scale(metab2017RawOnlyMetabNoNA, center = TRUE, scale = TRUE)
pcaResTotal2017PpcaNoNA <- pcaMethods::pca(metab2017Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE
                                   )
pcaR2Total2017PpcaNoNA <- pcaResTotal2017PpcaNoNA@R2
barplot(pcaR2Total2017PpcaNoNA)
pcaR2CumTotal2017 <- pcaResTotal2017PpcaNoNA@R2cum
pcaScoresTotal2017 <- pcaResTotal2017PpcaNoNA@scores

See(pcaScoresTotal2017)
table(is.na(pcaScoresTotal2017))

metabPCA2017 <- data.frame(metab2017RawNoNA[, 1:(metabStart - 1)],
                           pcaScoresTotal2017)
See(metabPCA2017, coln = 12)

fileName <- paste0(dirMidSTAMPCA, scriptID, "_pcaMethods_PCA_without_NA_metab_2017.csv")
write.csv(x = metabPCA2017, file = fileName)


# Figure #
nPC <- 30
pcaResTotal2017 <- pcaMethods::pca(metab2017Scaled,
                                   method = "ppca",
                                   nPcs = nPC,
                                   completeObs = FALSE)
pcaR2Total2017 <- pcaResTotal2017@R2
pdf(paste0(dirMidSTAMPCA, scriptID, "_PCA_Explained_variance_of_Metab_2017.pdf"))
barplot(pcaR2Total2017, xlab = "PC", ylab = "R^2")
dev.off()
















