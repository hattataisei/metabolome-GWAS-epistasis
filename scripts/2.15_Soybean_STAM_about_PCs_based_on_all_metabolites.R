##########################################################################################
######  Title: 2.15_Soybean_STAM_about_PCs                            ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/02 (Created), 2024/04/02 (Last Updated)                       ######
##########################################################################################





###### 1. Settings ######
##### 1.0. Reset workspace ######
# rm(list=ls())



##### 1.1. Setting working directory to the "projectName" directory #####
cropName <- "soybean"
project <- "STAM"
os <- osVersion

isRproject <- function(path = getwd()) {
  files <- list.files(path)

  if (length(grep(".Rproj", files)) >= 1) {
    out <- TRUE
  } else {
    out <-  FALSE
  }
  return(out)
}

scriptID <- "2.15"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"



dirMidSTAMAboutPCsAllmetabolites <- paste0(dirMidSTAMBase, scriptID,
                        "_About_PCs_of_all_metabolites/")
dir.create(dirMidSTAMAboutPCsAllmetabolites)
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





###### 2. Analysis of PCs from PCA for metabolomic data in 2017 ######
##### 2.1.1. Read PCA data in 2017 into R #####
# Metab2017PCA <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_2017.csv")
# class(Metab2017PCA)
# head(Metab2017PCA)
# plot(Metab2017PCA$PC1,Metab2017PCA$PC2)


gvMetab2017PCA <- read.csv("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_for_PC_Score_in_2017.csv", row.names = 1)
class(gvMetab2017PCA)
head(gvMetab2017PCA)
plot(gvMetab2017PCA$PC1,gvMetab2017PCA$PC2)

varietyNames <- rownames(gvMetab2017PCA)

# For Metab2017PCA #
# metabStart <- 11
# metabEnd <- ncol(Metab2017PCA)
# metabNames <- colnames(Metab2017PCA)[metabStart:metabEnd]
# varietyNames <- unique(Metab2017PCA$variety)
# blockNames <- unique(Metab2017PCA$block)
#
# nMetab <- metabEnd - metabStart + 1
# nVariety <- length(varietyNames)
# nRep <- length(unique(Metab2017PCA$ind))
# nBlock <- length(blockNames)



##### 2.1.2. Read marker genotype into R #####
gastonData0 <- gaston::read.vcf(file = "raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz")


chrNos <- gastonData0@snps$chr
chrNames <- sprintf(fmt = paste0("Chr",
                                 "%0", floor(log10(max(chrNos))) + 1, "i"),
                    chrNos)
pos <- gastonData0@snps$pos
mrkNames <- paste0(chrNames, "_", pos)
gastonData0@snps$id <- mrkNames

gastonData0Matrix <- as.matrix(gastonData0)
See(gastonData0Matrix)
markerInterest <- gastonData0Matrix[, "Chr06_18760995"]
See(markerInterest)

markerVarietyNames <- names(markerInterest)


##### 2.1.3. Plot PCs with different colors #####
CommonVarietyNames <- varietyNames[varietyNames %in% markerVarietyNames]
gvMetab2017PCAOnlyWithCommonVariety <- gvMetab2017PCA[CommonVarietyNames, ]
See(gvMetab2017PCAOnlyWithCommonVariety)

markerInterest <- as.data.frame(markerInterest)
markerInterest <- markerInterest[CommonVarietyNames, ]
markerInterest <- as.data.frame(markerInterest)
rownames(markerInterest) <- CommonVarietyNames
See(markerInterest)

gvMetab2017PCAOnlyWithCommonVarietyWithmarkerInterest <- cbind(gvMetab2017PCAOnlyWithCommonVariety,markerInterest)

pdf("midstream/2.15_About_PCs_of_All_Metabolites/2.15_plot_lmer_PC1_PC2.pdf")
ggplot(gvMetab2017PCAOnlyWithCommonVarietyWithmarkerInterest) + geom_point(mapping = aes(x = PC1, y = PC2, color = markerInterest))
dev.off()


