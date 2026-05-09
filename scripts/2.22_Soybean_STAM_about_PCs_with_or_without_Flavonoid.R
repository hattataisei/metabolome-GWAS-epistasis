##########################################################################################
######  Title: 2.22_Soybean_STAM_about_PCs_with_or_without_Flavonoid                ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/22 (Created), 2024/04/22 (Last Updated)                       ######
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

scriptID <- "2.22"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"



dirMidSTAMAboutPCsFlavonoidMetabolites <- paste0(dirMidSTAMBase, scriptID, "_About_PCs_with_or_without_flavonoid/")
dir.create(dirMidSTAMAboutPCsFlavonoidMetabolites)
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




###### 2. Analysis of PCs from PCA for flavonoid-related data in 2017 ######
##### 2.1. with Flavonoid #####
#### 2.1.1. Read genotipic value based on Flavonoid PCA data in 2017 into R ####
gvMetab2017FlavonoidPCA <- read.csv("midstream/2.18_BSH_for_PCA_with_or_without_flavonoid/2.18_lmer_genotypic_values_for_PC_Score_for_flavonoid_related_metab_in_2017.csv", row.names = 1)
class(gvMetab2017FlavonoidPCA)
head(gvMetab2017FlavonoidPCA)
plot(gvMetab2017FlavonoidPCA$PC1,gvMetab2017FlavonoidPCA$PC2, xlab = "PC1Score", ylab = "PC2Score")
abline(h = 0, v = 0)

varietyNames <- rownames(gvMetab2017FlavonoidPCA)

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



#### 2.1.2. Read marker genotype into R ####
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
markerInterestID <- "Chr06_18760995"
markerInterest <- gastonData0Matrix[, markerInterestID]
See(markerInterest)

markerVarietyNames <- names(markerInterest)


#### 2.1.3. Plot PCs with different colors ####
# CommonVarietyNames <- varietyNames[varietyNames %in% markerVarietyNames]
# gvMetab2017FlavonoidPCAOnlyWithCommonVariety <- gvMetab2017FlavonoidPCA[CommonVarietyNames, ]
# See(gvMetab2017FlavonoidPCAOnlyWithCommonVariety)
#
# markerInterest <- as.data.frame(markerInterest)
# markerInterest <- markerInterest[CommonVarietyNames, ]
# markerInterest <- as.data.frame(markerInterest)
# rownames(markerInterest) <- CommonVarietyNames
# See(markerInterest)
#
# gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest <- cbind(gvMetab2017FlavonoidPCAOnlyWithCommonVariety,markerInterest)
#
# pdf("midstream/2.21_Amount_of_metab_based_on_PCs_Flavonoid/2.21_plot_lmer_PC1_PC2_flavonoid.pdf")
# ggplot(gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest) + geom_point(mapping = aes(x = PC1, y = PC2, color = markerInterest))
# dev.off()




# PCNo1 <- "PC1"
# PCNo2 <- "PC2"
markerInterestID <- "Chr06_18760995"

dir.create(paste0(dirMidSTAMAboutPCsFlavonoidMetabolites, scriptID, "_plot_lmer_PCs_with_flavonoid"))
PCs <- colnames(gvMetab2017FlavonoidPCA)
markerInterestIDs <- c("Chr06_18760995", "Chr10_42562665", "Chr06_47426527")
for (PCNo1 in PCs){
  PCsWithoutUntilPCNo1 <- PCs[-(1:which(PCNo1 == PCs))]

  for (PCNo2 in PCsWithoutUntilPCNo1){

    for (markerInterestID in markerInterestIDs ){
      varietyNames <- rownames(gvMetab2017FlavonoidPCA)
      markerInterest <- gastonData0Matrix[, markerInterestID]
      markerVarietyNames <- names(markerInterest)
      CommonVarietyNames <- varietyNames[varietyNames %in% markerVarietyNames]
      gvMetab2017FlavonoidPCAOnlyWithCommonVariety <- gvMetab2017FlavonoidPCA[CommonVarietyNames, ]
      See(gvMetab2017FlavonoidPCAOnlyWithCommonVariety)

      markerInterest <- as.data.frame(markerInterest)
      markerInterest <- markerInterest[CommonVarietyNames, ]
      markerInterest <- as.data.frame(markerInterest)
      rownames(markerInterest) <- CommonVarietyNames
      See(markerInterest)

      gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest <- cbind(gvMetab2017FlavonoidPCAOnlyWithCommonVariety, markerInterest)

      pdf(paste0(dirMidSTAMAboutPCsFlavonoidMetabolites, scriptID,  "_plot_lmer_PCs_flavonoid/", scriptID, "_", markerInterestID, "_", PCNo1, "_", PCNo2, ".pdf"))
      p <- ggplot(gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest) + geom_point(mapping = aes(x = eval(parse(text = paste0(PCNo1))), y = eval(parse(text = paste0(PCNo2))), color = markerInterest))
      print(p)
      dev.off()
    }
  }
}



##### 2.2 Without Flavonoid #####



