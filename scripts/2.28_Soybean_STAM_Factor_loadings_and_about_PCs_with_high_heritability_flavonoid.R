##########################################################################################
######  Title: 2.28_Soybean_STAM_Factor_loadings_and_about_PCs_with_high_heritability_flavonoid         ######
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

scriptID <- "2.28"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"



dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites <- paste0(dirMidSTAMBase, scriptID, "_Factor_loadings_and_about_PCs_with_high_heritability_flavonoid/")
dir.create(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
install.packages('openxlsx')

require(openxlsx)
require(data.table)
require(RAINBOWR)
require(gaston)
require(ggplot2)
require(tidyverse)
require(plotly)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### 2. Analysis of PCs from PCA for high heritability Flavonoid in 2017 ######
#### 2.1. Read genotipic value based on Flavonoid PCA data in 2017 into R ####
gvMetab2017FlavonoidPCA <- read.csv("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_2017_pcaMethods_nPC=6.csv", row.names = 1)
class(gvMetab2017FlavonoidPCA)
head(gvMetab2017FlavonoidPCA)
plot(gvMetab2017FlavonoidPCA$PC1,gvMetab2017FlavonoidPCA$PC2, xlab = "PC1Score", ylab = "PC2Score")
abline(h = 0, v = 0)

varietyNames <- rownames(gvMetab2017FlavonoidPCA)



#### 2.2. Read marker genotype into R ####
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




#### 2.3. Plot PCs with different colors ####
# PCNo1 <- "PC1"
# PCNo2 <- "PC2"
# markerInterestID <- "Chr06_18760995"

dir.create(paste0(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites, scriptID, "_plot_lmer_PCScores/"))
PCs <- colnames(gvMetab2017FlavonoidPCA)
markerInterestIDs <- c("Chr06_18760995", "Chr10_42562665", "Chr06_47426527", "Chr17_16065902")
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
      # colnames(markerInterest) <- markerInterestID

      markerInterest$markerInterest <- factor(markerInterest$markerInterest)
      See(markerInterest)

      gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest <- cbind(gvMetab2017FlavonoidPCAOnlyWithCommonVariety, markerInterest)
      See(gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest, coln = 7)

      pdf(paste0(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites, scriptID,  "_plot_lmer_PCScores/", scriptID, "_", markerInterestID, "_", PCNo1, "_", PCNo2, ".pdf"))
      p <- ggplot(data = gvMetab2017FlavonoidPCAOnlyWithCommonVarietyWithmarkerInterest) + geom_point(mapping = aes(x = eval(parse(text = paste0(PCNo1))), y = eval(parse(text = paste0(PCNo2))), color = markerInterest)) + geom_hline( yintercept = 0 ) + geom_vline( xintercept = 0) + labs(title = "PCScore", x = PCNo1, y = PCNo2) + theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title.x  = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), legend.title = element_text(size = 15),legend.text = element_text(size = 15) ) + scale_color_hue(name = print(markerInterestID) )
      print(p)
      dev.off()
    }
  }
}


### PC Scores with C, D labels
See(pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9, coln = 15)
See(pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9$block)
pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9$block <- factor(pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9$block)
block <- pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9$block
ggplot(data = pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9, aes(x = PC4, y = PC6, color = block)) + geom_point() + geom_hline( yintercept = 0 ) + geom_vline( xintercept = 0) + labs(title = "PCScore", x = "PC4", y = "PC2") + theme(plot.title = element_text(hjust = 0.5) ) + scale_color_hue(name = "treatment" )








###### 3. Select metabolites with higher factor loading ######
##### 3.1. Metabolites data without row containing NA #####
Metab2017FlavonoidHeritabilityMoreThan0.9NoOutlier <- read.csv("data/phenotype/_2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(Metab2017FlavonoidHeritabilityMoreThan0.9NoOutlier, coln = 15)

MetabStart <- 11
MetabEnd <- ncol(Metab2017FlavonoidHeritabilityMoreThan0.9NoOutlier)
OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlier <- Metab2017FlavonoidHeritabilityMoreThan0.9NoOutlier[, MetabStart:MetabEnd]
# table(is.na(OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlier))
dim(OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlier)

missing <- apply(is.na(OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlier), 1, sum) > 0
table(missing)
OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlierNoNA <- OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlier[!missing, ]
table(is.na(OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlierNoNA))
See(OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlierNoNA)



##### 3.2. PCscores without rows containing NA based on OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlier #####
pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9 <- read.csv("midstream/2.24_PCA_with_high_heritability_flavonoid/2.24_pcaMethods_PCA_nPC=6_flavonoid_metab_>0.9_Heritability_2017.csv")
See(pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9, coln = 15)

PCAMetabStart <- 12
PCAMetabEnd <- ncol(pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9)
OnlyPCAscoreMetab2017Flavonoid <- pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9[, PCAMetabStart:PCAMetabEnd]
table(is.na(OnlyPCAscoreMetab2017Flavonoid))
See(OnlyPCAscoreMetab2017Flavonoid)

OnlyPCAscoreMetab2017Flavonoid <- OnlyPCAscoreMetab2017Flavonoid[!missing, ]
table(is.na(OnlyPCAscoreMetab2017Flavonoid))
See(OnlyPCAscoreMetab2017Flavonoid)

plot(pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9$PC1, pcaMethodsPCAMetab2017FlavonoidHeritabilityMoreThan0.9$PC4)



##### 3.3. Factor loadings for each PC #####
FactorLoadings <- cor(OnlyMetab2017FlavonoidHeritabilityMoreThan0.9NoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
See(FactorLoadings)

MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)

dir.create(paste0(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites, scriptID, "_factor_loading_for_PCs/"))
PCNames <- colnames(FactorLoadings)
PCNo <- 2
for (PCNo in 1:length(PCNames)){

  FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadings[, PCNo]
  FactorLoadingsBasedOnMetabPCNoScore <- as.data.frame(FactorLoadingsBasedOnMetabPCNoScore)
  FactorLoadingsBasedOnMetabPCNoScore <- arrange(FactorLoadingsBasedOnMetabPCNoScore, desc(FactorLoadingsBasedOnMetabPCNoScore))

  MetabAnnotationPCNoFactorLoadings <- MetabAnnotation[rownames(FactorLoadingsBasedOnMetabPCNoScore), ]
  MetabAnnotationPCNoFactorLoadings <- as.data.frame(MetabAnnotationPCNoFactorLoadings)
  rownames(MetabAnnotationPCNoFactorLoadings) <- rownames(FactorLoadingsBasedOnMetabPCNoScore)
  names(MetabAnnotationPCNoFactorLoadings) <- "Annotation"

  MetabAnnotationPCNoFactorLoadings[, 2] <- FactorLoadingsBasedOnMetabPCNoScore[, 1]
  colnames(MetabAnnotationPCNoFactorLoadings)[2] <- paste0("PC",PCNo,"Factorloading")
  See(MetabAnnotationPCNoFactorLoadings)

  fileMetabAnnotationPCNoFactorLoadings <- paste0(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites, scriptID, "_factor_loading_for_PCs/", scriptID, "_MetabInfo_of_PC", PCNo, "_factor_loadings_for_flavonoid_heritability_>0.9.csv")

  write.csv(x = MetabAnnotationPCNoFactorLoadings, file = fileMetabAnnotationPCNoFactorLoadings)
}


### plot factor loadings
PCNames <- c("PC1", "PC2", "PC3", "PC4", "PC6")
dir.create(paste0(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites, scriptID, "_plot_factor_loadings/"))
for (PCNo1 in PCNames){
  PCNamsesWithoutPCNo <- PCNames[-(1:which(PCNames == PCNo1))]

  for (PCNo2 in PCNamsesWithoutPCNo){
    pdf(paste0(dirMidSTAMFactorLoadingsAndAboutPCsHighHeritabilityFlavonoidMetabolites, scriptID, "_plot_factor_loadings/", "Factor_loading_", PCNo1, "_", PCNo2, ".pdf"))
    plot(FactorLoadings[, PCNo1], FactorLoadings[, PCNo2], main = "FactorLoading", xlab = paste0(PCNo1), ylab = paste0(PCNo2))
    abline( h = 0, v = 0 )
    dev.off()
  }
}



### for PC1, PC2
See(MetabAnnotationPCNoFactorLoadingsPC1)
MetabAnnotationPCNoFactorLoadingsPC1 <- MetabAnnotationPCNoFactorLoadings
MetabAnnotationPCNoFactorLoadingsPC2 <- MetabAnnotationPCNoFactorLoadings
MetabAnnotationPCNoFactorLoadingsPC1PC2 <- merge(x = MetabAnnotationPCNoFactorLoadingsPC1, y = MetabAnnotationPCNoFactorLoadingsPC2, by = "Annotation")
See(MetabAnnotationPCNoFactorLoadingsPC1PC2)

metabPC1PC2FactorloadingForPC1Positive <- MetabAnnotationPCNoFactorLoadingsPC1PC2[MetabAnnotationPCNoFactorLoadingsPC1PC2$PC1Factorloading > 0, ]
metabPC1PC2FactorloadingForPC1Negative <- MetabAnnotationPCNoFactorLoadingsPC1PC2[MetabAnnotationPCNoFactorLoadingsPC1PC2$PC1Factorloading < 0, ]
See(metabPC1PC2FactorloadingForPC1Positive)
See(metabPC1PC2FactorloadingForPC1Negative)

metabPC1PC2FactorloadingForPC1Positive <- arrange(metabPC1PC2FactorloadingForPC1Positive, desc = PC2Factorloading)
metabPC1PC2FactorloadingForPC1Negative <- arrange(metabPC1PC2FactorloadingForPC1Negative, desc = PC2Factorloading)



