##########################################################################################
######  Title: 2.20_Soybean_STAM_Factor_loading_for_PCA_based_on_metabolomic_data with or without flavonoid   ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/04 (Created), 2024/04/04 (Last Updated)                       ######
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

scriptID <- "2.20"


##### 1.2. Setting some parameters #####


dirMidSTAMBase <- "midstream/"

dirMidSTAMMetabFactorLBasedOnPCAFlavonoid <- paste0(dirMidSTAMBase, scriptID, "_Factor_loading_for_PCA_With_Or_Without_Flavonoid/")

dir.create(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
install.packages('openxlsx')

library(openxlsx)







###### 2. Select metabolites with higher factor loading ######
##### 2.1. Flavonoid Metabolites #####
##### ##### 2.1.1 Metabolites data without row containing NA #####
Metab2017FlavonoidNoOutlier <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
head(Metab2017FlavonoidNoOutlier)
See(Metab2017FlavonoidNoOutlier)

MetabStart <- 10
MetabEnd <- ncol(Metab2017FlavonoidNoOutlier)
OnlyMetab2017FlavonoidNoOutlier <- Metab2017FlavonoidNoOutlier[, MetabStart:MetabEnd]
# table(is.na(OnlyMetab2017FlavonoidNoOutlier))
dim(OnlyMetab2017FlavonoidNoOutlier)

missing <- apply(is.na(OnlyMetab2017FlavonoidNoOutlier), 1, sum) > 0
table(missing)
OnlyMetab2017FlavonoidNoOutlierNoNA <- OnlyMetab2017FlavonoidNoOutlier[!missing, ]

# table(is.na(OnlyMetab2017FlavonoidNoOutlierNoNA))
See(OnlyMetab2017FlavonoidNoOutlierNoNA)
dim(OnlyMetab2017FlavonoidNoOutlierNoNA)



##### 2.2. PCscores without rows containing NA based on OnlyMetab2017FlavonoidNoOutlier #####
pcaMethodsPCAMetab2017Flavonoid <- read.csv("midstream/2.16_PCA_with_or_without_flavonoid/2.16_pcaMethods_PCA_metab_related_to_flavonoid_pathway_2017.csv")
head(pcaMethodsPCAMetab2017Flavonoid)
See(pcaMethodsPCAMetab2017Flavonoid)

PCAMetabStart <- 11
PCAMetabEnd <- ncol(pcaMethodsPCAMetab2017Flavonoid)
OnlyPCAscoreMetab2017Flavonoid <- pcaMethodsPCAMetab2017Flavonoid[, PCAMetabStart:PCAMetabEnd]
table(is.na(OnlyPCAscoreMetab2017Flavonoid))
See(OnlyPCAscoreMetab2017Flavonoid)

Metab2017FlavonoidNoOutlier <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway.csv")
head(Metab2017FlavonoidNoOutlier)
MetabStart <- 10
MetabEnd <- ncol(Metab2017FlavonoidNoOutlier)
OnlyMetab2017FlavonoidNoOutlier <- Metab2017FlavonoidNoOutlier[, 10:MetabEnd]
missing <- apply(is.na(OnlyMetab2017FlavonoidNoOutlier), 1, sum) > 0
OnlyPCAscoreMetab2017Flavonoid <- OnlyPCAscoreMetab2017Flavonoid[!missing, ]
table(missing)
# table(is.na(OnlyPCAscoreMetab2017Flavonoid))
# See(OnlyPCAscoreMetab2017Flavonoid)





##### 2.3. Factor loadings of Flavonoid-related metabolites for each PC #####
FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
# FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid)
head(FactorLoadings)
class(FactorLoadings)
#sort(FactorLoadings)
#sort(FactorLoadings[,1])




FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
class(FactorLoadings)
MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)
PCNames <- colnames(FactorLoadings)
for (PCNo in 1:length(PCNames)){

  FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadings[, PCNo]
  MetabWithTop20PCNoFactorLoadings <- FactorLoadings[names(sort(abs(FactorLoadingsBasedOnMetabPCNoScore), decreasing = T)[1:20]), ]
  MetabWithTop20PCNoFactorLoadings <- as.data.frame(MetabWithTop20PCNoFactorLoadings)
  MetabWithTop20PCNoFactorLoadings <- select(MetabWithTop20PCNoFactorLoadings, PCNo)
  MetabWithTop20PCNoFactorLoadings <- arrange(MetabWithTop20PCNoFactorLoadings, desc(eval(parse(text = paste0("PC",PCNo)))))
  MetabNamesWithTop20PCNoFactorLoadings <- rownames(MetabWithTop20PCNoFactorLoadings)

  MetabAnnotationWithTop20PCNoFactorLoadings <- MetabAnnotation[MetabNamesWithTop20PCNoFactorLoadings, ]
  MetabAnnotationWithTop20PCNoFactorLoadings <- as.data.frame(MetabAnnotationWithTop20PCNoFactorLoadings)
  rownames(MetabAnnotationWithTop20PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
  names(MetabAnnotationWithTop20PCNoFactorLoadings) <- "Annotation"

  MetabAnnotationWithTop20PCNoFactorLoadings[, 2] <- MetabWithTop20PCNoFactorLoadings
  colnames(MetabAnnotationWithTop20PCNoFactorLoadings)[2] <- paste0("PC",PCNo,"Factorloading")

  fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")

  write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
}





# Plot factor loadings flavonoid related #
# FactorLoadingsBasedOnMetabFlavonoidPC1Score <- FactorLoadings[,1]
# FactorLoadingsBasedOnMetabFlavonoidPC2Score <- FactorLoadings[,2]
# FactorLoadingsBasedOnMetabFlavonoidPC3Score <- FactorLoadings[,3]
# FactorLoadingsBasedOnMetabFlavonoidPC4Score <- FactorLoadings[,4]
# FactorLoadingsBasedOnMetabFlavonoidPC6Score <- FactorLoadings[,6]
# head(FactorLoadings)
#


PCNames <- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo1 in PCNames){
  PCNamsesWithoutPCNo <- PCNames[-(1:which(PCNames == PCNo1))]

  for (PCNo2 in PCNamsesWithoutPCNo){
    pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_", PCNo1, "_", PCNo2, "_flavonoid_related.pdf"))
    plot(FactorLoadings[, PCNo1], FactorLoadings[, PCNo2], main = "FactorLoading", xlab = paste0(PCNo1), ylab = paste0(PCNo2))
    abline( h = 0, v = 0 )
    dev.off()
  }
}

# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC1_PC2_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC1Score, FactorLoadingsBasedOnMetabFlavonoidPC2Score, main = "FactorLoading", xlab = "PC1", ylab = "PC2")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC1_PC3_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC1Score, FactorLoadingsBasedOnMetabFlavonoidPC3Score, main = "FactorLoading", xlab = "PC1", ylab = "PC3")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC1_PC4_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC1Score, FactorLoadingsBasedOnMetabFlavonoidPC4Score, main = "FactorLoading", xlab = "PC1", ylab = "PC4")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC1_PC6_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC1Score, FactorLoadingsBasedOnMetabFlavonoidPC6Score, main = "FactorLoading", xlab = "PC1", ylab = "PC6")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC2_PC3_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC2Score, FactorLoadingsBasedOnMetabFlavonoidPC3Score, main = "FactorLoading", xlab = "PC2", ylab = "PC3")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC2_PC4_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC2Score, FactorLoadingsBasedOnMetabFlavonoidPC4Score, main = "FactorLoading", xlab = "PC2", ylab = "PC4")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC2_PC6_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC2Score, FactorLoadingsBasedOnMetabFlavonoidPC6Score, main = "FactorLoading", xlab = "PC2", ylab = "PC6")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC3_PC4_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC3Score, FactorLoadingsBasedOnMetabFlavonoidPC4Score, main = "FactorLoading", xlab = "PC3", ylab = "PC4")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "Factor_loading_PC3_PC6_flavonoid_related.pdf"))
# plot(FactorLoadingsBasedOnMetabFlavonoidPC3Score, FactorLoadingsBasedOnMetabFlavonoidPC6Score, main = "FactorLoading", xlab = "PC3", ylab = "PC6")
# abline( h = 0, v = 0 )
# dev.off()



# FactorLoadingsDF <- as.data.frame(FactorLoadings)
# FactorLoadingsDF[order(FactorLoadingsDF$PC1), ]
# FactorLoadingsDF[order(FactorLoadingsDF$PC1)[1:20], ]
# FactorLoadingsDF[order(FactorLoadingsDF$PC3), ]
# FactorLoadingsDF[order(FactorLoadingsDF$PC3)[1:20], ]


##### 2.4. Choosing metabolites with Top20 PC factor loadings #####
#### 2.4.1. Top20 of PC2 absolute factor loadings ####
sort(abs(FactorLoadingsBasedOnMetabPC2Score),decreasing = T)[1:20]
str(sort(abs(FactorLoadingsBasedOnMetabPC2Score),decreasing = T)[1:20])
class(sort(abs(FactorLoadingsBasedOnMetabPC2Score),decreasing = T)[1:20])
names(sort(abs(FactorLoadingsBasedOnMetabPC2Score),decreasing = T)[1:20])



# MetabWithTop20PC2FactorLoadings <- factor.loadings[names(sort(abs(FactorLoadingsBasedOnMetabPC2Score), decreasing = T)[1:20]), ]
# class(MetabWithTop20PC2FactorLoadings)
#
# MetabWithTop20PC2FactorLoadings <- as.data.frame(MetabWithTop20PC2FactorLoadings)
# MetabNamesWithTop20PC2FactorLoadings <- rownames(MetabWithTop20PC2FactorLoadings)
# MetabWithTop20PC2FactorLoadings <- MetabWithTop20PC2FactorLoadings[, "PC2"]
# MetabWithTop20PC2FactorLoadings <- as.data.frame(MetabWithTop20PC2FactorLoadings)
# rownames(MetabWithTop20PC2FactorLoadings) <- MetabNamesWithTop20PC2FactorLoadings
# class(MetabWithTop20PC2FactorLoadings)
# MetabWithTop20PC2FactorLoadings <- sort(MetabWithTop20PC2FactorLoadings[,1], decreasing = T)
# MetabWithTop20PC2FactorLoadings <- na.omit(MetabWithTop20PC2FactorLoadings)
#
# MetabNamesWithTop20PC2FactorLoadings <- rownames(MetabWithTop20PC2FactorLoadings)





#### 2.5. Writing files on factor loadings ####
MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)
head(MetabAnnotation)
class(MetabAnnotation)





# This code also correctly function #
# FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
# class(FactorLoadings)
# MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)
# PCNames <- colnames(FactorLoadings)
# for (PCNo in 1:length(PCNames)){
#
#   FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadings[, PCNo]
#   MetabWithTop20PCNoFactorLoadings <- FactorLoadings[names(sort(abs(FactorLoadingsBasedOnMetabPCNoScore), decreasing = T)[1:20]), ]
#   MetabWithTop20PCNoFactorLoadings <- as.data.frame(MetabWithTop20PCNoFactorLoadings)
#   MetabWithTop20PCNoFactorLoadings <- select(MetabWithTop20PCNoFactorLoadings, PCNo)
#   MetabWithTop20PCNoFactorLoadings <- arrange(MetabWithTop20PCNoFactorLoadings, desc(eval(parse(text = paste0("PC",PCNo)))))
#   MetabNamesWithTop20PCNoFactorLoadings <- rownames(MetabWithTop20PCNoFactorLoadings)
#
#   MetabAnnotationWithTop20PCNoFactorLoadings <- MetabAnnotation[MetabNamesWithTop20PCNoFactorLoadings, ]
#   MetabAnnotationWithTop20PCNoFactorLoadings <- as.data.frame(MetabAnnotationWithTop20PCNoFactorLoadings)
#   rownames(MetabAnnotationWithTop20PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
#   names(MetabAnnotationWithTop20PCNoFactorLoadings) <- "Annotation"
#
#   MetabAnnotationWithTop20PCNoFactorLoadings[, 2] <- MetabWithTop20PCNoFactorLoadings
#   colnames(MetabAnnotationWithTop20PCNoFactorLoadings)[2] <- paste0("PC",PCNo,"Factorloading")
#
#   fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")
#
#   write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
# }





# This code also correctly function #
# PCNo <- 2
# FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
# class(FactorLoadings)
# MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)
# PCNames <- colnames(FactorLoadings)
# for (PCNo in 1:length(PCNames)){
#   FactorLoadingsDF <- as.data.frame(FactorLoadings)
#   FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadingsDF[order(abs(FactorLoadingsDF[, paste0("PC",PCNo)]),decreasing = T)[1:20], ]
#
#   PCNoFactorLoadings <- FactorLoadingsBasedOnMetabPCNoScore[, paste0("PC",PCNo)]
#   MetabNamesWithTop20PCNoFactorLoadings <- rownames(FactorLoadingsBasedOnMetabPCNoScore)[1:20]
#   names(PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
#   PCNoFactorLoadings <- sort(PCNoFactorLoadings, decreasing = T)
#
#   PCNoFactorLoadings <- as.data.frame(PCNoFactorLoadings)
#
#
#   MetabAnnotationWithTop20PCNoFactorLoadings <- MetabAnnotation[MetabNamesWithTop20PCNoFactorLoadings, ]
#   MetabAnnotationWithTop20PCNoFactorLoadings <- as.data.frame(MetabAnnotationWithTop20PCNoFactorLoadings)
#   rownames(MetabAnnotationWithTop20PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
#   names(MetabAnnotationWithTop20PCNoFactorLoadings) <- "Annotation"
#
#   MetabAnnotationWithTop20PCNoFactorLoadings[, 2] <- PCNoFactorLoadings
#   colnames(MetabAnnotationWithTop20PCNoFactorLoadings)[2] <- paste0("PC",PCNo,"Factorloading")
#
#
#   fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")
#
#   write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
# }




# This code also correctly function #
# PCNo <- "PC2"
FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
class(FactorLoadings)
MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)
PCNames <- colnames(FactorLoadings)
for (PCNo in PCNames){
  FactorLoadingsDF <- as.data.frame(FactorLoadings)
  FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadingsDF[order(abs(FactorLoadingsDF[, PCNo]),decreasing = T)[1:20], ]

  PCNoFactorLoadings <- FactorLoadingsBasedOnMetabPCNoScore[, PCNo]
  MetabNamesWithTop20PCNoFactorLoadings <- rownames(FactorLoadingsBasedOnMetabPCNoScore)
  names(PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
  PCNoFactorLoadings <- sort(PCNoFactorLoadings, decreasing = T)

  PCNoFactorLoadings <- as.data.frame(PCNoFactorLoadings)


  MetabAnnotationWithTop20PCNoFactorLoadings <- MetabAnnotation[MetabNamesWithTop20PCNoFactorLoadings, ]
  MetabAnnotationWithTop20PCNoFactorLoadings <- as.data.frame(MetabAnnotationWithTop20PCNoFactorLoadings)
  rownames(MetabAnnotationWithTop20PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
  names(MetabAnnotationWithTop20PCNoFactorLoadings) <- "Annotation"

  MetabAnnotationWithTop20PCNoFactorLoadings[, 2] <- PCNoFactorLoadings
  colnames(MetabAnnotationWithTop20PCNoFactorLoadings)[2] <- paste0(PCNo,"Factorloading")


  fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, scriptID,                                                    "_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings", "_with_flavonoid", ".csv")

  write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
}





### Top 30 PC1 factor loadings
PCNo <- 1
FactorLoadings <- cor(OnlyMetab2017FlavonoidNoOutlierNoNA, OnlyPCAscoreMetab2017Flavonoid, use = "pair")
class(FactorLoadings)
MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)

FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadings[, PCNo]
MetabWithTop30PCNoFactorLoadings <- FactorLoadings[names(sort(abs(FactorLoadingsBasedOnMetabPCNoScore), decreasing = T)[1:30]), ]
MetabWithTop30PCNoFactorLoadings <- as.data.frame(MetabWithTop30PCNoFactorLoadings)
MetabWithTop30PCNoFactorLoadings <- select(MetabWithTop30PCNoFactorLoadings, PCNo)
MetabWithTop30PCNoFactorLoadings <- arrange(MetabWithTop30PCNoFactorLoadings, desc(eval(parse(text = paste0("PC",PCNo)))))
MetabNamesWithTop30PCNoFactorLoadings <- rownames(MetabWithTop30PCNoFactorLoadings)

MetabAnnotationWithTop30PCNoFactorLoadings <- MetabAnnotation[MetabNamesWithTop30PCNoFactorLoadings, ]
MetabAnnotationWithTop30PCNoFactorLoadings <- as.data.frame(MetabAnnotationWithTop30PCNoFactorLoadings)
rownames(MetabAnnotationWithTop30PCNoFactorLoadings) <- MetabNamesWithTop30PCNoFactorLoadings
names(MetabAnnotationWithTop30PCNoFactorLoadings) <- "Annotation"

MetabAnnotationWithTop30PCNoFactorLoadings[, 2] <- MetabWithTop30PCNoFactorLoadings
colnames(MetabAnnotationWithTop30PCNoFactorLoadings)[2] <- paste0("PC",PCNo,"Factorloading")

fileMetabAnnotationWithTop30PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCAFlavonoid, scriptID,                                                    "_MetabInfo_of_Top30_PC", PCNo, "_Factor_Loadings_with_flavonoid.csv")

write.csv(x = MetabAnnotationWithTop30PCNoFactorLoadings, file = fileMetabAnnotationWithTop30PCNoFactorLoadings)

