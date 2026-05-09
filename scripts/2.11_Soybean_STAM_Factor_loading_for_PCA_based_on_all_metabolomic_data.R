##########################################################################################
######  Title: 2.11_Soybean_STAM_Factor_loading_for_PCA_based_on_metabolomic_data   ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/07 (Created), 2024/03/13 (Last Updated)                       ######
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

scriptID <- "2.11"


##### 1.2. Setting some parameters #####
# nPC <- 6

dirMidSTAMBase <- "midstream/"

dirMidSTAMMetabFactorLBasedOnPCA <- paste0(dirMidSTAMBase, scriptID,
                                "_Factor_loading_for_PCA/")
dir.create(dirMidSTAMMetabFactorLBasedOnPCA)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
install.packages('openxlsx')

library(openxlsx)







###### 2. Selecting metabolites with higher loading factor ######
##### 2.1. Metabolites data without row containing NA #####
Metab2017NoOutlier <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
head(Metab2017NoOutlier)

MetabEnd <- ncol(Metab2017NoOutlier)
OnlyMetab2017NoOutlier <- Metab2017NoOutlier[, 10:MetabEnd]
# table(is.na(OnlyMetab2017NoOutlier))
dim(OnlyMetab2017NoOutlier)

missing <- apply(is.na(OnlyMetab2017NoOutlier), 1, sum) > 0
OnlyMetab2017NoOutlierNoNA <- OnlyMetab2017NoOutlier[!missing, ]

# table(is.na(OnlyMetab2017NoOutlier))
dim(OnlyMetab2017NoOutlierNoNA)



##### 2.2. PCscores without rows containing NA based on OnlyMetab2017NoOutlier #####


pcaMethodsPCAMetab2017 <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_2017.csv")
head(pcaMethodsPCAMetab2017)
str(pcaMethodsPCAMetab2017)

PCAMetabEnd <- ncol(pcaMethodsPCAMetab2017)
OnlyPCAscoreMetab2017 <- pcaMethodsPCAMetab2017[, 11:PCAMetabEnd]

Metab2017NoOutlier <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
MetabEnd <- ncol(Metab2017NoOutlier)
OnlyMetab2017NoOutlier <- Metab2017NoOutlier[, 10:MetabEnd]
missing <- apply(is.na(OnlyMetab2017NoOutlier), 1, sum) > 0
OnlyPCAscoreMetab2017 <- OnlyPCAscoreMetab2017[!missing, ]

# table(is.na(OnlyPCAscoreMetab2017))
# str(OnlyPCAscoreMetab2017)




##### 2.3. Factor loadings of all metabolites for each PC #####
dim(OnlyMetab2017NoOutlierNoNA)

FactorLoadings <- cor(OnlyMetab2017NoOutlierNoNA, OnlyPCAscoreMetab2017, use = "pair")
# FactorLoadings <- cor(OnlyMetab2017NoOutlierNoNA, OnlyPCAscoreMetab2017)


head(FactorLoadings)
# loadings based on ppca?
#sort(FactorLoadings)
class(FactorLoadings)

#sort(factor.loadings[,1])




FactorLoadings <- cor(OnlyMetab2017NoOutlierNoNA, OnlyPCAscoreMetab2017, use = "pair")
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

  fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")

  write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
}





### Plot factor loadings
# FactorLoadingsBasedOnMetabPC1Score <- FactorLoadings[,1]
# FactorLoadingsBasedOnMetabPC2Score <- FactorLoadings[,2]
# FactorLoadingsBasedOnMetabPC3Score <- FactorLoadings[,3]
# FactorLoadingsBasedOnMetabPC6Score <- FactorLoadings[,6]
# head(FactorLoadings)



PCNames <- c("PC1", "PC2", "PC3", "PC6")
for (PCNo1 in PCNames){
  PCNamsesWithoutPCNo <- PCNames[-(1:which(PCNames == PCNo1))]

  for (PCNo2 in PCNamsesWithoutPCNo){
    pdf(paste0(dirMidSTAMMetabFactorLBasedOnPCA, "Factor_loading_", PCNo1, "_", PCNo2, ".pdf"))
    plot(FactorLoadings[, PCNo1], FactorLoadings[, PCNo2], main = "FactorLoading", xlab = paste0(PCNo1), ylab = paste0(PCNo2))
    abline( h = 0, v = 0 )
    dev.off()
  }
}



# pdf("Factor_loading_PC1_PC2.pdf")
# plot(FactorLoadingsBasedOnMetabPC1Score, FactorLoadingsBasedOnMetabPC2Score, main = "FactorLoading", xlab = "PC1", ylab = "PC2")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf("Factor_loading_PC1_PC3.pdf")
# plot(FactorLoadingsBasedOnMetabPC1Score, FactorLoadingsBasedOnMetabPC3Score, main = "FactorLoading", xlab = "PC1", ylab = "PC3")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf("Factor_loading_PC1_PC6.pdf")
# plot(FactorLoadingsBasedOnMetabPC1Score, FactorLoadingsBasedOnMetabPC6Score, main = "FactorLoading", xlab = "PC1", ylab = "PC6")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf("Factor_loading_PC2_PC3.pdf")
# plot(FactorLoadingsBasedOnMetabPC2Score, FactorLoadingsBasedOnMetabPC3Score, main = "FactorLoading", xlab = "PC2", ylab = "PC3")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf("Factor_loading_PC2_PC6.pdf")
# plot(FactorLoadingsBasedOnMetabPC2Score, FactorLoadingsBasedOnMetabPC6Score, main = "FactorLoading", xlab = "PC2", ylab = "PC6")
# abline( h = 0, v = 0 )
# dev.off()
#
# pdf("Factor_loading_PC3_PC6.pdf")
# plot(FactorLoadingsBasedOnMetabPC3Score, FactorLoadingsBasedOnMetabPC6Score, main = "FactorLoading", xlab = "PC3", ylab = "PC6")
# abline( h = 0, v = 0 )
# dev.off()
#
# ?plot
#
# FactorLoadingsDF <- as.data.frame(FactorLoadings)
# FactorLoadingsDF[order(FactorLoadingsDF$PC1), ]
# FactorLoadingsDF[order(FactorLoadingsDF$PC1)[1:20], ]
# FactorLoadingsDF[order(FactorLoadingsDF$PC3), ]
# FactorLoadingsDF[order(FactorLoadingsDF$PC3)[1:20], ]


##### 2.4. Choosing metabolites with Top20 PC factor loadings #####
# #### 2.4.1. Factor loadings : PC1 > 0.6  &  -0.2 < PC2 < 0.2 ####
# MetabFactorLoadingsMoreRelatedToPC1FactorLoadings <- factor.loadings[FactorLoadingsBasedOnMetabPC1Score > 0.6 & abs(FactorLoadingsBasedOnMetabPC2Score) < 0.2, ]
# #plot(MetabFactorLoadingsMoreRelatedToPC1FactorLoadings[, 1], MetabFactorLoadingsMoreRelatedToPC1FactorLoadings[, 2], xlim = c(-1,1), ylim = c(-1, 1))
# MetabMoreRelatedToPC1FactorLoadings <- rownames(MetabFactorLoadingsMoreRelatedToPC1FactorLoadings)
#
#
#
# #### 2.4.2 Factor loadings : PC2 > 0.6  &  -0.2 < PC1 < 0.2 ####
# MetabFactorLoadingsMoreRelatedToPC2FactorLoadings <- factor.loadings[FactorLoadingsBasedOnMetabPC2Score > 0.6 & abs(FactorLoadingsBasedOnMetabPC1Score) < 0.2, ]
# #plot(MetabFactorLoadingsMoreRelatedToPC2FactorLoadings[, 1], MetabFactorLoadingsMoreRelatedToPC2FactorLoadings[, 2], xlim = c(-1,1), ylim = c(-1, 1))
# MetabMoreRelatedToPC2FactorLoadings <- rownames(MetabFactorLoadingsMoreRelatedToPC2FactorLoadings)
#
#
#
# #### 2.4.3. Factor loadings : PC1 < -0.6, 0.6 < PC1 ####
# MetabWithHigherPC1FactorLoadings <- factor.loadings[FactorLoadingsBasedOnMetabPC1Score > 0.6, ]
# MetabNamesWithHigherPC1FactorLoadings <- rownames(MetabWithHigherPC1FactorLoadings)
#
# MetabWithLowerPC1FactorLoadings <- factor.loadings[  FactorLoadingsBasedOnMetabPC1Score < -0.6, ]
# MetabNamesWithLowerPC1FactorLoadings <- rownames(MetabWithLowerPC1FactorLoadings)
#
#
# #### 2.4.4. Factor loadings : PC2 < -0.6, 0.6 < PC2 ####
# MetabWithHigherPC2FactorLoadings <- factor.loadings[FactorLoadingsBasedOnMetabPC2Score > 0.6, ]
# MetabNamesWithHigherPC2FactorLoadings <- rownames(MetabWithHigherPC2FactorLoadings)
#
# MetabWithLowerPC2FactorLoadings <- factor.loadings[  FactorLoadingsBasedOnMetabPC2Score < -0.6, ]
# MetabNamesWithLowerPC2FactorLoadings <- rownames(MetabWithLowerPC2FactorLoadings)
#
#
# #### 2.4.5. Factor loadings : PC3 < -0.5, 0.5 < PC3 ####
# PC3 <- sort(factor.loadings[, "PC3"])
#
# MetabWithHigherPC3FactorLoadings <- factor.loadings[FactorLoadingsBasedOnMetabPC3Score > 0.5, ]
# MetabNamesWithHigherPC3FactorLoadings <- rownames(MetabWithHigherPC3FactorLoadings)
#
# MetabWithLowerPC3FactorLoadings <- factor.loadings[  FactorLoadingsBasedOnMetabPC3Score < -0.5, ]
# MetabNamesWithLowerPC3FactorLoadings <- rownames(MetabWithLowerPC3FactorLoadings)
#
#
# #### 2.4.6. Factor loadings : PC6 < -0.4, 0.4 < PC6 ####
# PC6 <- sort(factor.loadings[, "PC6"])
#
# MetabWithHigherPC6FactorLoadings <- factor.loadings[FactorLoadingsBasedOnMetabPC6Score > 0.4, ]
# MetabNamesWithHigherPC6FactorLoadings <- rownames(MetabWithHigherPC6FactorLoadings)
#
# MetabWithLowerPC6FactorLoadings <- factor.loadings[  FactorLoadingsBasedOnMetabPC6Score < -0.4, ]
# MetabNamesWithLowerPC6FactorLoadings <- rownames(MetabWithLowerPC6FactorLoadings)




#### 2.4.7. Top20 of PC2 absolute factor loadings ####
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


# Mainly PC1 #
# MetabAnnotationMainlyRelatedToPC1FactorLoadings <- MetabAnnotation[MetabMoreRelatedToPC1FactorLoadings, ]
# MetabAnnotationMainlyRelatedToPC1FactorLoadings <- as.data.frame(MetabAnnotationMainlyRelatedToPC1FactorLoadings)
# rownames(MetabAnnotationMainlyRelatedToPC1FactorLoadings) <- MetabMoreRelatedToPC1FactorLoadings
# names(MetabAnnotationMainlyRelatedToPC1FactorLoadings) <- "Annotation"
#
#
# # Mainly PC2 #
# MetabAnnotationMainlyRelatedToPC2FactorLoadings <- MetabAnnotation[MetabMoreRelatedToPC2FactorLoadings, ]
# MetabAnnotationMainlyRelatedToPC2FactorLoadings <- as.data.frame(MetabAnnotationMainlyRelatedToPC2FactorLoadings)
# rownames(MetabAnnotationMainlyRelatedToPC2FactorLoadings) <- MetabMoreRelatedToPC2FactorLoadings
# names(MetabAnnotationMainlyRelatedToPC2FactorLoadings) <- "Annotation"
#
#
# # Higher PC1 #
# MetabAnnotationWithHigherPC1FactorLoadings <- MetabAnnotation[MetabNamesWithHigherPC1FactorLoadings, ]
# MetabAnnotationWithHigherPC1FactorLoadings <- as.data.frame(MetabAnnotationWithHigherPC1FactorLoadings)
# rownames(MetabAnnotationWithHigherPC1FactorLoadings) <- MetabNamesWithHigherPC1FactorLoadings
# names(MetabAnnotationWithHigherPC1FactorLoadings) <- "Annotation"
#
# # Lower PC1 #
# MetabAnnotationWithLowerPC1FactorLoadings <- MetabAnnotation[MetabNamesWithLowerPC1FactorLoadings, ]
# MetabAnnotationWithLowerPC1FactorLoadings <- as.data.frame(MetabAnnotationWithLowerPC1FactorLoadings)
# rownames(MetabAnnotationWithLowerPC1FactorLoadings) <- MetabNamesWithLowerPC1FactorLoadings
# names(MetabAnnotationWithLowerPC1FactorLoadings) <- "Annotation"
#
#
# # Higher PC2 #
# MetabAnnotationWithHigherPC2FactorLoadings <- MetabAnnotation[MetabNamesWithHigherPC2FactorLoadings, ]
# MetabAnnotationWithHigherPC2FactorLoadings <- as.data.frame(MetabAnnotationWithHigherPC2FactorLoadings)
# rownames(MetabAnnotationWithHigherPC2FactorLoadings) <- MetabNamesWithHigherPC2FactorLoadings
# names(MetabAnnotationWithHigherPC2FactorLoadings) <- "Annotation"
#
# # Lower PC2 #
# MetabAnnotationWithLowerPC2FactorLoadings <- MetabAnnotation[MetabNamesWithLowerPC2FactorLoadings, ]
# MetabAnnotationWithLowerPC2FactorLoadings <- as.data.frame(MetabAnnotationWithLowerPC2FactorLoadings)
# rownames(MetabAnnotationWithLowerPC2FactorLoadings) <- MetabNamesWithLowerPC2FactorLoadings
# names(MetabAnnotationWithLowerPC2FactorLoadings) <- "Annotation"
#
#
# # Higher PC3 #
# MetabAnnotationWithHigherPC3FactorLoadings <- MetabAnnotation[MetabNamesWithHigherPC3FactorLoadings, ]
# MetabAnnotationWithHigherPC3FactorLoadings <- as.data.frame(MetabAnnotationWithHigherPC3FactorLoadings)
# rownames(MetabAnnotationWithHigherPC3FactorLoadings) <- MetabNamesWithHigherPC3FactorLoadings
# names(MetabAnnotationWithHigherPC3FactorLoadings) <- "Annotation"
#
# # Lower PC3 #
# MetabAnnotationWithLowerPC3FactorLoadings <- MetabAnnotation[MetabNamesWithLowerPC3FactorLoadings, ]
# MetabAnnotationWithLowerPC3FactorLoadings <- as.data.frame(MetabAnnotationWithLowerPC3FactorLoadings)
# rownames(MetabAnnotationWithLowerPC3FactorLoadings) <- MetabNamesWithLowerPC3FactorLoadings
# names(MetabAnnotationWithLowerPC3FactorLoadings) <- "Annotation"
#
#
# # Higher PC6 #
# MetabAnnotationWithHigherPC6FactorLoadings <- MetabAnnotation[MetabNamesWithHigherPC6FactorLoadings, ]
# MetabAnnotationWithHigherPC6FactorLoadings <- as.data.frame(MetabAnnotationWithHigherPC6FactorLoadings)
# rownames(MetabAnnotationWithHigherPC6FactorLoadings) <- MetabNamesWithHigherPC6FactorLoadings
# names(MetabAnnotationWithHigherPC6FactorLoadings) <- "Annotation"
#
# # Lower PC6 #
# MetabAnnotationWithLowerPC6FactorLoadings <- MetabAnnotation[MetabNamesWithLowerPC6FactorLoadings, ]
# MetabAnnotationWithLowerPC6FactorLoadings <- as.data.frame(MetabAnnotationWithLowerPC6FactorLoadings)
# rownames(MetabAnnotationWithLowerPC6FactorLoadings) <- MetabNamesWithLowerPC6FactorLoadings
# names(MetabAnnotationWithLowerPC6FactorLoadings) <- "Annotation"




## file names ##
# fileMetabAnnotationMainlyRelatedToPC1FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC1>0.6_and_-0.2<PC2<0.2.csv")
#
# fileMetabAnnotationMainlyRelatedToPC2FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC2>0.6_and_-0.2<PC1<0.2.csv")
#
#
# fileMetabAnnotationWithHigherPC1FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC1>0.6.csv")
# fileMetabAnnotationWithLowerPC1FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC1<-0.6.csv")
#
# fileMetabAnnotationWithHigherPC2FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC2>0.6.csv")
# fileMetabAnnotationWithLowerPC2FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC2<-0.6.csv")
#
# fileMetabAnnotationWithHigherPC3FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC3>0.5.csv")
# fileMetabAnnotationWithLowerPC3FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC3<-0.5.csv")
#
# fileMetabAnnotationWithHigherPC6FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC6>0.4,.csv")
# fileMetabAnnotationWithLowerPC6FactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, scriptID,
#                 "_MetabInfo_of_Factor_Loadings_for_PC6<-0.4,.csv")




## write files ##
# write.csv(x = MetabAnnotationMainlyRelatedToPC1FactorLoadings, file = fileMetabAnnotationMainlyRelatedToPC1FactorLoadings)
# write.csv(x = MetabAnnotationMainlyRelatedToPC2FactorLoadings, file = fileMetabAnnotationMainlyRelatedToPC2FactorLoadings)
#
# write.csv(x = MetabAnnotationWithHigherPC1FactorLoadings, file = fileMetabAnnotationWithHigherPC1FactorLoadings)
# write.csv(x = MetabAnnotationWithLowerPC1FactorLoadings, file = fileMetabAnnotationWithLowerPC1FactorLoadings)
#
# write.csv(x = MetabAnnotationWithHigherPC2FactorLoadings, file = fileMetabAnnotationWithHigherPC2FactorLoadings)
# write.csv(x = MetabAnnotationWithLowerPC2FactorLoadings, file = fileMetabAnnotationWithLowerPC2FactorLoadings)
#
# write.csv(x = MetabAnnotationWithHigherPC3FactorLoadings, file = fileMetabAnnotationWithHigherPC3FactorLoadings)
# write.csv(x = MetabAnnotationWithLowerPC3FactorLoadings, file = fileMetabAnnotationWithLowerPC3FactorLoadings)
#
# write.csv(x = MetabAnnotationWithHigherPC6FactorLoadings, file = fileMetabAnnotationWithHigherPC6FactorLoadings)
# write.csv(x = MetabAnnotationWithLowerPC6FactorLoadings, file = fileMetabAnnotationWithLowerPC6FactorLoadings)





# This code also correctly function #
# FactorLoadings <- cor(OnlyMetab2017NoOutlierNoNA, OnlyPCAscoreMetab2017, use = "pair")
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
#   fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")
#
#   write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
# }




# This code also correctly function #
# PCNo <- 2
FactorLoadings <- cor(OnlyMetab2017NoOutlierNoNA, OnlyPCAscoreMetab2017, use = "pair")
class(FactorLoadings)
MetabAnnotation <- read.xlsx("raw_data/extra/Metabolome_README.xlsx", rowNames = T)
PCNames <- colnames(FactorLoadings)
for (PCNo in 1:length(PCNames)){
  FactorLoadingsDF <- as.data.frame(FactorLoadings)
  FactorLoadingsBasedOnMetabPCNoScore <- FactorLoadingsDF[order(abs(FactorLoadingsDF[, paste0("PC",PCNo)]),decreasing = T)[1:20], ]

  PCNoFactorLoadings <- FactorLoadingsBasedOnMetabPCNoScore[, paste0("PC",PCNo)]
  MetabNamesWithTop20PCNoFactorLoadings <- rownames(FactorLoadingsBasedOnMetabPCNoScore)[1:20]
  names(PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
  PCNoFactorLoadings <- sort(PCNoFactorLoadings, decreasing = T)

  PCNoFactorLoadings <- as.data.frame(PCNoFactorLoadings)


  MetabAnnotationWithTop20PCNoFactorLoadings <- MetabAnnotation[MetabNamesWithTop20PCNoFactorLoadings, ]
  MetabAnnotationWithTop20PCNoFactorLoadings <- as.data.frame(MetabAnnotationWithTop20PCNoFactorLoadings)
  rownames(MetabAnnotationWithTop20PCNoFactorLoadings) <- MetabNamesWithTop20PCNoFactorLoadings
  names(MetabAnnotationWithTop20PCNoFactorLoadings) <- "Annotation"

  MetabAnnotationWithTop20PCNoFactorLoadings[, 2] <- PCNoFactorLoadings
  colnames(MetabAnnotationWithTop20PCNoFactorLoadings)[2] <- paste0("PC",PCNo,"Factorloading")

  fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")

  write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
}




# This code also correctly function #
# PCNo <- "PC2"
FactorLoadings <- cor(OnlyMetab2017NoOutlierNoNA, OnlyPCAscoreMetab2017, use = "pair")
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

  fileMetabAnnotationWithTop20PCNoFactorLoadings <- paste0(dirMidSTAMMetabFactorLBasedOnPCA, "2.11_Metab_with_Top20_Factor_loading/", scriptID,                                                    "_MetabInfo_of_Top20_PC", PCNo, "_Factor_Loadings.csv")

  write.csv(x = MetabAnnotationWithTop20PCNoFactorLoadings, file = fileMetabAnnotationWithTop20PCNoFactorLoadings)
}
