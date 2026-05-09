##########################################################################################
######  Title: 2.21_Soybean_STAM_amount_of_metab_based_on_PCs_factor_loadings_with_or_without_flavonoid  ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/05 (Created), 2024/04/05 (Last Updated)                       ######
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

scriptID <- "2.21"


##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
# sigLevel <- 0.05

dirMidSTAMMetabAmountBasedOnPCsFlavonoid <- paste0(dirMidSTAMBase, scriptID,
                                "_Amount_of_metab_based_on_PCs_Flavonoid/")
dir.create(dirMidSTAMMetabAmountBasedOnPCsFlavonoid)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)






###### 2. Check amount of several metabolites in flavonoid pathway ######
##### 2.1. Read data with all metabolites in 2017 into R #####
gvMetab2017 <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv")

See(gvMetab2017)

#gvMetab2017$X200014[1:10]
#gvMetab2017$X250001[1:10]
#gvMetab2017$X500080[1:10]
#gvMetab2017$X00422[1:10]




##### 2.2. Read marker genotype into R #####
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
# markerInterest <- gastonData0Matrix[, "Chr06_47426527"]

# aho <- c(29, 14, 53, 102, 431, 4, 329, 213)
# baka <- c(0, 1, 0, 2, 2, 0, 2, 2)
# boxplot(aho ? baka)
# #baka:genotype





# ### flavonoid-related PC6
# markerInterestID <- "Chr06_47426527"
#
# markerInterest <- gastonData0Matrix[, markerInterestID]
# See(markerInterest)
#
# Metab1ID <- "X00867"
# Metab2ID <- "X500009"
#
# Metab1 <- as.data.frame(gvMetab2017[, Metab1ID])   # Saponarin
# Metab2 <- as.data.frame(gvMetab2017[, Metab2ID])   # Synephrine
# LineNames <- gvMetab2017[, "X"]
#
# rownames(Metab1) <- LineNames
# colnames(Metab1) <- "X250001"
# rownames(Metab2) <- LineNames
# colnames(Metab2) <- "X00422"
#
# length(markerInterest)
# See(markerInterest)
# # markerInterest <- as.data.frame(markerInterest)
# # markerVarietyNames <- rownames(markerInterest)
# markerVarietyNames <- names(markerInterest)
# LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
# CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]
# table(LineNames%in%markerVarietyNames)
# table(markerVarietyNames%in%LineNames)
#
# # df_mergeNames%in%CommonNames
# # table(df_mergeNames%in%CommonNames)
#
# Metab1 <- Metab1[CommonNames, ]
# Metab2 <- Metab2[CommonNames, ]
# See(Metab1)
# See(Metab2)
# Metab1 <- as.data.frame(Metab1)
# Metab2 <- as.data.frame(Metab2)
# rownames(Metab1) <- CommonNames
# rownames(Metab2) <- CommonNames
#
# markerInterest <- as.data.frame(markerInterest)
# colnames(markerInterest) <- markerInterestID
# markerInterest <- markerInterest[CommonNames, ]
# markerInterest <- as.data.frame(markerInterest)
# rownames(markerInterest) <- CommonNames
#
# # Metab1 <- as.vector(Metab1)
# # Metab2 <- as.vector(Metab2)
# # markerInterest <- as.vector(markerInterest)
# # class(Metab1)
# # boxplot(Metab1 ? markerInterest)
#
#
# boxplot(Metab1$Metab1 ? markerInterest$markerInterest)
# boxplot(Metab2$Metab2 ? markerInterest$markerInterest)
# head(Metab1)
# head(markerInterest)
# ?boxplot
# class(Metab1$Metab1)



### flavonoid-related, highest marker for PC1
markerInterestID <- "Chr06_18760995"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC1/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC1/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC1/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



### flavonoid-related, higher marker for PC2
## left
markerInterestID <- "Chr06_21382495"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC2/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC2/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC2/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}


## middle
markerInterestID <- "Chr06_28785446"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC2/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC2/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC2/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}


## right
markerInterestID <- "Chr06_37826501"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC2/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC2/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC2/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}




### flavonoid-related, higher marker for PC3
## 6th chromosome
markerInterestID <- "Chr06_47606869"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



## 18th chromosome
markerInterestID <- "Chr18_6592763"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_", markerInterestID, "_mainly_for_PC3/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



## 5th chromosome
markerInterestID <- "Chr05_39845221"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



## 2th chromosome
markerInterestID <- "Chr02_5754920"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC3/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}




### flavonoid-related, highest marker for PC4
# metabFlavonoidPC4 <- read.csv("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_PC4_Factor_Loadings_with_flavonoid.csv")
# See(metabFlavonoidPC4)
# metabNamesFlavonoidPC4 <- metabFlavonoidPC4[, "X"]


markerInterestID <- "Chr10_42562665"
See(markerInterest)
# PCNo <- "PC1"
# MetabNo <- "X500134"
PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC4/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC4/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC4/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}




### flavonoid-related, higher marker for PC6
## 6th chromosome in the latter
markerInterestID <- "Chr06_47426527"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



## 6th chromosome in the former
markerInterestID <- "Chr06_15443311"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



## 15th chromosome
markerInterestID <- "Chr15_45643650"
See(markerInterest)

PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"))
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}



## 3rd chromosome
markerInterestID <- "Chr03_41369779"
See(markerInterest)

PCNo <- "PC1"
MetabNo <-
PCs<- c("PC1", "PC2", "PC3", "PC4", "PC6")
for (PCNo in PCs){
  metabFlavonoidPCNo <- read.csv(paste0("midstream/2.20_Factor_loading_for_PCA_With_Or_Without_Flavonoid/2.20_MetabInfo_of_Top20_", PCNo, "_Factor_Loadings_with_flavonoid.csv"))
  metabNamesFlavonoidPCNo <- metabFlavonoidPCNo[, "X"]


  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/"))
  dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo))

  for (MetabNo in metabNamesFlavonoidPCNo){
    LineNames <- gvMetab2017[, "X"]
    markerInterest <- gastonData0Matrix[, markerInterestID]

    Metab <- as.data.frame(gvMetab2017[, MetabNo])

    rownames(Metab) <- LineNames
    colnames(Metab) <- MetabNo

    markerVarietyNames <- names(markerInterest)
    LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
    CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

    Metab <- Metab[CommonNames, ]
    Metab <- as.data.frame(Metab)
    rownames(Metab) <- CommonNames

    markerInterest <- as.data.frame(markerInterest)
    colnames(markerInterest) <- markerInterestID
    markerInterest <- markerInterest[CommonNames, ]
    markerInterest <- as.data.frame(markerInterest)
    rownames(markerInterest) <- CommonNames

    pdf(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID,  "_",markerInterestID, "_mainly_for_PC6/", scriptID, "_", PCNo, "/", MetabNo, ".pdf"), )
    boxplot(Metab$Metab ? markerInterest$markerInterest, main = MetabNo, xlab = "Marker value", ylab = "Amount")
    dev.off()
  }
}




