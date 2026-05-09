##########################################################################################
######  Title: 2.23_Soybean_STAM_Boxplot_of_Amount_of_Metab_with_or_without_Flavonoid_for_different_Markers ######
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

scriptID <- "2.23"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"



dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites <- paste0(dirMidSTAMBase, scriptID, "_Boxplot_Amount_of_flavonoid_for_each_marker_or_different_markers/")
dir.create(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites)
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




###### 2. Read genotypic values of Flavonoid-related Metabolites and read marker genotype #####
##### 2.1. Read genotypic values of all metabolites in 2017 into R #####
gvMetab2017 <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
See(gvMetab2017)



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




###### 3. Boxplot of amount of Flavonoid-related Metabolites #####
##### 3.1. For each marker #####
#### Chr06_18760995(pcaMethods, nPC = 6, PC1),
Chr10_42562665(pcaMethods, nPC = 6, PC4)




metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]

MetabNo <- "X200014"
markerInterestID <- "Chr15_45643039"
markers <- c("Chr06_18760995", "Chr10_42562665", "Chr06_47486371","Chr17_16065902", "Chr17_17152757", "Chr02_29956148", "Chr11_4271460", "Chr05_39845221", "Chr14_7215537", "Chr19_33074985", "Chr14_4306528", "Chr20_37715838", "Chr18_6213704", "Chr05_1716802", "Chr04_6507472", "Chr15_45643039", "Chr18_6592763", "Chr05_1724003", "Chr03_35410010", "Chr20_43338174", "Chr07_6053910", "Chr07_19701087", "Chr12_31046936")

for (markerInterestID in markers){

  dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites, scriptID, "_For_each_marker/"))
  dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_", markerInterestID, "/"))

  for (MetabNo in metabNamesFlavonoid){
    LineNames <- rownames(gvMetab2017)
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

    markerInterest[, 2] <- CommonNames
    markerInterest_0 <- markerInterest[markerInterest$markerInterest == 0, ]

    markerInterest_2 <- markerInterest[markerInterest$markerInterest == 2, ]
    See(markerInterest_0)

    Metab[, 2] <- CommonNames

    gvMarkerInterest_0 <- Metab[rownames(markerInterest_0), ]
    gvMarkerInterest_2 <- Metab[rownames(markerInterest_2), ]

    See(gvMarkerInterest_0)
    See(gvMarkerInterest_2)

    pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_", markerInterestID,"/", scriptID, "_", MetabNo, ".pdf"))
    boxplot(gvMarkerInterest_0$Metab, gvMarkerInterest_2$Metab, names = c("0",  "2"))
    dev.off()

  }
}

hist(gvMarkerInterest_0$Metab)
hist(gvMarkerInterest_2$Metab)




##### 3.2. For different markers #####
#### 3.2.1. For combination of marker values, (0,0),(0,2),(2,0),(2,2) ####
MetabNo <- "X00015"
markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr10_42562665"
metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites, scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", "marker_value_combination/"))
# dir.create(paste0(dirMidSTAMMetabAmountBasedOnPCsFlavonoid, scriptID, "_",  markerInterestID, "_mainly_for_PC1/", scriptID, "_", PCNo))

for (MetabNo in metabNamesFlavonoid){
  LineNames <- rownames(gvMetab2017)
  markerInterest1 <- gastonData0Matrix[, markerInterestID1]
  markerInterest2 <- gastonData0Matrix[, markerInterestID2]

  markerInterest1DF <- as.data.frame(markerInterest1)
  markerInterest2DF <- as.data.frame(markerInterest2)

  markerVarietyNames <- names(markerInterest1)
  LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
  CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

  markerInterest1DF <- markerInterest1DF[CommonNames, ]
  names(markerInterest1DF) <- CommonNames
  markerInterest1DF <- as.data.frame(markerInterest1DF)
  markerInterest2DF <- markerInterest2DF[CommonNames, ]
  names(markerInterest2DF) <- CommonNames
  markerInterest2DF <- as.data.frame(markerInterest2DF)

  markerInterest1And2DF <- cbind(markerInterest1DF, markerInterest2DF)
  See(markerInterest1And2DF)

  marker0_0 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 0, ]
  marker0_2 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 2, ]
  marker2_0 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 0, ]
  marker2_2 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 2, ]
  See(marker0_0)
  See(marker0_2)
  See(marker2_0)
  See(marker2_2)
  varietyNames0_0 <- rownames(marker0_0)
  varietyNames0_2 <- rownames(marker0_2)
  varietyNames2_0 <- rownames(marker2_0)
  varietyNames2_2 <- rownames(marker2_2)

  Metab <- as.data.frame(gvMetab2017[, MetabNo])

  rownames(Metab) <- rownames(gvMetab2017)
  colnames(Metab) <- MetabNo

  Metab0_0 <- Metab[varietyNames0_0, ]
  Metab0_0 <- as.data.frame(Metab0_0)
  rownames(Metab0_0) <- varietyNames0_0
  Metab0_2 <- Metab[varietyNames0_2, ]
  Metab0_2 <- as.data.frame(Metab0_2)
  rownames(Metab0_2) <- varietyNames0_2
  Metab2_0 <- Metab[varietyNames2_0, ]
  Metab2_0 <- as.data.frame(Metab2_0)
  rownames(Metab2_0) <- varietyNames2_0
  Metab2_2 <- Metab[varietyNames2_2, ]
  Metab2_2 <- as.data.frame(Metab2_2)
  rownames(Metab2_2) <- varietyNames2_2

  pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersFlavonoidMetabolites, "2.23_Chr06_18760995_Chr10_42562665_marker_value_combination/", scriptID, "_", MetabNo, "_marker_value_combination.pdf"))
  boxplot(Metab0_0$Metab0_0, Metab0_2$Metab0_2, Metab2_0$Metab2_0, Metab2_2$Metab2_2, names = c("0_0",  "0_2", "2_0", "2_2"))
  dev.off()
}



##### 4.4. Characterize Lines of each marker based on group information #####
markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr10_42562665"

groupInfo <- read.csv("midstream/2.30_Haplotype_based_PhyloTree_and_Network/0.3_group_information.csv", row.names = 1)
# groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"
See(groupInfo)


LineNames <- rownames(gvMetab2017)
markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]

markerInterest1DF <- as.data.frame(markerInterest1)
markerInterest2DF <- as.data.frame(markerInterest2)

markerVarietyNames <- names(markerInterest1)
LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

markerInterest1DF <- markerInterest1DF[CommonNames, ]
names(markerInterest1DF) <- CommonNames
markerInterest1DF <- as.data.frame(markerInterest1DF)
markerInterest2DF <- markerInterest2DF[CommonNames, ]
names(markerInterest2DF) <- CommonNames
markerInterest2DF <- as.data.frame(markerInterest2DF)

markerInterest1And2DF <- cbind(markerInterest1DF, markerInterest2DF)
See(markerInterest1And2DF)

marker0_0 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 0, ]
marker0_2 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 2, ]
marker2_0 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 0, ]
marker2_2 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 2, ]
See(marker0_0)
See(marker0_2)
See(marker2_0)
See(marker2_2)
varietyNames0_0 <- rownames(marker0_0)
varietyNames0_2 <- rownames(marker0_2)
varietyNames2_0 <- rownames(marker2_0)
varietyNames2_2 <- rownames(marker2_2)


# groupInfoLineNames <- groupInfo$lineNames
# rownames(groupInfo) <- groupInfoLineNames
#
# groupInfo0_0 <- groupInfo[varietyNames0_0, ]
# groupInfo0_2 <- groupInfo[varietyNames0_2, ]
# groupInfo2_0 <- groupInfo[varietyNames2_0, ]
# groupInfo2_2 <- groupInfo[varietyNames2_2, ]
# table(groupInfo0_0$group)
# table(groupInfo0_2$group)
# table(groupInfo2_0$group)
# table(groupInfo2_2$group)

groupInfo0_0 <- groupInfo[groupInfo$lineNames %in% varietyNames0_0, ]
groupInfo0_2 <- groupInfo[groupInfo$lineNames %in% varietyNames0_2, ]
groupInfo2_0 <- groupInfo[groupInfo$lineNames %in% varietyNames2_0, ]
groupInfo2_2 <- groupInfo[groupInfo$lineNames %in% varietyNames2_2, ]
table(groupInfo0_0$group)
table(groupInfo0_2$group)
table(groupInfo2_0$group)
table(groupInfo2_2$group)





###### 3. Boxplot of amount of Flavonoid Non-related Metabolites for different markers #####








