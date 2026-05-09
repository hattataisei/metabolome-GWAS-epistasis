##########################################################################################
######  Title: 2.31_Soybean_STAM_Genotypes_and_phenotypes_with_marker_and_area_information   ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/05/24 (Created), 2025/03/15 (Last Updated)                       ######
##########################################################################################





###### 1. Settings ######
##### 1.0. Reset workspace #####
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

scriptID <- "2.31"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMColorPhenotypeAndMarkers <- paste0(dirMidSTAMBase, scriptID, "_Color_phenotype_and_area_infomation_with_markers/")
dir.create(dirMidSTAMColorPhenotypeAndMarkers)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
require(openxlsx)
require(viridis)
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)
require(plotly)
require(manhattanly)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)



###### 2. Read information ######
##### 2.1. Genotype information #####
statusOriginInformationRaw <- read.xlsx("raw_data/extra/Supplementary_Tables.xlsx", startRow = 2)
See(statusOriginInformationRaw)
tail(statusOriginInformationRaw)

statusOriginInformation <- statusOriginInformationRaw[3:(nrow(statusOriginInformationRaw) - 1), 1:5]
See(statusOriginInformation)
tail(statusOriginInformation)


rownames(statusOriginInformation) <- statusOriginInformation$Accession.ID
See(statusOriginInformation)

unique(statusOriginInformation$Status)
unique(statusOriginInformation$Origin)

SouthEastAsia <- c("Malaysia", "Taiwan", "Philippines", "Myanmar", "Viet Nam", "Thailand", "Indonesia", "Nepal", "Laos", "East Timor", "Cambodia")
originArea <- sapply(statusOriginInformation$Origin, function(x){
  if (x %in% SouthEastAsia){
    return(paste("South East Asia"))
  }
    else {
      return(x)
  }
}
       )
See(originArea)
table(originArea)

statusOriginInformation <- cbind(statusOriginInformation, originArea)
See(statusOriginInformation)
statusOriginInformation$originArea
rownames(statusOriginInformation)

# "HOUJAKU_KUWAZU"
rownames(statusOriginInformation)[rownames(statusOriginInformation) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
rownames(statusOriginInformation)



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
rownames(gastonData0Matrix)



###### 3. Lines with status, origin and group information, and marker information ######
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/"))

# "group information" is in "2.29_"
# dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_group_information/"))
# groupInfo <- read.csv("data/extra/0.3_group_information.csv", row.names = 1)



##### 3.1. For each marker #####
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_each_marker/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_each_marker/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_each_marker/", scriptID, "_Ratio/"))

dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_each_marker/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_each_marker/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_each_marker/", scriptID, "_Ratio/"))

markerInterestIDs <- c("Chr06_18760995",  "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# markerInterestID <- "Chr06_18760995"
for (markerInterestID in markerInterestIDs){

  lineNames <- rownames(statusOriginInformation)
  markerVarietyNames <- rownames(gastonData0Matrix)
  commonNames <- markerVarietyNames[(markerVarietyNames %in% lineNames)]

  markerInterest <- gastonData0Matrix[, markerInterestID]
  markerInterest <- markerInterest[commonNames]

  statusOriginInformation <- statusOriginInformation[commonNames, ]
  See(statusOriginInformation)

  markerStatusOriginInformation <- cbind(statusOriginInformation, markerInterest)
  See(markerStatusOriginInformation, coln = 7)

  markerStatusOriginInformation$markerInterest <- factor(markerStatusOriginInformation$markerInterest)
  markerStatusOriginInformation$Status <- factor(markerStatusOriginInformation$Status)
  markerStatusOriginInformation$Origin <- factor(markerStatusOriginInformation$Origin)
  markerStatusOriginInformation$originArea <- factor(markerStatusOriginInformation$originArea)


  pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_each_marker/", scriptID, "_Number/", scriptID, "_", markerInterestID, ".pdf"))
  g <- ggplot(markerStatusOriginInformation, aes(x = markerInterest, fill = Status)) + geom_bar(stat = "count", position = "dodge") + labs(x = markerInterestID, y ="") + theme(text = element_text(size = 24))
  plot(g)
  dev.off()


  pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_each_marker/", scriptID, "_Number/", scriptID, "_", markerInterestID, ".pdf"))
  g <- ggplot(markerStatusOriginInformation, aes(x = markerInterest, fill = originArea)) + geom_bar(stat = "count", position = "dodge") + labs(x = markerInterestID, y ="") + theme(text = element_text(size = 24))
  plot(g)
  dev.off()

}




##### 3.2. For two markers combination #####
#### 3.2.1. Chr06_18760995, Chr10_42562665, Chr06_47426527, Chr17_16065902 ####
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_two_markers/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_two_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_two_markers/", scriptID, "_Ratio/"))

dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_two_markers/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_two_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_two_markers/", scriptID, "_Ratio/"))


markerInterestIDs <- c("Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# markerInterestID1 <- "Chr06_18760995"
# markerInterestID2 <- "Chr10_42562665"
for (markerInterestID1 in markerInterestIDs){
  markerInterestIDswithoutUntilMarkerInterestID1 <- markerInterestIDs[-(1:which(markerInterestID1 == markerInterestIDs))]

  for (markerInterestID2 in markerInterestIDswithoutUntilMarkerInterestID1){

    lineNames <- rownames(statusOriginInformation)
    markerInterest1 <- gastonData0Matrix[, markerInterestID1]
    markerInterest2 <- gastonData0Matrix[, markerInterestID2]

    # markerInterest1DF <- as.data.frame(markerInterest1)
    # markerInterest2DF <- as.data.frame(markerInterest2)

    markerVarietyNames <- names(markerInterest1)
    commonVarietyNames <- lineNames[(lineNames %in% markerVarietyNames)]


    markerInterest1 <- markerInterest1[commonVarietyNames]
    markerInterest2 <- markerInterest2[commonVarietyNames]

    markerInterest1And2DF <- cbind(markerInterest1, markerInterest2)
    See(markerInterest1And2DF)

    statusOriginInformation <- statusOriginInformation[commonNames, ]
    See(statusOriginInformation)

    markerStatusOriginInformation <- cbind(statusOriginInformation, markerInterest1And2DF)
    See(markerStatusOriginInformation, coln = 8)

    markerStatusOriginInformation <- mutate(markerStatusOriginInformation, group = paste(markerInterest1, markerInterest2, sep = "_"))

    markerStatusOriginInformation$markerInterest1 <- factor(markerStatusOriginInformation$markerInterest1)
    markerStatusOriginInformation$markerInterest2 <- factor(markerStatusOriginInformation$markerInterest2)
    markerStatusOriginInformation$group <- factor(markerStatusOriginInformation$group)
    markerStatusOriginInformation$Status <- factor(markerStatusOriginInformation$Status)
    markerStatusOriginInformation$Origin <- factor(markerStatusOriginInformation$Origin)
    markerStatusOriginInformation$originArea <- factor(markerStatusOriginInformation$originArea)


    a <- count(x = markerStatusOriginInformation, group, Status)
    b <- a %>% complete(group, Status, fill = list(n = 0))

    pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_two_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, ".pdf"))
    g <- ggplot(b, aes(x = group, y = n, fill = Status)) + geom_bar(stat = "identity", position = "dodge") + labs(x = paste0(markerInterestID1, "_", markerInterestID2), y ="") + theme(text = element_text(size = 15))
    plot(g)
    dev.off()


    a <- count(x = markerStatusOriginInformation, group, originArea)
    b <- a %>% complete(group, originArea, fill = list(n = 0))

    pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_two_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, ".pdf"))
    g <- ggplot(b, aes(x = group, y = n, fill = originArea)) + geom_bar(stat = "identity", position = "dodge") + labs(x = paste0(markerInterestID1, "_", markerInterestID2), y ="") + theme(text = element_text(size = 15))
    plot(g)
    dev.off()




  }
}



##### 3.3 For three markers combination #####
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_number_of_lines_with_group_information/", scriptID, "_For_three_markers/"))






##### 3.4 For four markers combination #####
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_number_of_lines_with_origin_area_information/", scriptID, "_For_four_markers/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_number_of_lines_with_status_information/", scriptID, "_For_four_markers/"))

#### 3.4.1 "Chr06_18760995", "Chr10_42562665", "Chr06_47426527", "Chr17_16065902" ####
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_four_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_four_markers/", scriptID, "_Ratio/"))

dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_four_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_four_markers/", scriptID, "_Ratio/"))

statusOriginInformationRaw <- read.xlsx("raw_data/extra/Supplementary_Tables.xlsx", startRow = 2)
See(statusOriginInformationRaw)

statusOriginInformation <- statusOriginInformationRaw[3:(nrow(statusOriginInformationRaw) - 1), 1:5]
See(statusOriginInformation)
tail(statusOriginInformation)

rownames(statusOriginInformation) <- statusOriginInformation$Accession.ID
See(statusOriginInformation)
rownames(statusOriginInformation)[rownames(statusOriginInformation) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"

unique(statusOriginInformation$Status)
unique(statusOriginInformation$Origin)

SouthEastAsia <- c("Malaysia", "Taiwan", "Philippines", "Myanmar", "Viet Nam", "Thailand", "Indonesia", "Nepal", "Laos", "East Timor", "Cambodia")
originArea <- sapply(statusOriginInformation$Origin, function(x){
  if (x %in% SouthEastAsia){
    return(paste("South East Asia"))
  } else {
    return(x)
  }
}
)
See(originArea)
table(originArea)

statusOriginInformation <- cbind(statusOriginInformation, originArea)
See(statusOriginInformation)
# statusOriginInformation$originArea


markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47426527"
markerInterestID3 <- "Chr10_42562665"
markerInterestID4 <- "Chr17_16065902"

LineNames <- rownames(statusOriginInformation)
markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]
markerInterest3 <- gastonData0Matrix[, markerInterestID3]
markerInterest4 <- gastonData0Matrix[, markerInterestID4]

markerInterest1DF <- as.data.frame(markerInterest1)
markerInterest2DF <- as.data.frame(markerInterest2)
markerInterest3DF <- as.data.frame(markerInterest3)
markerInterest4DF <- as.data.frame(markerInterest4)

markerVarietyNames <- names(markerInterest1)
LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

statusOriginInformation <- statusOriginInformation[CommonNames, ]
See(statusOriginInformation)

markerInterest1DF <- markerInterest1DF[CommonNames, ]
markerInterest2DF <- markerInterest2DF[CommonNames, ]
markerInterest3DF <- markerInterest3DF[CommonNames, ]
markerInterest4DF <- markerInterest4DF[CommonNames, ]
See(markerInterest1DF)

markerValueCombination <- paste(markerInterest1DF, markerInterest2DF, markerInterest3DF, markerInterest4DF, sep = "_")
names(markerValueCombination) <- CommonNames
markerValueCombination <- as.data.frame(markerValueCombination)
See(markerValueCombination)

statusOriginMarkerCombinationInformation <- cbind(statusOriginInformation, markerValueCombination)
See(statusOriginInformation)
See(markerValueCombination)
See(statusOriginMarkerCombinationInformation, coln = 10)

# statusOriginMarkerCombinationInformation <- mutate(statusOriginMarkerCombinationInformation, group = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest2, sep = "_"))

# statusOriginMarkerCombinationInformation$markerInterest1 <- factor(statusOriginMarkerCombinationInformation$markerInterest1)
# statusOriginMarkerCombinationInformation$markerInterest2 <- factor(statusOriginMarkerCombinationInformation$markerInterest2)
# statusOriginMarkerCombinationInformation$markerInterest3 <- factor(statusOriginMarkerCombinationInformation$markerInterest3)
# statusOriginMarkerCombinationInformation$markerInterest4 <- factor(statusOriginMarkerCombinationInformation$markerInterest4)
statusOriginMarkerCombinationInformation$markerValueCombination <- factor(statusOriginMarkerCombinationInformation$markerValueCombination)
statusOriginMarkerCombinationInformation$Status <- factor(statusOriginMarkerCombinationInformation$Status)
statusOriginMarkerCombinationInformation$Origin <- factor(statusOriginMarkerCombinationInformation$Origin)
statusOriginMarkerCombinationInformation$originArea <- factor(statusOriginMarkerCombinationInformation$originArea)
See(statusOriginMarkerCombinationInformation)


a <- count(x = statusOriginMarkerCombinationInformation, markerValueCombination, Status)
b <- a %>% complete(markerValueCombination, Status, fill = list(n = 0))


pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_status_information/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4,".pdf"))
g <- ggplot(b, aes(x = markerValueCombination, y =  n, fill = Status)) + geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
plot(g)
dev.off()

pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_origin_area_information/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, ".pdf"))
g <- ggplot(statusOriginMarkerCombinationInformation, aes(x = markerValueCombination, y =  originArea, fill = originArea)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
plot(g)
dev.off()


g <- ggplot(b, aes(x = group, y = n, fill = Status)) + geom_bar(stat = "identity", position = "dodge") + labs(x = paste0(markerInterestID1, "_", markerInterestID2), y ="") + theme(text = element_text(size = 15))
plot(g)
dev.off()


### check
statusOriginMarkerCombinationInformation[statusOriginMarkerCombinationInformation$markerValueCombination == "2_2_0_0", ]
# 2_2_0_0 = "MisuzuDaizu", breeders line


### Check breeders lines
breedersLine <- statusOriginMarkerCombinationInformation[statusOriginMarkerCombinationInformation$Status == "Breeders line", ]
breedersLine020N <- breedersLine[breedersLine$markerValueCombination == "0_2_0_0" | breedersLine$markerValueCombination == "0_2_0_2", ]
breedersLine220N <- breedersLine[breedersLine$markerValueCombination == "2_2_0_0" | breedersLine$markerValueCombination == "2_2_0_2", ]
breedersLine000N <- breedersLine[breedersLine$markerValueCombination == "0_0_0_0" | breedersLine$markerValueCombination == "0_0_0_2", ]

breedersLine0000 <- breedersLine[breedersLine$markerValueCombination == "0_0_0_0", ]
breedersLine0002 <- breedersLine[breedersLine$markerValueCombination == "0_0_0_2", ]


###### 4. Lines with color phenotype and marker information ######
dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/"))

# Pubescence, Flower, Hypocotyl, Cotyldon and Mature leaf color
colorPhenotypeNames <- c("Pubescence_color", "Flower_color", "Hypocotyl_color", "Cotyldon_seed_color", "Leaf_color_at_maturity")
colorPhenotypeName <- "Pubescence_color"
markerInterestIDs <- c("Chr06_18760995","Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47426527"
markerInterestID3 <- "Chr10_42562665"
markerInterestID4 <- "Chr17_16065902"
for (colorPhenotypeName in colorPhenotypeNames){
  dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/"))
  dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/"))

  for (markerInterestID in markerInterestIDs){
    dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_", markerInterestID, "/"))
    dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_", markerInterestID, "/", scriptID, "_Number/"))
  }

  dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_four_markers/"))
  dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_four_markers/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/"))
  dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_four_markers/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/",scriptID, "_Number/"))


  colorPhenotypeRaw <- read.xlsx("raw_data/extra/Supplementary_Tables.xlsx", sheet = "Supplementary Table S5")
  See(colorPhenotypeRaw)
  colorPhenotypeRaw <- colorPhenotypeRaw[2:nrow(colorPhenotypeRaw), ]
  See(colorPhenotypeRaw)
  colnames(colorPhenotypeRaw) <- NULL
  colnames(colorPhenotypeRaw) <- colorPhenotypeRaw[1, ]
  colorPhenotypeRaw <- colorPhenotypeRaw[-1,]
  tail(colorPhenotypeRaw, 10)
  colorPhenotypeRaw <- colorPhenotypeRaw[1:197, ]
  rownames(colorPhenotypeRaw) <- colorPhenotypeRaw[, 1]
  colorPhenotypeRaw <- colorPhenotypeRaw[, -1]
  # finish reshaping
  See(colorPhenotypeRaw)



  LineNames <- rownames(colorPhenotypeRaw)
  markerInterest1 <- gastonData0Matrix[, markerInterestID1]
  markerInterest2 <- gastonData0Matrix[, markerInterestID2]
  markerInterest3 <- gastonData0Matrix[, markerInterestID3]
  markerInterest4 <- gastonData0Matrix[, markerInterestID4]

  # markerInterest1DF <- as.data.frame(markerInterest1)
  # markerInterest2DF <- as.data.frame(markerInterest2)
  # markerInterest3DF <- as.data.frame(markerInterest3)
  # markerInterest4DF <- as.data.frame(markerInterest4)

  markerVarietyNames <- names(markerInterest1)
  LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
  CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

  colorPhenotype <- colorPhenotypeRaw[CommonNames, ]
  See(colorPhenotype)


  markerInterest1 <- markerInterest1[CommonNames]
  markerInterest2 <- markerInterest2[CommonNames]
  markerInterest3 <- markerInterest3[CommonNames]
  markerInterest4 <- markerInterest4[CommonNames]
  See(markerInterest1)

  # markerInterest1DF <- markerInterest1DF[CommonNames, ]
  # markerInterest2DF <- markerInterest2DF[CommonNames, ]
  # markerInterest3DF <- markerInterest3DF[CommonNames, ]
  # markerInterest4DF <- markerInterest4DF[CommonNames, ]
  # See(markerInterest1DF)


  markerValueCombination <- paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_")
  names(markerValueCombination) <- CommonNames
  # markerValueCombination <- as.data.frame(markerValueCombination)
  See(markerValueCombination)

  # markerInterest1DF <- factor(markerInterest1DF)
  # markerInterest2DF <- factor(markerInterest2DF)
  # markerInterest3DF <- factor(markerInterest3DF)
  # markerInterest4DF <- factor(markerInterest4DF)

  # colorPhenotype$`Pubescence color` <- factor(colorPhenotype$`Pubescence color`)

  allInformation <- cbind(colorPhenotype, markerValueCombination)
  allInformation <- cbind(allInformation, markerInterest1)
  allInformation <- cbind(allInformation, markerInterest2)
  allInformation <- cbind(allInformation, markerInterest3)
  allInformation <- cbind(allInformation, markerInterest4)
  See(allInformation, coln = 30)

  # allInformation <- allInformation[!is.na(allInformation$`Pubescence color`), ]

  # markerValueCombination <- paste(markerInterest1DF, markerInterest2DF, markerInterest3DF, markerInterest4DF, sep = "_")
  # names(markerValueCombination) <- CommonNames
  # pubescenceLocusMarkers <- paste(markerValueCombination, colorPhenotype$T, sep = "_")
  # markerValueCombination <- as.data.frame(markerValueCombination)
  # See(markerValueCombination)

  # See(pubescenceLocusMarkers)
  # pubescenceColorAndLocusMarkers <- cbind(statusOriginInformation, markerValueCombination)

  # dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_Number/", scriptID, "_", markerInterestID1))
  # dir.create(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4))


  colnames(allInformation)[colnames(allInformation) == "Pubescence color"] <- "Pubescence_color"
  colnames(allInformation)[colnames(allInformation) == "Flower color"] <- "Flower_color"
  colnames(allInformation)[colnames(allInformation) == "Hypocotyl color"] <- "Hypocotyl_color"
  colnames(allInformation)[colnames(allInformation) == "Cotyledon seed color"] <- "Cotyldon_seed_color"
  colnames(allInformation)[colnames(allInformation) == "Leaf color at maturity"] <- "Leaf_color_at_maturity"
  table(allInformation$`Pubescence color`)
  See(allInformation, coln = 30)






# ### to edit
#   pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_", markerInterestID1, "/", scriptID, "_Number/", scriptID, "_", markerInterestID1, ".pdf"))
#   g <- ggplot(allInformation, aes(x = markerInterest1DF, y = colorPhenotypeName, fill = eval(parse(text = paste(colorPhenotypeName))))) + geom_bar(stat = "identity") + labs(x = markerInterestID1, fill = colorPhenotypeName)
#   plot(g)
#   dev.off()
#
#   pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_", markerInterestID2, "/", scriptID, "_Number/", scriptID, "_", markerInterestID2, ".pdf"))
#   g <- ggplot(allInformation, aes(x = markerInterest2DF, y = colorPhenotypeName, fill = eval(parse(text = paste(colorPhenotypeName))))) + geom_bar(stat = "identity")+ labs(x = markerInterestID2, fill = colorPhenotypeName)
#   plot(g)
#   dev.off()
#
#   pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_", markerInterestID3, "/", scriptID, "_Number/", scriptID, "_", markerInterestID3, ".pdf"))
#   g <- ggplot(allInformation, aes(x = markerInterest3DF, y = colorPhenotypeName, fill = eval(parse(text = paste(colorPhenotypeName))))) + geom_bar(stat = "identity")+ labs(x = markerInterestID3, fill = colorPhenotypeName)
#   plot(g)
#   dev.off()
#
#   pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_each_marker/", scriptID, "_", markerInterestID4, "/", scriptID, "_Number/", scriptID, "_", markerInterestID4, ".pdf"))
#   g <- ggplot(allInformation, aes(x = markerInterest4DF, y = eval(parse(text = paste(colorPhenotypeName))), fill = eval(parse(text = paste(colorPhenotypeName))))) + geom_bar(stat = "identity")+ labs(x = markerInterestID4, fill = colorPhenotypeName)
#   plot(g)
#   dev.off()
# ###







  a <- count(x = allInformation, markerValueCombination, !!as.name(colorPhenotypeName))
  colnames(a)[2] <- colorPhenotypeName
  b <- a %>% complete(markerValueCombination, !!as.name(colorPhenotypeName), fill = list(n = 0))

  pdf(paste0(dirMidSTAMColorPhenotypeAndMarkers, scriptID, "_Number_of_lines_with_color_phenotype/", scriptID, "_", colorPhenotypeName, "/", scriptID, "_For_four_markers/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/",  scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4,".pdf"))
  g <- ggplot(b, aes(x = markerValueCombination, y = n, fill = eval(parse(text = paste(colorPhenotypeName))))) + geom_bar(stat = "identity", position = "dodge") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))+ labs(fill = colorPhenotypeName)
  plot(g)
  dev.off()

}

### check, F3'H and Glyma.06G202300( marker:Chr06_18760995 )
allInformation[allInformation$Pubescence_color == "Gray" & allInformation$markerInterest1DF == 0, ]






