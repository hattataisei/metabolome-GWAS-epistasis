##########################################################################################
######  Title: 2.29_Soybean_STAM_Boxplot_of_amount_of_metabolites_with_high_heritability_flavonoid_for_different_Markers ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/25 (Created), 2025/04/08 (Last Updated)                       ######
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

scriptID <- "2.29"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites <- paste0(dirMidSTAMBase, scriptID, "_Boxplot_Amount_and_Line_Number_of_high_heritability_flavonoid/")
dir.create(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
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




###### 2. Read genotypic values of Flavonoid-related Metabolites and read marker genotype ######
##### 2.1. Read genotypic values of all metabolites in 2017 into R #####
gvMetab2017 <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
See(gvMetab2017)
rownames(gvMetab2017)[rownames(gvMetab2017) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"


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




###### 3. Boxplot of amount of Flavonoid-related Metabolites ######
##### 3.1. For each marker #####
#### Chr06_18760995(pcaMethods, nPC = 6, PC1), Chr10_42562665(pcaMethods, nPC = 6, PC4) etc.
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_each_marker/"))

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])



### boxplot
## alt code
markerInterestIDs <- c("Chr06_18760995", "Chr06_47486371", "Chr10_42562665", "Chr17_16065902")
# markerInterestID <- "Chr06_18760995"
# metabNo <- "X00015"
for (markerInterestID in markerInterestIDs){
  dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_each_marker/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID))
  dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_each_marker/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability>0.9/"))
  dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_each_marker/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability<0.9/"))

  for (metabNo in metabNamesFlavonoid){
    lineNames <- rownames(gvMetab2017)
    markerVarietyNames <- rownames(gastonData0Matrix)
    commonVarietyNames <- lineNames[(lineNames %in% markerVarietyNames)]

    markerInterest <- gastonData0Matrix[, markerInterestID]

    markerInterest <- markerInterest[commonVarietyNames]
    markerInterestDF <- as.data.frame(markerInterest)

    metabNow <- as.matrix(gvMetab2017)[commonVarietyNames, metabNo]

    markerInterestgvMetabNow <- cbind(markerInterestDF, metabNow)
    colnames(markerInterestgvMetabNow)[colnames(markerInterestgvMetabNow) == "metabNow"] <- metabNo
    markerInterestgvMetabNow$markerInterest <- factor(markerInterestgvMetabNow$markerInterest)

    # markerInterestgvMetabNowComplete <- markerInterestgvMetabNow %>%
    #   right_join(allMarkerValueComb, by = "markerValueCombination") %>%
    #   mutate(metabNo = ifelse(is.na(metabNo), NA, !!as.name(metabNo)))


    if (metabNo %in% metabNamesFlavonoidHeritabilityMoreThan0.9){
      pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_each_marker/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability>0.9/",  scriptID, "_", metabNo, ".pdf"))
      a <- ggplot(markerInterestgvMetabNow, aes(markerInterest, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
      plot(a)
      dev.off()

    } else {
      pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_each_marker/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability<0.9/", scriptID, "_", metabNo, ".pdf"))
      a <- ggplot(markerInterestgvMetabNow, aes(markerInterest, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
      plot(a)
      dev.off()

    }
  }
}



### violin plot
# dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_Violinplot/"))
#
# metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
# metabNamesFlavonoid <- metabFlavonoid[, "Name"]
# metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
# See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
# metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])
#
#
# # MetabNo <- "X200014"
# # markerInterestID <- "Chr06_18760995"
# markers <- c("Chr06_18760995", "Chr06_47486371", "Chr10_42562665", "Chr17_16065902")
# for (markerInterestID in markers){
#
#   dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID, "/"))
#   dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability<0.9/"))
#   dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability>0.9/"))
#   for (MetabNo in metabNamesFlavonoid){
#     LineNames <- rownames(gvMetab2017)
#     markerInterest <- gastonData0Matrix[, markerInterestID]
#
#     Metab <- as.data.frame(gvMetab2017[, MetabNo])
#
#     rownames(Metab) <- LineNames
#     colnames(Metab) <- MetabNo
#
#     markerVarietyNames <- names(markerInterest)
#     LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
#     CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]
#
#     Metab <- Metab[CommonNames, ]
#     Metab <- as.data.frame(Metab)
#     rownames(Metab) <- CommonNames
#
#     markerInterest <- as.data.frame(markerInterest)
#     colnames(markerInterest) <- markerInterestID
#     markerInterest <- markerInterest[CommonNames, ]
#     markerInterest <- as.data.frame(markerInterest)
#     rownames(markerInterest) <- CommonNames
#
#     markerInterest[, 2] <- CommonNames
#     markerInterest$markerInterest <- factor(markerInterest$markerInterest)
#     # See(markerInterest$markerInterest)
#     gvMetabMarkerInterest <- cbind(Metab, markerInterest)
#     # See(gvMetabMarkerInterest)
#
#     if (MetabNo %in% metabNamesFlavonoidHeritabilityMoreThan0.9){
#       pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability>0.9/", scriptID, "_", MetabNo, ".pdf"))
#
#       a <- ggplot(gvMetabMarkerInterest, aes(x = markerInterest, y = Metab, fill = markerInterest))+
#         geom_violin(alpha = 0.6)+
#         scale_fill_viridis(option = "D", discrete = T)+
#         geom_dotplot(binaxis = "y", stackdir = "down", dotsize = 0.3,
#                      position = position_nudge(-0.05))+
#         geom_boxplot(width = 0.1, color = "black", alpha=0.7,
#                      position = position_nudge(+0.05))+
#         theme_bw(base_size = 10)+
#         # labs(title="Sepal.Length")+
#         labs(x = NULL, y = NULL)+
#         theme(legend.position="none")
#       print(a)
#
#       dev.off()
#
#     } else {
#       pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_each_marker/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID, "/", scriptID, "_heritability<0.9/", scriptID, "_", MetabNo, ".pdf"))
#
#       b <- ggplot(gvMetabMarkerInterest, aes(x = markerInterest, y = Metab, fill = markerInterest))+
#         geom_violin(alpha = 0.6)+
#         scale_fill_viridis(option = "D", discrete = T)+
#         geom_dotplot(binaxis = "y", stackdir = "down", dotsize = 0.3,
#                      position = position_nudge(-0.05))+
#         geom_boxplot(width = 0.1, color = "black", alpha=0.7,
#                      position = position_nudge(+0.05))+
#         theme_bw(base_size = 10)+
#         # labs(title="Sepal.Length")+
#         labs(x = NULL, y = NULL)+
#         theme(legend.position="none")
#       print(b)
#
#       dev.off()
#
#     }
#   }
# }




# ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species))+
#   geom_violin(alpha = 0.6)+
#   scale_fill_viridis(option = "D", discrete = T)+
#   geom_dotplot(binaxis = "y", stackdir = "down", dotsize = 0.5,
#                position = position_nudge(-0.05))+
#   geom_boxplot(width = 0.01, color = "black", alpha=0.7,
#                position = position_nudge(+0.05))+
#   theme_bw(base_size = 10)+
#   labs(title="Sepal.Length")+
#   #labs(x = NULL, y = NULL)+
#   theme(legend.position="none")




##### 3.2. For combination of marker values #####
#### 3.2.1. For combination of marker values, (0,0),(0,2),(2,0),(2,2) ####



##### 3.3. For combination of two marker values, (0,0),(0,2),(2,0),(2,2) #####
#### 3.3.1. For "Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902" ####
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/"))

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/"))
metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])



markerInterestIDs <- c("Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# metabNo <- "X00015"
# markerInterestID1 <- "Chr06_18760995"
# markerInterestID2 <- "Chr06_47426527"


### boxplot
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/", scriptID, "_Boxplot/"))
for (markerInterestID1 in markerInterestIDs){
  markerInterestIDswithoutUntilMarkerInterestID1 <- markerInterestIDs[-(1:which(markerInterestID1 == markerInterestIDs))]

  for (markerInterestID2 in markerInterestIDswithoutUntilMarkerInterestID1){
    dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/"))
    dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability>0.9/"))
    dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability<0.9/"))

    for (metabNo in metabNamesFlavonoid){
      lineNames <- rownames(gvMetab2017)
      markerVarietyNames <- rownames(gastonData0Matrix)
      commonVarietyNames <- lineNames[(lineNames %in% markerVarietyNames)]

      markerInterest1 <- gastonData0Matrix[, markerInterestID1]
      markerInterest2 <- gastonData0Matrix[, markerInterestID2]

      markerInterest1 <- markerInterest1[commonVarietyNames]
      markerInterest2 <- markerInterest2[commonVarietyNames]

      markerInterest12Mat <- cbind(markerInterest1, markerInterest2)
      See(markerInterest12Mat)
      markerInterest12DF <- as.data.frame(markerInterest12Mat)

      markerInterest12DF <- mutate(.data = markerInterest12DF, markerValueCombination = paste(markerInterest1, markerInterest2, sep = "_"))

      metabNow <- as.matrix(gvMetab2017)[commonVarietyNames, metabNo]

      markerInterest12gvMetabNow <- cbind(markerInterest12DF, metabNow)
      colnames(markerInterest12gvMetabNow)[colnames(markerInterest12gvMetabNow) == "metabNow"] <- metabNo

      allMarkerValueComb <- expand.grid(
        markerInterest1 = c(0, 2),
        markerInterest2 = c(0, 2)
      ) %>%
        mutate(markerValueCombination = paste(markerInterest1, markerInterest2, sep = "_"))

      markerInterest12gvMetabNowComplete <- markerInterest12gvMetabNow %>%
        right_join(allMarkerValueComb, by = "markerValueCombination") %>%
        mutate(metabNo = ifelse(is.na(metabNo), NA, !!as.name(metabNo)))


      if (metabNo %in% metabNamesFlavonoidHeritabilityMoreThan0.9){
        pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability>0.9/",  scriptID, "_", metabNo, ".pdf"))
        a <- ggplot(markerInterest12gvMetabNowComplete, aes(markerValueCombination, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
        plot(a)
        dev.off()

      } else {
        pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_two_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability<0.9/", scriptID, "_", metabNo, ".pdf"))
        a <- ggplot(markerInterest12gvMetabNowComplete, aes(markerValueCombination, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
        plot(a)
        dev.off()

      }
    }
  }
}



### violin plot
# # MetabNo <- "X00015"
# # markerInterestID1 <- "Chr06_18760995"
# # markerInterestID2 <- "Chr06_47426527"
# dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/"))
# dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/", scriptID, "_Violinplot/"))
# for (markerInterestID1 in markerInterestIDs){
#   markerInterestIDswithoutUntilMarkerInterestID1 <- markerInterestIDs[-(1:which(markerInterestID1 == markerInterestIDs))]
#
#   for (markerInterestID2 in markerInterestIDswithoutUntilMarkerInterestID1){
#     dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/"))
#     dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability>0.9/"))
#     dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/", scriptID, "_Violinplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability<0.9/"))
#
#
#     for (MetabNo in metabNamesFlavonoid){
#       LineNames <- rownames(gvMetab2017)
#       markerInterest1 <- gastonData0Matrix[, markerInterestID1]
#       markerInterest2 <- gastonData0Matrix[, markerInterestID2]
#
#       markerInterest1DF <- as.data.frame(markerInterest1)
#       markerInterest2DF <- as.data.frame(markerInterest2)
#
#       markerVarietyNames <- names(markerInterest1)
#       LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
#       CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]
#
#       markerInterest1DF <- markerInterest1DF[CommonNames, ]
#       names(markerInterest1DF) <- CommonNames
#       markerInterest1DF <- as.data.frame(markerInterest1DF)
#       markerInterest2DF <- markerInterest2DF[CommonNames, ]
#       names(markerInterest2DF) <- CommonNames
#       markerInterest2DF <- as.data.frame(markerInterest2DF)
#
#       markerInterest1And2DF <- cbind(markerInterest1DF, markerInterest2DF)
#       # See(markerInterest1And2DF)
#
#       markerInterest1And2DF <- markerInterest1And2DF %>%
#         unite(
#           col = "markerInterest1And2",
#           c("markerInterest1DF", "markerInterest2DF"),
#           sep = " ",
#           remove = TRUE,
#           na.rm = TRUE
#         )
#
#       Metab <- as.data.frame(gvMetab2017[, MetabNo])
#       rownames(Metab) <- rownames(gvMetab2017)
#       colnames(Metab) <- MetabNo
#       Metab <- Metab[CommonNames, ]
#       Metab <- as.data.frame(Metab)
#       rownames(Metab) <- CommonNames
#       See(Metab)
#
#       MetabAndMarkerInterest1And2DF <- cbind(Metab, markerInterest1And2DF)
#       See(MetabAndMarkerInterest1And2DF)
#
#       if (MetabNo %in% metabNamesFlavonoidHeritabilityMoreThan0.9){
#         pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/", scriptID, "_Violinplot/",  scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability>0.9/", scriptID, "_", MetabNo, ".pdf"))
#
#         a <- ggplot(MetabAndMarkerInterest1And2DF, aes(x = markerInterest1And2, y = Metab, fill = markerInterest1And2))+
#           geom_violin(alpha = 0.6)+
#           scale_fill_viridis(option = "D", discrete = T)+
#           geom_dotplot(binaxis = "y", stackdir = "down", dotsize = 0.1,
#                        position = position_nudge(-0.05))+
#           geom_boxplot(width = 0.1, color = "black", alpha=0.7,
#                        position = position_nudge(+0.05))+
#           theme_bw(base_size = 10)+
#           # labs(title="Sepal.Length")+
#           labs(x = NULL, y = NULL)+
#           theme(legend.position="none")
#         print(a)
#
#         dev.off()
#
#       } else {
#         pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_two_markers/", scriptID, "_Violinplot/",  scriptID, "_", markerInterestID1, "_", markerInterestID2, "/", scriptID, "_heritability<0.9/", scriptID, "_", MetabNo, ".pdf"))
#
#         b <- ggplot(MetabAndMarkerInterest1And2DF, aes(x = markerInterest1And2, y = Metab, fill = markerInterest1And2))+
#           geom_violin(alpha = 0.6)+
#           scale_fill_viridis(option = "D", discrete = T)+
#           geom_dotplot(binaxis = "y", stackdir = "down", dotsize = 0.1,
#                        position = position_nudge(-0.05))+
#           geom_boxplot(width = 0.1, color = "black", alpha=0.7,
#                        position = position_nudge(+0.05))+
#           theme_bw(base_size = 10)+
#           # labs(title="Sepal.Length")+
#           labs(x = NULL, y = NULL)+
#           theme(legend.position="none")
#         print(b)
#
#         dev.off()
#
#       }
#     }
#   }
# }





##### 3.4. For combination of three marker values, (0,0,0),... #####
#### 3.4.1. For "Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902" ####
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/"))

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_For_three_markers/"))
metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])



markerInterestIDs <- c("Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# metabNo <- "X00015"
# markerInterestID1 <- "Chr06_18760995"
# markerInterestID2 <- "Chr06_47426527"
# markerInterestID3 <- "Chr10_42562665"


### boxplot
## alt code
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/", scriptID, "_Boxplot/"))

for (markerInterestID1 in markerInterestIDs){
  markerInterestIDswithoutUntilMarkerInterestID1 <- markerInterestIDs[-(1:which(markerInterestID1 == markerInterestIDs))]

  for (markerInterestID2 in markerInterestIDswithoutUntilMarkerInterestID1){
    markerInterestIDswithoutUntilMarkerInterestID1And2 <- markerInterestIDs[-(1:which(markerInterestID2 == markerInterestIDs))]
    for (markerInterestID3 in markerInterestIDswithoutUntilMarkerInterestID1And2){

      dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "/"))
      dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "/", scriptID, "_heritability>0.9/"))
      dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "/", scriptID, "_heritability<0.9/"))

      for (metabNo in metabNamesFlavonoid){
        lineNames <- rownames(gvMetab2017)
        markerVarietyNames <- rownames(gastonData0Matrix)
        commonVarietyNames <- lineNames[(lineNames %in% markerVarietyNames)]

        markerInterest1 <- gastonData0Matrix[, markerInterestID1]
        markerInterest2 <- gastonData0Matrix[, markerInterestID2]
        markerInterest3 <- gastonData0Matrix[, markerInterestID3]

        markerInterest1 <- markerInterest1[commonVarietyNames]
        markerInterest2 <- markerInterest2[commonVarietyNames]
        markerInterest3 <- markerInterest3[commonVarietyNames]

        markerInterest123Mat <- cbind(markerInterest1, markerInterest2, markerInterest3)
        See(markerInterest123Mat)
        markerInterest123DF <- as.data.frame(markerInterest123Mat)

        markerInterest123DF <- mutate(.data = markerInterest123DF, markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, sep = "_"))

        metabNow <- as.matrix(gvMetab2017)[commonVarietyNames, metabNo]

        markerInterest123gvMetabNow <- cbind(markerInterest123DF, metabNow)
        colnames(markerInterest123gvMetabNow)[colnames(markerInterest123gvMetabNow) == "metabNow"] <- metabNo

        allMarkerValueComb <- expand.grid(
          markerInterest1 = c(0, 2),
          markerInterest2 = c(0, 2),
          markerInterest3 = c(0, 2)
        ) %>%
          mutate(markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, sep = "_"))

        markerInterest123gvMetabNowComplete <- markerInterest123gvMetabNow %>%
          right_join(allMarkerValueComb, by = "markerValueCombination") %>%
          mutate(metabNo = ifelse(is.na(metabNo), NA, !!as.name(metabNo)))


        if (metabNo %in% metabNamesFlavonoidHeritabilityMoreThan0.9){
          pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "/", scriptID, "_heritability>0.9/",  scriptID, "_", metabNo, ".pdf"))
          a <- ggplot(markerInterest123gvMetabNowComplete, aes(markerValueCombination, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
          plot(a)
          dev.off()

        } else {
          pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_three_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "/", scriptID, "_heritability<0.9/", scriptID, "_", metabNo, ".pdf"))
          a <- ggplot(markerInterest123gvMetabNowComplete, aes(markerValueCombination, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
          plot(a)
          dev.off()

        }
      }
    }
  }
}




##### 3.5. For combination of four marker values, (0,0,0,0),... #####
#### 3.5.1. For "Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902" ####

### boxplot
## alt code (using ggplot)
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/", scriptID, "_Boxplot/"))

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_marker_value_combination/"))
metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metabFlavonoidHeritabilityMoreThan0.9, coln = 12)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metabFlavonoidHeritabilityMoreThan0.9[, 11:ncol(metabFlavonoidHeritabilityMoreThan0.9)])

# markerInterestIDs <- c("Chr06_18760995", "Chr10_42562665", "Chr06_47426527", "Chr17_16065902")
# MetabNo <- "X00015"
markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47426527"
markerInterestID3 <- "Chr10_42562665"
markerInterestID4 <- "Chr17_16065902"

dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/", scriptID, "_heritability>0.9/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/", scriptID, "_heritability<0.9/"))

# metabNo <- "X00015"
for (metabNo in metabNamesFlavonoid){
  lineNames <- rownames(gvMetab2017)
  markerVarietyNames <- rownames(gastonData0Matrix)
  commonVarietyNames <- lineNames[(lineNames %in% markerVarietyNames)]

  markerInterest1 <- gastonData0Matrix[, markerInterestID1]
  markerInterest2 <- gastonData0Matrix[, markerInterestID2]
  markerInterest3 <- gastonData0Matrix[, markerInterestID3]
  markerInterest4 <- gastonData0Matrix[, markerInterestID4]

  markerInterest1 <- markerInterest1[commonVarietyNames]
  markerInterest2 <- markerInterest2[commonVarietyNames]
  markerInterest3 <- markerInterest3[commonVarietyNames]
  markerInterest4 <- markerInterest4[commonVarietyNames]

  markerInterest1234Mat <- cbind(markerInterest1, markerInterest2, markerInterest3, markerInterest4)
  See(markerInterest1234Mat)
  markerInterest1234DF <- as.data.frame(markerInterest1234Mat)

  markerInterest1234DF <- mutate(.data = markerInterest1234DF, markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_"))

  metabNow <- as.matrix(gvMetab2017)[commonVarietyNames, metabNo]

  markerInterest1234gvMetabNow <- cbind(markerInterest1234DF, metabNow)
  colnames(markerInterest1234gvMetabNow)[colnames(markerInterest1234gvMetabNow) == "metabNow"] <- metabNo

  allMarkerValueComb <- expand.grid(
    markerInterest1 = c(0, 2),
    markerInterest2 = c(0, 2),
    markerInterest3 = c(0, 2),
    markerInterest4 = c(0, 2)
  ) %>%
    mutate(markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_"))

  markerInterest1234gvMetabNowComplete <- markerInterest1234gvMetabNow %>%
    right_join(allMarkerValueComb, by = "markerValueCombination") %>%
    mutate(metabNo = ifelse(is.na(metabNo), NA, !!as.name(metabNo)))


  if (metabNo %in% metabNamesFlavonoidHeritabilityMoreThan0.9){
    pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/", scriptID, "_heritability>0.9/",  scriptID, "_", metabNo, ".pdf"))
    a <- ggplot(markerInterest1234gvMetabNowComplete, aes(markerValueCombination, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
    plot(a)
    dev.off()

  } else {
    pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Metabolite_accumulation/", scriptID, "_For_four_markers/", scriptID, "_Boxplot/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "/", scriptID, "_heritability<0.9/", scriptID, "_", metabNo, ".pdf"))
    a <- ggplot(markerInterest1234gvMetabNowComplete, aes(markerValueCombination, !!as.name(metabNo))) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15), axis.text.y = element_text(size = 15)) + labs(x = NULL, y = NULL)
    plot(a)
    dev.off()

  }
}




#### 3.5.2. For "Chr06_18760995", "Chr06_47426527", "Chr15_45643039", "Chr18_6592763" ####



###### 4. Number of Lines with group information for flavonoid-related markers ######
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/"))


### Check "Houjaku Kuwazu"
groupInfo <- read.csv("data/extra/0.3_group_information.csv", row.names = 1)
groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"
See(groupInfo)

rownames(groupInfo) <- groupInfo$lineNames
See(groupInfo)



##### 4.1. For each marker #####
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_each_marker/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_each_marker/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_each_marker/", scriptID, "_Ratio/"))


### alt code
markerInterestIDs <- c("Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# markerInterestID <- "Chr06_18760995"
for (markerInterestID in markerInterestIDs){

  lineNames <- rownames(gvMetab2017)
  markerVarietyNames <- rownames(gastonData0Matrix)
  commonVarietyNames <- lineNames[(lineNames %in% markerVarietyNames)]

  markerInterest <- gastonData0Matrix[, markerInterestID]
  markerInterest <- markerInterest[commonVarietyNames]

  groupInfo <- groupInfo[commonVarietyNames, ]

  group <- groupInfo$group

  groupMarkerInterest <- cbind(group, markerInterest)
  groupMarkerInterestDF <- as.data.frame(groupMarkerInterest)

  groupMarkerInterestDF$group <- as.factor(groupMarkerInterestDF$group)
  groupMarkerInterestDF$markerInterest <- as.factor(groupMarkerInterestDF$markerInterest)


  a <- count(x = groupMarkerInterestDF, markerInterest, group)
  b <- a %>% complete(markerInterest, group, fill = list(n = 0))

  pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_each_marker/", scriptID, "_Number/", scriptID,  "_", markerInterestID, ".pdf"))
  g <- ggplot(b, aes(x = markerInterest, y = n, fill = group)) + geom_bar(stat = "identity", position = "dodge") + labs(x = paste0(markerInterestID), y ="") + theme(text = element_text(size = 15))
  plot(g)
  dev.off()


  a <- groupMarkerInterestDF %>%
    count(markerInterest, group) %>%
    group_by(markerInterest) %>%
    mutate(percentage = n / sum(n) * 100)

  pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_each_marker/", scriptID, "_Ratio/", scriptID,  "_", markerInterestID, ".pdf"))
  g <- ggplot(a, aes(x = markerInterest, y = percentage, fill = group)) + geom_bar(stat = "identity") + labs(x = paste0(markerInterestID), y ="") + theme(text = element_text(size = 15))
  plot(g)
  dev.off()

}



##### 4.2. For two markers combination #####
# to edit

dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_number_of_lines_with_group_information/", scriptID, "_For_two_markers/"))


#### 4.2.1. Chr06_18760995, Chr10_42562665, Chr06_47426527, Chr17_16065902 ####
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_two_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_two_markers/", scriptID, "_Ratio/"))

markerInterestIDs <- c("Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# markerInterestID1 <- "Chr06_47426527"
# markerInterestID2 <- "Chr17_16065902"
for (markerInterestID1 in markerInterestIDs){
  markerInterestIDswithoutUntilMarkerInterestID1 <- markerInterestIDs[-(1:which(markerInterestID1 == markerInterestIDs))]

  for (markerInterestID2 in markerInterestIDswithoutUntilMarkerInterestID1){

    lineNames <- rownames(gvMetab2017)
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
    # See(markerInterest1And2DF)

    marker0_0 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1DF == 0 & markerInterest1And2DF$markerInterest2DF == 0, ]
    marker0_2 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1DF == 0 & markerInterest1And2DF$markerInterest2DF == 2, ]
    marker2_0 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1DF == 2 & markerInterest1And2DF$markerInterest2DF == 0, ]
    marker2_2 <- markerInterest1And2DF[markerInterest1And2DF$markerInterest1DF == 2 & markerInterest1And2DF$markerInterest2DF == 2, ]
    # See(marker0_0)
    # See(marker0_2)
    # See(marker2_0)
    # See(marker2_2)
    varietyNames0_0 <- rownames(marker0_0)
    varietyNames0_2 <- rownames(marker0_2)
    varietyNames2_0 <- rownames(marker2_0)
    varietyNames2_2 <- rownames(marker2_2)

    groupInfo0_0 <- groupInfo[groupInfo$lineNames %in% varietyNames0_0, ]
    groupInfo0_2 <- groupInfo[groupInfo$lineNames %in% varietyNames0_2, ]
    groupInfo2_0 <- groupInfo[groupInfo$lineNames %in% varietyNames2_0, ]
    groupInfo2_2 <- groupInfo[groupInfo$lineNames %in% varietyNames2_2, ]
    table(groupInfo0_0$group)
    table(groupInfo0_2$group)
    table(groupInfo2_0$group)
    table(groupInfo2_2$group)


    groupMat <- bind_rows(table(groupInfo0_0$group), table(groupInfo0_2$group))
    groupMat <- bind_rows(groupMat, table(groupInfo2_0$group))
    groupMat <- bind_rows(groupMat, table(groupInfo2_2$group))
    groupMat <- as.data.frame(groupMat)
    groupMat[is.na(groupMat)] <- 0
    rownames(groupMat) <- c("0_0", "0_2", "2_0", "2_2")
    groupNames <- c("Japan", "Primitive", "World")
    groupMat <- groupMat[, groupNames]

    pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_two_markers/", scriptID, "_Number/", scriptID,  "_", markerInterestID1, "_", markerInterestID2, ".pdf"))
    barplot(t(groupMat),
            las = 2,
            legend.text = groupNames,
            main = paste0(markerInterestID1, "_", markerInterestID2),
            cex.names = 1.4,
            cex.axis = 1.4)
    dev.off()

    groupMat <- NA
    groupMat <- bind_rows(table(groupInfo0_0$group) / sum(table(groupInfo0_0$group)), table(groupInfo0_2$group) / sum(table(groupInfo0_2$group)))
    groupMat <- bind_rows(groupMat, table(groupInfo2_0$group) / sum(table(groupInfo2_0$group)))
    groupMat <- bind_rows(groupMat, table(groupInfo2_2$group) / sum(table(groupInfo2_2$group)))
    groupMat <- as.data.frame(groupMat)
    groupMat[is.na(groupMat)] <- 0
    rownames(groupMat) <- c("0_0", "0_2", "2_0", "2_2")
    groupMat <- groupMat[, groupNames]
    groupNames <- c("Japan", "Primitive", "World")

    pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_two_markers/", scriptID, "_Ratio/", scriptID,  "_", markerInterestID1, "_", markerInterestID2, ".pdf"))
    barplot(t(groupMat),
            las = 2,
            legend.text = groupNames,
            main = paste0(markerInterestID1, "_", markerInterestID2),
            cex.names = 1.4,
            cex.axis = 1.4)
    dev.off()


  }
}



##### 4.3 For three markers combination #####
# to edit

dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_number_of_lines_with_group_information/", scriptID, "_For_three_markers/"))

dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_three_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_three_markers/", scriptID, "_Ratio/"))

markerInterestIDs <- c("Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902")
# markerInterestID1 <- "Chr06_47426527"
# markerInterestID2 <- "Chr17_16065902"
for (markerInterestID1 in markerInterestIDs){
  markerInterestIDswithoutUntilMarkerInterestID1 <- markerInterestIDs[-(1:which(markerInterestID1 == markerInterestIDs))]

  for (markerInterestID2 in markerInterestIDswithoutUntilMarkerInterestID1){
    markerInterestIDswithoutUntilMarkerInterestID2 <- markerInterestIDs[-(1:which(markerInterestID2 == markerInterestIDs))]

    for (markerInterestID3 in markerInterestIDswithoutUntilMarkerInterestID2){

      LineNames <- rownames(gvMetab2017)
      markerInterest1 <- gastonData0Matrix[, markerInterestID1]
      markerInterest2 <- gastonData0Matrix[, markerInterestID2]
      markerInterest3 <- gastonData0Matrix[, markerInterestID3]

      markerInterest1DF <- as.data.frame(markerInterest1)
      markerInterest2DF <- as.data.frame(markerInterest2)
      markerInterest3DF <- as.data.frame(markerInterest3)

      markerVarietyNames <- names(markerInterest1)
      LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
      CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]

      markerInterest1DF <- markerInterest1DF[CommonNames, ]
      names(markerInterest1DF) <- CommonNames
      markerInterest1DF <- as.data.frame(markerInterest1DF)
      markerInterest2DF <- markerInterest2DF[CommonNames, ]
      names(markerInterest2DF) <- CommonNames
      markerInterest2DF <- as.data.frame(markerInterest2DF)
      markerInterest3DF <- markerInterest3DF[CommonNames, ]
      names(markerInterest3DF) <- CommonNames
      markerInterest3DF <- as.data.frame(markerInterest3DF)

      markerInterest1And2And3DF <- cbind(markerInterest1DF, markerInterest2DF, markerInterest3DF)
      # See(markerInterest1And2DF)

      marker0_0_0 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 0 & markerInterest1And2And3DF$markerInterest2DF == 0 & markerInterest1And2And3DF$markerInterest3DF == 0, ]
      marker0_0_2 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 0 & markerInterest1And2And3DF$markerInterest2DF == 0 & markerInterest1And2And3DF$markerInterest3DF == 2, ]
      marker0_2_0 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 0 & markerInterest1And2And3DF$markerInterest2DF == 2 & markerInterest1And2And3DF$markerInterest3DF == 0, ]
      marker0_2_2 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 0 & markerInterest1And2And3DF$markerInterest2DF == 2 & markerInterest1And2And3DF$markerInterest3DF == 2, ]
      marker2_0_0 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 2 & markerInterest1And2And3DF$markerInterest2DF == 0 & markerInterest1And2And3DF$markerInterest3DF == 0, ]
      marker2_0_2 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 2 & markerInterest1And2And3DF$markerInterest2DF == 0 & markerInterest1And2And3DF$markerInterest3DF == 2, ]
      marker2_2_0 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 2 & markerInterest1And2And3DF$markerInterest2DF == 2 & markerInterest1And2And3DF$markerInterest3DF == 0, ]
      marker2_2_2 <- markerInterest1And2And3DF[markerInterest1And2And3DF$markerInterest1DF == 2 & markerInterest1And2And3DF$markerInterest2DF == 2 & markerInterest1And2And3DF$markerInterest3DF == 2, ]


      varietyNames0_0_0 <- rownames(marker0_0_0)
      varietyNames0_0_2 <- rownames(marker0_0_2)
      varietyNames0_2_0 <- rownames(marker0_2_0)
      varietyNames0_2_2 <- rownames(marker0_2_2)
      varietyNames2_0_0 <- rownames(marker2_0_0)
      varietyNames2_0_2 <- rownames(marker2_0_2)
      varietyNames2_2_0 <- rownames(marker2_2_0)
      varietyNames2_2_2 <- rownames(marker2_2_2)

      groupInfo0_0_0 <- groupInfo[groupInfo$lineNames %in% varietyNames0_0_0, ]
      groupInfo0_0_2 <- groupInfo[groupInfo$lineNames %in% varietyNames0_0_2, ]
      groupInfo0_2_0 <- groupInfo[groupInfo$lineNames %in% varietyNames0_2_0, ]
      groupInfo0_2_2 <- groupInfo[groupInfo$lineNames %in% varietyNames0_2_2, ]
      groupInfo2_0_0 <- groupInfo[groupInfo$lineNames %in% varietyNames2_0_0, ]
      groupInfo2_0_2 <- groupInfo[groupInfo$lineNames %in% varietyNames2_0_2, ]
      groupInfo2_2_0 <- groupInfo[groupInfo$lineNames %in% varietyNames2_2_0, ]
      groupInfo2_2_2 <- groupInfo[groupInfo$lineNames %in% varietyNames2_2_2, ]
      table(groupInfo0_0_0$group)
      table(groupInfo0_0_2$group)
      table(groupInfo0_2_0$group)
      table(groupInfo0_2_2$group)
      table(groupInfo2_0_0$group)
      table(groupInfo2_0_2$group)
      table(groupInfo2_2_0$group)
      table(groupInfo2_2_2$group)


      groupMat <- bind_rows(table(groupInfo0_0_0$group), table(groupInfo0_0_2$group))
      groupMat <- bind_rows(groupMat, table(groupInfo0_2_0$group))
      groupMat <- bind_rows(groupMat, table(groupInfo0_2_2$group))
      groupMat <- bind_rows(groupMat, table(groupInfo2_0_0$group))
      groupMat <- bind_rows(groupMat, table(groupInfo2_0_2$group))
      groupMat <- bind_rows(groupMat, table(groupInfo2_2_0$group))
      groupMat <- bind_rows(groupMat, table(groupInfo2_2_2$group))
      groupMat <- as.data.frame(groupMat)
      groupMat[is.na(groupMat)] <- 0

      markerBarNames <- c("0_0_0", "0_0_2", "0_2_0", "0_2_2", "2_0_0", "2_0_2", "2_2_0", "2_2_2")
      rownames(groupMat) <- markerBarNames
      groupNames <- c("Japan", "Primitive", "World")
      groupMat <- groupMat[, groupNames]

      pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_three_markers/", scriptID, "_Number/", scriptID,  "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, ".pdf"))
      barplot(t(groupMat),
              las = 2,
              legend.text = groupNames,
              main = paste0(markerInterestID1, "_", markerInterestID2, "_", markerInterestID3),
              cex.names = 1.2
              )
      dev.off()

      groupMat <- NA
      groupMat <- bind_rows(table(groupInfo0_0_0$group) / sum(table(groupInfo0_0_0$group)), table(groupInfo0_0_2$group) / sum(table(groupInfo0_0_2$group)))
      groupMat <- bind_rows(groupMat, table(groupInfo0_2_0$group) / sum(table(groupInfo0_2_0$group)))
      groupMat <- bind_rows(groupMat, table(groupInfo0_2_2$group) / sum(table(groupInfo0_2_2$group)))
      groupMat <- bind_rows(groupMat, table(groupInfo2_0_0$group) / sum(table(groupInfo2_0_0$group)))
      groupMat <- bind_rows(groupMat, table(groupInfo2_0_2$group) / sum(table(groupInfo2_0_2$group)))
      groupMat <- bind_rows(groupMat, table(groupInfo2_2_0$group) / sum(table(groupInfo2_2_0$group)))
      groupMat <- bind_rows(groupMat, table(groupInfo2_2_2$group) / sum(table(groupInfo2_2_2$group)))
      groupMat <- as.data.frame(groupMat)

      groupMat[is.na(groupMat)] <- 0
      rownames(groupMat) <- markerBarNames
      groupMat <- groupMat[, groupNames]
      groupNames <- c("Japan", "Primitive", "World")

      pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_three_markers/", scriptID, "_Ratio/", scriptID,  "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, ".pdf"))
      barplot(t(groupMat), las = 2,
              legend.text = groupNames,
              main = paste0(markerInterestID1, "_", markerInterestID2, "_", markerInterestID3),
              cex.names = 1.2
              )
      dev.off()

    }
  }
}






##### 4.4 For four markers combination #####
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_number_of_lines_with_group_information/", scriptID, "_For_four_markers/"))

dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_four_markers/", scriptID, "_Number/"))
dir.create(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_four_markers/", scriptID, "_Ratio/"))



#### 4.4.1 "Chr06_18760995", "Chr06_47426527", "Chr10_42562665", "Chr17_16065902" ####

### alt code
### Read "group data", and check "Houjaku Kuwazu"
groupInfo <- read.csv("data/extra/0.3_group_information.csv", row.names = 1)
groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"
See(groupInfo)

rownames(groupInfo) <- groupInfo$lineNames
See(groupInfo)


# genome, 198 varieties
markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47426527"
markerInterestID3 <- "Chr10_42562665"
markerInterestID4 <- "Chr17_16065902"

markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]
markerInterest3 <- gastonData0Matrix[, markerInterestID3]
markerInterest4 <- gastonData0Matrix[, markerInterestID4]

groupVarietyNames <- rownames(groupInfo)
markerVarietyNames <- rownames(gastonData0Matrix)
commonVarietyNames <- markerVarietyNames[(markerVarietyNames %in% groupVarietyNames)]

markerInterest1 <- markerInterest1[commonVarietyNames]
markerInterest2 <- markerInterest2[commonVarietyNames]
markerInterest3 <- markerInterest3[commonVarietyNames]
markerInterest4 <- markerInterest4[commonVarietyNames]


markerInterest1234Mat <- cbind(markerInterest1, markerInterest2, markerInterest3, markerInterest4)
See(markerInterest1234Mat)

markerInterest1234DF <- as.data.frame(markerInterest1234Mat)
markerInterest1234DF <- mutate(.data = markerInterest1234DF, markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_"))
See(markerInterest1234DF)
length(unique(markerInterest1234DF$markerValueCombination))


allMarkerValueCombinations <- expand.grid(
  markerInterest1 = c(0, 2),
  markerInterest2 = c(0, 2),
  markerInterest3 = c(0, 2),
  markerInterest4 = c(0, 2)
  ) %>%
  mutate(markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_"))
See(allMarkerValueCombinations)

groupInfo <- groupInfo[commonVarietyNames, ]
group <- groupInfo$group
markerInterest1234GroupDF <- cbind(markerInterest1234DF, group)
See(markerInterest1234GroupDF)


allMarkerValueGroupCombinations <- expand.grid(
  markerValueCombination = unique(allMarkerValueCombinations$markerValueCombination),
  group = unique(group)
)


markerInterest1234GroupDFCount <- markerInterest1234GroupDF %>%
  count(markerValueCombination, group)

markerInterest1234GroupDFCompleteCombination <- full_join(markerInterest1234GroupDFCount, allMarkerValueGroupCombinations, by = c("markerValueCombination", "group")) %>%
  replace_na(list(n = 0))
markerInterest1234GroupDFCompleteCombination
See(markerInterest1234GroupDFCompleteCombination)


markerInterest1234GroupDFCompleteCombination <- markerInterest1234GroupDFCompleteCombination %>%
  mutate(group = factor(group, levels = c("Japan", "World", "Primitive")))


pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "_genome_198_varieties_dodge", ".pdf"))
g <- ggplot(markerInterest1234GroupDFCompleteCombination, aes(x = markerValueCombination, y = n, fill = group)) + geom_bar(stat = "identity", position = "dodge") + labs(x = paste0(markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4), y ="") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot(g)
dev.off()

pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "_genome_198_varieties", ".pdf"))
g <- ggplot(markerInterest1234GroupDFCompleteCombination, aes(x = markerValueCombination, y = n, fill = group)) + geom_bar(stat = "identity") + labs(x = paste0(markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4), y ="") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot(g)
dev.off()



# genome and metabolome, 195 varieties
markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47426527"
markerInterestID3 <- "Chr10_42562665"
markerInterestID4 <- "Chr17_16065902"

markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]
markerInterest3 <- gastonData0Matrix[, markerInterestID3]
markerInterest4 <- gastonData0Matrix[, markerInterestID4]

metabVarietyNames <- rownames(gvMetab2017)
markerVarietyNames <- rownames(gastonData0Matrix)
commonVarietyNames <- markerVarietyNames[(markerVarietyNames %in% metabVarietyNames)]

markerInterest1 <- markerInterest1[commonVarietyNames]
markerInterest2 <- markerInterest2[commonVarietyNames]
markerInterest3 <- markerInterest3[commonVarietyNames]
markerInterest4 <- markerInterest4[commonVarietyNames]


markerInterest1234Mat <- cbind(markerInterest1, markerInterest2, markerInterest3, markerInterest4)
See(markerInterest1234Mat)

markerInterest1234DF <- as.data.frame(markerInterest1234Mat)
markerInterest1234DF <- mutate(.data = markerInterest1234DF, markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_"))
See(markerInterest1234DF)
length(unique(markerInterest1234DF$markerValueCombination))


allMarkerValueCombinations <- expand.grid(
  markerInterest1 = c(0, 2),
  markerInterest2 = c(0, 2),
  markerInterest3 = c(0, 2),
  markerInterest4 = c(0, 2)
) %>%
  mutate(markerValueCombination = paste(markerInterest1, markerInterest2, markerInterest3, markerInterest4, sep = "_"))
See(allMarkerValueCombinations)

groupInfo <- groupInfo[commonVarietyNames, ]
group <- groupInfo$group
markerInterest1234GroupDF <- cbind(markerInterest1234DF, group)
See(markerInterest1234GroupDF)


allMarkerValueGroupCombinations <- expand.grid(
  markerValueCombination = unique(allMarkerValueCombinations$markerValueCombination),
  group = unique(group)
)


markerInterest1234GroupDFCount <- markerInterest1234GroupDF %>%
  count(markerValueCombination, group)

markerInterest1234GroupDFCompleteCombination <- full_join(markerInterest1234GroupDFCount, allMarkerValueGroupCombinations, by = c("markerValueCombination", "group")) %>%
  replace_na(list(n = 0))
markerInterest1234GroupDFCompleteCombination
See(markerInterest1234GroupDFCompleteCombination)


markerInterest1234GroupDFCompleteCombination <- markerInterest1234GroupDFCompleteCombination %>%
  mutate(group = factor(group, levels = c("Japan", "World", "Primitive")))

pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "_genome_and_metabolome_195_varieties_dodge", ".pdf"))
g <- ggplot(markerInterest1234GroupDFCompleteCombination, aes(x = markerValueCombination, y = n, fill = group)) + geom_bar(stat = "identity", position = "dodge") + labs(x = paste0(markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4), y ="") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot(g)
dev.off()

pdf(paste0(dirMidSTAMBoxplotOfFlavonoidMetabolitesAmountForDifferentMarkersHighHeritabilityFlavonoidMetabolites, scriptID, "_Number_of_lines_with_group_information/", scriptID, "_For_four_markers/", scriptID, "_Number/", scriptID, "_", markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4, "_genome_and_metabolome_195_varieties", ".pdf"))
g <- ggplot(markerInterest1234GroupDFCompleteCombination, aes(x = markerValueCombination, y = n, fill = group)) + geom_bar(stat = "identity") + labs(x = paste0(markerInterestID1, "_", markerInterestID2, "_", markerInterestID3, "_", markerInterestID4), y ="") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
plot(g)
dev.off()

