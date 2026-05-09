##########################################################################################
######  Title: 2.32_Soybean_STAM_GWAS_for_flavonoid_related_three_markers_for_metabolomic_data          ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/06/07 (Created), 2025/04/4 (Last Updated)                       ######
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

scriptID <- "2.32"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMGWAS <- paste0(dirMidSTAMBase, scriptID,
                         "_GWAS_for_three_markers/")
dir.create(dirMidSTAMGWAS)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)




thresLD <- 0.95


##### 1.3. Import packages #####
install.packages("GGally")
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)
require(plotly)
require(manhattanly)

source("scripts/1.1_Soybean_STAM_qqPlotly.R")

##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### 2. Perform gaston for Metabolomic data in 2017 ######
##### 2.1. Read marker genotype into R #####
gastonData0 <- gaston::read.vcf(file = "raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz")
lineNamesGeno <- gastonData0@ped$id

haploBlockList <- data.frame(fread("raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2._haplotype_block_list.csv"))

chrNos <- gastonData0@snps$chr
chrNames <- sprintf(fmt = paste0("Chr",
                                 "%0", floor(log10(max(chrNos))) + 1, "i"),
                    chrNos)
pos <- gastonData0@snps$pos
mrkNames <- paste0(chrNames, "_", pos)
gastonData0@snps$id <- mrkNames

blockNames0 <- haploBlockList$block
names(blockNames0) <- haploBlockList$marker
blockNames <- blockNames0[mrkNames]
gastonData0@snps$block <- blockNames

gastonDataSmall <- LD.thin(gastonData0, threshold = thresLD)




##### 2.2. Read genotypic values of Metabolomic data in 2017 of 0 or 2 for Chr06_18760995 into R #####
#### 2.2.1 For each maker ####
cultivationInfo <- "2017"
targetInfo <- "Metabolome"

targetType <- c("0", "2")

markerInterestID <- "Chr06_18760995"



gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
See(gvMetab2017Total)
rownames(gvMetab2017Total)
rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"


gvMetab2017Total <- scale(gvMetab2017Total, center = TRUE, scale = TRUE)
See(gvMetab2017Total)
gvMetab2017Total <- as.data.frame(gvMetab2017Total)

gastonData0Matrix <- as.matrix(gastonData0)

lineNames <- rownames(gvMetab2017Total)
markerInterest <- gastonData0Matrix[, markerInterestID]


markerInterestDF <- as.data.frame(markerInterest)


# markerVarietyNames <- names(markerInterest)

table(markerInterestDF$markerInterest)
marker0VarietyNames <- rownames(markerInterestDF)[markerInterestDF$markerInterest == 0]
marker2VarietyNames <- rownames(markerInterestDF)[markerInterestDF$markerInterest == 2]
See(marker0VarietyNames)
See(marker2VarietyNames)


lineNames0 <- lineNames[(lineNames%in%marker0VarietyNames)]
commonNames0 <- marker0VarietyNames[(marker0VarietyNames%in%lineNames0)]
lineNames2 <- lineNames[(lineNames%in%marker2VarietyNames)]
commonNames2 <- marker2VarietyNames[(marker2VarietyNames%in%lineNames2)]

gvMetab2017Total0 <- gvMetab2017Total[commonNames0, ]
gvMetab2017Total2 <- gvMetab2017Total[commonNames2, ]
# gvMetab2017Total0 <- na.omit(gvMetab2017Total0)
See(gvMetab2017Total0, rown = 20)
See(gvMetab2017Total2)




### LD1
markerNames <- c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")
betaMat0LD1 <- matrix(NA, nrow = length(gvMetab2017Total0), ncol = length(markerNames))
betaMat2LD1 <- matrix(NA, nrow = length(gvMetab2017Total2), ncol = length(markerNames))
rownames(betaMat0LD1) <- colnames(gvMetab2017Total0)
rownames(betaMat2LD1) <- colnames(gvMetab2017Total2)
colnames(betaMat0LD1) <- markerNames
colnames(betaMat2LD1) <- markerNames
See(betaMat0LD1)
See(betaMat2LD1)

# targetTypeNo <- 1
# traitNo <- 1

thresLD <- 2
for (targetTypeNo in 1:length(targetType)) {

  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- eval(parse( text = paste0("gvMetab2017Total", targetTypeNow)))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]

  See(gvMetab)

  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))



  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf >= 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)


  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)

  dirMarkerTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             markerInterestID,
                             "_results/")
  dir.create(dirMarkerTypeNow)

  dirTargetTypeNow <- paste0(dirMarkerTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  dirLDNow <- paste0(dirTargetTypeNow,
                     scriptID, "_",
                     "LD=1",
                     "_results/")
  dir.create(dirLDNow)


  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
      # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)

      KNow <- GRM(gastonDataNow)
      eigenKNow <- eigen(KNow)
    } else {
      gastonDataNow <- gastonData

      KNow <- K
      eigenKNow <- eigenK
    }
    st <- Sys.time()
    gastonRes <- association.test(x = gastonDataNow, method = "lmm",
                                  response = "quantitative",
                                  test = "wald", eigenK = eigenKNow,
                                  p = 2)
    end <- Sys.time()
    print(end - st)

    if (targetTypeNow == "0"){
      betaMat0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else {
      betaMat2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
      gastonRes$id[which(gastonRes$id %in% markerNames)]
      # table(gastonRes$id %in% "Chr17_16065902")
      # table(gastonRes$id %in% markerNames)
      # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]

    }


    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirLDNow, scriptID,
                       "_", traitNow, "/")
    dir.create(dirSave0)

    dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
    dir.create(dirSave)

    fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
                           thresLD, "_gaston_wald")

    fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
    write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)

    fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
    save(gastonRes, file = fileNameSaveRData)

    fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
    fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")

    png(fileNameSaveManhattan, width = 1200, height = 900)
    gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    dev.off()

    png(fileNameSaveQq, width = 900, height = 900)
    RAINBOWR::qq(- log10(gastonRes$p))
    dev.off()


    if (thresLD <= 0.4) {
      gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
                                            BP = gastonRes$pos,
                                            P = gastonRes$p,
                                            SNP = gastonRes$id,
                                            BLOCK = gastonDataNow@snps$block)
      sigSNPs <- gastonRes$id[pAdj < 0.05]

      fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
      fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")

      plotlyManhattan <- manhattanly(gastonResForManhattanly,
                                     snp = "SNP", gene = "BLOCK",
                                     highlight = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
                                               basename(fileNameSavePlotlyManhattan)))

      plotlyQq <- qqPlotly(data = gastonResForManhattanly,
                           highlightSNPNames = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
                                               basename(fileNameSavePlotlyQq)))
    }
    print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
  }
}


dir.create(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/"))


fileName0 <- paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_0.csv")
fileName2 <- paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID ,"_Chr06_18760995_2.csv")


write.csv(x = betaMat0LD1, file = fileName0)
write.csv(x = betaMat2LD1, file = fileName2)











#### comparison of coefficients
dirMidSTAMBase <- "midstream/"

dirMidSTAMGWAS <- paste0(dirMidSTAMBase, scriptID,
                         "_GWAS_for_three_markers/")


betaMat0 <- read.csv(paste0(dirMidSTAMGWAS, scriptID, "_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_0.csv"), row.names = 1)
betaMat2 <- read.csv(paste0(dirMidSTAMGWAS, scriptID, "_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_2.csv"), row.names = 1)
See(betaMat0)
See(betaMat2)


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]

table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("green4", nrow(betaMat0))
names(colVec) <- rownames(betaMat0)
colVec[metabNamesFlavonoid] <- "orange1"


#### Without Sd line ####
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_47490224.pdf"))
plot(betaMat0[, 1], betaMat2[, 1],
     xlim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     ylim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4,
     cex.axis = 1.3)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr10_42562665.pdf"))
plot(betaMat0[, 2], betaMat2[, 2],
     xlim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     ylim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr17_16065902.pdf"))
plot(betaMat0[, 3], betaMat2[, 3],
     xlim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     ylim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()


rownames(betaMat2)[betaMat2[, 2] >= 60]


pairs(cbind(betaMat0[, 2], betaMat2[, 2]))



#### With SD line, of significant value of flavonoid, considering Non-flavonoid as null distribution, and extract significant metabolites ####
betaMat0 <- read.csv(paste0(dirMidSTAMGWAS, scriptID, "_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_0.csv"), row.names = 1)
betaMat2 <- read.csv(paste0(dirMidSTAMGWAS, scriptID, "_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_2.csv"), row.names = 1)
See(betaMat0)
See(betaMat2)

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
metabNamesNonFlavonoid <- rownames(betaMat0)[!(rownames(betaMat0) %in% metabNamesFlavonoid)]
See(metabNamesNonFlavonoid)

betaMat0NonFlavonoid <- betaMat0[metabNamesNonFlavonoid, ]
betaMat0NonFlavonoid <- na.omit(betaMat0NonFlavonoid)
betaMat2NonFlavonoid <- betaMat2[metabNamesNonFlavonoid, ]
betaMat2NonFlavonoid <- na.omit(betaMat2NonFlavonoid)

betaMat0NonFlavonoidChr06_47490224Sd <- sqrt(mean((betaMat0NonFlavonoid[, "Chr06_47490224"] - mean(betaMat0NonFlavonoid[, "Chr06_47490224"]))^2))
betaMat0NonFlavonoidChr10_42562665Sd <- sqrt(mean((betaMat0NonFlavonoid[, "Chr10_42562665"] - mean(betaMat0NonFlavonoid[, "Chr10_42562665"]))^2))
betaMat0NonFlavonoidChr17_16065902Sd <- sqrt(mean((betaMat0NonFlavonoid[, "Chr17_16065902"] - mean(betaMat0NonFlavonoid[, "Chr17_16065902"]))^2))

betaMat2NonFlavonoidChr06_47490224Sd <- sqrt(mean((betaMat2NonFlavonoid[, "Chr06_47490224"] - mean(betaMat2NonFlavonoid[, "Chr06_47490224"]))^2))
betaMat2NonFlavonoidChr10_42562665Sd <- sqrt(mean((betaMat2NonFlavonoid[, "Chr10_42562665"] - mean(betaMat2NonFlavonoid[, "Chr10_42562665"]))^2))
betaMat2NonFlavonoidChr17_16065902Sd <- sqrt(mean((betaMat2NonFlavonoid[, "Chr17_16065902"] - mean(betaMat2NonFlavonoid[, "Chr17_16065902"]))^2))


### thresSd = 2
thresSd <- 2

betaMat0NonFlavonoidChr06_47490224ThresSd2<- thresSd * betaMat0NonFlavonoidChr06_47490224Sd
betaMat0NonFlavonoidChr10_42562665ThresSd2<- thresSd * betaMat0NonFlavonoidChr10_42562665Sd
betaMat0NonFlavonoidChr17_16065902ThresSd2<- thresSd * betaMat0NonFlavonoidChr17_16065902Sd

betaMat2NonFlavonoidChr06_47490224ThresSd2<- thresSd * betaMat2NonFlavonoidChr06_47490224Sd
betaMat2NonFlavonoidChr10_42562665ThresSd2<- thresSd * betaMat2NonFlavonoidChr10_42562665Sd
betaMat2NonFlavonoidChr17_16065902ThresSd2<- thresSd * betaMat2NonFlavonoidChr17_16065902Sd

intervalThresSd2BetaMat0FlavonoidChr06_47490224 <- c( mean(betaMat0NonFlavonoid[, "Chr06_47490224"])- betaMat0NonFlavonoidChr06_47490224ThresSd2, mean(betaMat0NonFlavonoid[, "Chr06_47490224"])+ betaMat0NonFlavonoidChr06_47490224ThresSd2)
intervalThresSd2BetaMat0FlavonoidChr10_42562665 <- c( mean(betaMat0NonFlavonoid[, "Chr10_42562665"])- betaMat0NonFlavonoidChr10_42562665ThresSd2, mean(betaMat0NonFlavonoid[, "Chr10_42562665"])+ betaMat0NonFlavonoidChr10_42562665ThresSd2)
intervalThresSd2BetaMat0FlavonoidChr17_16065902 <- c( mean(betaMat0NonFlavonoid[, "Chr17_16065902"])- betaMat0NonFlavonoidChr17_16065902ThresSd2, mean(betaMat0NonFlavonoid[, "Chr17_16065902"])+ betaMat0NonFlavonoidChr17_16065902ThresSd2)

intervalThresSd2BetaMat2FlavonoidChr06_47490224 <- c( mean(betaMat2NonFlavonoid[, "Chr06_47490224"])- betaMat2NonFlavonoidChr06_47490224ThresSd2, mean(betaMat2NonFlavonoid[, "Chr06_47490224"])+ betaMat2NonFlavonoidChr06_47490224ThresSd2)
intervalThresSd2BetaMat2FlavonoidChr10_42562665 <- c( mean(betaMat2NonFlavonoid[, "Chr10_42562665"])- betaMat2NonFlavonoidChr10_42562665ThresSd2, mean(betaMat2NonFlavonoid[, "Chr10_42562665"])+ betaMat2NonFlavonoidChr10_42562665ThresSd2)
intervalThresSd2BetaMat2FlavonoidChr17_16065902 <- c( mean(betaMat2NonFlavonoid[, "Chr17_16065902"])- betaMat2NonFlavonoidChr17_16065902ThresSd2, mean(betaMat2NonFlavonoid[, "Chr17_16065902"])+ betaMat2NonFlavonoidChr17_16065902ThresSd2)

intervalThresSd2BetaMat0FlavonoidChr06_47490224
intervalThresSd2BetaMat0FlavonoidChr10_42562665
intervalThresSd2BetaMat0FlavonoidChr17_16065902
intervalThresSd2BetaMat2FlavonoidChr06_47490224
intervalThresSd2BetaMat2FlavonoidChr10_42562665
intervalThresSd2BetaMat2FlavonoidChr17_16065902

### thresSd = 3
thresSd <- 3

betaMat0NonFlavonoidChr06_47490224ThresSd3<- thresSd * betaMat0NonFlavonoidChr06_47490224Sd
betaMat0NonFlavonoidChr10_42562665ThresSd3<- thresSd * betaMat0NonFlavonoidChr10_42562665Sd
betaMat0NonFlavonoidChr17_16065902ThresSd3<- thresSd * betaMat0NonFlavonoidChr17_16065902Sd

betaMat2NonFlavonoidChr06_47490224ThresSd3<- thresSd * betaMat2NonFlavonoidChr06_47490224Sd
betaMat2NonFlavonoidChr10_42562665ThresSd3<- thresSd * betaMat2NonFlavonoidChr10_42562665Sd
betaMat2NonFlavonoidChr17_16065902ThresSd3<- thresSd * betaMat2NonFlavonoidChr17_16065902Sd

intervalThresSd3BetaMat0FlavonoidChr06_47490224 <- c( mean(betaMat0NonFlavonoid[, "Chr06_47490224"])- betaMat0NonFlavonoidChr06_47490224ThresSd3, mean(betaMat0NonFlavonoid[, "Chr06_47490224"])+ betaMat0NonFlavonoidChr06_47490224ThresSd3)
intervalThresSd3BetaMat0FlavonoidChr10_42562665 <- c( mean(betaMat0NonFlavonoid[, "Chr10_42562665"])- betaMat0NonFlavonoidChr10_42562665ThresSd3, mean(betaMat0NonFlavonoid[, "Chr10_42562665"])+ betaMat0NonFlavonoidChr10_42562665ThresSd3)
intervalThresSd3BetaMat0FlavonoidChr17_16065902 <- c( mean(betaMat0NonFlavonoid[, "Chr17_16065902"])- betaMat0NonFlavonoidChr17_16065902ThresSd3, mean(betaMat0NonFlavonoid[, "Chr17_16065902"])+ betaMat0NonFlavonoidChr17_16065902ThresSd3)

intervalThresSd3BetaMat2FlavonoidChr06_47490224 <- c( mean(betaMat2NonFlavonoid[, "Chr06_47490224"])- betaMat2NonFlavonoidChr06_47490224ThresSd3, mean(betaMat2NonFlavonoid[, "Chr06_47490224"])+ betaMat2NonFlavonoidChr06_47490224ThresSd3)
intervalThresSd3BetaMat2FlavonoidChr10_42562665 <- c( mean(betaMat2NonFlavonoid[, "Chr10_42562665"])- betaMat2NonFlavonoidChr10_42562665ThresSd3, mean(betaMat2NonFlavonoid[, "Chr10_42562665"])+ betaMat2NonFlavonoidChr10_42562665ThresSd3)
intervalThresSd3BetaMat2FlavonoidChr17_16065902 <- c( mean(betaMat2NonFlavonoid[, "Chr17_16065902"])- betaMat2NonFlavonoidChr17_16065902ThresSd3, mean(betaMat2NonFlavonoid[, "Chr17_16065902"])+ betaMat2NonFlavonoidChr17_16065902ThresSd3)

intervalThresSd3BetaMat0FlavonoidChr06_47490224
intervalThresSd3BetaMat0FlavonoidChr10_42562665
intervalThresSd3BetaMat0FlavonoidChr17_16065902
intervalThresSd3BetaMat2FlavonoidChr06_47490224
intervalThresSd3BetaMat2FlavonoidChr10_42562665
intervalThresSd3BetaMat2FlavonoidChr17_16065902


#PDF files with SD lines
dirMidSTAMBase <- "midstream/"

dirMidSTAMGWAS <- paste0(dirMidSTAMBase, scriptID,
                         "_GWAS_for_three_markers/")


betaMat0 <- read.csv(paste0(dirMidSTAMGWAS, scriptID, "_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_0.csv"), row.names = 1)
betaMat2 <- read.csv(paste0(dirMidSTAMGWAS, scriptID, "_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_2.csv"), row.names = 1)
See(betaMat0)
See(betaMat2)


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"), row.names = 1)
metabNamesFlavonoid <- rownames(metabFlavonoid)

table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("blue", nrow(betaMat0))
names(colVec) <- rownames(betaMat0)
colVec[metabNamesFlavonoid] <- "orange1"

betaMat0Flavonoid <- betaMat0[metabNamesFlavonoid, ]
betaMat2Flavonoid <- betaMat2[metabNamesFlavonoid, ]



# thresSd = 2
# Chr06_47490224; thresSd = 2, 3, 2&3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_47490224_thresSD_2.pdf"))
plot(betaMat0[, 1], betaMat2[, 1],
     xlim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     ylim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4,
     cex.axis = 1.3)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0FlavonoidChr06_47490224, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2FlavonoidChr06_47490224, col = 4, lty = 4, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_47490224_thresSD_3.pdf"))
plot(betaMat0[, 1], betaMat2[, 1],
     xlim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     ylim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4,
     cex.axis = 1.3)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0FlavonoidChr06_47490224, col = 1, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2FlavonoidChr06_47490224, col = 1, lty = 4, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_47490224_thresSD_2_3.pdf"))
plot(betaMat0[, 1], betaMat2[, 1],
     xlim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     ylim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4,
     cex.axis = 1.3)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0FlavonoidChr06_47490224, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0FlavonoidChr06_47490224, col = 1, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2FlavonoidChr06_47490224, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2FlavonoidChr06_47490224, col = 1, lty = 4, lwd = 1.5)
dev.off()

# Chr10_42562665; thresSd = 2, 3, 2&3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr10_42562665_thresSD_2.pdf"))
plot(betaMat0[, 2], betaMat2[, 2],
     xlim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     ylim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr10_42562665_thresSD_thresSD_3.pdf"))
plot(betaMat0[, 2], betaMat2[, 2],
     xlim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     ylim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0FlavonoidChr10_42562665, col = 1, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2FlavonoidChr10_42562665, col = 1, lty = 4, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr10_42562665_thresSD_2_3.pdf"))
plot(betaMat0[, 2], betaMat2[, 2],
     xlim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     ylim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0FlavonoidChr10_42562665, col = 1, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2FlavonoidChr10_42562665, col = 1, lty = 4, lwd = 1.5)
dev.off()

# Chr17_16065902; thresSd = 2, 3, 2&3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr17_16065902_thresSD_2.pdf"))
plot(betaMat0[, 3], betaMat2[, 3],
     xlim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     ylim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr17_16065902_thresSD_3.pdf"))
plot(betaMat0[, 3], betaMat2[, 3],
     xlim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     ylim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0FlavonoidChr17_16065902, col = 1, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2FlavonoidChr17_16065902, col = 1, lty = 4, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr17_16065902_thresSD_2_3.pdf"))
plot(betaMat0[, 3], betaMat2[, 3],
     xlim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     ylim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0FlavonoidChr17_16065902, col = 1, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2FlavonoidChr17_16065902, col = 1, lty = 4, lwd = 1.5)
dev.off()



#### Extract metabolites
betaMat0SortedBy0 <- arrange(betaMat0, Chr06_47490224)
betaMat0Flavonoid <- betaMat0[metabNamesFlavonoid, ]
betaMat2Flavonoid <- betaMat2[metabNamesFlavonoid, ]
See(betaMat0Flavonoid)
See(betaMat2Flavonoid)

### Chr06_47490224

rownames(betaMat0)[betaMat0[, 1] < -0.4 & betaMat2[, 1] < 0.1]
rownames(betaMat0[betaMat0[, 1] > 0.5, ])
rownames(betaMat2)[betaMat2[, 1] < -1.0 & betaMat0[, 1] > -0.3 ]
betaMat0[rownames(betaMat2)[betaMat2[, 1] < -1.0 & betaMat0[, 1] > -0.3 ], ]
rownames(betaMat0)[betaMat0[, 1] < -0.3 & abs(betaMat2[, 1]) < 0.1]

rownames(betaMat0Flavonoid)[betaMat0Flavonoid[, 1] < -0.4 & abs(betaMat2Flavonoid[, 1]) <0.1]
rownames(betaMat0)[betaMat0[, 1] < 0 & betaMat0[, 1] > -0.5 & betaMat2[, 1] < -1.0]


betaMat0[order(betaMat0$Chr06_47490224), ]
betaMat2[order(betaMat2$Chr06_47490224), ]


metabNames20 <- rownames(betaMat0Flavonoid)[betaMat0Flavonoid[, 1] < -0.3 & abs(betaMat2Flavonoid[, 1]) <0.1]
# metabNames20 <- rownames(betaMat0)[betaMat0[, 1] < -0.3 & abs(betaMat2[, 1]) <0.1]
metabNames00 <- rownames(betaMat0)[betaMat0[, 1] < 0 & betaMat0[, 1] > -0.5 & betaMat2[, 1] < -0.7]
metabNames00And20 <- c(metabNames20, metabNames00)






metabNamesFlavonoidLessThanThresSd3BetaMat0Chr06_47490224 <- rownames(betaMat0Flavonoid)[
  betaMat0Flavonoid[, 1] < intervalThresSd3BetaMat0FlavonoidChr06_47490224[1] &
    betaMat2Flavonoid[, 1] > intervalThresSd3BetaMat2FlavonoidChr06_47490224[1] &
    betaMat2Flavonoid[, 1] < intervalThresSd3BetaMat2FlavonoidChr06_47490224[2]
  ]
metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr06_47490224 <- rownames(betaMat0Flavonoid)[
  betaMat0Flavonoid[, 1] > intervalThresSd3BetaMat0FlavonoidChr06_47490224[2] &
    betaMat2Flavonoid[, 1] > intervalThresSd3BetaMat2FlavonoidChr06_47490224[1] &
    betaMat2Flavonoid[, 1] < intervalThresSd3BetaMat2FlavonoidChr06_47490224[2]
  ]

metabNamesFlavonoidLessThanThresSd3BetaMat2Chr06_47490224 <- rownames(betaMat2Flavonoid)[
   # betaMat0Flavonoid[, 1] < intervalThresSd3BetaMat0FlavonoidChr06_47490224[2] &
   #  betaMat0Flavonoid[, 1] > intervalThresSd3BetaMat0FlavonoidChr06_47490224[1] &
    betaMat2Flavonoid[, 1] < intervalThresSd3BetaMat2FlavonoidChr06_47490224[1]
  ]
metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr06_47490224 <- rownames(betaMat2Flavonoid)[
  # betaMat0Flavonoid[, 1] < intervalThresSd3BetaMat0FlavonoidChr06_47490224[2] &
    # betaMat0Flavonoid[, 1] > intervalThresSd3BetaMat0FlavonoidChr06_47490224[1] &
    betaMat2Flavonoid[, 1] > intervalThresSd3BetaMat2FlavonoidChr06_47490224[2]
  ]

### Non-flavonoid related
rownames(betaMat2)[
  # betaMat0[, 1] < intervalThresSd3BetaMat0FlavonoidChr06_47490224[2] &
  #   betaMat0Flavonoid[, 1] > intervalThresSd3BetaMat0FlavonoidChr06_47490224[1] &
    betaMat2[, 1] < intervalThresSd3BetaMat2FlavonoidChr06_47490224[1]
]


### Chr10_42562665
rownames(betaMat0[betaMat0[, 2] > 0.6, ])
rownames(betaMat0[betaMat0[, 2] > 0.4, ])
rownames(betaMat2[betaMat2[, 2] > 1, ])
rownames(betaMat0[betaMat0[, 2] < -0.7, ])
rownames(betaMat0[betaMat0[, 2] < -0.5, ])
rownames(betaMat2[betaMat2[, 2] < -1, ])


betaMat0[order(betaMat0$Chr10_42562665), ]
betaMat2[order(betaMat2$Chr10_42562665), ]





metabNamesFlavonoidLessThanThresSd3BetaMat0Chr10_42562665 <- rownames(betaMat0Flavonoid)[
  betaMat0Flavonoid[, 2] < intervalThresSd3BetaMat0FlavonoidChr10_42562665[1] &
    betaMat2Flavonoid[, 2] > intervalThresSd3BetaMat2FlavonoidChr10_42562665[1] &
    betaMat2Flavonoid[, 2] < intervalThresSd3BetaMat2FlavonoidChr10_42562665[2]
]
metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr10_42562665 <- rownames(betaMat0Flavonoid)[
  betaMat0Flavonoid[, 2] > intervalThresSd3BetaMat0FlavonoidChr10_42562665[2] &
    betaMat2Flavonoid[, 2] > intervalThresSd3BetaMat2FlavonoidChr10_42562665[1] &
    betaMat2Flavonoid[, 2] < intervalThresSd3BetaMat2FlavonoidChr10_42562665[2]
]

metabNamesFlavonoidLessThanThresSd3BetaMat2Chr10_42562665 <- rownames(betaMat2Flavonoid)[
  # betaMat0Flavonoid[, 2] < intervalThresSd3BetaMat0FlavonoidChr10_42562665[2] &
    # betaMat0Flavonoid[, 2] > intervalThresSd3BetaMat0FlavonoidChr10_42562665[1] &
    betaMat2Flavonoid[, 2] < intervalThresSd3BetaMat2FlavonoidChr10_42562665[1]
]
metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr10_42562665 <- rownames(betaMat2Flavonoid)[
  # betaMat0Flavonoid[, 2] < intervalThresSd3BetaMat0FlavonoidChr10_42562665[2] &
    # betaMat0Flavonoid[, 2] > intervalThresSd3BetaMat0FlavonoidChr10_42562665[1] &
    betaMat2Flavonoid[, 2] > intervalThresSd3BetaMat2FlavonoidChr10_42562665[2]
]

### Non-Flavonoid related
rownames(betaMat2)[
  # betaMat0Flavonoid[, 2] < intervalThresSd3BetaMat0FlavonoidChr10_42562665[2] &
  # betaMat0Flavonoid[, 2] > intervalThresSd3BetaMat0FlavonoidChr10_42562665[1] &
  betaMat2[, 2] > intervalThresSd3BetaMat2FlavonoidChr10_42562665[2]
]



### Chr17_16065902
rownames(betaMat0[betaMat0[, 3] > 0.25, ])
rownames(betaMat2)[betaMat2[, 3] < -0.3]


betaMat0[order(betaMat0$Chr17_16065902), ]
betaMat2[order(betaMat2$Chr17_16065902), ]



metabNamesFlavonoidLessThanThresSd3BetaMat0Chr17_16065902 <- rownames(betaMat0Flavonoid)[
  betaMat0Flavonoid[, 3] < intervalThresSd3BetaMat0FlavonoidChr17_16065902[1] &
    betaMat2Flavonoid[, 3] > intervalThresSd3BetaMat2FlavonoidChr17_16065902[1] &
    betaMat2Flavonoid[, 3] < intervalThresSd3BetaMat2FlavonoidChr17_16065902[2]
]
metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr17_16065902 <- rownames(betaMat0Flavonoid)[
  betaMat0Flavonoid[, 3] > intervalThresSd3BetaMat0FlavonoidChr17_16065902[2] &
    betaMat2Flavonoid[, 3] > intervalThresSd3BetaMat2FlavonoidChr17_16065902[1] &
    betaMat2Flavonoid[, 3] < intervalThresSd3BetaMat2FlavonoidChr17_16065902[2]
]

metabNamesFlavonoidLessThanThresSd3BetaMat2Chr17_16065902 <- rownames(betaMat2Flavonoid)[
  # betaMat0Flavonoid[, 3] < intervalThresSd3BetaMat0FlavonoidChr17_16065902[2] &
  # betaMat0Flavonoid[, 3] > intervalThresSd3BetaMat0FlavonoidChr17_16065902[1] &
  betaMat2Flavonoid[, 3] < intervalThresSd3BetaMat2FlavonoidChr17_16065902[1]
]
metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr17_16065902 <- rownames(betaMat2Flavonoid)[
  # betaMat0Flavonoid[, 3] < intervalThresSd3BetaMat0FlavonoidChr17_16065902[2] &
  # betaMat0Flavonoid[, 3] > intervalThresSd3BetaMat0FlavonoidChr17_16065902[1] &
  betaMat2Flavonoid[, 3] > intervalThresSd3BetaMat2FlavonoidChr17_16065902[2]
]


metabNamesFlavonoidLessThanThresSd3BetaMat0Chr06_47490224
metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr06_47490224
metabNamesFlavonoidLessThanThresSd3BetaMat2Chr06_47490224
metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr06_47490224

metabNamesFlavonoidLessThanThresSd3BetaMat0Chr10_42562665
metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr10_42562665
metabNamesFlavonoidLessThanThresSd3BetaMat2Chr10_42562665
metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr10_42562665

metabNamesFlavonoidLessThanThresSd3BetaMat0Chr17_16065902
metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr17_16065902
metabNamesFlavonoidLessThanThresSd3BetaMat2Chr17_16065902
metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr17_16065902


metabNamesFlavonoidAllDetectedThresSd3BetaMat0And2Chr06_47490224 <- c(metabNamesFlavonoidLessThanThresSd3BetaMat0Chr06_47490224, metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr06_47490224, metabNamesFlavonoidLessThanThresSd3BetaMat2Chr06_47490224, metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr06_47490224)
metabNamesFlavonoidAllDetectedThresSd3BetaMat0And2Chr06_47490224 <- na.omit(metabNamesFlavonoidAllDetectedThresSd3BetaMat0And2Chr06_47490224)

# write csv
dirMetabCoefficient <- paste("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/")
write.csv(x = metabNamesFlavonoidAllDetectedThresSd3BetaMat0And2Chr06_47490224, file = paste0(dirMetabCoefficient, "2.32_metabNames_Flavonoid_All_Detected_ThresSd3_BetaMat_0_And_2_Chr06_47490224.csv"))




#### Check annotation of extracted metabolites ####
metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"), row.names = 1)
# metabNamesFlavonoid <- rownames(metabFlavonoid)

annotationMetabFlavonoidLessThanThresSd3BetaMat0Chr06_47490224 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0Chr06_47490224, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0Chr06_47490224 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr06_47490224, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat2Chr06_47490224 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2Chr06_47490224, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat2Chr06_47490224 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr06_47490224, ]

annotationMetabFlavonoidLessThanThresSd3BetaMat0Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0Chr10_42562665, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr10_42562665, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat2Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2Chr10_42562665, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat2Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr10_42562665, ]

annotationMetabFlavonoidLessThanThresSd3BetaMat0Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0Chr17_16065902, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0Chr17_16065902, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat2Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2Chr17_16065902, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat2Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2Chr17_16065902, ]


annotationMetabFlavonoidLessThanThresSd3BetaMat0Chr06_47490224
annotationMetabFlavonoidMoreThanThresSd3BetaMat0Chr06_47490224
annotationMetabFlavonoidLessThanThresSd3BetaMat2Chr06_47490224
annotationMetabFlavonoidMoreThanThresSd3BetaMat2Chr06_47490224

annotationMetabFlavonoidLessThanThresSd3BetaMat0Chr10_42562665
annotationMetabFlavonoidMoreThanThresSd3BetaMat0Chr10_42562665
annotationMetabFlavonoidLessThanThresSd3BetaMat2Chr10_42562665
annotationMetabFlavonoidMoreThanThresSd3BetaMat2Chr10_42562665

annotationMetabFlavonoidLessThanThresSd3BetaMat0Chr17_16065902
annotationMetabFlavonoidMoreThanThresSd3BetaMat0Chr17_16065902
annotationMetabFlavonoidLessThanThresSd3BetaMat2Chr17_16065902
annotationMetabFlavonoidMoreThanThresSd3BetaMat2Chr17_16065902


# memo
# intervalThresSd2BetaMat0FlavonoidChr06_47490224
# intervalThresSd2BetaMat0FlavonoidChr10_42562665
# intervalThresSd2BetaMat0FlavonoidChr17_16065902
# intervalThresSd2BetaMat2FlavonoidChr06_47490224
# intervalThresSd2BetaMat2FlavonoidChr10_42562665
# intervalThresSd2BetaMat2FlavonoidChr17_16065902
#
# intervalThresSd3BetaMat0FlavonoidChr06_47490224
# intervalThresSd3BetaMat0FlavonoidChr10_42562665
# intervalThresSd3BetaMat0FlavonoidChr17_16065902
# intervalThresSd3BetaMat2FlavonoidChr06_47490224
# intervalThresSd3BetaMat2FlavonoidChr10_42562665
# intervalThresSd3BetaMat2FlavonoidChr17_16065902




#### 2.2.2. For two marker combinations ( Chr06_18760995, Chr06_47490224 )
markerNames <- c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")
betaMat0LD1 <- matrix(NA, nrow = length(gvMetab2017Total0), ncol = length(markerNames))
betaMat2LD1 <- matrix(NA, nrow = length(gvMetab2017Total2), ncol = length(markerNames))
rownames(betaMat0LD1) <- colnames(gvMetab2017Total0)
rownames(betaMat2LD1) <- colnames(gvMetab2017Total2)
colnames(betaMat0LD1) <- markerNames
colnames(betaMat2LD1) <- markerNames
See(betaMat0LD1)
See(betaMat2LD1)

# targetTypeNo <- 1
# traitNo <- 1

thresLD <- 2
for (targetTypeNo in 1:length(targetType)) {

  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- eval(parse( text = paste0("gvMetab2017Total", targetTypeNow)))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]

  See(gvMetab)

  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))



  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf >= 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)


  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)

  dirMarkerTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             markerInterestID,
                             "_results/")
  dir.create(dirMarkerTypeNow)

  dirTargetTypeNow <- paste0(dirMarkerTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  dirLDNow <- paste0(dirTargetTypeNow,
                     scriptID, "_",
                     "LD=1",
                     "_results/")
  dir.create(dirLDNow)


  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
      # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)

      KNow <- GRM(gastonDataNow)
      eigenKNow <- eigen(KNow)
    } else {
      gastonDataNow <- gastonData

      KNow <- K
      eigenKNow <- eigenK
    }
    st <- Sys.time()
    gastonRes <- association.test(x = gastonDataNow, method = "lmm",
                                  response = "quantitative",
                                  test = "wald", eigenK = eigenKNow,
                                  p = 2)
    end <- Sys.time()
    print(end - st)

    if (targetTypeNow == "0"){
      betaMat0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else {
      betaMat2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
      gastonRes$id[which(gastonRes$id %in% markerNames)]
      # table(gastonRes$id %in% "Chr17_16065902")
      # table(gastonRes$id %in% markerNames)
      # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]

    }


    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirLDNow, scriptID,
                       "_", traitNow, "/")
    dir.create(dirSave0)

    dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
    dir.create(dirSave)

    fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
                           thresLD, "_gaston_wald")

    fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
    write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)

    fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
    save(gastonRes, file = fileNameSaveRData)

    # fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
    # fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")
    #
    # png(fileNameSaveManhattan, width = 1200, height = 900)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    # dev.off()
    #
    # png(fileNameSaveQq, width = 900, height = 900)
    # RAINBOWR::qq(- log10(gastonRes$p))
    # dev.off()
    #
    #
    # if (thresLD <= 0.4) {
    #   gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
    #                                         BP = gastonRes$pos,
    #                                         P = gastonRes$p,
    #                                         SNP = gastonRes$id,
    #                                         BLOCK = gastonDataNow@snps$block)
    #   sigSNPs <- gastonRes$id[pAdj < 0.05]
    #
    #   fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
    #   fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")
    #
    #   plotlyManhattan <- manhattanly(gastonResForManhattanly,
    #                                  snp = "SNP", gene = "BLOCK",
    #                                  highlight = sigSNPs)
    #   htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
    #                           file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
    #                                            basename(fileNameSavePlotlyManhattan)))
    #
    #   plotlyQq <- qqPlotly(data = gastonResForManhattanly,
    #                        highlightSNPNames = sigSNPs)
    #   htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
    #                           file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
    #                                            basename(fileNameSavePlotlyQq)))
    # }
    print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
  }
}


# dir.create(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/"))



fileName <- paste0(dirMidSTAMGWAS, scriptID, "_Chr06_18760995_2_coefficient_LD1.csv")
write.csv(x = betaMat2LD1, file = fileName)






##### 2.3. Read genotypic values of Metabolomic data in 2017 of (0,0),(0,2),(2,0),(2,2) for Chr06_18760995 and Chr06_47490224 into R #####
#### 2.3.1 For each maker ####
cultivationInfo <- "2017"
targetInfo <- "Metabolome"

targetType <- c("0_0", "0_2", "2_0", "2_2")
# targetType <- c("0_2", "2_0", "2_2")

markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47490224"



gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
See(gvMetab2017Total)


gvMetab2017Total <- scale(gvMetab2017Total, center = TRUE, scale = TRUE)
See(gvMetab2017Total)
gvMetab2017Total <- as.data.frame(gvMetab2017Total)

lineNames <- rownames(gvMetab2017Total)
markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]


markerInterest1DF <- as.data.frame(markerInterest1)
markerInterest2DF <- as.data.frame(markerInterest2)
markerInterest1And2DF <- cbind(markerInterest1DF, markerInterest2DF)
See(markerInterest1And2DF)

# markerVarietyNames <- names(markerInterest)

table(markerInterest1And2DF$markerInterest1)
table(markerInterest1And2DF$markerInterest2)
marker0_0VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 0]
marker0_2VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 2]
marker2_0VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 0]
marker2_2VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 2]

See(marker0_0VarietyNames)
See(marker0_2VarietyNames)
See(marker2_0VarietyNames)
See(marker2_2VarietyNames)


lineNames0_0 <- lineNames[(lineNames %in% marker0_0VarietyNames)]
commonNames0_0 <- marker0_0VarietyNames[(marker0_0VarietyNames %in% lineNames0_0)]
lineNames0_2 <- lineNames[(lineNames %in% marker0_2VarietyNames)]
commonNames0_2 <- marker0_2VarietyNames[(marker0_2VarietyNames %in% lineNames0_2)]
lineNames2_0 <- lineNames[(lineNames %in% marker2_0VarietyNames)]
commonNames2_0 <- marker2_0VarietyNames[(marker2_0VarietyNames %in% lineNames2_0)]
lineNames2_2 <- lineNames[(lineNames %in% marker2_2VarietyNames)]
commonNames2_2 <- marker2_2VarietyNames[(marker2_2VarietyNames %in% lineNames2_2)]

gvMetab2017Total0_0 <- gvMetab2017Total[commonNames0_0, ]
gvMetab2017Total0_2 <- gvMetab2017Total[commonNames0_2, ]
gvMetab2017Total2_0 <- gvMetab2017Total[commonNames2_0, ]
gvMetab2017Total2_2 <- gvMetab2017Total[commonNames2_2, ]
# gvMetab2017Total0 <- na.omit(gvMetab2017Total0)
See(gvMetab2017Total0_0, rown = 20)
See(gvMetab2017Total0_2)
See(gvMetab2017Total2_0)
See(gvMetab2017Total2_2)




### LD1
markerNames <- c("Chr10_42562665", "Chr17_16065902")
betaMat0_0LD1 <- matrix(NA, nrow = length(gvMetab2017Total0_0), ncol = length(markerNames))
betaMat0_2LD1 <- matrix(NA, nrow = length(gvMetab2017Total0_2), ncol = length(markerNames))
betaMat2_0LD1 <- matrix(NA, nrow = length(gvMetab2017Total2_0), ncol = length(markerNames))
betaMat2_2LD1 <- matrix(NA, nrow = length(gvMetab2017Total2_2), ncol = length(markerNames))
rownames(betaMat0_0LD1) <- colnames(gvMetab2017Total0_0)
rownames(betaMat0_2LD1) <- colnames(gvMetab2017Total0_2)
rownames(betaMat2_0LD1) <- colnames(gvMetab2017Total2_0)
rownames(betaMat2_2LD1) <- colnames(gvMetab2017Total2_2)
colnames(betaMat0_0LD1) <- markerNames
colnames(betaMat0_2LD1) <- markerNames
colnames(betaMat2_0LD1) <- markerNames
colnames(betaMat2_2LD1) <- markerNames
See(betaMat0_0LD1)
See(betaMat0_2LD1)
See(betaMat2_0LD1)
See(betaMat2_2LD1)

targetTypeNo <- 1
traitNo <- 1

thresLD <- 2
for (targetTypeNo in 1:length(targetType)) {

  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- eval(parse( text = paste0("gvMetab2017Total", targetTypeNow)))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]

  See(gvMetab)

  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))



  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf >= 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)


  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)

  dirMarkerTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             markerInterestID1,
                             "_",
                             markerInterestID2,
                             "_results/")
  dir.create(dirMarkerTypeNow)

  dirTargetTypeNow <- paste0(dirMarkerTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  dirLDNow <- paste0(dirTargetTypeNow,
                     scriptID, "_",
                     "LD=1",
                     "_results/")
  dir.create(dirLDNow)


  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
      # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)

      KNow <- GRM(gastonDataNow)
      eigenKNow <- eigen(KNow)
    } else {
      gastonDataNow <- gastonData

      KNow <- K
      eigenKNow <- eigenK
    }
    st <- Sys.time()
    gastonRes <- association.test(x = gastonDataNow, method = "lmm",
                                  response = "quantitative",
                                  test = "wald", eigenK = eigenKNow,
                                  p = 2)
    end <- Sys.time()
    print(end - st)

    if (targetTypeNow == "0_0"){
      betaMat0_0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "0_2") {
      betaMat0_2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "2_0") {
      betaMat2_0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "2_2") {
      betaMat2_2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
      # gastonRes$id[which(gastonRes$id %in% markerNames)]
      # table(gastonRes$id %in% "Chr17_16065902")
      # table(gastonRes$id %in% markerNames)
      # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]

    }


    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirLDNow, scriptID,
                       "_", traitNow, "/")
    dir.create(dirSave0)

    dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
    dir.create(dirSave)

    fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
                           thresLD, "_gaston_wald")

    fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
    write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)

    fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
    save(gastonRes, file = fileNameSaveRData)

    fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
    fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")

    png(fileNameSaveManhattan, width = 1200, height = 900)
    gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    dev.off()

    png(fileNameSaveQq, width = 900, height = 900)
    RAINBOWR::qq(- log10(gastonRes$p))
    dev.off()


    if (thresLD <= 0.4) {
      gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
                                            BP = gastonRes$pos,
                                            P = gastonRes$p,
                                            SNP = gastonRes$id,
                                            BLOCK = gastonDataNow@snps$block)
      sigSNPs <- gastonRes$id[pAdj < 0.05]

      fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
      fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")

      plotlyManhattan <- manhattanly(gastonResForManhattanly,
                                     snp = "SNP", gene = "BLOCK",
                                     highlight = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
                                               basename(fileNameSavePlotlyManhattan)))

      plotlyQq <- qqPlotly(data = gastonResForManhattanly,
                           highlightSNPNames = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
                                               basename(fileNameSavePlotlyQq)))
    }
    print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
  }
}


dirTwoMarkers2017 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/")
dir.create(dirTwoMarkers2017)

fileName0_0 <- paste0(dirTwoMarkers2017, scriptID, "_(0,0).csv")
fileName0_2 <- paste0(dirTwoMarkers2017, scriptID, "_(0,2).csv")
fileName2_0 <- paste0(dirTwoMarkers2017, scriptID, "_(2,0).csv")
fileName2_2 <- paste0(dirTwoMarkers2017, scriptID, "_(2,2).csv")

write.csv(x = betaMat0_0LD1, file = fileName0_0)
write.csv(x = betaMat0_2LD1, file = fileName0_2)
write.csv(x = betaMat2_0LD1, file = fileName2_0)
write.csv(x = betaMat2_2LD1, file = fileName2_2)



##### comparison of coefficients for upstream two markers #####
#### Without SD line ####
dirTwoMarkers2017 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/")
dir.create(dirTwoMarkers2017)


betaMat0_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,0).csv"), row.names = 1)
betaMat0_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,2).csv"), row.names = 1)
betaMat2_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,0).csv"), row.names = 1)
betaMat2_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,2).csv"), row.names = 1)


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]

table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("blue1", nrow(betaMat0))
names(colVec) <- rownames(betaMat0)
colVec[metabNamesFlavonoid] <- "orange1"



dir.create(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/"))
dir.create(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/"))

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(0,2).pdf"))
plot(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(2,0).pdf"))
plot(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(2,2).pdf"))
plot(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_(0,2)_(2,0).pdf"))
plot(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1],
     xlim = range(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_(0,2)_(2,2).pdf"))
plot(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_(2,0)_(2,2).pdf"))
plot(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()



pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,0)_(0,2).pdf"))
plot(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,0)_(2,0).pdf"))
plot(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,0)_(2,2).pdf"))
plot(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,2)_(2,0).pdf"))
plot(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2],
     xlim = range(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,2)_(2,2).pdf"))
plot(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(2,0)_(2,2).pdf"))
plot(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()



### pairs plot
Chr10_42562665DF <- data.frame("0_0" = betaMat0_0LD1[, 1], "0_2" = betaMat0_2LD1[, 1], "2_0" = betaMat2_0LD1[, 1], "2_2" = betaMat2_2LD1[, 1])
colnames(Chr10_42562665DF) <- c("0_0", "0_2", "2_0", "2_2")
Chr17_16065902DF <- data.frame("0_0" = betaMat0_0LD1[, 2], "0_2" = betaMat0_2LD1[, 2], "2_0" = betaMat2_0LD1[, 2], "2_2" = betaMat2_2LD1[, 2])
colnames(Chr17_16065902DF) <- c("0_0", "0_2", "2_0", "2_2")


# upper <- function(x, y, ...){
#   oldpar <- par(usr = c(0, 1, 0, 1))
#   v <- abs(cor(x, y))
#   sz <- 1 + v * 2
#   text(0.5, 0.5, sprintf("%.3f", v), cex = sz)
#   par(oldpar)
# }

lower <- function(x, y, ...){
  points(x, y, pch = 21, bg = c("blue1", "orange1")[as.numeric(as.factor(colVec))])
  if( 1==1 ){
    abline(0, 1, col = 2, lty = 2, lwd = 1.5)
  }
}

# pairs(Chr10_42562665DF, upper.panel=upper, lower.panel=lower)
pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr10_42562665/", scriptID, "_pairs_plot.pdf"))
pairs(Chr10_42562665DF, lower.panel=lower, oma = c(3,3,3,15))
par(xpd = TRUE)
legend(x = "bottomright", legend = c("Non-flavonoid", "Flavonoid"), col = c("blue1", "orange1"), pch = 19, cex = 1)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_pairs_plot.pdf"))
pairs(Chr17_16065902DF, lower.panel=lower, oma = c(3,3,3,15))
par(xpd = TRUE)
legend(x = "bottomright", legend = c("Non-flavonoid", "Flavonoid"), col = c("blue1", "orange1"), pch = 19)
dev.off()



# Chr10_42562665DF <- data.frame("0_0" = betaMat0_0LD1[, 1], "0_2" = betaMat0_2LD1[, 1], "2_0" = betaMat2_0LD1[, 1], "2_2" = betaMat2_2LD1[, 1])
# colnames(Chr10_42562665DF) <- c("0_0", "0_2", "2_0", "2_2")
# Chr17_16065902DF <- data.frame("0_0" = betaMat0_0LD1[, 2], "0_2" = betaMat0_2LD1[, 2], "2_0" = betaMat2_0LD1[, 2], "2_2" = betaMat2_2LD1[, 2])
# colnames(Chr17_16065902DF) <- c("0_0", "0_2", "2_0", "2_2")
#
# See(Chr10_42562665DF)
# See(Chr17_16065902DF)


# ggpairs(Chr10_42562665DF, aes_string(colour="colVec", alpha=0.5)) + geom_abline(slope = 1, intercept = 0)
# ggpairs(Chr17_16065902DF, aes_string(colour="colVec", alpha=0.5))
#
#
# ggpairs(Chr10_42562665DF, mapping = aes(color = colVec), lower = list(continuous = "smooth"))
# ggpairs(Chr17_16065902DF, mapping = aes(color = colVec), lower = list(continuous = "smooth"))











####  With SD line, of significant value of flavonoid, considering Non-flavonoid as null distribution, and extract significant metabolites ####
betaMat0_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(0,0).csv", row.names = 1)
betaMat0_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(0,2).csv", row.names = 1)
betaMat2_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(2,0).csv", row.names = 1)
betaMat2_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(2,2).csv", row.names = 1)

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
metabNamesNonFlavonoid <- rownames(betaMat0_0LD1)[!(rownames(betaMat0_0LD1) %in% metabNamesFlavonoid)]
See(metabNamesNonFlavonoid)

betaMat0_0NonFlavonoid <- betaMat0_0LD1[metabNamesNonFlavonoid, ]
betaMat0_0NonFlavonoid <- na.omit(betaMat0_0NonFlavonoid)
betaMat0_2NonFlavonoid <- betaMat0_2LD1[metabNamesNonFlavonoid, ]
betaMat0_2NonFlavonoid <- na.omit(betaMat0_2NonFlavonoid)
betaMat2_0NonFlavonoid <- betaMat2_0LD1[metabNamesNonFlavonoid, ]
betaMat2_0NonFlavonoid <- na.omit(betaMat2_0NonFlavonoid)
betaMat2_2NonFlavonoid <- betaMat2_2LD1[metabNamesNonFlavonoid, ]
betaMat2_2NonFlavonoid <- na.omit(betaMat2_2NonFlavonoid)


betaMat0_0NonFlavonoidChr10_42562665Sd <- sqrt(mean((betaMat0_0NonFlavonoid[, "Chr10_42562665"] - mean(betaMat0_0NonFlavonoid[, "Chr10_42562665"]))^2))
betaMat0_0NonFlavonoidChr17_16065902Sd <- sqrt(mean((betaMat0_0NonFlavonoid[, "Chr17_16065902"] - mean(betaMat0_0NonFlavonoid[, "Chr17_16065902"]))^2))

betaMat0_2NonFlavonoidChr10_42562665Sd <- sqrt(mean((betaMat0_2NonFlavonoid[, "Chr10_42562665"] - mean(betaMat0_2NonFlavonoid[, "Chr10_42562665"]))^2))
betaMat0_2NonFlavonoidChr17_16065902Sd <- sqrt(mean((betaMat0_2NonFlavonoid[, "Chr17_16065902"] - mean(betaMat0_2NonFlavonoid[, "Chr17_16065902"]))^2))

betaMat2_0NonFlavonoidChr10_42562665Sd <- sqrt(mean((betaMat2_0NonFlavonoid[, "Chr10_42562665"] - mean(betaMat2_0NonFlavonoid[, "Chr10_42562665"]))^2))
betaMat2_0NonFlavonoidChr17_16065902Sd <- sqrt(mean((betaMat2_0NonFlavonoid[, "Chr17_16065902"] - mean(betaMat2_0NonFlavonoid[, "Chr17_16065902"]))^2))

betaMat2_2NonFlavonoidChr10_42562665Sd <- sqrt(mean((betaMat2_2NonFlavonoid[, "Chr10_42562665"] - mean(betaMat2_2NonFlavonoid[, "Chr10_42562665"]))^2))
betaMat2_2NonFlavonoidChr17_16065902Sd <- sqrt(mean((betaMat2_2NonFlavonoid[, "Chr17_16065902"] - mean(betaMat2_2NonFlavonoid[, "Chr17_16065902"]))^2))

### thresSd = 2
thresSd <- 2

betaMat0_0NonFlavonoidChr10_42562665ThresSd2 <- thresSd * betaMat0_0NonFlavonoidChr10_42562665Sd
betaMat0_2NonFlavonoidChr10_42562665ThresSd2 <- thresSd * betaMat0_2NonFlavonoidChr10_42562665Sd
betaMat2_0NonFlavonoidChr10_42562665ThresSd2 <- thresSd * betaMat2_0NonFlavonoidChr10_42562665Sd
betaMat2_2NonFlavonoidChr10_42562665ThresSd2 <- thresSd * betaMat2_2NonFlavonoidChr10_42562665Sd

betaMat0_0NonFlavonoidChr17_16065902ThresSd2 <- thresSd * betaMat0_0NonFlavonoidChr17_16065902Sd
betaMat0_2NonFlavonoidChr17_16065902ThresSd2 <- thresSd * betaMat0_2NonFlavonoidChr17_16065902Sd
betaMat2_0NonFlavonoidChr17_16065902ThresSd2 <- thresSd * betaMat2_0NonFlavonoidChr17_16065902Sd
betaMat2_2NonFlavonoidChr17_16065902ThresSd2 <- thresSd * betaMat2_2NonFlavonoidChr17_16065902Sd


intervalThresSd2BetaMat0_0FlavonoidChr10_42562665 <- c( mean(betaMat0_0NonFlavonoid[, "Chr10_42562665"])- betaMat0_0NonFlavonoidChr10_42562665ThresSd2, mean(betaMat0_0NonFlavonoid[, "Chr10_42562665"])+ betaMat0_0NonFlavonoidChr10_42562665ThresSd2)
intervalThresSd2BetaMat0_2FlavonoidChr10_42562665 <- c( mean(betaMat0_2NonFlavonoid[, "Chr10_42562665"])- betaMat0_2NonFlavonoidChr10_42562665ThresSd2, mean(betaMat0_2NonFlavonoid[, "Chr10_42562665"])+ betaMat0_2NonFlavonoidChr10_42562665ThresSd2)
intervalThresSd2BetaMat2_0FlavonoidChr10_42562665 <- c( mean(betaMat2_0NonFlavonoid[, "Chr10_42562665"])- betaMat2_0NonFlavonoidChr10_42562665ThresSd2, mean(betaMat2_0NonFlavonoid[, "Chr10_42562665"])+ betaMat2_0NonFlavonoidChr10_42562665ThresSd2)
intervalThresSd2BetaMat2_2FlavonoidChr10_42562665 <- c( mean(betaMat2_2NonFlavonoid[, "Chr10_42562665"])- betaMat2_2NonFlavonoidChr10_42562665ThresSd2, mean(betaMat2_2NonFlavonoid[, "Chr10_42562665"])+ betaMat2_2NonFlavonoidChr10_42562665ThresSd2)

intervalThresSd2BetaMat0_0FlavonoidChr17_16065902 <- c( mean(betaMat0_0NonFlavonoid[, "Chr17_16065902"])- betaMat0_0NonFlavonoidChr17_16065902ThresSd2, mean(betaMat0_0NonFlavonoid[, "Chr17_16065902"])+ betaMat0_0NonFlavonoidChr17_16065902ThresSd2)
intervalThresSd2BetaMat0_2FlavonoidChr17_16065902 <- c( mean(betaMat0_2NonFlavonoid[, "Chr17_16065902"])- betaMat0_2NonFlavonoidChr17_16065902ThresSd2, mean(betaMat0_2NonFlavonoid[, "Chr17_16065902"])+ betaMat0_2NonFlavonoidChr17_16065902ThresSd2)
intervalThresSd2BetaMat2_0FlavonoidChr17_16065902 <- c( mean(betaMat2_0NonFlavonoid[, "Chr17_16065902"])- betaMat2_0NonFlavonoidChr17_16065902ThresSd2, mean(betaMat2_0NonFlavonoid[, "Chr17_16065902"])+ betaMat2_0NonFlavonoidChr17_16065902ThresSd2)
intervalThresSd2BetaMat2_2FlavonoidChr17_16065902 <- c( mean(betaMat2_2NonFlavonoid[, "Chr17_16065902"])- betaMat2_2NonFlavonoidChr17_16065902ThresSd2, mean(betaMat2_2NonFlavonoid[, "Chr17_16065902"])+ betaMat2_2NonFlavonoidChr17_16065902ThresSd2)


intervalThresSd2BetaMat0_0FlavonoidChr10_42562665
intervalThresSd2BetaMat0_2FlavonoidChr10_42562665
intervalThresSd2BetaMat2_0FlavonoidChr10_42562665
intervalThresSd2BetaMat2_2FlavonoidChr10_42562665
intervalThresSd2BetaMat0_0FlavonoidChr17_16065902
intervalThresSd2BetaMat0_2FlavonoidChr17_16065902
intervalThresSd2BetaMat2_0FlavonoidChr17_16065902
intervalThresSd2BetaMat2_2FlavonoidChr17_16065902


### thresSd = 3
thresSd <- 3

betaMat0_0NonFlavonoidChr10_42562665ThresSd3 <- thresSd * betaMat0_0NonFlavonoidChr10_42562665Sd
betaMat0_2NonFlavonoidChr10_42562665ThresSd3 <- thresSd * betaMat0_2NonFlavonoidChr10_42562665Sd
betaMat2_0NonFlavonoidChr10_42562665ThresSd3 <- thresSd * betaMat2_0NonFlavonoidChr10_42562665Sd
betaMat2_2NonFlavonoidChr10_42562665ThresSd3 <- thresSd * betaMat2_2NonFlavonoidChr10_42562665Sd

betaMat0_0NonFlavonoidChr17_16065902ThresSd3 <- thresSd * betaMat0_0NonFlavonoidChr17_16065902Sd
betaMat0_2NonFlavonoidChr17_16065902ThresSd3 <- thresSd * betaMat0_2NonFlavonoidChr17_16065902Sd
betaMat2_0NonFlavonoidChr17_16065902ThresSd3 <- thresSd * betaMat2_0NonFlavonoidChr17_16065902Sd
betaMat2_2NonFlavonoidChr17_16065902ThresSd3 <- thresSd * betaMat2_2NonFlavonoidChr17_16065902Sd


intervalThresSd3BetaMat0_0FlavonoidChr10_42562665 <- c( mean(betaMat0_0NonFlavonoid[, "Chr10_42562665"])- betaMat0_0NonFlavonoidChr10_42562665ThresSd3, mean(betaMat0_0NonFlavonoid[, "Chr10_42562665"])+ betaMat0_0NonFlavonoidChr10_42562665ThresSd3)
intervalThresSd3BetaMat0_2FlavonoidChr10_42562665 <- c( mean(betaMat0_2NonFlavonoid[, "Chr10_42562665"])- betaMat0_2NonFlavonoidChr10_42562665ThresSd3, mean(betaMat0_2NonFlavonoid[, "Chr10_42562665"])+ betaMat0_2NonFlavonoidChr10_42562665ThresSd3)
intervalThresSd3BetaMat2_0FlavonoidChr10_42562665 <- c( mean(betaMat2_0NonFlavonoid[, "Chr10_42562665"])- betaMat2_0NonFlavonoidChr10_42562665ThresSd3, mean(betaMat2_0NonFlavonoid[, "Chr10_42562665"])+ betaMat2_0NonFlavonoidChr10_42562665ThresSd3)
intervalThresSd3BetaMat2_2FlavonoidChr10_42562665 <- c( mean(betaMat2_2NonFlavonoid[, "Chr10_42562665"])- betaMat2_2NonFlavonoidChr10_42562665ThresSd3, mean(betaMat2_2NonFlavonoid[, "Chr10_42562665"])+ betaMat2_2NonFlavonoidChr10_42562665ThresSd3)

intervalThresSd3BetaMat0_0FlavonoidChr17_16065902 <- c( mean(betaMat0_0NonFlavonoid[, "Chr17_16065902"])- betaMat0_0NonFlavonoidChr17_16065902ThresSd3, mean(betaMat0_0NonFlavonoid[, "Chr17_16065902"])+ betaMat0_0NonFlavonoidChr17_16065902ThresSd3)
intervalThresSd3BetaMat0_2FlavonoidChr17_16065902 <- c( mean(betaMat0_2NonFlavonoid[, "Chr17_16065902"])- betaMat0_2NonFlavonoidChr17_16065902ThresSd3, mean(betaMat0_2NonFlavonoid[, "Chr17_16065902"])+ betaMat0_2NonFlavonoidChr17_16065902ThresSd3)
intervalThresSd3BetaMat2_0FlavonoidChr17_16065902 <- c( mean(betaMat2_0NonFlavonoid[, "Chr17_16065902"])- betaMat2_0NonFlavonoidChr17_16065902ThresSd3, mean(betaMat2_0NonFlavonoid[, "Chr17_16065902"])+ betaMat2_0NonFlavonoidChr17_16065902ThresSd3)
intervalThresSd3BetaMat2_2FlavonoidChr17_16065902 <- c( mean(betaMat2_2NonFlavonoid[, "Chr17_16065902"])- betaMat2_2NonFlavonoidChr17_16065902ThresSd3, mean(betaMat2_2NonFlavonoid[, "Chr17_16065902"])+ betaMat2_2NonFlavonoidChr17_16065902ThresSd3)


intervalThresSd3BetaMat0_0FlavonoidChr10_42562665
intervalThresSd3BetaMat0_2FlavonoidChr10_42562665
intervalThresSd3BetaMat2_0FlavonoidChr10_42562665
intervalThresSd3BetaMat2_2FlavonoidChr10_42562665
intervalThresSd3BetaMat0_0FlavonoidChr17_16065902
intervalThresSd3BetaMat0_2FlavonoidChr17_16065902
intervalThresSd3BetaMat2_0FlavonoidChr17_16065902
intervalThresSd3BetaMat2_2FlavonoidChr17_16065902



#PDF files with border lines
dirMidSTAMBase <- "midstream/"
dirMidSTAMGWAS <- paste0(dirMidSTAMBase, scriptID,
                         "_GWAS_for_three_markers/")


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"), row.names = 1)
metabNamesFlavonoid <- rownames(metabFlavonoid)

# table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("blue", nrow(betaMat0_0LD1))
names(colVec) <- rownames(betaMat0_0LD1)
colVec[metabNamesFlavonoid] <- "orange1"

betaMat0_0LD1Flavonoid <- betaMat0_0LD1[metabNamesFlavonoid, ]
betaMat0_2LD1Flavonoid <- betaMat0_2LD1[metabNamesFlavonoid, ]
betaMat2_0LD1Flavonoid <- betaMat2_0LD1[metabNamesFlavonoid, ]
betaMat2_2LD1Flavonoid <- betaMat2_2LD1[metabNamesFlavonoid, ]



### Chr10_42562665
markerID <- "Chr10_42562665"

# (0,0)_(0,2), thresSd = 2, 3, 2&3

thresSd <- 2
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(0,2)_thresSD_", thresSd, ".pdf"))
plot(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID],
     xlim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat0_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- 3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(0,2)_thresSD_", thresSd, ".pdf"))
plot(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID],
     xlim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat0_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- "2_3"
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(0,2)_thresSD_", thresSd, ".pdf"))
plot(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID],
     xlim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat0_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat0_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()


# (2,0)_(2,2), thresSd = 2, 3, 2&3
thresSd <- 2
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr10_42562665/", scriptID, "_(2,0)_(2,2)_thresSD_", thresSd, ".pdf"))
plot(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID],
     xlim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat2_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- 3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr10_42562665/", scriptID, "_(2,0)_(2,2)_thresSD_", thresSd, ".pdf"))
plot(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID],
     xlim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat2_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- "2_3"
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr10_42562665/", scriptID, "_(2,0)_(2,2)_thresSD_", thresSd, ".pdf"))
plot(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID],
     xlim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat2_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat2_0FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2_2FlavonoidChr10_42562665, col = 4, lty = 4, lwd = 1.5)
dev.off()




### Chr17_16065902; thresSd = 2, 3, 2&3
markerID <- "Chr17_16065902"

thresSd <- 2
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr17_16065902/", scriptID, "_(0_0)_(0_2)_", "_thresSD_", thresSd, ".pdf"))
plot(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID],
     xlim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat0_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- 3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr17_16065902/", scriptID, "_(0_0)_(0_2)_", "_thresSD_", thresSd, ".pdf"))
plot(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID],
     xlim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat0_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- "2_3"
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr17_16065902/", scriptID, "_(0_0)_(0_2)_", "_thresSD_", thresSd, ".pdf"))
plot(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID],
     xlim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, markerID], betaMat0_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 0, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat0_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat0_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat0_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat0_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()


thresSd <- 2
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr17_16065902/", scriptID, "_(2_0)_(2_2)_", "_thresSD_", thresSd, ".pdf"))
plot(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID],
     xlim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat2_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- 3
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr17_16065902/", scriptID, "_(2_0)_(2_2)_", "_thresSD_", thresSd, ".pdf"))
plot(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID],
     xlim = range(intervalThresSd3BetaMat2_0FlavonoidChr17_16065902),
     ylim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd3BetaMat2_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()

thresSd <- "2_3"
pdf(paste0(dirMidSTAMGWAS, "2.32_2017_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/", scriptID, "_Chr17_16065902/", scriptID, "_(2_0)_(2_2)_", "_thresSD_", thresSd, ".pdf"))
plot(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID],
     xlim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, markerID], betaMat2_2LD1[, markerID], na.rm = TRUE),
     main = paste0("Coefficient of ", markerID),
     xlab = "Chr06_18760995 = 2, Chr06_47490224 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47490224 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
abline(v = intervalThresSd2BetaMat2_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd2BetaMat2_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(v = intervalThresSd3BetaMat2_0FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
abline(h = intervalThresSd3BetaMat2_2FlavonoidChr17_16065902, col = 4, lty = 4, lwd = 1.5)
dev.off()




#### 2.3.2. Extract metabolites ####
betaMat0_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(0,0).csv", row.names = 1)
betaMat0_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(0,2).csv", row.names = 1)
betaMat2_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(2,0).csv", row.names = 1)
betaMat2_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/2.32_(2,2).csv", row.names = 1)

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
# metabNamesNonFlavonoid <- rownames(betaMat0_0LD1)[!(rownames(betaMat0_0LD1) %in% metabNamesFlavonoid)]

betaMat0_0LD1Flavonoid <- betaMat0_0LD1[metabNamesFlavonoid, ]
betaMat0_2LD1Flavonoid <- betaMat0_2LD1[metabNamesFlavonoid, ]
betaMat2_0LD1Flavonoid <- betaMat2_0LD1[metabNamesFlavonoid, ]
betaMat2_2LD1Flavonoid <- betaMat2_2LD1[metabNamesFlavonoid, ]


# #### Chr10_42562665
# ### (0,0),(0,2)
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 2, ])
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 1 & betaMat0_2LD1[, 1] < 0.5, ])
# rownames(betaMat0_2LD1)[betaMat0_0LD1[, 1] < -1 & betaMat0_2LD1[, 1] > -0.7]
#
# rownames(betaMat0_2LD1)[betaMat0_0LD1[, 1] > 1 & betaMat0_2LD1[, 1] > 1]
#
# ### (2,0),(2,2)
# rownames(betaMat2_0LD1[betaMat2_0LD1[, 1] > 2, ])
#
# ### (0,0),(2,0)
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 2, ])
#
# ### (0,2),(2,2)
# rownames(betaMat2_2LD1[betaMat2_2LD1[, 1] > 1, ])
# rownames(betaMat0_2LD1[betaMat0_2LD1[, 1] < -0.6, ])
#
#
#
# #### Chr17_16065902
# ### (0,0),(0,2)
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] < -0.4, ])
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] < -0.4 & betaMat0_2LD1[, 2] < -0.1, ])
# betaMat0_0LD1[betaMat0_0LD1[, 2] < -0.4 & betaMat0_2LD1[, 2] < -0.1, ]
# betaMat0_2LD1[betaMat0_0LD1[, 2] < -0.4 & betaMat0_2LD1[, 2] < -0.1, ]
#
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.6, ])
#
#
# ### (2,0),(2,2)
#
#
# ### (0,0),(2,0)
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.4 & betaMat2_0LD1[, 2] > -0.2, ])
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] < -0.4, ])
#
# ### (0,2),(2,2)
# rownames(betaMat0_2LD1[betaMat0_2LD1[, 2] > 0.24, ])
# rownames(betaMat2_2LD1[betaMat2_2LD1[, 2] > 0.4, ])



### Chr10_42562665
betaMat0_0LD1[order(betaMat0_0LD1$Chr10_42562665), ]
betaMat0_2LD1[order(betaMat0_2LD1$Chr10_42562665), ]
betaMat2_0LD1[order(betaMat2_0LD1$Chr10_42562665), ]
betaMat2_2LD1[order(betaMat2_2LD1$Chr10_42562665), ]


## (0,0),(0,2)
metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr10_42562665 <- rownames(betaMat0_0LD1Flavonoid)[
  betaMat0_0LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat0_0FlavonoidChr10_42562665[1]
  # & betaMat0_2LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[1]
  # & betaMat0_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[2]
]
metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr10_42562665 <- rownames(betaMat0_0LD1Flavonoid)[
  betaMat0_0LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat0_0FlavonoidChr10_42562665[2]
  # & betaMat0_2LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[1]
  # & betaMat0_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[2]
]

metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr10_42562665 <- rownames(betaMat0_2LD1Flavonoid)[
  betaMat0_0LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat0_0FlavonoidChr10_42562665[2] &
  betaMat0_0LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[1] &
  betaMat0_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[1]
]

metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr10_42562665 <- rownames(betaMat0_2LD1Flavonoid)[
  # betaMat0_0LD1Flavonoid[, 2] < intervalThresSd3BetaMat0_0FlavonoidChr10_42562665[2] &
  # betaMat0_0LD1Flavonoid[, 2] > intervalThresSd3BetaMat0_0FlavonoidChr10_42562665[1] &
  betaMat0_2LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat0_2FlavonoidChr10_42562665[2]
]

# Around y = x line, significant flavonoid
rownames(betaMat0_2LD1Flavonoid)[betaMat0_0LD1Flavonoid[, 1] > 1 & betaMat0_2LD1Flavonoid[, 1] > 1]

### Non-Flavonoid related
rownames(betaMat0_0LD1)[
  betaMat0_0LD1[, "Chr10_42562665"] > intervalThresSd3BetaMat0_0FlavonoidChr10_42562665[2]
]



## (2,0),(2,2)
metabNamesFlavonoidLessThanThresSd3BetaMat2_0Chr10_42562665 <- rownames(betaMat2_0LD1Flavonoid)[
  betaMat2_0LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat2_0FlavonoidChr10_42562665[1]
  # & betaMat2_2LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[1]
  # & betaMat2_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[2]
]
metabNamesFlavonoidMoreThanThresSd3BetaMat2_0Chr10_42562665 <- rownames(betaMat2_0LD1Flavonoid)[
  betaMat2_0LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat2_0FlavonoidChr10_42562665[2]
  # & betaMat2_2LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[1]
  # & betaMat2_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[2]
]

metabNamesFlavonoidLessThanThresSd3BetaMat2_2Chr10_42562665 <- rownames(betaMat2_2LD1Flavonoid)[
  betaMat2_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[1]
  # & betaMat2_0LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[1]
  # & betaMat2_2LD1Flavonoid[, "Chr10_42562665"] < intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[1]
]

metabNamesFlavonoidMoreThanThresSd3BetaMat2_2Chr10_42562665 <- rownames(betaMat2_2LD1Flavonoid)[
  # betaMat2_0LD1Flavonoid[, 2] < intervalThresSd3BetaMat2_0FlavonoidChr10_42562665[2] &
  # betaMat2_0LD1Flavonoid[, 2] > intervalThresSd3BetaMat2_0FlavonoidChr10_42562665[1] &
  betaMat2_2LD1Flavonoid[, "Chr10_42562665"] > intervalThresSd3BetaMat2_2FlavonoidChr10_42562665[2]
]



### Chr17_16065902
# rownames(betaMat0[betaMat0[, 3] > 0.25, ])
# rownames(betaMat2)[betaMat2[, 3] < -0.3]


betaMat0_0LD1[order(betaMat0_0LD1$Chr17_16065902), ]
betaMat0_2LD1[order(betaMat0_2LD1$Chr17_16065902), ]
betaMat2_0LD1[order(betaMat2_0LD1$Chr17_16065902), ]
betaMat2_2LD1[order(betaMat2_2LD1$Chr17_16065902), ]



metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr17_16065902 <- rownames(betaMat0_0LD1Flavonoid)[
  betaMat0_0LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[1] &
    betaMat0_2LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[1] &
    betaMat0_2LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[2]
]

metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr17_16065902 <- rownames(betaMat0_0LD1Flavonoid)[
  betaMat0_0LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[2]
  # & betaMat0_2LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[1]
  # & betaMat0_2LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[2]
]

metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr17_16065902 <- rownames(betaMat0_2LD1Flavonoid)[
  # betaMat0_0LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[2] &
  # betaMat0_0LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[1] &
  betaMat0_2LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[1]
]

metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902 <- rownames(betaMat0_2LD1Flavonoid)[
  (betaMat0_2LD1Flavonoid[, "Chr17_16065902"] > betaMat0_0LD1Flavonoid[, "Chr17_16065902"]) &
  betaMat0_2LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[2]
]
# metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902 <- rownames(betaMat0_2LD1)[
#   betaMat0_0LD1[, "Chr17_16065902"] < intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[2] &
#   betaMat0_0LD1[, "Chr17_16065902"] > intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[1] &
#   betaMat0_2LD1[, "Chr17_16065902"] > intervalThresSd3BetaMat0_2FlavonoidChr17_16065902[2]
# ]


metabNamesFlavonoidLessThanThresSd3BetaMat2_0Chr17_16065902 <- rownames(betaMat2_0LD1Flavonoid)[
  betaMat2_0LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat2_0FlavonoidChr17_16065902[1]
]

metabNamesFlavonoidMoreThanThresSd3BetaMat2_0Chr17_16065902 <- rownames(betaMat2_0LD1Flavonoid)[
  betaMat2_0LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat2_0FlavonoidChr17_16065902[2]
]

metabNamesFlavonoidLessThanThresSd3BetaMat2_2Chr17_16065902 <- rownames(betaMat2_2LD1Flavonoid)[
  betaMat2_2LD1Flavonoid[, "Chr17_16065902"] < intervalThresSd3BetaMat2_2FlavonoidChr17_16065902[1]
]

metabNamesFlavonoidMoreThanThresSd3BetaMat2_2Chr17_16065902 <- rownames(betaMat2_2LD1Flavonoid)[
    betaMat2_2LD1Flavonoid[, "Chr17_16065902"] > intervalThresSd3BetaMat2_2FlavonoidChr17_16065902[2]
]

### Non-flavonoid related
rownames(betaMat0_0LD1)[
  betaMat0_0LD1[, "Chr17_16065902"] > intervalThresSd3BetaMat0_0FlavonoidChr17_16065902[2]
]

rownames(betaMat2_0LD1)[
  betaMat2_0LD1[, "Chr17_16065902"] < intervalThresSd3BetaMat2_0FlavonoidChr17_16065902[1]
]



#
metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr10_42562665
metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr10_42562665
metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr10_42562665
metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr10_42562665

metabNamesFlavonoidLessThanThresSd3BetaMat2_0Chr10_42562665
metabNamesFlavonoidMoreThanThresSd3BetaMat2_0Chr10_42562665
metabNamesFlavonoidLessThanThresSd3BetaMat2_2Chr10_42562665
metabNamesFlavonoidMoreThanThresSd3BetaMat2_2Chr10_42562665


metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr17_16065902
metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr17_16065902
metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr17_16065902
metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902

metabNamesFlavonoidLessThanThresSd3BetaMat2_0Chr17_16065902
metabNamesFlavonoidMoreThanThresSd3BetaMat2_0Chr17_16065902
metabNamesFlavonoidLessThanThresSd3BetaMat2_2Chr17_16065902
metabNamesFlavonoidMoreThanThresSd3BetaMat2_2Chr17_16065902



metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr10_42562665 <- c(metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr10_42562665, metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr10_42562665, metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr10_42562665, metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr10_42562665)
metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr10_42562665 <- na.omit(metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr10_42562665)


metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr17_16065902 <- c(metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr17_16065902, metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr17_16065902, metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr17_16065902, metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902)
metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr17_16065902 <- na.omit(metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr17_16065902)


dirMetabCoefficient <- paste("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/")
write.csv(x = metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr10_42562665, file = paste0(dirMetabCoefficient, "2.32_metabNames_Flavonoid_All_Detected_ThresSd3_given_0_0_And_0_2_Chr10_42562665.csv"))
write.csv(x = metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr17_16065902, file = paste0(dirMetabCoefficient, "2.32_metabNames_Flavonoid_All_Detected_ThresSd3_given_0_0_And_0_2_Chr17_16065902.csv"))



#### Check annotation of extracted metabolites ####
metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"), row.names = 1)
# metabNamesFlavonoid <- rownames(metabFlavonoid)


annotationMetabFlavonoidLessThanThresSd3BetaMat0_0Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr10_42562665, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_0Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr10_42562665, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat0_2Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr10_42562665, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_2Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr10_42562665, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat2_0Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2_0Chr10_42562665, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat2_0Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2_0Chr10_42562665, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat2_2Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2_2Chr10_42562665, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat2_2Chr10_42562665 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2_2Chr10_42562665, ]

annotationMetabFlavonoidLessThanThresSd3BetaMat0_0Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0_0Chr17_16065902, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_0Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0_0Chr17_16065902, ]
annotationMetabFlavonoidLessThanThresSd3BetaMat0_2Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat0_2Chr17_16065902, ]
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902, ]
# annotationMetabFlavonoidLessThanThresSd3BetaMat2_0Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2_0Chr17_16065902, ]
# annotationMetabFlavonoidMoreThanThresSd3BetaMat2_0Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2_0Chr17_16065902, ]
# annotationMetabFlavonoidLessThanThresSd3BetaMat2_2Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidLessThanThresSd3BetaMat2_2Chr17_16065902, ]
# annotationMetabFlavonoidMoreThanThresSd3BetaMat2_2Chr17_16065902 <- metabFlavonoid[metabNamesFlavonoidMoreThanThresSd3BetaMat2_2Chr17_16065902, ]


annotationMetabFlavonoidLessThanThresSd3BetaMat0_0Chr10_42562665
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_0Chr10_42562665
annotationMetabFlavonoidLessThanThresSd3BetaMat0_2Chr10_42562665
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_2Chr10_42562665
annotationMetabFlavonoidLessThanThresSd3BetaMat2_0Chr10_42562665
annotationMetabFlavonoidMoreThanThresSd3BetaMat2_0Chr10_42562665
annotationMetabFlavonoidLessThanThresSd3BetaMat2_2Chr10_42562665
annotationMetabFlavonoidMoreThanThresSd3BetaMat2_2Chr10_42562665

annotationMetabFlavonoidLessThanThresSd3BetaMat0_0Chr17_16065902
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_0Chr17_16065902
annotationMetabFlavonoidLessThanThresSd3BetaMat0_2Chr17_16065902
annotationMetabFlavonoidMoreThanThresSd3BetaMat0_2Chr17_16065902
# annotationMetabFlavonoidLessThanThresSd3BetaMat2_0Chr17_16065902
# annotationMetabFlavonoidMoreThanThresSd3BetaMat2_0Chr17_16065902
# annotationMetabFlavonoidLessThanThresSd3BetaMat2_2Chr17_16065902
# annotationMetabFlavonoidMoreThanThresSd3BetaMat2_2Chr17_16065902





##### 2.4. Read genotypic values of Metabolomic data in 2017 of (0,0),(0,2),(2,0),(2,2) for Chr06_18760995 and Chr10_42562665 into R #####
#### 2.4.1 For each maker ####
cultivationInfo <- "2017"
targetInfo <- "Metabolome"

targetType <- c("0_0", "0_2", "2_0", "2_2")
# targetType <- c("0_2", "2_0", "2_2")

markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr10_42562665"



gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
See(gvMetab2017Total)


gvMetab2017Total <- scale(gvMetab2017Total, center = TRUE, scale = TRUE)
See(gvMetab2017Total)
gvMetab2017Total <- as.data.frame(gvMetab2017Total)

lineNames <- rownames(gvMetab2017Total)
markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]


markerInterest1DF <- as.data.frame(markerInterest1)
markerInterest2DF <- as.data.frame(markerInterest2)
markerInterest1And2DF <- cbind(markerInterest1DF, markerInterest2DF)
See(markerInterest1And2DF)

# markerVarietyNames <- names(markerInterest)

table(markerInterest1And2DF$markerInterest1)
table(markerInterest1And2DF$markerInterest2)
marker0_0VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 0]
marker0_2VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 2]
marker2_0VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 0]
marker2_2VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 2]

See(marker0_0VarietyNames)
See(marker0_2VarietyNames)
See(marker2_0VarietyNames)
See(marker2_2VarietyNames)


lineNames0_0 <- lineNames[(lineNames %in% marker0_0VarietyNames)]
commonNames0_0 <- marker0_0VarietyNames[(marker0_0VarietyNames %in% lineNames0_0)]
lineNames0_2 <- lineNames[(lineNames %in% marker0_2VarietyNames)]
commonNames0_2 <- marker0_2VarietyNames[(marker0_2VarietyNames %in% lineNames0_2)]
lineNames2_0 <- lineNames[(lineNames %in% marker2_0VarietyNames)]
commonNames2_0 <- marker2_0VarietyNames[(marker2_0VarietyNames %in% lineNames2_0)]
lineNames2_2 <- lineNames[(lineNames %in% marker2_2VarietyNames)]
commonNames2_2 <- marker2_2VarietyNames[(marker2_2VarietyNames %in% lineNames2_2)]

gvMetab2017Total0_0 <- gvMetab2017Total[commonNames0_0, ]
gvMetab2017Total0_2 <- gvMetab2017Total[commonNames0_2, ]
gvMetab2017Total2_0 <- gvMetab2017Total[commonNames2_0, ]
gvMetab2017Total2_2 <- gvMetab2017Total[commonNames2_2, ]
# gvMetab2017Total0 <- na.omit(gvMetab2017Total0)
See(gvMetab2017Total0_0, rown = 20)
See(gvMetab2017Total0_2)
See(gvMetab2017Total2_0)
See(gvMetab2017Total2_2)




### LD1
markerNames <- c("Chr06_47490224", "Chr17_16065902")
betaMat0_0LD1 <- matrix(NA, nrow = length(gvMetab2017Total0_0), ncol = length(markerNames))
betaMat0_2LD1 <- matrix(NA, nrow = length(gvMetab2017Total0_2), ncol = length(markerNames))
betaMat2_0LD1 <- matrix(NA, nrow = length(gvMetab2017Total2_0), ncol = length(markerNames))
betaMat2_2LD1 <- matrix(NA, nrow = length(gvMetab2017Total2_2), ncol = length(markerNames))
rownames(betaMat0_0LD1) <- colnames(gvMetab2017Total0_0)
rownames(betaMat0_2LD1) <- colnames(gvMetab2017Total0_2)
rownames(betaMat2_0LD1) <- colnames(gvMetab2017Total2_0)
rownames(betaMat2_2LD1) <- colnames(gvMetab2017Total2_2)
colnames(betaMat0_0LD1) <- markerNames
colnames(betaMat0_2LD1) <- markerNames
colnames(betaMat2_0LD1) <- markerNames
colnames(betaMat2_2LD1) <- markerNames
See(betaMat0_0LD1)
See(betaMat0_2LD1)
See(betaMat2_0LD1)
See(betaMat2_2LD1)

targetTypeNo <- 1
traitNo <- 1

thresLD <- 2
for (targetTypeNo in 1:length(targetType)) {

  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- eval(parse( text = paste0("gvMetab2017Total", targetTypeNow)))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]

  See(gvMetab)

  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))



  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf >= 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)


  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)

  dirMarkerTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             markerInterestID1,
                             "_",
                             markerInterestID2,
                             "_results/")
  dir.create(dirMarkerTypeNow)

  dirTargetTypeNow <- paste0(dirMarkerTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  dirLDNow <- paste0(dirTargetTypeNow,
                     scriptID, "_",
                     "LD=1",
                     "_results/")
  dir.create(dirLDNow)


  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
      # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)

      KNow <- GRM(gastonDataNow)
      eigenKNow <- eigen(KNow)
    } else {
      gastonDataNow <- gastonData

      KNow <- K
      eigenKNow <- eigenK
    }
    st <- Sys.time()
    gastonRes <- association.test(x = gastonDataNow, method = "lmm",
                                  response = "quantitative",
                                  test = "wald", eigenK = eigenKNow,
                                  p = 2)
    end <- Sys.time()
    print(end - st)

    if (targetTypeNow == "0_0"){
      betaMat0_0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "0_2") {
      betaMat0_2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "2_0") {
      betaMat2_0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "2_2") {
      betaMat2_2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
      # gastonRes$id[which(gastonRes$id %in% markerNames)]
      # table(gastonRes$id %in% "Chr17_16065902")
      # table(gastonRes$id %in% markerNames)
      # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]

    }


    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirLDNow, scriptID,
                       "_", traitNow, "/")
    dir.create(dirSave0)

    dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
    dir.create(dirSave)

    fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
                           thresLD, "_gaston_wald")

    fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
    write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)

    fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
    save(gastonRes, file = fileNameSaveRData)

    fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
    fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")

    png(fileNameSaveManhattan, width = 1200, height = 900)
    gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    dev.off()

    png(fileNameSaveQq, width = 900, height = 900)
    RAINBOWR::qq(- log10(gastonRes$p))
    dev.off()


    if (thresLD <= 0.4) {
      gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
                                            BP = gastonRes$pos,
                                            P = gastonRes$p,
                                            SNP = gastonRes$id,
                                            BLOCK = gastonDataNow@snps$block)
      sigSNPs <- gastonRes$id[pAdj < 0.05]

      fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
      fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")

      plotlyManhattan <- manhattanly(gastonResForManhattanly,
                                     snp = "SNP", gene = "BLOCK",
                                     highlight = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
                                               basename(fileNameSavePlotlyManhattan)))

      plotlyQq <- qqPlotly(data = gastonResForManhattanly,
                           highlightSNPNames = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
                                               basename(fileNameSavePlotlyQq)))
    }
    print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
  }
}


dirTwoMarkers2017 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr10_42562665/")
dir.create(dirTwoMarkers2017)

fileName0_0 <- paste0(dirTwoMarkers2017, scriptID, "_(0,0).csv")
fileName0_2 <- paste0(dirTwoMarkers2017, scriptID, "_(0,2).csv")
fileName2_0 <- paste0(dirTwoMarkers2017, scriptID, "_(2,0).csv")
fileName2_2 <- paste0(dirTwoMarkers2017, scriptID, "_(2,2).csv")

write.csv(x = betaMat0_0LD1, file = fileName0_0)
write.csv(x = betaMat0_2LD1, file = fileName0_2)
write.csv(x = betaMat2_0LD1, file = fileName2_0)
write.csv(x = betaMat2_2LD1, file = fileName2_2)


### comparison of coefficients
dirTwoMarkers2017 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr10_42562665/")
dir.create(dirTwoMarkers2017)

betaMat0_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,0).csv"), row.names = 1)
betaMat0_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,2).csv"), row.names = 1)
betaMat2_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,0).csv"), row.names = 1)
betaMat2_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,2).csv"), row.names = 1)


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]

table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("blue1", nrow(betaMat0))
names(colVec) <- rownames(betaMat0)
colVec[metabNamesFlavonoid] <- "orange1"



dir.create(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/"))
dir.create(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/"))

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_(0,0)_(0,2).pdf"))
plot(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 0, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_(0,0)_(2,0).pdf"))
plot(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_(0,0)_(2,2).pdf"))
plot(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_(0,2)_(2,0).pdf"))
plot(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1],
     xlim = range(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 2",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_(0,2)_(2,2).pdf"))
plot(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 2",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_(2,0)_(2,2).pdf"))
plot(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 2, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()



pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,0)_(0,2).pdf"))
plot(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 0, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,0)_(2,0).pdf"))
plot(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,0)_(2,2).pdf"))
plot(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,2)_(2,0).pdf"))
plot(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2],
     xlim = range(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 2",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(0,2)_(2,2).pdf"))
plot(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr10_42562665 = 2",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_(2,0)_(2,2).pdf"))
plot(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 2, Chr10_42562665 = 0",
     ylab = "Chr06_18760995 = 2, Chr10_42562665 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()



### pairs plot
Chr06_47490224DF <- data.frame("0_0" = betaMat0_0LD1[, 1], "0_2" = betaMat0_2LD1[, 1], "2_0" = betaMat2_0LD1[, 1], "2_2" = betaMat2_2LD1[, 1])
colnames(Chr06_47490224DF) <- c("0_0", "0_2", "2_0", "2_2")
Chr17_16065902DF <- data.frame("0_0" = betaMat0_0LD1[, 2], "0_2" = betaMat0_2LD1[, 2], "2_0" = betaMat2_0LD1[, 2], "2_2" = betaMat2_2LD1[, 2])
colnames(Chr17_16065902DF) <- c("0_0", "0_2", "2_0", "2_2")


# upper <- function(x, y, ...){
#   oldpar <- par(usr = c(0, 1, 0, 1))
#   v <- abs(cor(x, y))
#   sz <- 1 + v * 2
#   text(0.5, 0.5, sprintf("%.3f", v), cex = sz)
#   par(oldpar)
# }

lower <- function(x, y, ...){
  points(x, y, pch = 21, bg = c("blue1", "orange1")[as.numeric(as.factor(colVec))])
  if( 1==1 ){
    abline(0, 1, col = 2, lty = 2, lwd = 1.5)
  }
}

# pairs(Chr06_47490224DF, upper.panel=upper, lower.panel=lower)
pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr06_47490224/", scriptID, "_pairs_plot.pdf"))
pairs(Chr06_47490224DF, lower.panel=lower, oma = c(3,3,3,15))
par(xpd = TRUE)
legend(x = "bottomright", legend = c("Non-flavonoid", "Flavonoid"), col = c("blue1", "orange1"), pch = 19, cex = 1)
dev.off()

pdf(paste0(dirTwoMarkers2017, scriptID, "_Chr17_16065902/", scriptID, "_pairs_plot.pdf"))
pairs(Chr17_16065902DF, lower.panel=lower, oma = c(3,3,3,15))
par(xpd = TRUE)
legend(x = "bottomright", legend = c("Non-flavonoid", "Flavonoid"), col = c("blue1", "orange1"), pch = 19)
dev.off()



# Chr06_47490224DF <- data.frame("0_0" = betaMat0_0LD1[, 1], "0_2" = betaMat0_2LD1[, 1], "2_0" = betaMat2_0LD1[, 1], "2_2" = betaMat2_2LD1[, 1])
# colnames(Chr06_47490224DF) <- c("0_0", "0_2", "2_0", "2_2")
# Chr17_16065902DF <- data.frame("0_0" = betaMat0_0LD1[, 2], "0_2" = betaMat0_2LD1[, 2], "2_0" = betaMat2_0LD1[, 2], "2_2" = betaMat2_2LD1[, 2])
# colnames(Chr17_16065902DF) <- c("0_0", "0_2", "2_0", "2_2")
#
# See(Chr06_47490224DF)
# See(Chr17_16065902DF)


# ggpairs(Chr06_47490224DF, aes_string(colour="colVec", alpha=0.5)) + geom_abline(slope = 1, intercept = 0)
# ggpairs(Chr17_16065902DF, aes_string(colour="colVec", alpha=0.5))
#
#
# ggpairs(Chr06_47490224DF, mapping = aes(color = colVec), lower = list(continuous = "smooth"))
# ggpairs(Chr17_16065902DF, mapping = aes(color = colVec), lower = list(continuous = "smooth"))



### extract metabolites
betaMat0_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr10_42562665/2.32_(0,0).csv", row.names = 1)
betaMat0_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr10_42562665/2.32_(0,2).csv", row.names = 1)
betaMat2_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr10_42562665/2.32_(2,0).csv", row.names = 1)
betaMat2_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr10_42562665/2.32_(2,2).csv", row.names = 1)


#### Chr06_47490224
### (0,0),(0,2)
chr06 <- arrange(betaMat0_2LD1, Chr06_47490224)
rownames(betaMat0_2LD1[betaMat0_2LD1[, 1] < -1, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 2, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 1, ])
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 0.8 & betaMat0_2LD1[, 1] < 0.8 , ])

### (2,0),(2,2)
rownames(betaMat2_0LD1[betaMat2_0LD1[, 1] > 2, ])
rownames(betaMat2_0LD1[betaMat2_0LD1[, 1] < -2, ])

### (0,0),(2,0)
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 2, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] < -0.6, ])

### (0,2),(2,2)
rownames(betaMat0_2LD1[betaMat0_2LD1[, 1] > 0.45, ])
rownames(betaMat2_2LD1[betaMat2_2LD1[, 1] > 1, ])
rownames(betaMat0_2LD1[betaMat0_2LD1[, 1] < -0.6, ])
rownames(betaMat2_2LD1[betaMat2_2LD1[, 1] < -1, ])


#### Chr17_16065902
chr17beta00 <- betaMat0_0LD1[, 2]
chr17beta02 <- betaMat0_2LD1[, 2]
chr17beta00beta02 <- cbind(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2])
rownames(chr17beta00beta02) <- rownames(betaMat0_0LD1)
colnames(chr17beta00beta02) <- c("beta0_0", "beta0_2")
See(chr17beta00beta02)
chr17beta00beta02DF <- as.data.frame(chr17beta00beta02)
arrange(chr17beta00beta02DF, beta0_0)


### (0,0),(0,2)
chr17beta00 <- arrange(betaMat0_0LD1, Chr17_16065902)
rownames(betaMat0_2LD1[betaMat0_2LD1[, 2] < -0.6, ])
rownames(betaMat0_2LD1)[betaMat0_2LD1[, 2] < -0.6]
rownames(betaMat0_2LD1)[betaMat0_2LD1[, 2] < 0.1 & betaMat0_0LD1[, 2] > 0.39]
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.6, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.4, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.5, ])

### (2,0),(2,2)


### (0,0),(2,0)
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.4 & betaMat2_0LD1[, 2] > -0.2, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] < -0.4, ])

### (0,2),(2,2)
rownames(betaMat0_2LD1[betaMat0_2LD1[, 2] > 0.24, ])
rownames(betaMat2_2LD1[betaMat2_2LD1[, 2] > 0.4, ])










###### 3. Perform gaston for Metabolomic data in 2018 ######
##### 3.1. Read marker genotype into R #####
gastonData0 <- gaston::read.vcf(file = "raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz")
lineNamesGeno <- gastonData0@ped$id

haploBlockList <- data.frame(fread("raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2._haplotype_block_list.csv"))

chrNos <- gastonData0@snps$chr
chrNames <- sprintf(fmt = paste0("Chr",
                                 "%0", floor(log10(max(chrNos))) + 1, "i"),
                    chrNos)
pos <- gastonData0@snps$pos
mrkNames <- paste0(chrNames, "_", pos)
gastonData0@snps$id <- mrkNames

blockNames0 <- haploBlockList$block
names(blockNames0) <- haploBlockList$marker
blockNames <- blockNames0[mrkNames]
gastonData0@snps$block <- blockNames

gastonDataSmall <- LD.thin(gastonData0, threshold = thresLD)




##### 3.2. Read genotypic values of Metabolomic data in 2018 of 0 or 2 for Chr06_18760995 into R #####
#### 3.2.1 For each maker ####
cultivationInfo <- "2018"
targetInfo <- "Metabolome"

targetType <- c("0", "2")

markerInterestID <- "Chr06_18760995"



gvMetab2018Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2018.csv", row.names = 1)
See(gvMetab2018Total)
# rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"


gvMetab2018Total <- scale(gvMetab2018Total, center = TRUE, scale = TRUE)
See(gvMetab2018Total)
gvMetab2018Total <- as.data.frame(gvMetab2018Total)

lineNames <- rownames(gvMetab2018Total)
markerInterest <- gastonData0Matrix[, markerInterestID]


markerInterestDF <- as.data.frame(markerInterest)


# markerVarietyNames <- names(markerInterest)

table(markerInterestDF$markerInterest)
marker0VarietyNames <- rownames(markerInterestDF)[markerInterestDF$markerInterest == 0]
marker2VarietyNames <- rownames(markerInterestDF)[markerInterestDF$markerInterest == 2]
See(marker0VarietyNames)
See(marker2VarietyNames)


lineNames0 <- lineNames[(lineNames%in%marker0VarietyNames)]
commonNames0 <- marker0VarietyNames[(marker0VarietyNames%in%lineNames0)]
lineNames2 <- lineNames[(lineNames%in%marker2VarietyNames)]
commonNames2 <- marker2VarietyNames[(marker2VarietyNames%in%lineNames2)]

gvMetab2018Total0 <- gvMetab2018Total[commonNames0, ]
gvMetab2018Total2 <- gvMetab2018Total[commonNames2, ]
# gvMetab2018Total0 <- na.omit(gvMetab2018Total0)
See(gvMetab2018Total0, rown = 20)
See(gvMetab2018Total2)




### LD1
markerNames <- c("Chr06_47428824", "Chr10_42562665", "Chr17_15930807")
betaMat0LD1 <- matrix(NA, nrow = length(gvMetab2018Total0), ncol = length(markerNames))
betaMat2LD1 <- matrix(NA, nrow = length(gvMetab2018Total2), ncol = length(markerNames))
rownames(betaMat0LD1) <- colnames(gvMetab2018Total0)
rownames(betaMat2LD1) <- colnames(gvMetab2018Total2)
colnames(betaMat0LD1) <- markerNames
colnames(betaMat2LD1) <- markerNames
See(betaMat0LD1)
See(betaMat2LD1)

# targetTypeNo <- 1
# traitNo <- 1

thresLD <- 2
for (targetTypeNo in 1:length(targetType)) {

  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- eval(parse( text = paste0("gvMetab2018Total", targetTypeNow)))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]

  See(gvMetab)

  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))



  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf >= 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)


  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)

  dirMarkerTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             markerInterestID,
                             "_results/")
  dir.create(dirMarkerTypeNow)

  dirTargetTypeNow <- paste0(dirMarkerTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  dirLDNow <- paste0(dirTargetTypeNow,
                     scriptID, "_",
                     "LD=1",
                     "_results/")
  dir.create(dirLDNow)


  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
      # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)

      KNow <- GRM(gastonDataNow)
      eigenKNow <- eigen(KNow)
    } else {
      gastonDataNow <- gastonData

      KNow <- K
      eigenKNow <- eigenK
    }
    st <- Sys.time()
    gastonRes <- association.test(x = gastonDataNow, method = "lmm",
                                  response = "quantitative",
                                  test = "wald", eigenK = eigenKNow,
                                  p = 2)
    end <- Sys.time()
    print(end - st)

    if (targetTypeNow == "0"){
      betaMat0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else {
      betaMat2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
      gastonRes$id[which(gastonRes$id %in% markerNames)]
      # table(gastonRes$id %in% "Chr17_16065902")
      # table(gastonRes$id %in% markerNames)
      # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]

    }


    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirLDNow, scriptID,
                       "_", traitNow, "/")
    dir.create(dirSave0)

    dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
    dir.create(dirSave)

    fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
                           thresLD, "_gaston_wald")

    fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
    write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)

    fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
    save(gastonRes, file = fileNameSaveRData)

    fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
    fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")

    png(fileNameSaveManhattan, width = 1200, height = 900)
    gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    dev.off()

    png(fileNameSaveQq, width = 900, height = 900)
    RAINBOWR::qq(- log10(gastonRes$p))
    dev.off()


    if (thresLD <= 0.4) {
      gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
                                            BP = gastonRes$pos,
                                            P = gastonRes$p,
                                            SNP = gastonRes$id,
                                            BLOCK = gastonDataNow@snps$block)
      sigSNPs <- gastonRes$id[pAdj < 0.05]

      fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
      fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")

      plotlyManhattan <- manhattanly(gastonResForManhattanly,
                                     snp = "SNP", gene = "BLOCK",
                                     highlight = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
                                               basename(fileNameSavePlotlyManhattan)))

      plotlyQq <- qqPlotly(data = gastonResForManhattanly,
                           highlightSNPNames = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
                                               basename(fileNameSavePlotlyQq)))
    }
    print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
  }
}



dir.create(paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/"))


fileName0 <- paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_0.csv")
fileName2 <- paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID ,"_Chr06_18760995_2.csv")


write.csv(x = betaMat0LD1, file = fileName0)
write.csv(x = betaMat2LD1, file = fileName2)











#### comparison of coefficients
betaMat0 <- read.csv(paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_0.csv"), row.names = 1)
betaMat2 <- read.csv(paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_18760995_2.csv"), row.names = 1)
See(betaMat0)
See(betaMat2)


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]

table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("green4", nrow(betaMat0))
names(colVec) <- rownames(betaMat0)
colVec[metabNamesFlavonoid] <- "orange1"



pdf(paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr06_47428824.pdf"))
plot(betaMat0[, 1], betaMat2[, 1],
     xlim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     ylim = range(betaMat0[, 1], betaMat2[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4,
     cex.axis = 1.3)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr10_42562665.pdf"))
plot(betaMat0[, 2], betaMat2[, 2],
     xlim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     ylim = range(betaMat0[, 2], betaMat2[, 2], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirMidSTAMGWAS, "2.32_2018_Metabolome_results/", scriptID, "_Coefficient_LD1_scaled_background_Chr06_18760995/", scriptID, "_Chr17_15930807.pdf"))
plot(betaMat0[, 3], betaMat2[, 3],
     xlim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     ylim = range(betaMat0[, 3], betaMat2[, 3], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()


rownames(betaMat2)[betaMat2[, 2] >= 60]


pairs(cbind(betaMat0[, 2], betaMat2[, 2]))


#### Extract metabolites
### Chr06_47428824
rownames(betaMat0[betaMat0[, 1] < -0.4 & betaMat2[, 1] < 0.1, ])
rownames(betaMat0[betaMat0[, 1] > 0.5, ])
rownames(betaMat2[betaMat2[, 1] < -1.2, ])

### Chr10_42562665
rownames(betaMat0[betaMat0[, 2] > 0.6, ])
rownames(betaMat2[betaMat2[, 2] > 1, ])
rownames(betaMat0[betaMat0[, 2] < -0.7, ])
rownames(betaMat2[betaMat2[, 2] < -1, ])

### Chr17_15930807
rownames(betaMat0[betaMat0[, 3] > 0.25, ])
rownames(betaMat2)[betaMat2[, 3] < -0.3]




#### 3.2.2. For two marker combinations ( Chr06_18760995, Chr06_47428824 )
# markerNames <- c("Chr06_47428824", "Chr10_42562665", "Chr17_15930807")
# betaMat0LD1 <- matrix(NA, nrow = length(gvMetab2018Total0), ncol = length(markerNames))
# betaMat2LD1 <- matrix(NA, nrow = length(gvMetab2018Total2), ncol = length(markerNames))
# rownames(betaMat0LD1) <- colnames(gvMetab2018Total0)
# rownames(betaMat2LD1) <- colnames(gvMetab2018Total2)
# colnames(betaMat0LD1) <- markerNames
# colnames(betaMat2LD1) <- markerNames
# See(betaMat0LD1)
# See(betaMat2LD1)
#
# targetTypeNo <- 1
# traitNo <- 1
#
# thresLD <- 2
# for (targetTypeNo in 1:length(targetType)) {
#
#   targetTypeNow <- targetType[targetTypeNo]
#   gvMetab <- eval(parse( text = paste0("gvMetab2018Total", targetTypeNow)))
#
#   gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
#   gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]
#
#   See(gvMetab)
#
#   traitNames <- colnames(gvMetab)
#   lineNamesMetab <- rownames(gvMetab)
#   nTrait <- ncol(gvMetab)
#
#   lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))
#
#
#
#   gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
#   gastonData <- select.snps(gastonData, maf >= 0.025)
#   gastonData <- LD.thin(gastonData, threshold = thresLD)
#
#
#   K <- GRM(gastonData)
#   eigenK <- eigen(K)
#
#
#
#   dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
#                               "_", cultivationInfo, "_",
#                               targetInfo, "_results/")
#   dir.create(dirCultivationNow)
#
#   dirMarkerTypeNow <- paste0(dirCultivationNow,
#                              scriptID, "_",
#                              markerInterestID,
#                              "_results/")
#   dir.create(dirMarkerTypeNow)
#
#   dirTargetTypeNow <- paste0(dirMarkerTypeNow,
#                              scriptID, "_",
#                              targetTypeNow,
#                              "_results/")
#   dir.create(dirTargetTypeNow)
#
#   dirLDNow <- paste0(dirTargetTypeNow,
#                      scriptID, "_",
#                      "LD=1",
#                      "_results/")
#   dir.create(dirLDNow)
#
#
#   for (traitNo in 1:nTrait) {
#     traitNow <- traitNames[traitNo]
#     print(paste0(targetTypeNow, "_", traitNow))
#
#     gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]
#
#
#     isNAPheno <- is.na(gastonData@ped$pheno)
#
#     if (sum(isNAPheno) > 0) {
#       gastonDataNow <- gastonData[-which(isNAPheno), ]
#       gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
#       gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
#       # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)
#
#       KNow <- GRM(gastonDataNow)
#       eigenKNow <- eigen(KNow)
#     } else {
#       gastonDataNow <- gastonData
#
#       KNow <- K
#       eigenKNow <- eigenK
#     }
#     st <- Sys.time()
#     gastonRes <- association.test(x = gastonDataNow, method = "lmm",
#                                   response = "quantitative",
#                                   test = "wald", eigenK = eigenKNow,
#                                   p = 2)
#     end <- Sys.time()
#     print(end - st)
#
#     if (targetTypeNow == "0"){
#       betaMat0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
#     } else {
#       betaMat2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
#       gastonRes$id[which(gastonRes$id %in% markerNames)]
#       # table(gastonRes$id %in% "Chr17_16065902")
#       # table(gastonRes$id %in% markerNames)
#       # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]
#
#     }
#
#
#     gastonResOrd <- gastonRes[order(gastonRes$p), ]
#
#     pAdj <- p.adjust(gastonRes$p, method = "BH")
#
#     colGaston <- rep("black", nrow(gastonRes))
#     colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
#     colGaston[pAdj < 0.05] <- "blue"
#
#
#     dirSave0 <- paste0(dirLDNow, scriptID,
#                        "_", traitNow, "/")
#     dir.create(dirSave0)
#
#     dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
#     dir.create(dirSave)
#
#     fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
#                            thresLD, "_gaston_wald")
#
#     fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
#     write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)
#
#     fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
#     save(gastonRes, file = fileNameSaveRData)
#
#     # fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
#     # fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")
#     #
#     # png(fileNameSaveManhattan, width = 1200, height = 900)
#     # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
#     # dev.off()
#     #
#     # png(fileNameSaveQq, width = 900, height = 900)
#     # RAINBOWR::qq(- log10(gastonRes$p))
#     # dev.off()
#     #
#     #
#     # if (thresLD <= 0.4) {
#     #   gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
#     #                                         BP = gastonRes$pos,
#     #                                         P = gastonRes$p,
#     #                                         SNP = gastonRes$id,
#     #                                         BLOCK = gastonDataNow@snps$block)
#     #   sigSNPs <- gastonRes$id[pAdj < 0.05]
#     #
#     #   fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
#     #   fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")
#     #
#     #   plotlyManhattan <- manhattanly(gastonResForManhattanly,
#     #                                  snp = "SNP", gene = "BLOCK",
#     #                                  highlight = sigSNPs)
#     #   htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
#     #                           file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
#     #                                            basename(fileNameSavePlotlyManhattan)))
#     #
#     #   plotlyQq <- qqPlotly(data = gastonResForManhattanly,
#     #                        highlightSNPNames = sigSNPs)
#     #   htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
#     #                           file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
#     #                                            basename(fileNameSavePlotlyQq)))
#     # }
#     print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
#   }
# }
#
#
# fileName <- paste0(dirMidSTAMGWAS, scriptID, "_Chr06_18760995_2_coefficient_LD1.csv")
# write.csv(x = betaMat2LD1, file = fileName)






##### 3.3. Read genotypic values of Metabolomic data in 2018 of (0,0),(0,2),(2,0),(2,2) for Chr06_18760995 and Chr06_47428824 into R #####
#### 3.2.1 For each maker ####
cultivationInfo <- "2018"
targetInfo <- "Metabolome"

targetType <- c("0_0", "0_2", "2_0", "2_2")
# targetType <- c("0_2", "2_0", "2_2")

markerInterestID1 <- "Chr06_18760995"
markerInterestID2 <- "Chr06_47428824"



gvMetab2018Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2018.csv", row.names = 1)
See(gvMetab2018Total)
# rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"


gvMetab2018Total <- scale(gvMetab2018Total, center = TRUE, scale = TRUE)
See(gvMetab2018Total)
gvMetab2018Total <- as.data.frame(gvMetab2018Total)

lineNames <- rownames(gvMetab2018Total)
markerInterest1 <- gastonData0Matrix[, markerInterestID1]
markerInterest2 <- gastonData0Matrix[, markerInterestID2]


markerInterest1DF <- as.data.frame(markerInterest1)
markerInterest2DF <- as.data.frame(markerInterest2)
markerInterest1And2DF <- cbind(markerInterest1DF, markerInterest2DF)
See(markerInterest1And2DF)

# markerVarietyNames <- names(markerInterest)

table(markerInterest1And2DF$markerInterest1)
table(markerInterest1And2DF$markerInterest2)
marker0_0VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 0]
marker0_2VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 0 & markerInterest1And2DF$markerInterest2 == 2]
marker2_0VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 0]
marker2_2VarietyNames <- rownames(markerInterest1And2DF)[markerInterest1And2DF$markerInterest1 == 2 & markerInterest1And2DF$markerInterest2 == 2]

See(marker0_0VarietyNames)
See(marker0_2VarietyNames)
See(marker2_0VarietyNames)
See(marker2_2VarietyNames)


lineNames0_0 <- lineNames[(lineNames %in% marker0_0VarietyNames)]
commonNames0_0 <- marker0_0VarietyNames[(marker0_0VarietyNames %in% lineNames0_0)]
lineNames0_2 <- lineNames[(lineNames %in% marker0_2VarietyNames)]
commonNames0_2 <- marker0_2VarietyNames[(marker0_2VarietyNames %in% lineNames0_2)]
lineNames2_0 <- lineNames[(lineNames %in% marker2_0VarietyNames)]
commonNames2_0 <- marker2_0VarietyNames[(marker2_0VarietyNames %in% lineNames2_0)]
lineNames2_2 <- lineNames[(lineNames %in% marker2_2VarietyNames)]
commonNames2_2 <- marker2_2VarietyNames[(marker2_2VarietyNames %in% lineNames2_2)]

gvMetab2018Total0_0 <- gvMetab2018Total[commonNames0_0, ]
gvMetab2018Total0_2 <- gvMetab2018Total[commonNames0_2, ]
gvMetab2018Total2_0 <- gvMetab2018Total[commonNames2_0, ]
gvMetab2018Total2_2 <- gvMetab2018Total[commonNames2_2, ]
# gvMetab2018Total0 <- na.omit(gvMetab2018Total0)
See(gvMetab2018Total0_0, rown = 20)
See(gvMetab2018Total0_2)
See(gvMetab2018Total2_0)
See(gvMetab2018Total2_2)
# table(is.na(gvMetab2018Total2_0))
# table(is.na(gvMetab2018Total2_2))
# table(is.na(gvMetab2018Total0_0))
table(is.na(gvMetab2018Total2_0$X00006))



### LD1
markerNames <- c("Chr10_42562665", "Chr17_15930807")
betaMat0_0LD1 <- matrix(NA, nrow = length(gvMetab2018Total0_0), ncol = length(markerNames))
betaMat0_2LD1 <- matrix(NA, nrow = length(gvMetab2018Total0_2), ncol = length(markerNames))
betaMat2_0LD1 <- matrix(NA, nrow = length(gvMetab2018Total2_0), ncol = length(markerNames))
betaMat2_2LD1 <- matrix(NA, nrow = length(gvMetab2018Total2_2), ncol = length(markerNames))
rownames(betaMat0_0LD1) <- colnames(gvMetab2018Total0_0)
rownames(betaMat0_2LD1) <- colnames(gvMetab2018Total0_2)
rownames(betaMat2_0LD1) <- colnames(gvMetab2018Total2_0)
rownames(betaMat2_2LD1) <- colnames(gvMetab2018Total2_2)
colnames(betaMat0_0LD1) <- markerNames
colnames(betaMat0_2LD1) <- markerNames
colnames(betaMat2_0LD1) <- markerNames
colnames(betaMat2_2LD1) <- markerNames
See(betaMat0_0LD1)
See(betaMat0_2LD1)
See(betaMat2_0LD1)
See(betaMat2_2LD1)

# targetTypeNo <- 1
# traitNo <- 1

# targetType <- c("2_0", "2_2")
# targetType <- c("2_2")
thresLD <- 2
for (targetTypeNo in 1:length(targetType)) {

  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- eval(parse( text = paste0("gvMetab2018Total", targetTypeNow)))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]

  See(gvMetab)

  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))



  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf >= 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)


  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)

  dirMarkerTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             markerInterestID1,
                             "_",
                             markerInterestID2,
                             "_results/")
  dir.create(dirMarkerTypeNow)

  dirTargetTypeNow <- paste0(dirMarkerTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  dirLDNow <- paste0(dirTargetTypeNow,
                     scriptID, "_",
                     "LD=1",
                     "_results/")
  dir.create(dirLDNow)


  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf >= 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)
      # gastonDataNow <- LD.thin(gastonDataNow, threshold = 1)

      KNow <- GRM(gastonDataNow)
      eigenKNow <- eigen(KNow)
    } else {
      gastonDataNow <- gastonData

      KNow <- K
      eigenKNow <- eigenK
    }
    st <- Sys.time()
    gastonRes <- association.test(x = gastonDataNow, method = "lmm",
                                  response = "quantitative",
                                  test = "wald", eigenK = eigenKNow,
                                  p = 2)
    end <- Sys.time()
    print(end - st)

    if (targetTypeNow == "0_0"){
      betaMat0_0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "0_2") {
      betaMat0_2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "2_0") {
      betaMat2_0LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
    } else if (targetTypeNow == "2_2") {
      betaMat2_2LD1[traitNow, markerNames] <- gastonRes$beta[gastonRes$id %in% markerNames]
      # gastonRes$id[which(gastonRes$id %in% markerNames)]
      # table(gastonRes$id %in% "Chr17_16065902")
      # table(gastonRes$id %in% markerNames)
      # gastonRes$id[grepl(pattern = "Chr06", x = gastonRes$id)]

    }


    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirLDNow, scriptID,
                       "_", traitNow, "/")
    dir.create(dirSave0)

    dirSave <- paste0(dirSave0, scriptID, "_gaston_results/")
    dir.create(dirSave)

    fileNameSave <- paste0(dirSave, scriptID, "_", traitNow, "_thresLD=",
                           thresLD, "_gaston_wald")

    fileNameSaveOrderedRes <- paste0(fileNameSave, "_ordered_results.csv")
    write.csv(gastonResOrd, file = fileNameSaveOrderedRes, quote = FALSE)

    fileNameSaveRData <- paste0(fileNameSave, "_raw_results.RData")
    save(gastonRes, file = fileNameSaveRData)

    fileNameSaveManhattan <- paste0(fileNameSave, "_manhattan_plot.png")
    fileNameSaveQq <- paste0(fileNameSave, "_qq_plot.png")

    png(fileNameSaveManhattan, width = 1200, height = 900)
    gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    dev.off()

    png(fileNameSaveQq, width = 900, height = 900)
    RAINBOWR::qq(- log10(gastonRes$p))
    dev.off()


    if (thresLD <= 0.4) {
      gastonResForManhattanly <- data.frame(CHR = gastonRes$chr,
                                            BP = gastonRes$pos,
                                            P = gastonRes$p,
                                            SNP = gastonRes$id,
                                            BLOCK = gastonDataNow@snps$block)
      sigSNPs <- gastonRes$id[pAdj < 0.05]

      fileNameSavePlotlyManhattan <- paste0(fileNameSave, "_manhattan_plotly.html")
      fileNameSavePlotlyQq <- paste0(fileNameSave, "_qq_plotly.html")

      plotlyManhattan <- manhattanly(gastonResForManhattanly,
                                     snp = "SNP", gene = "BLOCK",
                                     highlight = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyManhattan),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyManhattan)),
                                               basename(fileNameSavePlotlyManhattan)))

      plotlyQq <- qqPlotly(data = gastonResForManhattanly,
                           highlightSNPNames = sigSNPs)
      htmlwidgets::saveWidget(widget = partial_bundle(plotlyQq),
                              file = file.path(normalizePath(dirname(fileNameSavePlotlyQq)),
                                               basename(fileNameSavePlotlyQq)))
    }
    print(paste0("There are ", sum(pAdj < 0.05), " peaks!"))
  }
}




### comparison of coefficients
dirTwoMarkers2018 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2018_Metabolome_results/2.32_Coefficient_LD1_scaled_background_Chr06_18760995_Chr06_47428824/")
dir.create(dirTwoMarkers2018)


metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]

table(abs(betaMat0 - betaMat2) > 1)

colVec <- rep("blue1", nrow(betaMat0))
names(colVec) <- rownames(betaMat0)
colVec[metabNamesFlavonoid] <- "orange1"



dir.create(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/"))
dir.create(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/"))

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(0,2).pdf"))
plot(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat0_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(2,0).pdf"))
plot(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_(0,0)_(2,2).pdf"))
plot(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_(0,2)_(2,0).pdf"))
plot(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1],
     xlim = range(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 1], betaMat2_0LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_(0,2)_(2,2).pdf"))
plot(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_(2,0)_(2,2).pdf"))
plot(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1],
     xlim = range(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, 1], betaMat2_2LD1[, 1], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 2, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()



pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_(0,0)_(0,2).pdf"))
plot(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat0_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 0, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_(0,0)_(2,0).pdf"))
plot(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_(0,0)_(2,2).pdf"))
plot(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_(0,2)_(2,0).pdf"))
plot(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2],
     xlim = range(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 2], betaMat2_0LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 0",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_(0,2)_(2,2).pdf"))
plot(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat0_2LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0, Chr06_47428824 = 2",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_(2,0)_(2,2).pdf"))
plot(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2],
     xlim = range(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     ylim = range(betaMat2_0LD1[, 2], betaMat2_2LD1[, 2], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 2, Chr06_47428824 = 0",
     ylab = "Chr06_18760995 = 2, Chr06_47428824 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()



### pairs plot
Chr10_42562665DF <- data.frame("0_0" = betaMat0_0LD1[, 1], "0_2" = betaMat0_2LD1[, 1], "2_0" = betaMat2_0LD1[, 1], "2_2" = betaMat2_2LD1[, 1])
colnames(Chr10_42562665DF) <- c("0_0", "0_2", "2_0", "2_2")
Chr17_15930807DF <- data.frame("0_0" = betaMat0_0LD1[, 2], "0_2" = betaMat0_2LD1[, 2], "2_0" = betaMat2_0LD1[, 2], "2_2" = betaMat2_2LD1[, 2])
colnames(Chr17_15930807DF) <- c("0_0", "0_2", "2_0", "2_2")


# upper <- function(x, y, ...){
#   oldpar <- par(usr = c(0, 1, 0, 1))
#   v <- abs(cor(x, y))
#   sz <- 1 + v * 2
#   text(0.5, 0.5, sprintf("%.3f", v), cex = sz)
#   par(oldpar)
# }

lower <- function(x, y, ...){
  points(x, y, pch = 21, bg = c("blue1", "orange1")[as.numeric(as.factor(colVec))])
  if( 1==1 ){
    abline(0, 1, col = 2, lty = 2, lwd = 1.5)
  }
}

# pairs(Chr10_42562665DF, upper.panel=upper, lower.panel=lower)
pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr10_42562665/", scriptID, "_pairs_plot.pdf"))
pairs(Chr10_42562665DF, lower.panel=lower, oma = c(3,3,3,15))
par(xpd = TRUE)
legend(x = "bottomright", legend = c("Non-flavonoid", "Flavonoid"), col = c("blue1", "orange1"), pch = 19, cex = 1)
dev.off()

pdf(paste0(dirTwoMarkers2018, scriptID, "_Chr17_15930807/", scriptID, "_pairs_plot.pdf"))
pairs(Chr17_15930807DF, lower.panel=lower, oma = c(3,3,3,15))
par(xpd = TRUE)
legend(x = "bottomright", legend = c("Non-flavonoid", "Flavonoid"), col = c("blue1", "orange1"), pch = 19)
dev.off()



# Chr10_42562665DF <- data.frame("0_0" = betaMat0_0LD1[, 1], "0_2" = betaMat0_2LD1[, 1], "2_0" = betaMat2_0LD1[, 1], "2_2" = betaMat2_2LD1[, 1])
# colnames(Chr10_42562665DF) <- c("0_0", "0_2", "2_0", "2_2")
# Chr17_15930807DF <- data.frame("0_0" = betaMat0_0LD1[, 2], "0_2" = betaMat0_2LD1[, 2], "2_0" = betaMat2_0LD1[, 2], "2_2" = betaMat2_2LD1[, 2])
# colnames(Chr17_15930807DF) <- c("0_0", "0_2", "2_0", "2_2")
#
# See(Chr10_42562665DF)
# See(Chr17_15930807DF)


# ggpairs(Chr10_42562665DF, aes_string(colour="colVec", alpha=0.5)) + geom_abline(slope = 1, intercept = 0)
# ggpairs(Chr17_15930807DF, aes_string(colour="colVec", alpha=0.5))
#
#
# ggpairs(Chr10_42562665DF, mapping = aes(color = colVec), lower = list(continuous = "smooth"))
# ggpairs(Chr17_15930807DF, mapping = aes(color = colVec), lower = list(continuous = "smooth"))




fileName0_0 <- paste0(dirTwoMarkers2018, scriptID, "_(0,0).csv")
fileName0_2 <- paste0(dirTwoMarkers2018, scriptID, "_(0,2).csv")
fileName2_0 <- paste0(dirTwoMarkers2018, scriptID, "_(2,0).csv")
fileName2_2 <- paste0(dirTwoMarkers2018, scriptID, "_(2,2).csv")


write.csv(x = betaMat0_0LD1, file = fileName0_0)
write.csv(x = betaMat0_2LD1, file = fileName0_2)
write.csv(x = betaMat2_0LD1, file = fileName2_0)
write.csv(x = betaMat2_2LD1, file = fileName2_2)



### extract metabolites
betaMat0_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2018_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47428824/2.32_(0,0).csv", row.names = 1)
betaMat0_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2018_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47428824/2.32_(0,2).csv", row.names = 1)
betaMat2_0LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2018_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47428824/2.32_(2,0).csv", row.names = 1)
betaMat2_2LD1 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2018_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47428824/2.32_(2,2).csv", row.names = 1)


#### Chr10_42562665
### (0,0),(0,2)
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 2, ])
# rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 0.8 & betaMat0_2LD1[, 1] < 0.8 , ])

### (2,0),(2,2)
rownames(betaMat2_0LD1[betaMat2_0LD1[, 1] > 2, ])
rownames(betaMat2_0LD1[betaMat2_0LD1[, 1] < -2, ])

### (0,0),(2,0)
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] > 2, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 1] < -0.6, ])

### (0,2),(2,2)
rownames(betaMat0_2LD1[betaMat0_2LD1[, 1] > 0.45, ])
rownames(betaMat2_2LD1[betaMat2_2LD1[, 1] > 1, ])
rownames(betaMat0_2LD1[betaMat0_2LD1[, 1] < -0.6, ])
rownames(betaMat2_2LD1[betaMat2_2LD1[, 1] < -1, ])


#### Chr17_15930807
### (0,0),(0,2)


### (2,0),(2,2)


### (0,0),(2,0)
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] > 0.4 & betaMat2_0LD1[, 2] > -0.2, ])
rownames(betaMat0_0LD1[betaMat0_0LD1[, 2] < -0.4, ])

### (0,2),(2,2)
rownames(betaMat0_2LD1[betaMat0_2LD1[, 2] > 0.24, ])
rownames(betaMat2_2LD1[betaMat2_2LD1[, 2] > 0.4, ])








