##########################################################################################
######  Title: 2.6_Soybean_STAM_GWAS_for_metabolomic_data                           ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/13 (Created), 2024/08/14, 2026/01/13 (Last Updated)                       ######
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

scriptID <- "2.6"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMGWAS <- paste0(dirMidSTAMBase, scriptID,
                         "_GWAS/")
dir.create(dirMidSTAMGWAS)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)




thresLD <- 0.95


##### 1.3. Import packages #####
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)
require(plotly)
require(manhattanly)

source("scripts/1.1_Soybean_STAM_qqPlotly.R")
source("scripts/manhattanTaisei.R")

##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### 2. Perform gaston for Metabolomic BLUP, from raw data without outlier ######
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
See(gastonDataSmall)




##### 2.2. 2017 #####
#### 2.2.1. Read genotypic values of Metabolomic data in 2017 into R ####
cultivationInfo <- "2017"
targetInfo <- "Metabolome"

targetType <- c("Total", "Control", "Drought",
                "CPlusD", "CMinusD")
# targetType <- c("CMinusD")
# targetTypeNo <- 1
# traitNo <- 57

for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.2_BSH/2.2_lmer_genotypic_values_",
                             targetTypeNow, "_", cultivationInfo,
                             ".csv"),
                      row.names = 1)

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]


  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))


  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf > 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)

  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)
  dirTargetTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf > 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)

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

    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirTargetTypeNow, scriptID,
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

    png(fileNameSaveManhattan, width = 1500, height = 900)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    manhattanTaisei(gastonRes, pch = 20, col = colGaston,
                      cex = 2.5, cex.axis = 2, cex.axis.x = 2, cex.lab = 1)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
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




##### 2.3. 2018 #####
#### 2.3.1. Read genotypic values of Metabolomic data in 2018 into R ####
cultivationInfo <- "2018"
targetInfo <- "Metabolome"

targetType <- c("Total", "Control", "Drought",
                "CPlusD", "CMinusD")


for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.2_BSH/2.2_lmer_genotypic_values_",
                             targetTypeNow, "_", cultivationInfo,
                             ".csv"),
                      row.names = 1)

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]


  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))


  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf > 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)

  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)
  dirTargetTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf > 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)

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

    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirTargetTypeNow, scriptID,
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







###### 3. Perform gaston for Metabolomic BLUP, from Box-Cox transformed data without outlier ######
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
See(gastonDataSmall)




##### 3.2. 2017 #####
#### 3.2.1. Read genotypic values of Metabolomic data in 2017 into R ####
cultivationInfo <- "2017"
targetInfo <- "Metabolome"
dataType <- "BoxCox"

targetType <- c("Total", "Control", "Drought",
                "CPlusD", "CMinusD")

traitNamesOrder <- c("X00006", "X00255", "X00698", "X00735", "X00799", "X00867", "X200076", "X200079", "X500113")

# gvMetab2017BoxCox <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017_BoxCox.csv", row.names = 1)
# See(gvMetab2017BoxCox)
#
# traitNos <- match(traitNamesOrder, colnames(gvMetab2017BoxCox))

# targetType <- c("CMinusD")
# targetTypeNo <- 1
# traitNo <- 57
# gvMetab2017BoxCoxD <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Drought_2017_BoxCox.csv", row.names = 1)
# colnames(gvMetab2017BoxCoxD)
# gvMetab2017BoxCoxD[, 73]

for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.2_BSH/2.2_lmer_genotypic_values_",
                             targetTypeNow, "_", cultivationInfo, "_", dataType,
                             ".csv"),
                      row.names = 1)

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]


  traitNames <- colnames(gvMetab)
  traitNos <- match(traitNamesOrder, traitNames)

  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))


  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf > 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)

  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)
  dirDataTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             dataType,
                             "_results/")
  dir.create(dirDataTypeNow)
  dirTargetTypeNow <- paste0(dirDataTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  # for (traitNo in 1:nTrait) {
  for(traitNo in traitNos){
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf > 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)

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

    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirTargetTypeNow, scriptID,
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

    png(fileNameSaveManhattan, width = 1500, height = 900)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    manhattanTaisei(gastonRes, pch = 20, col = colGaston,
                    cex = 2.5, cex.axis = 2, cex.axis.x = 2, cex.lab = 1)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
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




##### 3.3. 2018 #####
#### 3.3.1. Read genotypic values of Metabolomic data in 2018 into R ####
cultivationInfo <- "2018"
targetInfo <- "Metabolome"
dataType <- "BoxCox"

targetType <- c("Total", "Control", "Drought",
                "CPlusD", "CMinusD")

traitNamesOrder <- c("X00006", "X00255", "X00698", "X00735", "X00799", "X00867", "X200076", "X200079", "X500113")

# gvMetab2017BoxCox <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017_BoxCox.csv", row.names = 1)
# See(gvMetab2017BoxCox)
#
# traitNos <- match(traitNamesOrder, colnames(gvMetab2017BoxCox))

# targetType <- c("CMinusD")
# targetTypeNo <- 1
# traitNo <- 57
# gvMetab2017BoxCoxD <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Drought_2017_BoxCox.csv", row.names = 1)
# colnames(gvMetab2017BoxCoxD)
# gvMetab2017BoxCoxD[, 73]

for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.2_BSH/2.2_lmer_genotypic_values_",
                             targetTypeNow, "_", cultivationInfo, "_", dataType,
                             ".csv"),
                      row.names = 1)

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]


  traitNames <- colnames(gvMetab)
  traitNos <- match(traitNamesOrder, traitNames)

  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))


  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf > 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)

  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)
  dirDataTypeNow <- paste0(dirCultivationNow,
                           scriptID, "_",
                           dataType,
                           "_results/")
  dir.create(dirDataTypeNow)
  dirTargetTypeNow <- paste0(dirDataTypeNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  # for (traitNo in 1:nTrait) {
  for(traitNo in traitNos){
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf > 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)

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

    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirTargetTypeNow, scriptID,
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

    png(fileNameSaveManhattan, width = 1500, height = 900)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
    manhattanTaisei(gastonRes, pch = 20, col = colGaston,
                    cex = 2.5, cex.axis = 2, cex.axis.x = 2, cex.lab = 1)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston, cex = 2.5)
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










###### 10. Root Metabolomic data in 2020 ######
cultivationInfo <- "2020"
targetInfo <- "Root_Metabolome"

targetType <- c("Control", "Drought")

gvMetab <- read.csv(paste0("raw_data/phenotype/2020_Tottori_Main_RootMetabolome.csv"))
See(gvMetab, coln = 12)

gvMetabControl <- gvMetab[gvMetab$block == "W1", ]
gvMetabDrought <- gvMetab[gvMetab$block == "W4", ]
rownames(gvMetabControl) <- gvMetabControl$variety
rownames(gvMetabDrought) <- gvMetabDrought$variety
See(gvMetabControl)
See(gvMetabDrought)

gvMetabControl <- gvMetabControl[, 10:ncol(gvMetab)]
gvMetabDrought <- gvMetabDrought[, 10:ncol(gvMetab)]

targetTypeNo <- 1
traitNo <- 1
for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]

  assign(paste0("gvMetab"), eval(parse(text = paste0("gvMetab", targetType[targetTypeNo]))))

  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(is.na(x))))]
  gvMetab <- gvMetab[, which(!apply(gvMetab, 2, function(x) all(x == 0)))]


  traitNames <- colnames(gvMetab)
  lineNamesMetab <- rownames(gvMetab)
  nTrait <- ncol(gvMetab)

  lineNamesMatch <- Reduce(intersect, list(lineNamesGeno, lineNamesMetab))


  gastonData <- gastonDataSmall[lineNamesGeno %in% lineNamesMatch, ]
  gastonData <- select.snps(gastonData, maf > 0.025)
  gastonData <- LD.thin(gastonData, threshold = thresLD)

  K <- GRM(gastonData)
  eigenK <- eigen(K)



  dirCultivationNow <- paste0(dirMidSTAMGWAS, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results/")
  dir.create(dirCultivationNow)
  dirTargetTypeNow <- paste0(dirCultivationNow,
                             scriptID, "_",
                             targetTypeNow,
                             "_results/")
  dir.create(dirTargetTypeNow)

  for (traitNo in 1:nTrait) {
    traitNow <- traitNames[traitNo]
    print(paste0(targetTypeNow, "_", traitNow))

    gastonData@ped$pheno <- gvMetab[lineNamesMatch, traitNow]


    isNAPheno <- is.na(gastonData@ped$pheno)

    if (sum(isNAPheno) > 0) {
      gastonDataNow <- gastonData[-which(isNAPheno), ]
      gastonDataNow <- select.snps(gastonDataNow, maf > 0.025)
      gastonDataNow <- LD.thin(gastonDataNow, threshold = thresLD)

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

    gastonResOrd <- gastonRes[order(gastonRes$p), ]

    pAdj <- p.adjust(gastonRes$p, method = "BH")

    colGaston <- rep("black", nrow(gastonRes))
    colGaston[gastonRes$chr %% 2 == 0] <- "gray50"
    colGaston[pAdj < 0.05] <- "blue"


    dirSave0 <- paste0(dirTargetTypeNow, scriptID,
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



