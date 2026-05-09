##########################################################################################
######  Title: 2.27_Soybean_STAM_GWAS_for_PCA_with_high_heritabiliity_flavonoid     ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/04 (Created), 2024/08/18 (Last Updated)                       ######
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

scriptID <- "2.27"



##### 1.2. Setting some parameters #####
dirMidSTAMBasePCAFlavonoid <- "midstream/"

dirMidSTAMGWASPCAFlavonoidMoreThan0.9Heritability <- paste0(dirMidSTAMBasePCAFlavonoid, scriptID, "_GWAS_for_PCA_with_high_heritability_flavonoid/")
dir.create(dirMidSTAMGWASPCAFlavonoidMoreThan0.9Heritability)
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



###### 3. Read genotypic values for PCA with flavonoid metab >0.9 heritability in 2017 into R and GWAS ######
##### 3.1. pcaMethods, nPC = 6 #####
# targetTypeNo <- 1
nPC <- 6
PCAMethod <- "pcaMethods"

cultivationInfo <- "2017"
targetInfo <- "Metabolome"

targetType <- "flavonoid_>0.9_heritability"
# targetType <- c("flavonoid_related", "flavonoid_non_related")


for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_2017_", PCAMethod, "_nPC=", nPC, ".csv"), row.names = 1)
  #gvMetab <- read.csv(paste0("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_",
  #                           targetTypeNow, "_of_PCA_", cultivationInfo,
  #                           ".csv"),
  #                    row.names = 1)

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



  dirCultivationNow <- paste0(dirMidSTAMGWASPCAFlavonoidMoreThan0.9Heritability, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results_", PCAMethod, "_nPC=", nPC, "/")
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

    png(fileNameSaveManhattan, width = 1400, height = 900)
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston)
    manhattanTaisei(gastonRes, pch = 20, col = colGaston,
                    cex = 2.5, cex.axis = 2, cex.axis.x = 2)
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



##### 3.2. pcaMethods, nPC = 20 #####
# targetTypeNo <- 1
nPC <- 20
PCAMethod <- "pcaMethods"

for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_2017_", PCAMethod, "_nPC=", nPC, ".csv"), row.names = 1)
  #gvMetab <- read.csv(paste0("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_",
  #                           targetTypeNow, "_of_PCA_", cultivationInfo,
  #                           ".csv"),
  #                    row.names = 1)

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



  dirCultivationNow <- paste0(dirMidSTAMGWASPCAFlavonoidMoreThan0.9Heritability, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results_", PCAMethod, "_nPC=", nPC, "/")
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
    colGaston[pAdj < 0.05] <- "green"


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
    gaston::manhattan(gastonRes, pch = 20, col = colGaston)
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




###### 4. Read genotypic values for PCA with flavonoid metab >0.9 heritability in 2018 into R and GWAS ######
##### 4.1. pcaMethods, nPC = 6 #####
# targetTypeNo <- 1
nPC <- 6
PCAMethod <- "pcaMethods"

cultivationInfo <- "2018"
targetInfo <- "Metabolome"

targetType <- "flavonoid_>0.9_heritability"
# targetType <- c("flavonoid_related", "flavonoid_non_related")


for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_",  cultivationInfo, "_", PCAMethod, "_nPC=", nPC, ".csv"), row.names = 1)
  #gvMetab <- read.csv(paste0("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_",
  #                           targetTypeNow, "_of_PCA_", cultivationInfo,
  #                           ".csv"),
  #                    row.names = 1)

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



  dirCultivationNow <- paste0(dirMidSTAMGWASPCAFlavonoidMoreThan0.9Heritability, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results_", PCAMethod, "_nPC=", nPC, "/")
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
    # gaston::manhattan(gastonRes, pch = 20, col = colGaston)
    manhattanTaisei(gastonRes, pch = 20, col = colGaston,
                    cex = 2.5, cex.axis = 2, cex.axis.x = 2)
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



##### 4.2. pcaMethods, nPC = 20 #####
# targetTypeNo <- 1
nPC <- 20
PCAMethod <- "pcaMethods"

cultivationInfo <- "2018"
targetInfo <- "Metabolome"

targetType <- "flavonoid_>0.9_heritability"
# targetType <- c("flavonoid_related", "flavonoid_non_related")


for (targetTypeNo in 1:length(targetType)) {
  targetTypeNow <- targetType[targetTypeNo]
  gvMetab <- read.csv(paste0("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_", cultivationInfo, "_", PCAMethod, "_nPC=", nPC, ".csv"), row.names = 1)
  #gvMetab <- read.csv(paste0("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_",
  #                           targetTypeNow, "_of_PCA_", cultivationInfo,
  #                           ".csv"),
  #                    row.names = 1)

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



  dirCultivationNow <- paste0(dirMidSTAMGWASPCAFlavonoidMoreThan0.9Heritability, scriptID,
                              "_", cultivationInfo, "_",
                              targetInfo, "_results_", PCAMethod, "_nPC=", nPC, "/")
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
    colGaston[pAdj < 0.05] <- "green"


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
    gaston::manhattan(gastonRes, pch = 20, col = colGaston)
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







