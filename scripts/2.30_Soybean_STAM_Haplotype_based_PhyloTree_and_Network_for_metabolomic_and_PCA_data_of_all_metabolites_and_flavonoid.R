##########################################################################################
######  Title: 2.35_Soybean_STAM_Haplotype_based_PhyloTree_and_Network_for_metabolomic_and_PCA_data_of_all_metabolites_and_flavonoid  ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/04/06 (Created), 2024/04/25 (Last Updated)                       ######
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

scriptID <- "2.35"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMHaplotypeBasedPhyloTreeAndNetwork <- paste0(dirMidSTAMBase, scriptID, "_Haplotype_based_PhyloTree_and_Network/")
dir.create(dirMidSTAMHaplotypeBasedPhyloTreeAndNetwork)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)


# cultivationInfo <- "2017"
# targetInfo <- "Metabolome"

# targetType <- c("Total", "Control", "Drought", "CPlusD", "CMinusD")


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
source("scripts/10.0_Soybean_STAM_myPlotHaploNetwork.R")

##### 1.4. Project options #####
options(stringAsFactors = FALSE)




###### 2. Perform Phylotree and Network analysis of PCs for all Metabolomic data in 2017 ######
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
genoMat <- as.matrix(gastonDataSmall)


##### 2.2. Read lmer PCA score and Group infomation into R #####
lmerGvMetabPCA <- read.csv("midstream/2.5_BSH_for_PCA/2.5_lmer_genotypic_values_for_PC_Score_in_2017.csv", row.names = 1)
lmerGvMetabPCA <- as.matrix(lmerGvMetabPCA)
rownames(lmerGvMetabPCA)[rownames(lmerGvMetabPCA) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
See(lmerGvMetabPCA)

groupInfo <- read.csv("midstream/2.30_Haplotype_based_PhyloTree_and_Network/0.3_group_information.csv", row.names = 1)
groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"


##### 2.3. Perform Phylotree and Network analysis of PC2 #####
detectedSnp <- "Chr06_18760995"
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc2 <- data.frame(line = rownames(phenoMat),
                       PC2 = phenoMat[, 2])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc2,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes$haplotypeInfo$haploCluster
estPhyloRes



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc2,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype




##### 2.4. Check marker values for each group #####
groupInfo <- read.csv("midstream/2.22_Haplotype_based_PhyloTree_and_Network/0.3_group_information.csv", row.names = 1)
groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"
class(groupInfo)
head(groupInfo)
tail(groupInfo)
lineNamesGroupInfo <- groupInfo$lineNames

gastonData0Matrix <- as.matrix(gastonData0)

markerInterestID <- "Chr06_18760995"
markerInterest <- gastonData0Matrix[, markerInterestID]
class(markerInterest)
head(markerInterest)
tail(markerInterest)
markerLineNames <- names(markerInterest)

commonLineNames <- markerLineNames[markerLineNames %in% lineNamesGroupInfo]

rownames(groupInfo) <- groupInfo$lineNames
groupInfo <- groupInfo[commonLineNames, ]
class(groupInfo)
class(markerInterest)
markerInterest <- as.data.frame(markerInterest)
lineNames <- rownames(markerInterest)
markerInterest <- transform(markerInterest, lineNames = lineNames)

markerGroup <- inner_join(markerInterest, groupInfo, by = "lineNames")
head(markerGroup)
tail(markerGroup)

barplot(markerGroup$markerInterest, markerGroup$group)

ggplot(markerGroup, aes(x = factor(markerInterest), y = group , fill = group)) +
  geom_col(position = "dodge", colour = "black") +
  #scale_x_discrete(limit = c("M", "F")) +
  scale_fill_discrete(limit = c("WTM", "PdM", "WTF", "PdM")) +
  scale_fill_manual(values = c("#2980B9", "#154360", "#CB4335", "#641E16"))



markerGroup
table(markerGroup$group)
sum(table(markerGroup$group))

marker0Group <- markerGroup[markerGroup$markerInterest == 0, ]
marker2Group <- markerGroup[markerGroup$markerInterest == 2, ]

table(marker0Group$group)
table(marker2Group$group)

table(marker0Group$group) / 198
table(marker2Group$group) / 198

table(marker0Group$group) / sum(table(marker0Group$group))
table(marker2Group$group) / sum(table(marker2Group$group))





###### 3. Perform Phylotree and Network analysis of PCs for Flavonoid Metabolomic data in 2017 ######
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
genoMat <- as.matrix(gastonDataSmall)


##### 3.2. Read lmer PCA score and Group infomation into R #####
lmerGvMetabFlavonoidPCA <- read.csv("midstream/2.18_BSH_for_PCA_with_or_without_flavonoid/2.18_lmer_genotypic_values_for_PC_Score_for_flavonoid_related_metab_in_2017.csv", row.names = 1)
lmerGvMetabFlavonoidPCA <- as.matrix(lmerGvMetabFlavonoidPCA)
rownames(lmerGvMetabFlavonoidPCA)[rownames(lmerGvMetabFlavonoidPCA) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
See(lmerGvMetabFlavonoidPCA)

groupInfo <- read.csv("midstream/2.22_Haplotype_based_PhyloTree_and_Network/0.3_group_information.csv", row.names = 1)
groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"


##### 3.3. Perform Phylotree and Network analysis #####
#### 3.3.1 PC1 ####
detectedSnp <- "Chr06_18760995"
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabFlavonoidPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc1 <- data.frame(line = rownames(phenoMat),
                       PC1 = phenoMat[, 1])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc1,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes$haplotypeInfo$haploCluster
estPhyloRes$pValChi2Test
estPhyloRes$clusterNosForHaplotype



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc1,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype



#### 2.3.2 PC2 ####
detectedSnp <- "Chr06_21382495"
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabFlavonoidPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc2 <- data.frame(line = rownames(phenoMat),
                       PC2 = phenoMat[, 2])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc2,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes$haplotypeInfo$haploCluster



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc2,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype



#### 2.3.3 PC3 ####
detectedSnp <- "Chr06_47606869"
# detectedSnp <- "Chr18_6592763"
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabFlavonoidPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc2 <- data.frame(line = rownames(phenoMat),
                       PC2 = phenoMat[, 2])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc2,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes$haplotypeInfo$haploCluster



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc1,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype



#### 2.3.4.1 PC4 ####
detectedSnp <- "Chr10_42562665"
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabFlavonoidPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc4 <- data.frame(line = rownames(phenoMat),
                       PC4 = phenoMat[, 4])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc4,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes$clusterNosForHaplotype
estPhyloRes$haplotypeInfo$haploCluster
estPhyloRes$subpopInfo
estPhyloRes$haplotypeInfo
estPhyloRes$pValChi2Test



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc4,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype



##### 2.3.4.2 Check marker values for each group #####
groupInfo <- read.csv("midstream/2.22_Haplotype_based_PhyloTree_and_Network/0.3_group_information.csv", row.names = 1)
groupInfo$lineNames[groupInfo$lineNames == "Houjaku-Kuwazu"] <- "HOUJAKU_KUWAZU"
class(groupInfo)
head(groupInfo)
tail(groupInfo)
lineNamesGroupInfo <- groupInfo$lineNames

gastonData0Matrix <- as.matrix(gastonData0)

markerInterestID <- "Chr10_42562665"
markerInterest <- gastonData0Matrix[, markerInterestID]
class(markerInterest)
head(markerInterest)
tail(markerInterest)
markerLineNames <- names(markerInterest)

commonLineNames <- markerLineNames[markerLineNames %in% lineNamesGroupInfo]

rownames(groupInfo) <- groupInfo$lineNames
groupInfo <- groupInfo[commonLineNames, ]
class(groupInfo)
class(markerInterest)
markerInterest <- as.data.frame(markerInterest)
lineNames <- rownames(markerInterest)
markerInterest <- transform(markerInterest, lineNames = lineNames)

markerGroup <- inner_join(markerInterest, groupInfo, by = "lineNames")
head(markerGroup)
tail(markerGroup)

# barplot(markerGroup$markerInterest, markerGroup$group)

# ggplot(markerGroup, aes(x = factor(markerInterest), y = group , fill = group)) +
#   geom_col(position = "dodge", colour = "black") +
#   #scale_x_discrete(limit = c("M", "F")) +
#   scale_fill_discrete(limit = c("WTM", "PdM", "WTF", "PdM")) +
#   scale_fill_manual(values = c("#2980B9", "#154360", "#CB4335", "#641E16"))



markerGroup
table(markerGroup$group)
sum(table(markerGroup$group))

marker0Group <- markerGroup[markerGroup$markerInterest == 0, ]
marker2Group <- markerGroup[markerGroup$markerInterest == 2, ]

table(marker0Group$group)
table(marker2Group$group)

table(marker0Group$group) / 198
table(marker2Group$group) / 198

table(marker0Group$group) / sum(table(marker0Group$group))
table(marker2Group$group) / sum(table(marker2Group$group))




#### 2.3.5 PC5 ####
detectedSnp <- ""
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabFlavonoidPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc2 <- data.frame(line = rownames(phenoMat),
                       PC2 = phenoMat[, 2])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc2,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc1,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype



#### 2.3.6 PC6 ####
detectedSnp <- "Chr06_47426527"
detectedBlockName <- blockNames[mrkNames == detectedSnp]
detectedMarkers0 <- names(blockNames0[blockNames0 == detectedBlockName])
detectedMarkers <- detectedMarkers0[detectedMarkers0 %in% mrkNames]

rownames(genoMat)

modifyRes <- modify.data(pheno.mat = lmerGvMetabFlavonoidPCA,
                         geno.mat = genoMat,
                         return.ZETA = TRUE,
                         return.GWAS.format = TRUE)
phenoMat <- modifyRes$pheno.modi
genoModi <- modifyRes$geno.modi
ZETA <- modifyRes$ZETA

detectedMakrersLd <- colnames(genoModi)[colnames(genoModi) %in% detectedMarkers]
genoMatInBlock <- genoModi[, detectedMakrersLd]

phenoPc6 <- data.frame(line = rownames(phenoMat),
                       PC6 = phenoMat[, 6])
# phenoPc2 <- data.frame(line = rownames(phenoModi),
#                        PC2 = phenoModi[, 2])

grpInfoVec <- groupInfo$group
names(grpInfoVec) <- groupInfo$lineNames
grpNow <- grpInfoVec[rownames(phenoMat)]


estPhyloRes <- estPhylo(blockInterest = genoMatInBlock,
                        blockName = detectedBlockName,
                        pheno = phenoPc6,
                        subpopInfo = factor(grpNow),
                        chi2Test = TRUE,
                        kernelTypes = "gaussian",
                        rangeHStart = 10^c(-1:4),
                        plotTree = TRUE,
                        ggPlotTree = FALSE,
                        ZETA = ZETA)
# plotPhyloTree(estPhyloRes = estPhyloRes,
#               plotTree = FALSE,
#               ggPlotTree = TRUE)
class(estPhyloRes)
estPhyloRes$clusterNosForHaplotype
estPhyloRes$haplotypeInfo$haploCluster
estPhyloRes$subpopInfo
estPhyloRes$haplotypeInfo
estPhyloRes$pValChi2Test



estNetworkRes <- estNetwork(blockInterest = genoMatInBlock,
                            blockName = detectedBlockName,
                            pheno = phenoPc6,
                            subpopInfo = factor(grpNow),
                            nGrp = 3,
                            chi2Test = TRUE,
                            ZETA = ZETA,
                            kernelTypes = "diffusion",
                            plotNetwork = FALSE,
                            ggPlotNetwork = FALSE)
myPlotHaploNetwork(estNetworkRes = estNetworkRes,
                   colHaploBase = c(3, 5, 6),
                   plotNetwork = FALSE,
                   ggPlotNetwork = TRUE)
# estNetworkRes$clusterNosForHaplotype
