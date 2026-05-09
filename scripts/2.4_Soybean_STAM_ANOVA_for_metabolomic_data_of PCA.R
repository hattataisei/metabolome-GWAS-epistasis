##########################################################################################
######  Title: 2.4_Soybean_STAM_ANOVA_for_metab_PCA_data                          ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/07 (Created), 2024/08/14 (Last Updated)                       ######
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

scriptID <- "2.4"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
sigLevel <- 0.05

dirMidSTAMANOVAForPCA <- paste0(dirMidSTAMBase, scriptID,
                          "_ANOVA_for_PCA/")
dir.create(dirMidSTAMANOVAForPCA)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)




##### 1.4. Project options #####
options(stringAsFactors = FALSE)






###### 2. Perform ANOVA and adjust p-values for Metabolomic data in 2017 ######
##### 2.1. Read Metabolomic data in 2017 into R #####
metabPCA2017 <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_2017.csv")
See(metabPCA2017, rown = 6, coln = 12)



##### 2.2. Perform ANOVA for X00006 in Metabolomic data in 2017 #####
#lmRes <- lm(formula = X00006 ? variety + block + variety * block, data = metabPCA2017)
#summary(lmRes)

#anovaRes <- anova(object = lmRes)
#anovaRes$`Pr(>F)`

### `anova(lm())` can be replaced by `aov()`.



##### 2.3. Perform ANOVA for all Metabolomic PCA data in 2017 #####
metabPCAStart <- 11
metabPCAEnd <- ncol(metabPCA2017)
metabPCANames <- colnames(metabPCA2017)[metabPCAStart:metabPCAEnd]
levelNames <- c("variety", "block", "variety x block")

nPCAMetab <- metabPCAEnd - metabPCAStart + 1
nLevel <- 3

pValMetabPCA2017 <- matrix(NA, nrow = nPCAMetab, ncol = nLevel)

for (metabPCANo in metabPCAStart:metabPCAEnd) {
  metabPCANow <- metabPCA2017[, metabPCANo]

  lmRes <- lm(formula = metabPCANow ? variety + block + variety * block,
              data = metabPCA2017)
  anovaRes <- anova(object = lmRes)
  pValMetabPCA2017[metabPCANo - metabPCAStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetabPCA2017) <- metabPCANames
colnames(pValMetabPCA2017) <- levelNames


### `apply`/ `sapply` can be used instead of `for` loop.
pValMetab2017Apply <- t(apply(
  X = metabPCA2017[, metabStart:metabEnd],
  MARGIN = 2,
  FUN = function (metabNow) {
    lmRes <- lm(formula = metabNow ? variety + block + variety * block,
                data = metabPCA2017)
    anovaRes <- anova(object = lmRes)

    return(anovaRes$`Pr(>F)`[1:nLevel])
  }
))

colnames(pValMetab2017Apply) <- levelNames


pValMetab2017Sapply <- t(sapply(
  X = metabStart:metabEnd,
  FUN = function (metabNo) {
    metabNow <- metabPCA2017[, metabNo]

    lmRes <- lm(formula = metabNow ? variety + block + variety * block,
                data = metabPCA2017)
    anovaRes <- anova(object = lmRes)

    return(anovaRes$`Pr(>F)`[1:nLevel])
  }
))

rownames(pValMetab2017Sapply) <- metabNames
colnames(pValMetab2017Sapply) <- levelNames


##### 2.4. Adjust p-values of ANOVA by BH method #####
pAdjMetabPCA2017 <- apply(pValMetabPCA2017, 2,
                       FUN = function(pVals) {
                         pAdjs <- p.adjust(p = pVals, method = "BH")

                         return(pAdjs)
                       })
# pAdjMetab2017 <- apply(pValMetab2017, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetabPCA2017 <- paste0(dirMidSTAMANOVAForPCA, scriptID,
                            "_adjusted_p_values_by_ANOVA_for_PCA_in_2017.csv")
write.csv(x = pAdjMetabPCA2017, file = filePAdjMetabPCA2017)


nDetectedMetabPCA2017 <- apply(pAdjMetabPCA2017 <= sigLevel, 2, sum)
propDetectedMetabPCA2017 <- nDetectedMetabPCA2017 / nPCAMetab

pAdjMetabPCA2017SortedPvalues <- pAdjMetabPCA2017[order(pAdjMetabPCA2017[, 1], decreasing = F), ]










###### 3. Perform ANOVA and adjust p-values for Metabolomic data in 2018 ######
##### 3.1. Read Metabolomic data in 2018 into R #####
metabPCA2018 <- read.csv("midstream/2.3_PCA/2.3_pcaMethods_PCA_metab_2018.csv")
See(metabPCA2018, rown = 6, coln = 12)



##### 3.2. Perform ANOVA for X00006 in Metabolomic data in 2018 #####
# lmRes <- lm(formula = X00006 ? variety + block + variety * block,
#             data = metabPCA2018)
# summary(lmRes)
#
# anovaRes <- anova(object = lmRes)
# anovaRes$`Pr(>F)`
#
# ### `anova(lm())` can be replaced by `aov()`.



##### 3.3. Perform ANOVA for all Metabolomic PCA data in 2018 #####
metabPCAStart <- 11
metabPCAEnd <- ncol(metabPCA2018)
metabPCANames <- colnames(metabPCA2018)[metabPCAStart:metabPCAEnd]
levelNames <- c("variety", "block", "variety x block")

nPCAMetab <- metabPCAEnd - metabPCAStart + 1
nLevel <- 3

pValMetabPCA2018 <- matrix(NA, nrow = nPCAMetab, ncol = nLevel)

for (metabPCANo in metabPCAStart:metabPCAEnd) {
  metabPCANow <- metabPCA2018[, metabPCANo]

  lmRes <- lm(formula = metabPCANow ? variety + block + variety * block,
              data = metabPCA2018)
  anovaRes <- anova(object = lmRes)
  pValMetabPCA2018[metabPCANo - metabPCAStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetabPCA2018) <- metabPCANames
colnames(pValMetabPCA2018) <- levelNames


### `apply`/ `sapply` can be used instead of `for` loop.
pValMetab2018Apply <- t(apply(
  X = metabPCA2018[, metabStart:metabEnd],
  MARGIN = 2,
  FUN = function (metabNow) {
    lmRes <- lm(formula = metabNow ? variety + block + variety * block,
                data = metabPCA2018)
    anovaRes <- anova(object = lmRes)

    return(anovaRes$`Pr(>F)`[1:nLevel])
  }
))

colnames(pValMetab2018Apply) <- levelNames


pValMetab2018Sapply <- t(sapply(
  X = metabStart:metabEnd,
  FUN = function (metabNo) {
    metabNow <- metabPCA2018[, metabNo]

    lmRes <- lm(formula = metabNow ? variety + block + variety * block,
                data = metabPCA2018)
    anovaRes <- anova(object = lmRes)

    return(anovaRes$`Pr(>F)`[1:nLevel])
  }
))

rownames(pValMetab2018Sapply) <- metabNames
colnames(pValMetab2018Sapply) <- levelNames


##### 3.4. Adjust p-values of ANOVA by BH method #####
pAdjMetabPCA2018 <- apply(pValMetabPCA2018, 2,
                          FUN = function(pVals) {
                            pAdjs <- p.adjust(p = pVals, method = "BH")

                            return(pAdjs)
                          })
# pAdjMetab2018 <- apply(pValMetab2018, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetabPCA2018 <- paste0(dirMidSTAMANOVAForPCA, scriptID,
                               "_adjusted_p_values_by_ANOVA_for_PCA_in_2018.csv")
write.csv(x = pAdjMetabPCA2018, file = filePAdjMetabPCA2018)


nDetectedMetabPCA2018 <- apply(pAdjMetabPCA2018 <= sigLevel, 2, sum)
propDetectedMetabPCA2018 <- nDetectedMetabPCA2018 / nPCAMetab

pAdjMetabPCA2018SortedPvalues <- pAdjMetabPCA2018[order(pAdjMetabPCA2018[, 1], decreasing = F), ]


### `order()` can be used for sorting p-values.
pAdjMetabPCA2018SortedPvalues <- pAdjMetab2018[order(pAdjMetab2018[, 1], decreasing = F), ]




### plot pAdjMetab in 2018 and 2018 ###
# pAdjMetab2018 <- read.csv("midstream/2.1_ANOVA/2.1_adjusted_p_values_by_ANOVA_in_2018.csv")
# pAdjMetab2018 <- read.csv("midstream/2.1_ANOVA/2.1_adjusted_p_values_by_ANOVA_in_2018.csv")

# pAdjMetab2018 <- as.matrix(pAdjMetab2018)
# pAdjMetab2018 <- as.matrix(pAdjMetab2018) ?


plot(-log10(pAdjMetabPCA2018[, 1]),
     -log10(pAdjMetabPCA2018[, 1]))
plot(-log10(pAdjMetabPCA2018[, 2]),
     -log10(pAdjMetabPCA2018[, 2]))
plot(-log10(pAdjMetabPCA2018[, 3]),
     -log10(pAdjMetabPCA2018[, 3]))






###### 4. Perform ANOVA and adjust p-values for Metabolomic data using both in 2017 and 2018 ######
##### 4.1. Read Metabolomic data in 2017 and 2018, and integrate them into R #####
metabPCA2017and2018 <- rbind(metabPCA2017, metabPCA2018)



##### 4.2. Perform ANOVA for X00006 in Metabolomic data in 2017 and 2018 #####
# lmRes <- lm(formula = X00006 ? variety + block + variety * block + year,
#             data = metabPCA2017and2018)
# summary(lmRes)
#
# anovaRes <- anova(object = lmRes)
# anovaRes$`Pr(>F)`
#
# ### `anova(lm())` can be replaced by `aov()`.



##### 4.3. Perform ANOVA for all Metabolomic data in 2017 and 2018 #####
metabStart <- 11
metabEnd <- ncol(metabPCA2017and2018)
metabNames <- colnames(metabPCA2017and2018)[metabStart:metabEnd]
levelNames <- c("variety", "block", "variety x block", "year")

nMetab <- metabEnd - metabStart + 1
nLevel <- 4

pValMetabPCA2017and2018 <- matrix(NA, nrow = nMetab, ncol = nLevel)

for (metabNo in metabStart:metabEnd) {
  metabNow <- metabPCA2017and2018[, metabNo]

  lmRes <- lm(formula = metabNow ? variety + block + variety * block + year,
              data = metabPCA2017and2018)
  anovaRes <- anova(object = lmRes)
  pValMetabPCA2017and2018[metabNo - metabStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetabPCA2017and2018) <- metabNames
colnames(pValMetabPCA2017and2018) <- levelNames




##### 4.4. Adjust p-values of ANOVA by BH method #####
pAdjMetabPCA2017and2018 <- apply(pValMetabPCA2017and2018, 2,
                            FUN = function(pVals) {
                              pAdjs <- p.adjust(p = pVals, method = "BH")

                              return(pAdjs)
                            })
# pAdjMetabPCA2017and2018 <- apply(pValMetabPCA2017and2018, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetabPCA2017and2018 <- paste0(dirMidSTAMANOVA, scriptID,
                                   "_adjusted_p_values_by_ANOVA_For_PCA_in_2017_&_2018.csv")
write.csv(x = pAdjMetabPCA2017and2018, file = filePAdjMetabPCA2017and2018)


nDetectedMetabPCA2017and2018 <- apply(pAdjMetabPCA2017and2018 <= sigLevel, 2, sum)
propDetectedMetabPCA2017and2018 <- nDetectedMetabPCA2017and2018 / nMetab



### `order()` can be used for sorting p-values.
pAdjMetabPCA2017and2018SortedPvalues <- pAdjMetabPCA2017and2018[order(pAdjMetabPCA2017and2018[, 1], decreasing = F), ]






