##########################################################################################
######  Title: 2.1_Soybean_STAM_ANOVA_for_metabolomic_data                          ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                            ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/07 (Created), 2026/01/09 (Last Updated)                       ######
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

scriptID <- "2.1"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
sigLevel <- 0.05

dirMidSTAMANOVA <- paste0(dirMidSTAMBase, scriptID,
                          "_ANOVA/")
dir.create(dirMidSTAMANOVA)
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




###### 2. Perform ANOVA and adjust p-values for No Outlier, raw metabolite data ######
##### 2.1. ANOVA for 2017  #####
#### 2.1.1. Read Metabolomic data in 2017 into R ####
metab2017Raw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
See(metab2017Raw, rown = 6, coln = 12)



#### 2.1.2. Perform ANOVA for X00006 in Metabolomic data in 2017 ####
lmRes <- lm(formula = X00006 ? variety + block + variety * block,
            data = metab2017Raw)
summary(lmRes)

anovaRes <- anova(object = lmRes)
anovaRes$`Pr(>F)`

### `anova(lm())` can be replaced by `aov()`.



#### 2.1.3. Perform ANOVA for all Metabolomic data in 2017 ####
metabStart <- 10
metabEnd <- ncol(metab2017Raw)
metabNames <- colnames(metab2017Raw)[metabStart:metabEnd]
levelNames <- c("variety", "block", "variety x block")

nMetab <- metabEnd - metabStart + 1
nLevel <- 3

pValMetab2017 <- matrix(NA, nrow = nMetab, ncol = nLevel)

for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2017Raw[, metabNo]

  lmRes <- lm(formula = metabNow ? variety + block + variety * block,
              data = metab2017Raw)
  anovaRes <- anova(object = lmRes)
  pValMetab2017[metabNo - metabStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetab2017) <- metabNames
colnames(pValMetab2017) <- levelNames


### `apply`/ `sapply` can be used instead of `for` loop.
pValMetab2017Apply <- t(apply(
  X = metab2017Raw[, metabStart:metabEnd],
  MARGIN = 2,
  FUN = function (metabNow) {
    lmRes <- lm(formula = metabNow ? variety + block + variety * block,
                data = metab2017Raw)
    anovaRes <- anova(object = lmRes)

    return(anovaRes$`Pr(>F)`[1:nLevel])
  }
))

colnames(pValMetab2017Apply) <- levelNames


pValMetab2017Sapply <- t(sapply(
  X = metabStart:metabEnd,
  FUN = function (metabNo) {
    metabNow <- metab2017Raw[, metabNo]

    lmRes <- lm(formula = metabNow ? variety + block + variety * block,
                data = metab2017Raw)
    anovaRes <- anova(object = lmRes)

    return(anovaRes$`Pr(>F)`[1:nLevel])
  }
))

rownames(pValMetab2017Sapply) <- metabNames
colnames(pValMetab2017Sapply) <- levelNames


#### 2.1.4. Adjust p-values of ANOVA by BH method ####
pAdjMetab2017 <- apply(pValMetab2017, 2,
                       FUN = function(pVals) {
                         pAdjs <- p.adjust(p = pVals, method = "BH")

                         return(pAdjs)
                       })
# pAdjMetab2017 <- apply(pValMetab2017, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetab2017 <- paste0(dirMidSTAMANOVA, scriptID,
                            "_adjusted_p_values_by_ANOVA_in_2017.csv")
write.csv(x = pAdjMetab2017, file = filePAdjMetab2017)


nDetectedMetab2017 <- apply(pAdjMetab2017 <= sigLevel, 2, sum)
propDetectedMetab2017 <- nDetectedMetab2017 / nMetab

pAdjMetab2017SortedPvalues <- pAdjMetab2017[order(pAdjMetab2017[, 1], decreasing = F), ]





##### 2.2. ANOVA for 2017 #####
#### 2.2.1. Read Metabolomic data in 2018 into R ####
metab2018Raw <- read.csv("raw_data/phenotype/2018_Tottori_May_Metabolome.csv")
See(metab2018Raw, rown = 6, coln = 12)



##### 2.2.2. Perform ANOVA for X00006 in Metabolomic data in 2018 #####
lmRes <- lm(formula = X00006 ? variety + block + variety * block,
            data = metab2018Raw)
summary(lmRes)

anovaRes <- anova(object = lmRes)
anovaRes$`Pr(>F)`

### `anova(lm())` can be replaced by `aov()`.



##### 2.2.3. Perform ANOVA for all Metabolomic data in 2018 #####
metabStart <- 10
metabEnd <- ncol(metab2018Raw)
metabNames <- colnames(metab2018Raw)[metabStart:metabEnd]
levelNames <- c("variety", "block", "variety x block")

nMetab <- metabEnd - metabStart + 1
nLevel <- 3

pValMetab2018 <- matrix(NA, nrow = nMetab, ncol = nLevel)

for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2018Raw[, metabNo]

  lmRes <- lm(formula = metabNow ? variety + block + variety * block,
              data = metab2018Raw)
  anovaRes <- anova(object = lmRes)
  pValMetab2018[metabNo - metabStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetab2018) <- metabNames
colnames(pValMetab2018) <- levelNames




##### 2.2.4. Adjust p-values of ANOVA by BH method #####
pAdjMetab2018 <- apply(pValMetab2018, 2,
                       FUN = function(pVals) {
                         pAdjs <- p.adjust(p = pVals, method = "BH")

                         return(pAdjs)
                       })
# pAdjMetab2018 <- apply(pValMetab2018, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetab2018 <- paste0(dirMidSTAMANOVA, scriptID,
                            "_adjusted_p_values_by_ANOVA_in_2018.csv")
write.csv(x = pAdjMetab2018, file = filePAdjMetab2018)


nDetectedMetab2018 <- apply(pAdjMetab2018 <= sigLevel, 2, sum)
propDetectedMetab2018 <- nDetectedMetab2018 / nMetab



### `order()` can be used for sorting p-values.
pAdjMetab2018SortedPvalues <- pAdjMetab2018[order(pAdjMetab2018[, 1], decreasing = F), ]





##### 2.3. Perform ANOVA and adjust p-values for Metabolomic data using both in 2017 and 2018 #####
#### 2.3.1. Read Metabolomic data in 2018 into R ####
# metab2017Raw <- read.csv("raw_data/phenotype/2017_Tottori_May_Metabolome.csv")
# metab2018Raw <- read.csv("raw_data/phenotype/2018_Tottori_May_Metabolome.csv")
# See(metab2017Raw, rown = 6, coln = 12)
# See(metab2018Raw, rown = 6, coln = 12)

metab2017_2018Raw <- rbind(metab2017Raw, metab2018Raw)



#### 2.3.2. Perform ANOVA for X00006 in Metabolomic data in 2017 and 2018 #####
lmRes <- lm(formula = X00006 ? variety + block + variety * block + year,
            data = metab2017_2018Raw)
summary(lmRes)

anovaRes <- anova(object = lmRes)
anovaRes$`Pr(>F)`

### `anova(lm())` can be replaced by `aov()`.



##### 2.3.3. Perform ANOVA for all Metabolomic data in 2017 and 2018 #####
metabStart <- 10
metabEnd <- ncol(metab2017_2018Raw)
metabNames <- colnames(metab2017_2018Raw)[metabStart:metabEnd]
levelNames <- c("variety", "block", "variety x block", "year")

nMetab <- metabEnd - metabStart + 1
nLevel <- 4

pValMetab2017_2018 <- matrix(NA, nrow = nMetab, ncol = nLevel)

for (metabNo in metabStart:metabEnd) {

  metabNow <- metab2017_2018Raw[, metabNo]
  lmRes <- lm(formula = metabNow ? variety + block + variety * block + year,
              data = metab2017_2018Raw)
  anovaRes <- anova(object = lmRes)
  pValMetab2017_2018[metabNo - metabStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetab2017_2018) <- metabNames
colnames(pValMetab2017_2018) <- levelNames



##### 2.3.4. Adjust p-values of ANOVA by BH method #####
pAdjMetab2017_2018 <- apply(pValMetab2017_2018, 2,
                            FUN = function(pVals) {
                              pAdjs <- p.adjust(p = pVals, method = "BH")

                              return(pAdjs)
                            })
# pAdjMetab2017_2018 <- apply(pValMetab2017_2018, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetab2017and2018 <- paste0(dirMidSTAMANOVA, scriptID,
                                   "_adjusted_p_values_by_ANOVA_in_2017_&_2018.csv")
write.csv(x = pAdjMetab2017_2018, file = filePAdjMetab2017and2018)


nDetectedMetab2017_2018 <- apply(pAdjMetab2017_2018 <= sigLevel, 2, sum)
propDetectedMetab2017_2018 <- nDetectedMetab2017_2018 / nMetab



### `order()` can be used for sorting p-values.
pAdjMetab2017_2018SortedPvalues <- pAdjMetab2017_2018[order(pAdjMetab2017_2018[, 1], decreasing = F), ]





###### 3. Perform ANOVA and adjust p-values for BoX-Cox metabolite data ######
##### 3.1. ANOVA for Box-Cox metabolite data in 2017 #####
#### 3.1.1. Read Metabolomic data in 2017 into R ####
metab2017Raw <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
See(metab2017Raw, rown = 6, coln = 12)



#### 3.1.2. Perform ANOVA for X00006 in Metabolomic data in 2017 ####
# lmRes <- lm(formula = X00006 ? variety + block + variety * block,
#             data = metab2017Raw)
# summary(lmRes)
#
# anovaRes <- anova(object = lmRes)
# anovaRes$`Pr(>F)`

### `anova(lm())` can be replaced by `aov()`.



#### 3.1.3. Perform ANOVA for all Metabolomic data in 2017 ####
metabStart <- 10
metabEnd <- ncol(metab2017Raw)
metabNames <- colnames(metab2017Raw)[metabStart:metabEnd]
levelNames <- c("variety", "block", "variety x block")

nMetab <- metabEnd - metabStart + 1
nLevel <- 3

pValMetab2017BoxCox <- matrix(NA, nrow = nMetab, ncol = nLevel)
metab2017BoxCox <- metab2017Raw[, 1:metabStart - 1]

for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2017Raw[, metabNo]

  lmResBoxCox <- boxcox(lm(formula = metabNow ? variety + block + variety * block,
              data = metab2017Raw))
  lambda <- lmResBoxCox$x[which.max(lmResBoxCox$y)]
  metabNowBoxCox <- (metabNow ^ lambda - 1) / lambda
  metab2017BoxCox <- cbind(metab2017BoxCox, metabNowBoxCox)
}

colnames(metab2017BoxCox)[metabStart:metabEnd] <- metabNames
See(metab2017BoxCox, coln = 15)

qqnorm(metab2017BoxCox[, 11])

fileNameMetab2017RawBoxCox <- paste0("2017_Tottori_May_Metabolome_BoxCox.csv")
write.csv(x = metab2017BoxCox, file = paste0("data/phenotype/", fileNameMetab2017RawBoxCox))


metab2017BoxCox <- read.csv(paste0("data/phenotype/2017_Tottori_May_Metabolome_BoxCox.csv"), row.names = 1)
See(metab2017BoxCox, coln = 12)


for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2017BoxCox[, metabNo]
  lmResMetabNowBoxCox <- lm(formula = metabNow ? variety + block + variety * block,
                           data = metab2017BoxCox)

  anovaRes <- anova(object = lmResMetabNowBoxCox)
  pValMetab2017BoxCox[metabNo - metabStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetab2017BoxCox) <- metabNames
colnames(pValMetab2017BoxCox) <- levelNames




#### 3.1.4. Adjust p-values of ANOVA by BH method ####
pAdjMetab2017BoxCox <- apply(pValMetab2017BoxCox, 2,
                       FUN = function(pVals) {
                         pAdjs <- p.adjust(p = pVals, method = "BH")

                         return(pAdjs)
                       })
# pAdjMetab2017BoxCox <- apply(pValMetab2017BoxCox, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetab2017BoxCox <- paste0(dirMidSTAMANOVA, scriptID, "_metabolome_BoxCox_adjusted_p_values_by_ANOVA_in_2017.csv")
write.csv(x = pAdjMetab2017BoxCox, file = filePAdjMetab2017BoxCox)


nDetectedMetab2017BoxCox <- apply(pAdjMetab2017BoxCox <= sigLevel, 2, sum)
propDetectedMetab2017BoxCox <- nDetectedMetab2017BoxCox / nMetab

pAdjMetab2017BoxCoxSortedPvalues <- pAdjMetab2017BoxCox[order(pAdjMetab2017BoxCox[, 1], decreasing = F), ]




##### 3.2. ANOVA for Box-Cox metabolite data in 2018 #####
#### 3.2.1. Read Metabolomic data in 2018 into R ####
metab2018Raw <- read.csv("data/phenotype/2018_Tottori_May_Metabolome_No_Outlier.csv")
See(metab2018Raw, rown = 6, coln = 12)



#### 3.2.2. Perform ANOVA for X00006 in Metabolomic data in 2018 ####
# lmRes <- lm(formula = X00006 ? variety + block + variety * block,
#             data = metab2018Raw)
# summary(lmRes)
#
# anovaRes <- anova(object = lmRes)
# anovaRes$`Pr(>F)`

### `anova(lm())` can be replaced by `aov()`.



#### 3.2.3. Perform ANOVA for all Metabolomic data in 2018 ####
metabStart <- 10
metabEnd <- ncol(metab2018Raw)
metabNames <- colnames(metab2018Raw)[metabStart:metabEnd]
levelNames <- c("variety", "block", "variety x block")

nMetab <- metabEnd - metabStart + 1
nLevel <- 3

pValMetab2018BoxCox <- matrix(NA, nrow = nMetab, ncol = nLevel)
metab2018BoxCox <- metab2018Raw[, 1:metabStart - 1]

for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2018Raw[, metabNo]

  lmResBoxCox <- boxcox(lm(formula = metabNow ? variety + block + variety * block,
                           data = metab2018Raw))
  lambda <- lmResBoxCox$x[which.max(lmResBoxCox$y)]
  metabNowBoxCox <- (metabNow ^ lambda - 1) / lambda
  metab2018BoxCox <- cbind(metab2018BoxCox, metabNowBoxCox)
}

colnames(metab2018BoxCox)[metabStart:metabEnd] <- metabNames
See(metab2018BoxCox, coln = 15)

qqnorm(metab2018BoxCox[, 11])

fileNameMetab2018RawBoxCox <- paste0("2018_Tottori_May_Metabolome_BoxCox.csv")
write.csv(x = metab2018BoxCox, file = paste0("data/phenotype/", fileNameMetab2018RawBoxCox))


metab2018BoxCox <- read.csv(paste0("data/phenotype/2018_Tottori_May_Metabolome_BoxCox.csv"), row.names = 1)
See(metab2018BoxCox, coln = 12)


for (metabNo in metabStart:metabEnd) {
  metabNow <- metab2018BoxCox[, metabNo]
  lmResMetabNowBoxCox <- lm(formula = metabNow ? variety + block + variety * block,
                            data = metab2018BoxCox)

  anovaRes <- anova(object = lmResMetabNowBoxCox)
  pValMetab2018BoxCox[metabNo - metabStart + 1, ] <- anovaRes$`Pr(>F)`[1:nLevel]
}
rownames(pValMetab2018BoxCox) <- metabNames
colnames(pValMetab2018BoxCox) <- levelNames




#### 3.2.4. Adjust p-values of ANOVA by BH method ####
pAdjMetab2018BoxCox <- apply(pValMetab2018BoxCox, 2,
                             FUN = function(pVals) {
                               pAdjs <- p.adjust(p = pVals, method = "BH")

                               return(pAdjs)
                             })
# pAdjMetab2018BoxCox <- apply(pValMetab2018BoxCox, 2,
#                        FUN = p.adjust, "BH")

filePAdjMetab2018BoxCox <- paste0(dirMidSTAMANOVA, scriptID, "_metabolome_BoxCox_adjusted_p_values_by_ANOVA_in_2018.csv")
write.csv(x = pAdjMetab2018BoxCox, file = filePAdjMetab2018BoxCox)


nDetectedMetab2018BoxCox <- apply(pAdjMetab2018BoxCox <= sigLevel, 2, sum)
propDetectedMetab2018BoxCox <- nDetectedMetab2018BoxCox / nMetab

pAdjMetab2018BoxCoxSortedPvalues <- pAdjMetab2018BoxCox[order(pAdjMetab2018BoxCox[, 1], decreasing = F), ]








###### 4. plot adjusted p-values in 2017 and 2018 ######
## No outlier data
pAdjMetab2017 <- read.csv("midstream/2.1_ANOVA/2.1_metabolome_BoxCox_adjusted_p_values_by_ANOVA_in_2017.csv", row.names = 1)
pAdjMetab2018 <- read.csv("midstream/2.1_ANOVA/2.1_metabolome_BoxCox_adjusted_p_values_by_ANOVA_in_2018.csv", row.names = 1)
See(pAdjMetab2017)
See(pAdjMetab2018)

pdf(paste0(dirMidSTAMANOVA, scriptID,"_adjusted_p_values_by_ANOVA_for_variety_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017[, 1], 10), -log(pAdjMetab2018[, 1], 10))
dev.off()

pdf(paste0(dirMidSTAMANOVA, scriptID,"_adjusted_p_values_by_ANOVA_for_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017[, 2], 10), -log(pAdjMetab2018[, 2], 10))
dev.off()

pdf(paste0(dirMidSTAMANOVA, scriptID,"_adjusted_p_values_by_ANOVA_for_variety_*_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017[, 3], 10), -log(pAdjMetab2018[, 3], 10))
dev.off()



## BoxCox
pAdjMetab2017BoxCox <- read.csv(paste0("midstream/2.1_ANOVA/2.1_metabolome_BoxCox_adjusted_p_values_by_ANOVA_in_2017.csv"), row.names = 1)
pAdjMetab2018BoxCox <- read.csv(paste0("midstream/2.1_ANOVA/2.1_metabolome_BoxCox_adjusted_p_values_by_ANOVA_in_2018.csv"), row.names = 1)

table(is.na(pAdjMetab2017BoxCox))
table(is.na(pAdjMetab2018BoxCox))

table(pAdjMetab2017BoxCox[, 1] == 0)
table(pAdjMetab2017BoxCox[, 2] == 0)
table(pAdjMetab2017BoxCox[, 3] == 0)
table(pAdjMetab2018BoxCox[, 1] == 0)
table(pAdjMetab2018BoxCox[, 2] == 0)
table(pAdjMetab2018BoxCox[, 3] == 0)

which(pAdjMetab2017BoxCox[, 1] == 0)
which(pAdjMetab2018BoxCox[, 1] == 0)
rownames(pAdjMetab2017BoxCox) == rownames(pAdjMetab2018BoxCox)

pAdjMetab2017BoxCox <- !(pAdjMetab2017BoxCox[, 1] == 0 | pAdjMetab2018BoxCox[, 1] == 0)

min(pAdjMetab2017BoxCox[, 1][!pAdjMetab2017BoxCox[, 1] == 0])
min(pAdjMetab2017BoxCox[, 2][!pAdjMetab2017BoxCox[, 2] == 0])
min(pAdjMetab2017BoxCox[, 3][!pAdjMetab2017BoxCox[, 3] == 0])
min(pAdjMetab2018BoxCox[, 1][!pAdjMetab2018BoxCox[, 1] == 0])
min(pAdjMetab2018BoxCox[, 2][!pAdjMetab2018BoxCox[, 2] == 0])
min(pAdjMetab2018BoxCox[, 3][!pAdjMetab2018BoxCox[, 3] == 0])



# Zero excluded
pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_variety_in_2017_&_2018_zero_excluded.pdf"))
plot(-log(pAdjMetab2017BoxCox[, 1][!(pAdjMetab2017BoxCox[, 1] == 0 | pAdjMetab2018BoxCox[, 1] == 0)], 10), -log(pAdjMetab2018BoxCox[, 1][!(pAdjMetab2017BoxCox[, 1] == 0 | pAdjMetab2018BoxCox[, 1] == 0)], 10))
dev.off()

pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017BoxCoxNoZero[, 2], 10), -log(pAdjMetab2018BoxCoxNoZero[, 2], 10))
dev.off()

pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_variety_*_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017BoxCoxNoZero[, 3], 10), -log(pAdjMetab2018BoxCoxNoZero[, 3], 10))
dev.off()

# not substituted
pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_variety_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017BoxCox[, 1][!pAdjMetab2017BoxCox[, 1] == 0], 10), -log(pAdjMetab2018BoxCox[, 1][!pAdjMetab2018BoxCox[, 1] == 0], 10))
dev.off()

pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017BoxCoxNoZero[, 2], 10), -log(pAdjMetab2018BoxCoxNoZero[, 2], 10))
dev.off()

pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_variety_*_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017BoxCoxNoZero[, 3], 10), -log(pAdjMetab2018BoxCoxNoZero[, 3], 10))
dev.off()




### Extract metabolites with higher adjusted p-values
## variety
metabNamespAdjMoreThan50ForVariety2017 <- rownames(pAdjMetab2017[-log(pAdjMetab2017[, 1], 10) > 50, ])
metabNamespAdjMoreThan50ForVariety2018 <- rownames(pAdjMetab2018[-log(pAdjMetab2018[, 1], 10) > 50, ])
metabNamespAdjMoreThan50ForVariety2017And2018 <- metabNamespAdjMoreThan50ForVariety2017[metabNamespAdjMoreThan50ForVariety2017 %in% metabNamespAdjMoreThan50ForVariety2018]

fileMetabNamesPAdj50ForVariety2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_50_for_variety_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThan50ForVariety2017And2018, file = fileMetabNamesPAdj50ForVariety2017And2018)


metabNamespAdjMoreThan30ForVariety2017 <- rownames(pAdjMetab2017[-log(pAdjMetab2017[, 1], 10) > 30, ])
metabNamespAdjMoreThan30ForVariety2018 <- rownames(pAdjMetab2018[-log(pAdjMetab2018[, 1], 10) > 30, ])
metabNamespAdjMoreThan30ForVariety2017And2018 <- metabNamespAdjMoreThan30ForVariety2017[metabNamespAdjMoreThan30ForVariety2017 %in% metabNamespAdjMoreThan30ForVariety2018]

fileMetabNamesPAdj30ForVariety2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_30_for_variety_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThan30ForVariety2017And2018, file = fileMetabNamesPAdj30ForVariety2017And2018)


## treatment
metabNamespAdjMoreThan50ForTreatment2017 <- rownames(pAdjMetab2017[-log(pAdjMetab2017[, 2], 10) > 50, ])
metabNamespAdjMoreThan50ForTreatment2018 <- rownames(pAdjMetab2018[-log(pAdjMetab2018[, 2], 10) > 50, ])
metabNamespAdjMoreThan50ForTreatment2017And2018 <- metabNamespAdjMoreThan50ForTreatment2017[metabNamespAdjMoreThan50ForTreatment2017 %in% metabNamespAdjMoreThan50ForTreatment2018]

fileMetabNamesPAdj50ForTreatment2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_50_for_treatment_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThan50ForTreatment2017And2018, file = fileMetabNamesPAdj50ForTreatment2017And2018)


metabNamespAdjMoreThan30ForTreatment2017 <- rownames(pAdjMetab2017[-log(pAdjMetab2017[, 2], 10) > 30, ])
metabNamespAdjMoreThan30ForTreatment2018 <- rownames(pAdjMetab2018[-log(pAdjMetab2018[, 2], 10) > 30, ])
metabNamespAdjMoreThan30ForTreatment2017And2018 <- metabNamespAdjMoreThan30ForTreatment2017[metabNamespAdjMoreThan30ForTreatment2017 %in% metabNamespAdjMoreThan30ForTreatment2018]

fileMetabNamesPAdj30ForTreatment2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_30_for_treatment_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThan30ForTreatment2017And2018, file = fileMetabNamesPAdj30ForTreatment2017And2018)


## variety * treatment
metabNamespAdjMoreThan50ForVarietyAndTreatment2017 <- rownames(pAdjMetab2017[-log(pAdjMetab2017[, 3], 10) > 50, ])
metabNamespAdjMoreThan50ForVarietyAndTreatment2018 <- rownames(pAdjMetab2018[-log(pAdjMetab2018[, 3], 10) > 50, ])
metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018 <- metabNamespAdjMoreThan50ForVarietyAndTreatment2017[metabNamespAdjMoreThan50ForVarietyAndTreatment2017 %in% metabNamespAdjMoreThan50ForVarietyAndTreatment2018]

fileMetabNamesPAdj50ForVarietyAndTreatment2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_50_for_variety_&_treatment_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018, file = fileMetabNamesPAdj50ForVarietyAndTreatment2017And2018)


metabNamespAdjMoreThan30ForVarietyAndTreatment2017 <- rownames(pAdjMetab2017[-log(pAdjMetab2017[, 3], 10) > 30, ])
metabNamespAdjMoreThan30ForVarietyAndTreatment2018 <- rownames(pAdjMetab2018[-log(pAdjMetab2018[, 3], 10) > 30, ])
metabNamespAdjMoreThan30ForVarietyAndTreatment2017And2018 <- metabNamespAdjMoreThan30ForVarietyAndTreatment2017[metabNamespAdjMoreThan30ForVarietyAndTreatment2017 %in% metabNamespAdjMoreThan30ForVarietyAndTreatment2018]

fileMetabNamesPAdj30ForVarietyAndTreatment2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_30_for_variety_&_treatment_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThan30ForVarietyAndTreatment2017And2018, file = fileMetabNamesPAdj30ForVarietyAndTreatment2017And2018)


## Box-Cox
# variety



# treatment



# variety * treatment
pAdjMetab2017BoxCoxNoZeroSubstituted <- pAdjMetab2017BoxCox
pAdjMetab2017BoxCoxNoZeroSubstituted[pAdjMetab2017BoxCox == 0] <- 5.035606e-318
pAdjMetab2018BoxCoxNoZeroSubstituted <- pAdjMetab2018BoxCox
pAdjMetab2018BoxCoxNoZeroSubstituted[pAdjMetab2018BoxCox == 0] <- 5.035606e-318


hist(pAdjMetab2017BoxCoxNoZeroSubstituted[, 1], breaks = 1000)
hist(-log(pAdjMetab2017BoxCoxNoZeroSubstituted[, 3]), 10, breaks = 100)
hist(-log(pAdjMetab2018BoxCoxNoZeroSubstituted[, 3]), 10, breaks = 100)



See(pAdjMetab2017BoxCoxNoZeroSubstituted)
pAdjMetab2017BoxCoxNoZeroSubstituted <- as.matrix(pAdjMetab2017BoxCoxNoZeroSubstituted)

logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted <- -log(pAdjMetab2017BoxCoxNoZeroSubstituted[, 3])
q1In2017<- quantile(logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted, 0.25)
q3In2017 <- quantile(logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted, 0.75)
iqrValVarietyTreatment2017 <- IQR(logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted)
significantLevel2017 <-  q3In2017 + 1.5 * iqrValVarietyTreatment2017

significantMetabVarietyTreatment2017BoxCox <- logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted[logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted > significantLevel2017]


See(pAdjMetab2018BoxCoxNoZeroSubstituted)
pAdjMetab2018BoxCoxNoZeroSubstituted <- as.matrix(pAdjMetab2018BoxCoxNoZeroSubstituted)

logpAdjMetabVarietyTreatment2018BoxCoxNoZeroSubstituted <- -log(pAdjMetab2018BoxCoxNoZeroSubstituted[, 3])
q1In2018<- quantile(logpAdjMetabVarietyTreatment2018BoxCoxNoZeroSubstituted, 0.25)
q3In2018 <- quantile(logpAdjMetabVarietyTreatment2018BoxCoxNoZeroSubstituted, 0.75)
iqrValVarietyTreatment2018 <- IQR(logpAdjMetabVarietyTreatment2018BoxCoxNoZeroSubstituted)
significantLevel2018 <-  q3In2018 + 1.5 * iqrValVarietyTreatment2018

significantMetabVarietyTreatment2017BoxCox <- logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted[logpAdjMetabVarietyTreatment2017BoxCoxNoZeroSubstituted > significantLevel2018]


metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2017 <- rownames(pAdjMetab2017BoxCox[-log(pAdjMetab2017BoxCox[, 3], 10) > significantLevel2017, ])
metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2018 <- rownames(pAdjMetab2018BoxCox[-log(pAdjMetab2018BoxCox[, 3], 10) > significantLevel2018, ])
metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2017And2018 <- metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2017[metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2017 %in% metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2018]

fileMetabNamesPAdjSigLevelForVarietyAndTreatment2017And2018 <- paste0(dirMidSTAMANOVA, scriptID, "_Metab_names_padj_over_1.5_IQR_for_variety_&_treatment_2017_&_2018.csv")
write.csv(x = metabNamespAdjMoreThanSigLevelForVarietyAndTreatment2017And2018, file = fileMetabNamesPAdjSigLevelForVarietyAndTreatment2017And2018)



pdf(paste0(dirMidSTAMANOVA, scriptID,"_metabolome_BoxCox_adjusted_p_values_by_ANOVA_for_variety_*_treatment_in_2017_&_2018.pdf"))
plot(-log(pAdjMetab2017BoxCoxNoZeroSubstituted[, 3], 10), -log(pAdjMetab2018BoxCoxNoZeroSubstituted[, 3], 10))
abline(v=significantLevel2017, col="red", lty=2)
abline(h=significantLevel2018, col="red", lty=2)
dev.off()




# -log10P > 50
metabNamespAdjMoreThan50ForVarietyAndTreatment2017BoxCox <- rownames(pAdjMetab2017BoxCoxNoZeroSubstituted[-log(pAdjMetab2017BoxCoxNoZeroSubstituted[, 3], 10) > 50, ])
metabNamespAdjMoreThan50ForVarietyAndTreatment2018BoxCox <- rownames(pAdjMetab2018BoxCoxNoZeroSubstituted[-log(pAdjMetab2018BoxCoxNoZeroSubstituted[, 3], 10) > 50, ])
metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018 <- metabNamespAdjMoreThan50ForVarietyAndTreatment2017BoxCox[metabNamespAdjMoreThan50ForVarietyAndTreatment2017BoxCox %in% metabNamespAdjMoreThan50ForVarietyAndTreatment2018BoxCox]

metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018

# -log10P > 30
metabNamespAdjMoreThan30ForVarietyAndTreatment2017BoxCox <- rownames(pAdjMetab2017BoxCoxNoZeroSubstituted[-log(pAdjMetab2017BoxCoxNoZeroSubstituted[, 3], 10) > 30, ])
metabNamespAdjMoreThan30ForVarietyAndTreatment2018BoxCox <- rownames(pAdjMetab2018BoxCoxNoZeroSubstituted[-log(pAdjMetab2018BoxCoxNoZeroSubstituted[, 3], 10) > 30, ])
metabNamespAdjMoreThan30ForVarietyAndTreatment2017And2018 <- metabNamespAdjMoreThan30ForVarietyAndTreatment2017BoxCox[metabNamespAdjMoreThan30ForVarietyAndTreatment2017BoxCox %in% metabNamespAdjMoreThan30ForVarietyAndTreatment2018BoxCox]

metabNamespAdjMoreThan30ForVarietyAndTreatment2017And2018





### metab for variety, and flavonoid and heritability > 0.9
metabFlavonoid<- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
See(metabFlavonoid)
metabNamesFlavonoid <- metabFlavonoid$Name

table(metabNamespAdjMoreThan50ForVariety2017And2018 %in% metabNamesFlavonoid)
table(metabNamespAdjMoreThan30ForVariety2017And2018 %in% metabNamesFlavonoid)

table(metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018 %in% metabNamesFlavonoid)
table(metabNamespAdjMoreThan30ForVarietyAndTreatment2017And2018 %in% metabNamesFlavonoid)



metabFlavonoidHeritabilityMoreThan0.9 <- read.csv("midstream/2.2_BSH/2.2_Flavonoid_>0.9_heritability_2017.csv")
metabNamesFlavonoidHeritabilityMoreThan0.9 <- metabFlavonoidHeritabilityMoreThan0.9$X
table(metabNamespAdjMoreThan50ForVariety2017And2018 %in% metabNamesFlavonoidHeritabilityMoreThan0.9)


### Treatment and variety * treatment
metabNamespAdjMoreThan50ForTreatmentAndVarietyTreatment2017And2018 <- metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018[metabNamespAdjMoreThan50ForVarietyAndTreatment2017And2018 %in% metabNamespAdjMoreThan50ForTreatment2017And2018]












