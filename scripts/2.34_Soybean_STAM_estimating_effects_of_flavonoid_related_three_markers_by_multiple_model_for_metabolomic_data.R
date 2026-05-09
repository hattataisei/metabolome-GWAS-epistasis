##########################################################################################
######  Title: 2.34_Soybean_STAM_estimating_effects_of_flavonoid_related_three_markers_by_multiple_model_for_metabolomic_data          ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/11/13 (Created), 2024/11/13 (Last Updated)                       ######
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

scriptID <- "2.34"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel <- paste0(dirMidSTAMBase, scriptID,
                                          "_Flavonoid_coefficients_of_three_markers_based_on_multiple_model/")
dir.create(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)




thresLD <- 2



##### 1.3. Import packages #####
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)
require(plotly)
# require(manhattanly)




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

# gastonDataSmall <- LD.thin(gastonData0, threshold = thresLD)
gastonDataSmall <- gastonData0
See(gastonDataSmall)




##### 2.2. Read  Metabolomic data in 2017 into R #####
#### 2.2.1. Genotypic values ####
gvMetab2017TotalInterceptPlusRanef <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_intercept_plus_ranef_2017.csv", row.names = 1)
See(gvMetab2017TotalInterceptPlusRanef)

rownames(gvMetab2017TotalInterceptPlusRanef)[rownames(gvMetab2017TotalInterceptPlusRanef) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"


# gvMetab2017TotalInterceptPlusRanefScaled <- scale(gvMetab2017Total, center = TRUE, scale = TRUE)

# metabLogicalMoreThan0orNA <- apply(gvMetab2017TotalInterceptPlusRanef, 2, function(x)any(x < 0))
# metabNamesNA <- names(metabLogicalMoreThan0orNA[is.na(metabLogicalMoreThan0orNA)])
# metabLogicalMoreThan0 <- metabLogicalMoreThan0orNA
# metabLogicalMoreThan0[metabNamesNA] <- FALSE
# table(metabLogicalMoreThan0)
# gvMetab2017TotalInterceptPlusRanef <- gvMetab2017TotalInterceptPlusRanef[, !metabLogicalMoreThan0]
# gvMetab2017TotalInterceptPlusRanef$X500098
# gvMetab2017TotalInterceptPlusRanef$X00853



####  2.2.2. Raw data ####
metabRaw2017 <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier.csv")
rownames(metabRaw2017)[rownames(metabRaw2017) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
See(metabRaw2017, coln = 10)
table(is.na(metabRaw2017))

metabStart <- 10
metabEnd <- ncol(metabRaw2017)

metabRaw2017MetabOnly <- metabRaw2017[metabStart:metabEnd]
See(metabRaw2017MetabOnly)

# metabRaw2017Scaled <- scale(metabRaw2017MetabOnly, center = TRUE, scale = TRUE)
# metabRaw2017Scaled <- scale(metabRaw2017MetabOnly, center = FALSE, scale = TRUE)
# See(metabRaw2017Scaled)
# apply(metabRaw2017Scaled, 2, min)




##### 2.3. Estimating effects of three markers based on multiple model, using intercept + ranef #####
#### 2.3.1. dividing group into F3H = 0 or 2, and extracting > 0 metabolite names ####
genoMat <- as.matrix(gastonDataSmall)
See(genoMat)

### dividing group by Chr06_18760995 = 0 or 2
## Chr06_18760995 = 0
varietyNamesF3H0 <- rownames(genoMat)[genoMat[, "Chr06_18760995"] == 0]
varietyNamesMetab <- rownames(gvMetab2017TotalInterceptPlusRanef)
varietyNamesCommonF3H0 <- varietyNamesF3H0[varietyNamesF3H0 %in% varietyNamesMetab]

gvMetab2017TotalInterceptPlusRanefF3H0 <- gvMetab2017TotalInterceptPlusRanef[varietyNamesCommonF3H0, ]
See(gvMetab2017TotalInterceptPlusRanefF3H0)

## Chr06_18760995 = 2
varietyNamesF3H2 <- rownames(genoMat)[genoMat[, "Chr06_18760995"] == 2]
varietyNamesMetab <- rownames(gvMetab2017TotalInterceptPlusRanef)
varietyNamesCommonF3H2 <- varietyNamesF3H2[varietyNamesF3H2 %in% varietyNamesMetab]

gvMetab2017TotalInterceptPlusRanefF3H2 <- gvMetab2017TotalInterceptPlusRanef[varietyNamesCommonF3H2, ]
See(gvMetab2017TotalInterceptPlusRanefF3H2)

###


### check whether > 0 or < 0
## Chr06_18760995 = 0
apply(gvMetab2017TotalInterceptPlusRanefF3H0, 2, function(x)any(x < 0))
apply(gvMetab2017TotalInterceptPlusRanefF3H0, 2, function(x)any(is.na(x < 0)))
gvMetab2017TotalInterceptPlusRanefF3H0[, "X500098"]

gvMetab2017TotalInterceptPlusRanefF3H0Logical <- apply(gvMetab2017TotalInterceptPlusRanefF3H0, 2, function(x)any(x < 0))
table(is.na(gvMetab2017TotalInterceptPlusRanefF3H0Logical))
gvMetab2017TotalInterceptPlusRanefF3H0Logical[is.na(gvMetab2017TotalInterceptPlusRanefF3H0Logical)]
gvMetab2017TotalInterceptPlusRanefF3H0Logical["X500098"] <- FALSE

gvMetab2017TotalInterceptPlusRanefF3H0[, gvMetab2017TotalInterceptPlusRanefF3H0Logical]
## X00750, X00959, X200079, X500018 include  many < 0 elements in F3H = 0 varieties ##

gvMetab2017TotalInterceptPlusRanefF3H0[, !gvMetab2017TotalInterceptPlusRanefF3H0Logical]

## Chr06_18760995 = 2
apply(gvMetab2017TotalInterceptPlusRanefF3H2, 2, function(x)any(x < 0))
apply(gvMetab2017TotalInterceptPlusRanefF3H2, 2, function(x)any(is.na(x < 0)))

gvMetab2017TotalInterceptPlusRanefF3H2Logical <- apply(gvMetab2017TotalInterceptPlusRanefF3H2, 2, function(x)any(x < 0))
table(is.na(gvMetab2017TotalInterceptPlusRanefF3H2Logical))
gvMetab2017TotalInterceptPlusRanefF3H2Logical[is.na(gvMetab2017TotalInterceptPlusRanefF3H2Logical)]
gvMetab2017TotalInterceptPlusRanefF3H2Logical[c("X00317", "X00853", "X200005")] <- FALSE

gvMetab2017TotalInterceptPlusRanefF3H2[, gvMetab2017TotalInterceptPlusRanefF3H2Logical]
## X00750, X00959, X200079, X500018 include  many < 0 elements in F3H = 0 varieties ##

gvMetab2017TotalInterceptPlusRanefF3H2[, !gvMetab2017TotalInterceptPlusRanefF3H2Logical]

###



#### 2.3.2. linear regression ####
## Chr06_18760995 = 0
markersMat <- genoMat[, c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")]
See(markersMat)
markersCommonMat <- markersMat[varietyNamesCommonF3H0, ]
See(markersCommonMat)

markersCommonMat[markersCommonMat == 2] <- 1
See(markersCommonMat)

gvMetab2017TotalInterceptPlusRanefF3H0ThreeMarkers <- modify.data(pheno.mat = gvMetab2017TotalInterceptPlusRanefF3H0,
                                          geno.mat = markersCommonMat,
                                          return.ZETA = FALSE,
                                          return.GWAS.format = FALSE)
See(gvMetab2017TotalInterceptPlusRanefF3H0ThreeMarkers$geno.modi)
See(gvMetab2017TotalInterceptPlusRanefF3H0ThreeMarkers$pheno.modi)
genoModiThreeMarkers <- gvMetab2017TotalInterceptPlusRanefF3H0ThreeMarkers$geno.modi
phenoModiThreeMarkers <- gvMetab2017TotalInterceptPlusRanefF3H0ThreeMarkers$pheno.modi
See(phenoModiThreeMarkers)

log(phenoModiThreeMarkers)
min(phenoModiThreeMarkers[, 1])

# phenoModiThreeMarkersPlusMin <- apply(phenoModiThreeMarkers, 2, function(x){x - min(x) + 1})
# apply(phenoModiThreeMarkersPlusMin, 2, min)
# table(apply(phenoModiThreeMarkersPlusMin, 2, min) < 0)

phenoModiThreeMarkersLog <- log(phenoModiThreeMarkers)
table(is.nan(phenoModiThreeMarkersLog))
genoModiThreeMarkersNumeric <- apply(genoModiThreeMarkers, 2, as.numeric)
rownames(genoModiThreeMarkersNumeric) <- rownames(genoModiThreeMarkers)


metabNames <- colnames(phenoModiThreeMarkersLog)
nMetab <- length(metabNames)
markerNames <- c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")
nMarkers <- length(markerNames)

coefMatF3H0 <- matrix(NA, nrow = nMetab, ncol = nMarkers)
rownames(coefMatF3H0) <- metabNames
colnames(coefMatF3H0) <- markerNames


# metabNo <- 1

for( metabNo in 1:nMetab ){
  metabNow <- phenoModiThreeMarkersLog[, metabNo]

  if( all(is.na(metabNow))) {
    coefMatF3H0[metabNo, ] <- NA
  } else {
    lmRes <- lm( metabNow ? genoModiThreeMarkersNumeric)
    # summary(lmRes)
    # See(lmRes)
    # lmRes$coefficients
    # lmRes$coefficients[2:4]
    # See(lmRes$coefficients)
    # names(lmRes$coefficients)

    coefMatF3H0[metabNo, ] <- lmRes$coefficients[2:4]
    # return(coefMat)

  }
}

# coefMatZero <- coefMat  # result with using all varieties

See(coefMatF3H0)
coefMatF3H0


write.csv(coefMat, paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_0.csv"))


## Chr06_18760995 = 2
markersMat <- genoMat[, c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")]
See(markersMat)
markersCommonMat <- markersMat[varietyNamesCommonF3H2, ]
See(markersCommonMat)

markersCommonMat[markersCommonMat == 2] <- 1
See(markersCommonMat)

gvMetab2017TotalInterceptPlusRanefF3H2ThreeMarkers <- modify.data(pheno.mat = gvMetab2017TotalInterceptPlusRanefF3H2,
                                                                  geno.mat = markersCommonMat,
                                                                  return.ZETA = FALSE,
                                                                  return.GWAS.format = FALSE)
See(gvMetab2017TotalInterceptPlusRanefF3H2ThreeMarkers$geno.modi)
See(gvMetab2017TotalInterceptPlusRanefF3H2ThreeMarkers$pheno.modi)
genoModiThreeMarkers <- gvMetab2017TotalInterceptPlusRanefF3H2ThreeMarkers$geno.modi
phenoModiThreeMarkers <- gvMetab2017TotalInterceptPlusRanefF3H2ThreeMarkers$pheno.modi
See(phenoModiThreeMarkers)

log(phenoModiThreeMarkers)
min(phenoModiThreeMarkers[, 1])

# phenoModiThreeMarkersPlusMin <- apply(phenoModiThreeMarkers, 2, function(x){x - min(x) + 1})
# apply(phenoModiThreeMarkersPlusMin, 2, min)
# table(apply(phenoModiThreeMarkersPlusMin, 2, min) < 0)

phenoModiThreeMarkersLog <- log(phenoModiThreeMarkers)
table(is.nan(phenoModiThreeMarkersLog))
genoModiThreeMarkersNumeric <- apply(genoModiThreeMarkers, 2, as.numeric)
rownames(genoModiThreeMarkersNumeric) <- rownames(genoModiThreeMarkers)


metabNames <- colnames(phenoModiThreeMarkersLog)
nMetab <- length(metabNames)
markerNames <- c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")
nMarkers <- length(markerNames)

coefMatF3H2 <- matrix(NA, nrow = nMetab, ncol = nMarkers)
rownames(coefMatF3H2) <- metabNames
colnames(coefMatF3H2) <- markerNames


# metabNo <- 1

for( metabNo in 1:nMetab ){
  metabNow <- phenoModiThreeMarkersLog[, metabNo]

  if( all(is.na(metabNow))) {
    coefMatF3H2[metabNo, ] <- NA
  } else {
    lmRes <- lm( metabNow ? genoModiThreeMarkersNumeric)
    # summary(lmRes)
    # See(lmRes)
    # lmRes$coefficients
    # lmRes$coefficients[2:4]
    # See(lmRes$coefficients)
    # names(lmRes$coefficients)

    coefMatF3H2[metabNo, ] <- lmRes$coefficients[2:4]
    # return(coefMat)

  }
}

# coefMatZero <- coefMat  # result with using all varieties

See(coefMatF3H2)
coefMatF3H2


write.csv(coefMatF3H2, paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_2.csv"))




#### 2.3.3. histgram of all of 188 metabolites ####
### All of 188 metabolites, including 4 metabolites < 0
coefMatF3H0 <- read.csv(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_0.csv"))

coefMatExp <- exp(coefMatF3H0)
coefMatAllMetab <- coefMatExp
coefMatAllMetabNoNA <- na.omit(coefMatAllMetab)
See(coefMatAllMetabNoNA)

# min(coefMatAllMetabNoNA[, 1])
# max(coefMatAllMetabNoNA[, 1])
# min(coefMatAllMetabNoNA[, 2])
# max(coefMatAllMetabNoNA[, 2])
# min(coefMatAllMetabNoNA[, 3])
# max(coefMatAllMetabNoNA[, 3])
# hist(coefMatAllMetabNoNA[, 1], main = "Chr06_47490224", xlab = "Coefficient")
# hist(coefMatAllMetabNoNA[, 2], main = "Chr10_42562665", xlab = "Coefficient")
# hist(coefMatAllMetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")
hist(coefMatAllMetabNoNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient")

hist(coefMatAllMetabNoNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMatAllMetabNoNA[coefMatAllMetabNoNA[, 2] < 10, 2], breaks = seq(0, 10, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMatAllMetabNoNA[coefMatAllMetabNoNA[, 2] < 2, 2], breaks = seq(0, 2, 0.1), main = "Chr10_42562665", xlab = "Coefficient")

hist(coefMatAllMetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/"))
dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/"))

dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_All_188_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_All_188_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatAllMetabNoNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_All_188_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatAllMetabNoNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_All_188_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatAllMetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")
dev.off()



### flavonoid metabolites
metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]

# table(metabNamesCoefRound %in% metabNamesFlavonoid)


coefMatFlavonoid <- coefMatExp[metabNamesFlavonoid, ]
# coefMatFlavonoid <- round(coefMatFlavonoid, 2)
See(coefMatFlavonoid)
coefMatFlavonoidNonNA <- na.omit(coefMatFlavonoid)
See(coefMatFlavonoidNonNA)


hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")
# hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0.4, 2, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
# hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 3, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
# hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Flavonoid_83_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Flavonoid_83_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Flavonoid_83_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Flavonoid_83_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()



### Non-Flavonoid metabolites
metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]
metabNamesNonFlavonoid <- rownames(coefMat)[!(rownames(coefMat) %in% metabNamesFlavonoid)]

# table(metabNamesCoefRound %in% metabNamesFlavonoid)

coefMatNonFlavonoid <- coefMatExp[metabNamesNonFlavonoid, ]
# coefMatNonFlavonoid <- round(coefMatNonFlavonoid, 2)
See(coefMatNonFlavonoid)
coefMatNonFlavonoidNonNA <- na.omit(coefMatNonFlavonoid)
See(coefMatNonFlavonoidNonNA)


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Non_Flavonoid_105_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Non_Flavonoid_105_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatNonFlavonoidNonNA[, 1], main = "Chr06_47490224", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Non_Flavonoid_105_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatNonFlavonoidNonNA[, 2], main = "Chr10_42562665", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_All_188_metabolites/", scriptID, "_Non_Flavonoid_105_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatNonFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()


### all flavonoid metabolites
UUU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
UUD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
UDU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
UDD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
DUU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
DUD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
DDU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
DDD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]

UUU <- round(UUU, 2)
UUD <- round(UUD, 2)
UDU <- round(UDU, 2)
UDD <- round(UDD, 2)
DUU <- round(DUU, 2)
DUD <- round(DUD, 2)
DDU <- round(DDU, 2)
DDD <- round(DDD, 2)

UUU
UUD
UDU
UDD
DUU
DUD
DDU
DDD

See(UUU)
See(UUD)
See(UDU)
See(UDD)
See(DUU)
See(DUD)
See(DDU)
See(DDD)



### excluding " = 1 " flavonoid metabolites
coefMatFlavonoidNonNAReplacedBy1 <- coefMatFlavonoidNonNA
coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.15 & coefMatFlavonoidNonNAReplacedBy1 > 0.85] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.1 & coefMatFlavonoidNonNAReplacedBy1 > 0.9] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.11 & coefMatFlavonoidNonNAReplacedBy1 > 0.89] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.05 & coefMatFlavonoidNonNAReplacedBy1 > 0.95] <- 1
coefMatFlavonoidNonNAReplacedBy1

# apply(coefMatExpRound, 1, function(x)all(x == 1))
# table(apply(coefMatExpRound, 1, function(x)all(x == 1)))
# table(!apply(coefMatExpRound, 1, function(x)all(x == 1)))



# metabNamesCoefRound <- rownames(coefMatExpRound)[!apply(coefMatExpRound, 1,


metabNamesAll1 <-  rownames(coefMatFlavonoidNonNAReplacedBy1)[coefMatFlavonoidNonNAReplacedBy1[, 1] == 1 & coefMatFlavonoidNonNAReplacedBy1[, 2] == 1 & coefMatFlavonoidNonNAReplacedBy1[, 3] == 1]
metabNamesNotAll1 <- rownames(coefMatFlavonoidNonNAReplacedBy1)[!(rownames(coefMatFlavonoidNonNAReplacedBy1) %in% metabNamesAll1)]

coefMatFlavonoidNotAll1 <- coefMatFlavonoidNonNAReplacedBy1[metabNamesNotAll1, ]
See(coefMatFlavonoidNotAll1)
# coefMatFlavonoid1NA <- coefMatFlavonoid
# coefMatFlavonoid1NA[coefMatFlavonoid1NA == 1] <- NA
# NA > 1
# NaN > 1



UUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
# See(UUD1)
# UUD1 <- t(as.matrix(UUD1))
# rownames(UUD1) <- "X00431"
UUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
UDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
UDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
UEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
UEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]

DUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
DUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
DDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
DDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(DEU1)
# DEU1 <- t(as.matrix(DEU1))
# rownames(DEU1) <- "X500118"
DDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
DED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]

EUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
EUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
EUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
# See(EUE1)
# EUE1 <- t(as.matrix(EUE1))
# rownames(EUE1) <- "X00853"
EDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(EDU1)
# EDU1 <- t(as.matrix(EDU1))
# rownames(EDU1) <- "X500127"
EDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
EDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
EEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
EED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
# See(EED1)
# EED1 <- t(as.matrix(EED1))
# rownames(EED1) <- "X100039"
EEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]



# apply(coefMatFlavonoid, 2, function(x){
#   x
# })

groupNames <- c("UUU1", "UUD1", "UUE1", "UDU1", "UDD1", "UDE1", "UEU1", "UED1", "UEE1", "DUU1", "DUD1", "DUE1", "DDU1", "DDD1", "DDE1", "DEU1", "DED1", "DEE1", "EUU1", "EUD1", "EUE1", "EDU1", "EDD1", "EDE1", "EEU1", "EED1", "EEE1")

# i <- 1
for(i in 1:length(groupNames)){
  assign(paste(groupNames[i]), round(eval(parse(text = groupNames[i])), 2))
}

See(UUU1)
See(UUD1)
See(UUE1)
# See(UDU1)
# See(UDD1)
# See(UDE1)
# See(UEU1)
# See(UED1)
See(UEE1)
See(DUU1)
See(DUD1)
See(DUE1)
See(DDU1)
See(DDD1)
See(DDE1)
See(DEU1)
# See(DED1)
See(DEE1)
# See(EUU1)
# See(EUD1)
See(EUE1)
See(EDU1)
# See(EDD1)
See(EDE1)
# See(EEU1)
# See(EED1)
# See(EEE1)



# listMetabNotAll1 <- list(UUU, UUD, UDU, UDD, DUU, DUD, DDU, DDD)
# See(listMetabNotAll1[1])
# See(listMetabNotAll1[[1]])
# as.matrix(listMetabNotAll1[1])
#
# i <- 1
# metabNos <- 0
# for(i in length(listMetabNotAll1)){
#   metabNos <- metabNos + nrow(listMetabNotAll1[[i]])
# }






#### 2.3.4. histgram of 184 metabolites, each of which is > 0 ####
### read "coefMat"
coefMatF3H0 <- read.csv(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_0.csv"), row.names = 1)
See(coefMatF3H0)

coefMatF3H2 <- read.csv(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_2.csv"), row.names = 1)
See(coefMatF3H2)
###


### 184 metabolites, BLUP + ranef > 0
coefMatExp <- exp(coefMatF3H0)
coefMatAllMetab <- coefMatExp
See(coefMatAllMetab)
table(apply(coefMatAllMetab, 1, function(x)any(is.na(x))))
coefMatAllMetabNoNA <- na.omit(coefMatAllMetab)
See(coefMatAllMetabNoNA)
coefMat184MetabNoNA <- coefMatAllMetabNoNA[!gvMetab2017TotalInterceptPlusRanefF3H0Logical, ]
See(coefMat184MetabNoNA)


min(coefMat184MetabNoNA[, 1])
max(coefMat184MetabNoNA[, 1])
min(coefMat184MetabNoNA[, 2])
max(coefMat184MetabNoNA[, 2])
min(coefMat184MetabNoNA[, 3])
max(coefMat184MetabNoNA[, 3])
# hist(coefMat184MetabNoNA[, 1], main = "Chr06_47490224", xlab = "Coefficient")
# hist(coefMat184MetabNoNA[, 2], main = "Chr10_42562665", xlab = "Coefficient")
# hist(coefMat184MetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")
hist(coefMat184MetabNoNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient")

hist(coefMat184MetabNoNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMat184MetabNoNA[coefMat184MetabNoNA[, 2] < 10, 2], breaks = seq(0, 10, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMat184MetabNoNA[coefMat184MetabNoNA[, 2] < 2, 2], breaks = seq(0, 2, 0.1), main = "Chr10_42562665", xlab = "Coefficient")

hist(coefMat184MetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/"))
dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0"))

dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_All_184_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_All_184_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMat184MetabNoNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_All_184_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMat184MetabNoNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_All_184_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMat184MetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")
dev.off()


### flavonoid metabolites
metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]

# table(metabNamesCoefRound %in% metabNamesFlavonoid)


coefMatAllMetab <- coefMatExp
See(coefMatAllMetab)
table(apply(coefMatAllMetab, 1, function(x)any(is.na(x))))
coefMatAllMetabNoNA <- na.omit(coefMatAllMetab)
See(coefMatAllMetabNoNA)
coefMat184MetabNoNA <- coefMatAllMetabNoNA[!gvMetab2017TotalInterceptPlusRanefF3H0Logical, ]
See(coefMat184MetabNoNA)


coefMatFlavonoidNonNA <- coefMat184MetabNoNA[rownames(coefMat184MetabNoNA) %in% metabNamesFlavonoid, ]
# coefMatFlavonoid <- round(coefMatFlavonoid, 2)


See(coefMatFlavonoidNonNA)     # 82 Flaonoid metabolites



hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")
# hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0.4, 2, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
# hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 3, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
# hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Flavonoid_82_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Flavonoid_82_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0, 5, 0.1), main = "Chr06_47490224", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Flavonoid_82_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 45, 0.1), main = "Chr10_42562665", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Flavonoid_82_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()


## Extracting metabolites
rownames(coefMatFlavonoidNonNA)[coefMatFlavonoidNonNA[, 2] > 40] # X00434
coefMatFlavonoidNonNA["X00434", ]

rownames(coefMatFlavonoidNonNA)[coefMatFlavonoidNonNA[, 2] > 10] # X00434, X01219, X500081
coefMatFlavonoidNonNA[c("X00434", "X01219", "X500081"), ]

rownames(coefMatFlavonoidNonNA)[coefMatFlavonoidNonNA[, 2] > 5] # X00434, X00863, X01219, X500081, X500103
coefMatFlavonoidNonNA[c("X00434", "X00863", "X01219", "X500081", "X500103"), ]

rownames(coefMatFlavonoidNonNA)[coefMatFlavonoidNonNA[, 2] > 3] # X00434, X00863, X01219, X500081, X500098, X500099, X500103
coefMatFlavonoidNonNA[c("X00434", "X00863", "X01219", "X500081", "X500098", "X500099","X500103"), ]


### Non-Flavonoid metabolites
coefMatAllMetab <- coefMatExp
See(coefMatAllMetab)
table(apply(coefMatAllMetab, 1, function(x)any(is.na(x))))
coefMatAllMetabNoNA <- na.omit(coefMatAllMetab)
See(coefMatAllMetabNoNA)
coefMat184MetabNoNA <- coefMatAllMetabNoNA[!gvMetab2017TotalInterceptPlusRanefF3H0Logical, ]
See(coefMat184MetabNoNA)


metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]
metabNamesNonFlavonoid <- rownames(coefMat184MetabNoNA)[!(rownames(coefMat184MetabNoNA) %in% metabNamesFlavonoid)]

# table(metabNamesCoefRound %in% metabNamesFlavonoid)


coefMatNonFlavonoidNonNA <- coefMat184MetabNoNA[metabNamesNonFlavonoid, ]
# coefMatNonFlavonoid <- round(coefMatNonFlavonoid, 2)
See(coefMatNonFlavonoidNonNA) # 102 Non-Flavnoid metabolites


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Non_Flavonoid_102_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Non_Flavonoid_102_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatNonFlavonoidNonNA[, 1], main = "Chr06_47490224", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Non_Flavonoid_102_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatNonFlavonoidNonNA[, 2], main = "Chr10_42562665", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_intercept_plus_ranef/", scriptID, "_184_metabolites_>0/", scriptID, "_Non_Flavonoid_102_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatNonFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()


### all flavonoid metabolites
UUU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
UUD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
UDU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
UDD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
DUU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
DUD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
DDU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
DDD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]

UUU <- round(UUU, 2)
UUD <- round(UUD, 2)
UDU <- round(UDU, 2)
UDD <- round(UDD, 2)
DUU <- round(DUU, 2)
DUD <- round(DUD, 2)
DDU <- round(DDU, 2)
DDD <- round(DDD, 2)

UUU
UUD
UDU
UDD
DUU
DUD
DDU
DDD

See(UUU)
See(UUD)
See(UDU)
See(UDD)
See(DUU)
See(DUD)
See(DDU)
See(DDD)



### excluding " = 1 " flavonoid metabolites
coefMatFlavonoidNonNAReplacedBy1 <- coefMatFlavonoidNonNA
coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.15 & coefMatFlavonoidNonNAReplacedBy1 > 0.85] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.1 & coefMatFlavonoidNonNAReplacedBy1 > 0.9] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.11 & coefMatFlavonoidNonNAReplacedBy1 > 0.89] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.05 & coefMatFlavonoidNonNAReplacedBy1 > 0.95] <- 1
coefMatFlavonoidNonNAReplacedBy1

# apply(coefMatExpRound, 1, function(x)all(x == 1))
# table(apply(coefMatExpRound, 1, function(x)all(x == 1)))
# table(!apply(coefMatExpRound, 1, function(x)all(x == 1)))



# metabNamesCoefRound <- rownames(coefMatExpRound)[!apply(coefMatExpRound, 1,


metabNamesAll1 <-  rownames(coefMatFlavonoidNonNAReplacedBy1)[coefMatFlavonoidNonNAReplacedBy1[, 1] == 1 & coefMatFlavonoidNonNAReplacedBy1[, 2] == 1 & coefMatFlavonoidNonNAReplacedBy1[, 3] == 1]
metabNamesNotAll1 <- rownames(coefMatFlavonoidNonNAReplacedBy1)[!(rownames(coefMatFlavonoidNonNAReplacedBy1) %in% metabNamesAll1)]

coefMatFlavonoidNotAll1 <- coefMatFlavonoidNonNAReplacedBy1[metabNamesNotAll1, ]
See(coefMatFlavonoidNotAll1)
# coefMatFlavonoid1NA <- coefMatFlavonoid
# coefMatFlavonoid1NA[coefMatFlavonoid1NA == 1] <- NA
# NA > 1
# NaN > 1



UUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
# See(UUD1)
# UUD1 <- t(as.matrix(UUD1))
# rownames(UUD1) <- "X00431"
UUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
UDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
UDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
UEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
UEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]

DUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
DUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
DDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
DDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(DEU1)
# DEU1 <- t(as.matrix(DEU1))
# rownames(DEU1) <- "X500118"
DDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
DED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]

EUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
EUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
EUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
# See(EUE1)
# EUE1 <- t(as.matrix(EUE1))
# rownames(EUE1) <- "X00853"
EDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(EDU1)
# EDU1 <- t(as.matrix(EDU1))
# rownames(EDU1) <- "X500127"
EDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
EDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
EEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
EED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
# See(EED1)
# EED1 <- t(as.matrix(EED1))
# rownames(EED1) <- "X100039"
EEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]



# apply(coefMatFlavonoid, 2, function(x){
#   x
# })

groupNames <- c("UUU1", "UUD1", "UUE1", "UDU1", "UDD1", "UDE1", "UEU1", "UED1", "UEE1", "DUU1", "DUD1", "DUE1", "DDU1", "DDD1", "DDE1", "DEU1", "DED1", "DEE1", "EUU1", "EUD1", "EUE1", "EDU1", "EDD1", "EDE1", "EEU1", "EED1", "EEE1")

# i <- 1
for(i in 1:length(groupNames)){
  assign(paste(groupNames[i]), round(eval(parse(text = groupNames[i])), 2))
}

See(UUU1)
See(UUD1)
See(UUE1)
# See(UDU1)
# See(UDD1)
# See(UDE1)
# See(UEU1)
# See(UED1)
See(UEE1)
See(DUU1)
See(DUD1)
See(DUE1)
See(DDU1)
See(DDD1)
See(DDE1)
See(DEU1)
# See(DED1)
See(DEE1)
# See(EUU1)
# See(EUD1)
See(EUE1)
See(EDU1)
# See(EDD1)
See(EDE1)
# See(EEU1)
# See(EED1)
# See(EEE1)



# listMetabNotAll1 <- list(UUU, UUD, UDU, UDD, DUU, DUD, DDU, DDD)
# See(listMetabNotAll1[1])
# See(listMetabNotAll1[[1]])
# as.matrix(listMetabNotAll1[1])
#
# i <- 1
# metabNos <- 0
# for(i in length(listMetabNotAll1)){
#   metabNos <- metabNos + nrow(listMetabNotAll1[[i]])
# }


#### 2.3.5. compare effects of three SNPs between F3H = 0 and = 2 for 184 metabolites ####
coefMatF3H0 <- read.csv(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_0.csv"), row.names = 1)
See(coefMatF3H0)

coefMatF3H0Exp <- exp(coefMatF3H0)
coefMatF3H0AllMetab <- coefMatF3H0Exp
See(coefMatF3H0AllMetab)
table(apply(coefMatF3H0AllMetab, 1, function(x)any(is.na(x))))
coefMatF3H0AllMetabNoNA <- na.omit(coefMatF3H0AllMetab)
See(coefMatF3H0AllMetabNoNA)
coefMatF3H0184MetabNoNA <- coefMatF3H0AllMetabNoNA[!gvMetab2017TotalInterceptPlusRanefF3H0Logical, ]
See(coefMatF3H0184MetabNoNA)



coefMatF3H2 <- read.csv(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_coefficient_of_188_metabolites_for_three_markers_given_F3H_2.csv"), row.names = 1)
See(coefMatF3H2)

coefMatF3H2Exp <- exp(coefMatF3H2)
coefMatF3H2AllMetab <- coefMatF3H2Exp
See(coefMatF3H2AllMetab)
table(apply(coefMatF3H2AllMetab, 1, function(x)any(is.na(x))))
coefMatF3H2AllMetabNoNA <- na.omit(coefMatF3H2AllMetab)
See(coefMatF3H2AllMetabNoNA)
coefMatF3H2184MetabNoNA <- coefMatF3H2AllMetabNoNA[!gvMetab2017TotalInterceptPlusRanefF3H2Logical, ]
See(coefMatF3H2184MetabNoNA)



metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]

colVec <- rep("blue1", nrow(coefMatF3H0))
names(colVec) <- rownames(coefMatF3H0)
colVec[metabNamesFlavonoid] <- "orange1"

pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_Chr06_47490224_given_F3H.pdf"))
plot(coefMatF3H0184MetabNoNA[, 1], coefMatF3H2184MetabNoNA[, 1],
     xlim = range(coefMatF3H0184MetabNoNA[, 1], coefMatF3H0184MetabNoNA[, 1], na.rm = TRUE),
     ylim = range(coefMatF3H2184MetabNoNA[, 1], coefMatF3H2184MetabNoNA[, 1], na.rm = TRUE),
     main = "Coefficient of Chr06_47490224",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()


pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_Chr10_42562665_given_F3H.pdf"))
plot(coefMatF3H0184MetabNoNA[, 2], coefMatF3H2184MetabNoNA[, 2],
     xlim = range(coefMatF3H0184MetabNoNA[, 2], coefMatF3H0184MetabNoNA[, 2], na.rm = TRUE),
     ylim = range(coefMatF3H2184MetabNoNA[, 2], coefMatF3H2184MetabNoNA[, 2], na.rm = TRUE),
     main = "Coefficient of Chr10_42562665",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()


pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_Chr17_16065902_given_F3H.pdf"))
plot(coefMatF3H0184MetabNoNA[, 3], coefMatF3H2184MetabNoNA[, 3],
     xlim = range(coefMatF3H0184MetabNoNA[, 3], coefMatF3H0184MetabNoNA[, 3], na.rm = TRUE),
     ylim = range(coefMatF3H2184MetabNoNA[, 3], coefMatF3H2184MetabNoNA[, 3], na.rm = TRUE),
     main = "Coefficient of Chr17_16065902",
     xlab = "Chr06_18760995 = 0",
     ylab = "Chr06_18760995 = 2",
     col = colVec,
     pch = 19,
     cex.main = 1.5,
     cex.lab = 1.4)
abline(0, 1, col = 2, lty = 2, lwd = 1.5)
dev.off()












#### 2.4. Estimating effects of three markers based on multiple model, using raw data ####
genoMat <- as.matrix(gastonDataSmall)
See(genoMat)

### dividing group by Chr06_18760995 = 0 or 2
varietyNamesF3H0 <- rownames(genoMat)[genoMat[, "Chr06_18760995"] == 0]
varietyNamesMetab <- rownames(gvMetab2017TotalScaled)
varietyNamesCommon <- varietyNamesF3H0[varietyNamesF3H0 %in% varietyNamesMetab]

gvMetab2017TotalScaledF3H0 <- gvMetab2017TotalScaled[varietyNamesCommon, ]
See(gvMetab2017TotalScaledF3H0)
###


markersMat <- genoMat[, c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")]
See(markersMat)
markersCommonMat <- markersMat[varietyNamesCommon, ]
See(markersCommonMat)

markersCommonMat[markersCommonMat == 2] <- 1
See(markersCommonMat)

gvMetab2017TotalScaledF3H0ThreeMarkers <- modify.data(pheno.mat = gvMetab2017TotalScaledF3H0,
                                                      geno.mat = markersCommonMat,
                                                      return.ZETA = FALSE,
                                                      return.GWAS.format = FALSE)
See(gvMetab2017TotalScaledF3H0ThreeMarkers$geno.modi)
See(gvMetab2017TotalScaledF3H0ThreeMarkers$pheno.modi)
genoModiThreeMarkers <- gvMetab2017TotalScaledF3H0ThreeMarkers$geno.modi
phenoModiThreeMarkers <- gvMetab2017TotalScaledF3H0ThreeMarkers$pheno.modi
See(phenoModiThreeMarkers)

log(phenoModiThreeMarkers)
min(phenoModiThreeMarkers[, 1])

phenoModiThreeMarkersPlusMin <- apply(phenoModiThreeMarkers, 2, function(x){x - min(x) + 1})
apply(phenoModiThreeMarkersPlusMin, 2, min)
table(apply(phenoModiThreeMarkersPlusMin, 2, min) < 0)

phenoModiThreeMarkersPlusMinLog <- log(phenoModiThreeMarkersPlusMin)

genoModiThreeMarkersNumeric <- apply(genoModiThreeMarkers, 2, as.numeric)
rownames(genoModiThreeMarkersNumeric) <- rownames(genoModiThreeMarkers)


metabNames <- colnames(phenoModiThreeMarkersPlusMinLog)
nMetab <- length(metabNames)
markerNames <- c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")
nMarkers <- length(markerNames)

coefMat <- matrix(NA, nrow = nMetab, ncol = nMarkers)
rownames(coefMat) <- metabNames
colnames(coefMat) <- markerNames


# metabNo <- 1

for( metabNo in 1:nMetab ){
  metabNow <- phenoModiThreeMarkersPlusMinLog[, metabNo]

  if( all(is.na(metabNow))) {
    coefMat[metabNo, ] <- NA
  } else {
    lmRes <- lm( metabNow ? genoModiThreeMarkersNumeric)
    # summary(lmRes)
    # See(lmRes)
    # lmRes$coefficients
    # lmRes$coefficients[2:4]
    # See(lmRes$coefficients)
    # names(lmRes$coefficients)

    coefMat[metabNo, ] <- lmRes$coefficients[2:4]
    # return(coefMat)

  }
}

# coefMatZero <- coefMat  # result with using all varieties

coefMat
coefMatExp <- exp(coefMat)
# function(x)all(x == 1))]



### All of 188 metabolites
coefMatAllMetab <- coefMatExp
coefMatAllMetabNoNA <- na.omit(coefMatAllMetab)

min(coefMatAllMetabNoNA[, 1])
max(coefMatAllMetabNoNA[, 1])
min(coefMatAllMetabNoNA[, 2])
max(coefMatAllMetabNoNA[, 2])
min(coefMatAllMetabNoNA[, 3])
max(coefMatAllMetabNoNA[, 3])
hist(coefMatAllMetabNoNA[, 1], breaks = seq(0.5, 1.9, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
hist(coefMatAllMetabNoNA[, 2], breaks = seq(0.4, 2.9, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMatAllMetabNoNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/"))

dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_All_188_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_All_188_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatAllMetabNoNA[, 1], breaks = seq(0, 3, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_All_188_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatAllMetabNoNA[, 2], breaks = seq(0, 3, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_All_188_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatAllMetabNoNA[, 3], breaks = seq(0, 3, 0.1), main = "Chr17_16065902", xlab = "Coefficient")
dev.off()


### flavonoid metabolites
metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]

# table(metabNamesCoefRound %in% metabNamesFlavonoid)


coefMatFlavonoid <- coefMatExp[metabNamesFlavonoid, ]
# coefMatFlavonoid <- round(coefMatFlavonoid, 2)
See(coefMatFlavonoid)
coefMatFlavonoidNonNA <- na.omit(coefMatFlavonoid)
See(coefMatFlavonoidNonNA)


hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0.4, 2, 0.1), main = "Chr06_47490224", xlab = "Coefficient")
hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 3, 0.1), main = "Chr10_42562665", xlab = "Coefficient")
hist(coefMatFlavonoidNonNA[, 3], main = "Chr17_16065902", xlab = "Coefficient")


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Flavonoid_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Flavonoid_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatFlavonoidNonNA[, 1], breaks = seq(0, 3, 0.1), main = "Chr06_47490224", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Flavonoid_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatFlavonoidNonNA[, 2], breaks = seq(0, 3, 0.1), main = "Chr10_42562665", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Flavonoid_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatFlavonoidNonNA[, 3], breaks = seq(0, 3, 0.1), main = "Chr17_16065902", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()



### Non-Flavonoid metabolites
metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]
metabNamesNonFlavonoid <- rownames(coefMat)[!(rownames(coefMat) %in% metabNamesFlavonoid)]

# table(metabNamesCoefRound %in% metabNamesFlavonoid)

coefMatNonFlavonoid <- coefMatExp[metabNamesNonFlavonoid, ]
# coefMatNonFlavonoid <- round(coefMatNonFlavonoid, 2)
See(coefMatNonFlavonoid)
coefMatNonFlavonoidNonNA <- na.omit(coefMatNonFlavonoid)
See(coefMatNonFlavonoidNonNA)


dir.create(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Non_Flavonoid_metabolites/"))
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Non_Flavonoid_metabolites/", scriptID, "_histgram_of_Chr06_47490224.pdf"))
hist(coefMatNonFlavonoidNonNA[, 1], breaks = seq(0, 3, 0.1), main = "Chr06_47490224", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Non_Flavonoid_metabolites/", scriptID, "_histgram_of_Chr10_42562665.pdf"))
hist(coefMatNonFlavonoidNonNA[, 2], breaks = seq(0, 3, 0.1), main = "Chr10_42562665", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()
pdf(paste0(dirMidSTAMFlavonoidThreeMarkersEffectMultipleModel, scriptID, "_BLUP_min_+_1/", scriptID, "_Non_Flavonoid_metabolites/", scriptID, "_histgram_of_Chr17_16065902.pdf"))
hist(coefMatNonFlavonoidNonNA[, 3], breaks = seq(0, 3, 0.1), main = "Chr17_16065902", xlab = "Coefficient", cex.main = 2, cex.lab = 1.5, cex.axis = 2)
dev.off()


### all flavonoid metabolites
UUU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
UUD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
UDU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
UDD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] > 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
DUU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
DUD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] > 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]
DDU <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] > 1, ]
DDD <- coefMatFlavonoidNonNA[coefMatFlavonoidNonNA[, "Chr06_47490224"] < 1 & coefMatFlavonoidNonNA[, "Chr10_42562665"] < 1 & coefMatFlavonoidNonNA[, "Chr17_16065902"] < 1, ]

UUU <- round(UUU, 2)
UUD <- round(UUD, 2)
UDU <- round(UDU, 2)
UDD <- round(UDD, 2)
DUU <- round(DUU, 2)
DUD <- round(DUD, 2)
DDU <- round(DDU, 2)
DDD <- round(DDD, 2)

UUU
UUD
UDU
UDD
DUU
DUD
DDU
DDD

See(UUU)
See(UUD)
See(UDU)
See(UDD)
See(DUU)
See(DUD)
See(DDU)
See(DDD)



### excluding " = 1 " flavonoid metabolites
coefMatFlavonoidNonNAReplacedBy1 <- coefMatFlavonoidNonNA
coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.1 & coefMatFlavonoidNonNAReplacedBy1 > 0.9] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.11 & coefMatFlavonoidNonNAReplacedBy1 > 0.89] <- 1
# coefMatFlavonoidNonNAReplacedBy1[coefMatFlavonoidNonNAReplacedBy1 < 1.05 & coefMatFlavonoidNonNAReplacedBy1 > 0.95] <- 1
coefMatFlavonoidNonNAReplacedBy1

# apply(coefMatExpRound, 1, function(x)all(x == 1))
# table(apply(coefMatExpRound, 1, function(x)all(x == 1)))
# table(!apply(coefMatExpRound, 1, function(x)all(x == 1)))



# metabNamesCoefRound <- rownames(coefMatExpRound)[!apply(coefMatExpRound, 1,


metabNamesAll1 <-  rownames(coefMatFlavonoidNonNAReplacedBy1)[coefMatFlavonoidNonNAReplacedBy1[, 1] == 1 & coefMatFlavonoidNonNAReplacedBy1[, 2] == 1 & coefMatFlavonoidNonNAReplacedBy1[, 3] == 1]
metabNamesNotAll1 <- rownames(coefMatFlavonoidNonNAReplacedBy1)[!(rownames(coefMatFlavonoidNonNAReplacedBy1) %in% metabNamesAll1)]

coefMatFlavonoidNotAll1 <- coefMatFlavonoidNonNAReplacedBy1[metabNamesNotAll1, ]
See(coefMatFlavonoidNotAll1)
# coefMatFlavonoid1NA <- coefMatFlavonoid
# coefMatFlavonoid1NA[coefMatFlavonoid1NA == 1] <- NA
# NA > 1
# NaN > 1



UUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(UUU1)
# UUU1 <- t(as.matrix(UUU1))
# rownames(UUU1) <- "X01222"
UUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
# See(UUD1)
# UUD1 <- t(as.matrix(UUD1))
# rownames(UUD1) <- "X00431"
UUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
UDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
UDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
UEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
UED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
UEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] > 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]

DUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
DUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
DDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
DDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(DEU1)
# DEU1 <- t(as.matrix(DEU1))
# rownames(DEU1) <- "X500126"
DDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
DED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
DEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] < 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]

EUU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
# See(EUU1)
# EUU1 <- t(as.matrix(EUU1))
# rownames(EUU1) <- "X00853"
EUD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
EUE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] > 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
EDU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
EDD1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
EDE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] < 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]
EEU1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] > 1, ]
EED1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] < 1, ]
# See(EED1)
# EED1 <- t(as.matrix(EED1))
# rownames(EED1) <- "X100039"
EEE1 <- coefMatFlavonoidNotAll1[coefMatFlavonoidNotAll1[, "Chr06_47490224"] == 1 & coefMatFlavonoidNotAll1[, "Chr10_42562665"] == 1 & coefMatFlavonoidNotAll1[, "Chr17_16065902"] == 1, ]



# apply(coefMatFlavonoid, 2, function(x){
#   x
# })

groupNames <- c("UUU1", "UUD1", "UUE1", "UDU1", "UDD1", "UDE1", "UEU1", "UED1", "UEE1", "DUU1", "DUD1", "DUE1", "DDU1", "DDD1", "DDE1", "DEU1", "DED1", "DEE1", "EUU1", "EUD1", "EUE1", "EDU1", "EDD1", "EDE1", "EEU1", "EED1", "EEE1")

# i <- 1
for(i in 1:length(groupNames)){
  assign(paste(groupNames[i]), round(eval(parse(text = groupNames[i])), 2))
}

See(UUU1)
See(UUD1)
# See(UUE1)
# See(UDU1)
# See(UDD1)
See(UDE1)
# See(UEU1)
# See(UED1)
See(UEE1)
See(DUU1)
See(DUD1)
See(DUE1, rown = 7)
See(DDU1)
See(DDD1)
See(DDE1)
See(DEU1)
# See(DED1)
See(DEE1)
See(EUU1)
# See(EUD1)
See(EUE1)
See(EDU1)
# See(EDD1)
See(EDE1)
# See(EEU1)
See(EED1)
# See(EEE1)



# listMetabNotAll1 <- list(UUU, UUD, UDU, UDD, DUU, DUD, DDU, DDD)
# See(listMetabNotAll1[1])
# See(listMetabNotAll1[[1]])
# as.matrix(listMetabNotAll1[1])
#
# i <- 1
# metabNos <- 0
# for(i in length(listMetabNotAll1)){
#   metabNos <- metabNos + nrow(listMetabNotAll1[[i]])
# }



