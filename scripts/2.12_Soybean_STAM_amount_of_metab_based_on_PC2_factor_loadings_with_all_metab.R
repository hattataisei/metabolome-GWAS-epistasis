##########################################################################################
######  Title: 2.12_Soybean_STAM_amount_of_metab_based_on_PC2_factor_loadings_with_all_metab                    ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2020/06/07 (Created), 2024/04/05 (Last Updated)                       ######
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

scriptID <- "2.12"


##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"
# sigLevel <- 0.05

dirMidSTAMMetabAmountBasedOnPC2FactorLoadings <- paste0(dirMidSTAMBase, scriptID,
                                "_Amount_of_metab_based_on_PC2/")
dir.create(dirMidSTAMMetabAmountBasedOnPC2FactorLoadings)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)



##### 1.3. Import packages #####
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)



##### 1.4. Project options #####
options(stringAsFactors = FALSE)






###### 2. Check certain metabolites of flavonoid pathway ######
##### 2.1. Read Metabolomic data in 2017 into R #####
gvMetab2017 <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv")

See(gvMetab2017)

#gvMetab2017$X200014[1:10]
#gvMetab2017$X250001[1:10]
#gvMetab2017$X500080[1:10]
#gvMetab2017$X00422[1:10]

Quercetin <- gvMetab2017[, "X250001"]
Kaempferol <- gvMetab2017[, "X00422"]
line_name <- gvMetab2017[, "X"]

#names(Quercetin) <- line_name

See(Quercetin)

#names(Quercetin) <-
#colnamesVariety = gvMetab2017[, 1]
#Quercetin <- cbind(colnamesVariety, Quercetin)

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
markerInterest <- gastonData0Matrix[, "Chr06_18760995"]
See(markerInterest)

#varietyLogical <- markerInterest %in% Quercetin
#varietyLogical

#boxplot(Quercetin ? markerInterest)


# ## Quercetin ##
# length(markerInterest)
# length(Quercetin)
#
# head(Quercetin)
# df_Q <- data.frame(Quercetin, row.names = line_name)
# df_m <- data.frame(markerInterest)
#
# length(df_m[row.names(df_Q), ])
# head(df_m)
#
# length(unique(c(row.names(df_Q), row.names(df_m))))
#
# unique(markerInterest)
# unique(df_m)
# #unique(gastonData0Matrix)
#
# val_mI <- data.frame(df_m[row.names(df_Q), ], row.names = row.names(df_Q))
# val_Q <- data.frame(df_Q[row.names(df_m), ], row.names = row.names(df_m))
#
# omitted_val_mI <- na.omit(val_mI)
# omitted_val_Q <- na.omit(val_Q)
#
# df_merge <- merge(omitted_val_Q, omitted_val_mI, by = 0)
# head(df_merge)
#
# boxplot(df_merge$df_Q.row.names.df_m.... ? df_merge$df_m.row.names.df_Q....)
# summary(df_merge[, 1])
#
# ## Kaempferol ##
# length(markerInterest)
# length(Kaempferol)
#
# head(Kaempferol)
# df_K <- data.frame(Kaempferol, row.names = line_name)
# df_m <- data.frame(markerInterest)
#
# length(df_m[row.names(df_K), ])
# head(df_m)
#
# length(unique(c(row.names(df_K), row.names(df_m))))
#
# unique(markerInterest)
# unique(df_m)
# #unique(gastonData0Matrix)
#
# val_mI <- data.frame(df_m[row.names(df_K), ], row.names = row.names(df_K))
# val_K <- data.frame(df_K[row.names(df_m), ], row.names = row.names(df_m))
#
# omitted_val_mI <- na.omit(val_mI)
# omitted_val_K <- na.omit(val_K)
#
# df_merge <- merge(omitted_val_K, omitted_val_mI, by = 0)
# head(df_merge)
# See(df_merge)
# # df_mergeNames <- df_merge$Row.names
#
# boxplot(df_merge$df_K.row.names.df_m.... ? df_merge$df_m.row.names.df_K....)
#
#
#
# aho <- c(29, 14, 53, 102, 431, 4, 329, 213)
# baka <- c(0, 1, 0, 2, 2, 0, 2, 2)
# boxplot(aho ? baka)
# #baka:genotype







markerInterest <- gastonData0Matrix[, "Chr06_18760995"]
See(markerInterest)
Quercetin <- as.data.frame(gvMetab2017[, "X250001"])
Kaempferol <- as.data.frame(gvMetab2017[, "X00422"])
LineNames <- gvMetab2017[, "X"]

rownames(Quercetin) <- LineNames
colnames(Quercetin) <- "X250001"
rownames(Kaempferol) <- LineNames
colnames(Kaempferol) <- "X00422"

length(markerInterest)
See(markerInterest)
# markerInterest <- as.data.frame(markerInterest)
# markerVarietyNames <- rownames(markerInterest)
markerVarietyNames <- names(markerInterest)
LineNames <- LineNames[(LineNames%in%markerVarietyNames)]
CommonNames <- markerVarietyNames[(markerVarietyNames%in%LineNames)]
table(LineNames%in%markerVarietyNames)
table(markerVarietyNames%in%LineNames)

# df_mergeNames%in%CommonNames
# table(df_mergeNames%in%CommonNames)

Quercetin <- Quercetin[CommonNames, ]
Kaempferol <- Kaempferol[CommonNames, ]
See(Quercetin)
See(Kaempferol)
Quercetin <- as.data.frame(Quercetin)
Kaempferol <- as.data.frame(Kaempferol)
rownames(Quercetin) <- CommonNames
rownames(Kaempferol) <- CommonNames

markerInterest <- as.data.frame(markerInterest)
colnames(markerInterest) <- "Chr06_18760995"
markerInterest <- markerInterest[CommonNames, ]
markerInterest <- as.data.frame(markerInterest)
rownames(markerInterest) <- CommonNames

# Quercetin <- as.vector(Quercetin)
# Kaempferol <- as.vector(Kaempferol)
# markerInterest <- as.vector(markerInterest)
# class(Quercetin)
# boxplot(Quercetin ? markerInterest)


boxplot(Quercetin$Quercetin ? markerInterest$markerInterest)
boxplot(Kaempferol$Kaempferol ? markerInterest$markerInterest)
head(Quercetin)
head(markerInterest)
?boxplot
class(Quercetin$Quercetin)

