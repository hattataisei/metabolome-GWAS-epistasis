##########################################################################################
######  Title: 2.36_Soybean_STAM_bayesian_network_including_flavonoid_related_four_markers_and_metabolomic_data          ######
######  Author: Taisei Hatta (hatta@ut-biomet.org)                                  ######
######  Affiliation: Lab. of Biometry and Bioinformatics, The University of Tokyo   ######
######  Date: 2024/12/2 (Created), 2026/02/22 (Last Updated)                       ######
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

scriptID <- "2.36"



##### 1.2. Setting some parameters #####
dirMidSTAMBase <- "midstream/"

dirMidSTAMBayesianNetwork <- paste0(dirMidSTAMBase, scriptID,
                                                             "_Bayesian_Network/")
dir.create(dirMidSTAMBayesianNetwork)
# fileParamsSTAM <- paste0(dirMidSTAMBase, scriptID,
#                                   "_", project, "_all_parameters.RData")
# save.image(fileParamsprojectName)




thresLD <- 2



##### 1.3. Import packages #####
# install.packages("bnlearn")
# install.packages("igraph")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Rgraphviz")

require(bnlearn)
require(igraph)
require(Rgraphviz)
require(data.table)
require(RAINBOWR)
require(ggplot2)
require(tidyverse)
require(gaston)
require(plotly)
# require(manhattanly)




##### 1.4. Project options #####
options(stringAsFactors = FALSE)





###### 2. Perform Bayesian Network for Metabolomic data and 3 SNPs scores in 2017 ######
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




##### 2.2. Read Metabolomic data in 2017 into R #####
#### 2.2.1. Genotypic values ####
### All metab
gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
See(gvMetab2017Total)

rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"



### flavonoid, 83 metabolites
gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"

metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]

gvMetab2017TotalFlavonoid <- gvMetab2017Total[, metabNamesFlavonoid]
See(gvMetab2017TotalFlavonoid)


### flavonoid, 40 metabolites, heritability > 0.9
gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"

metab2017RawHeritabilityMoreThan0.9All <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metab2017RawHeritabilityMoreThan0.9All, coln = 11)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metab2017RawHeritabilityMoreThan0.9All[, 11:50])

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9 <- gvMetab2017Total[, metabNamesFlavonoidHeritabilityMoreThan0.9]
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9)


### PCs, 40 metabolites, heritability > 0.9
gvPCScore2017Flavonoid40Metab <- read.csv("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_2017_pcaMethods_nPC=6.csv", row.names = 1)
See(gvPCScore2017Flavonoid40Metab)

rownames(gvPCScore2017Flavonoid40Metab)[rownames(gvPCScore2017Flavonoid40Metab) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"


table(is.na(gvPCScore2017Flavonoid40Metab))



##### 2.3. Bayesian Network using bnlearn #####
#### 2.3.1. Binding genotypic values and 4 SNPs scores ####
genoMat <- as.matrix(gastonDataSmall)


gvMetab2017TotalNoNA <- na.omit(gvMetab2017Total)
See(gvMetab2017TotalNoNA)
varietyNamesMetab <- rownames(gvMetab2017TotalNoNA)
varietyNamesGeno <- rownames(genoMat)
varietyNamesCommon <- varietyNamesMetab[varietyNamesMetab %in% varietyNamesGeno]
See(varietyNamesCommon)

gvMetab2017TotalCommon <- gvMetab2017Total[varietyNamesCommon, ]
gvMetab2017TotalCommon <- apply(gvMetab2017TotalCommon, 2, as.numeric)
See(gvMetab2017TotalCommon)

genoMatCommon <- genoMat[varietyNamesCommon, ]
fourSNPsScoresCommon <- genoMatCommon[, c("Chr06_18760995", "Chr06_47490224", "Chr10_42562665", "Chr17_16065902")]
See(fourSNPsScoresCommon)



# 188 metabolites
gvMetab2017TotalFourSNPs <- cbind(gvMetab2017TotalCommon, fourSNPsScoresCommon)
See(gvMetab2017TotalFourSNPs)

# PC Score, based on 40 metabolites
gvPCScore2017Flavonoid40MetabCommon <- gvPCScore2017Flavonoid40Metab[varietyNamesCommon, ]
gvPCScore2017Flavonoid40MetabFourSNPs <- cbind(gvPCScore2017Flavonoid40MetabCommon, fourSNPsScoresCommon)




# 188 metabolites, given F3H
fourSNPsScores0Common <- fourSNPsScoresCommon[fourSNPsScoresCommon[, "Chr06_18760995"] == 0, ]
fourSNPsScores2Common <- fourSNPsScoresCommon[fourSNPsScoresCommon[, "Chr06_18760995"] == 2, ]
See(fourSNPsScores0Common)
See(fourSNPsScores2Common)
threeSNPSScoresGivenF3H0Common <- fourSNPsScores0Common[, 2:4]
threeSNPSScoresGivenF3H2Common <- fourSNPsScores2Common[, 2:4]

gvMetab2017TotalF3H0Common <- gvMetab2017Total[rownames(threeSNPSScoresGivenF3H0Common), ]
gvMetab2017TotalF3H2Common <- gvMetab2017Total[rownames(threeSNPSScoresGivenF3H2Common), ]
See(gvMetab2017TotalF3H0Common)
See(gvMetab2017TotalF3H2Common)


# 40 metabolites, given F3H
fourSNPsScoresF3H0Common <- fourSNPsScoresCommon[fourSNPsScoresCommon[, "Chr06_18760995"] == 0, ]
fourSNPsScoresF3H2Common <- fourSNPsScoresCommon[fourSNPsScoresCommon[, "Chr06_18760995"] == 2, ]
See(fourSNPsScores0Common)
See(fourSNPsScores2Common)
threeSNPSScoresGivenF3H0Common <- fourSNPsScoresF3H0Common[, 2:4]
threeSNPSScoresGivenF3H2Common <- fourSNPsScoresF3H2Common[, 2:4]

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H0Common <- gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9[rownames(threeSNPSScoresGivenF3H0Common), ]
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H2Common <- gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9[rownames(threeSNPSScoresGivenF3H2Common), ]
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H0Common)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H2Common)



# 40 metabolites, given F3H and Chr06_47490224
fourSNPsScoresF3H_0_Chr06_47490224_0_Common <- fourSNPsScoresCommon[fourSNPsScoresCommon[, "Chr06_18760995"] == 0 & fourSNPsScoresCommon[, "Chr06_47490224"] == 0, ]
# fourSNPsScores2Common <- fourSNPsScoresCommon[fourSNPsScoresCommon[, "Chr06_18760995"] == 2, ]

See(fourSNPsScoresF3H_0_Chr06_47490224_0_Common)
# See(fourSNPsScoresF3H_0_Chr06_47490224_2_Common)

twoSNPsScoresGivenF3H_0_Chr06_47490224_0_Common <- fourSNPsScoresF3H_0_Chr06_47490224_0_Common[, 3:4]
# threeSNPSScoresGivenF3H_0_Chr06_47490224_2_Common <- fourSNPsScores2Common[, 2:4]
See(twoSNPsScoresGivenF3H_0_Chr06_47490224_0_Common)

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H_0_Chr06_47490224_0_Common <- gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9[rownames(twoSNPsScoresGivenF3H_0_Chr06_47490224_0_Common), ]
# gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H2Common <- gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9[rownames(threeSNPSScoresGivenF3H2Common), ]
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H_0_Chr06_47490224_0_Common)
# See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H2Common)






See(gvMetab2017Total)
table(is.na(gvMetab2017Total))
apply(gvMetab2017Total, 2, function(x)any(is.na(x)))
rownames(gvMetab2017Total)[apply(gvMetab2017Total, 2, function(x)any(is.na(x)))]
colnames(gvMetab2017Total)[apply(gvMetab2017Total, 2, function(x)any(is.na(x)))]
gvMetab2017Total[, "X500098"]
gvMetab2017Total[, "X200005"]

gvMetab2017Total[, "X00682"] <- as.numeric(gvMetab2017Total[, "X00682"])
gvMetab2017Total[, "X500019"] <- as.numeric(gvMetab2017Total[, "X500019"])
gvMetab2017Total <- na.omit(gvMetab2017Total)
See(gvMetab2017Total)


gvMetab2017TotalFourSNPs <- na.omit(gvMetab2017TotalFourSNPs)
gvMetab2017TotalFourSNPs[, "X00682"] <- as.numeric(gvMetab2017TotalFourSNPs[, "X00682"])
gvMetab2017TotalFourSNPs[, "X500019"] <- as.numeric(gvMetab2017TotalFourSNPs[, "X500019"])
gvMetab2017TotalFourSNPs[, "Chr06_18760995"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr06_18760995"])
gvMetab2017TotalFourSNPs[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr06_47490224"])
gvMetab2017TotalFourSNPs[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr10_42562665"])
gvMetab2017TotalFourSNPs[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr17_16065902"])

gvMetab2017Total188MetabFourSNPs <- gvMetab2017TotalFourSNPs






#### 2.3.2. Metabolites ####
betaMat0 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_background_Chr06_18760995/2.32_Chr06_18760995_0.csv", row.names = 1)
betaMat2 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_background_Chr06_18760995/2.32_Chr06_18760995_2.csv", row.names = 1)

dirTwoMarkers2017 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/")
betaMat0_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,0).csv"), row.names = 1)
betaMat0_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,2).csv"), row.names = 1)
betaMat2_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,0).csv"), row.names = 1)
betaMat2_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,2).csv"), row.names = 1)

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
betaMat0Flavonoid <- betaMat0[metabNamesFlavonoid, ]
betaMat2Flavonoid <- betaMat2[metabNamesFlavonoid, ]


## only metab
# all of 188 metabolites
dagGvMetab2017Total <- hc(gvMetab2017Total, score = "bic-g")


# Flavonoid 83 metabolites
dagGvMetab2017Total <- hc(gvMetab2017TotalFlavonoid, score = "bic-g")




# Flavonoid 40 metabolites, heritability > 0.9




#### metabolites, with 4 SNPs ####
# all of 188 metabolites
gvMetab2017Total188MetabFourSNPs <- cbind(gvMetab2017TotalCommon, fourSNPsScoresCommon)
See(gvMetab2017Total188MetabFourSNPs[, 188:192])

gvMetab2017TotalFourSNPs[, "Chr06_18760995"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr06_18760995"])
gvMetab2017TotalFourSNPs[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr06_47490224"])
gvMetab2017TotalFourSNPs[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr10_42562665"])
gvMetab2017TotalFourSNPs[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFourSNPs[, "Chr17_16065902"])
# i <- 189
# for(i in 189:192){
#   gvMetab2017Total188MetabFourSNPs[, i] <- as.factor(gvMetab2017Total188MetabFourSNPs[, i])
# }

gvMetab2017Total188MetabFourSNPs <- as.data.frame(gvMetab2017Total188MetabFourSNPs)
?hc # if containing categorical variables, score = "bic-cg"
dagGvMetab2017Total188MetabFourSNPs <-  hc(gvMetab2017Total188MetabFourSNPs, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_All_of_188_metabolites.pdf"))
graphviz.plot(dagGvMetab2017Total188MetabFourSNPs, shape = "ellipse")
dev.off()


# Flavonoid 83 metabolites
gvMetab2017TotalFlavonoidCommon <- gvMetab2017TotalFlavonoid[varietyNamesCommon, ]
gvMetab2017TotalFlavonoidMetabFourSNPs <- cbind(gvMetab2017TotalFlavonoidCommon, fourSNPsScoresCommon)

gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_18760995"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_18760995"])
gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_47490224"])
gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr10_42562665"])
gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr17_16065902"])
### without bl or wl
See(gvMetab2017TotalFlavonoidMetabFourSNPs)
dagGvMetab2017TotalFlavonoidMetabFourSNPs <- hc(gvMetab2017TotalFlavonoidMetabFourSNPs, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_Flavonoid_83_metabolites.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidMetabFourSNPs, shape = "ellipse")
dev.off()
### with bl
bl
bl1 <- tiers2blacklist(list("Chr06_18760995", c("Chr06_47490224", "Chr10_42562665", "Chr17_16065902")))
bl2 <- tiers2blacklist(list("Chr06_47490224", c("Chr06_18760995", "Chr10_42562665", "Chr17_16065902")))
bl3 <- tiers2blacklist(list("Chr10_42562665", c("Chr06_18760995", "Chr06_47490224", "Chr17_16065902")))
bl4 <- tiers2blacklist(list("Chr17_16065902", c("Chr06_18760995", "Chr06_47490224", "Chr10_42562665")))

blFlavonoidMetabFourSNPs <- rbind(bl1, bl2, bl3, bl4)



dagGvMetab2017TotalFlavonoidMetabFourSNPs <- hc(gvMetab2017TotalFlavonoidMetabFourSNPs, score = "bic-cg", blacklist = blFlavonoidMetabFourSNPs)
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_Flavonoid_83_metabolites_with_blacklist_between_4_SNPs.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidMetabFourSNPs, shape = "ellipse")
dev.off()



# Flavonoid 40 metabolites, heritability > 0.9
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9Common <- gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9[varietyNamesCommon, ]
gvMetab2017TotalFlavonoidMetabFourSNPs <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9Common, fourSNPsScoresCommon)
See(gvMetab2017TotalFlavonoidMetabFourSNPs)

gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_18760995"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_18760995"])
gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr06_47490224"])
gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr10_42562665"])
gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFlavonoidMetabFourSNPs[, "Chr17_16065902"])
See(gvMetab2017TotalFlavonoidMetabFourSNPs)

dagGvMetab2017TotalFlavonoidMetabFourSNPs <- hc(gvMetab2017TotalFlavonoidMetabFourSNPs, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_Flavonoid_40_metabolites_heritability_>_0.9_without_blacklist_between_4SNPs.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidMetabFourSNPs, shape = "ellipse")
dev.off()


metabNamesChildrenChr06_18760995 <- children(dagGvMetab2017TotalFlavonoidMetabFourSNPs, "Chr06_18760995")
metabNamesChildrenChr06_47490224 <- children(dagGvMetab2017TotalFlavonoidMetabFourSNPs, "Chr06_47490224")
metabNamesChildrenChr10_42562665 <- children(dagGvMetab2017TotalFlavonoidMetabFourSNPs, "Chr10_42562665")
metabNamesChildrenChr17_16065902 <- children(dagGvMetab2017TotalFlavonoidMetabFourSNPs, "Chr17_16065902")

metabNamesChildrenChr06_18760995
metabNamesChildrenChr06_47490224
metabNamesChildrenChr10_42562665
metabNamesChildrenChr17_16065902


### comparison children metabolites with the metabolites detected in (given) GWAS
# metabNamesDetectedInGWASAroundChr06_18760995 <- c("X00414", "X00419", "X00422", "X00425", "X00431", "X00432", "X00854", "X00928", "X00931", "X01128", "X200004", "X200012", "X200014", "X200086", "X200087","X210004", "X210010", "X250001", "X250010", "X500080", "X500097", "X500134")

# table(metabNamesChildrenChr06_18760995 %in% metabNamesDetectedInGWASAroundChr06_18760995)


metabNamesFlavonoidAllDetectedThresSd3BetaMat0And2Chr06_47490224 <- read.csv(paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_metabNames_Flavonoid_All_Detected_ThresSd3_BetaMat_0_And_2_Chr06_47490224.csv"), row.names = 1)

metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr10_42562665 <- read.csv(paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_metabNames_Flavonoid_All_Detected_ThresSd3_given_0_0_And_0_2_Chr10_42562665.csv"), row.names = 1)

metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr17_16065902 <- read.csv(paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_metabNames_Flavonoid_All_Detected_ThresSd3_given_0_0_And_0_2_Chr17_16065902.csv"), row.names = 1)



table(metabNamesChildrenChr06_47490224 %in% metabNamesFlavonoidAllDetectedThresSd3BetaMat0And2Chr06_47490224[, 1])
table(metabNamesChildrenChr10_42562665 %in% metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr10_42562665[, 1])
table(metabNamesChildrenChr17_16065902 %in% metabNamesFlavonoidAllDetectedThresSd3BetaMat0_0And0_2Chr17_16065902[, 1])


###


betaMat0 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_Coefficient_LD1_scaled_background_Chr06_18760995/2.32_Chr06_18760995_0.csv", row.names = 1)
betaMat2 <- read.csv("midstream/2.32_GWAS_for_three_markers/2.32_Coefficient_LD1_scaled_background_Chr06_18760995/2.32_Chr06_18760995_2.csv", row.names = 1)
# metabNamesChr06_18760995Score20LinkedChr06_47490224 <- rownames(betaMat0)[betaMat0[, 1] < -0.3 & abs(betaMat2[, 1]) <0.1]
metabNamesChr06_18760995Score20LinkedChr06_47490224 <- rownames(betaMat0)[betaMat0[, 1] < -0.5 & abs(betaMat2[, 1]) <0.1]
metabNamesChr06_18760995Score20LinkedChr06_47490224 <- c(metabNamesChr06_18760995Score20LinkedChr06_47490224, rownames(betaMat0)[betaMat0[, 1] > 0.5])

metabNamesChr06_18760995Score00LinkedChr06_47490224 <- rownames(betaMat0)[betaMat0[, 1] < 0 & betaMat0[, 1] > -0.5 & betaMat2[, 1] < -0.7]
metabNamesChr06_18760995Score20And00LinkedChr06_47490224 <- c(metabNamesChr06_18760995Score20LinkedChr06_47490224, metabNamesChr06_18760995Score00LinkedChr06_47490224)
metabNamesChr06_18760995Score20And00LinkedChr06_47490224 <- na.omit(metabNamesChr06_18760995Score20And00LinkedChr06_47490224)

table(metabNamesChildrenChr06_47490224 %in% metabNamesChr06_18760995Score20And00LinkedChr06_47490224)
table(metabNamesChr06_18760995Score20And00LinkedChr06_47490224 %in% metabNamesChildrenChr06_47490224)
table(metabNamesChr06_18760995Score20LinkedChr06_47490224 %in% metabNamesChildrenChr06_47490224)




children(dagGvMetab2017TotalFlavonoidMetabFourSNPs, "Chr06_47490224")
metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
betaMat0Flavonoid <- betaMat0[metabNamesFlavonoid, ]
betaMat2Flavonoid <- betaMat2[metabNamesFlavonoid, ]

# metabNames20 <- rownames(betaMat0)[betaMat0[, 1] < -0.3 & abs(betaMat2[, 1]) <0.1]
metabNames00 <- rownames(betaMat0)[betaMat0[, 1] < 0 & betaMat0[, 1] > -0.5 & betaMat2[, 1] < -0.7]
metabNames20 <- rownames(betaMat0Flavonoid)[betaMat0Flavonoid[, 1] < -0.5 & abs(betaMat2Flavonoid[, 1]) <0.1]
metabNames00And20 <- c(metabNames20, metabNames00, "X00431")
metabNames00And20 <- na.omit(metabNames00And20)
metabNames00And20

metabNames40Metabolites4SNPsChildrenChr06_47490224 <- children(dagGvMetab2017TotalFlavonoidMetabFourSNPs, "Chr06_47490224")
metabNames40Metabolites4SNPsChildrenChr06_47490224 <- which(metabNames40Metabolites4SNPsChildrenChr06_47490224 == c("Chr10_42562665", "Chr17_16065902"))
metabNames40Metabolites4SNPsChildrenChr06_47490224 <- metabNames40Metabolites4SNPsChildrenChr06_47490224[!metabNames40Metabolites4SNPsChildrenChr06_47490224 %in% c("Chr10_42562665", "Chr17_16065902")]

table(metabNames40Metabolites4SNPsChildrenChr06_47490224 %in% metabNames00And20)
table(metabNames00And20 %in% metabNames40Metabolites4SNPsChildrenChr06_47490224)


## PC loadings
metabNamesChildrenChr06_18760995
metabNamesChildrenChr06_47490224
metabNamesChildrenChr10_42562665
metabNamesChildrenChr17_16065902

table(metabNamesChildrenChr10_42562665 %in% metabNamesChildrenChr17_16065902)
table(metabNamesChildrenChr17_16065902 %in% metabNamesChildrenChr10_42562665)

PC1FactorLoadings <- read.csv("midstream/2.28_Factor_loadings_and_about_PCs_with_high_heritability_flavonoid/2.28_factor_loading_for_PCs/2.28_MetabInfo_of_PC1_factor_loadings_for_flavonoid_heritability_>0.9.csv", row.names = 1)
PC2FactorLoadings <- read.csv("midstream/2.28_Factor_loadings_and_about_PCs_with_high_heritability_flavonoid/2.28_factor_loading_for_PCs/2.28_MetabInfo_of_PC2_factor_loadings_for_flavonoid_heritability_>0.9.csv", row.names = 1)
PC3FactorLoadings <- read.csv("midstream/2.28_Factor_loadings_and_about_PCs_with_high_heritability_flavonoid/2.28_factor_loading_for_PCs/2.28_MetabInfo_of_PC3_factor_loadings_for_flavonoid_heritability_>0.9.csv", row.names = 1)
PC4FactorLoadings <- read.csv("midstream/2.28_Factor_loadings_and_about_PCs_with_high_heritability_flavonoid/2.28_factor_loading_for_PCs/2.28_MetabInfo_of_PC4_factor_loadings_for_flavonoid_heritability_>0.9.csv", row.names = 1)
See(PC1FactorLoadings)
# PC1FactorLoadings <- as.matrix(PC1FactorLoadings)
PC1OnlyFactorLoadings <- PC1FactorLoadings[, 2]
PC2OnlyFactorLoadings <- PC2FactorLoadings[, 2]
PC3OnlyFactorLoadings <- PC3FactorLoadings[, 2]
PC4OnlyFactorLoadings <- PC4FactorLoadings[, 2]
names(PC1OnlyFactorLoadings) <- rownames(PC1FactorLoadings)
names(PC2OnlyFactorLoadings) <- rownames(PC2FactorLoadings)
names(PC3OnlyFactorLoadings) <- rownames(PC3FactorLoadings)
names(PC4OnlyFactorLoadings) <- rownames(PC4FactorLoadings)
PC1OnlyFactorLoadings <- as.matrix(PC1OnlyFactorLoadings)
PC2OnlyFactorLoadings <- as.matrix(PC2OnlyFactorLoadings)
PC3OnlyFactorLoadings <- as.matrix(PC3OnlyFactorLoadings)
PC4OnlyFactorLoadings <- as.matrix(PC4OnlyFactorLoadings)
metabNamesPC1HigherFactorLoadings <- names(PC1OnlyFactorLoadings[order(abs(PC1OnlyFactorLoadings), decreasing = T), ])
metabNamesPC2HigherFactorLoadings <- names(PC2OnlyFactorLoadings[order(abs(PC2OnlyFactorLoadings), decreasing = T), ])
metabNamesPC3HigherFactorLoadings <- names(PC3OnlyFactorLoadings[order(abs(PC3OnlyFactorLoadings), decreasing = T), ])
metabNamesPC4HigherFactorLoadings <- names(PC4OnlyFactorLoadings[order(abs(PC4OnlyFactorLoadings), decreasing = T), ])

metabNamesPC1Top27FactorLoadings <- metabNamesPC1HigherFactorLoadings[1:27]
metabNamesPC2Top17FactorLoadings <- metabNamesPC2HigherFactorLoadings[1:17]
metabNamesPC3Top21FactorLoadings <- metabNamesPC3HigherFactorLoadings[1:21]
metabNamesPC4Top6FactorLoadings <- metabNamesPC4HigherFactorLoadings[1:6]
metabNamesPC1Top27FactorLoadings
metabNamesPC2Top17FactorLoadings
metabNamesPC3Top21FactorLoadings
metabNamesPC4Top6FactorLoadings

table(metabNamesChildrenChr06_18760995 %in% metabNamesPC1Top27FactorLoadings)
table(metabNamesChildrenChr06_47490224 %in% metabNamesPC2Top17FactorLoadings)
table(metabNamesChildrenChr10_42562665 %in% metabNamesPC3Top21FactorLoadings)
table(metabNamesChildrenChr17_16065902 %in% metabNamesPC4Top6FactorLoadings)




start <- Sys.time()
bootStrGvMetab2017TotalFlavonoidgvMetab2017TotalFlavonoidMetabFourSNPs <-  boot.strength(gvMetab2017TotalFlavonoidMetabFourSNPs,
                                                                                                   R = 200,
                                                                                                   algorithm = "hc",
                                                                                                   algorithm.args = list(score="bic-cg"))
end <- Sys.time()

# ?boot.strength
head(bootStrGvMetab2017TotalFlavonoidgvMetab2017TotalFlavonoidMetabFourSNPs)
plot(bootStrGvMetab2017TotalFlavonoidgvMetab2017TotalFlavonoidMetabFourSNPs)
avg.diff = averaged.network(bootStrGvMetab2017TotalFlavonoidgvMetab2017TotalFlavonoidMetabFourSNPs)

pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_strength_40_metabolites_and_4SNPs.pdf"))
strength.plot(avg.diff, bootStrGvMetab2017TotalFlavonoidgvMetab2017TotalFlavonoidMetabFourSNPs, shape = "ellipse")
dev.off()








#### 2.3.3. PC scores ####
## only PCs
# based on all of 188 metabolites




# based on 83 flavonoid metabolites




# based on 40 flavonoid metabolites, > 0.9 heritability









## PCs, with 4 SNPs
# based on all of 188 metabolites




# based on 83 flavonoid metabolites




# based on 40 flavonoid metabolites, > 0.9 heritability
gvPCScore2017Flavonoid40MetabCommon <- gvPCScore2017Flavonoid40Metab[varietyNamesCommon, ]
gvPCScore2017Flavonoid40MetabFourSNPs <- cbind(gvPCScore2017Flavonoid40MetabCommon, fourSNPsScoresCommon)

gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr06_18760995"] <- as.factor(gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr06_18760995"])
gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr06_47490224"] <- as.factor(gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr06_47490224"])
gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr10_42562665"] <- as.factor(gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr10_42562665"])
gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr17_16065902"] <- as.factor(gvPCScore2017Flavonoid40MetabFourSNPs[, "Chr17_16065902"])
See(gvPCScore2017Flavonoid40MetabFourSNPs)
table(is.na(gvPCScore2017Flavonoid40MetabFourSNPs))
# without blacklist
dagGvPCScore2017Flavonoid40MetabFourSNPs <- hc(gvPCScore2017Flavonoid40MetabFourSNPs[, -c(5, 6)], score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_PC_based_on_40_Flavonoid_metabolites_heritability_>0.9_without_blacklist.pdf"))
graphviz.plot(dagGvPCScore2017Flavonoid40MetabFourSNPs, shape = "ellipse")
dev.off()
# with blacklist
dagGvPCScore2017Flavonoid40MetabFourSNPs <- hc(gvPCScore2017Flavonoid40MetabFourSNPs, score = "bic-cg", blacklist = blFlavonoidMetabFourSNPs)
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_Flavonoid_40_metabolites_heritability_>_0.9_with_blacklist.pdf"))
graphviz.plot(dagGvPCScore2017Flavonoid40MetabFourSNPs, shape = "ellipse")
dev.off()

children(dagGvPCScore2017Flavonoid40MetabFourSNPs, "Chr06_18760995")





#### 2.3.4. Only 4 SNPs ####
See(fourSNPsScoresCommon)
fourSNPsScoresCommon <- as.data.frame(fourSNPsScoresCommon)

fourSNPsScoresCommon[, "Chr06_18760995"] <- as.factor(fourSNPsScoresCommon[, "Chr06_18760995"])
fourSNPsScoresCommon[, "Chr06_47490224"] <- as.factor(fourSNPsScoresCommon[, "Chr06_47490224"])
fourSNPsScoresCommon[, "Chr10_42562665"] <- as.factor(fourSNPsScoresCommon[, "Chr10_42562665"])
fourSNPsScoresCommon[, "Chr17_16065902"] <- as.factor(fourSNPsScoresCommon[, "Chr17_16065902"])

?hc

dagFourSNPsScoresCommon <- hc(fourSNPsScoresCommon, score = "bic")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_Four_SNPs_without_blacklist.pdf"))
graphviz.plot(dagFourSNPsScoresCommon, shape = "ellipse")
dev.off()





#### 2.3.5 Given SNPs ####
### 188 metab, given F3H
gvMetab2017TotalthreeSNPsScoresGivenF3H0Common <- cbind(gvMetab2017TotalF3H0Common, threeSNPSScoresGivenF3H0Common)
gvMetab2017TotalthreeSNPsScoresGivenF3H2Common <- cbind(gvMetab2017TotalF3H2Common, threeSNPSScoresGivenF3H2Common)
See(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common)
See(gvMetab2017TotalthreeSNPsScoresGivenF3H2Common)

gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "X00682"] <- as.numeric(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "X00682"])
gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "X500019"] <- as.numeric(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "X500019"])
gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "Chr06_47490224"])
gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "Chr10_42562665"])
gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common[, "Chr17_16065902"])


dagGvMetab2017TotalthreeSNPsScoresGivenF3H0Common <- hc(gvMetab2017TotalthreeSNPsScoresGivenF3H0Common, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_Four_SNPs_without_blacklist.pdf"))
graphviz.plot(dagGvMetab2017TotalthreeSNPsScoresGivenF3H0Common, shape = "ellipse")
dev.off()



### 83 flavonoid metabolites




### 40 flavonoid metabolites, heritability > 0.9, given F3H
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H0Common, threeSNPSScoresGivenF3H0Common)
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H2Common <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H2Common, threeSNPSScoresGivenF3H2Common)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H2Common)

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common[, "Chr06_47490224"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common[, "Chr10_42562665"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common[, "Chr17_16065902"])

dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common <- hc(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_40_metabolites_and_3SNPs_given_F3H_without_blacklist.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common, shape = "ellipse")
dev.off()


metabNamesChildrenChr06_47490224

metabFlavonoid <- read.csv(paste0("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv"))
metabNamesFlavonoid <- metabFlavonoid[, "Name"]
betaMat0Flavonoid <- betaMat0[metabNamesFlavonoid, ]
betaMat2Flavonoid <- betaMat2[metabNamesFlavonoid, ]

# metabNames20 <- rownames(betaMat0)[betaMat0[, 1] < -0.3 & abs(betaMat2[, 1]) <0.1]
metabNames00 <- rownames(betaMat0)[betaMat0[, 1] < 0 & betaMat0[, 1] > -0.5 & betaMat2[, 1] < -0.7]
metabNames20 <- rownames(betaMat0Flavonoid)[betaMat0Flavonoid[, 1] < -0.5 & abs(betaMat2Flavonoid[, 1]) <0.1]
metabNames00And20 <- c(metabNames20, metabNames00)
metabNames00And20 <- na.omit(metabNames00And20)
metabNames00And20




?Sys.time
start <- Sys.time()
bootStrGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common = boot.strength(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common,
                                        R = 1000,
                                        algorithm = "hc",
                                        algorithm.args = list(score="bic-cg"))
end <- Sys.time()

# ?boot.strength
head(bootStrGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common)
plot(bootStrGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common)
avg.diff = averaged.network(bootStrGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common)

pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_strength_40_metabolites_and_3SNPs_given_F3H.pdf"))
strength.plot(avg.diff, bootStrGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H0Common, shape = "ellipse")
dev.off()






### 40 flavonoid metabolites, heritability > 0.9, given F3H and Chr06_47490224
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H_0_Chr06_47490224_0_Common, twoSNPsScoresGivenF3H_0_Chr06_47490224_0_Common)
# gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H2Common <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9F3H2Common, threeSNPSScoresGivenF3H2Common)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common)
# See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9threeSNPSScoresF3H2Common)

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common[, "Chr10_42562665"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common[, "Chr17_16065902"])

dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common <- hc(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_40_metabolites_and_2SNPs_given_F3H_and_Chr06_47490224_without_blacklist.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common, shape = "ellipse")
dev.off()



metabNamesChildrenChr10_42562665GivenChr06_18760995AndChr06_47490224 <- children(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common, "Chr10_42562665")
metabNamesChildrenChr10_42562665GivenChr06_18760995AndChr06_47490224

metabNamesChildrenChr17_16065902GivenChr06_18760995AndChr06_47490224 <- children(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common, "Chr17_16065902")
metabNamesChildrenChr17_16065902GivenChr06_18760995AndChr06_47490224


dirTwoMarkers2017 <- paste0("midstream/2.32_GWAS_for_three_markers/2.32_2017_Metabolome_results/2.32_Coefficient_LD1_scaled_Chr06_18760995_Chr06_47490224/")
betaMat0_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,0).csv"), row.names = 1)
betaMat0_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(0,2).csv"), row.names = 1)
betaMat2_0LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,0).csv"), row.names = 1)
betaMat2_2LD1 <- read.csv(paste0(dirTwoMarkers2017, "2.32_(2,2).csv"), row.names = 1)

rownames(betaMat0_0LD1)[betaMat0_0LD1[, 1] > 0.9 & betaMat0_2LD1[, 1] < 0.5 & betaMat0_2LD1[, 1] > 0]
rownames(betaMat0_2LD1)[betaMat0_0LD1[, 1] < -1 & betaMat0_2LD1[, 1] > -0.7]
metabNamesF3HChr06_47490224ForChr10_42562665 <- rownames(betaMat0_0LD1)[betaMat0_0LD1[, 1] > 0.9 & betaMat0_2LD1[, 1] < 0.5 & betaMat0_2LD1[, 1] > 0]
metabNamesF3HChr06_47490224ForChr10_42562665 <- c(metabNamesF3HChr06_47490224ForChr10_42562665, rownames(betaMat0_2LD1)[betaMat0_0LD1[, 1] < -0.75 & betaMat0_2LD1[, 1] > -0.7])
metabNamesF3HChr06_47490224ForChr10_42562665 <- na.omit(metabNamesF3HChr06_47490224ForChr10_42562665)

table(metabNamesF3HChr06_47490224ForChr10_42562665 %in% metabNamesChildrenChr10_42562665GivenChr06_18760995AndChr06_47490224)
table(metabNamesChildrenChr10_42562665GivenChr06_18760995AndChr06_47490224 %in% metabNamesF3HChr06_47490224ForChr10_42562665)



rownames(betaMat0_0LD1)[betaMat0_0LD1[, 2] > 0.6]
rownames(betaMat0_0LD1)[betaMat0_0LD1[, 2] < -0.35 & betaMat0_2LD1[, 2] > -0.15]


bootStrGvMetab2017Total = boot.strength(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9twoSNPsScoresF3H_0_Chr06_47490224_0_Common[1:50], R = 200,
                                        algorithm = "hc",
                                        algorithm.args = list(score="bic-g"))







bootStrGvMetab2017Total = boot.strength(gvMetab2017Total[1:50], R = 200,
                                        algorithm = "hc",
                                        algorithm.args = list(score="bic-g"))

bootStrGvMetab2017Total





##### 2.4. Metabolites and SNP and SNP interaction #####
gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
See(gvMetab2017Total)

rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"



# #### 2.4.1. flavonoid, 83 metabolites ####
# gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
# rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"
#
# metabNamesAnnotationFlavonoid <- read.csv("data/extra/2017_Tottori_May_Metabolome_Flavonoid_Metab_Names_And_Annotation.csv")
# metabNamesFlavonoid <- metabNamesAnnotationFlavonoid[, "Name"]
#
# gvMetab2017TotalFlavonoid <- gvMetab2017Total[, metabNamesFlavonoid]
# See(gvMetab2017TotalFlavonoid)


#### 2.4.2. flavonoid, 40 metabolites, heritability > 0.9 ####
gvMetab2017Total <- read.csv("midstream/2.2_BSH/2.2_lmer_genotypic_values_Total_2017.csv", row.names = 1)
rownames(gvMetab2017Total)[rownames(gvMetab2017Total) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"

metab2017RawHeritabilityMoreThan0.9All <- read.csv("data/phenotype/2017_Tottori_May_Metabolome_No_Outlier_Related_To_Flavonoid_Pathway_>0.9_heritability.csv")
See(metab2017RawHeritabilityMoreThan0.9All, coln = 11)
metabNamesFlavonoidHeritabilityMoreThan0.9 <- colnames(metab2017RawHeritabilityMoreThan0.9All[, 11:50])

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9 <- gvMetab2017Total[, metabNamesFlavonoidHeritabilityMoreThan0.9]
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9)





genoMat <- as.matrix(gastonDataSmall)
rownames(genoMat)[rownames(genoMat) == "HOUJAKU_KUWAZU"] <- "HOUJAKU_KUWAZU"

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNA <- na.omit(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNA)
varietyNamesMetab <- rownames(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNA)
varietyNamesGeno <- rownames(genoMat)
varietyNamesCommon <- varietyNamesMetab[varietyNamesMetab %in% varietyNamesGeno]
See(varietyNamesCommon)

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon <- gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNA[varietyNamesCommon, ]
# gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon <- apply(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon, 2, as.numeric)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon)


genoMatCommon <- genoMat[varietyNamesCommon, ]
See(genoMatCommon)
fourSNPsScoresCommon <- genoMatCommon[, c("Chr06_18760995", "Chr06_47490224", "Chr10_42562665", "Chr17_16065902")]
See(fourSNPsScoresCommon)


SNPChr06_18760995 <- genoMatCommon[, "Chr06_18760995"]
SNPChr06_47490224 <- genoMatCommon[, "Chr06_47490224"]
SNPChr10_42562665 <- genoMatCommon[, "Chr10_42562665"]
SNPChr17_16065902 <- genoMatCommon[, "Chr17_16065902"]

Chr06_18760995_Chr06_47490224 <- interaction(SNPChr06_18760995, SNPChr06_47490224, drop = TRUE)
Chr06_18760995_Chr10_42562665 <- interaction(SNPChr06_18760995, SNPChr10_42562665, drop = TRUE)
Chr06_18760995_Chr17_16065902 <- interaction(SNPChr06_18760995, SNPChr17_16065902, drop = TRUE)

Chr06_47490224_Chr10_42562665 <- interaction(SNPChr06_47490224, SNPChr10_42562665, drop = TRUE)
Chr06_47490224_Chr17_16065902 <- interaction(SNPChr06_47490224, SNPChr17_16065902, drop = TRUE)
Chr10_42562665_Chr17_16065902 <- interaction(SNPChr10_42562665, SNPChr17_16065902, drop = TRUE)

Chr06_18760995_Chr06_47490224_Chr10_42562665 <- interaction(SNPChr06_18760995, SNPChr06_47490224, SNPChr10_42562665)


gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon, fourSNPsScoresCommon, Chr06_18760995_Chr06_47490224, Chr06_18760995_Chr10_42562665, Chr06_18760995_Chr17_16065902, Chr06_47490224_Chr10_42562665, Chr06_47490224_Chr17_16065902, Chr10_42562665_Chr17_16065902)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, 40:50], coln = 10, rown= 10)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction)
# gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction <- as.data.frame(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction)

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_18760995"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_18760995"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_47490224"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr10_42562665"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr17_16065902"])



dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction <-  hc(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_40_metabolites_2SNP_Interaction.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction, shape = "ellipse")
dev.off()



gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon, fourSNPsScoresCommon, Chr06_18760995_Chr06_47490224, Chr06_18760995_Chr10_42562665, Chr06_18760995_Chr17_16065902, Chr06_47490224_Chr10_42562665, Chr06_47490224_Chr17_16065902, Chr10_42562665_Chr17_16065902, Chr06_18760995_Chr06_47490224_Chr10_42562665)

gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_18760995"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_18760995"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_47490224"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr06_47490224"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr10_42562665"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr10_42562665"])
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr17_16065902"] <- as.factor(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, "Chr17_16065902"])

See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction[, 40:51], coln = 14)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction)

dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction <-  hc(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_40_metabolites_2&3SNP_Interaction.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction, shape = "ellipse")
dev.off()


# only 2SNP interactions
gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPInteraction <- cbind(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon, Chr06_18760995_Chr06_47490224, Chr06_18760995_Chr10_42562665, Chr06_18760995_Chr17_16065902, Chr06_47490224_Chr10_42562665, Chr06_47490224_Chr17_16065902, Chr10_42562665_Chr17_16065902)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPInteraction[, 40:46], coln = 10, rown= 10)
See(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPInteraction)

dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPInteraction <-  hc(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPInteraction, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_40_metabolites_only_2SNP_Interaction.pdf"))
graphviz.plot(dagGvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPInteraction, shape = "ellipse")
dev.off()



### PCs, 40 metabolites, heritability > 0.9
gvPCScore2017Flavonoid40Metab <- read.csv("midstream/2.26_BSH_for_PCA_with_high_heritability_flavonoid/2.26_lmer_genotypic_values_for_PC_Score_for_flavonoid_metab_>0.9_heritability_in_2017_pcaMethods_nPC=6.csv", row.names = 1)
See(gvPCScore2017Flavonoid40Metab)

gvPCScore2017Flavonoid40Metab <- gvPCScore2017Flavonoid40Metab[, 1:4]

rownames(gvPCScore2017Flavonoid40Metab)[rownames(gvPCScore2017Flavonoid40Metab) == "Houjaku Kuwazu"] <- "HOUJAKU_KUWAZU"

table(is.na(gvPCScore2017Flavonoid40Metab))



genoMat <- as.matrix(gastonDataSmall)
rownames(genoMat)[rownames(genoMat) == "HOUJAKU_KUWAZU"] <- "HOUJAKU_KUWAZU"

gvPCScore2017Flavonoid40Metab <- na.omit(gvPCScore2017Flavonoid40Metab)
See(gvPCScore2017Flavonoid40Metab)
varietyNamesMetab <- rownames(gvPCScore2017Flavonoid40Metab)
varietyNamesGeno <- rownames(genoMat)
varietyNamesCommon <- varietyNamesMetab[varietyNamesMetab %in% varietyNamesGeno]
See(varietyNamesCommon)

gvPCScore2017Flavonoid40MetabCommon <- gvPCScore2017Flavonoid40Metab[varietyNamesCommon, ]
# gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon <- apply(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9NoNACommon, 2, as.numeric)
See(gvPCScore2017Flavonoid40MetabCommon)


genoMatCommon <- genoMat[varietyNamesCommon, ]
fourSNPsScoresCommon <- genoMatCommon[, c("Chr06_18760995", "Chr06_47490224", "Chr10_42562665", "Chr17_16065902")]
See(fourSNPsScoresCommon)

SNPChr06_18760995 <- genoMatCommon[, "Chr06_18760995"]
SNPChr06_47490224 <- genoMatCommon[, "Chr06_47490224"]
SNPChr10_42562665 <- genoMatCommon[, "Chr10_42562665"]
SNPChr17_16065902 <- genoMatCommon[, "Chr17_16065902"]

Chr06_18760995_Chr06_47490224 <- interaction(SNPChr06_18760995, SNPChr06_47490224, drop = TRUE)
Chr06_18760995_Chr10_42562665 <- interaction(SNPChr06_18760995, SNPChr10_42562665, drop = TRUE)
Chr06_18760995_Chr17_16065902 <- interaction(SNPChr06_18760995, SNPChr17_16065902, drop = TRUE)

Chr06_47490224_Chr10_42562665 <- interaction(SNPChr06_47490224, SNPChr10_42562665, drop = TRUE)
Chr06_47490224_Chr17_16065902 <- interaction(SNPChr06_47490224, SNPChr17_16065902, drop = TRUE)
Chr10_42562665_Chr17_16065902 <- interaction(SNPChr10_42562665, SNPChr17_16065902, drop = TRUE)

Chr06_18760995_Chr06_47490224_Chr10_42562665 <- interaction(SNPChr06_18760995, SNPChr06_47490224, SNPChr10_42562665)


# blacklist
targetNodes <- c("Chr06_18760995",
                 "Chr06_47490224",
                 "Chr10_42562665",
                 "Chr17_16065902",
                 "Chr06_18760995_Chr06_47490224",
                 "Chr06_18760995_Chr10_42562665",
                 "Chr06_18760995_Chr17_16065902",
                 "Chr06_47490224_Chr10_42562665",
                 "Chr06_47490224_Chr17_16065902",
                 "Chr10_42562665_Chr17_16065902"
                 )
blacklist4SNPand4SNPInteractions <- expand.grid(from = targetNodes, to = targetNodes, stringsAsFactors = FALSE)
blacklist4SNPand4SNPInteractions <- blacklist4SNPand4SNPInteractions[blacklist4SNPand4SNPInteractions$from != blacklist4SNPand4SNPInteractions$to, ]
rownames(blacklist4SNPand4SNPInteractions) <- NULL
print(blacklist4SNPand4SNPInteractions)



targetNodes <- c("Chr06_18760995_Chr06_47490224",
                 "Chr06_18760995_Chr10_42562665",
                 "Chr06_18760995_Chr17_16065902",
                 "Chr06_47490224_Chr10_42562665",
                 "Chr06_47490224_Chr17_16065902",
                 "Chr10_42562665_Chr17_16065902"
                 )
blacklistOnly4SNPInteractions <- expand.grid(from = targetNodes, to = targetNodes, stringsAsFactors = FALSE)
blacklistOnly4SNPInteractions <- blacklistOnly4SNPInteractions[blacklistOnly4SNPInteractions$from != blacklistOnly4SNPInteractions$to, ]
rownames(blacklistOnly4SNPInteractions) <- NULL
print(blacklistOnly4SNPInteractions)



# 4SNPs and 2SNP-interactions
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction <- cbind(gvPCScore2017Flavonoid40MetabCommon, fourSNPsScoresCommon, Chr06_18760995_Chr06_47490224, Chr06_18760995_Chr10_42562665, Chr06_18760995_Chr17_16065902, Chr06_47490224_Chr10_42562665, Chr06_47490224_Chr17_16065902, Chr10_42562665_Chr17_16065902)
See(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, coln = 15, rown= 10)
See(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction)
# gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction <- as.data.frame(gvMetab2017TotalFlavonoidHeritabilityMoreThan0.9FourSNPsSNPInteraction)

gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_18760995"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_18760995"])
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_47490224"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_47490224"])
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr10_42562665"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr10_42562665"])
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr17_16065902"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr17_16065902"])



dagGvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction <-  hc(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_PCs_from_40_metabolites_2SNP_Interaction.pdf"))
graphviz.plot(dagGvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, shape = "ellipse")
dev.off()


# 4SNPs and 2&3 SNP-interactions
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction <- cbind(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, fourSNPsScoresCommon, Chr06_18760995_Chr06_47490224, Chr06_18760995_Chr10_42562665, Chr06_18760995_Chr17_16065902, Chr06_47490224_Chr10_42562665, Chr06_47490224_Chr17_16065902, Chr10_42562665_Chr17_16065902, Chr06_18760995_Chr06_47490224_Chr10_42562665)

gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_18760995"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_18760995"])
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_47490224"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr06_47490224"])
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr10_42562665"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr10_42562665"])
gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr17_16065902"] <- as.factor(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction[, "Chr17_16065902"])

See(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, coln = 14)
See(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction)

dagGvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction <-  hc(gvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_PCs_from_40_metabolites_2&3SNP_Interaction.pdf"))
graphviz.plot(dagGvPCScore2017Flavonoid40MetabCommonFourSNPsSNPInteraction, shape = "ellipse")
dev.off()



# only 2SNP-interactions of 4SNP
gvPCScore2017Flavonoid40MetabCommonFourSNPInteraction <- cbind(gvPCScore2017Flavonoid40MetabCommon, Chr06_18760995_Chr06_47490224, Chr06_18760995_Chr10_42562665, Chr06_18760995_Chr17_16065902, Chr06_47490224_Chr10_42562665, Chr06_47490224_Chr17_16065902, Chr10_42562665_Chr17_16065902)
See(gvPCScore2017Flavonoid40MetabCommonFourSNPInteraction, coln = 15, rown= 10)
See(gvPCScore2017Flavonoid40MetabCommonFourSNPInteraction)

dagGvPCScore2017Flavonoid40MetabCommonFourSNPInteraction <-  hc(gvPCScore2017Flavonoid40MetabCommonFourSNPInteraction, score = "bic-cg")
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_PCs_from_40_metabolites_only_2SNP_Interaction.pdf"))
graphviz.plot(dagGvPCScore2017Flavonoid40MetabCommonFourSNPInteraction, shape = "ellipse")
dev.off()


dagGvPCScore2017Flavonoid40MetabCommonFourSNPInteraction <-  hc(gvPCScore2017Flavonoid40MetabCommonFourSNPInteraction, score = "bic-cg", bl = blacklistOnly4SNPInteractions)
pdf(paste0(dirMidSTAMBayesianNetwork, scriptID, "_PCs_from_40_metabolites_only_2SNP_Interaction_with_bl.pdf"))
graphviz.plot(dagGvPCScore2017Flavonoid40MetabCommonFourSNPInteraction, shape = "ellipse")
dev.off()











