# Characterize methylation profile of each individual - get percent methylation per base and depth of coverage
# bwameth aligned reads

#load pacakges
library(methylKit)
library(ggplot2)
library(viridis)
library(dplyr)
library(car)
library(tidyverse)

#import metadata
setwd("/minDepth10")

# Set file.list of samples
file.list.10=list("CG_22_GOS_001_CpG.methylKit", "CG_22_GOS_002_CpG.methylKit", "CG_GOS_35_04_CpG.methylKit", "CG_GOS_2210_01_CpG.methylKit", "FD_GOS_27_CpG.methylKit", "FD_GOS_28_CpG.methylKit", "FD_GOS_33_CpG.methylKit", "CG_22_ROB_001_CpG.methylKit", "CG_22_ROB_002_CpG.methylKit", "CG_ROB_003_CpG.methylKit", "CG_ROB_31_01_CpG.methylKit", "FD_ROB_98_CpG.methylKit", "FD_ROB_102_CpG.methylKit", "FD_ROB_106_CpG.methylKit", "CG_22_SAY_001_CpG.methylKit", "CG_22_SAY_002_CpG.methylKit", "CG_SAY_4_01_CpG.methylKit", "CG_SAY_2343_003_CpG.methylKit", "CG_22_WK_001_CpG.methylKit", "CG_22_WK_002_CpG.methylKit", "CG_22_WK_003_CpG.methylKit", "CG_22_WK_004_CpG.methylKit", "FD_22_WK_002_CpG.methylKit", "FD_22_WK_005_CpG.methylKit", "FD_22_WK_011_CpG.methylKit", "FD_22_WK_012_CpG.methylKit", "CG_22_WT_001_CpG.methylKit", "CG_22_WT_002_CpG.methylKit", "CG_22_WT_003_CpG.methylKit", "CG_22_WT_004_CpG.methylKit", "FD_22_WT_004_CpG.methylKit", "FD_22_WT_007_CpG.methylKit", "FD_22_WT_016_CpG.methylKit", "FD_22_WT_017_CpG.methylKit")

# read the files to a methylRawList object: myobj, min coverage of 10
myobj.10=methRead(file.list.10,
                  sample.id=list("CG_GOS_001", "CG_GOS_002", "CG_GOS_003", "CG_GOS_004", "FD_GOS_027", "FD_GOS_028", "FD_GOS_033", "CG_ROB_001", "CG_ROB_002", "CG_ROB_003", "CG_ROB_004", "FD_ROB_098", "FD_ROB_102", "FD_ROB_106", "CG_SAY_001", "CG_SAY_002", "CG_SAY_003", "CG_SAY_004", "CG_WK_001", "CG_WK_002", "CG_WK_003", "CG_WK_004", "FD_WK_002", "FD_WK_005", "FD_WK_011", "FD_WK_012","CG_WT_001", "CG_WT_002", "CG_WT_003", "CG_WT_004", "FD_WT_004", "FD_WT_007", "FD_WT_016", "FD_WT_017"),
                  pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                  assembly="StickleGeneAnnotations",
                  treatment=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5),
                  context="CpG",
                  mincov = 10
)

# Methylation stats and plotting
getMethylationStats(myobj.10[[1]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[1]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[2]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[2]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[3]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[3]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[4]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[4]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[5]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[5]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[6]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[6]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[7]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[7]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[8]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[8]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[9]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[9]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[10]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[10]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[11]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[11]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[12]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[12]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[13]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[13]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[14]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[14]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[15]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[15]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[16]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[16]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[17]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[17]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[18]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[18]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[19]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[19]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[20]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[20]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[21]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[21]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[22]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[22]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[23]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[23]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[24]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[24]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[25]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[25]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[26]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[26]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[27]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[27]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[28]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[28]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[29]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[29]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[30]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[30]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[31]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[31]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[32]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[32]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[33]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[33]],plot=TRUE,both.strands=FALSE)

getMethylationStats(myobj.10[[34]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[34]],plot=TRUE,both.strands=FALSE)

#plot based on read coverage per base
getCoverageStats(myobj.10[[1]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[1]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[2]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[3]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[3]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[4]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[4]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[5]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[5]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[6]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[6]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[7]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[7]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[8]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[8]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[9]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[9]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[10]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[10]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[11]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[11]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[12]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[12]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[13]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[13]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[14]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[14]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[15]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[15]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[16]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[16]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[17]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[17]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[18]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[18]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[19]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[19]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[20]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[20]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[21]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[21]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[22]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[22]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[23]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[23]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[24]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[24]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[25]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[25]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[26]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[26]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[27]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[27]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[28]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[28]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[29]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[29]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[30]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[30]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[31]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[31]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[32]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[32]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[33]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[33]],plot=FALSE,both.strands=FALSE)

getCoverageStats(myobj.10[[34]],plot=TRUE,both.strands=FALSE)
getCoverageStats(myobj.10[[34]],plot=FALSE,both.strands=FALSE)

# Filter reads with <10x coverage and in the 99.9th percentile (to account for PCR bias)
filtered.myobjDepth10=filterByCoverage(myobj.10,lo.count=10,lo.perc=NULL,
                                       hi.count=NULL,hi.perc=99.9)

# normalize read coverage between samples to avoid bias introduced by systematically more sequenced sameples
normalized_myobjDepth10=normalizeCoverage(filtered.myobjDepth10, method="median")

#merge files to include only base pairs covered in every sample. destrand=TRUE merges both strands of a CpG dinucleotide to provide greater depth of coverage
meth.10=methylKit::unite(normalized_myobjDepth10, destrand=TRUE)
head(meth.10)
#creates dataframe with 252 rows - so 252 CpG sites were covered in all individuals. This parameter can be relaxed to allow for sites with missing data in some individuals to be included.

# Separate each sample into individual methylRaw objects - contains info about mehtylation and chr location
CG_GOS_001=myobj.10[[1]] 
CG_GOS_002=myobj.10[[2]] 
CG_GOS_003=myobj.10[[3]] 
CG_GOS_004=myobj.10[[4]] 
FD_GOS_027=myobj.10[[5]] 
FD_GOS_028=myobj.10[[6]] 
FD_GOS_033=myobj.10[[7]] 
CG_ROB_001=myobj.10[[8]] 
CG_ROB_002=myobj.10[[9]] 
CG_ROB_003=myobj.10[[10]] 
CG_ROB_004=myobj.10[[11]] 
FD_ROB_098=myobj.10[[12]] 
FD_ROB_102=myobj.10[[13]] 
FD_ROB_106=myobj.10[[14]] 
CG_SAY_001=myobj.10[[15]] 
CG_SAY_002=myobj.10[[16]] 
CG_SAY_003=myobj.10[[17]] 
CG_SAY_004=myobj.10[[18]] 
CG_WK_001=myobj.10[[19]] 
CG_WK_002=myobj.10[[20]] 
CG_WK_003=myobj.10[[21]] 
CG_WK_004=myobj.10[[22]] 
FD_WK_002=myobj.10[[23]] 
FD_WK_005=myobj.10[[24]] 
FD_WK_011=myobj.10[[25]] 
FD_WK_012=myobj.10[[26]] 
CG_WT_001=myobj.10[[27]] 
CG_WT_002=myobj.10[[28]] 
CG_WT_003=myobj.10[[29]] 
CG_WT_004=myobj.10[[30]] 
FD_WT_004=myobj.10[[31]] 
FD_WT_007=myobj.10[[32]] 
FD_WT_016=myobj.10[[33]] 
FD_WT_017=myobj.10[[34]] 



CG_GOS_001 <- getData(CG_GOS_001)
CG_GOS_002 <- getData(CG_GOS_002)
CG_GOS_003 <- getData(CG_GOS_003)
CG_GOS_004 <- getData(CG_GOS_004)
FD_GOS_027 <- getData(FD_GOS_027)
FD_GOS_028 <- getData(FD_GOS_028)
FD_GOS_033 <- getData(FD_GOS_033)
CG_ROB_001 <- getData(CG_ROB_001)
CG_ROB_002 <- getData(CG_ROB_002)
CG_ROB_003 <- getData(CG_ROB_003)
CG_ROB_004 <- getData(CG_ROB_004)
FD_ROB_098 <- getData(FD_ROB_098)
FD_ROB_102 <- getData(FD_ROB_102)
FD_ROB_106 <- getData(FD_ROB_106)
CG_SAY_001 <- getData(CG_SAY_001)
CG_SAY_002 <- getData(CG_SAY_002)
CG_SAY_003 <- getData(CG_SAY_003)
CG_SAY_004 <- getData(CG_SAY_004)
CG_WK_001 <- getData(CG_WK_001)
CG_WK_002 <- getData(CG_WK_002)
CG_WK_003 <- getData(CG_WK_003)
CG_WK_004 <- getData(CG_WK_004)
FD_WK_002 <- getData(FD_WK_002)
FD_WK_005 <- getData(FD_WK_005)
FD_WK_011 <- getData(FD_WK_011)
FD_WK_012 <- getData(FD_WK_012)
CG_WT_001 <- getData(CG_WT_001)
CG_WT_002 <- getData(CG_WT_002)
CG_WT_003 <- getData(CG_WT_003)
CG_WT_004 <- getData(CG_WT_004)
FD_WT_004 <- getData(FD_WT_004)
FD_WT_007 <- getData(FD_WT_007)
FD_WT_016 <- getData(FD_WT_016)
FD_WT_017 <- getData(FD_WT_017)

CG_GOS_001$SampleID <- "CG_GOS_001"
CG_GOS_002$SampleID <- "CG_GOS_002"
CG_GOS_003$SampleID <- "CG_GOS_003"
CG_GOS_004$SampleID <- "CG_GOS_004"
FD_GOS_027$SampleID <- "FD_GOS_027"
FD_GOS_028$SampleID <- "FD_GOS_028"
FD_GOS_033$SampleID <- "FD_GOS_033"
CG_ROB_001$SampleID <- "CG_ROB_001"
CG_ROB_002$SampleID <- "CG_ROB_002"
CG_ROB_003$SampleID <- "CG_ROB_003"
CG_ROB_004$SampleID <- "CG_ROB_004"
FD_ROB_098$SampleID <- "FD_ROB_098"
FD_ROB_102$SampleID <- "FD_ROB_102"
FD_ROB_106$SampleID <- "FD_ROB_106"
CG_SAY_001$SampleID <- "CG_SAY_001"
CG_SAY_002$SampleID <- "CG_SAY_002"
CG_SAY_003$SampleID <- "CG_SAY_003"
CG_SAY_004$SampleID <- "CG_SAY_004"
CG_WK_001$SampleID <- "CG_WK_001"
CG_WK_002$SampleID <- "CG_WK_002"
CG_WK_003$SampleID <- "CG_WK_003"
CG_WK_004$SampleID <- "CG_WK_004"
FD_WK_002$SampleID <- "FD_WK_002"
FD_WK_005$SampleID <- "FD_WK_005"
FD_WK_011$SampleID <- "FD_WK_011"
FD_WK_012$SampleID <- "FD_WK_012"
CG_WT_001$SampleID <- "CG_WT_001"
CG_WT_002$SampleID <- "CG_WT_002"
CG_WT_003$SampleID <- "CG_WT_003"
CG_WT_004$SampleID <- "CG_WT_004"
FD_WT_004$SampleID <- "FD_WT_004"
FD_WT_007$SampleID <- "FD_WT_007"
FD_WT_016$SampleID <- "FD_WT_016"
FD_WT_017$SampleID <- "FD_WT_017"

CG_GOS_001$Population <- "Gosling"
CG_GOS_002$Population <- "Gosling"
CG_GOS_003$Population <- "Gosling"
CG_GOS_004$Population <- "Gosling"
FD_GOS_027$Population <- "Gosling"
FD_GOS_028$Population <- "Gosling"
FD_GOS_033$Population <- "Gosling"
CG_ROB_001$Population <- "Roberts"
CG_ROB_002$Population <- "Roberts"
CG_ROB_003$Population <- "Roberts"
CG_ROB_004$Population <- "Roberts"
FD_ROB_098$Population <- "Roberts"
FD_ROB_102$Population <- "Roberts"
FD_ROB_106$Population <- "Roberts"
CG_SAY_001$Population <- "Sayward"
CG_SAY_002$Population <- "Sayward"
CG_SAY_003$Population <- "Sayward"
CG_SAY_004$Population <- "Sayward"
CG_WK_001$Population <- "Wik"
CG_WK_002$Population <- "Wik"
CG_WK_003$Population <- "Wik"
CG_WK_004$Population <- "Wik"
FD_WK_002$Population <- "Wik"
FD_WK_005$Population <- "Wik"
FD_WK_011$Population <- "Wik"
FD_WK_012$Population <- "Wik"
CG_WT_001$Population <- "Watson"
CG_WT_002$Population <- "Watson"
CG_WT_003$Population <- "Watson"
CG_WT_004$Population <- "Watson"
FD_WT_004$Population <- "Watson"
FD_WT_007$Population <- "Watson"
FD_WT_016$Population <- "Watson"
FD_WT_017$Population <- "Watson"

CG_GOS_001$Environment <- "CommonGarden"
CG_GOS_002$Environment <- "CommonGarden"
CG_GOS_003$Environment <- "CommonGarden"
CG_GOS_004$Environment <- "CommonGarden"
FD_GOS_027$Environment <- "Field"
FD_GOS_028$Environment <- "Field"
FD_GOS_033$Environment <- "Field"
CG_ROB_001$Environment <- "CommonGarden"
CG_ROB_002$Environment <- "CommonGarden"
CG_ROB_003$Environment <- "CommonGarden"
CG_ROB_004$Environment <- "CommonGarden"
FD_ROB_098$Environment <- "Field"
FD_ROB_102$Environment <- "Field"
FD_ROB_106$Environment <- "Field"
CG_SAY_001$Environment <- "CommonGarden"
CG_SAY_002$Environment <- "CommonGarden"
CG_SAY_003$Environment <- "CommonGarden"
CG_SAY_004$Environment <- "CommonGarden"
CG_WK_001$Environment <- "CommonGarden"
CG_WK_002$Environment <- "CommonGarden"
CG_WK_003$Environment <- "CommonGarden"
CG_WK_004$Environment <- "CommonGarden"
FD_WK_002$Environment <- "Field"
FD_WK_005$Environment <- "Field"
FD_WK_011$Environment <- "Field"
FD_WK_012$Environment <- "Field"
CG_WT_001$Environment <- "CommonGarden"
CG_WT_002$Environment <- "CommonGarden"
CG_WT_003$Environment <- "CommonGarden"
CG_WT_004$Environment <- "CommonGarden"
FD_WT_004$Environment <- "Field"
FD_WT_007$Environment <- "Field"
FD_WT_016$Environment <- "Field"
FD_WT_017$Environment <- "Field"

CG_GOS_001$Batch <- "1"
CG_GOS_002$Batch <- "1"
CG_GOS_003$Batch <- "2"
CG_GOS_004$Batch <- "2"
FD_GOS_027$Batch <- "2"
FD_GOS_028$Batch <- "2"
FD_GOS_033$Batch <- "2"
CG_ROB_001$Batch <- "1"
CG_ROB_002$Batch <- "1"
CG_ROB_003$Batch <- "2"
CG_ROB_004$Batch <- "2"
FD_ROB_098$Batch <- "2"
FD_ROB_102$Batch <- "2"
FD_ROB_106$Batch <- "2"
CG_SAY_001$Batch <- "1"
CG_SAY_002$Batch <- "1"
CG_SAY_003$Batch <- "2"
CG_SAY_004$Batch <- "2"
CG_WK_001$Batch <- "1"
CG_WK_002$Batch <- "1"
CG_WK_003$Batch <- "1"
CG_WK_004$Batch <- "1"
FD_WK_002$Batch <- "1"
FD_WK_005$Batch <- "1"
FD_WK_011$Batch <- "1"
FD_WK_012$Batch <- "1"
CG_WT_001$Batch <- "1"
CG_WT_002$Batch <- "1"
CG_WT_003$Batch <- "1"
CG_WT_004$Batch <- "1"
FD_WT_004$Batch <- "1"
FD_WT_007$Batch <- "1"
FD_WT_016$Batch <- "1"
FD_WT_017$Batch <- "1"

CG_GOS_001$PercMeth <- (CG_GOS_001$numCs/CG_GOS_001$coverage)*100
CG_GOS_002$PercMeth <- (CG_GOS_002$numCs/CG_GOS_002$coverage)*100
CG_GOS_003$PercMeth <- (CG_GOS_003$numCs/CG_GOS_003$coverage)*100
CG_GOS_004$PercMeth <- (CG_GOS_004$numCs/CG_GOS_004$coverage)*100
FD_GOS_027$PercMeth <- (FD_GOS_027$numCs/FD_GOS_027$coverage)*100
FD_GOS_028$PercMeth <- (FD_GOS_028$numCs/FD_GOS_028$coverage)*100
FD_GOS_033$PercMeth <- (FD_GOS_033$numCs/FD_GOS_033$coverage)*100
CG_ROB_001$PercMeth <- (CG_ROB_001$numCs/CG_ROB_001$coverage)*100
CG_ROB_002$PercMeth <- (CG_ROB_002$numCs/CG_ROB_002$coverage)*100
CG_ROB_003$PercMeth <- (CG_ROB_003$numCs/CG_ROB_003$coverage)*100
CG_ROB_004$PercMeth <- (CG_ROB_004$numCs/CG_ROB_004$coverage)*100
FD_ROB_098$PercMeth <- (FD_ROB_098$numCs/FD_ROB_098$coverage)*100
FD_ROB_102$PercMeth <- (FD_ROB_102$numCs/FD_ROB_102$coverage)*100
FD_ROB_106$PercMeth <- (FD_ROB_106$numCs/FD_ROB_106$coverage)*100
CG_SAY_001$PercMeth <- (CG_SAY_001$numCs/CG_SAY_001$coverage)*100
CG_SAY_002$PercMeth <- (CG_SAY_002$numCs/CG_SAY_002$coverage)*100
CG_SAY_003$PercMeth <- (CG_SAY_003$numCs/CG_SAY_003$coverage)*100
CG_SAY_004$PercMeth <- (CG_SAY_004$numCs/CG_SAY_004$coverage)*100
CG_WK_001$PercMeth <- (CG_WK_001$numCs/CG_WK_001$coverage)*100
CG_WK_002$PercMeth <- (CG_WK_002$numCs/CG_WK_002$coverage)*100
CG_WK_003$PercMeth <- (CG_WK_003$numCs/CG_WK_003$coverage)*100
CG_WK_004$PercMeth <- (CG_WK_004$numCs/CG_WK_004$coverage)*100
FD_WK_002$PercMeth <- (FD_WK_002$numCs/FD_WK_002$coverage)*100
FD_WK_005$PercMeth <- (FD_WK_005$numCs/FD_WK_005$coverage)*100
FD_WK_011$PercMeth <- (FD_WK_011$numCs/FD_WK_011$coverage)*100
FD_WK_012$PercMeth <- (FD_WK_012$numCs/FD_WK_012$coverage)*100
CG_WT_001$PercMeth <- (CG_WT_001$numCs/CG_WT_001$coverage)*100
CG_WT_002$PercMeth <- (CG_WT_002$numCs/CG_WT_002$coverage)*100
CG_WT_003$PercMeth <- (CG_WT_003$numCs/CG_WT_003$coverage)*100
CG_WT_004$PercMeth <- (CG_WT_004$numCs/CG_WT_004$coverage)*100
FD_WT_004$PercMeth <- (FD_WT_004$numCs/FD_WT_004$coverage)*100
FD_WT_007$PercMeth <- (FD_WT_007$numCs/FD_WT_007$coverage)*100
FD_WT_016$PercMeth <- (FD_WT_016$numCs/FD_WT_016$coverage)*100
FD_WT_017$PercMeth <- (FD_WT_017$numCs/FD_WT_017$coverage)*100

AllSamples <- rbind(CG_GOS_001, CG_GOS_002)
AllSamples <- rbind(AllSamples, CG_GOS_003)
AllSamples <- rbind(AllSamples, CG_GOS_004)
AllSamples <- rbind(AllSamples, FD_GOS_027)
AllSamples <- rbind(AllSamples, FD_GOS_028)
AllSamples <- rbind(AllSamples, FD_GOS_033)
AllSamples <- rbind(AllSamples, CG_ROB_001)
AllSamples <- rbind(AllSamples, CG_ROB_002)
AllSamples <- rbind(AllSamples, CG_ROB_003)
AllSamples <- rbind(AllSamples, CG_ROB_004)
AllSamples <- rbind(AllSamples, FD_ROB_098)
AllSamples <- rbind(AllSamples, FD_ROB_102)
AllSamples <- rbind(AllSamples, FD_ROB_106)
AllSamples <- rbind(AllSamples, CG_SAY_001)
AllSamples <- rbind(AllSamples, CG_SAY_002)
AllSamples <- rbind(AllSamples, CG_SAY_003)
AllSamples <- rbind(AllSamples, CG_SAY_004)
AllSamples <- rbind(AllSamples, CG_WK_001)
AllSamples <- rbind(AllSamples, CG_WK_002)
AllSamples <- rbind(AllSamples, CG_WK_003)
AllSamples <- rbind(AllSamples, CG_WK_004)
AllSamples <- rbind(AllSamples, FD_WK_002)
AllSamples <- rbind(AllSamples, FD_WK_005)
AllSamples <- rbind(AllSamples, FD_WK_011)
AllSamples <- rbind(AllSamples, FD_WK_012)
AllSamples <- rbind(AllSamples, CG_WT_001)
AllSamples <- rbind(AllSamples, CG_WT_002)
AllSamples <- rbind(AllSamples, CG_WT_003)
AllSamples <- rbind(AllSamples, CG_WT_004)
AllSamples <- rbind(AllSamples, FD_WT_004)
AllSamples <- rbind(AllSamples, FD_WT_007)
AllSamples <- rbind(AllSamples, FD_WT_016)
AllSamples <- rbind(AllSamples, FD_WT_017)


BatchPlot <- ggplot(AllSamples)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  #theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("Percent Methylation per Base") +
  facet_wrap(~Batch)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)) +
  labs(tag = "A")
BatchPlot

PopPlot <- ggplot(AllSamples)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  #theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("Percent Methylation per Base") +
  facet_wrap(~Population)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)) +
  labs(tag = "B")
PopPlot

EnvPlot <- ggplot(AllSamples)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  #theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("Percent Methylation per Base") +
  facet_wrap(~Environment)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)) +
  labs(tag = "B")
EnvPlot

SamPlot <- ggplot(AllSamples)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  #theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("Percent Methylation per Base") +
  facet_wrap(~SampleID)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)) +
  labs(tag = "B")
SamPlot


BatchPlot
