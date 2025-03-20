#Methylation methods comparison
# Bismark aligned reads
#MethylKit analysis

#load pacakges
library(methylKit)
library(ggplot2)
library(viridis)
library(dplyr)
library(car)
library(tidyverse)

setwd("/PATH/cov")

# Set file.list of samples
file.list.10=list("CG_22_GOS_001_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_GOS_002_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_GOS_35_04_R1_001_val_1_bismark_bt2_pe.bismark.cov", "CG_GOS_2210_01_R1_001_val_1_bismark_bt2_pe.bismark.cov", "FD_GOS_27_R1_001_val_1_bismark_bt2_pe.bismark.cov", "FD_GOS_28_R1_001_val_1_bismark_bt2_pe.bismark.cov", "FD_GOS_33_R1_001_val_1_bismark_bt2_pe.bismark.cov",  "CG_22_ROB_001_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_ROB_002_R1_val_1_bismark_bt2_pe.bismark.cov",  "CG_ROB_003_R1_001_val_1_bismark_bt2_pe.bismark.cov", "CG_ROB_31_01_R1_001_val_1_bismark_bt2_pe.bismark.cov", "FD_ROB_98_R1_001_val_1_bismark_bt2_pe.bismark.cov", "FD_ROB_102_R1_001_val_1_bismark_bt2_pe.bismark.cov", "FD_ROB_106_R1_001_val_1_bismark_bt2_pe.bismark.cov", "CG_22_SAY_001_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_SAY_002_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_SAY_4_01_R1_001_val_1_bismark_bt2_pe.bismark.cov", "CG_SAY_2343_003_R1_001_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WK_001_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WK_002_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WK_003_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WK_004_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WK_002_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WK_005_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WK_011_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WK_012_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WT_001_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WT_002_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WT_003_R1_val_1_bismark_bt2_pe.bismark.cov", "CG_22_WT_004_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WT_004_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WT_007_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WT_016_R1_val_1_bismark_bt2_pe.bismark.cov", "FD_22_WT_017_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
myobj.10 <- methRead(file.list.10,
                     sample.id=list("CG_GOS_001", "CG_GOS_002", "CG_GOS_003", "CG_GOS_004", "FD_GOS_027", "FD_GOS_028", "FD_GOS_033", "CG_ROB_001", "CG_ROB_002", "CG_ROB_003", "CG_ROB_004", "FD_ROB_098", "FD_ROB_102", "FD_ROB_106", "CG_SAY_001", "CG_SAY_002", "CG_SAY_003", "CG_SAY_004", "CG_WK_001", "CG_WK_002", "CG_WK_003","CG_WK_004", "FD_WK_002", "FD_WK_005", "FD_WK_011", "FD_WK_012", "CG_WT_001", "CG_WT_002", "CG_WT_003", "CG_WT_004", "FD_WT_004", "FD_WT_007", "FD_WT_016", "FD_WT_017"),
                     assembly="StickleGeneAnnotations",
                     treatment = c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5),
                     pipeline = "bismarkCoverage",
                     mincov = 10)

# Methylation stats and plotting
getMethylationStats(myobj.10[[1]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[1]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[1]])

getMethylationStats(myobj.10[[2]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[2]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[2]])

getMethylationStats(myobj.10[[3]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[3]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[3]])

getMethylationStats(myobj.10[[4]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[4]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[4]])

getMethylationStats(myobj.10[[5]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[5]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[5]])

getMethylationStats(myobj.10[[6]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[6]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[6]])

getMethylationStats(myobj.10[[7]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[7]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[7]])

getMethylationStats(myobj.10[[8]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[8]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[8]])

getMethylationStats(myobj.10[[9]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[9]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[9]])

getMethylationStats(myobj.10[[10]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[10]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[10]])

getMethylationStats(myobj.10[[11]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[11]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[11]])

getMethylationStats(myobj.10[[12]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[12]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[12]])

getMethylationStats(myobj.10[[13]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[13]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[13]])

getMethylationStats(myobj.10[[14]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[14]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[14]])

getMethylationStats(myobj.10[[15]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[15]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[15]])

getMethylationStats(myobj.10[[16]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[16]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[16]])

getMethylationStats(myobj.10[[17]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[17]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[17]])

getMethylationStats(myobj.10[[18]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[18]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[18]])

getMethylationStats(myobj.10[[19]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[19]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[19]])

getMethylationStats(myobj.10[[20]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[20]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[20]])

getMethylationStats(myobj.10[[21]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[21]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[21]])

getMethylationStats(myobj.10[[22]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[22]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[22]])

getMethylationStats(myobj.10[[23]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[23]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[23]])

getMethylationStats(myobj.10[[24]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[24]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[24]])

getMethylationStats(myobj.10[[25]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[25]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[25]])

getMethylationStats(myobj.10[[26]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[26]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[26]])

getMethylationStats(myobj.10[[27]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[27]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[27]])

getMethylationStats(myobj.10[[28]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[28]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[28]])

getMethylationStats(myobj.10[[29]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[29]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[29]])

getMethylationStats(myobj.10[[30]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[30]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[30]])

getMethylationStats(myobj.10[[31]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[31]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[31]])

getMethylationStats(myobj.10[[32]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[32]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[32]])

getMethylationStats(myobj.10[[33]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[33]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[33]])

getMethylationStats(myobj.10[[34]],plot=FALSE,both.strands=FALSE)
getMethylationStats(myobj.10[[34]],plot=TRUE,both.strands=FALSE)
count(myobj.10[[34]])

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
