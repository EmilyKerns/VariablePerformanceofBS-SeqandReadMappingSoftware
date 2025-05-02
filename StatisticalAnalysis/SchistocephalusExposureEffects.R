library(methylKit)
library(ggplot2)
library(viridis)
library(dplyr)
library(magrittr)
library(car)
library(tidyverse)

setwd('PATH/RRBS/methylKit/MaxVarFrac8MinDepth5')

info = read.delim("R:/Genohub_seq/methylation_methods_comparison/Metadata/Metadata.txt")
labels = read.delim("R:/Genohub_seq/methylation_methods_comparison/Metadata/Label.txt")
info <- info %>% 
  mutate("Label" = (names_from = labels))

file.list.10 <- list(
  "CG_22_GOS_001_CpG.methylKit", "CG_22_GOS_002_CpG.methylKit", "CG_GOS_35_04_CpG.methylKit", "CG_GOS_2210_01_CpG.methylKit", "FD_GOS_27_CpG.methylKit", "FD_GOS_28_CpG.methylKit", "FD_GOS_33_CpG.methylKit", "CG_22_ROB_001_CpG.methylKit", "CG_22_ROB_002_CpG.methylKit", "CG_ROB_003_CpG.methylKit", "CG_ROB_31_01_CpG.methylKit", "FD_ROB_98_CpG.methylKit", "FD_ROB_102_CpG.methylKit", "FD_ROB_106_CpG.methylKit", "CG_22_SAY_001_CpG.methylKit", "CG_22_SAY_002_CpG.methylKit", "CG_SAY_4_01_CpG.methylKit", "CG_SAY_2343_003_CpG.methylKit", "CG_22_WK_001_CpG.methylKit", "CG_22_WK_002_CpG.methylKit", "CG_22_WK_003_CpG.methylKit", "CG_22_WK_004_CpG.methylKit", "FD_22_WK_002_CpG.methylKit", "FD_22_WK_005_CpG.methylKit", "FD_22_WK_011_CpG.methylKit", "FD_22_WK_012_CpG.methylKit", "CG_22_WT_001_CpG.methylKit", "CG_22_WT_002_CpG.methylKit", "CG_22_WT_003_CpG.methylKit", "CG_22_WT_004_CpG.methylKit", "FD_22_WT_004_CpG.methylKit", "FD_22_WT_007_CpG.methylKit", "FD_22_WT_016_CpG.methylKit", "FD_22_WT_017_CpG.methylKit"
)

myobj.10=methRead(file.list.10,
                  sample.id=list("CG_GOS_001", "CG_GOS_002", "CG_GOS_003", "CG_GOS_004", "FD_GOS_027", "FD_GOS_028", "FD_GOS_033", "CG_ROB_001", "CG_ROB_002", "CG_ROB_003", "CG_ROB_004", "FD_ROB_098", "FD_ROB_102", "FD_ROB_106", "CG_SAY_001", "CG_SAY_002", "CG_SAY_003", "CG_SAY_004", "CG_WK_001", "CG_WK_002", "CG_WK_003", "CG_WK_004", "FD_WK_002", "FD_WK_005", "FD_WK_011", "FD_WK_012","CG_WT_001", "CG_WT_002", "CG_WT_003", "CG_WT_004", "FD_WT_004", "FD_WT_007", "FD_WT_016", "FD_WT_017"),
                  pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                  assembly="StickleGeneAnnotations",
                  treatment=c(1,1,2,2,2,2,2,1,1,2,2,2,2,2,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                  context="CpG",
                  mincov = 10
)

# Filter reads with <10x coverage and in the 99.9th percentile (to account for PCR bias)
filtered.myobjDepth10=filterByCoverage(myobj.10,lo.count=10,lo.perc=NULL,
                                       hi.count=NULL,hi.perc=99.9)

# normalize read coverage between samples to avoid bias introduced by systematically more sequenced sameples
normalized_myobjDepth10=normalizeCoverage(filtered.myobjDepth10, method="median")

#merge files to include only base pairs covered in every sample. destrand=TRUE merges both strands of a CpG dinucleotide to provide greater depth of coverage
meth.10=methylKit::unite(normalized_myobjDepth10, destrand=TRUE)
head(meth.10)

# remove CpGs that aren't variable between samples to reduce dimensionality of the data
# get percent methylation matrix
pm=percMethylation(meth.10)
nrow(pm) #number of CpG sites 

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation to determine a suitable cutoff
hist(sds, breaks = 100)

# keep only CpG with standard deviations larger than 2%
meth_10 <- meth.10[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth_10) 

#sample correlation
getCorrelation(meth_10,plot=FALSE)

#clustering samples
clusterSamples(meth_10, dist="correlation", method="ward.D", plot=TRUE)

hc.10 = clusterSamples(meth_10, dist="correlation", method="ward.D", plot=FALSE)

PCASamples(meth_10, screeplot=TRUE)
PCASamples(meth_10)

### ggplot PCA

pca.10=prcomp(t(pm), center = T)
summary(pca.10)
# Importance of components:
#                            PC1     PC2      PC3      PC4      PC5
# Standard deviation     38.3089 31.1476 25.79283 24.87515 23.22819
# Proportion of Variance  0.1544  0.1021  0.07001  0.06512  0.05678
# Cumulative Proportion   0.1544  0.2565  0.32654  0.39166  0.44844

# Add population info to pca matrix
df=as.data.frame(pca.10$x)
info$Population=as.character(info$Population)

# PCA in all individuals-------
# calculate variance explained by each PC
var.all.10=as.data.frame(pca.10$x)
var.all.10=apply(var.all.10, 2, var)
var.all.10=var.all.10/sum(var.all.10)
var.all.10[1:2]
#      PC1       PC2 
# 0.1544378 0.1020947  

axes.10=as.data.frame(pca.10$x)
# only use PC1 and PC2
axes.10=axes.10[,c(1,2)]
# only use PC2 and PC3
#axes.10=axes.10[,c(1,2)]
axes.all.10=axes.10
axes.all.10$Label=rownames(axes.all.10)

updated_info <- dplyr::select(info, Environment, Population, FishID, Batch, Exposure, Label)

axes.all.10=merge(axes.all.10, updated_info, by = 'Label', all.x=TRUE)

axes.all.10 <- axes.all.10 %>% 
  dplyr::select(., Label, PC1, PC2, Environment, Population, Batch) %>% 
  mutate(Environment = c("Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden", "Field","Field","Field","Field","Field","Field","Field","Field", "Field", "Field", "Field", "Field", "Field", "Field"), Population = c("Gosling", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Roberts", "Sayward", "Sayward", "Sayward", "Sayward", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson"), Batch = c("1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2",  "1", "1", "1", "1", "1", "1", "1", "1"), Sex = c("Male", "Female", "Female" ,"Male", "Female", "Male", "Male", "Male", "Male", "Female", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Unknown", "Male", "Male", "Male", "Female", "Female", "Unknown"), KnownExposure = c("Y", "Y", "N", "N", "Y", "Y", "N", "N", "Y", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"))

axes.all.10$Population=factor(axes.all.10$Population, levels = c("Gosling", "Roberts", "Sayward", "Wik", "Watson"), ordered = F)

# Plot PC1 and PC2

g1=ggplot(axes.all.10, aes(PC1, PC2))
g2=g1+geom_point(aes(color=KnownExposure, fill=KnownExposure, shape = Batch), size=6)+
  labs(x="PC1 (15.4%)", y="PC2 (10.2%)")
g3=g2 + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.text = element_text(size = 14, face = "bold"),
                         axis.title = element_text(size = 16, face = "bold"),
                         legend.text = element_text(size = 14, face = "bold"),
                         legend.title = element_text(size = 16, face = "bold")
)
g3

# only use PC2 and PC3
var.all.10[2:3]
#        PC2        PC3 
# 0.10209473 0.07000861 

axes.11=as.data.frame(pca.10$x)
# only use PC1 and PC2
axes.11=axes.11[,c(2,3)]
axes.all.11=axes.11
axes.all.11$Label=rownames(axes.all.11)

axes.all.11=merge(axes.all.11, updated_info, by = 'Label', all.x=TRUE)

axes.all.11 <- axes.all.11 %>% 
  dplyr::select(., Label, PC2, PC3, Environment, Population, Batch) %>% 
  mutate(Environment = c("Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Field","Field","Field","Field","Field","Field","Field","Field", "Field", "Field", "Field", "Field", "Field", "Field"), Population = c("Gosling", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Roberts", "Sayward", "Sayward", "Sayward", "Sayward", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson"), Batch = c("1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2",  "1", "1", "1", "1", "1", "1", "1", "1"), Sex = c("Male", "Female", "Female" ,"Male", "Female", "Male", "Male", "Male", "Male", "Female", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Unknown", "Male", "Male", "Male", "Female", "Female", "Unknown"), KnownExposure = c("Y", "Y", "N", "N", "Y", "Y", "N", "N", "Y", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"))

axes.all.10$Population=factor(axes.all.10$Population, levels = c("Gosling", "Roberts", "Sayward", "Wik", "Watson"), ordered = F)

#Plot PC2 and PC3

g4=ggplot(axes.all.11, aes(PC2, PC3))
g5=g4+geom_point(aes(color=KnownExposure, fill=KnownExposure, shape=Batch), size=6)+
  labs(x="PC2 (10.2%)", y="PC3 (7.0%)")
g6=g5 + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.text = element_text(size = 14, face = "bold"),
                         axis.title = element_text(size = 16, face = "bold"),
                         legend.text = element_text(size = 14, face = "bold"),
                         legend.title = element_text(size = 16, face = "bold")
)
g6

# only use PC3 and PC4
var.all.10[3:4]
#        PC3        PC4 
# 0.07000861 0.06511561 

axes.12=as.data.frame(pca.10$x)
# only use PC3 and PC4
axes.12=axes.12[,c(3,4)]
axes.all.12=axes.12
axes.all.12$Label=rownames(axes.all.12)

axes.all.12=merge(axes.all.12, updated_info, by = 'Label', all.x=TRUE)

axes.all.12 <- axes.all.12 %>% 
  dplyr::select(., Label, PC3, PC4, Environment, Population, Batch) %>% 
  mutate(Environment = c("Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Field","Field","Field","Field","Field","Field","Field","Field", "Field", "Field", "Field", "Field", "Field", "Field"), Population = c("Gosling", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Roberts", "Sayward", "Sayward", "Sayward", "Sayward", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson"), Batch = c("1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2",  "1", "1", "1", "1", "1", "1", "1", "1"), Sex = c("Male", "Female", "Female" ,"Male", "Female", "Male", "Male", "Male", "Male", "Female", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Unknown", "Male", "Male", "Male", "Female", "Female", "Unknown"), KnownExposure = c("Y", "Y", "N", "N", "Y", "Y", "N", "N", "Y", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"))

axes.all.12$Population=factor(axes.all.12$Population, levels = c("Gosling", "Roberts", "Sayward", "Wik", "Watson"), ordered = F)

#Plot PC3 and PC4

g4=ggplot(axes.all.12, aes(PC3, PC4))
g5=g4+geom_point(aes(color=KnownExposure, fill=KnownExposure, shape=Batch), size=6)+
  labs(x="PC3 (7.0%)", y="PC4 (6.5%)")
g6=g5 + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.text = element_text(size = 14, face = "bold"),
                         axis.title = element_text(size = 16, face = "bold"),
                         legend.text = element_text(size = 14, face = "bold"),
                         legend.title = element_text(size = 16, face = "bold")
)
g6

# only use PC4 and PC5
var.all.10[4:5]
#        PC3        PC4 
# 0.06511561 0.05677855 

axes.13=as.data.frame(pca.10$x)
# only use PC3 and PC4
axes.13=axes.13[,c(4,5)]
axes.all.13=axes.13
axes.all.13$Label=rownames(axes.all.13)

axes.all.13=merge(axes.all.13, updated_info, by = 'Label', all.x=TRUE)

axes.all.13 <- axes.all.13 %>% 
  dplyr::select(., Label, PC4, PC5, Environment, Population, Batch) %>% 
  mutate(Environment = c("Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden","Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Common Garden", "Field","Field","Field","Field","Field","Field","Field","Field", "Field", "Field", "Field", "Field", "Field", "Field"), Population = c("Gosling", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Roberts", "Sayward", "Sayward", "Sayward", "Sayward", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson", "Gosling", "Gosling", "Gosling", "Roberts", "Roberts", "Roberts", "Wik", "Wik", "Wik", "Wik", "Watson", "Watson", "Watson", "Watson"), Batch = c("1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "2", "2", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2",  "1", "1", "1", "1", "1", "1", "1", "1"), Sex = c("Male", "Female", "Female" ,"Male", "Female", "Male", "Male", "Male", "Male", "Female", "Male", "Male", "Male", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Unknown", "Male", "Male", "Male", "Female", "Female", "Unknown"), KnownExposure = c("Y", "Y", "N", "N", "Y", "Y", "N", "N", "Y", "Y", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"))

axes.all.13$Population=factor(axes.all.13$Population, levels = c("Gosling", "Roberts", "Sayward", "Wik", "Watson"), ordered = F)

#Plot PC3 and PC4

g4=ggplot(axes.all.13, aes(PC4, PC5))
g5=g4+geom_point(aes(color=KnownExposure, fill=KnownExposure, shape=Batch), size=6)+
  labs(x="PC4 (6.5%)", y="PC5 (5.7%)")
g6=g5 + theme_bw()+theme(panel.grid=element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(),
                         axis.text = element_text(size = 14, face = "bold"),
                         axis.title = element_text(size = 16, face = "bold"),
                         legend.text = element_text(size = 14, face = "bold"),
                         legend.title = element_text(size = 16, face = "bold")
)
g6
