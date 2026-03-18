### Sequencing batch effects for threespine stickleback - Supplemental Figures 2 & 3


# Import methylation data -------------------------------------------------



setwd("/MethylDackel/Stickle")

file.list.10 <- list(
  "CG_22_GOS_001_CpG.methylKit", "CG_22_GOS_002_CpG.methylKit", "CG_GOS_35_04_CpG.methylKit", "CG_GOS_2210_01_CpG.methylKit", "FD_GOS_27_CpG.methylKit", "FD_GOS_28_CpG.methylKit", "FD_GOS_33_CpG.methylKit", "CG_22_ROB_001_CpG.methylKit", "CG_22_ROB_002_CpG.methylKit", "CG_ROB_003_001_CpG.methylKit", "CG_ROB_31_01_001_CpG.methylKit", "FD_ROB_98_001_CpG.methylKit", "FD_ROB_102_001_CpG.methylKit", "FD_ROB_106_001_CpG.methylKit", "CG_22_SAY_001_CpG.methylKit", "CG_22_SAY_002_CpG.methylKit", "CG_SAY_4_01_001_CpG.methylKit", "CG_SAY_2343_003_001_CpG.methylKit", "CG_22_WK_001_CpG.methylKit", "CG_22_WK_002_CpG.methylKit", "CG_22_WK_003_CpG.methylKit", "CG_22_WK_004_CpG.methylKit", "FD_22_WK_002_CpG.methylKit", "FD_22_WK_005_CpG.methylKit", "FD_22_WK_011_CpG.methylKit", "FD_22_WK_012_CpG.methylKit", "CG_22_WT_001_CpG.methylKit", "CG_22_WT_002_CpG.methylKit", "CG_22_WT_003_CpG.methylKit", "CG_22_WT_004_CpG.methylKit", "FD_22_WT_004_CpG.methylKit", "FD_22_WT_007_CpG.methylKit", "FD_22_WT_016_CpG.methylKit", "FD_22_WT_017_CpG.methylKit"
)

myobj.10=methRead(file.list.10,
                  sample.id=list("CG_GOS_001", "CG_GOS_002", "CG_GOS_003", "CG_GOS_004", "FD_GOS_027", "FD_GOS_028", "FD_GOS_033", "CG_ROB_001", "CG_ROB_002", "CG_ROB_003", "CG_ROB_004", "FD_ROB_098", "FD_ROB_102", "FD_ROB_106", "CG_SAY_001", "CG_SAY_002", "CG_SAY_003", "CG_SAY_004", "CG_WK_001", "CG_WK_002", "CG_WK_003", "CG_WK_004", "FD_WK_002", "FD_WK_005", "FD_WK_011", "FD_WK_012","CG_WT_001", "CG_WT_002", "CG_WT_003", "CG_WT_004", "FD_WT_004", "FD_WT_007", "FD_WT_016", "FD_WT_017"),
                  pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                  assembly="StickleGeneAnnotations",
                  treatment=c(1,1,2,2,2,2,2,1,1,2,2,2,2,2,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                  context="CpG",
                  mincov = 5
)


# PCA & Clusters based on DNAm data ---------------------------------------



filtered.myobjDepth10=filterByCoverage(myobj.10,lo.count=5,lo.perc=NULL,
                                       hi.count=NULL,hi.perc=99.9)

normalized_myobjDepth10=normalizeCoverage(filtered.myobjDepth10, method="median")

meth.10=methylKit::unite(normalized_myobjDepth10, min.per.group = 2L, destrand=TRUE)

pm=percMethylation(meth.10)
nrow(pm) # 1124550

clusterSamples(meth.10, dist="correlation", method="ward.D", plot=TRUE)

PCASamples(meth.10, screeplot=TRUE)
PCASamples(meth.10, screeplot=FALSE)
PCASamples(meth.10, screeplot=FALSE, comp = c(2,3))


# Mapping efficiency ------------------------------------------------------


## Import read mapping data ------------------------------------------------



setwd("/MethylMethods")

dat <- read_excel("QCStats_1.xlsx", sheet = "BamtoolsStats")
{
  dat <- filter(dat, !AlignmentMethod == "bwa mem")
  
  dat$PercentFailedQC[is.na(dat$PercentFailedQC)] <- 0
  
  dat$PercMappedErr <- dat$PercentMapped - dat$PercentFailedQC
  
  dat$AlignmentMethod <- factor(dat$AlignmentMethod, levels = c("Bismark","BismarkLocal","BisulfiteBolt", "Biscuit", "BWA meth"))


Bismark <- dat %>% 
  filter(AlignmentMethod == "Bismark" & Organism == "Stickleback")

BismarkL <- dat %>% 
  filter(AlignmentMethod == "BismarkLocal" & Organism == "Stickleback")

bwameth <- dat %>% 
  filter(AlignmentMethod == "BWA meth" & Organism == "Stickleback")

Biscuit <- dat %>% 
  filter(AlignmentMethod == "Biscuit" & Organism == "Stickleback")

BSBolt <- dat %>% 
  filter(AlignmentMethod == "BisulfiteBolt" & Organism == "Stickleback")
}

## Bismark -----------------------------------------------------------------



# Shapiro-Wilk normality test
shapiro.test(Bismark$PercentMapped)
# data:  Bismark$PercentMapped
# W = 0.91734, p-value = 0.01359 

# Welch Two Sample t-test
wilcox.test(PercentMapped ~ Batch, data = Bismark)
# Wilcoxon rank sum test with continuity correction
# 
# data:  PercentMapped by Batch
# W = 142, p-value = 0.732
# alternative hypothesis: true location shift is not equal to 0
# 
# Warning message:
#   In wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...) :
#   cannot compute exact p-value with ties

Bis.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = Bismark,
  method = "wilcox.test"
)

Bismark$Batch <- as.factor(Bismark$Batch)

Bism <- ggplot(Bismark) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "Bismark", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(Bis.wilcox.test, label = "p", y.position = 0.7, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")
Bism

## Bismark Local --------------------------------------------------------

# Welch Two Sample t-test
wilcox.test(PercentMapped ~ Batch, data = BismarkL)
# Wilcoxon rank sum test with continuity correction
# 
# data:  PercentMapped by Batch
# W = 193, p-value = 0.006319
# alternative hypothesis: true location shift is not equal to 0
# 
# Warning message:
#   In wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...) :
#   cannot compute exact p-value with ties

BisL.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = BismarkL,
  method = "wilcox.test"
)

BismarkL$Batch <- as.factor(BismarkL$Batch)

BismL <- ggplot(BismarkL) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "Bismark Local", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BisL.wilcox.test, label = "p", y.position = 0.9, step.increase = 1) +
  theme_cowplot() +
  theme(legend.position = "none")
BismL


## BWA meth ----------------------------------------------------------------

bwameth$PercentFailedQC[is.na(bwameth$PercentFailedQC)] <- 0
  
bwameth$PercMappedErr <- bwameth$PercentMapped - bwameth$PercentFailedQC

bwameth$Batch <- as.factor(bwameth$Batch)

# Shapiro-Wilk normality test
shapiro.test(bwameth$PercMappedErr)
# data:  bwameth$PercentMapped
#W = 0.92718, p-value = 0.02591

# Wilcoxon rank sum exact test
wilcox.test(PercMappedErr~Batch, data = bwameth) 
# data:  PercentMapped by Batch
# W = 0, p-value = 3.647e-09
# alternative hypothesis: true location shift is not equal to 0
# p<0.05, so there is a significant difference between batch 1 and batch 2 mapping efficiencies

BWAMeth.wilcox.test <- compare_means(
  PercMappedErr ~ Batch,
  data = bwameth,
  method = "wilcox.test"
)

bwam <- ggplot(bwameth) +
  geom_boxplot(aes(x = Batch, y = PercMappedErr, fill = Batch)) +
  labs(title = "BWA Meth", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMeth.wilcox.test, label = "p", y.position = 1.0, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")
bwam


## Biscuit --------------------------------------------------------------


# Wilcoxon rank sum exact test
wilcox.test(PercentMapped~Batch, data = Biscuit) 
# data:  PercentMapped by Batch
# W = 244, p-value = 9.461e-06
# alternative hypothesis: true location shift is not equal to 0
# p<0.05, so there is a significant difference between batch 1 and batch 2 mapping efficiencies

Biscuit.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = Biscuit,
  method = "wilcox.test"
)

Biscuit$Batch <- as.factor(Biscuit$Batch)

bisc <- ggplot(Biscuit) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "Biscuit", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMeth.wilcox.test, label = "p", y.position = 1.0, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")
bisc


## BSBolt ---------------------------------------------------------------


# Wilcoxon rank sum exact test
wilcox.test(PercentMapped~Batch, data = BSBolt) 
# data:  PercentMapped by Batch
# W = 260, p-value = 4.377e-08
# alternative hypothesis: true location shift is not equal to 0

BSBolt.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = BSBolt,
  method = "wilcox.test"
)

BSBolt$Batch <- as.factor(BSBolt$Batch)

bolt <- ggplot(BSBolt) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "BiSulfite Bolt", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMeth.wilcox.test, label = "p", y.position = 1.0, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")
bolt

Bism + BismL + bwam + bisc + bolt
