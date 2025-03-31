#### Test for differrences in mapping efficiency and number of mapped reads between Bismark, BWA meth, and BWA mem

# import data

library(readxl)
library(tidyverse)
library(ggplot2)
library(car)
library(ggpubr)
library(dplyr)
library(cowplot)
library(coin)
library(rstatix)


MappingEfficiency <- read_excel("R:/Genohub_seq/methylation_methods_comparison/Metadata/QC/QCStats_1.xlsx", sheet = "BamtoolsStats")

dat <- MappingEfficiency %>% 
  select("AlignmentMethod", "PercentMapped", "Batch", "Population", "Environment","SampleID")

dat <- na.omit(dat)

################## Mapping Efficiency ##################

######## Visualize data

ggplot(dat) +
  aes(x = AlignmentMethod, y = PercentMapped, color = Population) +
  geom_jitter()

ggplot(dat) +
  aes(x = AlignmentMethod, y = PercentMapped, color = Environment) +
  geom_jitter()

ggplot(dat) +
  aes(x = AlignmentMethod, y = PercentMapped, color = AlignmentMethod) +
  geom_jitter() +
  theme(legend.position = "none")

ggplot(dat) +
  aes(x = AlignmentMethod, y = PercentMapped, color = AlignmentMethod) +
  geom_boxplot() +
  theme(legend.position = "none")


####### Test for differences between read mapping methods

shapiro.test(dat$PercentMapped)
#W = 0.73905, p-value = 3.57e-12

# histogram
hist(dat$PercentMapped)

#### Does not meet assumption of normality, using Friedman Test

res.fried <- friedman_test(dat, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#    .y.               n      statistic    df        p          method       
#   1 PercentMapped    34      60.9         2       5.85e-14  Friedman test

dat %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n   effsize method      magnitude
#   1 PercentMapped    34   0.896   Kendall W     large 

# pairwise comparisons
pwc <- dat %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1  group2      n1    n   statistic    p          p.adj    p.adj.signif
# PercentMapped Bismark bwa mem     34    34   548      2.5 e- 6    7.5 e- 6      ****        
# PercentMapped Bismark bwa meth    34    34   0        1.16e-10    3.48e-10      ****        
# PercentMapped bwa mem bwa meth    34    34   0        1.16e-10    3.48e-10      ****     

dat %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.             group1  group2   effsize    n1    n2 magnitude
# 1 PercentMapped Bismark bwa mem    0.734    34    34 large    
# 2 PercentMapped Bismark bwa meth   0.872    34    34 large    
# 3 PercentMapped bwa mem bwa meth   0.872    34    34 large

######### Make a nice plot displaying significant differences between treatments

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "AlignmentMethod")
ggboxplot(dat, x = "AlignmentMethod", y = "PercentMapped", add = "point") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(x = "Alignment Method", y = "Mapping Efficiency"
    )
pwc

####### Mean and SD of mapping efficiency

by_group <- group_by(dat, AlignmentMethod)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# AlignmentMethod  mean   sd
# Bismark         0.544 0.0435 
# bwa mem         0.493 0.0446 
# bwa meth        0.992 0.00181

# Bismark batch 1 & batch 2
Bismark <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "Bismark")

by_group <- group_by(Bismark, Batch)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# Batch  mean     sd
# 1     0.545 0.0447
# 2     0.542 0.0431

# bwa meth batch 1 & batch 2
BwaMeth <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa meth")

by_group <- group_by(BwaMeth, Batch)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# Batch  mean       sd
# 1     0.993   0.000864
# 2     0.990   0.00127

# bwa mem batch 1 & batch 2
BwaMem <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa mem")

by_group <- group_by(BwaMem, Batch)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# Batch  mean     sd
#   1   0.517   0.0301
# 2     0.447   0.0279

###### Test for batch effects on mapping efficiency using each mapping tool

### Bismark

# Shapiro-Wilk normality test
shapiro.test(Bismark$PercentMapped)
# data:  Bismark$PercentMapped
# W = 0.91734, p-value = 0.01359 - p <0.05, not normally distributed so we need to use Wilcox test

# Wilcox ranks sum test
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

ggplot(Bismark) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "Bismark", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(Bis.wilcox.test, label = "p", y.position = 0.7, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")


### BWA meth
# Shapiro-Wilk normality test
shapiro.test(BwaMeth$PercentMapped)
# data:  bwameth$PercentMapped
# W = 0.9292, p-value = 0.02965 - p<0.05, data not normally distributed, may need to use a nonparametric test

# Wilcoxon rank sum exact test
wilcox.test(PercentMapped~Batch, data = BwaMeth) 
# data:  PercentMapped by Batch
# W = 261, p-value = 2.553e-08
# alternative hypothesis: true location shift is not equal to 0
# p<0.05, so there is a significant difference between batch 1 and batch 2 mapping efficiencies

BWAMeth.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = BwaMeth,
  method = "wilcox.test"
)

BwaMeth$Batch <- as.factor(BwaMeth$Batch)

ggplot(BwaMeth) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "BWA Meth", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMeth.wilcox.test, label = "p", y.position = 1.0, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")

### BWA mem
# Shapiro-Wilk normality test
shapiro.test(BwaMem$PercentMapped)
# data:  BwaMem$PercentMapped
# W = 0.95749, p-value = 0.2052 - p>0.05, data normally distributed, but will use Wilcox test for consistency

# Wilcox test
wilcox.test(PercentMapped~Batch, data = BwaMem) 
# W = 253, p-value = 7.112e-07

BWAMem.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = BwaMem,
  method = "wilcox.test"
)

BwaMem$Batch <- as.factor(BwaMem$Batch)

ggplot(BwaMem) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "BWA Mem", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMem.wilcox.test, label = "p", y.position = 0.7, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")


################## Mapped Reads ##################

####### Test for differences between read mapping methods

dat2 <- MappingEfficiency %>% 
  select("AlignmentMethod", "MappedReads", "Batch", "Population", "Environment", "SampleID")

dat2 <- na.omit(dat2)

shapiro.test(as.numeric(dat2$MappedReads))
#W = 0.68734, p-value = 1.989e-13

# histogram
hist(as.numeric(dat2$MappedReads))

#### Does not meet assumption of normality, using Friedman Test

dat2$MappedReads <- as.numeric(dat2$MappedReads)

res.fried <- dat2 %>% friedman_test(MappedReads ~ AlignmentMethod |SampleID)
res.fried
#    .y.               n      statistic    df        p          method       
#    MappedReads    34        68            2      1.71e-15     Friedman test

dat2 %>% friedman_effsize(MappedReads ~ AlignmentMethod |SampleID)
#     .y.               n   effsize method      magnitude
#   MappedReads         34       1  Kendall W     large  

# pairwise comparisons
pwc <- dat2 %>%
  wilcox_test(MappedReads ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1    group2      n1    n   statistic    p       p.adj    p.adj.signif
#   MappedReads Bismark   bwa mem     34    34         0  1.16e-10  3.48e-10  ****        
#   MappedReads Bismark   bwa meth    34    34         0  1.16e-10  3.48e-10  ****        
#   MappedReads bwa mem   bwa meth    34    34         0  1.16e-10  3.48e-10  ****     


######### Make a nice plot displaying significant differences between treatments

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "AlignmentMethod")
ggboxplot(dat2, x = "AlignmentMethod", y = "MappedReads", add = "point") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(x = "Alignment Method", y = "Mapped Reads"
  )
pwc


##### Mean & SD of Mapped Reads

by_group <- group_by(dat2, AlignmentMethod)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# AlignmentMethod    mean        sd
# Bismark          4875224.  3741376.
# bwa mem          8586430.  6110385.
# bwa meth        17790041. 13454374.

# Bismark batch 1 & 2
Bismark_reads <- dat2 %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "Bismark")

by_group <- group_by(Bismark_reads, Batch)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# Batch   mean       sd
#   1   2855790.  448334.
# 2     8577520. 4301236.

# bwa meth batch 1 & 2
BwaMeth_reads <- dat2 %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa meth")

by_group <- group_by(BwaMeth_reads, Batch)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# Batch   mean        sd
#   1   10428487.  1268689.
# 2     31286224. 15260964.

# bwa mem batch 1 & 2
BwaMem_reads <- dat2 %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa mem")

by_group <- group_by(BwaMem_reads, Batch)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# Batch   mean       sd
#   1   5443143.  844552.
# 2     14349122. 7393282.

###### Test for batch effects on mapped reads using each mapping tool

### Bismark

# Shapiro-Wilk normality test
shapiro.test(Bismark_reads$MappedReads)
# data:  Bismark_mappedreads$MappedReads
# W = 0.68505, p-value = 2.926e-07 - p<0.05, need to use non parametric test

# Wilcoxon rank sum exact test
Bismark.wilcox <- wilcox.test(MappedReads ~ Batch, data = Bismark_reads)
Bismark.wilcox
# data:  MappedReads by Batch
# W = 0, p-value = 3.647e-09
# alternative hypothesis: true location shift is not equal to 0

Bismark.wilcox.test <- compare_means(
  MappedReads ~ Batch,
  data = Bismark_reads,
  method = "wilcox.test"
)

Bismark_reads$Batch <- as.factor(Bismark_reads$Batch)

ggplot(Bismark_reads) +
  geom_boxplot(aes(x = Batch, y = MappedReads, fill = Batch)) +
  labs(title = "Bismark", x = "Sequencing Batch", y = "Mapped Reads") +
  stat_pvalue_manual(Bismark.wilcox.test, label = "p", y.position = 20000000, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")


### BWA meth
# Shapiro-Wilk normality test
shapiro.test(BwaMeth_reads$MappedReads)
# data:  BwaMeth_reads$MappedReads
# W = 0.67436, p-value = 2.036e-07 - p<0.05, need to use non parametric test

# Wilcoxon rank sum exact test
bwameth.wilcox.reads <- wilcox.test(MappedReads ~ Batch, data = BwaMeth_reads)
bwameth.wilcox.reads
# data:  MappedReads by Batch
# W = 0, p-value = 3.647e-09
# alternative hypothesis: true location shift is not equal to 0

bwameth.wilcox.reads <- compare_means(
  MappedReads ~ Batch,
  data = BwaMeth_reads,
  method = "wilcox.test"
)

BwaMeth_reads$Batch <- as.factor(BwaMeth_reads$Batch)

ggplot(BwaMeth_reads) +
  geom_boxplot(aes(x = Batch, y = MappedReads, fill = Batch)) +
  labs(title = "BWA meth", x = "Sequencing Batch", y = "Mapped Reads") +
  stat_pvalue_manual(bwameth.wilcox.reads, label = "p", y.position = 75000000, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")

### BWA mem
# Shapiro-Wilk normality test
shapiro.test(BwaMem_reads$MappedReads)
# data:  BwaMem_reads$MappedReads
# W = 0.66366, p-value = 1.426e-07 - p<0.05, need to use non parametric test

# Wilcoxon rank sum exact test
bwamem.wilcox.reads <- wilcox.test(MappedReads ~ Batch, data = BwaMem_reads)
bwamem.wilcox.reads
# data:  MappedReads by Batch
# W = 2, p-value = 1.459e-08
# alternative hypothesis: true location shift is not equal to 0

bwamem.wilcox.reads <- compare_means(
  MappedReads ~ Batch,
  data = BwaMem_reads,
  method = "wilcox.test"
)

BwaMem_reads$Batch <- as.factor(BwaMem_reads$Batch)

ggplot(BwaMem_reads) +
  geom_boxplot(aes(x = Batch, y = MappedReads, fill = Batch)) +
  labs(title = "BWA mem", x = "Sequencing Batch", y = "Mapped Reads") +
  stat_pvalue_manual(bwamem.wilcox.reads, label = "p", y.position = 36000000, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")


##################################################
################ WGBS ############################
##################################################

MappingEfficiency <- read_excel("R:/Genohub_seq/methylation_methods_comparison/Metadata/QC/QCStats_1.xlsx", sheet = "WBGSBamtoolsStats")

dat <- MappingEfficiency %>%
  select("Alignment Method", "PercentMapped", "Population")

dat <- na.omit(dat)

dat$AlignmentMethod <- dat$`Alignment Method`

# Visualize the data

ggplot(dat) +
  aes(x = AlignmentMethod, y = PercentMapped, color = Population) +
  geom_jitter()

ggplot(dat) +
  aes(x = AlignmentMethod, y = PercentMapped, color = AlignmentMethod) +
  geom_jitter() +
  theme(legend.position = "none")

ggboxplot(dat, x = "AlignmentMethod", y = "PercentMapped", add = "point") +
  labs(x = "Alignment Method", y = "Mapping Efficiency"
  )

####### Too few WGBS samples to do statistical tests, just report mean and SD 

by_group <- group_by(dat, `Alignment Method`)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# 
# `Alignment Method`  mean      sd
# <chr>              <dbl>   <dbl>
#   1 Bismark         0.387   0.0124 
# 2 bwa mem           0.746   0.0343 
# 3 bwa meth          0.996   0.00157

####### Mean and SD of mapped reads

dat <- MappingEfficiency %>%
  select("Alignment Method", "MappedReads", "Population")

dat <- na.omit(dat)

dat$AlignmentMethod <- dat$`Alignment Method`

by_group <- group_by(dat, `Alignment Method`)

by_group %>% summarise_each_(funs(mean(., na.rm = TRUE), sd(., na.rm=TRUE)), names(by_group)[2])
# `Alignment Method`      mean       sd
# <chr>                  <dbl>    <dbl>
# 1 Bismark            11321912.  502776.
# 2 bwa mem            44017088  2656871.
# 3 bwa meth           58620964. 1107988.