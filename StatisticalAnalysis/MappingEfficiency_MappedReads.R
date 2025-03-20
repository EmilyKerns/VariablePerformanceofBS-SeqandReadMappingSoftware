#### Test for differrences in mapping efficiency and number of mapped reads between Bismark, BWA meth, and BWA mem

# import data

library(readxl)
library(tidyverse)
library(ggplot2)
library(car)
library(FSA)
library(ggstatsplot)
library(ggpubr)
library(dplyr)
library(cowplot)


MappingEfficiency <- read_excel("QCStats_1.xlsx", sheet = "BamtoolsStats")

dat <- MappingEfficiency %>% 
  select("AlignmentMethod", "PercentMapped", "Batch", "Population", "Environment")

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


# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
    data = dat
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
# Terms:
#   AlignmentMethod Residuals
# Sum of Squares         5.140212  0.128049
# Deg. of Freedom               2        99
# 
# Residual standard error: 0.03596426
# Estimated effects may be unbalanced

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.94372, p-value = 0.0002821

#### Does not meet assumption of normality, moving to Kruskal Wallis rank sum test

## Kruskal Wallis Test
kruskal.test(PercentMapped ~ AlignmentMethod,
             data = dat
)

# Kruskal-Wallis rank sum test
# 
# data:  PercentMapped by AlignmentMethod
# Kruskal-Wallis chi-squared = 75.952, df = 2, p-value < 2.2e-16

# post hoc test
dunnTest(PercentMapped ~ AlignmentMethod,
         data = dat,
         method = "holm"
)
# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
#       Comparison         Z      P.unadj        P.adj
# 1 Bismark - bwa mem  2.934500 3.340851e-03 3.340851e-03
# 2 Bismark - bwa meth -5.639487 1.705577e-08 3.411154e-08
# 3 bwa mem - bwa meth -8.573987 9.996252e-18 2.998876e-17
# Warning message:
#   AlignmentMethod was coerced to a factor. 


######### Make a nice plot displaying p value


dat$Batch <- as.character(dat$Batch)

p <- ggboxplot(dat, x = "AlignmentMethod", y = "PercentMapped",
          color = "Batch", palette = "jco", xlab = "Alignment Method", ylab = "Mapping Efficiency", add = "jitter")+
  theme(text=element_text(size=14))+
  stat_compare_means(method = "kruskal.test") # Pairwise comparison against all
 
ggpar(p, ylim = c(0,1.1))
   
Figure <- ggboxplot(dat, x = "AlignmentMethod", y = "PercentMapped",
  color = "AlignmentMethod", palette = "jco", xlab = "Alignment Method", ylab = "Mapping Efficiency", add = "jitter")+
  theme(text=element_text(size=14), legend.position = "none")+
  stat_compare_means(method = "kruskal.test") # Pairwise comparison against all
   
ggpar(Figure, ylim = c(0,1.1))
   
q <- ggscatter(dat, x = "AlignmentMethod", y = "PercentMapped",
          color = as.character("Batch"), palette = "jco")+
  stat_compare_means(method = "kruskal.test")+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")+   # Pairwise comparison against all
  aes(xlab = "Alignment Method", ylab = "Mapping Efficiency")

q

####### Mean and SD of mapping efficiency

# bwa meth
bwameth <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa meth")

mean(bwameth$PercentMapped) 
# 0.9923139
sd(bwameth$PercentMapped)
# 0.0018085
sd(bwameth$PercentMapped)/sqrt(length((bwameth))) #standard error
# 0.001044138


# Bismark
Bismark <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "Bismark")

mean(Bismark$PercentMapped)
# 0.5437647
sd(Bismark$PercentMapped)
# 0.04348114
sd(Bismark$PercentMapped)/sqrt(length((Bismark))) #standard error
# 0.02510385

# bwa mem
bwamem <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa mem")

mean(bwamem$PercentMapped)
# 0.492582
sd(bwamem$PercentMapped)
# 0.0445691
sd(bwamem$PercentMapped)/sqrt(length((bwamem))) #standard error
# 0.02573198

# Bismark batch 1 & batch 2
Bismark_1 <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "Bismark") %>% 
  filter(Batch == "1")

mean(Bismark_1$PercentMapped)
# 0.5448636
sd(Bismark_1$PercentMapped)
# 0.04466307
sd(Bismark_1$PercentMapped)/sqrt(length((Bismark_1))) #standard error
# 0.02578624

Bismark_2 <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "Bismark") %>% 
  filter(Batch == "2")

mean(Bismark_2$PercentMapped)
# 0.54175
sd(Bismark_2$PercentMapped)
# 0.04309002
sd(Bismark_2$PercentMapped)/sqrt(length((Bismark_2))) #standard error
# 0.02487804

# bwa meth batch 1 & batch 2
BwaMeth_1 <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa meth") %>% 
  filter(Batch == "1")

mean(BwaMeth_1$PercentMapped)
# 0.9934065
sd(BwaMeth_1$PercentMapped)
# 0.0008636809
sd(BwaMeth_1$PercentMapped)/sqrt(length((BwaMeth_1))) #standard error
# 0.0004986464

BwaMeth_2 <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa meth") %>% 
  filter(Batch == "2")

mean(BwaMeth_2$PercentMapped)
# 0.990311
sd(BwaMeth_2$PercentMapped)
# 0.001274477
sd(BwaMeth_2$PercentMapped)/sqrt(length((BwaMeth_2))) #standard error
# 0.0007358194

# bwa mem batch 1 & batch 2
BwaMem_1 <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa mem") %>% 
  filter(Batch == "1")

mean(BwaMem_1$PercentMapped)
# 0.5172608
sd(BwaMem_1$PercentMapped)
# 0.03007624
sd(BwaMem_1$PercentMapped)/sqrt(length((BwaMem_1))) #standard error
# 0.01736452

BwaMem_2 <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped, Batch) %>% 
  filter(AlignmentMethod == "bwa mem") %>% 
  filter(Batch == "2")

mean(BwaMem_2$PercentMapped)
# 0.4473376
sd(BwaMem_2$PercentMapped)
# 0.02794727
sd(BwaMem_2$PercentMapped)/sqrt(length((BwaMem_2))) #standard error
# 0.01613537

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
shapiro.test(bwameth$PercentMapped)
# data:  bwameth$PercentMapped
# W = 0.9292, p-value = 0.02965 - p<0.05, data not normally distributed, may need to use a nonparametric test

# Wilcoxon rank sum exact test
wilcox.test(PercentMapped~Batch, data = bwameth) 
# data:  PercentMapped by Batch
# W = 261, p-value = 2.553e-08
# alternative hypothesis: true location shift is not equal to 0
# p<0.05, so there is a significant difference between batch 1 and batch 2 mapping efficiencies

BWAMeth.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = bwameth,
  method = "wilcox.test"
)

bwameth$Batch <- as.factor(bwameth$Batch)

ggplot(bwameth) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "BWA Meth", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMeth.wilcox.test, label = "p", y.position = 1.0, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")

### BWA mem
# Shapiro-Wilk normality test
shapiro.test(bwamem$PercentMapped)
# data:  bwamem$PercentMapped
# W = 0.95749, p-value = 0.2052 - p>0.05, data normally distributed, but will use Wilcox test for consistency

# Wilcox test
wilcox.test(PercentMapped~Batch, data = bwamem) 


BWAMem.wilcox.test <- compare_means(
  PercentMapped ~ Batch,
  data = bwamem,
  method = "wilcox.test"
)

bwamem$Batch <- as.factor(bwamem$Batch)

ggplot(bwamem) +
  geom_boxplot(aes(x = Batch, y = PercentMapped, fill = Batch)) +
  labs(title = "BWA Mem", x = "Sequencing Batch", y = "Mapping Efficiency") +
  stat_pvalue_manual(BWAMem.wilcox.test, label = "p", y.position = 0.7, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "none")


################## Mapped Reads ##################

####### Test for differences between read mapping methods

dat2 <- MappingEfficiency %>% 
  select("AlignmentMethod", "MappedReads", "Batch", "Population", "Environment")

# ANOVA
ANOVA <- aov(MappedReads ~ AlignmentMethod,
             data = dat2
)
ANOVA
# Call:
#   aov(formula = MappedReads ~ AlignmentMethod, data = dat2)
# 
# Terms:
# AlignmentMethod    Residuals
# Sum of Squares     3.006416e+15 7.667711e+15
# Deg. of Freedom               2           99
# 
# Residual standard error: 8800661
# Estimated effects may be unbalanced

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# data:  ANOVA$residuals
# W = 0.71788, p-value = 1.047e-12

#### Does not meet assumption of normality, moving to Kruskal Wallis rank sum test

## Kruskal-Wallis rank sum test
kruskal.test(MappedReads ~ AlignmentMethod,
             data = dat2
)
# data:  MappedReads by AlignmentMethod
# Kruskal-Wallis chi-squared = 50.388, df = 2, p-value = 1.144e-11

# post hoc test
dunnTest(MappedReads ~ AlignmentMethod,
         data = dat2,
         method = "holm"
)
# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Holm method.
# 
#   Comparison                Z      P.unadj        P.adj
# 1 Bismark - bwa mem -3.250055 1.153829e-03 1.153829e-03
# 2 Bismark - bwa meth -7.090283 1.338380e-12 4.015140e-12
# 3 bwa mem - bwa meth -3.840228 1.229199e-04 2.458398e-04
# Warning message:
#   AlignmentMethod was coerced to a factor. 


######### Make a nice plot displaying p value


# Visualize: Specify the comparisons you want
# my_comparisons <- list( c("Bismark", "bwa mem"), c("Bismark", "bwa meth"), c("bwa mem", "bwa meth") )

dat2$Batch <- as.character(dat2$Batch)

r <- ggboxplot(dat2, x = "AlignmentMethod", y = "MappedReads",
               color = "Batch", palette = "jco", xlab = "Alignment Method", ylab = "Mapped Reads", add = "jitter")+
  theme(text=element_text(size=14))+
  #stat_compare_means(comparisons = my_comparisons)+ #Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test") # Pairwise comparison against all

ggpar(r, ylim = c(0,70000000))

Figure <- ggboxplot(dat2, x = "AlignmentMethod", y = "MappedReads",
                    color = "AlignmentMethod", palette = "jco", xlab = "Alignment Method", ylab = "Mapped Reads", add = "jitter")+
  theme(text=element_text(size=14), legend.position = "none")+
  #stat_compare_means(comparisons = my_comparisons)+ #Add pairwise comparisons p-value
  stat_compare_means(method = "kruskal.test") # Pairwise comparison against all

ggpar(Figure, ylim = c(0, 70000000))


##### Mean & SD of Mapped Reads


# bwa meth
bwameth_mappedreads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa meth")

mean(bwameth_mappedreads$MappedReads)
# 17,790,041
sd(bwameth_mappedreads$MappedReads)
# 13,454,374

#bwa mem
bwamem_mappedreads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa mem")

mean(bwamem_mappedreads$MappedReads)
# 8,586,430
sd(bwamem_mappedreads$MappedReads)
# 6,110,385

# Bismark
Bismark_mappedreads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "Bismark")

mean(Bismark_mappedreads$MappedReads)
# 4,875,224
sd(Bismark_mappedreads$MappedReads)
# 3,741,376

# Bismark batch 1 & 2
Bismark_1_reads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "Bismark") %>% 
  filter(Batch == "1")
  
mean(Bismark_1_reads$MappedReads)
# 1,450,314
sd(Bismark_1_reads$MappedReads)
# 196,475.8
  
Bismark_2_reads <- MappingEfficiency %>%
    select(AlignmentMethod, MappedReads, Batch) %>% 
    filter(AlignmentMethod == "Bismark") %>% 
    filter(Batch == "2")

mean(Bismark_2_reads$MappedReads)
# 2,855,790
sd(Bismark_2_reads$MappedReads)
# 4,301,236

# bwa meth batch 1 & 2
BwaMeth_1_reads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa meth") %>% 
  filter(Batch == "1")

mean(BwaMeth_1_reads$MappedReads)
# 10,428,487
sd(BwaMeth_1_reads$MappedReads)
# 1,268,689

BwaMeth_2_reads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa meth") %>% 
  filter(Batch == "2")

mean(BwaMeth_2_reads $MappedReads)
# 31,286,224
sd(BwaMeth_2_reads $MappedReads)
# 15,260,964

# bwa mem batch 1 & 2
BwaMem_1_reads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa mem") %>% 
  filter(Batch == "1")

mean(BwaMem_1_reads$MappedReads)
# 5,443,143
sd(BwaMem_1_reads$MappedReads)
# 844,552.3

BwaMem_2_reads <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads, Batch) %>% 
  filter(AlignmentMethod == "bwa mem") %>% 
  filter(Batch == "2")

mean(BwaMem_2_reads$MappedReads)
# 14,349,123
sd(BwaMem_2_reads$MappedReads)
# 7,393,282

###### Test for batch effects on mapped reads using each mapping tool


### Bismark

# Shapiro-Wilk normality test
shapiro.test(Bismark_mappedreads$MappedReads)
# data:  Bismark_mappedreads$MappedReads
# W = 0.68505, p-value = 2.926e-07 - p<0.05, need to use non parametric test

# Wilcoxon rank sum exact test
Bismark.wilcox <- wilcox.test(MappedReads ~ Batch, data = Bismark_mappedreads)
Bismark.wilcox
# data:  MappedReads by Batch
# W = 0, p-value = 3.647e-09
# alternative hypothesis: true location shift is not equal to 0

Bismark.wilcox.test <- compare_means(
  MappedReads ~ Batch,
  data = Bismark_mappedreads,
  method = "wilcox.test"
)

Bismark_mappedreads$Batch <- as.factor(Bismark_mappedreads$Batch)

ggplot(Bismark_mappedreads) +
  geom_boxplot(aes(x = Batch, y = MappedReads, fill = Batch)) +
  labs(title = "Bismark", x = "Sequencing Batch", y = "Mapped Reads") +
  stat_pvalue_manual(Bismark.wilcox.test, label = "p", y.position = 20000000, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")


### BWA meth
# Shapiro-Wilk normality test
shapiro.test(bwameth_mappedreads$MappedReads)
# data:  bwameth_mappedreads$MappedReads
# W = 0.67436, p-value = 2.036e-07 - p<0.05, need to use non parametric test

# Wilcoxon rank sum exact test
bwameth.wilcox.reads <- wilcox.test(MappedReads ~ Batch, data = bwameth_mappedreads)
bwameth.wilcox.reads
# data:  MappedReads by Batch
# W = 0, p-value = 3.647e-09
# alternative hypothesis: true location shift is not equal to 0

bwameth.wilcox.reads <- compare_means(
  MappedReads ~ Batch,
  data = bwameth_mappedreads,
  method = "wilcox.test"
)

bwameth_mappedreads$Batch <- as.factor(bwameth_mappedreads$Batch)

ggplot(bwameth_mappedreads) +
  geom_boxplot(aes(x = Batch, y = MappedReads, fill = Batch)) +
  labs(title = "BWA meth", x = "Sequencing Batch", y = "Mapped Reads") +
  stat_pvalue_manual(bwameth.wilcox.reads, label = "p", y.position = 75000000, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")

### BWA mem
# Shapiro-Wilk normality test
shapiro.test(bwamem_mappedreads$MappedReads)
# data:  bwamem_mappedreads$MappedReads
# W = 0.66366, p-value = 1.426e-07 - p<0.05, need to use non parametric test

# Wilcoxon rank sum exact test
bwamem.wilcox.reads <- wilcox.test(MappedReads ~ Batch, data = bwamem_mappedreads)
bwamem.wilcox.reads
# data:  MappedReads by Batch
# W = 2, p-value = 1.459e-08
# alternative hypothesis: true location shift is not equal to 0

bwamem.wilcox.reads <- compare_means(
  MappedReads ~ Batch,
  data = bwamem_mappedreads,
  method = "wilcox.test"
)

bwamem_mappedreads$Batch <- as.factor(bwamem_mappedreads$Batch)

ggplot(bwamem_mappedreads) +
  geom_boxplot(aes(x = Batch, y = MappedReads, fill = Batch)) +
  labs(title = "BWA mem", x = "Sequencing Batch", y = "Mapped Reads") +
  stat_pvalue_manual(bwamem.wilcox.reads, label = "p", y.position = 36000000, step.increase = 0.2) +
  theme_cowplot() +
  theme(legend.position = "none")


##################################################
################ WGBS ############################
##################################################

MappingEfficiency <- read_excel("QCStats_1.xlsx", sheet = "WBGSBamtoolsStats")

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

ggplot(dat) +
  aes(x = "Alignment Method", y = "PercentMapped", color = "Alignment Method") +
  geom_boxplot() +
  theme(legend.position = "none")

####### Too few WGBS samples to do statistical tests, just report mean and SD 

bwameth <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped) %>% 
  filter(AlignmentMethod == "bwa meth")

mean(bwameth$PercentMapped) 
# 0.9958422
sd(bwameth$PercentMapped)
# 0.001569624

Bismark <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped) %>% 
  filter(AlignmentMethod == "Bismark")

mean(Bismark$PercentMapped)
# 0.387
sd(Bismark$PercentMapped)
# 0.01235584

bwamem <- MappingEfficiency %>%
  select(AlignmentMethod, PercentMapped) %>% 
  filter(AlignmentMethod == "bwa mem")

mean(bwamem$PercentMapped)
# 0.7464473
sd(bwamem$PercentMapped)
# 0.03429545

####### Mean and SD of mapped reads
bwameth <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads) %>% 
  filter(AlignmentMethod == "bwa meth")

mean(bwameth$MappedReads)
# 58,620,964
sd(bwameth$MappedReads)
# 1,107,988

bwamem <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads) %>% 
  filter(AlignmentMethod == "bwa mem")

mean(bwamem$MappedReads)
# 44,017,088
sd(bwamem$MappedReads)
# 2,656,871

Bismark <- MappingEfficiency %>%
  select(AlignmentMethod, MappedReads) %>% 
  filter(AlignmentMethod == "Bismark")

mean(Bismark$MappedReads)
sd(Bismark$MappedReads)