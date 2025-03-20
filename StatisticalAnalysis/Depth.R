# Depth of RRBS and WGBS

library(readxl)
library(tidyverse)
library(ggplot2)
library(car)
library(FSA)
library(ggstatsplot)
library(ggpubr)
library(dplyr)
library(cowplot)
library(ggpmisc)

#########################################
############## Depth ####################
#########################################

Depth <- read_excel("/PATH/QCStats_1.xlsx", sheet = "MethylKitDepth")
hist(Depth$Mean)
hist(Depth$Max)

###################### RRBS

#### BWA meth
bwameth <- Depth %>%
  filter(AlignmentMethod == "bwameth") %>% 
  filter(SequencingMethod == "RRBS")

mean(bwameth$Mean) 
# 16.39324
sd(bwameth$Mean)
# 3.455886
sd(bwameth$Mean)/sqrt(length((bwameth))) #standard error
# 0.9236245

mean(bwameth$Max)
# 18894.88

# bwa meth batch 1 & batch 2
BwaMeth_1 <- Depth %>%
  filter(AlignmentMethod == "bwameth") %>% 
  filter(SequencingMethod == "RRBS") %>% 
  filter(Batch == "1")

mean(BwaMeth_1$Mean)
# 14.72909
sd(BwaMeth_1$Mean)
# 1.095619
sd(BwaMeth_1$Mean)/sqrt(length((BwaMeth_1))) #standard error
# 0.2928164

BwaMeth_2 <- Depth %>%
  filter(AlignmentMethod == "bwameth") %>% 
  filter(SequencingMethod == "RRBS") %>% 
  filter(Batch == "2")

mean(BwaMeth_2$Mean)
# 19.44417
sd(BwaMeth_2$Mean)
# 4.224297
sd(BwaMeth_2$Mean)/sqrt(length((BwaMeth_2))) #standard error
# 1.128991

# Shapiro-Wilk normality test
shapiro.test(bwameth$Mean)
# data:  bwameth$Mean
# W = 0.73439, p-value = 1.734e-06

# Wilcoxon rank sum exact test
wilcox.test(Mean~Batch, data = bwameth) 
# data:  Mean by Batch
# W = 24.5, p-value = 0.0001151
# alternative hypothesis: true location shift is not equal to 0
# 
# Warning message:
#   In wilcox.test.default(x = DATA[[1L]], y = DATA[[2L]], ...) :
#   cannot compute exact p-value with ties

BWAMeth.wilcox.test <- compare_means(
  Mean~Batch,
  data = bwameth,
  method = "wilcox.test"
)

bwameth$Batch <- as.factor(bwameth$Batch)

ggplot(bwameth) +
  geom_boxplot(aes(x = Batch, y = Mean, fill = Batch)) +
  labs(title = "BWA Meth", x = "Sequencing Batch", y = "Mean Depth") +
  stat_pvalue_manual(BWAMeth.wilcox.test, label = "p", y.position = 30, step.increase = 0.001) +
  theme_cowplot() +
  theme(legend.position = "top")

###################### WGBS

# bwa meth
WGBSbwameth <- Depth %>%
  filter(AlignmentMethod == "bwameth") %>% 
  filter(SequencingMethod == "WGBS")

mean(WGBSbwameth$Mean) 
# 12.7325
sd(WGBSbwameth$Mean)
# 0.06238322
sd(WGBSbwameth$Mean)/sqrt(length((WGBSbwameth))) #standard error
# 0.01667262

bwameth_all <- Depth %>%
     filter(AlignmentMethod == "bwameth")

ggplot(Depth) +
  geom_boxplot(aes(x = SequencingMethod, y = Mean, colour = AlignmentMethod)) +
  labs(x = "Sequencing Method", y = "Mean Depth") +
  theme_cowplot() +
  theme(legend.position = "top")