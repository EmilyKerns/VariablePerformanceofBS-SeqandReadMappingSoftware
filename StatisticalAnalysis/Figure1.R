#### Test for differrences in mapping efficiency and number of mapped reads between Bismark, BWA meth, BWA mem, and Bisulfite Bolt


# import data

library(readxl)
library(tidyverse)
library(ggplot2)
library(car)
library(FSA)
library(ggstatsplot)
library(ggpubr)
library(dplyr)
library(patchwork)


# Create Figure 1 ---------------------------------------------------------

setwd("/MethylMethods")

MappingEfficiency <- read_excel("QCStats_1.xlsx", sheet = "BamtoolsStats")

MappingEfficiency <- filter(MappingEfficiency, !AlignmentMethod == "bwa mem")

MappingEfficiency$PercentFailedQC[is.na(MappingEfficiency$PercentFailedQC)] <- 0

MappingEfficiency$PercMappedErr <- MappingEfficiency$PercentMapped - MappingEfficiency$PercentFailedQC

MappingEfficiency$AlignmentMethod <- factor(MappingEfficiency$AlignmentMethod, levels = c("Bismark","BismarkLocal","BisulfiteBolt", "Biscuit", "BWA meth"))


# Ensure we're only including samples mapped with all 5 methods

filtered_ME <- MappingEfficiency %>%
  filter(!is.na(AlignmentMethod) & !is.na(PercentMapped)) %>%
  group_by(Organism, SampleID) %>%
  filter(n_distinct(AlignmentMethod) == 5) %>%
  ungroup()

filtered_ME %>%
  filter(Organism == "Stickleback") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )

filtered_ME %>%
  filter(Organism == "Cichlid") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )

filtered_ME %>%
  filter(Organism == "Great Tit") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )

filtered_ME %>%
  filter(Organism == "Sea Urchin") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )

RRBS <- ggplot(filtered_ME) +
  aes(x = AlignmentMethod, y = PercMappedErr, color = Organism) +
  geom_boxplot()+
  scale_color_manual(values=c('#88CCEE','#44AA99', '#117733','#332288'))+
  ylab("Mapping Efficiency")+
  xlab("Alignment Method")+
  ylim(c(0,1))+
  geom_jitter(aes(color = Organism, shape = Layout), position = position_jitterdodge(jitter.width = 0.15, dodge.width = .75), size = 3, alpha = 0.2) +
  scale_shape_manual(values = c(17, 18))+
  theme_classic()+
  theme(legend.position = "bottom", legend.spacing = unit(.25, 'cm'), text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
RRBS


WGBSMappingEfficiency <- read_excel("QCStats_1.xlsx", sheet = "WBGSBamtoolsStats")

WGBSMappingEfficiency <- filter(WGBSMappingEfficiency, !AlignmentMethod == "bwa mem")

WGBSMappingEfficiency$PercentFailedQC[is.na(WGBSMappingEfficiency$PercentFailedQC)] <- 0

WGBSMappingEfficiency$PercMappedErr <- WGBSMappingEfficiency$PercentMapped - WGBSMappingEfficiency$PercentFailedQC

WGBSMappingEfficiency$AlignmentMethod <- factor(WGBSMappingEfficiency$AlignmentMethod, levels = c("Bismark","BismarkLocal","BisulfiteBolt","Biscuit", "bwa meth"))

# Ensure we're only including samples mapped with all 5 methods

filtered_WGBS <- WGBSMappingEfficiency %>%
  filter(!is.na(AlignmentMethod) & !is.na(PercentMapped)) %>%
  group_by(Organism, SampleID) %>%
  filter(n_distinct(AlignmentMethod) == 5) %>%
  ungroup()

filtered_WGBS %>%
  filter(Organism == "Stickleback") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )

filtered_WGBS %>%
  filter(Organism == "Cichlid") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )

filtered_WGBS %>%
  filter(Organism == "Coral") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )


filtered_WGBS %>%
  filter(Organism == "StickInsect") %>%
  group_by(AlignmentMethod) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE)
  )


WGBS <- ggplot(WGBSMappingEfficiency) +
  aes(x = AlignmentMethod, y = PercMappedErr, color = Organism) +
  geom_boxplot()+
  scale_color_manual(values=c('#88CCEE','#0AAB3F','#37F375','#332288'))+
  scale_x_discrete(labels = c("Bismark", "BismarkLocal", "BisulfiteBolt", "Biscuit", "BWA meth")) +
  xlab("Alignment Method") +
  ylim(c(0,1))+
  geom_jitter(aes(color = Organism), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 3, alpha = 0.2) +
  theme_classic() +
  theme(axis.title.y=element_blank(), legend.position = "bottom", legend.spacing = unit(.25, 'cm'), legend.title = element_text(size = 14), text = element_text(size = 14), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size=14, face="bold", colour = "black")) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
WGBS

RRBS + WGBS


pdf(file = "/Figures/Figure1.pdf",
    width = 8.2, 
    height = 6)

RRBS + WGBS 

dev.off()

#   mean in group RRBS mean in group WGBS 
# 22.99964           12.73250 
