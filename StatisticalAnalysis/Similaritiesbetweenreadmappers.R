#### Test for differrences in mean percent methylation between read mapping methods

# load libraries and import data

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


PercMethyl <- read_excel("QCStats_1.xlsx", sheet = "MeanPercMeth")

PercMethyl <- PercMethyl %>%
  filter(SequencingMethod == "RRBS")

## Plot
ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bwamem.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bismark.mean, y = bwamem.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  geom_point()

## Model - should be the same as the model in the figures

bismark.bwameth.mean <- glm(bismark.mean ~ bwameth.mean, data = PercMethyl)
summary(bismark.bwameth.mean)

bwamem.bwameth.mean <- glm(bwamem.mean ~ bwameth.mean, data = PercMethyl)
summary(bwamem.bwameth.mean)

bwamem.bismark.mean <- glm(bwamem.mean ~ bismark.mean, data = PercMethyl)
summary(bwamem.bismark.mean)
