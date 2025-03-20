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


PercMethyl <- read_excel("R:/Genohub_seq/methylation_methods_comparison/Metadata/QC/QCStats_1.xlsx", sheet = "MeanPercMeth")

PercMethyl <- PercMethyl %>%
  filter(SequencingMethod == "RRBS")

## Plot
ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bwamem.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bismark.mean, y = bwamem.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

## Model - should be the same as the model in the figures

bismark.bwameth.mean <- glm(bismark.mean ~ bwameth.mean, data = PercMethyl)
summary(bismark.bwameth.mean)
# Calculate R-squared
deviance <- summary(bismark.bwameth.mean)$deviance
null_deviance <- summary(bismark.bwameth.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.9686299

bwamem.bwameth.mean <- glm(bwamem.mean ~ bwameth.mean, data = PercMethyl)
summary(bwamem.bwameth.mean)
# Calculate R-squared
deviance <- summary(bwamem.bwameth.mean)$deviance
null_deviance <- summary(bwamem.bwameth.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.6659768

bwamem.bismark.mean <- glm(bwamem.mean ~ bismark.mean, data = PercMethyl)
summary(bwamem.bismark.mean)
# Calculate R-squared
deviance <- summary(bwamem.bismark.mean)$deviance
null_deviance <- summary(bwamem.bismark.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.6260214