# Asssess how mean percent methylation and depth per sample varies due to read mapping method
# Figures 2 & 3


# Load libraries ----------------------------------------------------------

library(readxl)
library(tidyverse)
library(ggplot2)
library(car)
library(FSA)
library(ggstatsplot)
library(ggpubr)
library(dplyr)
library(ggpmisc)
library(cowplot)
library(DHARMa)


# Depth -------------------------------------------------------------------


setwd("/MethylMethods")


# Mean percent methylation ------------------------------------------------



PercMethyl <- read_excel("QCStats_1.xlsx", sheet = "MeanPercMeth")

PercMethyl <- PercMethyl %>%
  filter(SequencingMethod == "RRBS")

PercMethylClean <- PercMethyl %>% 
  filter(!is.na(bismarkL.mean) & !is.na(bismark.mean)) %>%
  mutate(Bismark.BismarkL.residuals = bismark.mean - bismarkL.mean)

summary <- PercMethylClean %>%
  group_by(AlignmentMethod, SequencingMethod) %>%
  summarise(
    n = n(),
    mean = mean(Mean, na.rm = TRUE),
    sd = sd(Mean, na.rm = TRUE),
    se = (sd(Mean, na.rm = TRUE))/(sqrt(n())),
    min = min(Mean, na.rm = TRUE),
    max = max(Mean, na.rm = TRUE)
  )
summary

## Plot
ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismarkL.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark Local (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bwameth.mean, y = biscuit.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  xlim("0", "100") +
  ylim("0", "100") +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bismark.mean, y = biscuit.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  xlim("0", "100") +
  ylim("0", "100") +
  geom_point()

ggplot(data = PercMethyl, aes(x = bismarkL.mean, y = biscuit.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Biscuit Local (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  xlim("0", "100") +
  ylim("0", "100") +
  geom_point()

ggplot(data = PercMethyl, aes(x = bsbolt.mean, y = biscuit.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Biscuit Local (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  xlim("0", "100") +
  ylim("0", "100") +
  geom_point()

ggplot(data = PercMethyl, aes(x = bismark.mean, y = bismarkL.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Bismark Local (mean % methylation)") +
  theme_cowplot() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  xlim("0", "100") +
  ylim("0", "100") +
  geom_point()


# Model -------------------------------------------------------------------


bismark.bwameth.mean <- glm(bismark.mean ~ bwameth.mean, data = PercMethyl)
summary(bismark.bwameth.mean)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   9.24244    2.10652   4.388 0.000117 ***
# bwameth.mean  0.93516    0.05342  17.507  < 2e-16 ***
# AIC: 194.2
# Calculate R-squared
deviance <- summary(bismark.bwameth.mean)$deviance
null_deviance <- summary(bismark.bwameth.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.9054638

bismarkL.bwameth.mean <- glm(bismarkL.mean ~ bwameth.mean, data = PercMethyl)
summary(bismarkL.bwameth.mean)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  13.09621    2.35468   5.562 4.29e-06 ***
# bwameth.mean  0.85790    0.05967  14.377 2.94e-15 ***
# AIC: 195.9
deviance <- summary(bismarkL.bwameth.mean)$deviance
null_deviance <- summary(bismarkL.bwameth.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.8695874

Biscuit.bwameth.mean <- glm(biscuit.mean ~ bwameth.mean, data = PercMethyl)
summary(Biscuit.bwameth.mean)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  57.78798    1.11916   51.63  < 2e-16 ***
# bwameth.mean  0.37709    0.02838   13.29 1.43e-14 ***
# AIC: 151.19
# Calculate R-squared
deviance <- summary(Biscuit.bwameth.mean)$deviance
null_deviance <- summary(Biscuit.bwameth.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.8465655

Biscuit.bismark.mean <- glm(biscuit.mean ~ bismark.mean, data = PercMethyl)
summary(Biscuit.bismark.mean)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  55.17148    1.42737   38.65  < 2e-16 ***
# bismark.mean  0.37808    0.03111   12.15 1.57e-13 ***
# AIC: 156.25
# Calculate R-squared
deviance <- summary(Biscuit.bismark.mean)$deviance
null_deviance <- summary(Biscuit.bismark.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.8219413

Biscuit.bismarkL.mean <- glm(biscuit.mean ~ bismarkL.mean, data = PercMethyl)
summary(Biscuit.bismarkL.mean)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   53.70333    1.08353   49.56   <2e-16 ***
# bismarkL.mean  0.40834    0.02324   17.57   <2e-16 ***
# AIC: 128.16
# Calculate R-squared
deviance <- summary(Biscuit.bismarkL.mean)$deviance
null_deviance <- summary(Biscuit.bismarkL.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.9087485

Bismark.bismarkL.mean <- glm(bismarkL.mean ~ bismark.mean, data = PercMethylClean)
summary(Bismark.bismarkL.mean)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.52236    1.46215   3.093  0.00417 ** 
# bismark.mean  0.92303    0.03198  28.866  < 2e-16 ***
# AIC: 153.3
# Calculate R-squared
deviance <- summary(Bismark.bismarkL.mean)$deviance
null_deviance <- summary(Bismark.bismarkL.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.9641314

Biscuit.bsbolt.mean <- glm(biscuit.mean ~ bsbolt.mean, data = PercMethyl)
summary(Biscuit.bsbolt.mean)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 14.22742    2.76315   5.149 1.29e-05 ***
# bsbolt.mean  0.79275    0.03786  20.937  < 2e-16 ***
# AIC: 123.54
# Calculate R-squared
deviance <- summary(Biscuit.bsbolt.mean)$deviance
null_deviance <- summary(Biscuit.bsbolt.mean)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
# 0.9319693


# Residuals ---------------------------------------------------------------


# BWA meth vs Bismark Residuals
meth.Bismark_predicted_values <- fitted(bismark.bwameth.mean)
head(meth.Bismark_predicted_values)

meth.Bismark.residuals <- residuals(bismark.bwameth.mean)
head(meth.Bismark.residuals)
hist(meth.Bismark.residuals)

PercMethyl$meth.Bis.Res <- meth.Bismark.residuals

# BWA meth vs Bismark Local Residuals
meth.BismarkL_predicted_values <- fitted(bismarkL.bwameth.mean)
head(meth.BismarkL_predicted_values)

meth.BismarkL.residuals <- residuals(bismarkL.bwameth.mean)
head(meth.BismarkL.residuals)
hist(meth.BismarkL.residuals)

PercMethyl$meth.BisL.Res <- meth.BismarkL.residuals

# BWA meth vs biscuit Residuals
biscuit.meth_predicted_values <- fitted(Biscuit.bwameth.mean)
head(biscuit.meth_predicted_values)

biscuit.meth.residuals <- residuals(Biscuit.bwameth.mean)
head(biscuit.meth.residuals)
hist(biscuit.meth.residuals)

PercMethyl$biscuit.meth.Res <- biscuit.meth.residuals

# BSBolt vs biscuit Residuals
biscuit.bsbolt_predicted_values <- fitted(Biscuit.bsbolt.mean)
head(biscuit.bsbolt_predicted_values)

biscuit.bsbolt.meth.residuals <- residuals(Biscuit.bsbolt.mean)
head(biscuit.bsbolt.meth.residuals)
hist(biscuit.bsbolt.meth.residuals)

PercMethyl$biscuit.bsbolt.Res <- biscuit.bsbolt.meth.residuals

# Bismark vs biscuit Residuals
biscuit.Bismark_predicted_values <- fitted(Biscuit.bwameth.mean)
head(biscuit.Bismark_predicted_values)

biscuit.Bismark.residuals <- residuals(Biscuit.bwameth.mean)
head(biscuit.Bismark.residuals)
hist(biscuit.Bismark.residuals)

PercMethyl$biscuit.Bis.Res <- biscuit.Bismark.residuals

# Bismark vs Bismark Local Residuals
Bismark.BismarkL_predicted_values <- fitted(Bismark.bismarkL.mean)
head(Bismark.BismarkL_predicted_values)

Bismark.BismarkL.residuals <- residuals(Bismark.bismarkL.mean)
head(Bismark.BismarkL.residuals)
hist(Bismark.BismarkL.residuals)

residuals_with_na <- rep(NA, 34)

missing_row <- 12
residuals_with_na[-missing_row] <- Bismark.BismarkL.residuals

PercMethyl$Bis.BisL.Res <- residuals_with_na





# Figure 2 ----------------------------------------------------------------

## Figure 2. Plot to see if each read mapper results in similar mean percent methylation per sample. If they agree, there should be a 1:1 relationship between mean percent methylation. Color by batch (p plots) and sequencing effort (q plots) to see if they had an effect on agreement 

PercMethyl$Batch <- as.factor(PercMethyl$Batch)
PercMethylClean$Batch <- as.factor(PercMethylClean$Batch)


p1<-ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 30, r = 10, b = 10, l = 10))+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(color=Batch), size = 3)+  
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  ylim(c(0,100))+
  xlim(c(0,100))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
p1

q1 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.1) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
q1

q1b <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismarkL.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark Local (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.1) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
q1b

p2 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = biscuit.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none",plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85)+
  geom_point(aes(color=Batch), size = 3)+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
p2

q2 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = biscuit.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.1)+
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
q2

p3 <- ggplot(data = PercMethyl, aes(x = bismark.mean, y = biscuit.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(color=Batch), size=3)+  
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
p3

q3 <- ggplot(data = PercMethyl, aes(x = bismark.mean, y = biscuit.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12))+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.1) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
q3

p3b <- ggplot(data = PercMethyl, aes(x = bismarkL.mean, y = biscuit.mean)) +
  labs(x = "Bismark Local (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(color=Batch), size=3)+  
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
p3b

q3b <- ggplot(data = PercMethyl, aes(x = bismarkL.mean, y = biscuit.mean)) +
  labs(x = "Bismark Local (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12))+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.1) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
q3b

q4 <- ggplot(data = PercMethyl, aes(x = bsbolt.mean, y = biscuit.mean)) +
  labs(x = "BiSulfite Bolt (mean % methylation)", y = "Biscuit (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.1) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
q4

Bis <- ggplot(data = PercMethyl, aes(x = bismark.mean, y = bismarkL.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "Bismark Local (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12))+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.1) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(0,100))+
  xlim(c(0,100))
Bis


q1 + Bis + q2 + q4


# Figure 3: Plot residuals against total sequences per sample --------------------


hist(PercMethyl$TotalSequences)
range(PercMethyl$TotalSequences)


BisR <- ggplot(data = PercMethylClean, aes(x = TotalSequences, y = Bismark.BismarkL.residuals)) +
  labs(x = "Total Sequences", y = "Bismark vs Bismark Local Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color=Batch), size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
BisR

Bis + BisR

p4 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = meth.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "BWA meth vs Bismark Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color=Batch), size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
p4

q4 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = meth.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "BWA meth vs Bismark Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
q4

p5 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = biscuit.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "Bismark vs Biscuit Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color = Batch), size=2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(-3, 3.5)
p5

q5 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = biscuit.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "Bismark vs Biscuit Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position="none")+
  stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 1) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(-3, 3.5)
q5

p6 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = biscuit.meth.residuals)) +
  labs(x = "Total Sequences", y = "BWA meth vs Biscuit Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color = Batch), size=2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(-2.5, 3.5))
p6

q6 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = biscuit.meth.residuals)) +
  labs(x = "Total Sequences", y = "BWA meth vs Biscuit Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(-2.5, 3.5))
q6

BisR <- ggplot(data = PercMethylClean, aes(x = TotalSequences, y = Bismark.BismarkL.residuals)) +
  labs(x = "Total Sequences", y = "Bismark vs Bismark Local Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -1.5, vjust = 1, size = 6, fontface = "bold")
BisR

BiscBSBR <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = biscuit.bsbolt.meth.residuals)) +
  labs(x = "Total Sequences", y = "BiSulfite Bolt vs Biscuit Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position="none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("<10M","10M – 20M",">20M"))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
BiscBSBR

q4 + BisR + q6 + BiscBSBR

p1 + p2 + p3 + plot_layout(ncol = 3)
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, nrow = 2)

q1 + q2 + q3 + plot_layout(ncol = 3)
q1 + q2 + q3 + q4 + q5 + q6 + Bis + BisR + plot_layout(ncol = 3)



# Model residuals as exponential decay ------------------------------------

## Bismark vs BWA meth residuals

MethBisExp <- glm(abs(meth.Bis.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(MethBisExp)
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.545e+00  2.210e-01   6.989 6.41e-08 ***
# TotalSequences -3.028e-08  9.891e-09  -3.061  0.00444 ** 
# 
# (Dispersion parameter for Gamma family taken to be 0.5945129)
# 
# Null deviance: 23.561  on 33  degrees of freedom
# Residual deviance: 19.636  on 32  degrees of freedom
# AIC: 135.6

simulationOutput <- simulateResiduals(fittedModel = MethBisExp)
plot(simulationOutput)

testDispersion(simulationOutput)

MethBisLinear <- glm(abs(meth.Bis.Res)~TotalSequences, data = PercMethyl)
summary(MethBisLinear)
# AIC: 163.14


# Bismark vs Bismark Local residuals

BisBisLExp <- glm(abs(Bis.BisL.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(BisBisLExp)
#                 Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     6.107e-01  5.184e-01   1.178   0.2478  
# TotalSequences -4.818e-08  2.312e-08  -2.084   0.0455 *
#
# (Dispersion parameter for Gamma family taken to be 3.246341)
# 
# Null deviance: 71.419  on 32  degrees of freedom
# Residual deviance: 58.700  on 31  degrees of freedom
# (1 observation deleted due to missingness)
# AIC: 52.612

simulationOutput <- simulateResiduals(fittedModel = BisBisLExp)
plot(simulationOutput)

testDispersion(simulationOutput)

BisBisLLinear <- glm(abs(Bis.BisL.Res)~TotalSequences, data = PercMethyl)
summary(BisBisLLinear)
# AIC: 145.37


# BWA meth vs Biscuit residuals

MethBiscExp <- glm(abs(biscuit.meth.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(MethBiscExp)
#                    Estimate Std. Error t value Pr(>|t|)  
# (Intercept)     5.989e-01  3.444e-01   1.739   0.0916 .
# TotalSequences -1.806e-08  1.541e-08  -1.172   0.2499  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for Gamma family taken to be 1.442961)
# 
# Null deviance: 30.863  on 33  degrees of freedom
# Residual deviance: 29.297  on 32  degrees of freedom
# AIC: 91.642

simulationOutput <- simulateResiduals(fittedModel = MethBiscExp)
plot(simulationOutput)

testDispersion(simulationOutput)

MethBiscLinear <- glm(abs(biscuit.meth.Res)~TotalSequences, data = PercMethyl)
summary(MethBiscLinear)
# AIC: 131.06


# BSBolt vs Biscuit residuals

BSBBiscExp <- glm(abs(biscuit.bsbolt.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(BSBBiscExp)
#                   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)     5.643e-01  2.470e-01   2.284  0.02913 * 
# TotalSequences -3.450e-08  1.105e-08  -3.121  0.00381 **
#   
# (Dispersion parameter for Gamma family taken to be 0.7425881)
# 
# Null deviance: 30.324  on 33  degrees of freedom
# Residual deviance: 25.804  on 32  degrees of freedom
# AIC: 67.899

simulationOutput <- simulateResiduals(fittedModel = BSBBiscExp)
plot(simulationOutput)

testDispersion(simulationOutput)

BSBBiscLinear <- glm(abs(biscuit.bsbolt.Res)~TotalSequences, data = PercMethyl)
summary(BSBBiscLinear)
# AIC: 93.301
