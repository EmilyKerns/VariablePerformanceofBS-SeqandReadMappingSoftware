# Analyze depth of RRBS and WGBS


# import data and load packages

library(readxl)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpmisc)
library(patchwork)
library(DHARMa)


PercMethyl <- read_excel("PATH/Metadata/QC/QCStats_1.xlsx", sheet = "MeanPercMeth")

PercMethyl <- PercMethyl %>%
  filter(SequencingMethod == "RRBS")

## Figure 2A-C. Plot to see if each read mapper results in similar mean percent methylation per sample. If they agree, there should be a 1:1 relationship between mean percent methylation. Color by batch (p plots) and sequencing effort (q plots) to see if they had an effect on agreement 

PercMethyl$Batch <- as.factor(PercMethyl$Batch)


p1<-ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 30, r = 10, b = 10, l = 10))+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(color=Batch), size = 3)+  
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
p1

q1 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "Bismark (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 3) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
q1

p2 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bwamem.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none",plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85)+
  geom_point(aes(color=Batch), size = 3)+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(90,100))
p2

q2 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bwamem.mean)) +
  labs(x = "BWA meth (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.85)+
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 3) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(90,100))
q2

p3 <- ggplot(data = PercMethyl, aes(x = bismark.mean, y = bwamem.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12), plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(color=Batch), size=3)+  
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(90,100))
p3

q3 <- ggplot(data = PercMethyl, aes(x = bismark.mean, y = bwamem.mean)) +
  labs(x = "Bismark (mean % methylation)", y = "BWA mem (mean % methylation)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.text=element_text(size=12))+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.85) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 3) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(90,100))
q3


## Linear model the mean percent methylation between each sample

bismark.bwameth.mean <- glm(bismark.mean ~ bwameth.mean, data = PercMethyl)
summary(bismark.bwameth.mean)

bwamem.bwameth.mean <- glm(bwamem.mean ~ bwameth.mean, data = PercMethyl)
summary(bwamem.bwameth.mean)

bwamem.bismark.mean <- glm(bwamem.mean ~ bismark.mean, data = PercMethyl)
summary(bwamem.bismark.mean)

## The relationship between BWA mem and the other read mappers doesn't appear to be linear. There may be differences in agreement of percent methylation due to sequencing effort for each sample. Compare the amount of sequencing effort for each sample to the residuals for the models (q plots). Also look at batch (p plots).

# BWA meth vs Bismark Residuals
meth.Bismark_predicted_values <- fitted(bismark.bwameth.mean)
head(meth.Bismark_predicted_values)

meth.Bismark.residuals <- residuals(bismark.bwameth.mean)
head(meth.Bismark.residuals)
hist(meth.Bismark.residuals)

PercMethyl$meth.Bis.Res <- meth.Bismark.residuals

# BWA meth vs BWA mem Residuals
mem.meth_predicted_values <- fitted(bwamem.bwameth.mean)
head(mem.meth_predicted_values)

mem.meth.residuals <- residuals(bwamem.bwameth.mean)
head(mem.meth.residuals)
hist(mem.meth.residuals)

PercMethyl$mem.meth.Res <- mem.meth.residuals

# Bismark vs BWA mem Residuals
mem.Bismark_predicted_values <- fitted(bwamem.bismark.mean)
head(mem.Bismark_predicted_values)

mem.Bismark.residuals <- residuals(bwamem.bismark.mean)
head(mem.Bismark.residuals)
hist(mem.Bismark.residuals)

PercMethyl$mem.Bis.Res <- mem.Bismark.residuals

## Plot residuals against total sequences per sample

hist(PercMethyl$TotalSequences)
range(PercMethyl$TotalSequences)

p4 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = meth.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "BWA meth vs Bismark Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color=Batch), size = 2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
p4

q4 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = meth.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "BWA meth vs Bismark Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  #stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")
q4

p5 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = mem.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "Bismark vs BWA mem Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color = Batch), size=2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(-3, 3.5)
p5

q5 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = mem.Bismark.residuals)) +
  labs(x = "Total Sequences", y = "Bismark vs BWA mem Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position="none")+
  #stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(-3, 3.5)
q5

p6 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = mem.meth.residuals)) +
  labs(x = "Total Sequences", y = "Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 10, l = 10) )+
  stat_poly_line() +
  geom_point(aes(color = Batch), size=2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  scale_color_manual(values =c("#56B1F7","#132B43"))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(-2.5, 3.5))
p6

q6 <- ggplot(data = PercMethyl, aes(x = TotalSequences, y = mem.meth.residuals)) +
  labs(x = "Total Sequences", y = "Residuals") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  #stat_poly_line() +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(-2.5, 3.5))
q6

p1 + p2 + p3 + plot_layout(ncol = 3)
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3, nrow = 2)

q1 + q2 + q3 + plot_layout(ncol = 3)
q1 + q2 + q3 + q4 + q5 + q6 + plot_layout(ncol = 3)



#### Find difference in mean percent methylation when samples have less than 20 million reads

PercMethyl$SeqCat <-ifelse(PercMethyl$TotalSequences<=10000000,"Low",
                           ifelse(PercMethyl$TotalSequences>=20000000,"High",
                                  ifelse(PercMethyl$TotalSequences>10000000&PercMethyl$TotalSequences<20000000,"Inter", "WHOOPS")))

PercMethyl$MethMem <- PercMethyl$bwameth.mean - PercMethyl$bwamem.mean
PercMethyl$BisMem <- PercMethyl$bismark.mean - PercMethyl$bwamem.mean
PercMethyl$MethBis <- PercMethyl$bwameth.mean - PercMethyl$bismark.mean

mean(abs(PercMethyl$MethBis))
#7.302912
max(abs(PercMethyl$MethBis))
# 12.671
min(abs(PercMethyl$MethBis))
# 0.202

mean(abs(PercMethyl$MethMem))
# 56.76697
mean(abs(PercMethyl$BisMem))
# 49.94418

meanDiff = PercMethyl %>% group_by(SeqCat)  %>%
  summarise(MethMem = mean(abs(MethMem)), 
            BismMem = mean(abs(BisMem)), 
            MethBis = mean(abs(MethBis)), 
            .groups = 'drop')
meanDiff

maxDiff = PercMethyl %>% group_by(SeqCat)  %>%
  summarise(MethMem = max(abs(MethMem)), 
            BismMem = max(abs(BisMem)), 
            MethBis = max(abs(MethBis)), 
            .groups = 'drop')
maxDiff

minDiff = PercMethyl %>% group_by(SeqCat)  %>%
  summarise(MethMem = min(abs(MethMem)), 
            BismMem = min(abs(BisMem)), 
            MethBis = min(abs(MethBis)), 
            .groups = 'drop')
minDiff

#### Test for improved precision of percent methylation above 20M reads per sample for BWA meth vs Bismark

MethBisExp <- glm(abs(meth.Bis.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(MethBisExp)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
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
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     3.881e+00  7.210e-01   5.382 6.53e-06 ***
# TotalSequences -5.546e-08  3.226e-08  -1.719   0.0953 .
#
# (Dispersion parameter for gaussian family taken to be 6.324918)
# 
# Null deviance: 221.09  on 33  degrees of freedom
# Residual deviance: 202.40  on 32  degrees of freedom
# AIC: 163.14

simulationOutput <- simulateResiduals(fittedModel = MethBisLinear)
plot(simulationOutput)

testDispersion(simulationOutput)

MethMemExp <- glm(abs(mem.meth.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(MethMemExp)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)
# (Intercept)    -9.009e-02  2.585e-01  -0.348    0.730
# TotalSequences -8.130e-09  1.157e-08  -0.703    0.487
# 
# (Dispersion parameter for Gamma family taken to be 0.8133221)
# 
# Null deviance: 30.710  on 33  degrees of freedom
# Residual deviance: 30.148  on 32  degrees of freedom
# AIC: 57.14

simulationOutput <- simulateResiduals(fittedModel = MethMemExp)
plot(simulationOutput)

testDispersion(simulationOutput)

MethMemLin <- glm(abs(mem.meth.Res)~TotalSequences, data = PercMethyl)
summary(MethMemLin)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     9.510e-01  2.098e-01   4.532 7.69e-05 ***
#   TotalSequences -8.622e-09  9.389e-09  -0.918    0.365    
# 
# (Dispersion parameter for gaussian family taken to be 0.5357192)
# 
# Null deviance: 17.595  on 33  degrees of freedom
# Residual deviance: 17.143  on 32  degrees of freedom
# AIC: 79.206

simulationOutput <- simulateResiduals(fittedModel = MethMemLin)
plot(simulationOutput)

testDispersion(simulationOutput)

BisMemExp <- glm(abs(mem.Bis.Res)~TotalSequences, data = PercMethyl, family = Gamma(link="log"))
summary(BisMemExp)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)     8.176e-02  2.053e-01   0.398    0.693
# TotalSequences -1.106e-08  9.187e-09  -1.204    0.237
# 
# (Dispersion parameter for Gamma family taken to be 0.5128677)
# 
# Null deviance: 21.563  on 33  degrees of freedom
# Residual deviance: 20.618  on 32  degrees of freedom
# AIC: 60.349

simulationOutput <- simulateResiduals(fittedModel = BisMemExp)
plot(simulationOutput)

testDispersion(simulationOutput)

BisMemLin <- glm(abs(mem.Bis.Res)~TotalSequences, data = PercMethyl)
summary(BisMemLin)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     1.116e+00  1.902e-01   5.865 1.61e-06 ***
#   TotalSequences -1.190e-08  8.513e-09  -1.398    0.172    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 0.4403583)
# 
# Null deviance: 14.952  on 33  degrees of freedom
# Residual deviance: 14.091  on 32  degrees of freedom
# AIC: 72.541

simulationOutput <- simulateResiduals(fittedModel = BisMemLin)
plot(simulationOutput)

testDispersion(simulationOutput)

# Add GLM model to the residual plots
newdata <- data.frame(TotalSequences = seq(min(PercMethyl$TotalSequences), 
                                           max(PercMethyl$TotalSequences), 
                                           length.out = 100))


# Make the figure

library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)


pdf(file = "PATH/Figure2.pdf",
    width = 6.6, 
    height = 7)

# Modify the figure
q1 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bismark.mean)) +
  labs(x = "BWA meth (% meth)", y = "Bismark (% meth)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), label.x = 0.05, label.y = 0.85) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")

q2 <- ggplot(data = PercMethyl, aes(x = bwameth.mean, y = bwamem.mean)) +
  labs(x = "BWA meth (% meth)", y = "BWA mem (% meth)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.85)+
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(90,100))

q3 <- ggplot(data = PercMethyl, aes(x = bismark.mean, y = bwamem.mean)) +
  labs(x = "Bismark (% meth)", y = "BWA mem (% meth)") +
  theme_classic() +
  theme(text = element_text(size = 12), legend.position = "none")+
  stat_poly_line() +
  stat_poly_eq(use_label("eq"), , label.x = 0.05, label.y = 0.85) +
  geom_point(aes(colour = cut(TotalSequences, c(-Inf,10000000,20000000,Inf))),
             size = 2) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 6, fontface = "bold")+
  ylim(c(90,100))

# Predict values using the model
newdata$predicted <- predict(MethBisExp, newdata, type = "response")

q4 <- q4 + 
  # Add upper curve (positive residuals)
  geom_line(data = newdata, aes(x = TotalSequences, y = predicted),
            color = "firebrick", linetype = "solid", size = 1) +
  # Add lower curve (negative residuals)
  geom_line(data = newdata, aes(x = TotalSequences, y = -predicted), 
            color = "firebrick", linetype = "solid", size = 1) +
  labs(x = "Total Sequences", y = "Residuals") +
  theme(legend.position = "none") +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))

# Predict GLM values using the model
newdata$predicted <- predict(BisMemExp, newdata, type = "response")

# Add the GLM curve
q5 <- q5 + 
  # Add upper curve (positive residuals)
  geom_line(data = newdata, aes(x = TotalSequences, y = predicted), 
            color = "firebrick", linetype = "solid", size = 1) +
  # Add lower curve (negative residuals)
  geom_line(data = newdata, aes(x = TotalSequences, y = -predicted), 
            color = "firebrick", linetype = "solid", size = 1) +
  labs(x = "Total Sequences", y = "Residuals") +
  theme(legend.text=element_text(size=12)) +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))

# Predict GLM values using the model
newdata$predicted <- predict(MethMemExp, newdata, type = "response")

# Add the GLM curve
q6 <- q6 + 
  # Add upper curve (positive residuals)
  geom_line(data = newdata, aes(x = TotalSequences, y = predicted), 
            color = "firebrick", linetype = "solid", size = 1) +
  labs(x = "Total Sequences", y = "Residuals") +
  # Add lower curve (negative residuals)
  geom_line(data = newdata, aes(x = TotalSequences, y = -predicted), 
            color = "firebrick", linetype = "solid", size = 1) +
  # Update the legend to include the GLM line
  theme(legend.position = "none") +
  scale_color_manual(name = "Total Sequences",
                     values = c(
                       "(-Inf,1e+07]" = "#56B1F7",
                       "(1e+07,2e+07]" = "#3670A0",
                       "(2e+07, Inf]" = "#132B43"),
                     labels = c("Minimum – 10M","10M – 20M","20M – Maximum"))

q1 + q4 + q2 + q5 + q3 + q6 + plot_layout(ncol = 2)

dev.off()

################# Plot median

ggplot(data = PercMethyl, aes(x = bwameth.median, y = bismark.median)) +
  labs(x = "BWA meth (median % methylation)", y = "Bismark (median % methylation)") +
  theme_classic() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bwameth.median, y = bwamem.median)) +
  labs(x = "BWA meth (median % methylation)", y = "BWA mem (median % methylation)") +
  theme_classic() +
  stat_poly_line() +
  xlim("0", "100") +
  ylim("0", "100") +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

ggplot(data = PercMethyl, aes(x = bismark.median, y = bwamem.median)) +
  labs(x = "Bismark (median % methylation)", y = "BWA mem (median % methylation)") +
  theme_classic() +
  stat_poly_line() +
  stat_poly_eq(use_label("eq")) +
  stat_poly_eq(label.y = 0.9) +
  geom_point()

## Model

bismark.bwameth.median <- glm(bismark.median ~ bwameth.median, data = PercMethyl)
summary(bismark.bwameth.median)
# Calculate R-squared
deviance <- summary(bismark.bwameth.median)$deviance
null_deviance <- summary(bismark.bwameth.median)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared)
