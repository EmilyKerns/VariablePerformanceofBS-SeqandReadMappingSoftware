### Assess overlapping CpG sites between sequencing methods within WT


# Load libraries ----------------------------------------------------------


{
library(methylKit)
library(ggvenn)
library(ggplot2)
library(patchwork)
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicPlot)
}

setwd("/Stickle/NoSNPFilters")



# Import and format data --------------------------------------------------

{
# import samples
WGBSfile.list=list("WGBS_CG_22_WT_002_CpG.methylKit", "WGBS_CG_22_WT_004_CpG.methylKit", "CG_22_WT_002_CpG.methylKit", "CG_22_WT_004_CpG.methylKit")

WGBSvsRRBS=methRead(WGBSfile.list,
                    sample.id=list("WGBS_CG_WT_002", "WGBS_CG_WT_004", "RRBS_CG_WT_002", "RRBS_CG_WT_004"),
                    pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                    assembly="StickleGeneAnnotations",
                    treatment=c(1,1,2,2),
                    context="CpG",
                    mincov = 5
)
}

{
  
  WGBSvsRRBS=filterByCoverage(WGBSvsRRBS,lo.count=5,lo.perc=NULL,
                                         hi.count=NULL,hi.perc=99.9)
# Separate each sample into individual methylRaw objects - contains info about mehtylation and chr location
WGBS_CG_WT_002=WGBSvsRRBS[[1]]
WGBS_CG_WT_004=WGBSvsRRBS[[2]]
RRBS_CG_WT_002=WGBSvsRRBS[[3]]
RRBS_CG_WT_004=WGBSvsRRBS[[4]]

WGBS_CG_WT_002$WGBS_PercMeth <- ((WGBS_CG_WT_002$numCs/WGBS_CG_WT_002$coverage)*100)
WGBS_CG_WT_004$WGBS_PercMeth <- ((WGBS_CG_WT_004$numCs/WGBS_CG_WT_004$coverage)*100)
RRBS_CG_WT_002$RRBS_PercMeth <- ((RRBS_CG_WT_002$numCs/RRBS_CG_WT_002$coverage)*100)
RRBS_CG_WT_004$RRBS_PercMeth <- ((RRBS_CG_WT_004$numCs/RRBS_CG_WT_004$coverage)*100)

WGBS_CG_WT_002 <- as.data.frame(WGBS_CG_WT_002)
WGBS_CG_WT_004 <- as.data.frame(WGBS_CG_WT_004)
RRBS_CG_WT_002 <- as.data.frame(RRBS_CG_WT_002)
RRBS_CG_WT_004 <- as.data.frame(RRBS_CG_WT_004)

WGBS_CG_WT_002 <- methylKit::getData(WGBS_CG_WT_002)
WGBS_CG_WT_004 <- methylKit::getData(WGBS_CG_WT_004)
RRBS_CG_WT_002 <- methylKit::getData(RRBS_CG_WT_002)
RRBS_CG_WT_004 <- methylKit::getData(RRBS_CG_WT_004)

## Overlapping CpG sites on biological replicates with RRBS and WGBS -------


WTWGBS_merged <- merge(WGBS_CG_WT_002, WGBS_CG_WT_004, by = "chr")
nrow(WTWGBS_merged)
# 13,130,739 CpG sites

WTRRBS_merged <- merge(RRBS_CG_WT_002, RRBS_CG_WT_004, by = "chr")
nrow(WTRRBS_merged)
# 255,991 CpG sites

WTWGBS_merged$MeanMeth <- rowMeans(WTWGBS_merged[, c("WGBS_PercMeth.x", "WGBS_PercMeth.y")], na.rm = TRUE)

WTRRBS_merged$MeanMeth <- rowMeans(WTRRBS_merged[, c("RRBS_PercMeth.x", "RRBS_PercMeth.y")], na.rm = TRUE)


## Overlapping CpG sites with RRBS and WGBS on technical replicates --------

WT_002_merged <- merge(WGBS_CG_WT_002, RRBS_CG_WT_002, by = "chr")
hist(WT_002_merged$WGBS_PercMeth)
hist(WT_002_merged$RRBS_PercMeth)


WT_004_merged <- merge(WGBS_CG_WT_004, RRBS_CG_WT_004, by = "chr")
hist(WT_004_merged$WGBS_PercMeth)
hist(WT_004_merged$RRBS_PercMeth)

}

WT_002_merged$Diff <- WT_002_merged$WGBS_PercMeth - WT_002_merged$RRBS_PercMeth
WT_004_merged$Diff <- WT_004_merged$WGBS_PercMeth - WT_004_merged$RRBS_PercMeth

WT_002_merged$SampleID <- "WT_02"
WT_004_merged$SampleID <- "WT_04"

AllWT <- rbind(WT_002_merged, WT_004_merged)

#### Figure 7

ggplot(AllWT, aes(x=coverage.x, y=Diff)) +
  geom_point()+
  theme_classic()+
  labs(x = "WGBS Depth", y = "Difference in Percent Methylation per Base")

ggplot(AllWT, aes(x=coverage.y, y=Diff)) +
  geom_point() +
  theme_classic()+
  labs(x = "RRBS Depth", y = "Difference in Percent Methylation per Base")

ggplot(data = AllWT, aes(x = Diff, color = SampleID)) +
  geom_density(alpha = 0.4) +  
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in Percent Methylation per Base", y = "Frequency")

# Need to run WKRRBSvsWGBS.R to create All dataframe 
proportion <- mean(All$Diff <= 10 & All$Diff >= -10)
proportion
#0.5922322
proportion2 <- mean(All$Diff > 10)
proportion2 #0.2070207
proportion3 <- mean(All$Diff < -10)
proportion3 #0.2007471

test <- glm(RRBS_PercMeth ~ WGBS_PercMeth, data = All)
summary(test)
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.2984303  0.0213267   201.6   <2e-16 ***
# WGBS_PercMeth 0.8854360  0.0004085  2167.4   <2e-16 ***
# 
# (Dispersion parameter for gaussian family taken to be 444.0301)
# 
# Null deviance: 2838053896  on 1693931  degrees of freedom
# Residual deviance:  752155961  on 1693930  degrees of freedom
# AIC: 15133202

n <- nrow(All)
p <- length(coef(test))
r_squared <- 1 - test$deviance / test$null.deviance
adjusted_r_squared <- 1 - (n - 1) / (n - p) * (1 - r_squared)
print(adjusted_r_squared) # 0.7349746

test2 <- glm(abs(Diff) ~ coverage.x, data = All)
summary(test2)
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 14.737956   0.037446   393.6   <2e-16 ***
# coverage.x  -0.272718   0.004126   -66.1   <2e-16 ***
# 
# (Dispersion parameter for gaussian family taken to be 309.515)
# 
# Null deviance: 525649164  on 1693931  degrees of freedom
# Residual deviance: 524296824  on 1693930  degrees of freedom
# AIC: 14521886

n <- nrow(All)
p <- length(coef(test2))
r_squared <- 1 - test2$deviance / test2$null.deviance
adjusted_r_squared <- 1 - (n - 1) / (n - p) * (1 - r_squared)
print(adjusted_r_squared) # 0.002572115

test3 <- glm(abs(Diff) ~ coverage.y, data = All)
summary(test3)
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 15.427461   0.033048  466.82   <2e-16 ***
# coverage.y  -0.381653   0.003841  -99.38   <2e-16 ***
# 
# (Dispersion parameter for gaussian family taken to be 308.5148)
# 
# Null deviance: 525649164  on 1693931  degrees of freedom
# Residual deviance: 522602453  on 1693930  degrees of freedom
# AIC: 14516403

n <- nrow(All)
p <- length(coef(test3))
r_squared <- 1 - test3$deviance / test3$null.deviance
adjusted_r_squared <- 1 - (n - 1) / (n - p) * (1 - r_squared)
print(adjusted_r_squared) #0.005795505


fill_colors <- c("WK_02" = "skyblue", "WK_03"= "darkblue", "WT_02"= "purple", "WT_04" = "darkgrey")

A <- ggplot(data = All, aes(x = Diff, color = SampleID)) +
  geom_density(alpha = 0.1, linewidth = .5) +  
  theme_classic() +
  scale_color_manual(name = "Sample ID", 
                    values = fill_colors,
                    labels = c("WK_02", 
                               "WK_03",
                               "WT_02",
                               "WT_04"))+
  theme(legend.position = "none", text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  labs(x = "Difference in Percent Methylation per Base", y = "Density")
A

B <- ggplot(All, aes(x=coverage.x, y=Diff, color = SampleID)) +
  geom_point(alpha = .3, size = .5, position = position_jitter(width = 0.3, height = 0))+
  theme_classic() +
  scale_color_manual(name = "Sample ID", 
                     values = fill_colors,
                     labels = c("WK_02", 
                                "WK_03",
                                "WT_02",
                                "WT_04"))+
  theme(legend.position = "none", text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  labs(x = "WGBS Depth", y = "Difference in Percent Methylation per Base")
B

C <- ggplot(All, aes(x=coverage.y, y=Diff, color = SampleID)) +
  geom_point(alpha = .3, size = .5, position = position_jitter(width = 0.3, height = 0))+
  theme_classic() +
  scale_color_manual(name = "Sample ID", 
                     values = fill_colors,
                     labels = c("WK_02", 
                                "WK_03",
                                "WT_02",
                                "WT_04"))+
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), legend.text = element_text(size=12), legend.title = element_text(face = "bold", size = 12))+
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4))) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  labs(x = "RRBS Depth", y = "Difference in Percent Methylation per Base")
C

A + B + C

# Percent Methylation histograms ------------------------------------------


## Overlapping sites -------------------------------------------------------

### Biological replicates ---------------------------------------------------


WT_WGBS <- ggplot(WTWGBS_merged)+
  geom_histogram(aes(x = MeanMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WT_WGBS


WT_RRBS <- ggplot(WTRRBS_merged)+
  geom_histogram(aes(x = MeanMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black"), plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = .5))+
  ylab("")+
  labs(x = "") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WT_RRBS


### Technical Replicates ---------------------------------------------


WT2_WGBS <- ggplot(WT_002_merged)+
  geom_histogram(aes(x = WGBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = .5, size = 5, fontface = "bold")+
  coord_cartesian(clip = "off")
WT2_WGBS


WT2_RRBS <- ggplot(WT_002_merged)+
  geom_histogram(aes(x = RRBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(axis.title.y = element_text(size=14, face="bold", colour = "black"), plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = .5), plot.title.position = "plot")+
  ylab("Frequency")+
  labs(x = "") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WT2_RRBS

WT2_RRBS +  WT_RRBS + WT2_WGBS  + WT_WGBS +  plot_layout(ncol = 2)


## Unique CpG Sites --------------------------------------------------------

{
colnames(RRBS_CG_WT_002)[8] <- "PercMeth"
colnames(WGBS_CG_WT_002)[8] <- "PercMeth"
WT2RRBS_Diff <- dplyr::setdiff(RRBS_CG_WT_002,WGBS_CG_WT_002)
# 1,245,354 CpG sites  unique to RRBS

WT2WGBS_Diff <- dplyr::setdiff(WGBS_CG_WT_002,RRBS_CG_WT_002)
# 18,345,686 CpG sites Unique to WGBS

colnames(RRBS_CG_WT_004)[8] <- "PercMeth"
colnames(WGBS_CG_WT_004)[8] <- "PercMeth"
WT4RRBS_Diff <- dplyr::setdiff(RRBS_CG_WT_004,WGBS_CG_WT_004)
WT4WGBS_Diff <- dplyr::setdiff(WGBS_CG_WT_004,RRBS_CG_WT_004)
}

RRBS <- ggplot(WT2RRBS_Diff)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.y = element_text(size=14, face="bold", colour = "black"), plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = .5))+
  ylab("")+
  labs(x = "") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
RRBS

WGBS <- ggplot(WT2WGBS_Diff)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  labs(x = "Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WGBS

WT2_RRBS +  RRBS + WT_RRBS +   WT2_WGBS  + WGBS +WT_WGBS +  plot_layout(ncol = 3)

RRBS4 <- ggplot(WT4RRBS_Diff)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.y = element_text(size=14, face="bold", colour = "black"), plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = .5))+
  ylab("")+
  labs(x = "") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
RRBS4


WGBS4 <- ggplot(WT4WGBS_Diff)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  labs(x = "Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WGBS4

WT4_WGBS <- ggplot(WT_004_merged)+
  geom_histogram(aes(x = WGBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = .5, size = 5, fontface = "bold")+
  coord_cartesian(clip = "off")
WT4_WGBS

WT4_RRBS <- ggplot(WT_004_merged)+
  geom_histogram(aes(x = RRBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = .5, size = 5, fontface = "bold")+
  coord_cartesian(clip = "off")
WT4_RRBS

WT4_RRBS + RRBS4 + WT4_WGBS + WGBS4 + plot_layout(ncol = 2)

# CpG Annotations -------------------------------------------------------------


## Format data -------------------------------------------------------------

{
WT_002_merged$start <- WT_002_merged$chr
WT_002_merged$end <- WT_002_merged$chr
WT_002_merged$start <- sub("^chr[^.]+\\.", "", WT_002_merged$start)
WT_002_merged$end <- sub("^chr[^.]+\\.", "", WT_002_merged$end)
WT_002_merged$chr <- sub("\\..*", "", WT_002_merged$chr)
head(WT_002_merged)
WT_002_merged$start <- as.numeric(WT_002_merged$start)
WT_002_merged$end <- as.numeric(WT_002_merged$end)

WT2RRBS_Diff$start <- WT2RRBS_Diff$chr
WT2RRBS_Diff$end <- WT2RRBS_Diff$chr
WT2RRBS_Diff$start <- sub("^chr[^.]+\\.", "", WT2RRBS_Diff$start)
WT2RRBS_Diff$end <- sub("^chr[^.]+\\.", "", WT2RRBS_Diff$end)
WT2RRBS_Diff$chr <- sub("\\..*", "", WT2RRBS_Diff$chr)
head(WT2RRBS_Diff)
WT2RRBS_Diff$start <- as.numeric(WT2RRBS_Diff$start)
WT2RRBS_Diff$end <- as.numeric(WT2RRBS_Diff$end)

WT2WGBS_Diff$start <- WT2WGBS_Diff$chr
WT2WGBS_Diff$end <- WT2WGBS_Diff$chr
WT2WGBS_Diff$start <- sub("^chr[^.]+\\.", "", WT2WGBS_Diff$start)
WT2WGBS_Diff$end <- sub("^chr[^.]+\\.", "", WT2WGBS_Diff$end)
WT2WGBS_Diff$chr <- sub("\\..*", "", WT2WGBS_Diff$chr)
head(WT2WGBS_Diff)
WT2WGBS_Diff$start <- as.numeric(WT2WGBS_Diff$start)
WT2WGBS_Diff$end <- as.numeric(WT2WGBS_Diff$end)

WT_004_merged$start <- WT_004_merged$chr
WT_004_merged$end <- WT_004_merged$chr
WT_004_merged$start <- sub("^chr[^.]+\\.", "", WT_004_merged$start)
WT_004_merged$end <- sub("^chr[^.]+\\.", "", WT_004_merged$end)
WT_004_merged$chr <- sub("\\..*", "", WT_004_merged$chr)
head(WT_004_merged)
WT_004_merged$start <- as.numeric(WT_004_merged$start)
WT_004_merged$end <- as.numeric(WT_004_merged$end)
sum(is.na(WT_004_merged))

WT4RRBS_Diff$start <- WT4RRBS_Diff$chr
WT4RRBS_Diff$end <- WT4RRBS_Diff$chr
WT4RRBS_Diff$start <- sub("^chr[^.]+\\.", "", WT4RRBS_Diff$start)
WT4RRBS_Diff$end <- sub("^chr[^.]+\\.", "", WT4RRBS_Diff$end)
WT4RRBS_Diff$chr <- sub("\\..*", "", WT4RRBS_Diff$chr)
head(WT4RRBS_Diff)
WT4RRBS_Diff$start <- as.numeric(WT4RRBS_Diff$start)
WT4RRBS_Diff$end <- as.numeric(WT4RRBS_Diff$end)
sum(is.na(WT4RRBS_Diff))

WT4WGBS_Diff$start <- WT4WGBS_Diff$chr
WT4WGBS_Diff$end <- WT4WGBS_Diff$chr
WT4WGBS_Diff$start <- sub("^chr[^.]+\\.", "", WT4WGBS_Diff$start)
WT4WGBS_Diff$end <- sub("^chr[^.]+\\.", "", WT4WGBS_Diff$end)
WT4WGBS_Diff$chr <- sub("\\..*", "", WT4WGBS_Diff$chr)
head(WT4WGBS_Diff)
WT4WGBS_Diff$start <- as.numeric(WT4WGBS_Diff$start)
WT4WGBS_Diff$end <- as.numeric(WT4WGBS_Diff$end)
sum(is.na(WT4WGBS_Diff))
}

## Convert Methylation DFs to GRanges objects ----------------------------------------------
{
# convert to a GRanges object
WT_002_merged_GR <- GRanges(
  seqnames = WT_002_merged$chr,
  ranges = IRanges(start = WT_002_merged$start, end = WT_002_merged$end),
  strand = "*"
)
head(WT_002_merged_GR)
seqlevels(WT_002_merged_GR)
class(WT_002_merged_GR)

WT2RRBS_Diff_GR <- GRanges(
  seqnames = WT2RRBS_Diff$chr,
  ranges = IRanges(start = WT2RRBS_Diff$start, end = WT2RRBS_Diff$end),
  strand = "*"
)
head(WT2RRBS_Diff_GR)
seqlevels(WT2RRBS_Diff_GR)
class(WT2RRBS_Diff_GR)

WT2WGBS_Diff_GR <- GRanges(
  seqnames = WT2WGBS_Diff$chr,
  ranges = IRanges(start = WT2WGBS_Diff$start, end = WT2WGBS_Diff$end),
  strand = "*"
)
head(WT2WGBS_Diff_GR)
seqlevels(WT2WGBS_Diff_GR)
class(WT2WGBS_Diff_GR)

# convert to a GRanges object
WT_004_merged_GR <- GRanges(
  seqnames = WT_004_merged$chr,
  ranges = IRanges(start = WT_004_merged$start.x, end = WT_004_merged$end.x),
  strand = "*"
)
head(WT_004_merged_GR)
seqlevels(WT_004_merged_GR)
class(WT_004_merged_GR)

WT4RRBS_Diff_GR <- GRanges(
  seqnames = WT4RRBS_Diff$chr,
  ranges = IRanges(start = WT4RRBS_Diff$start, end = WT4RRBS_Diff$end),
  strand = "*"
)
head(WT4RRBS_Diff_GR)
seqlevels(WT4RRBS_Diff_GR)
class(WT4RRBS_Diff_GR)

WT4WGBS_Diff_GR <- GRanges(
  seqnames = WT4WGBS_Diff$chr,
  ranges = IRanges(start = WT4WGBS_Diff$start, end = WT4WGBS_Diff$end),
  strand = "*"
)
head(WT4WGBS_Diff_GR)
seqlevels(WT4WGBS_Diff_GR)
class(WT4WGBS_Diff_GR)
}

## CpG Islands -------------------------------------------------------------


# Import CpG Islands based on reference genome. stickle_CGI.bed from running TaJoCGI on the stickleback reference genome to find CpG islands
CpGI_anot <- cpg_anot <- readFeatureFlank("/stickle_CGI.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
head(CpGI_anot)

# RRBS and WGBS merged sites
WT_002_merged_CpGI <- annotateWithFeatureFlank(as(WT_002_merged_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WT_002_merged_CpGI
# Rows in target set: 781846
# percentage of target elements overlapping with features:
#   CpGi shores  other 
# 43.59  21.14  35.27 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 50.45  33.46 

Merged <- plotTargetAnnotation(WT_002_merged_CpGI, main = "RRBS & WGBS Merged")


# RRBS only sites
WT2RRBS_Diff_CpGI <- annotateWithFeatureFlank(as(WT2RRBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WT2RRBS_Diff_CpGI
# Rows in target set: 1245354
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 48.52  20.07  31.41 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 61.60  38.74 

RRBS <- plotTargetAnnotation(WT2RRBS_Diff_CpGI, main = "Unique to RRBS")

# WGBS only sites
WT2WGBS_Diff_CpGI <- annotateWithFeatureFlank(as(WT2WGBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WT2WGBS_Diff_CpGI
# Rows in target set: 18345686
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 14.05  27.07  58.87 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 95.01  95.93 

WGBS <- plotTargetAnnotation(WT2WGBS_Diff_CpGI, main = "Unique to WGBS")


# RRBS and WGBS merged sites
WT_004_merged_CpGI <- annotateWithFeatureFlank(as(WT_004_merged_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WT_004_merged_CpGI
# Rows in target set: 803066
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 45.16  20.43  34.41 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 41.94  29.27 

plotTargetAnnotation(WT_002_merged_CpGI, main = "RRBS & WGBS Merged")

# RRBS only sites
WT4RRBS_Diff_CpGI <- annotateWithFeatureFlank(as(WT4RRBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WT4RRBS_Diff_CpGI
# Rows in target set: 1233607
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 48.64  19.81  31.55 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 51.59  34.06 

plotTargetAnnotation(WT4RRBS_Diff_CpGI, main = "Unique to RRBS")

# WGBS only sites
WT4WGBS_Diff_CpGI <- annotateWithFeatureFlank(as(WT4WGBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WT4WGBS_Diff_CpGI
# Rows in target set: 17826537
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 14.04  26.97  58.99 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 94.95  95.95  

plotTargetAnnotation(WT4WGBS_Diff_CpGI, main = "Unique to WGBS")

## Genome annotations ----------------------------------


### Change seq levels -------------------------------------------------------

{
### Load the genome annotation data; i.e the coordinates of promoters, TSS, intron and exons
refseq_anot <- readTranscriptFeatures("/GCF_016920845.1_GAculeatus_UGA_version5_genomic.bed",remove.unusual=FALSE)
head(refseq_anot)
seqlevels(refseq_anot)

# remove unmapped chromosomes
refseq_anot_filtered <- GRangesList(
  lapply(refseq_anot, function(gr) {
    # Get seqnames as character vector for this GRanges
    seqs <- as.character(seqnames(gr))
    # Filter out "NW"
    gr[!grepl("^NW", seqs)]
  })
)
seqlevels(refseq_anot_filtered)

# Drop empty elements
refseq_anot_filtered <- refseq_anot_filtered[elementNROWS(refseq_anot_filtered) > 0]

# Drop unused seqlevels in each GRanges
refseq_anot_filtered <- endoapply(refseq_anot_filtered, function(gr) {
  keepSeqlevels(gr, unique(seqnames(gr)), pruning.mode = "coarse")
})
seqlevels(refseq_anot_filtered)

seqlevels(WT_002_merged_GR)
seqlevels(WT2RRBS_Diff_GR)
seqlevels(WT2WGBS_Diff_GR)

# Make sure chromosome names match in both reference and data file 

chr_mapping <- c(
  "chrI" = "NC_053212.1",
  "chrII" = "NC_053213.1", 
  "chrIII" = "NC_053214.1",
  "chrIV" = "NC_053215.1",
  "chrV" = "NC_053216.1",
  "chrVI" = "NC_053217.1",
  "chrVII" = "NC_053218.1", 
  "chrVIII" = "NC_053219.1",
  "chrIX" = "NC_053220.1", 
  "chrX" = "NC_053221.1",
  "chrXI" = "NC_053222.1",
  "chrXII" = "NC_053223.1",
  "chrXIII" = "NC_053224.1",
  "chrXIV" = "NC_053225.1",
  "chrXV" = "NC_053226.1",
  "chrXVI" = "NC_053227.1",
  "chrXVII" = "NC_053228.1",
  "chrXVIII" = "NC_053229.1",
  "chrXIX" = "NC_053230.1",
  "chrXX" = "NC_053231.1",
  "chrXXI" = "NC_053232.1",
  "chrY" = "NC_053233.1"
)



WT_002_merged_GR <- renameSeqlevels(WT_002_merged_GR, chr_mapping)
WT2RRBS_Diff_GR <- renameSeqlevels(WT2RRBS_Diff_GR, chr_mapping)
WT2WGBS_Diff_GR <- renameSeqlevels(WT2WGBS_Diff_GR, chr_mapping)

WT_004_merged_GR <- renameSeqlevels(WT_004_merged_GR, chr_mapping)
WT4RRBS_Diff_GR <- renameSeqlevels(WT4RRBS_Diff_GR, chr_mapping)
WT4WGBS_Diff_GR <- renameSeqlevels(WT4WGBS_Diff_GR, chr_mapping)

seqlevels(WT_002_merged_GR)
seqlevels(WT2RRBS_Diff_GR)
seqlevels(WT2WGBS_Diff_GR)

seqlevels(WT_004_merged_GR)
seqlevels(WT4RRBS_Diff_GR)
seqlevels(WT4WGBS_Diff_GR)

seqlevels(WT_002_merged_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
seqlevels(WT_002_merged_GR)

seqlevels(WT2RRBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
seqlevels(WT2RRBS_Diff_GR)

seqlevels(WT2WGBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
seqlevels(WT2WGBS_Diff_GR)

seqlevels(WT_004_merged_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
seqlevels(WT_004_merged_GR)

seqlevels(WT4RRBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
seqlevels(WT4RRBS_Diff_GR)

seqlevels(WT4WGBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
seqlevels(WT4WGBS_Diff_GR)
}


### Annotate WGBS and RRBS  -------------------------------------------------


WT_002_merged.anot <- annotateWithGeneParts(target = as(WT_002_merged_GR,"GRanges"),
                                       feature = refseq_anot_filtered, strand = TRUE)
WT_002_merged.anot
# Rows in target set: 765054
#
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 35.69      42.84      39.03      17.59 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 35.69      19.99      26.73      17.59 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 55.36    12.32    13.13 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     488    2710    9313   10214  278984 

WT2RRBS_Diff.anot <- annotateWithGeneParts(target = as(WT2RRBS_Diff_GR,"GRanges"),
                                         feature = refseq_anot_filtered, strand = TRUE)
WT2RRBS_Diff.anot
# Summary of target set annotation with genic parts
# Rows in target set: 1217149

#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 30.31      43.67      37.89      18.10 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 30.31      24.06      27.53      18.10 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 58.67    14.97    15.16 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     663    3524    9975   11300  280838 

WT2WGBS_Diff.anot <- annotateWithGeneParts(target = as(WT2WGBS_Diff_GR,"GRanges"),
                                        feature = refseq_anot_filtered, strand = TRUE)
WT2WGBS_Diff.anot
# Rows in target set: 18004229
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 14.93      21.97      50.30      26.16 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 14.93      15.83      43.08      26.16 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 96.83    88.84    90.28 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    2090    6239   14109   16712  282513 

WT_004_merged.anot <- annotateWithGeneParts(target = as(WT_004_merged_GR,"GRanges"),
                                            feature = refseq_anot_filtered, strand = TRUE)
WT_004_merged.anot
# Rows in target set: 784737
#
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 45.37      43.42      39.35      16.24 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 45.37      15.12      23.27      16.24 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 56.49     9.81    11.74 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     329    1397    7874    7759  279037 

WT4RRBS_Diff.anot <- annotateWithGeneParts(target = as(WT4RRBS_Diff_GR,"GRanges"),
                                           feature = refseq_anot_filtered, strand = TRUE)
WT4RRBS_Diff.anot
# Rows in target set: 1204809
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 39.80      43.26      38.49      17.32 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 39.80      18.13      24.75      17.32 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 59.38    11.95    13.42 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     400    2148    8683    9317  281696 

WT4WGBS_Diff.anot <- annotateWithGeneParts(target = as(WT4WGBS_Diff_GR,"GRanges"),
                                           feature = refseq_anot_filtered, strand = TRUE)
WT4WGBS_Diff.anot
# Rows in target set: 17484051
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 15.08      21.99      50.24      26.17 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 15.08      15.78      42.96      26.17 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 96.83    88.40    89.99 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    2073    6213   14117   16682  282513

# Extract distance to closest TSS merged sites
WT2dist_tss <- getAssociationWithTSS(WT_002_merged.anot)
head(WT2dist_tss)

hist(WT2dist_tss$dist.to.feature)
plotTargetAnnotation(WT_002_merged.anot, main = "RRBS & WGBS Merged Sites")

WT4dist_tss <- getAssociationWithTSS(WT_004_merged.anot)
head(WT4dist_tss)

hist(WT4dist_tss$dist.to.feature)
plotTargetAnnotation(WT_004_merged.anot, main = "RRBS & WGBS Merged Sites")

# Extract distance to nearest TSS RRBS only sites
WT2dist_tss2 <- getAssociationWithTSS(WT2RRBS_Diff.anot)
head(WT2dist_tss2)

hist(WT2dist_tss2$dist.to.feature)
plotTargetAnnotation(WT2RRBS_Diff.anot, main = "Unique to RRBS")

WT4dist_tss2 <- getAssociationWithTSS(WT4RRBS_Diff.anot)
head(WT4dist_tss2)

hist(WT4dist_tss2$dist.to.feature)
plotTargetAnnotation(WT4RRBS_Diff.anot, main = "Unique to RRBS")

# Extract distance to nearest TSS WGBS only sites
WT2dist_tss3 <- getAssociationWithTSS(WT2WGBS_Diff.anot)
head(WT2dist_tss3)

hist(WT2dist_tss3$dist.to.feature)
plotTargetAnnotation(WT2WGBS_Diff.anot, main = "Unique to WGBS")

WT4dist_tss3 <- getAssociationWithTSS(WT4WGBS_Diff.anot)
head(WT4dist_tss3)

hist(WT4dist_tss3$dist.to.feature)
plotTargetAnnotation(WT4WGBS_Diff.anot, main = "Unique to WGBS")

# merge sites covered by both, RRBS, and WGBS, respectively
WT_Both_dist_tss <- rbind(WT2dist_tss, WT4dist_tss)
WT_Both_dist_tss$Method <- "Both"
WT_Both_dist_tss$Population <- "Watson"
head(WT_Both_dist_tss)

WT_RRBS_dist_tss <- rbind(WT2dist_tss2, WT4dist_tss2)
WT_RRBS_dist_tss$Method <- "RRBS"
WT_RRBS_dist_tss$Population <- "Watson"

WT_WGBS_dist_tss <- rbind(WT2dist_tss3, WT4dist_tss3)
WT_WGBS_dist_tss$Method <- "WGBS"
WT_WGBS_dist_tss$Population <- "Watson"
head(WT_WGBS_dist_tss)

WT_all_dist_tss <- rbind(WT_Both_dist_tss, WT_RRBS_dist_tss)
WT_all_dist_tss <- rbind(WT_all_dist_tss, WT_WGBS_dist_tss)

write.csv(WT_all_dist_tss, "/Figures//WT_dist_tss.csv")

