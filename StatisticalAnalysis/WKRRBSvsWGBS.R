### Assess overlapping CpG sites between sequencing methods

# load packages
{
library(methylKit)
library(ggvenn)
library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(GenomicPlot)
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(ggridges)
library(readxl)
}

# Import and format data -------------------------------------------------------



#### Min Depth 5x

setwd("/MethylDackel/Stickle/NoSNPFilters")
{
# import samples
WGBSfile.list=list("WGBS_CG_22_WK_002_CpG.methylKit", "WGBS_CG_22_WK_003_CpG.methylKit", "CG_22_WK_002_CpG.methylKit", "CG_22_WK_003_CpG.methylKit")

WGBSvsRRBS=methRead(WGBSfile.list,
                    sample.id=list("WGBS_CG_WK_002", "WGBS_CG_WK_003", "RRBS_CG_WK_002", "RRBS_CG_WK_003"),
                    pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                    assembly="StickleGeneAnnotations",
                    treatment=c(1,1,2,2),
                    context="CpG",
                    mincov = 5
)

WGBSvsRRBS=filterByCoverage(WGBSvsRRBS,lo.count=5,lo.perc=NULL,
                            hi.count=NULL,hi.perc=99.9)
}

{
# Separate each sample into individual methylRaw objects - contains info about mehtylation and chr location
WGBS_CG_WK_002=WGBSvsRRBS[[1]]
WGBS_CG_WK_003=WGBSvsRRBS[[2]]
RRBS_CG_WK_002=WGBSvsRRBS[[3]]
RRBS_CG_WK_003=WGBSvsRRBS[[4]]

WGBS_CG_WK_002$WGBS_PercMeth <- ((WGBS_CG_WK_002$numCs/WGBS_CG_WK_002$coverage)*100)
WGBS_CG_WK_003$WGBS_PercMeth <- ((WGBS_CG_WK_003$numCs/WGBS_CG_WK_003$coverage)*100)
RRBS_CG_WK_002$RRBS_PercMeth <- ((RRBS_CG_WK_002$numCs/RRBS_CG_WK_002$coverage)*100)
RRBS_CG_WK_003$RRBS_PercMeth <- ((RRBS_CG_WK_003$numCs/RRBS_CG_WK_003$coverage)*100)

WGBS_CG_WK_002 <- methylKit::getData(WGBS_CG_WK_002)
WGBS_CG_WK_003 <- methylKit::getData(WGBS_CG_WK_003)
RRBS_CG_WK_002 <- methylKit::getData(RRBS_CG_WK_002)
RRBS_CG_WK_003 <- methylKit::getData(RRBS_CG_WK_003)


## Overlapping CpG sites on biological replicates with RRBS and WGBS -------


WKWGBS_merged <- merge(WGBS_CG_WK_002, WGBS_CG_WK_003, by = "chr")
nrow(WKWGBS_merged)
# 13,937,683 CpG sites

WKRRBS_merged <- merge(RRBS_CG_WK_002, RRBS_CG_WK_003, by = "chr")
nrow(WKRRBS_merged)
# 386,965 CpG sites

WKWGBS_merged$MeanMeth <- rowMeans(WKWGBS_merged[, c("WGBS_PercMeth.x", "WGBS_PercMeth.y")], na.rm = TRUE)

WKRRBS_merged$MeanMeth <- rowMeans(WKRRBS_merged[, c("RRBS_PercMeth.x", "RRBS_PercMeth.y")], na.rm = TRUE)


## Overlapping CpG sites with RRBS and WGBS on technical replicates --------

WK_002_merged <- merge(WGBS_CG_WK_002, RRBS_CG_WK_002, by = "chr")
hist(WK_002_merged$WGBS_PercMeth)
hist(WK_002_merged$RRBS_PercMeth)
nrow(WK_002_merged) # 542,910


WK_003_merged <- merge(WGBS_CG_WK_003, RRBS_CG_WK_003, by = "chr")
hist(WK_003_merged$WGBS_PercMeth)
hist(WK_003_merged$RRBS_PercMeth)
nrow(WK_003_merged) # 399,648

}
  

WK_002_merged$Diff <- WK_002_merged$WGBS_PercMeth - WK_002_merged$RRBS_PercMeth
WK_003_merged$Diff <- WK_003_merged$WGBS_PercMeth - WK_003_merged$RRBS_PercMeth

WK_002_merged$SampleID <- "WK_02"
WK_003_merged$SampleID <- "WK_03"

AllWK <- rbind(WK_002_merged, WK_003_merged)

# Need to create AllWT from WTRRBSvsWGBS.R script before running this 
All <- rbind(AllWK, AllWT)

ggplot(AllWK, aes(x=coverage.x, y=Diff)) +
  geom_point()+
  theme_classic()+
  labs(x = "WGBS Depth", y = "Difference in Percent Methylation per Base")

ggplot(AllWK, aes(x=coverage.y, y=Diff)) +
  geom_point() +
  theme_classic()+
  labs(x = "RRBS Depth", y = "Difference in Percent Methylation per Base")

ggplot(data = AllWK, aes(x = Diff, color = SampleID)) +
  geom_density(alpha = 0.4) +  
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Difference in Percent Methylation per Base", y = "Frequency")



WK_002_sites <- WK_002_merged %>%
  mutate(WGBS_bin = cut(WGBS_PercMeth, breaks = 10),
         RRBS_bin = cut(RRBS_PercMeth, breaks = 10)) %>%
  group_by(WGBS_bin, RRBS_bin) %>%
  summarise(mean_diff = WGBS_PercMeth - RRBS_PercMeth, n_sites = n(), .groups = "drop")


ggplot(WK_002_sites, aes(x = WGBS_bin, y = RRBS_bin, fill = mean_diff)) +
  geom_point(aes(size = n_sites), 
             shape = 21, color = "black") + 
  scale_fill_gradient2(low = "deepskyblue3", mid = "lightgrey", high = "yellow2", midpoint = 0, name = "Methylation\nDifference") +
  scale_size_continuous(range = c(2, 15), name = "Number of\nSites") +
  labs(title = "WK_002", x = "WGBS (% Methylation Per Base)", y = "RRBS (% Methylation per Base)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggplot(WK_002_merged, aes(x=WGBS_PercMeth, y=RRBS_PercMeth)) +
  geom_point()

WK_002_sites <- WK_002_merged %>%
  mutate(WGBS_bin = cut(WGBS_PercMeth, breaks = 10),
         RRBS_bin = cut(RRBS_PercMeth, breaks = 10)) %>%
  group_by(WGBS_bin, RRBS_bin) %>%
  summarise(mean_diff = WGBS_PercMeth - RRBS_PercMeth, n_sites = n(), .groups = "drop")


ggplot(WK_002_sites, aes(x = WGBS_bin, y = RRBS_bin, fill = mean_diff)) +
  geom_point(aes(size = n_sites), 
             shape = 21, color = "black") + 
  scale_fill_gradient2(low = "deepskyblue3", mid = "lightgrey", high = "yellow2", midpoint = 0, name = "Methylation\nDifference") +
  scale_size_continuous(range = c(2, 15), name = "Number of\nSites") +
  labs(title = "WK_002", x = "WGBS (% Methylation Per Base)", y = "RRBS (% Methylation per Base)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  
# Percent Methylation histograms ------------------------------------------
  
  
## Overlapping sites -------------------------------------------------------
  
### Biological replicates ---------------------------------------------------

WGBS_binned <- WKWGBS_merged %>%
  mutate(bin = floor(MeanMeth / 5) * 5) %>%
  count(bin) %>%
  mutate(fill_group = case_when(
    bin < 15  ~ "<15%",
    bin >= 85 ~ ">85%",
    TRUE      ~ "15-85%"
  ))
  
WK_WGBS <- ggplot(WGBS_binned) +
  geom_col(aes(x = bin + 2.5, y = n, fill = fill_group),
           alpha = .8, color = "black", width = 5) +
  scale_fill_manual(values = c("<15%" = "pink", "15-85%" = "lightgreen", ">85%" = "lightblue")) +
  theme_classic() +
  theme(text = element_text(size = 14),
        axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
        legend.position = "none") +
  ylab("Frequency") +
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK_WGBS


WK_WGBS <- ggplot(WKWGBS_merged)+
    geom_histogram(aes(x = MeanMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
    #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
    ylab("Frequency")+
    xlab("") +
    annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK_WGBS
  
  
WK_RRBS <- ggplot(WKRRBS_merged)+
    geom_histogram(aes(x = MeanMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
    #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
    ylab("Frequency")+
    labs(x = "Percent Methylation per Base") +
    annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK_RRBS

WK_WGBS + WK_RRBS + plot_layout(ncol = 1)
  
### Technical Replicates ---------------------------------------------
  
  
WK2_WGBS <- ggplot(WGBS_CG_WK_002)+
    geom_histogram(aes(x = WGBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK2_WGBS
  
  
WK2_RRBS <- ggplot(RRBS_CG_WK_002)+
    geom_histogram(aes(x = RRBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  labs(x = "Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK2_RRBS

WK3_WGBS <- ggplot(WGBS_CG_WK_003)+
  geom_histogram(aes(x = WGBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK3_WGBS


WK3_RRBS <- ggplot(RRBS_CG_WK_003)+
  geom_histogram(aes(x = RRBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("Frequency")+
  labs(x = "Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WK3_RRBS
  
WK2_WGBS + WK2_RRBS + WK3_WGBS + WK3_RRBS + plot_layout(ncol = 2)
  
  
## Unique CpG Sites --------------------------------------------------------
  
{
  colnames(RRBS_CG_WK_002)[8] <- "PercMeth"
  colnames(WGBS_CG_WK_002)[8] <- "PercMeth"
  WK2RRBS_Diff <- dplyr::setdiff(RRBS_CG_WK_002,WGBS_CG_WK_002)
    # 1,245,354 CpG sites  unique to RRBS
    
  WK2WGBS_Diff <- dplyr::setdiff(WGBS_CG_WK_002,RRBS_CG_WK_002)
    # 18,345,686 CpG sites Unique to WGBS
    
  colnames(RRBS_CG_WK_003)[8] <- "PercMeth"
  colnames(WGBS_CG_WK_003)[8] <- "PercMeth"
  WK3RRBS_Diff <- dplyr::setdiff(RRBS_CG_WK_003,WGBS_CG_WK_003)
  WK3WGBS_Diff <- dplyr::setdiff(WGBS_CG_WK_003,RRBS_CG_WK_003)
  }
  
RRBS <- ggplot(WK2RRBS_Diff)+
    geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.y = element_text(size=14, face="bold", colour = "black"), plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = .5))+
    ylab("")+
    labs(x = "") +
    annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
RRBS
  
WGBS <- ggplot(WK2WGBS_Diff)+
    geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
    ylab("")+
    labs(x = "Percent Methylation per Base") +
    annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WGBS
  
WK2_RRBS +  RRBS + WK_RRBS +   WK2_WGBS  + WGBS +WK_WGBS +  plot_layout(ncol = 3)
  
RRBS4 <- ggplot(WK3RRBS_Diff)+
    geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.y = element_text(size=14, face="bold", colour = "black"), plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = .5))+
    ylab("")+
    labs(x = "") +
    annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
RRBS4
  
  
WGBS4 <- ggplot(WK3WGBS_Diff)+
    geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
    ylab("")+
    labs(x = "Percent Methylation per Base") +
    annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
WGBS4
  
WK3_WGBS <- ggplot(WK_003_merged)+
    geom_histogram(aes(x = WGBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
    ylab("Frequency")+
    xlab("Percent Methylation per Base") +
    annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = .5, size = 5, fontface = "bold")+
    coord_cartesian(clip = "off")
WK3_WGBS
  
WK3_RRBS <- ggplot(WK_003_merged)+
    geom_histogram(aes(x = RRBS_PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
    theme_classic() +
    theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
    ylab("Frequency")+
    xlab("") +
    annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = .5, size = 5, fontface = "bold")+
    coord_cartesian(clip = "off")
WK3_RRBS
  
WK3_RRBS + RRBS4 + WK3_WGBS + WGBS4 + plot_layout(ncol = 2)
  
# CpG Annotations -------------------------------------------------------------
  
  
## Format data -------------------------------------------------------------
  
  {
    WK_002_merged$start <- WK_002_merged$chr
    WK_002_merged$end <- WK_002_merged$chr
    WK_002_merged$start <- sub("^chr[^.]+\\.", "", WK_002_merged$start)
    WK_002_merged$end <- sub("^chr[^.]+\\.", "", WK_002_merged$end)
    WK_002_merged$chr <- sub("\\..*", "", WK_002_merged$chr)
    head(WK_002_merged)
    WK_002_merged$start <- as.numeric(WK_002_merged$start)
    WK_002_merged$end <- as.numeric(WK_002_merged$end)
    
    WK2RRBS_Diff$start <- WK2RRBS_Diff$chr
    WK2RRBS_Diff$end <- WK2RRBS_Diff$chr
    WK2RRBS_Diff$start <- sub("^chr[^.]+\\.", "", WK2RRBS_Diff$start)
    WK2RRBS_Diff$end <- sub("^chr[^.]+\\.", "", WK2RRBS_Diff$end)
    WK2RRBS_Diff$chr <- sub("\\..*", "", WK2RRBS_Diff$chr)
    head(WK2RRBS_Diff)
    WK2RRBS_Diff$start <- as.numeric(WK2RRBS_Diff$start)
    WK2RRBS_Diff$end <- as.numeric(WK2RRBS_Diff$end)
    
    WK2WGBS_Diff$start <- WK2WGBS_Diff$chr
    WK2WGBS_Diff$end <- WK2WGBS_Diff$chr
    WK2WGBS_Diff$start <- sub("^chr[^.]+\\.", "", WK2WGBS_Diff$start)
    WK2WGBS_Diff$end <- sub("^chr[^.]+\\.", "", WK2WGBS_Diff$end)
    WK2WGBS_Diff$chr <- sub("\\..*", "", WK2WGBS_Diff$chr)
    head(WK2WGBS_Diff)
    WK2WGBS_Diff$start <- as.numeric(WK2WGBS_Diff$start)
    WK2WGBS_Diff$end <- as.numeric(WK2WGBS_Diff$end)
    
    WK_003_merged$start <- WK_003_merged$chr
    WK_003_merged$end <- WK_003_merged$chr
    WK_003_merged$start <- sub("^chr[^.]+\\.", "", WK_003_merged$start)
    WK_003_merged$end <- sub("^chr[^.]+\\.", "", WK_003_merged$end)
    WK_003_merged$chr <- sub("\\..*", "", WK_003_merged$chr)
    head(WK_003_merged)
    WK_003_merged$start <- as.numeric(WK_003_merged$start)
    WK_003_merged$end <- as.numeric(WK_003_merged$end)
    sum(is.na(WK_003_merged))
    
    WK3RRBS_Diff$start <- WK3RRBS_Diff$chr
    WK3RRBS_Diff$end <- WK3RRBS_Diff$chr
    WK3RRBS_Diff$start <- sub("^chr[^.]+\\.", "", WK3RRBS_Diff$start)
    WK3RRBS_Diff$end <- sub("^chr[^.]+\\.", "", WK3RRBS_Diff$end)
    WK3RRBS_Diff$chr <- sub("\\..*", "", WK3RRBS_Diff$chr)
    head(WK3RRBS_Diff)
    WK3RRBS_Diff$start <- as.numeric(WK3RRBS_Diff$start)
    WK3RRBS_Diff$end <- as.numeric(WK3RRBS_Diff$end)
    sum(is.na(WK3RRBS_Diff))
    
    WK3WGBS_Diff$start <- WK3WGBS_Diff$chr
    WK3WGBS_Diff$end <- WK3WGBS_Diff$chr
    WK3WGBS_Diff$start <- sub("^chr[^.]+\\.", "", WK3WGBS_Diff$start)
    WK3WGBS_Diff$end <- sub("^chr[^.]+\\.", "", WK3WGBS_Diff$end)
    WK3WGBS_Diff$chr <- sub("\\..*", "", WK3WGBS_Diff$chr)
    head(WK3WGBS_Diff)
    WK3WGBS_Diff$start <- as.numeric(WK3WGBS_Diff$start)
    WK3WGBS_Diff$end <- as.numeric(WK3WGBS_Diff$end)
    sum(is.na(WK3WGBS_Diff))
  }
  
## Convert Methylation DFs to GRanges objects ----------------------------------------------
  {
    # convert to a GRanges object
    WK_002_merged_GR <- GRanges(
      seqnames = WK_002_merged$chr,
      ranges = IRanges(start = WK_002_merged$start, end = WK_002_merged$end),
      strand = "*"
    )
    head(WK_002_merged_GR)
    seqlevels(WK_002_merged_GR)
    class(WK_002_merged_GR)
    
    WK2RRBS_Diff_GR <- GRanges(
      seqnames = WK2RRBS_Diff$chr,
      ranges = IRanges(start = WK2RRBS_Diff$start, end = WK2RRBS_Diff$end),
      strand = "*"
    )
    head(WK2RRBS_Diff_GR)
    seqlevels(WK2RRBS_Diff_GR)
    class(WK2RRBS_Diff_GR)
    
    WK2WGBS_Diff_GR <- GRanges(
      seqnames = WK2WGBS_Diff$chr,
      ranges = IRanges(start = WK2WGBS_Diff$start, end = WK2WGBS_Diff$end),
      strand = "*"
    )
    head(WK2WGBS_Diff_GR)
    seqlevels(WK2WGBS_Diff_GR)
    class(WK2WGBS_Diff_GR)
    
    # convert to a GRanges object
    WK_003_merged_GR <- GRanges(
      seqnames = WK_003_merged$chr,
      ranges = IRanges(start = WK_003_merged$start.x, end = WK_003_merged$end.x),
      strand = "*"
    )
    head(WK_003_merged_GR)
    seqlevels(WK_003_merged_GR)
    class(WK_003_merged_GR)
    
    WK3RRBS_Diff_GR <- GRanges(
      seqnames = WK3RRBS_Diff$chr,
      ranges = IRanges(start = WK3RRBS_Diff$start, end = WK3RRBS_Diff$end),
      strand = "*"
    )
    head(WK3RRBS_Diff_GR)
    seqlevels(WK3RRBS_Diff_GR)
    class(WK3RRBS_Diff_GR)
    
    WK3WGBS_Diff_GR <- GRanges(
      seqnames = WK3WGBS_Diff$chr,
      ranges = IRanges(start = WK3WGBS_Diff$start, end = WK3WGBS_Diff$end),
      strand = "*"
    )
    head(WK3WGBS_Diff_GR)
    seqlevels(WK3WGBS_Diff_GR)
    class(WK3WGBS_Diff_GR)
  }
  
## CpG Islands -------------------------------------------------------------
  
  
# Import CpG Islands based on reference genome. TaJoCGI used on the reference genome to create CpG island bed file
CpGI_anot <- cpg_anot <- readFeatureFlank("/stickle_CGI.bed", feature.flank.name = c("CpGi", "shores"), flank=2000)
head(CpGI_anot)
  
# RRBS and WGBS merged sites
WK_002_merged_CpGI <- annotateWithFeatureFlank(as(WK_002_merged_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WK_002_merged_CpGI
# Rows in target set: 1000272
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 44.07  21.09  34.84 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 47.41  33.76
  
plotTargetAnnotation(WK_002_merged_CpGI, main = "RRBS & WGBS Merged")
  
  
# RRBS only sites
WK2RRBS_Diff_CpGI <- annotateWithFeatureFlank(as(WK2RRBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WK2RRBS_Diff_CpGI
# Rows in target set: 1438344
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 47.25  20.43  32.32 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 56.11  37.92
  
plotTargetAnnotation(WK2RRBS_Diff_CpGI, main = "Unique to RRBS")
  
# WGBS only sites
WK2WGBS_Diff_CpGI <- annotateWithFeatureFlank(as(WK2WGBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WK2WGBS_Diff_CpGI
# Rows in target set: 19152179
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 14.41  26.89  58.69 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 95.76  97.30
  
plotTargetAnnotation(WK2WGBS_Diff_CpGI, main = "Unique to WGBS")
  
  
# RRBS and WGBS merged sites
WK_003_merged_CpGI <- annotateWithFeatureFlank(as(WK_003_merged_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WK_003_merged_CpGI
# Rows in target set: 752438
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 47.08  20.15  32.77 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 41.49  27.93 
  
plotTargetAnnotation(WK_002_merged_CpGI, main = "RRBS & WGBS Merged")
  
# RRBS only sites
WK3RRBS_Diff_CpGI <- annotateWithFeatureFlank(as(WK3RRBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WK3RRBS_Diff_CpGI
# Rows in target set: 1124674
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 50.96  19.34  29.71 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 50.46  32.05 
  
plotTargetAnnotation(WK3RRBS_Diff_CpGI, main = "Unique to RRBS")
  
# WGBS only sites
WK3WGBS_Diff_CpGI <- annotateWithFeatureFlank(as(WK3WGBS_Diff_GR,"GRanges"), feature = CpGI_anot$CpGi, flank = CpGI_anot$shores, feature.name = "CpGi", flank.name = "shores")
WK3WGBS_Diff_CpGI
# Rows in target set: 18704536
# 
#   percentage of target elements overlapping with features:
#   CpGi shores  other 
# 14.25  27.09  58.65 
# 
# percentage of feature elements overlapping with target:
#   CpGi shores 
# 95.06  95.97  
  
plotTargetAnnotation(WK3WGBS_Diff_CpGI, main = "Unique to WGBS")
  
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
    
    seqlevels(WK_002_merged_GR)
    seqlevels(WK2RRBS_Diff_GR)
    seqlevels(WK2WGBS_Diff_GR)
    
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
    
    
    
    WK_002_merged_GR <- renameSeqlevels(WK_002_merged_GR, chr_mapping)
    WK2RRBS_Diff_GR <- renameSeqlevels(WK2RRBS_Diff_GR, chr_mapping)
    WK2WGBS_Diff_GR <- renameSeqlevels(WK2WGBS_Diff_GR, chr_mapping)
    
    WK_003_merged_GR <- renameSeqlevels(WK_003_merged_GR, chr_mapping)
    WK3RRBS_Diff_GR <- renameSeqlevels(WK3RRBS_Diff_GR, chr_mapping)
    WK3WGBS_Diff_GR <- renameSeqlevels(WK3WGBS_Diff_GR, chr_mapping)
    
    seqlevels(WK_002_merged_GR)
    seqlevels(WK2RRBS_Diff_GR)
    seqlevels(WK2WGBS_Diff_GR)
    
    seqlevels(WK_003_merged_GR)
    seqlevels(WK3RRBS_Diff_GR)
    seqlevels(WK3WGBS_Diff_GR)
    
    seqlevels(WK_002_merged_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
    seqlevels(WK_002_merged_GR)
    
    seqlevels(WK2RRBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
    seqlevels(WK2RRBS_Diff_GR)
    
    seqlevels(WK2WGBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
    seqlevels(WK2WGBS_Diff_GR)
    
    seqlevels(WK_003_merged_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
    seqlevels(WK_003_merged_GR)
    
    seqlevels(WK3RRBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
    seqlevels(WK3RRBS_Diff_GR)
    
    seqlevels(WK3WGBS_Diff_GR, pruning.mode="coarse") <- seqlevels(refseq_anot_filtered)
    seqlevels(WK3WGBS_Diff_GR)
  }
  
  
### Annotate WGBS and RRBS  -------------------------------------------------
  
  
WK_002_merged.anot <- annotateWithGeneParts(target = as(WK_002_merged_GR,"GRanges"), feature = refseq_anot_filtered, strand = TRUE)
WK_002_merged.anot
# Rows in target set: 977974
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 41.38      42.05      39.65      17.48 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 41.38      16.19      24.95      17.48 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 58.87    11.39    13.42 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     378    1921    8743    9167  281696 
  
WK2RRBS_Diff.anot <- annotateWithGeneParts(target = as(WK2RRBS_Diff_GR,"GRanges"), feature = refseq_anot_filtered, strand = TRUE)
WK2RRBS_Diff.anot
# Rows in target set: 1404911
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 36.17      41.54      39.14      18.56 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 36.17      18.84      26.43      18.56 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 61.65    13.33    14.95 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     470    2747    9492   10541  281696 
  
WK2WGBS_Diff.anot <- annotateWithGeneParts(target = as(WK2WGBS_Diff_GR,"GRanges"), feature = refseq_anot_filtered, strand = TRUE)
WK2WGBS_Diff.anot
# Rows in target set: 18781221
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 14.94      22.15      50.00      26.32 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 14.94      15.94      42.80      26.32 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 98.59    89.98    91.59 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    2093    6254   14136   16707  282513 

  
WK_003_merged.anot <- annotateWithGeneParts(target = as(WK_003_merged_GR,"GRanges"), feature = refseq_anot_filtered, strand = TRUE)
WK_003_merged.anot
# Rows in target set: 735151
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 45.69      44.68      38.57      15.84 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 45.69      15.62      22.85      15.84 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 54.95     9.39    11.03 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     322    1358    7632    7551  282116 
  
WK3RRBS_Diff.anot <- annotateWithGeneParts(target = as(WK3RRBS_Diff_GR,"GRanges"), feature = refseq_anot_filtered, strand = TRUE)
WK3RRBS_Diff.anot
# Rows in target set: 1098558
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 40.97      44.88      37.63      16.61 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 40.97      18.46      23.96      16.61 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 57.45    11.24    12.46 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0     381    1955    8331    8866  282116

  
WK3WGBS_Diff.anot <- annotateWithGeneParts(target = as(WK3WGBS_Diff_GR,"GRanges"), feature = refseq_anot_filtered, strand = TRUE)
WK3WGBS_Diff.anot
# Rows in target set: 18350129
# 
#   percentage of target features overlapping with annotation:
#   promoter       exon     intron intergenic 
# 14.74      22.00      50.30      26.16 
# 
# percentage of target features overlapping with annotation:
#   (with promoter > exon > intron precedence):
#   promoter       exon     intron intergenic 
# 14.74      15.94      43.16      26.16 
# 
# percentage of annotation boundaries with feature overlap:
#   promoter     exon   intron 
# 96.84    89.13    90.48 
# 
# summary of distances to the nearest TSS:
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    2120    6279   14137   16750  282513 
  
# Extract distance to closest TSS merged sites
dist_tss_WK2 <- getAssociationWithTSS(WK_002_merged.anot)
head(dist_tss_WK2)
  
hist(dist_tss_WK2$dist.to.feature)
plotTargetAnnotation(WK_002_merged.anot, main = "RRBS & WGBS Merged Sites")
  
dist_tss_WK3 <- getAssociationWithTSS(WK_003_merged.anot)
head(dist_tss_WK3)
  
hist(dist_tss_WK3$dist.to.feature)
plotTargetAnnotation(dist_tss_WK3, main = "RRBS & WGBS Merged Sites")
  
# Extract distance to nearest TSS RRBS only sites
dist_tss2_WK2 <- getAssociationWithTSS(WK2RRBS_Diff.anot)
head(dist_tss2_WK2)
  
hist(dist_tss2_WK2$dist.to.feature)
plotTargetAnnotation(WK2RRBS_Diff.anot, main = "Unique to RRBS")
  
dist_tss2_WK3 <- getAssociationWithTSS(WK3RRBS_Diff.anot)
head(dist_tss2_WK3)
  
hist(dist_tss2_WK3$dist.to.feature)
plotTargetAnnotation(WK3RRBS_Diff.anot, main = "Unique to RRBS")
  
# Extract distance to nearest TSS WGBS only sites
dist_tss3_WK2 <- getAssociationWithTSS(WK2WGBS_Diff.anot)
head(dist_tss3_WK2)
  
hist(dist_tss3_WK2$dist.to.feature)
plotTargetAnnotation(WK2WGBS_Diff.anot, main = "Unique to WGBS")
  
dist_tss3_WK3 <- getAssociationWithTSS(WK3WGBS_Diff.anot)
head(dist_tss3_WK3)
  
hist(dist_tss3_WK3$dist.to.feature)
plotTargetAnnotation(WK3WGBS_Diff.anot, main = "Unique to WGBS")

# merge sites covered by both, RRBS, and WGBS, respectively
Both_dist_tss <- rbind(dist_tss_WK2, dist_tss_WK3)
Both_dist_tss$Method <- "Both"
Both_dist_tss$Population <- "Wik"
head(Both_dist_tss)

RRBS_dist_tss <- rbind(dist_tss2_WK2, dist_tss2_WK3)
RRBS_dist_tss$Method <- "RRBS"
RRBS_dist_tss$Population <- "Wik"

WGBS_dist_tss <- rbind(dist_tss3_WK2, dist_tss3_WK3)
WGBS_dist_tss$Method <- "WGBS"
WGBS_dist_tss$Population <- "Wik"
head(WGBS_dist_tss)

all_dist_tss <- rbind(Both_dist_tss, RRBS_dist_tss)
all_dist_tss <- rbind(all_dist_tss, WGBS_dist_tss)

WT_dist_tss <- read.csv("/home/ekerns/MethylMethods/Figures/WT_dist_tss.csv")
WT_dist_tss <- WT_dist_tss[,-1]

all_dist_tss <- rbind(all_dist_tss, WT_dist_tss)
head(all_dist_tss)


# Visualize annotations ---------------------------------------------------


Annots <- read_xlsx("/WGBSRRBSAnnotations.xlsx")

####### Figure 6

Annots_long <- Annots %>% 
  pivot_longer(cols = !c(Population, FishID, Sequencing), names_to= 'GenomicFeature', values_to= 'Percent')

Annots_long$GenomicFeature

Annots_CpG <- Annots_long[Annots_long$GenomicFeature == 'CpGI' |
                          Annots_long$GenomicFeature == 'CpGShores' |
                          Annots_long$GenomicFeature == 'NonCpGI',]

Annots_CpG_avg <- Annots_CpG %>%
  group_by(Population, Sequencing, GenomicFeature) %>%
  summarise(Percent = mean(Percent), .groups = 'drop')
avg <- Annots_CpG_avg %>%
  group_by(Sequencing, GenomicFeature) %>%
  summarise(
    avg_percent = mean(Percent),
    sd_percent = sd(Percent),
    .groups = "drop"
  )
avg
# Sequencing GenomicFeature avg_percent sd_percent
# <chr>      <chr>                <dbl>      <dbl>
# 1 Both       CpGI                  45.0     0.849 
# 2 Both       CpGShores             20.7     0.117 
# 3 Both       NonCpGI               34.3     0.732 
# 4 RRBS       CpGI                  48.8     0.371 
# 5 RRBS       CpGShores             19.9     0.0389
# 6 RRBS       NonCpGI               31.2     0.329 
# 7 WGBS       CpGI                  14.2     0.202 
# 8 WGBS       CpGShores             27.0     0.0212
# 9 WGBS       NonCpGI               58.8     0.184 

fill_colors <- c("CpGI" = "skyblue", "CpGShores"= "darkblue", "NonCpGI"= "darkgrey")

CpGIPlot <- ggplot(data = Annots_CpG_avg, aes(fill = GenomicFeature, y = Percent, x = Sequencing)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = fill_colors) + 
  scale_x_discrete(labels=c("Both" = "Both", "RRBS" = "RRBS Only", "WGBS" = "WGBS Only"))+
  xlab("Sequencing Method") +
  ylab("Percent Overlap") +
  facet_wrap(~Population)+ 
  #annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold") + theme(axis.text.x = element_text(angle = 45, hjust = 1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)) +
  labs(tag = "A")
CpGIPlot

Annots_CpG$GenomicFeature

NumCpGs <- Annots_long[Annots_long$GenomicFeature == 'NumCpGsForCpgI', ]

NumCpGsPlot<- ggplot(data = NumCpGs, aes(shape = Population, y=Percent, x=Sequencing)) + 
  geom_point(size = 3, position = position_jitter(width = 0.2, height = 0)) +
  scale_x_discrete(labels=c("Both" = "Both", "RRBS" = "RRBS Only", "WGBS" = "WGBS Only"))+
  ylab("Number of CpG Sites")+
  xlab("Sequencing Method")+
  theme_classic()+
  theme(legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
NumCpGsPlot

Annots_CpG_PercCov <- Annots_long[Annots_long$GenomicFeature == 'PercCpGICovered' |
                            Annots_long$GenomicFeature == 'PercShoresCovered',]

PercCpGIsPlot<- ggplot(data = Annots_CpG_PercCov, aes(color = GenomicFeature, shape = Population, y=Percent, x=Sequencing)) + 
  geom_point(size = 3, position = position_jitter(width = 0.2, height = 0)) +
  scale_x_discrete(labels=c("Both" = "Both", "RRBS" = "RRBS Only", "WGBS" = "WGBS Only"))+
  ylab("Percent Covered")+
  xlab("Sequencing Method")+
  scale_color_manual(name = "Genomic Feature", 
                       labels = c("CpG Islands", 
                                  "CpG Shores"),
                       values = c("skyblue", "darkblue"))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  ylim(0,100)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PercCpGIsPlot


CpGIPlot + NumCpGsPlot + PercCpGIsPlot



Annots_long$GenomicFeature

Annots_Genes <- Annots_long[Annots_long$GenomicFeature == 'PercPromoterOverlap' |
                            Annots_long$GenomicFeature == 'PercExonOverlap' |
                            Annots_long$GenomicFeature == 'PercIntronOverlap' |
                            Annots_long$GenomicFeature == 'PercIntergenicOverlap',]

Annots_Genes_avg <- Annots_Genes %>%
  group_by(Population, Sequencing, GenomicFeature) %>%
  summarise(Percent = mean(Percent), .groups = 'drop')
Annots_Genes_avg
avg2 <- Annots_Genes_avg %>%
  group_by(Sequencing, GenomicFeature) %>%
  summarise(
    avg_percent = mean(Percent),
    sd_percent = sd(Percent),
    .groups = "drop"
  )
avg2
# Sequencing GenomicFeature        avg_percent sd_percent
# <chr>      <fct>                       <dbl>      <dbl>
# 1 Both       PercPromoterOverlap          42.0     2.12  
# 2 Both       PercExonOverlap              16.7     1.17  
# 3 Both       PercIntronOverlap            24.4     0.778 
# 4 Both       PercIntergenicOverlap        16.8     0.180 
# 5 RRBS       PercPromoterOverlap          36.8     2.49  
# 6 RRBS       PercExonOverlap              19.9     1.73  
# 7 RRBS       PercIntronOverlap            25.7     0.668 
# 8 RRBS       PercIntergenicOverlap        17.6     0.0884
# 9 WGBS       PercPromoterOverlap          14.9     0.117 
# 10 WGBS       PercExonOverlap              15.9     0.0955
# 11 WGBS       PercIntronOverlap            43       0.0283
# 12 WGBS       PercIntergenicOverlap        26.2     0.0530

fill_colors <- c("PercPromoterOverlap" = "purple", "PercExonOverlap"= "skyblue", "PercIntronOverlap"= "darkblue", "PercIntergenicOverlap" = "darkgrey")

Annots_Genes_avg$GenomicFeature <- factor(Annots_Genes_avg$GenomicFeature, levels = c("PercPromoterOverlap", "PercExonOverlap", "PercIntronOverlap", "PercIntergenicOverlap"))

GenePlot <- ggplot(data = Annots_Genes_avg, aes(fill = GenomicFeature, y = Percent, x = Sequencing)) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_classic() +
  scale_x_discrete(labels=c("Both" = "Both", "RRBS" = "RRBS Only", "WGBS" = "WGBS Only"))+
  xlab("Sequencing Method") +
  ylab("Percent Overlap") +
  scale_fill_manual(name = "Genomic Feature", 
                    values = fill_colors,
                     labels = c("Promoters", 
                                "Exons",
                                "Introns",
                                "Intergenic Regions"))+
  #annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  facet_wrap(~Population)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 16, face = "bold"),
        plot.tag.position = c(0, 1)) +
  labs(tag = "D")
GenePlot


Genes_NumCpGs <- Annots_long[Annots_long$GenomicFeature == 'NumCpGs', ]

NumCsGenePlot<- ggplot(data = Genes_NumCpGs, aes(shape = Population, y=Percent, x=Sequencing)) + 
  geom_point(size = 3, position = position_jitter(width = 0.2, height = 0)) +
  scale_x_discrete(labels=c("Both" = "Both", "RRBS" = "RRBS Only", "WGBS" = "WGBS Only"))+
  ylab("Number of CpG Sites")+
  xlab("Sequencing Method")+
  theme_classic()+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
NumCsGenePlot

Annots_Genes_PercCov <- Annots_long[Annots_long$GenomicFeature == 'PromoterCovered' |
                                    Annots_long$GenomicFeature == 'ExonCovered' |
                                    Annots_long$GenomicFeature == 'IntronCovered',]

PercGenesPlot<- ggplot(data = Annots_Genes_PercCov, aes(color = GenomicFeature, shape = Population, y=Percent, x=Sequencing)) + 
  geom_point(size = 3, position = position_jitter(width = 0.2, height = 0)) +
  scale_x_discrete(labels=c("Both" = "Both", "RRBS" = "RRBS Only", "WGBS" = "WGBS Only"))+
  ylab("Percent Covered")+
  xlab("Sequencing Method")+
  scale_color_manual(name = "Genomic Feature", 
                     labels = c("Exons", 
                                "Introns",
                                "Promoters"),
                     values = c("skyblue", "darkblue", "purple"))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -.5, vjust = 1, size = 5, fontface = "bold")+
  ylim(0,100)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PercGenesPlot


GenePlot + NumCsGenePlot + PercGenesPlot

(CpGIPlot + NumCpGsPlot + PercCpGIsPlot) / (GenePlot + NumCsGenePlot + PercGenesPlot)

(CpGIPlot + GenePlot) / (NumCpGsPlot + NumCsGenePlot) / (PercCpGIsPlot + PercGenesPlot)

