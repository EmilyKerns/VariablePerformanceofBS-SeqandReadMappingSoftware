# Figure 4

# Biscuit -----------------------------------------------------------------

setwd("/Biscuit/Stickle/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WK_002_R1_CpG.methylKit")

Biscuit=methRead(file.list.10,
                  sample.id=list("CG_22_WK_002"),
                  pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                  assembly="StickleGeneAnnotations",
                  treatment=c(1),
                  context="CpG",
                  mincov = 10
)

# BSBolt -----------------------------------------------------------------

setwd("/BSBolt/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WK_002_output_sorted_CpG.methylKit")

BSB=methRead(file.list.10,
                  sample.id=list("CG_22_WK_002"),
                  pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                  assembly="StickleGeneAnnotations",
                  treatment=c(1),
                  context="CpG",
                  mincov = 10
)


# Bismark Local -----------------------------------------------------------


setwd("/Bismark/Local")

# Set file.list of samples
file.list.10=list("CG_22_WK_002_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
BL <- methRead(file.list.10,
                     sample.id=list("CG_WK_002"),
                     assembly="StickleGeneAnnotations",
                     treatment = c(1),
                     pipeline = "bismarkCoverage",
                     mincov = 10)


# Bismark -----------------------------------------------------------------



setwd("/Bismark/NondirectionalFlag/cov")

# Set file.list of samples
file.list.10=list("CG_22_WK_002_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
Bism <- methRead(file.list.10,
                     sample.id=list("CG_WK_002"),
                     assembly="StickleGeneAnnotations",
                     treatment = c(1),
                     pipeline = "bismarkCoverage",
                     mincov = 10)


# BWA Meth ----------------------------------------------------------------

setwd("/MethylMethods")

# import samples
file.list=list("CG_22_WK_002_CpG.methylKit")
  
BWAMeth=methRead(file.list,
                      sample.id=list("RRBS_CG_WK_002"),
                      pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                      assembly="StickleGeneAnnotations",
                      treatment=c(1),
                      context="CpG",
                      mincov = 10
  )

getMethylationStats(BWAMeth[[1]], plot=FALSE, both.strands=FALSE)


# Extract Percent Methylation -------------------------------------------------

Bisc=Biscuit[[1]]
Bolt=BSB[[1]]
BismL=BL[[1]]
Bismark=Bism[[1]]
BWA=BWAMeth[[1]]

Bisc$PercMeth <- ((Bisc$numCs/Bisc$coverage)*100)
Bolt$PercMeth <- ((Bolt$numCs/Bolt$coverage)*100)
BismL$PercMeth <- ((BismL$numCs/BismL$coverage)*100)
Bismark$PercMeth <- ((Bismark$numCs/Bismark$coverage)*100)
BWA$PercMeth <- ((BWA$numCs/BWA$coverage)*100)

Bisc <- getData(Bisc)
Bolt <- getData(Bolt)
BismL <- getData(BismL)
Bismark <- getData(Bismark)
BWA <- getData(BWA)

nrow(Bisc)    #615,994
nrow(Bolt)    #549,976
nrow(BismL)   #598,238
nrow(Bismark) #476,717
nrow(BWA)     #188,608

Bisc_summary <- Bisc %>%
  summarise(mean_depth = mean(coverage), 
            n_sites = n(), 
            sd_depth = sd(coverage))
Bisc_summary
# mean_depth n_sites sd_depth
# 1   16.49521  615994 19.46032

Bolt_summary <- Bolt %>%
  summarise(mean_depth = mean(coverage), 
            n_sites = n(), 
            sd_depth = sd(coverage))
Bolt_summary
# mean_depth n_sites sd_depth
# 1   15.89105  549976 20.53162

BismL_summary <- BismL %>%
  summarise(mean_depth = mean(coverage), 
            n_sites = n(), 
            sd_depth = sd(coverage))
BismL_summary
# mean_depth n_sites sd_depth
# 1   16.26131  598238  21.8889

Bismark_summary <- Bismark %>%
  summarise(mean_depth = mean(coverage), 
            n_sites = n(), 
            sd_depth = sd(coverage))
Bismark_summary
# mean_depth n_sites sd_depth
# 1   16.11017  476717 21.48056

BWA_summary <- BWA %>%
  summarise(mean_depth = mean(coverage), 
            n_sites = n(), 
            sd_depth = sd(coverage))
BWA_summary
# mean_depth n_sites sd_depth
# 1   14.04099  188608 27.61839

# Percent Methylation Histograms ------------------------------------------

BiscuitPlot <- ggplot(Bisc)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BiscuitPlot

BoltPlot <- ggplot(Bolt)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 16), axis.title.x = element_text(size=18, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BoltPlot

BismLPlot <- ggplot(BismL)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismLPlot

BismarkPlot <- ggplot(Bismark)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismarkPlot

BWAPlot <- ggplot(BWA)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BWAPlot

BWAPlot + BismarkPlot + BismLPlot + BiscuitPlot + BoltPlot


# WK_003 ------------------------------------------------------------------

# Supplemental Figure 4

# Biscuit -----------------------------------------------------------------

setwd("/Biscuit/Stickle/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WK_003_R1_CpG.methylKit")

Biscuit=methRead(file.list.10,
                 sample.id=list("CG_22_WK_003"),
                 pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                 assembly="StickleGeneAnnotations",
                 treatment=c(1),
                 context="CpG",
                 mincov = 10
)

# BSBolt -----------------------------------------------------------------

setwd("/BSBolt/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WK_003_output_sorted_CpG.methylKit")

BSB=methRead(file.list.10,
             sample.id=list("CG_22_WK_003"),
             pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
             assembly="StickleGeneAnnotations",
             treatment=c(1),
             context="CpG",
             mincov = 10
)


# Bismark Local -----------------------------------------------------------


setwd("/Bismark/Local")

# Set file.list of samples
file.list.10=list("CG_22_WK_003_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
BL <- methRead(file.list.10,
               sample.id=list("CG_WK_003"),
               assembly="StickleGeneAnnotations",
               treatment = c(1),
               pipeline = "bismarkCoverage",
               mincov = 10)


# Bismark -----------------------------------------------------------------



setwd("/Bismark/NondirectionalFlag/cov")

# Set file.list of samples
file.list.10=list("CG_22_WK_003_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
Bism <- methRead(file.list.10,
                 sample.id=list("CG_WK_003"),
                 assembly="StickleGeneAnnotations",
                 treatment = c(1),
                 pipeline = "bismarkCoverage",
                 mincov = 10)


# BWA Meth ----------------------------------------------------------------

setwd("/MethylMethods")

# import samples
file.list=list("CG_22_WK_003_CpG.methylKit")

BWAMeth=methRead(file.list,
                 sample.id=list("RRBS_CG_WK_003"),
                 pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                 assembly="StickleGeneAnnotations",
                 treatment=c(1),
                 context="CpG",
                 mincov = 10
)

getMethylationStats(BWAMeth[[1]], plot=FALSE, both.strands=FALSE)


# Extract Percent Methylation -------------------------------------------------

Bisc=Biscuit[[1]]
Bolt=BSB[[1]]
BismL=BL[[1]]
Bismark=Bism[[1]]
BWA=BWAMeth[[1]]

Bisc$PercMeth <- ((Bisc$numCs/Bisc$coverage)*100)
Bolt$PercMeth <- ((Bolt$numCs/Bolt$coverage)*100)
BismL$PercMeth <- ((BismL$numCs/BismL$coverage)*100)
Bismark$PercMeth <- ((Bismark$numCs/Bismark$coverage)*100)
BWA$PercMeth <- ((BWA$numCs/BWA$coverage)*100)

Bisc <- getData(Bisc)
Bolt <- getData(Bolt)
BismL <- getData(BismL)
Bismark <- getData(Bismark)
BWA <- getData(BWA)

nrow(Bisc)    #487,314
nrow(Bolt)    #428,240
nrow(BismL)   #473,743
nrow(Bismark) #384,860
nrow(BWA)     #152,023


# Percent Methylation Histograms ------------------------------------------

BiscuitPlot <- ggplot(Bisc)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BiscuitPlot

BoltPlot <- ggplot(Bolt)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 16), axis.title.x = element_text(size=18, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BoltPlot

BismLPlot <- ggplot(BismL)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismLPlot

BismarkPlot <- ggplot(Bismark)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismarkPlot

BWAPlot <- ggplot(BWA)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BWAPlot

BWAPlot + BismarkPlot + BismLPlot + BiscuitPlot + BoltPlot


# WT_002 ------------------------------------------------------------------

# Supplemental Figure 5

# Biscuit -----------------------------------------------------------------

setwd("/Biscuit/Stickle/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WT_002_R1_CpG.methylKit")

Biscuit=methRead(file.list.10,
                 sample.id=list("CG_22_WT_002"),
                 pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                 assembly="StickleGeneAnnotations",
                 treatment=c(1),
                 context="CpG",
                 mincov = 10
)

# BSBolt -----------------------------------------------------------------

setwd("/BSBolt/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WT_002_output_sorted_CpG.methylKit")

BSB=methRead(file.list.10,
             sample.id=list("CG_22_WT_002"),
             pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
             assembly="StickleGeneAnnotations",
             treatment=c(1),
             context="CpG",
             mincov = 10
)


# Bismark Local -----------------------------------------------------------


setwd("/Bismark/Local")

# Set file.list of samples
file.list.10=list("CG_22_WT_002_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
BL <- methRead(file.list.10,
               sample.id=list("CG_WT_002"),
               assembly="StickleGeneAnnotations",
               treatment = c(1),
               pipeline = "bismarkCoverage",
               mincov = 10)


# Bismark -----------------------------------------------------------------



setwd("/Bismark/NondirectionalFlag/cov")

# Set file.list of samples
file.list.10=list("CG_22_WT_002_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
Bism <- methRead(file.list.10,
                 sample.id=list("CG_WT_002"),
                 assembly="StickleGeneAnnotations",
                 treatment = c(1),
                 pipeline = "bismarkCoverage",
                 mincov = 10)


# BWA Meth ----------------------------------------------------------------

setwd("/MethylMethods")

# import samples
file.list=list("CG_22_WT_002_CpG.methylKit")

BWAMeth=methRead(file.list,
                 sample.id=list("RRBS_CG_WT_002"),
                 pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                 assembly="StickleGeneAnnotations",
                 treatment=c(1),
                 context="CpG",
                 mincov = 10
)

getMethylationStats(BWAMeth[[1]], plot=FALSE, both.strands=FALSE)


# Extract Percent Methylation -------------------------------------------------

Bisc=Biscuit[[1]]
Bolt=BSB[[1]]
BismL=BL[[1]]
Bismark=Bism[[1]]
BWA=BWAMeth[[1]]

Bisc$PercMeth <- ((Bisc$numCs/Bisc$coverage)*100)
Bolt$PercMeth <- ((Bolt$numCs/Bolt$coverage)*100)
BismL$PercMeth <- ((BismL$numCs/BismL$coverage)*100)
Bismark$PercMeth <- ((Bismark$numCs/Bismark$coverage)*100)
BWA$PercMeth <- ((BWA$numCs/BWA$coverage)*100)

Bisc <- getData(Bisc)
Bolt <- getData(Bolt)
BismL <- getData(BismL)
Bismark <- getData(Bismark)
BWA <- getData(BWA)

nrow(Bisc)    #376,159
nrow(Bolt)    #310,700
nrow(BismL)   #310,700
nrow(Bismark) #384,860
nrow(BWA)     #60,313


# Percent Methylation Histograms ------------------------------------------

BiscuitPlot <- ggplot(Bisc)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BiscuitPlot

BoltPlot <- ggplot(Bolt)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 16), axis.title.x = element_text(size=18, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BoltPlot

BismLPlot <- ggplot(BismL)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismLPlot

BismarkPlot <- ggplot(Bismark)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismarkPlot

BWAPlot <- ggplot(BWA)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BWAPlot

BWAPlot + BismarkPlot + BismLPlot + BiscuitPlot + BoltPlot


# WT_004 ------------------------------------------------------------------

# Supplemental Figure 6

# Biscuit -----------------------------------------------------------------

setwd("/Biscuit/Stickle/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WT_004_R1_CpG.methylKit")

Biscuit=methRead(file.list.10,
                 sample.id=list("CG_22_WT_004"),
                 pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                 assembly="StickleGeneAnnotations",
                 treatment=c(1),
                 context="CpG",
                 mincov = 10
)

# BSBolt -----------------------------------------------------------------

setwd("/BSBolt/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_WT_004_output_sorted_CpG.methylKit")

BSB=methRead(file.list.10,
             sample.id=list("CG_22_WT_004"),
             pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
             assembly="StickleGeneAnnotations",
             treatment=c(1),
             context="CpG",
             mincov = 10
)


# Bismark Local -----------------------------------------------------------


setwd("/Bismark/Local")

# Set file.list of samples
file.list.10=list("CG_22_WT_004_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
BL <- methRead(file.list.10,
               sample.id=list("CG_WT_004"),
               assembly="StickleGeneAnnotations",
               treatment = c(1),
               pipeline = "bismarkCoverage",
               mincov = 10)


# Bismark -----------------------------------------------------------------



setwd("/Bismark/NondirectionalFlag")

# Set file.list of samples
file.list.10=list("CG_22_WT_002_R1_val_1_bismark_bt2_pe.bismark.cov")

# read the files to a methylRawList object: myobj, min coverage of 10
Bism <- methRead(file.list.10,
                 sample.id=list("CG_WT_004"),
                 assembly="StickleGeneAnnotations",
                 treatment = c(1),
                 pipeline = "bismarkCoverage",
                 mincov = 10)


# BWA Meth ----------------------------------------------------------------

setwd("/MethylMethods")

# import samples
file.list=list("CG_22_WT_004_CpG.methylKit")

BWAMeth=methRead(file.list,
                 sample.id=list("RRBS_CG_WT_004"),
                 pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                 assembly="StickleGeneAnnotations",
                 treatment=c(1),
                 context="CpG",
                 mincov = 10
)

getMethylationStats(BWAMeth[[1]], plot=FALSE, both.strands=FALSE)


# Extract Percent Methylation -------------------------------------------------

Bisc=Biscuit[[1]]
Bolt=BSB[[1]]
BismL=BL[[1]]
Bismark=Bism[[1]]
BWA=BWAMeth[[1]]

Bisc$PercMeth <- ((Bisc$numCs/Bisc$coverage)*100)
Bolt$PercMeth <- ((Bolt$numCs/Bolt$coverage)*100)
BismL$PercMeth <- ((BismL$numCs/BismL$coverage)*100)
Bismark$PercMeth <- ((Bismark$numCs/Bismark$coverage)*100)
BWA$PercMeth <- ((BWA$numCs/BWA$coverage)*100)

Bisc <- getData(Bisc)
Bolt <- getData(Bolt)
BismL <- getData(BismL)
Bismark <- getData(Bismark)
BWA <- getData(BWA)

nrow(Bisc)    #495,536
nrow(Bolt)    #433,616
nrow(BismL)   #481,172
nrow(Bismark) #266,950
nrow(BWA)     #132,541


# Percent Methylation Histograms ------------------------------------------

BiscuitPlot <- ggplot(Bisc)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BiscuitPlot

BoltPlot <- ggplot(Bolt)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 16), axis.title.x = element_text(size=18, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("Percent Methylation per Base") +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BoltPlot

BismLPlot <- ggplot(BismL)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismLPlot

BismarkPlot <- ggplot(Bismark)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))+
  ylab("")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BismarkPlot

BWAPlot <- ggplot(BWA)+
  geom_histogram(aes(x = PercMeth),fill = "skyblue", alpha = .8, color = "black", binwidth = 5)+
  #geom_histogram(aes(x = RRBS_PercMeth, fill = "RRBS"), color = "black", alpha = .8)+
  #scale_fill_manual(values = c("WGBS" = "skyblue", "RRBS" = "darkblue"), name = "Method") +
  theme_classic() +
  theme(text = element_text(size = 14), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=18, face="bold", colour = "black"))+
  ylab("Frequency")+
  xlab("") +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -.5, vjust = 1, size = 5, fontface = "bold")
BWAPlot

BWAPlot + BismarkPlot + BismLPlot + BiscuitPlot + BoltPlot

