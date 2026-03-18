# BSBolt -----------------------------------------------------------------

setwd("/BSBolt/MethylDackel")

# List of methylation files. Min depth 5, SNP filtered with Biscuit
file.list.10 <- list("CG_22_GOS_001_output_sorted_CpG.methylKit", "CG_22_GOS_002_output_sorted_CpG.methylKit", "CG_GOS_35_04_output_sorted_CpG.methylKit", "CG_GOS_2210_01_output_sorted_CpG.methylKit", "FD_GOS_27_output_sorted_CpG.methylKit", "FD_GOS_28_output_sorted_CpG.methylKit", "FD_GOS_33_output_sorted_CpG.methylKit", "CG_22_ROB_001_output_sorted_CpG.methylKit","CG_22_ROB_002_output_sorted_CpG.methylKit", "CG_ROB_003_output_sorted_CpG.methylKit", "CG_ROB_31_01_output_sorted_CpG.methylKit", "FD_ROB_98_output_sorted_CpG.methylKit", "FD_ROB_102_output_sorted_CpG.methylKit", "FD_ROB_106_output_sorted_CpG.methylKit", "CG_22_SAY_001_output_sorted_CpG.methylKit", "CG_22_SAY_002_output_sorted_CpG.methylKit", "CG_SAY_4_01_output_sorted_CpG.methylKit", "CG_SAY_2343_003_output_sorted_CpG.methylKit", "CG_22_WK_001_output_sorted_CpG.methylKit", "CG_22_WK_002_output_sorted_CpG.methylKit", "CG_22_WK_003_output_sorted_CpG.methylKit", "CG_22_WK_004_output_sorted_CpG.methylKit", "FD_22_WK_002_output_sorted_CpG.methylKit", "FD_22_WK_005_output_sorted_CpG.methylKit", "FD_22_WK_011_output_sorted_CpG.methylKit", "FD_22_WK_012_output_sorted_CpG.methylKit", "CG_22_WT_001_output_sorted_CpG.methylKit", "CG_22_WT_002_output_sorted_CpG.methylKit", "CG_22_WT_003_output_sorted_CpG.methylKit", "CG_22_WT_004_output_sorted_CpG.methylKit", "FD_22_WT_004_output_sorted_CpG.methylKit", "FD_22_WT_007_output_sorted_CpG.methylKit", "FD_22_WT_016_output_sorted_CpG.methylKit", "FD_22_WT_017_output_sorted_CpG.methylKit")

myobj.10=methRead(file.list.10,
                  sample.id=list("CG_22_GOS_001", "CG_22_GOS_002", "CG_GOS_35_04","CG_GOS_2210","FD_GOS_27", "FD_GOS_28", "FD_GOS_33", "CG_22_ROB_001","CG_22_ROB_002","CG_ROB_003","CG_ROB_31","FD_ROB_98","FD_ROB_102","FD_ROB_106","CG_22_SAY_001","CG_22_SAY_002","CG_SAY_4","CG_SAY_2343","CG_22_WK_001","CG_22_WK_002", "CG_22_WK_003", "CG_22_WK_004", "FD_22_WK_002", "FD_22_WK_005", "FD_22_WK_011","FD_22_WK_012","CG_22_WT_001", "CG_22_WT_002", "CG_22_WT_003", "CG_22_WT_004", "FD_22_WT_004", "FD_22_WT_007", "FD_22_WT_016", "FD_22_WT_017"),
                  pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                  assembly="StickleGeneAnnotations",
                  treatment=c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5),
                  context="CpG",
                  mincov = 10
)

# Separate each sample into individual methylRaw objects - contains info about mehtylation and chr location

getMethylationStats(myobj.10[[1]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[2]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[3]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[4]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[5]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[6]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[7]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[8]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[9]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[10]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[11]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[12]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[13]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[14]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[15]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[16]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[17]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[18]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[19]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[20]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[21]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[22]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[23]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[24]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[25]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[26]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[27]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[28]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[29]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[30]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[31]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[32]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[33]], plot=FALSE, both.strands=FALSE)
getMethylationStats(myobj.10[[34]], plot=FALSE, both.strands=FALSE)


