### Assess overlapping CpG sites between sequencing methods

# load packages
library(methylKit)
library(ggvenn)

#### Min Depth 10x

# import samples
WGBSfile.list=list("minDepth10/NoSNPFiltering/WGBS_CG_22_WK_002_CpG.methylKit", "minDepth10/NoSNPFiltering/WGBS_CG_22_WK_003_CpG.methylKit", "minDepth10/NoSNPFiltering/CG_22_WK_002_CpG.methylKit", "minDepth10/NoSNPFiltering/CG_22_WK_003_CpG.methylKit")

WGBSvsRRBS=methRead(WGBSfile.list,
                    sample.id=list("WGBS_CG_WK_002", "WGBS_CG_WK_003", "RRBS_CG_WK_002", "RRBS_CG_WK_003"),
                    pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                    assembly="StickleGeneAnnotations",
                    treatment=c(1,1,2,2),
                    context="CpG",
                    mincov = 10
)

# Separate each sample into individual methylRaw objects - contains info about mehtylation and chr location
WGBS_CG_WK_002=WGBSvsRRBS[[1]]
WGBS_CG_WK_003=WGBSvsRRBS[[2]]
RRBS_CG_WK_002=WGBSvsRRBS[[3]]
RRBS_CG_WK_003=WGBSvsRRBS[[4]]

chrs_WGBS_CG_WK_002=WGBS_CG_WK_002[[1]]
chrs_WGBS_CG_WK_003=WGBS_CG_WK_003[[1]]
chrs_RRBS_CG_WK_002=RRBS_CG_WK_002[[1]]
chrs_RRBS_CG_WK_003=RRBS_CG_WK_003[[1]]

unions <- list(WGBS_CG_WK_002 = chrs_WGBS_CG_WK_002, WGBS_CG_WK_003 = chrs_WGBS_CG_WK_003, RRBS_CG_WK_002 = chrs_RRBS_CG_WK_002, RRBS_CG_WK_003 = chrs_RRBS_CG_WK_003)
ggvenn(unions)

# Overlapping CpG sites on biological replicates with RRBS and WGBS
WKWGBS_unions <- list(WGBS_CG_WK_002 = chrs_WGBS_CG_WK_002, WGBS_CG_WK_003 = chrs_WGBS_CG_WK_003)
ggvenn(WKWGBS_unions)

WKRRBS_unions <- list(RRBS_CG_WK_002 = chrs_RRBS_CG_WK_002, RRBS_CG_WK_003 = chrs_RRBS_CG_WK_003)
ggvenn(WKRRBS_unions)

# Overlapping CpG sites with RRBS and WGBS on technical replicates
WK_002_Union <- list(WGBS_CG_WK_002 = chrs_WGBS_CG_WK_002, RRBS_CG_WK_002 = chrs_RRBS_CG_WK_002)
ggvenn(WK_002_Union)

WK_003_Union <- list(WGBS_CG_WK_003 = chrs_WGBS_CG_WK_003, RRBS_CG_WK_003 = chrs_RRBS_CG_WK_003)
ggvenn(WK_003_Union)

#### Min Depth 5x

# import samples
WGBSfile.list=list("minDepth5/NoSNPFiltering/WGBS_CG_22_WK_002_CpG.methylKit", "minDepth5/NoSNPFiltering/WGBS_CG_22_WK_003_CpG.methylKit", "minDepth5/NoSNPFiltering/CG_22_WK_002_CpG.methylKit", "minDepth5/NoSNPFiltering/CG_22_WK_003_CpG.methylKit")

WGBSvsRRBS=methRead(WGBSfile.list,
                    sample.id=list("WGBS_CG_WK_002", "WGBS_CG_WK_003", "RRBS_CG_WK_002", "RRBS_CG_WK_003"),
                    pipeline = list(fraction = FALSE, chr.col = 1, start.col = 3, end.col = 3, coverage.col = 5, strand.col = 4, freqC.col=6),
                    assembly="StickleGeneAnnotations",
                    treatment=c(1,1,2,2),
                    context="CpG",
                    mincov = 5
)

# Separate each sample into individual methylRaw objects - contains info about mehtylation and chr location
WGBS_CG_WK_002=WGBSvsRRBS[[1]]
WGBS_CG_WK_003=WGBSvsRRBS[[2]]
RRBS_CG_WK_002=WGBSvsRRBS[[3]]
RRBS_CG_WK_003=WGBSvsRRBS[[4]]

chrs_WGBS_CG_WK_002=WGBS_CG_WK_002[[1]]
chrs_WGBS_CG_WK_003=WGBS_CG_WK_003[[1]]
chrs_RRBS_CG_WK_002=RRBS_CG_WK_002[[1]]
chrs_RRBS_CG_WK_003=RRBS_CG_WK_003[[1]]

unions <- list(WGBS_CG_WK_002 = chrs_WGBS_CG_WK_002, WGBS_CG_WK_003 = chrs_WGBS_CG_WK_003, RRBS_CG_WK_002 = chrs_RRBS_CG_WK_002, RRBS_CG_WK_003 = chrs_RRBS_CG_WK_003)

ggvenn(unions)

# Overlapping CpG sites on biological replicates with RRBS and WGBS
WKWGBS_unions <- list(WGBS_CG_WK_002 = chrs_WGBS_CG_WK_002, WGBS_CG_WK_003 = chrs_WGBS_CG_WK_003)
ggvenn(WKWGBS_unions)

WKRRBS_unions <- list(RRBS_CG_WK_002 = chrs_RRBS_CG_WK_002, RRBS_CG_WK_003 = chrs_RRBS_CG_WK_003)
ggvenn(WKRRBS_unions)

# Overlapping CpG sites with RRBS and WGBS on technical replicates
WK_002_Union <- list(WGBS_CG_WK_002 = chrs_WGBS_CG_WK_002, RRBS_CG_WK_002 = chrs_RRBS_CG_WK_002)
ggvenn(WK_002_Union)

WK_003_Union <- list(WGBS_CG_WK_003 = chrs_WGBS_CG_WK_003, RRBS_CG_WK_003 = chrs_RRBS_CG_WK_003)
ggvenn(WK_003_Union)