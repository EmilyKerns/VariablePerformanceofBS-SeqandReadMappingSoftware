#!/bin/bash

#run MethylDackel with a min depth of 10, --maxVariantFrac --minOppositeDepth to differentiate between SNPs and converted cytosines, --keepDupes for RRBS. Min depth of 5 or 10. Also tested --maxVariantFrac 0.25, 0.5, 0.75, and 0.8 with --minOppositeDepth 5, or with no SNP calling at all. 
for sample_id in *.bam
 do
        base_id=$(basename $sample_id ".bam")
	MethylDackel extract --minDepth 10 --methylKit --keepDupes --maxVariantFrac 0.8 --minOppositeDepth 5 stickleback_v5_assembly.fa ${base_id}.bam
 done
