#!/bin/bash

for sample_id in *_bismark_bt2_pe.bam
 do
	base_id=$(basename $sample_id "_bismark_bt2_pe.bam")
	bismark_methylation_extractor --comprehensive --merge_non_CpG --gzip --bedGraph ${base_id}_bismark_bt2_pe.bam
 done
