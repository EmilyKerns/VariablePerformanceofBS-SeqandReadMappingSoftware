#!/bin/bash

for sample_id in *_bismark_bt2_pe.bam
 do
        base_id=$(basename $sample_id "_bismark_bt2_pe.bam")
        bismark2bedGraph --ucsc -o ${base_id} ${base_id}_bismark_bt2_pe.bam
 done
