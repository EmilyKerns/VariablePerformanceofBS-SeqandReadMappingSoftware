#!/bin/bash

#Using Trim Galore! for QC (Fastq) and adapater trimming (cutadapt) with -rrbs -paired
for sample_id in *_1.fq.gz
 do
        base_id=$(basename $sample_id "_1.fq.gz")
        trim_galore --paired --rrbs --fastqc -q 20 \
        ${base_id}_1.fq.gz ${base_id}_2.fq.gz
 done
