#!/usr/bin/bash

# align each read to the reference genome and output a sorted .bam file
sample_id=${1}_R1_val_1.fq.gz
        base_id=$(basename $sample_id "_R1_val_1.fq.gz")
        bwa mem stickleback_v5_assembly_1.fa ${base_id}_R1_val_1.fq.gz ${base_id}_R2_val_2.fq.gz | samtools view -b | samtools sort > ${base_id}.bam
