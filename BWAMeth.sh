#!/bin/bash

for sample_id in *_val_1.fq.gz
 do
        base_id=$(basename $sample_id "_val_1.fq.gz")
        /PATH/bwameth.py --reference /PATH/stickleback_v5_assembly.fa ${base_id}_val_1.fq.gz ${base_id}_val_2.fq.gz | /PATH/samtools view -b | /PATH/samtools sort --output-fmt BAM > ${base_id}.bam
 done
