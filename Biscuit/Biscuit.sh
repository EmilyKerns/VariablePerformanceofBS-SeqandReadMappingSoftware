#!/bin/bash

export PATH=/opt/conda/bin:$PATH

# prior to running this script, move the indexed genome files to one folder and create a tar.gz 
tar -xvzf biscuit_ref.tar.gz --strip-components=1

fq=$1

echo "Aligning WGBS reads with Biscuit on $(hostname) at $(date)"
echo "  Processing sample: ${fq}"
echo "  Running Biscuit with 16 threads..."

# Replace $REFERENCE_GENOME with your reference genome file and adjust -@ 16 to however many threads you want to use for alignment
biscuit align -@ 16 $REFERENCE_GENOME ${fq}_1_val_1.fq.gz ${fq}_2_val_2.fq.gz |
        /opt/conda/bin/samtools sort -o ${fq}_biscuit.bam -O BAM


echo "Read alignment completed at $(date)"
