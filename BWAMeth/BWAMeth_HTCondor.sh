#!/bin/bash

base_id=$1

echo "Aligning RRBS reads with BWA meth on $(hostname) at $(date)"
echo "  Processing sample: ${base_id}"
echo "  Running BWA meth with 8 threads..."

# Make sure reference file exists and show its timestamp
ls -l RFERENCE_GENOME.fna

# Update timestamps on index files FIRST, with a small delay. Avoids error where the computer thinks the reference genome is newer than the indexed reference genome
echo "  Updating index file timestamps..."
sleep 2
touch REFERENCE_GENOME.fna.bwameth.c2t*
ls -l REFERENCE_GENOME.fna.bwameth.c2t*

bwameth.py --reference REFERENCE_GENOME.fna --threads 8 ${base_id}_trimmed.fq.gz | \
        samtools sort --threads 6 -m 2G -o ${base_id}.bam -

echo "Read alignment completed at $(date)"
