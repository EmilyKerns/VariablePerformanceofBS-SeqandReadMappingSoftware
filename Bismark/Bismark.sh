#!/bin/bash

# Before running this script, create a tarball of the reference genome indexed with Bismark
tar -xvzf REFERENCE_Bisulfite_Genome.tar.gz

fq=$1

echo "Aligning BS-seq reads with Bismark on $(hostname) at $(date)"
echo "  Processing sample: ${base_id}"
echo "  Running Bismark..."

# replace ./ with the path to your indexed reference genome. For global alignment, remove the --local flag
bismark --non_directional --local --genome ./ -1 ${fq}_1_val_1.fq.gz -2 ${fq}_2_val_2.fq.gz

echo "Read alignment completed at $(date)"
