# decompress indexed reference genome
tar -xvzf BSBolt_ref.tar.gz

fq=$1

# Bisulfite bolt alignment for paired end reads
bsbolt Align -t 10 -F1 ${fq}_1_val_1.fq.gz -F2 ${fq}_2_val_2.fq.gz -DB ./BSBolt_ref/BSB_ref.fa -O ${fq}_BSBolt

