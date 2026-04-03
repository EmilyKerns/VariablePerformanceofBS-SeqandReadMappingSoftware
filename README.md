# VariablePerformanceofBS-SeqandReadMappingSoftware

These scripts accomppany the article titled "Variable performance of widely used bisulfite sequencing methods and read mapping software for DNA methylation." (https://doi.org/10.1101/2025.03.14.643302). The yml files can be used to download all of the necessary conda packages for read mapping and methylation calling. 
They can be used to analyze WGBS or RRBS data. 

**Analysis pipeline:**
1. Adapter trimming & FastQC: Trim Galore.
2. Index Reference Genome & Read Mapping: You can use either BWA meth, Biscuit, BisulfiteBolt (BSBolt), or Bismark for read mapping. Index the reference before proceeding with mapping your sample reads.
3. Call methylation: Scripts for calling methylation with Bismark are included. Methylation calling for all other read mapping methods was done using MethylDackel.

Some analysis was completed on the HTCondor Software Suite for Highthroughput Computing. Def & submission files are also included if you have access to this computing system. All scripts for statistical analyses are compatible with R. Raw seqeunce data is available on Dryad and NCBI SRA BioProject.
