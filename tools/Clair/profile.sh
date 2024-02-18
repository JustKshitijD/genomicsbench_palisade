#!/bin/bash
# Download and extract the testing dataset
wget 'http://www.bio8.cs.hku.hk/testingData.tar'
tar -xf testingData.tar

# Download the Illumina model
wget http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar
tar -xf 12345.tar

mkdir training

# Use cProfile to find the most time consuming functions
python3 -m cProfile -o cProfile.txt clair.py callVarBam --chkpnt_fn ./model --bam_fn testingData/chr21/chr21.bam --ref_fn testingData/chr21/chr21.fa --call_fn ./training/chr21.vcf --sampleName HG001 --pysam_for_all_indel_bases --threads 8 --qual 100 --ctgName chr21 --ctgStart 10269870 --ctgEnd 46672937 --pipe_line

# Store intermediate data with --store_loaded_mini_match
python3 clair.py callVarBam --chkpnt_fn ./model --bam_fn testingData/chr21/chr21.bam --ref_fn testingData/chr21/chr21.fa --call_fn ./training/chr21.vcf --sampleName HG001 --pysam_for_all_indel_bases --threads 8 --qual 100 --ctgName chr21 --ctgStart 10269870 --ctgEnd 46672937 --pipe_line --store_loaded_mini_match

# Only run prediction function in clair
python3 clair.py callVarBam --chkpnt_fn ./model --bam_fn testingData/chr21/chr21.bam --ref_fn testingData/chr21/chr21.fa --call_fn ./training/chr21.vcf --sampleName HG001 --pysam_for_all_indel_bases --threads 1 --qual 100 --ctgName chr21 --ctgStart 10269870 --ctgEnd 46672937 --only_prediction --time_counter_file_name time_counter_non_pipeline_1_thread.h5

python3 read_h5.py time_counter_non_pipeline_1_thread.h5

# Run my prediction function
python3 my_prediction.py --chkpnt_fn ./model --bam_fn testingData/chr21/chr21.bam --ref_fn testingData/chr21/chr21.fa --call_fn ./training/chr21.vcf --sampleName HG001 --pysam_for_all_indel_bases --threads 1 --qual 100

python3 read_h5.py time_counter_only_prediction.h5

python3 clair.py callVarBam --chkpnt_fn ./model --bam_fn /z/scratch7/yufenggu/dataset/HG002_prom_R941_guppy360_2_GRCh38_TEST.bam.bai --ref_fn /z/scratch7/yufenggu/dataset/GRCh38_no_alt_chr20.fa  --call_fn ./training/chr20.vcf --sampleName HG002 --pysam_for_all_indel_bases --threads 8 --qual 100 --ctgName chr20 --ctgStart 10269870 --ctgEnd 46672937 --pipe_line --store_loaded_mini_match
