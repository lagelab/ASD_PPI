#!/bin/bash

### UGER submission script to run hisat2 alignment for each of the 6 samples

#$ -cwd
#$ -N uger.hisat2

#$ -l h_vmem=8g
#$ -l h_rt=02:00:00

#$ -t 1-6
#$ -tc 6

### cores
#$ -pe smp 8
#$ -binding linear:8
#$ -R y

# activate conda env with RNA-seq software installed
source activate genoppiEnv

# sample names
SAMP_LONG=$(awk "NR==$SGE_TASK_ID {print \$1}" ASD_RNA_2_Hisat2.param)
SAMP_SHORT=$(awk "NR==$SGE_TASK_ID {print \$2}" ASD_RNA_2_Hisat2.param)
SAMP_TIME=$(awk "NR==$SGE_TASK_ID {print \$3}" ASD_RNA_2_Hisat2.param)

echo $SAMP_LONG
echo $SAMP_SHORT
echo $SAMP_TIME

# create output directory (if necessary)
mkdir -p 2_Hisat2

# hisat2 alignment command
hisat2 -p 8 --rg-id=${SAMP_SHORT} --rg SM:${SAMP_TIME} --rg PL:ILLUMINA --dta --rna-strandness FR -x refs/grch38_tran/genome_tran -1 1_TrimGalore/${SAMP_LONG}_L002_R1_001_val_1.fq.gz -2 1_TrimGalore/${SAMP_LONG}_L002_R2_001_val_2.fq.gz -S 2_Hisat2/${SAMP_SHORT}.sam

# sort and convert to bam
samtools sort -@ 8 -o 2_Hisat2/${SAMP_SHORT}.bam 2_Hisat2/${SAMP_SHORT}.sam

# bam index file
samtools index 2_Hisat2/${SAMP_SHORT}.bam

# delete sam file
if [ -f 2_Hisat2/${SAMP_SHORT}.bam ]; then
	rm 2_Hisat2/${SAMP_SHORT}.sam
fi

echo "TASK COMPLETED"
