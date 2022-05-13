##########################################################################################
## UGER job submission script
## bulk RNA-seq analysis for 6 samples (d0, d21, d51 of NGN2 iNs in duplicate)
## Step 1 of 3: run trim_galore to trim reads for each of the 6 samples
##
## Author: Yu-Han Hsu, Kalliopi Tsafou
##########################################################################################

#!/bin/bash

#$ -cwd
#$ -N uger.trim
#$ -l h_vmem=4g
#$ -l h_rt=02:00:00
#$ -t 1-6
#$ -tc 6

### cores 
#$ -pe smp 4
#$ -binding linear:4
#$ -R y


# activate conda env with RNA-seq software installed
source activate genoppiEnv

# folder or file name parameters for each sample
SAMP_DIR=$(awk "NR==$SGE_TASK_ID {print \$1}" ASD_RNA_1_TrimGalore.param)
SAMP_NAME=$(awk "NR==$SGE_TASK_ID {print \$2}" ASD_RNA_1_TrimGalore.param)

# path for raw data
RNA_DATA_DIR=${SAMP_DIR}

# create output directory (if necessary)
mkdir -p 1_TrimGalore

# trim_galore command
trim_galore --illumina --paired --retain_unpaired \
--stringency 7 --phred33 --length 20 -e 0.1 \
--fastqc --output_dir 1_TrimGalore --cores 4 \
${RNA_DATA_DIR}/${SAMP_NAME}_L002_R1_001.fastq.gz \
${RNA_DATA_DIR}/${SAMP_NAME}_L002_R2_001.fastq.gz
	
# settings note:
# --phred33 --length 20 -e 0.1 are all default cutoffs
# only --stringency 7 is different (less stringent) than default

echo "TASK COMPLETED"
