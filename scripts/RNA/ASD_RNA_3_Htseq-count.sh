##########################################################################################
## UGER job submission script
## bulk RNA-seq analysis for 6 samples (d0, d21, d51 of NGN2 iNs in duplicate)
## Step 3 of 3: run htseq-count to generate gene counts for each of the 6 samples
##
## Author: Yu-Han Hsu, Kalliopi Tsafou
##########################################################################################

#!/bin/bash

#$ -cwd
#$ -N uger.htseq-count
#$ -l h_vmem=8g
#$ -l h_rt=04:00:00
#$ -t 1-6
#$ -tc 6


# activate conda env with RNA-seq software installed
source activate genoppiEnv

# sample names
SAMP_SHORT=$(awk "NR==$SGE_TASK_ID {print \$2}" ASD_RNA_2_Hisat2.param)

echo $SAMP_SHORT

# create output directory (if necessary)
mkdir -p 3_Htseq-count 

# htseq-count command
htseq-count --format bam --order pos --mode union --stranded yes \
--minaqual 1 --type exon --idattr gene_id \
2_Hisat2/${SAMP_SHORT}.bam refs/Homo_sapiens.GRCh38.84.gtf > 3_Htseq-count/${SAMP_SHORT}.geneCounts.tsv 

echo "TASK COMPLETED"
