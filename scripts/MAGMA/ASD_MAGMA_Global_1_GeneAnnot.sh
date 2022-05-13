##########################################################################################
## UGER job submission script
## run global MAGMA analysis for each ASD network (ints vs. genomic background)
## Step 1 of 2: run gene annotation step for each network (same across all traits)
##
## Author: Yu-Han Hsu
##########################################################################################

#$ -cwd
#$ -N uger.magma.global.annot
#$ -l h_vmem=4g
#$ -l h_rt=08:00:00
#$ -t 1-29


listName=$(awk "NR==$SGE_TASK_ID {print \$1}" ASD_MAGMA_Conditional.param)
echo $listName

# create gene.loc file for each network
mkdir -p temp

awk "NR==$SGE_TASK_ID {print \$2}" ASD_MAGMA_Conditional.param | tr ',' '\n' > \
temp/${listName}.baits.txt

cat temp/${listName}.baits.txt

grep -vwf temp/${listName}.baits.txt Ensembl_BioMart_GRCh37.gene.loc > \
temp/${listName}.GLOBAL.gene.loc


# MAGMA directory
magmaDir=Software/MAGMA

# directory to store MAGMA output files
mkdir -p magma_output

# Annotation Step (SNP to gene mapping)
${magmaDir}/magma_v1.09_static/magma \
--annotate window=50 \
filter=${magmaDir}/ReferenceData/g1000_eur/g1000_eur.noMHC.bim \
--snp-loc ${magmaDir}/ReferenceData/g1000_eur/g1000_eur.bim \
--gene-loc temp/${listName}.GLOBAL.gene.loc \
--out magma_output/${listName}.GLOBAL


# delete temp files
rm temp/${listName}.baits.txt temp/${listName}.GLOBAL.gene.loc
