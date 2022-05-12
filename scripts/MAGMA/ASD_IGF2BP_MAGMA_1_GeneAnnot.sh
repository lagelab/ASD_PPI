# UGER job submission script
# run global MAGMA analysis for each IGF2BP target list
# run gene annotation step for each list (same across all traits)

#$ -cwd
#$ -N uger.magma.igf.annot
#$ -l h_vmem=4g
#$ -l h_rt=08:00:00
#$ -t 1-4


listName=$(awk "NR==$SGE_TASK_ID {print \$1}" ASD_IGF2BP_MAGMA_1_GeneAnnot.param)
echo $listName

# create gene.loc file for each network
mkdir -p temp

awk "NR==$SGE_TASK_ID {print \$2}" ASD_IGF2BP_MAGMA_1_GeneAnnot.param | tr ',' '\n' > \
temp/${listName}.baits.txt

cat temp/${listName}.baits.txt

grep -vwf temp/${listName}.baits.txt Ensembl_BioMart_GRCh37.gene.loc > \
temp/IGF2BP.${listName}.GLOBAL.gene.loc


# MAGMA directory
magmaDir=Software/MAGMA

# directory to store MAGMA output files
mkdir -p magma_output

# Annotation Step (SNP to gene mapping)
${magmaDir}/magma_v1.09_static/magma \
--annotate window=50 \
filter=${magmaDir}/ReferenceData/g1000_eur/g1000_eur.noMHC.bim \
--snp-loc ${magmaDir}/ReferenceData/g1000_eur/g1000_eur.bim \
--gene-loc temp/IGF2BP.${listName}.GLOBAL.gene.loc \
--out magma_output/IGF2BP.${listName}.GLOBAL


# delete temp files
rm temp/${listName}.baits.txt temp/IGF2BP.${listName}.GLOBAL.gene.loc 
