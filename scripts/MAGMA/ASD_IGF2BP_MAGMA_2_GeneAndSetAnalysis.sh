# UGER job submission script
# run global MAGMA analysis for each IGF2BP target list
# run gene and gene set annalysis steps for each network and each trait
# for GWAS traits: ASD ADHD BIP MDD SCZ height

#$ -cwd
#$ -N uger.magma.igf.gene
#$ -l h_vmem=4g
#$ -l h_rt=24:00:00
#$ -t 1-24
#$ -tc 24

listName=$(awk -F "\t" "NR==$SGE_TASK_ID {print \$1}" ASD_IGF2BP_MAGMA_2_GeneAndSetAnalysis.param)
trait=$(awk -F "\t" "NR==$SGE_TASK_ID {print \$2}" ASD_IGF2BP_MAGMA_2_GeneAndSetAnalysis.param) 
gwasParam=$(awk -F "\t" "NR==$SGE_TASK_ID {print \$3}" ASD_IGF2BP_MAGMA_2_GeneAndSetAnalysis.param) 
echo $listName
echo $trait

# create gene set file from each interactor list
mkdir -p temp

if [ $listName != 'COMBINED' ]
then

tail -n +2 Huang2018_IGF2BP1-3_RNAtargets.txt | \
awk -v var="$listName" 'BEGIN{print var} {if ($1==var) print $2}' | tr '\n' '\t' \
> temp/IGF2BP.${listName}.GLOBAL.${trait}.TargetSet.txt

else

tail -n +2 Huang2018_IGF2BP1-3_RNAtargets.txt | cut -f 2 | sort | uniq | \
awk -v var="$listName" 'BEGIN{print var} {print $1}' | tr '\n' '\t' \
> temp/IGF2BP.${listName}.GLOBAL.${trait}.TargetSet.txt

fi

# directory to store MAGMA output files
mkdir -p magma_output

# MAGMA directory
magmaDir=Software/MAGMA


# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/IGF2BP.${listName}.GLOBAL.genes.annot \
--pval ${gwasParam} \
--out magma_output/IGF2BP.${listName}.GLOBAL.${trait}

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/IGF2BP.${listName}.GLOBAL.${trait}.genes.raw \
--set-annot temp/IGF2BP.${listName}.GLOBAL.${trait}.TargetSet.txt \
--out magma_output/IGF2BP.${listName}.GLOBAL.${trait}

# delete temp files
rm temp/IGF2BP.${listName}.GLOBAL.${trait}.TargetSet.txt
 
