##########################################################################################
## UGER job submission script
## run conditional MAGMA analysis for each ASD network (ints vs. non-ints)
## for GWAS traits: ASD ADHD BIP MDD SCZ height
##
## Author: Yu-Han Hsu
##########################################################################################

#$ -cwd
#$ -N uger.magma
#$ -l h_vmem=4g
#$ -l h_rt=48:00:00
#$ -t 1-29
#$ -tc 29


listName=$(awk "NR==$SGE_TASK_ID {print \$1}" ASD_MAGMA_Conditional.param)
echo $listName

# create gene.loc and gene set files from each interactor list
mkdir -p temp

tail -n +2 ASD_MasterInteractorTable_withInWeb.txt | \
awk -v var="$listName" '{if ($1==var) {print $3,$5,$6,$7}}' > \
temp/${listName}.gene.loc

tail -n +2 ASD_MasterInteractorTable_withInWeb.txt | \
awk -v var="$listName" 'BEGIN{print var} {if ($1==var && $4=="TRUE") print $3}' | \
tr '\n' '\t' > temp/${listName}.InteractorGeneSet.txt


# directory to store MAGMA output files
mkdir -p magma_output

# MAGMA directory
magmaDir=Software/MAGMA

# Annotation Step (SNP to gene mapping)
${magmaDir}/magma_v1.09_static/magma \
--annotate window=50 \
filter=${magmaDir}/ReferenceData/g1000_eur/g1000_eur.noMHC.bim \
--snp-loc ${magmaDir}/ReferenceData/g1000_eur/g1000_eur.bim \
--gene-loc temp/${listName}.gene.loc \
--out magma_output/${listName}


### ASD
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.genes.annot \
--pval ${magmaDir}/GwasData/ASD/iPSYCH-PGC_ASD_Nov2017.txt N=46351 \
--out magma_output/${listName}.ASD

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.ASD.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.ASD


### ADHD
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.genes.annot \
--pval ${magmaDir}/GwasData/ADHD/adhd_eur_jun2017 N=53293 \
--out magma_output/${listName}.ADHD

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.ADHD.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.ADHD


### BIP
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.genes.annot \
--pval ${magmaDir}/GwasData/BIP/daner_PGC_BIP32b_mds7a_0416a_Ntot.txt ncol=Ntot \
--out magma_output/${listName}.BIP

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.BIP.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.BIP


### MDD
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.genes.annot \
--pval ${magmaDir}/GwasData/MDD/PGC_UKB_depression_genome-wide.txt snp-id=MarkerName pval=P N=500199 \
--out magma_output/${listName}.MDD

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.MDD.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.MDD


### SCZ
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.genes.annot \
--pval ${magmaDir}/GwasData/SCZ/daner_natgen_pgc_eur_Ntot.txt ncol=Ntot \
--out magma_output/${listName}.SCZ

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.SCZ.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.SCZ


### height
# Gene Analysis Step (calculate gene p-values + other gene-level metrics)
${magmaDir}/magma_v1.09_static/magma \
--bfile ${magmaDir}/ReferenceData/g1000_eur/g1000_eur \
--gene-annot magma_output/${listName}.genes.annot \
--pval ${magmaDir}/GwasData/Height/Meta-analysis_Wood_et_al+UKBiobank_2018.txt ncol=N \
--out magma_output/${listName}.height

# Gene Set Analysis Step (calculate enrichment of interactors compared to non-interactors)
${magmaDir}/magma_v1.09_static/magma \
--gene-results magma_output/${listName}.height.genes.raw \
--set-annot temp/${listName}.InteractorGeneSet.txt \
--out magma_output/${listName}.height


# delete temp files
rm temp/${listName}.gene.loc temp/${listName}.InteractorGeneSet.txt
