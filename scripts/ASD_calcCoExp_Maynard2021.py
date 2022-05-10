##########################################################################################
## Create gene x gene co-expression matrix containing -log10(Fisher's exact p-value)
## using human dorsolateral prefrontal cortex matrix from Maynard 2021 (Visium)
##
## Author: Yu-Han Hsu
##########################################################################################

import sys, gzip, time, numpy
from scipy.stats import fisher_exact

# time script
start  = time.time()

# input parameters (for parallelization)
# generate co-exp stats for genes with indices form iStart to (iEnd-1)
iStart = int(sys.argv[1])
iEnd = int(sys.argv[2])


# ----------------------------------------------------------------------------------------
# read in row index -> gene name mapping info
rowDict = {} # row index as key, protein-coding gene name as value
rowGenes = set()

# NOTE: input file below not included in GitHub
with open('../data/Maynard2021_rowData.txt','r') as inFile:
	next(inFile)

	rowInd = 1
	for line in inFile:
		fields = line.strip().split('\t')

		# skip duplciated (21) and non-protein-coding gene names
		if fields[5] not in rowGenes and fields[7]=='protein_coding':
			rowDict[str(rowInd)] = fields[5]
			rowGenes.add(fields[5])

		rowInd += 1

print('# of row indices mapped to pc gene names: ' + str(len(rowDict))) # 19893


# ----------------------------------------------------------------------------------------
# read in human perfrontal coretex gene x barcode expression data (Maynard 2021)

# dict with gene name as key, set of barcodes expressing the gene as value
expDict = {}
barcodes = set()

# NOTE: input file below not included in GitHub
with gzip.open('../data/Maynard2021_counts.mtx.gz','rt') as inFile:
	next(inFile) # skip 2 header lines
	next(inFile)

	for line in inFile:
		fields = line.strip().split()
		
		# skip non-protein-coding gene names
		if fields[0] in rowDict:
			expDict.setdefault(rowDict[fields[0]],set()).add(fields[1])
			barcodes.add(fields[1])

print('# of genes in raw exp matrix: ' + str(len(expDict))) # 17934 genes
print('# of barcodes in raw exp matrix: ' + str(len(barcodes))) # 47681 barcodes

# checked that no genes are expressed in all barcodes
#badCount = 0
#for g in expDict:
#	if len(expDict[g]) == (len(barcodes)):
#		badCount += 1
#print(badCount)


# ----------------------------------------------------------------------------------------
# generate gene x gene matrix storing -log10(Fisher's exact p-value) for each gene pair
# 2x2 contingency table: # barcodes detecting both genes, 1st gene, 2nd gene, or neither

pDict = {} # nested dict, 1st key = gene 1, 2nd key = gene 2, value = p-value
geneList = list(expDict.keys())
geneList.sort() # ensure genes are always in same order

for g1 in geneList[iStart:iEnd]:
	for g2 in geneList:
	
		M = len(barcodes) # total # of barcodes
		N = len(expDict[g1]) # # of barcodes expressing g1
		n = len(expDict[g2]) # # of barcodes expressing g2
		x = len(expDict[g1] & expDict[g2]) # # of barcodes expressing both 
		table = numpy.array([[x,n-x],[N-x,M-n-N+x]]) # 2 x 2 contingency table
		oddsRatio,p = fisher_exact(table,alternative='greater')
		#equivalent to: p = hypergeom.sf(x-1, M, n, N)

		pDict.setdefault(g1,{})[g2] = -numpy.log10(p)  # -log10 p range between 0 and inf

print('# of genes in co-exp matrix: ' + str(len(pDict)))


# output co-exp matrix
with open('../output/Maynard2021_CoExpFisherNegLogP.' + str(iStart) + '-' + str(iEnd) + '.txt','w') as outFile:
	outFile.write('Gene\t' + '\t'.join(geneList) + '\n')
	
	for g1 in geneList[iStart:iEnd]:
		outFields = [g1] + ['%.2f' % pDict[g1][g2] if pDict[g1][g2] > 0 else str(0) for g2 in geneList]
		outFile.write('\t'.join(outFields) + '\n')

		
# time script
elapsed = time.time() - start
print('SCRIPT RUN TIME: ' + str(elapsed))

