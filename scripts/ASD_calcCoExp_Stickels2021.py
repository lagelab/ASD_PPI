##########################################################################################
## Create gene x gene co-expression matrix containing -log10(Fisher's exact p-value)
## using mouse neocortex gene expression matrix from Stickels 2021 (Slide-seqV2)
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
# map mouse genes to human genes

# Ensembl protein-coding human genes
protGenes = set()
with open('../data/Homo_sapiens.GRCh38.84.GeneAnnotations.txt','r') as inFile:
	next(inFile)
	
	for line in inFile:
		fields = line.strip().split('\t')
		if fields[6]=='protein_coding':
			protGenes.add(fields[1])

# dict with  # mouse gene symbol as key, 'DB Class Key' as value
mouseToHomDict = {}

# dict with 'DB Class Key' as key, list of human protein-coding gene symbols as value
homToHumanDict = {} 

# MGI mouse-human gene homology mapping table
with open('../data/HOM_MouseHumanSequence.rpt.txt','r') as inFile:
	next(inFile)
	
	for line in inFile:
		fields = line.strip().split('\t')
	
		if fields[1] == 'mouse, laboratory':
			mouseToHomDict[fields[3]] = fields[0]
		elif fields[1] == 'human' and fields[3] in protGenes:
			homToHumanDict.setdefault(fields[0],[]).append(fields[3])

print('# of mouse genes mapped to homlogy IDs: ' + str(len(mouseToHomDict))) # 20597
print('# of homology IDs mapped to human pc genes: ' + str(len(homToHumanDict))) # 18008


# ----------------------------------------------------------------------------------------
# read in mouse neocoretx gene x barcode expression matrix (Stickels 2021)

# dict with gene name as key, set of barcodes expressing the gene as value
expDict = {}
barcodes = set()

# NOTE: input file below not included in GitHub
with gzip.open('../data/Puck_190921_19.digital_expression.txt.gz','rt') as inFile:
	headers = next(inFile).strip().split('\t')

	for line in inFile:
		fields = line.strip().split('\t')
		
		mouseGene = fields[0]
		if (mouseGene in mouseToHomDict) and (mouseToHomDict[mouseGene] in homToHumanDict):
			humanGenes = homToHumanDict[mouseToHomDict[mouseGene]]
		
			# only save genes that are not missing or detected across all barcodes:
			counts = list(map(int,fields[1:]))
			nDetected = sum([x > 0 for x in counts])
			
			if nDetected > 0 and nDetected < len(counts):
				for g in humanGenes: # save each unique gene symbol as separate entry
					if g not in expDict: # skip duplicated gene names
						for i in range(len(counts)):
							if counts[i] > 0: # only save if non-zero count
								expDict.setdefault(g,set()).add(headers[i+1])
								barcodes.add(headers[i+1])

print('# of human genes in raw exp matrix: ' + str(len(expDict))) # 15152 genes
print('# of barcodes in raw exp matrix: ' + str(len(barcodes))) # 33611 barcodes


# ----------------------------------------------------------------------------------------
# generate gene x gene matrix storing -log10(Fisher's exact p-value) for each gene pair
# 2x2 contingency table: # barcodes detecting both genes, 1st gene, 2nd gene, or neither

pDict = {} # nested dict, 1st key = gene1, 2nd key = gene2, value = p-value
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
		# equivalent to: hypergeom.sf(x-1, M, n, N)

		pDict.setdefault(g1,{})[g2] = -numpy.log10(p) # -log10 p range between 0 and inf

print('# of human genes in co-exp matrix: ' + str(len(pDict)))


# output co-exp matrix
with open('../output/Stickels2021_CoExpFisherNegLogP.' + str(iStart) + '-' + str(iEnd) + '.txt','w') as pFile:
	pFile.write('Gene\t' + '\t'.join(geneList) + '\n')
	
	for g1 in geneList[iStart:iEnd]:
		pFields = [g1] + ['%.2f' % pDict[g1][g2] if pDict[g1][g2] > 0 else str(0) for g2 in geneList]
		pFile.write('\t'.join(pFields) + '\n')
		
		
# time script
elapsed = time.time() - start
print('SCRIPT RUN TIME: ' + str(elapsed))

