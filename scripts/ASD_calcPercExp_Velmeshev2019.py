##########################################################################################
## Create gene x cell type matrix containing % of cells in each cell type expressing gene
## using human cortex gene expression matrix from Velmeshev 2019 (scRNA-seq)
##
## Author: Yu-Han Hsu
##########################################################################################

import gzip
import pandas as pd

# read in cell to cell type mapping info
celltypeDict = {} # cell ID -> cell type dict
cellcountDict = {} # cell type -> # of cells in cell type

# NOTE: input file below not included in GitHub 
with open('../data/Velmeshev2019_meta.tsv','r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.split()
		celltypeDict[fields[0]] = fields[1]

		if fields[1] in cellcountDict: cellcountDict[fields[1]] += 1
		else: cellcountDict[fields[1]] = 0

print(len(celltypeDict))
print(len(cellcountDict))


print('### Processing Expression Matrix ###')
# nested dicts: 1st key = gene, 2nd key = cell type, value = # cells detecting genes
expDict = {}
lineCount = 0

# NOTE: input file below not included in GitHub
with gzip.open('../data/Velmeshev2019_exprMatrix.tsv.gz','rt') as inFile:
	headers = next(inFile).split() # gene (ENSG ID) followed by cell IDs
	
	for line in inFile:
		fields = line.strip().split()
		gene = fields[0]

		expDict[gene] = {}
		for ct in cellcountDict:
			expDict[gene][ct] = 0

		for i in range(1,len(fields)):
			if float(fields[i]) > 0: #non-zero expression
				cell = headers[i]
				ct = celltypeDict[cell]
				expDict[gene][ct] += 1
		
		lineCount = lineCount + 1
		if lineCount % 1000 == 0:
			print(lineCount)

print('### Done Processing Expression Matrix ###')
print(len(expDict))


# output gene x cell type matrix containing % of cells in cell type expressing gene
with gzip.open('../data/Velmeshev2019_PercExpMatrix.txt.gz','wt') as outFile:
	celltypes = cellcountDict.keys()
	outFile.write('ENSG_ID\t' + '\t'.join(celltypes) + '\n')

	for gene in expDict:
		percs = []
		for ct in celltypes:
			percs.append(float(expDict[gene][ct])/cellcountDict[ct]*100)

		outFile.write(gene + '\t' +  '\t'.join(['%.2f' % x for x in percs]) + '\n')

