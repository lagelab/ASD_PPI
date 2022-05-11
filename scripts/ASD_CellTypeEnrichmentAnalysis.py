##########################################################################################
## Cell type enrichment analysis of ASD PPI networks using Velmeshev2019 data
## (1) assess overlap between networks vs. genes expressed in >50% of cells in cell type
## (2) assess overlap between networks vs. ASD-DEGs in each cell type
##
## Author: Yu-Han Hsu
##########################################################################################

import gzip
import pandas as pd
from scipy.stats import hypergeom


# gene symbol <-> ENSG ID mapping dicts
ensgDict = {} # gene symbol as key, set of ENSG IDs as value
symDict = {} # ENSG ID as key, gene symbol as value
with open('../data/Velmeshev2019_genes.tsv','r') as inFile:
	for line in inFile:
		fields = line.strip().split()
		ensgDict.setdefault(fields[1],set()).add(fields[0])
		symDict[fields[0]] = fields[1]


# dict containing interactor and non-interactor genes in each network
intDict = {'ADNP':[set(),set()],'ANK2':[set(),set()],'ARID1B':[set(),set()],
	'CHD8':[set(),set()],'CTNNB1':[set(),set()],'DYRK1A':[set(),set()],
	'GIGYF1':[set(),set()],'MED13L':[set(),set()],'PTEN':[set(),set()],
	'SCN2A':[set(),set()],'SHANK3':[set(),set()],'SYNGAP1':[set(),set()],
	'TLK2':[set(),set()],'COMBINED':[set(),set()]} # store [ints,non-ints] gene symbols
	
with open('../data/ASD_MasterInteractorTable.txt','r') as inFile:
	next(inFile)
	for line in inFile:
		fields = line.split()
		if fields[0] in intDict:
			if fields[3]=='TRUE': intDict[fields[0]][0].add(fields[2])
			else: intDict[fields[0]][1].add(fields[2])

# check int and non-int counts
#for d in intDict:
#	print(d)
#	print(len(intDict[d][0]))
#	print(len(intDict[d][1]))


# nested dict: cell type -> gene symbol -> % expressed 
prcDict = {}
allGenes = set()
with gzip.open('../data/Velmeshev2019_PercExpMatrix.txt.gz','rt') as inFile:
	headers = next(inFile).strip().split() # gene (ENSG ID) followed by cell types

	for ct in headers[1:]:
		prcDict[ct] = {}

	for line in inFile:
		fields = line.strip().split()
		gene  = symDict[fields[0]]
		allGenes.add(gene)
		
		for i in range(1,len(fields)):
			ct = headers[i]
			if gene not in prcDict[ct] or float(fields[i]) > prcDict[ct][gene]:
				prcDict[ct][gene] = float(fields[i])
		
print(len(prcDict))
print(len(allGenes))


# nested dicts: DEG type -> cell type -> DEGs
degDict = {}
degDict['All'] = {}
degDict['Up'] = {}
degDict['Down'] = {}

degTable = pd.read_excel('../data/Velmeshev2019_DataS4.xls')

for i in range(degTable.shape[0]): # iterate through each row
	
	degDict['All'].setdefault(degTable['Cell type'][i],set()).add(degTable['Gene name'][i])

	if float(degTable['Fold change'][i]) > 0:
		degDict['Up'].setdefault(degTable['Cell type'][i],set()).add(degTable['Gene name'][i])
	else:
		degDict['Down'].setdefault(degTable['Cell type'][i],set()).add(degTable['Gene name'][i])


# ----------------------------------------------------------------------------------------
# (1) overlap enrichment of ints vs. genes expressed in >50% of cells in each cell type

with open('../output/ASD_scExpEnrichment.txt','w') as outFile:

	outHeaders = ['Test','Dataset','CellType','PopCount',
		'IntCount','ExpCount','OverlapCount','P']
	
	outFile.write('\t'.join(outHeaders)+'\n')
	
	# iterate through network
	for d in intDict:
		ints = intDict[d][0]
		nonInts = intDict[d][1]
		
		# iterate through cell type
		for ct in prcDict:
		
			# global test
			# population = genes found in scRNA dataset
			popGenes = allGenes
			M = len(popGenes) # population count
			N = len(ints & popGenes) # int count
			ctExpressed = set([g for g in popGenes if prcDict[ct][g] > 50])
			n = len(ctExpressed) # expressed count
			overlapGenes = list(ctExpressed & ints)
			x = len(overlapGenes) # overlap count
			p = hypergeom.sf(x-1, M, n, N)
			
			outInfo = ['Global',d,ct,M,N,n,x,p]
			outFile.write('\t'.join(map(str,outInfo)) + '\n')
			
			# conditional test
			# population = (ints|non-ints) & genes found in scRNA dataset
			popGenes = (ints | nonInts) & allGenes
			M = len(popGenes) # population count
			N = len(ints & popGenes) # int count
			ctExpressed = set([g for g in popGenes if prcDict[ct][g] > 50])
			n = len(ctExpressed) # expressed count
			overlapGenes = list(ctExpressed & ints)
			x = len(overlapGenes) # overlap count
			p = hypergeom.sf(x-1, M, n, N)
			
			outInfo = ['Conditional',d,ct,M,N,n,x,p]
			outFile.write('\t'.join(map(str,outInfo)) + '\n')


# ----------------------------------------------------------------------------------------
# (2) overlap enrichment of ints vs. ASD-DEGs in each cell type

with open('../output/ASD_scDegEnrichment.txt','w') as outFile:

	outHeaders = ['Test','Dataset','DegType','CellType','PopCount',
		'IntCount','DegCount','OverlapCount','P','OverlapGenes']
	
	outFile.write('\t'.join(outHeaders)+'\n')

	# iterate through network
	for d in intDict:
		ints = intDict[d][0]
		nonInts = intDict[d][1]

		# iterate through DEG type (All, Up, Down)
		for degType in degDict:
		
			# iterate through cell type
			for ct in degDict[degType]:
				
				# global test
				# population = genes found in scRNA dataset
				popGenes = allGenes
				M = len(popGenes) # population count
				N = len(ints & popGenes) # int count
				ctDegs = degDict[degType][ct] & popGenes
				n = len(ctDegs) # DEG count
				overlapGenes = list(ctDegs & ints)
				x = len(overlapGenes) # overlap count
				p = hypergeom.sf(x-1, M, n, N)

				outInfo = ['Global',d,degType,ct,M,N,n,x,p,','.join(overlapGenes)]
				outFile.write('\t'.join(map(str,outInfo)) + '\n')

				# conditional test
				# population = (ints|non-ints) & genes found in scRNA dataset
				popGenes = (ints | nonInts) & allGenes
				M = len(popGenes) # population count
				N = len(ints & popGenes) # int count
				ctDegs = degDict[degType][ct] & popGenes
				n = len(ctDegs) # DEG count
				overlapGenes = list(ctDegs & ints)
				x = len(overlapGenes) # overlap count
				p = hypergeom.sf(x-1, M, n, N)

				outInfo = ['Conditional',d,degType,ct,M,N,n,x,p,','.join(overlapGenes)]
				outFile.write('\t'.join(map(str,outInfo)) + '\n')

