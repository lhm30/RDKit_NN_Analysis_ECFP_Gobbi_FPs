### K-Near Neighbour ECFP or Gobbi analysis of SMILES or SDF files
### Lewis Mervin
### lhm30@cam.ac.uk

## Script to perform near neighbour analysis of compound sets
## Can be used to calculate intra- or inter- similarity of compound sets
## Intra-similarity is automatically performed if no second file provided
## Can take multiple files using * wildcards

#libraries
import rdkit
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import *
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
import cPickle
import os
import sys
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing
from functools import partial
import operator
from optparse import OptionParser
import os
import glob
multiprocessing.freeze_support()

#optionparser options
parser = OptionParser()
parser.add_option("--i1", dest="inf1", help="First inputfile", metavar="FILE")
parser.add_option("--i2", dest="inf2", help="Second inputfile (optional)", metavar="FILE")
parser.add_option("-o", "--out", dest="outfile", help="Outputfile (optional)", metavar="FILE")
parser.add_option("--d1", "--delim1", default=' ', type=str, dest="delim1", help="Delimiter for first file")
parser.add_option("--d2", "--delim2", default=' ', type=str, dest="delim2", help="Delimiter for second file")
parser.add_option("--delimcol1", default=0, type=int, dest="delimcol1", help="Delimiter column of smiles in first file")
parser.add_option("--delimcol2", default=0, type=int, dest="delimcol2", help="Delimiter column of smiles in second file")
parser.add_option("-b", "--bits", default=2048, type=int, dest="bits", help="No. bits in morgan fingerprint")
parser.add_option("-r", "--radius", default=2, type=int, dest="radius", help="Radius of morgan fingerprint")
parser.add_option("-t", "--topn", default=1, type=int, dest="topn", help="Top n neighbors")
parser.add_option("-n", "--ncores", default=1, type=int, dest="ncores", help="No. cores")
parser.add_option("--pf1", action="store_true", default=False, dest="pf1", help="Parallel process list of first files, not lines per file")
parser.add_option("--pf2", action="store_true", default=False, dest="pf2", help="Parallel process list of second files, not lines per file")
parser.add_option("--gobbi", action="store_true", default=False, dest="gobbifp", help="Use Gobbi 2D pharmacophore FPs")
parser.add_option("--mw", default=None, type=int, dest="mw", help="Max. Mw filter")
(options, args) = parser.parse_args()
#if tab delimited then set to \t
if options.delim1 == 'tab': options.delim1 = '\t'
if options.delim2 == 'tab': options.delim2 = '\t'

N_cores = options.ncores

#generate rdkit ECFP/GOBBI fingerprint per smile
def molfp(smile):
	m = Chem.MolFromSmiles(smile)
	if m is None: raise ValueError('None mol in function "molfp" (line 58)')
	if options.mw is not None and Descriptors.MolWt(m) > options.mw: raise ValueError('Mol too small')
	if not options.gobbifp: return AllChem.GetMorganFingerprintAsBitVect(m,options.radius, nBits=options.bits)
	else: return Generate.Gen2DFingerprint(m,Gobbi_Pharm2D.factory)

#generate rdkit ECFP/GOBBI fingerprints for list of smiles
def genmol(ms, delim=' ', delimcol=0):
	ret = []
	for m in ms:
		m = m.split(delim)[delimcol]
		try:
			ret.append(molfp(m))
		except: pass
	return ret

#generate rdkit ECFP fingerprints for smile file
def genmol_wopen(inf, delim=' ', delimcol=0):
	ms = open(inf).read().splitlines()
	ret = []
	for m in ms:
		m = m.split(delim)[options.delimcol]
		try:
			ret.append(molfp(m))
		except: pass
	return ret

#generate rdkit ECFP fingerprints for one sdf file
def genmol_sdf(ms):
	suppl = Chem.SDMolSupplier(ms)
	ret = []
	for m in suppl:
		if m is None: continue
		try:
			if not options.gobbifp: ret.append(AllChem.GetMorganFingerprintAsBitVect(m,options.radius, nBits=options.bits))
			else: ret.append(Generate.Gen2DFingerprint(m,Gobbi_Pharm2D.factory))
		except: pass
	return ret

#parallel process list of smiles files (file per core)
def multiple_smi_parallel_file(infiles, d, c):
	pooler = Pool(processes=N_cores)  # set up resources
	genmol_wopen_worker=partial(genmol_wopen, delim=d, col=c)
	jobs = pooler.imap(genmol_wopen_worker, infiles)
	processed_smi = []
	for i, result in enumerate(jobs):
		processed_smi += result
	pooler.close()
	pooler.join()
	return processed_smi

#parallel process list of smiles files (subsets of each file per core)
def multiple_smi(infiles, d, c):
	processed_smi = []
	for f in infiles:
		query = open(f).read().splitlines()
		smiles_per_core = int(math.ceil(len(query) / N_cores)+1)/2
		print 'Smiles per core: ' + str(smiles_per_core)
		chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
		pooler = Pool(processes=N_cores)  # set up resources
		genmol_worker=partial(genmol, delim=d, col=c)
		jobs = pooler.imap(genmol_worker, chunked_smiles)
		for i, result in enumerate(jobs):
			processed_smi += result
		pooler.close()
		pooler.join()
	return processed_smi	

#parallel process one smile file (subsets of file per core)
def single_smi(inf, d, c):
	if d == 'tab': d = '\t'
	print d
	query = open(inf).read().splitlines()
	smiles_per_core = int(math.ceil(len(query) / N_cores)+1)/2
	print 'Smiles per core: ' + str(smiles_per_core)
	chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
	pooler = Pool(processes=N_cores)  # set up resources
	genmol_worker=partial(genmol, delim=d, delimcol=c)
	jobs = pooler.imap(genmol_worker, chunked_smiles)
	processed_smi = []
	for i, result in enumerate(jobs):
		processed_smi += result
	pooler.close()
	pooler.join()
	return processed_smi
	
#process one sdf file
def single_sdf(ms):
	processed_smi = genmol_sdf(ms)
	return processed_smi

#parallel process sdf files (file per core)
def multiple_sdf_parallel_file(infiles):
	pooler = Pool(processes=N_cores)  # set up resources
	jobs = pooler.imap(genmol_sdf, infiles)
	processed_smi = []
	for i, result in enumerate(jobs):
		processed_smi += result
	pooler.close()
	pooler.join()
	return processed_smi

#check smiles or sdf (if not quit)
def check_extension(inf):
	extension = inf.split('.')[-1]
	if extension not in ['smi','smiles','sdf']:
		print 'Provide "smi" or "sdf" filetype'
		quit()
	return extension

#check extension, number of files for processing, and method of parallel processing
#(if applicable) and perform suitable processing
def process_inf(inf,pf,delim,col):
	extension = check_extension(inf)
	if extension in ['smi','smiles']:
		#check if multiple files (using * wildcard)
		if '*' in inf:
			infiles = glob.glob(inf)
			#check if parallel file
			if pf: processed_smi = multiple_smi_parallel_file(infiles,delim,col)
			else: processed_smi = multiple_smi(infiles,delim,col)
		else: processed_smi = single_smi(inf,delim,col)
	#else input is sdf file
	else:
		if '*' in inf:
			infiles = glob.glob(inf)
			processed_smi = multiple_sdf_parallel_file(infiles)
		else: processed_smi = single_sdf(inf)
	return processed_smi

#near neighbour similarity worker
def do_sim(i1,i2,intra=False):
	ret = []
	for i in i1:
		if intra: sims = DataStructs.BulkTanimotoSimilarity(i[1],i2[:i[0]] + i2[i[0]+1:])
		else: sims = DataStructs.BulkTanimotoSimilarity(i[1],i2)
		if options.topn ==1: ret.append(max(sims))
		else: ret.append(np.average(sorted(sims,reverse=True)[:options.topn]))
	return ret

#perform either intra- or inter- near neighbour analysis 
def nn_sim(inp1,inp2,intra=False):
	sims_per_core = int(math.ceil(len(inp1) / N_cores)+1)/2
	print 'Sim calculations per core: ' + str(sims_per_core)
	chunked_sims = [zip(range(x,x+sims_per_core,1),inp1[x:x+sims_per_core]) for x in xrange(0, len(inp1), sims_per_core)]
	pooler = Pool(processes=N_cores)  # set up resources
	sim_worker=partial(do_sim, i2=inp2, intra=intra)
	jobs = pooler.imap(sim_worker, chunked_sims)
	processed_sim = []
	for idx, result in enumerate(jobs):
		percent = (float(idx)/float(len(inp1)))*100 + 1
		sys.stdout.write('Calculating sims: %3d%%\r' % percent)
		sys.stdout.flush()
		processed_sim += result
	pooler.close()
	pooler.join()
	print
	return processed_sim

#write out specified name or create suitable naming
def write_out(results):
	if options.gobbifp: fp = 'gobbi'
	else: fp = 'ecfp'
	if options.outfile: of = open(options.outfile,'w')
	else:
		if options.inf2:
			sf = options.inf2.split('/')[-1]
			of = open(options.inf1 + '_vs_' + sf + '_' + fp + 'sims.txt','w')
		else: of = open(options.inf1 + '_' + fp + '_sims.txt','w')
	for line in results: of.write(str(round(line,3)) + '\n')
	of.close()
	return

#calc inter- similarity between two input files
def process_sims_two_files(inp1,inp2):
	print 'Calculating inter-similarities between two files'
	results = nn_sim(inp1,inp2)
	write_out(results)
	return

#calc intra- similarity between two input files
def process_sims(inp1):
	print 'Calculating intra-similarities within file'
	results = nn_sim(inp1,inp1,intra=True)
	write_out(results)
	return

if __name__ == "__main__":
	#process first file
	processed_smi1 = process_inf(options.inf1, options.pf1, options.delim1, options.delimcol1)
	print 'Processed first file...'
	#if second file then process and do inter-sim
	if options.inf2:
		print 'Processing second file...',
		processed_smi2 = process_inf(options.inf2, options.pf2, options.delim2, options.delimcol2)
		print 'processed'
		process_sims_two_files(processed_smi1,processed_smi2)
	#if not second, do intra-sim
	else: process_sims(processed_smi1)

