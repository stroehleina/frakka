import sys
from datetime import datetime
import os


# Log a message to stderr
def msg(*args, **kwargs):
	now = f'[{datetime.now()}]: '
    print(now, *args, file=sys.stderr, **kwargs)


# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
    msg(*args, **kwargs)
    sys.exit(1)


def checkPath(path):
	if not os.path.exists(path):
		err(f'File {path} does not exist! Exiting.')


def readTabSep(file, type):
	'''
	checks if tab-separated file exists and reads it into a list of lists
	which it returns
	'''

	checkPath(file)

	with open(file, 'r') as f:
		msg(f'Reading {type} file {file}')
		l_of_l = []
		first = True
		col = None
		for line in f:
			cols = line.rstrip('\n').split('\t')
			if first:
				first = False
				col = len(cols) # Return number of columns present in tab-separated file
			l_of_l.append(cols)
		return l_of_l, col 


def readKraken(file):
	'''
	Reads Kraken standard output file (STDOUT / --output) 
	and returns a dict of all classified reads with read id as key
	and a dict as value (with keys taxid:, species: and kmerstr:)	
	'''

	krak = {}

	f = readTabSep(file, 'Kraken output')[0]
	for line in f:
			uc, read, taxid, kmerstr = line[0], line[1], line[2], line[4]
			conf = getConfidence(taxid, kmerstr)
			if uc != 'U':
				krak[read]={'taxid' : taxid, 'species' : '', 'conf' : conf}
    return krak


def readKReport(file):
	'''
	Reads in a Kraken report file created via the --report option
	and returns a dict of taxids as keys matched to a dict of name: and count:
	'''

	report = {}

	f,form  = readTabSep(file, 'Kraken report')
	level = form - 3
	tax_col = form - 2
	spec_col = form - 1
	for line in f:
		if line[level] == 'S':
		count, taxid, name = line[2], line[tax_col], line[spec_col]
		report[taxid] = {'count' : count, 'name' : name}

	return report


def fileOfFiles(file):
	'''
	Reads a tab-separated file of report and output-file pairs specified by --fof / -f (and optionally, a third column with species identifiers)
	and returns a list of lists to iterate over
	'''
	f = readTabSep(file, 'tab-separated input')[0]

	return f


def getConfidence(taxid, kmerstr):
	'''
	Reads the kmerstring and called taxid from column 5 and 3 of the Kraken standard output file (STDOUT / --output) 
	and returns the confidence score for the read
	see: https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring

	'''

	# Format: 562:13 561:4 A:31 0:1 562:3
	# paired read data will contain a '|:|' token in this list to indicate the end of one read and the beginning of another
	r_lst = kmerstr.split(' |:| ')
	conf_lst = []
	for r in r_lst:
		allkmers = 0
		taxkmers = 0
		cons_kmers = r.split(' ')
			for k in cons_kmers:
				run = k.split(':')
				if run[0] != 'A':
					allkmers+=run[1]
					if run[0] == taxid:
						taxkmers+=run[1]				
		conf_lst.append(taxkmers/allkmers)
		
	conf = sum(conf_lst)/len(conf_lst)

	return conf


def getCounts(krak_f, report_f):
	'''
	Implements the --counts option to return counts for each species
	instead of confidence scores for individual reads
	'''

	species = {}

	krak = readKraken(krak_f)
	report = readKReport(report_f)

	for read in krak:


	# krak[read]={'taxid' : taxid, 'species' : '', 'conf' : conf}
	# report[taxid] = {'count' : count, 'name' : name}


def output(record, name, header=True, sep='\t', counts=False, taxid=False):
	msg(f"Writing output for file {name}")
	if header:
		if counts:
			# TODO add taxid functionality for headers
			print sep.join(['file', 'K_spec_taxid', 'read_count', 'median_score'])
		else:
			print sep.join(['file', 'read_id', 'true_spec_taxid', 'K_spec_taxid', 'score'])
	for r in record:
		print sep.join(record)
