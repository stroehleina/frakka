import os
import sys
from datetime import datetime
from statistics import median as median

# Log a message to stderr
def msg(*args, **kwargs):
	now = f'[{datetime.now().isoformat(sep=" ", timespec="seconds")}]: '
	print(now, *args, file=sys.stderr, **kwargs)


# Log an error to stderr and quit with non-zero error code
def err(*args, **kwargs):
    msg('ERROR', *args, **kwargs)
    sys.exit(1)


def _checkPath(path):
	if not os.path.exists(path):
		err(f'File {path} does not exist! Exiting.')


def _readTabSep(file, type):
	'''
	checks if tab-separated file exists and reads it into a list of lists
	which it returns
	'''

	_checkPath(file)

	with open(file, 'r') as f:
		msg(f'Reading {type} file {file}')
		file_lst = []
		col = None
		coldict = {}
		for line in f:
			cols = line.rstrip('\n').split('\t')
			cols = [c.strip() for c in cols]
			col = len(cols)

			try:
				coldict[col] += 1				
			except KeyError:
				coldict[col] = 1

			file_lst.append(cols)

		if len(coldict) == 1:
			msg(f'Finished reading file {file}. All entries (n = {coldict[list(coldict.keys())[0]]}) have the same number of columns (n = {list(coldict.keys())[0]})')
		else:
			msg('File of files {file} contains different number of columns per line.')
			msg('The following number of lines with different number of columns have been detected:')
			for k in sorted(coldict, reverse=True):
				print(f'\t\tcoldict[k] lines with {k} columns', file=sys.stderr)
			err('Check your file is in the correct format and try again.')

		return file_lst, col 


def readKraken(file, score=0):
	'''
	Reads Kraken standard output file (STDOUT / --output) 
	and returns a dict of all classified reads with read id as key
	and a dict as value (with keys taxid:, species: and kmerstr:)	
	'''

	krak = {}

	f = _readTabSep(file, 'Kraken output')[0]
	for line in f:
			uc, read, taxid, kmerstr = line[0], line[1], line[2], line[4]
			# msg(f'Calculating confidence score for read {read} (kmerstring: {kmerstr}):')
			conf = _getConfidence(taxid, kmerstr)
			if uc != 'U' and conf >= score:
				krak[read]={'taxid' : taxid, 'conf' : conf}

	return krak


def readKReport(file):
	'''
	Reads in a Kraken report file created via the --report option
	and returns a dict of taxids as keys matched to a dict of name: and count:
	'''

	report = {}

	f,form  = _readTabSep(file, 'Kraken report')
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
	f = _readTabSep(file, 'tab-separated input')[0]

	return f


def _getConfidence(taxid, kmerstr):
	'''
	Reads the kmerstring and called taxid from column 5 and 3 of the Kraken standard output file (STDOUT / --output) 
	and returns the confidence score for the read
	see: https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring

	'''

	# Format: 562:13 561:4 A:31 0:1 562:3
	# paired read data will contain a '|:|' token in this list to indicate the end of one read and the beginning of another
	r_lst = kmerstr.split(' |:| ')
	assert 1 <= len(r_lst) <= 2
	conf_lst = []
	for r in r_lst:
		conf = 0
		allkmers = 0
		taxkmers = 0
		cons_kmers = r.strip().split(' ')

		for k in cons_kmers:
			run = k.split(':')
			if run[0] != 'A':
				allkmers += int(run[1])

				if run[0] == taxid:
					taxkmers += int(run[1])

		if allkmers != 0:
			conf_lst.append(taxkmers/allkmers)

	if 1 <= len(conf_lst) <= 2:
		conf = round(sum(conf_lst)/len(conf_lst), 3)

	return conf


def getCounts(report_f, krak_f, score=0):
	'''
	Implements the --counts option to return counts for each species
	instead of confidence scores for individual reads
	'''

	species = {}
	res_lst = []

	report = readKReport(report_f)
	krak = readKraken(krak_f, score=score)

	for read in krak:
		taxid=krak[read]['taxid']
		conf = krak[read]['conf']
		try:
			species[taxid]['read_count'] += 1
			species[taxid]['confs'].append(conf)
		except KeyError:
			species[taxid]['read_count'] = 1
			species[taxid]['confs'] = [conf]
			species[taxid]['name'] = report[krak[read]['taxid']]['name']

	for taxid in species:
		species[taxid]['median_score'] = median(species[taxid]['confs'])

	return species


def output(record, header=True, sep='\t', counts=False, taxid=False):
	msg(f"Writing output for file {record['name']}")

	if header:
		kspec = "K_spec"
		truespec = "true_spec"

		if taxid:
			kspec = kspec + "_taxid"
			truespec = truespec + "_taxid"

		hdlst = []

		if counts:
			hdlst = ['file', truespec, kspec, 'read_count', 'median_score']
		else:
			hdlst = ['file', truespec, kspec, 'read_id', 'score']

		print(sep.join(hdlst))

	for r in record['records']:
		print(sep.join([str(x) for x in r]))
