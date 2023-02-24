import os
import sys
from datetime import datetime
from statistics import median

class Logger:
	'''A class that provides basic logging functions'''

	@classmethod
	def msg(cls, *args, **kwargs):
		'''Log a message to stderr'''
		now = f'[{datetime.now().isoformat(sep=" ", timespec="seconds")}]: '
		print(now, *args, file=sys.stderr, **kwargs)

	@classmethod
	def err(cls, *args, **kwargs):
		'''Log an error to stderr and quit with non-zero error code'''
		cls.msg('ERROR', *args, **kwargs)
		sys.exit(1)

msg = Logger.msg
err = Logger.err

class DirHandler:
	'''A class for checking / creating output directories'''

	def __init__(self, argsdir):
		if argsdir[0] == "/":
			self.dir = argsdir
		else:
			self.dir = '/'.join([os.getcwd(), argsdir])
		# TODO make this a list of dirs so we can create a separate directory for each pair of files?

	def makeOutputDir(self):
		'''Checks if plot output directory exists and creates one if it doesn't'''
		if os.path.exists(self.dir) and os.path.isdir(self.dir):
			msg(f'Found output directory {self.dir}. Will save plots there.')
			return self.dir

		else:
			if os.path.exists(self.dir) and not os.path.isdir(self.dir):
				err(f'The path {self.dir} for the output directory you have provided exists but is a file, not a directory. Exiting')
			try:
				msg(f'Creating output directory {self.dir}.')
				os.mkdir(self.dir)
				return self.dir
			except PermissionError:
				err(f'You do not have permission to create directory {self.dir}. Exiting.') 

class FileReader:
	'''A class that provides functions for reading Kraken files'''
	def __init__(self, cwd, file, ftype, **kwargs):
		self.cwd = cwd
		self.file = file

		f_types = ['rep', 'krak', 'fof']
		if ftype not in f_types:
			raise ValueError('INTERNAL ERROR: Invalid file type. Expected one of: %s' % f_types)

		if ftype == 'krak':
			if 'score' not in kwargs.keys():
				raise ValueError(f'INTERNAL ERROR: Kraken file {self.file} provided but no cut-off score specified!')
			self.score = kwargs['score']

		self.ftype = ftype

	def _checkPath(self):
		'''A function to check that a FileReader object has a correct path'''
		if self.file[0] == "/":
			path = self.file
		else:
			path = '/'.join([self.cwd, self.file])

		if not os.path.exists(path):
			err(f'File {path} does not exist! Exiting.')
		return path

	def _readTabSep(self, message):
		'''
		checks if tab-separated file exists and reads it into a list of lists
		which it returns
		'''

		path = self._checkPath()

		with open(path, 'r') as f:
			msg(f'Reading {message} file {self.file}')
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
				msg(f'Finished reading file {path}. All entries (n = {coldict[list(coldict.keys())[0]]}) have the same number of columns (n = {list(coldict.keys())[0]})')
			else:
				msg(f'File of files {path} contains different number of columns per line.')
				msg('The following number of lines with different number of columns have been detected:')

				for k in sorted(coldict, reverse=True):
					print(f'\t\tcoldict[k] lines with {k} columns', file=sys.stderr)
				err('Check your file is in the correct format and try again.')

			return file_lst, col

	def fileOfFiles(self):
		'''
		Reads a tab-separated file of report and output-file pairs specified by --fof / -f (and optionally, a third column with species identifiers)
		and returns a list of lists to iterate over
		'''
		f = self._readTabSep(message='file-of-files')[0]

		return f

	def readKraken(self):
		'''
		Reads a read-level Kraken output file 
		and returns a list of KrakenLine objects	
		'''
		if self.ftype != 'krak':
			raise TypeError(f'INTERNAL ERROR: Cannot call method readKraken() on a file with ftype {self.ftype}!')

		krak = []

		f = self._readTabSep(message='Kraken output')[0]
		for line in f:
			kl = KrakenLine(uc=line[0], read_id=line[1], taxid=line[2], kmerstr=line[4]) 
			if kl.uc != 'U' and kl.score >= self.score:
				krak.append(kl)

		# Returns a list of KrakenLine objects
		return krak
		
	def readKReport(self):
		'''
		Reads in a Kraken report file created via the --report option
		and returns a dict of taxids as keys matched to a dict of name: and count:
		'''
		if self.ftype != 'rep':
			raise TypeError(f'INTERNAL ERROR: Cannot call method readKReport() on a file with ftype {self.ftype}! ftype "rep" expected!')

		report = {}
		gtdb = {} # We fill this with S1 lines instead of S lines in case we are dealing with a GTDB input file
		f, form  = self._readTabSep(message='Kraken report')
		level = form - 3
		tax_col = form - 2
		spec_col = form - 1
		check_S_gtdb = set()
		for line in f:
			if line[level] == 'S':
				taxid, species = line[tax_col], line[spec_col]
				report[taxid] = species
				# If we are dealing with a GTDB file the assumption is that the unique counts at S level are all 0 and the species names are in S1 lines
				# We store these in a set and take the union everytime we add one then later check that there is only one element in it and that is 0
				check_S_gtdb = check_S_gtdb.union({line[2]})
			elif line[level] == 'S1':
				taxid, species = line[tax_col], line[spec_col]
				gtdb[taxid] = species

		if {'0'} == check_S_gtdb and len(gtdb) == len(report):
			msg(f'WARNING: It appears that your input files are derived from the GTDB database, using S1 lines of the Kraken report to match taxids to species names.')
			# keep the keys from gtdb but keep the values from report.
			# We can do this because dicts are ordered since Python 3.7 (and assume S1 lines always follow S line in the Kraken report)
			gtdb_keys = list(gtdb)
			rep_values = list(report.values())
			report = {}
			for i,k in enumerate(gtdb_keys):
				report[k] = rep_values[i]

		return report

class Counter:

	@classmethod
	def getCounts(cls, specmap, kll, truespec, file):
		'''
		Implements the --counts option to return counts for each species
		instead of confidence scores for individual reads
		specmap is a return object from readKReport() call, mapping taxid to species
		kll is a return object from readKraken() call
		'''

		species = {}
		result = []

		for kl in kll:

			try:
				species[kl.taxid]['read_count'] += 1
				species[kl.taxid]['scores'].append(kl.score)
			except KeyError:
				try:
					species[kl.taxid] = {
										'read_count' : 1,
										'scores' : [kl.score],
										'name' : specmap[kl.taxid]
										}
				except KeyError as e:
					if specmap.get(kl.taxid, None):
							err(f'Another KeyError {e.args[0]} while discarding non-species reads has occurred that should not occur. Exiting.')

		for tx in species:
			name = species[tx]['name']
			read_count = species[tx]['read_count']
			median_score = round(median(species[tx]['scores']), 3)
			# NOTE do we need species name and kspec? kspec is taxid so we should never have to override it if we have species
			cr = CountRecord(file=file, truespec=truespec, kspec=tx, species=name, read_count=read_count, median_score=median_score)
			result.append(cr)

		return result

class KrakenLine:
	'''A class representing a single line of Kraken output file'''
	def __init__(self, uc, read_id, taxid, kmerstr):
		self.uc = uc
		self.read_id = read_id
		self.taxid = taxid
		self.kmerstr = kmerstr
		self.score = self.getConfidence()

	def getConfidence(self):
		'''
		Reads the kmerstring and called taxid from column 5 and 3 of the Kraken standard output file (STDOUT / --output) 
		and returns the confidence score for the read
		see: https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring
		'''

		# Format: 562:13 561:4 A:31 0:1 562:3
		# paired read data will contain a '|:|' token in this list to indicate the end of one read and the beginning of another
		r_lst = self.kmerstr.split(' |:| ')
		assert 1 <= len(r_lst) <= 2
		conf_lst = []
		score = 0
		for r in r_lst:
			allkmers = 0
			taxkmers = 0
			cons_kmers = r.strip().split(' ')

			for k in cons_kmers:
				run = k.split(':')
				if run[0] != 'A':
					allkmers += int(run[1])

					if run[0] == self.taxid:
						taxkmers += int(run[1])

			if allkmers != 0:
				conf_lst.append(taxkmers/allkmers)

		if 1 <= len(conf_lst) <= 2:
			score = round(sum(conf_lst)/len(conf_lst), 3)

		return score

class Record:
	'''A generic super-class for frakka record objects (CountRecord and ReadRecord)'''
	def __init__(self, file, truespec, kspec):
		self.file = file
		self.truespec = truespec
		self.kspec = kspec

class CountRecord(Record):
	'''A class that represents a species-read_count object'''
	def __init__(self, file, truespec, kspec, species, read_count, median_score=0):
		super().__init__(file, truespec, kspec)
		self.species = species
		self.read_count = read_count
		self.median_score = median_score # NOTE Set this as default 0 

	def join(self, sep):
		return sep.join([self.file, self.truespec, self.kspec, str(self.read_count), str(self.median_score)])

class ReadRecord(Record):
	'''A class that represents a species-read-kmerstring-confidence score object'''
	def __init__(self, file, truespec, kspec, read_id, score):
		super().__init__(file, truespec, kspec)
		self.read_id = read_id
		self.score = score

	def join(self, sep):
		return sep.join([self.file, self.truespec, self.kspec, str(self.read_id), str(self.score)])

class Output:
	'''A class for output objects (lines to be printed)'''

	def __init__(self, record, fh, isHeader=False, sep='\t', counts=False, useTaxid=False):
		self.fh = fh
		self.isHeader = isHeader
		self.sep = sep
		self.counts = counts
		self.useTaxid = useTaxid

		if self.isHeader:
			self.kspec = "K_spec"
			self.truespec = "true_spec"

			if self.useTaxid:
				self.kspec = self.kspec + "_taxid"
				self.truespec = self.truespec + "_taxid"

			if self.counts:
				self.record = self.sep.join(['file', self.truespec, self.kspec, 'read_count', 'median_score'])
			else:
				self.record = self.sep.join(['file', self.truespec, self.kspec, 'read_id', 'score'])
		else:
			self.record = record

	def printRecord(self):
		print(self.record, file=self.fh)
