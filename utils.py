import os
import sys
from datetime import datetime
from statistics import median as median

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

class DirHandler:
	'''A class for checking / creating output directories'''

	def __init__(self, argsdir):
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
				raise ValueError(f'INTERNAL ERROR: Kraken file {self.file} provided but not cut-off score specified!')
			self.score = score

		self.type = ftype

	def _checkPath(self):
		'''A function to check that a FileReader object has a correct path'''
		path = '/'.join([self.cwd, self.file])
		if not os.path.exists(path):
			Logger.err(f'File {path} does not exist! Exiting.')
		return path

	def _readTabSep(self, msg):
		'''
		checks if tab-separated file exists and reads it into a list of lists
		which it returns
		'''

		path = self._checkPath()

		with open(path, 'r') as f:
			Logger.msg(f'Reading {msg} file {self.file}')
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
				Logger.msg(f'Finished reading file {file}. All entries (n = {coldict[list(coldict.keys())[0]]}) have the same number of columns (n = {list(coldict.keys())[0]})')
			else:
				Logger.msg('File of files {file} contains different number of columns per line.')
				Logger.msg('The following number of lines with different number of columns have been detected:')

				for k in sorted(coldict, reverse=True):
					print(f'\t\tcoldict[k] lines with {k} columns', file=sys.stderr)
				Logger.err('Check your file is in the correct format and try again.')

			return file_lst, col

	def fileOfFiles(self):
		'''
		Reads a tab-separated file of report and output-file pairs specified by --fof / -f (and optionally, a third column with species identifiers)
		and returns a list of lists to iterate over
		'''
		f = self._readTabSep('file-of-files')[0]

		return f

	def readKraken(self):
		'''
		Reads Kraken standard output file (STDOUT / --output) 
		and returns a dict of all classified reads with read id as key
		and a dict as value (with keys taxid:, species: and kmerstr:)	
		'''
		if self.ftype != 'krak':
			raise TypeError(f'INTERNAL ERROR: Cannot call method readKraken() on a file with ftype {self.ftype}!')

		krak = []

		f = self._readTabSep('Kraken output')[0]
		for line in f:
			kl = KrakenLine(uc=line[0], read_id=line[1], taxid=line[2], kmerstr=line[4]) 
			if kl.uc != 'U' and kl.score >= self.score:
				krak.append(kl)

		# Return a list of KrakenLine objects
		return krak
		
	def readKReport(self):
		'''
		Reads in a Kraken report file created via the --report option
		and returns a dict of taxids as keys matched to a dict of name: and count:
		'''
		if self.ftype != 'rep':
			raise TypeError(f'INTERNAL ERROR: Cannot call method readKReport() on a file with ftype {self.ftype}! ftype "rep" expected!')

		report = {}
		f, form  = self._readTabSep('Kraken report')
		level = form - 3
		tax_col = form - 2
		spec_col = form - 1
		for line in f:
			if line[level] == 'S':
				taxid, species = line[tax_col], line[spec_col]
				report[taxid] = species

		return report

class Counter:

	@classmethod
	def getCounts(cls, specmap, kll, truespec, file):
		'''
		Implements the --counts option to return counts for each species
		instead of confidence scores for individual reads
		specmap is a return object from readKReport() call, mapping taxid to species
		kll is a return object from readKraken() call
		FIXME how do we pass file name?
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
				except KeyError:
					if specmap.get(kl.taxid, None):
							Logger.err('Another KeyError while discarding non-species reads has occurred that should not occur. Exiting.')

		for tx in species:
			read_count = species[tx]['read_count']
			median_score = round(median(species[tx]['scores']), 3)
			cr = CountRecord(file=file, truespec=truespec, kspec=tx, read_count=read_count, median_score=median_score)
			result.append(cr)

		return result

class KrakenLine:
	'''A class representing a single line of Kraken output file'''
	def __init__(self, uc, read_id, taxid, kmerstr)
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
		# TODO we must be able to get the truespec from the user input
		self.kspec = kspec

class CountRecord(Record):
	'''A class that represents a species-read_count object'''
	def __init__(self, file, truespec, kspec, read_count, median_score):
		super().__init__(file, truespec, kspec)
		self.read_count = read_count
		self.median_score = median_score

	def join(self, sep):
		return sep.join([self.file, self.truespec, self.kspec, self.read_count, self.median_score])

class ReadRecord(Record):
	'''A class that represents a species-read-kmerstring-confidence score object'''
	def __init__(self, file, truespec, kspec, read_id, score):
		super().__init__(file, truespec, kspec)
		self.read_id = read_id
		self.score = score

	def join(self, sep):
		return sep.join([self.file, self.truespec, self.kspec, self.read_id, self.score])

class Output:
	'''A class for output objects (lines to be printed)'''

	def __init__(self, record, header=True, sep='\t', counts=False, taxid=False):
		self.isHeader = isHeader
		self.useTaxid = useTaxid
		self.sep = sep

		if self.isHeader:
			self.kspec = "K_spec"
			self.truespec = "true_spec"

			if self.useTaxid:
				self.kspec = self.kspec + "_taxid"
				self.truespec = self.truespec + "_taxid"

			if counts:
				self.record = self.sep.join(['file', self.truespec, self.kspec, 'read_count', 'median_score'])
			else:
				self.record = self.sep.join(['file', self.truespec, self.kspec, 'read_id', 'score'])
		else:
			self.record = record.join(self.sep)

	def printRecord(self):
		print(self.record)
