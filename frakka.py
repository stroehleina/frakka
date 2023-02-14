from utils import Logger, DirHandler, FileReader, Counter, Output
from plot import CountPlotter
from datetime import datetime
import argparse
import sys
import re
import os

VERSION = '0.0.1'

def set_parsers():
	'''Sets command line argument options and parses them'''
	parser = argparse.ArgumentParser(description='frakka - a tool to filter Kraken output files and calculate read-level\
		and summary confidence score metrics per classified species',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)
	parser.add_argument('--kreport', '-k', help='Kraken2 report file(s) (comma-separated if multiple, required)')
	parser.add_argument('--kout', '-o', help='Kraken2 output file(s) (comma-separated if multiple, required, must be in same order as --kreport / -k)')
	parser.add_argument('--species', '-sp', help='List of corresponding true / known species names or taxids\
		(comma-separated if multiple, optional, must be in same order as -k and -o files,\
		single species provided assumes all files are same true / known species)')
	parser.add_argument('--sp_only', '-spo', help='NOT IMPLEMENTED Report only matches to true/known species (--species) option') # TODO
	parser.add_argument('--fof', '-f', help='Tab-separated file of one (-k, -o, -sp)-tuple per line')
	parser.add_argument('--score', '-s', help='Confidence score threshold, only reads / counts higher than this score are reported', default=0)
	parser.add_argument('--counts', '-c', action='store_true', default=False, help='Report total counts per species with score > --score / -s instead of per-read reporting')
	parser.add_argument('--taxid', '-t', action='store_true', default=False, help='Species input and output are NCBI taxIDs instead of species names')
	parser.add_argument('--plot', '-p', action='store_true', default=False, help='Plot distribution of score per species')
	parser.add_argument('--directory', '-d', default=f'frakka_{re.sub(":", "_", datetime.now().isoformat(sep="_", timespec="seconds"))}', help='Specify output directory')
	parser.add_argument('--prefix', '-x', default='.', help='NOT IMPLEMENTED Specify output file prefix') # TODO 
	parser.add_argument('--delim', '-del', default='\t', help='Specify output file delimiter')

	args = parser.parse_args()
	return args, parser

def main():

	msg = Logger.msg
	err = Logger.err

	args, parser = set_parsers()
	file_s = []
	k_files = {args.kreport, args.kout}

	# Creating or checking output directory
	# TODO catch outdir == "-" --> textfiles to STDOUT and plots to cwd
	outdir = DirHandler(args.directory).makeOutputDir()

	try:
		k_files.remove(None)
	except KeyError:
		pass

	if (len(k_files) < 2 and not args.fof) or (len(k_files) == 2 and args.fof):
		msg('You need to specify either both --kreport and --kout or alternatively, --fof. Exiting.')
		parser.print_help(sys.stderr)
		
	if args.fof:
		fof = FileReader(cwd=os.getcwd(), file=args.fof, ftype='fof')
		file_s = fof.fileOfFiles()
		if args.species and len(file_s[0]) == 3:
			err('Species provided in column 3 of file specified via --fof but also via --species. Provide one or the other. Exiting')
	else:
		f_kreport = args.kreport.split(',')
		f_kout = args.kout.split(',')
		f_species = ['N/A']

		if args.species:
			f_species = args.species.split(',')

		same_spec = False
		if len(f_species) == 1:
			same_spec = True

		if same_spec and len(f_kreport) > 1 and f_species[0] != 'N/A':
			msg(f'Only one species provided but {len(f_kreport)} report files. Assuming true / known species is {f_species[0]} for all reports.')
			
		if (len(f_kreport) != len(f_kout)) or (all(l != len(f_species) for l in (len(f_kout), len(f_kreport))) and (not same_spec)):
			err('The number of provided comma-separated files for --kreport and --kout do not match')

		for i,f in enumerate(f_kreport):
			try:
				file_lst = [f, f_kout[i], f_species[i]]
			except IndexError:
				file_lst = [f, f_kout[i], f_species[0]]

			file_s.append(file_lst)

	outlst = []

	# process individual file records
	for f in file_s:
		rec_lst = []
		truespec = f[2]

		if ((args.species and re.match('[^0-9]', args.species)) or re.match('[^0-9]', truespec)) and truespec != 'N/A' and args.taxid:
			err('--taxid has been provided but input via --species or --fof (column 3) contain non-numerical characters.\
			Did you accidentally provide species names instead?\
			Check your file or remove the --taxid argument.')	

		report = FileReader(cwd=os.getcwd(), file=f[0], ftype='rep')
		specmap = report.readKReport()

		krak = FileReader(cwd=os.getcwd(), file=f[1], ftype='krak', score=float(args.score))
		kll = krak.readKraken()

		if args.counts:
			# Counter.getCounts() returns a list of CountRecord objects
			counts = Counter.getCounts(specmap=specmap, kll=kll, truespec=truespec, file=f[1])

			for cr in counts:
				if not args.taxid:
					cr.kspec = cr.species

				rec_lst.append(cr.join(args.delim))

		else:

			# Create a ReadRecord object for each KrakenLine in kll
			for kl in kll:
				kspec = specmap[kl.taxid]

				if args.taxid:
					kspec = kl.taxid

				rec = ReadRecord(file=f[1], truespec=truespec, kspec=kspec, read_id=kl.read_id, score=kl.score)
				rec_lst.append(rec.join(args.delim))

		outlst += rec_lst


	for c,o in enumerate(outlst):
		if c == 0:
			header = Output(record=None, isHeader=True, sep=args.delim, counts=args.counts, useTaxid=args.taxid)
			header.printRecord()

		out = Output(record=o, sep=args.delim, useTaxid=args.taxid)
		out.printRecord()
		# TODO define output directory for printing output lines

	# TODO
	# if args.plot:

	# 	if args.counts:
	# 		# Create plots for summarised counts
	# 		# FIXME create object
	# 		CountPlotter.plot()
	# 		# FIXME when supplying multiple files we must plot this once per file
	# 		# iterate over files (we have a file name for this)
	# 		# iterate over line objects (for printing output)
	# 		# pass list of line objects to this function (we don't need to iterate over it in here)


	# 	else:
	# 		# Create distribution plot for reads


	msg('Done. Thank you for using frakka. Please cite https://github.com/stroehleina/frakka')

if __name__ == "__main__":
    main()
