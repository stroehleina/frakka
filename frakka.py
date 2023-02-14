from utils import Logger, DirHandler, FileReader, Counter, Output
from plot import CountPlotter
from datetime import datetime
import argparse
import sys
import re

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
	parser.add_argument('--directory', '-d', default=f'frakka_{re.sub(':', '_', datetime.now().isoformat(sep="_", timespec="seconds"))}', help='Specify output directory')
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

	outlst = [] # FIXME 

	# process individual file records
	for f in file_s:
		rec_lst = [] # FIXME

		truespec = f[2]

		if re.match('[^0-9]', args.truespec) and truespec != 'N/A' and args.taxid:
			err('--taxid has been provided but input via --species or --fof (column 3) contain non-numerical characters.\
			Did you accidentally provide species names instead?\
			Check your file or remove the --taxid argument.')	

		report = FileReader(cwd=os.getcwd(), file=f[0], ftype='rep')
		report_objs = report.readKReport() # TODO additional arguments?

		krak = FileReader(cwd=os.getcwd(), file=f[1], ftype='krak')
		krak_objs = krak.readKraken(score=float(args.score)) # TODO additional arguments?	

		if args.counts:
			counts = Counter.getCounts(report_f=report, krak_f=out, score=float(args.score)) # FIXME
			

			# One line per species
			for taxid in counts:
				if args.taxid:
					kspec = taxid
				else:
					kspec = counts[taxid]['name']

				rec = [report, truespec, kspec, counts[taxid]['read_count'], counts[taxid]['median_score']]
				rec_lst.append(rec) # FIXME
		else:


			for read in krak:
				try:
					kspec = report[krak[read]['taxid']]['name']
					conf = krak[read]['conf']

					if args.taxid:
						kspec = krak[read]['taxid']

					rec = [f[0], truespec, kspec, read, conf]
					rec_lst.append(rec) # FIXME
				except KeyError:
					# msg(f'Read {read} classified higher than species level (taxid: {krak[read]["taxid"]}). Skipping.')
					if report.get(krak[read]['taxid'], None):
						err('Another KeyError has occurred that should not occur. Exiting')

		outrec = {'name' : f[0].split('/')[-1:][0], 'records' : rec_lst} # FIXME
		outlst.append(outrec) # FIXME


	for c,o in enumerate(outlst): # FIXME
		if c == 0:
			fu.output(o, header=True, sep=args.delim, counts=args.counts, taxid=args.taxid) # FIXME 
			# FIXME does this not print the first data line?
			# FIXME Or does it print two lines at once (expected behaviour)?
		else:
			fu.output(o, header=False, sep=args.delim, counts=args.counts, taxid=args.taxid) # FIXME 

		# FIXME define output directory for printing output lines
		# 	# FIXME output() was called on a list of records whereas now it should be called on a single element
		# 	for r in record['records']:
		# 		print(sep.join([str(x) for x in r]))


	if args.plot:

		if args.counts:
			# Create plots for summarised counts
			# FIXME create object
			CountPlotter.plot()
			# FIXME when supplying multiple files we must plot this once per file
			# iterate over files (we have a file name for this)
			# iterate over line objects (for printing output)
			# pass list of line objects to this function (we don't need to iterate over it in here)


		else:
			# Create distribution plot for reads


	msg('Done. Thank you for using frakka. Please cite https://github.com/stroehleina/frakka')

if __name__ == "__main__":
    main()
