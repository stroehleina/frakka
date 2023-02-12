from utils import FrakkaUtils as fu
import plot
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
	parser.add_argument('--plot', '-p', action='store_true', default=False, help='NOT IMPLEMENTED Plot distribution of score per species') # TODO 
	parser.add_argument('--directory', '-d', default='.', help='NOT IMPLEMENTED Specify output directory') # TODO 
	parser.add_argument('--prefix', '-x', default='.', help='NOT IMPLEMENTED Specify output file prefix') # TODO 
	parser.add_argument('--delim', '-del', default='\t', help='Specify output file delimiter')

	args = parser.parse_args()
	return args, parser

def main():

	msg = fu.msg
	err = fu.err

	args, parser = set_parsers()
	file_s = []
	k_files = {args.kreport, args.kout}
	
	try:
		k_files.remove(None)
	except KeyError:
		pass

	if (len(k_files) < 2 and not args.fof) or (len(k_files) == 2 and args.fof):
		parser.print_help(sys.stderr)
		err('You need to specify either both --kreport and --kout or alteratively, --fof. Exiting.')
		
	if args.fof:
		file_s = fu.fileOfFiles(args.fof)
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
		report, out, truespec = f

		if re.match('[^0-9]', truespec) and truespec != 'N/A':
			err('--taxid has been provided but input via --species or --fof (column 3) contain non-numerical characters.\
			Did you provide species names instead?\
			Check your file or remove the --taxid argument.')		

		if args.counts:
			counts = fu.getCounts(report_f=report, krak_f=out, score=float(args.score))

			# One line per species
			for taxid in counts:
				if args.taxid:
					kspec = taxid
				else:
					kspec = counts[taxid]['name']

				rec = [report, truespec, kspec, counts[taxid]['read_count'], counts[taxid]['median_score']]
				rec_lst.append(rec)
		else:
			report = fu.readKReport(report)
			krak = fu.readKraken(out, score=float(args.score))

			for read in krak:
				try:
					kspec = report[krak[read]['taxid']]['name']
					conf = krak[read]['conf']

					if args.taxid:
						kspec = krak[read]['taxid']

					rec = [f[0], truespec, kspec, read, conf]
					rec_lst.append(rec)
				except KeyError:
					# msg(f'Read {read} classified higher than species level (taxid: {krak[read]["taxid"]}). Skipping.')
					if report.get(krak[read]['taxid'], None):
						err('Another KeyError has occurred that should not occur. Exiting')

		outrec = {'name' : f[0].split('/')[-1:][0], 'records' : rec_lst}
		outlst.append(outrec)


	for c,o in enumerate(outlst):
		if c == 0:
			fu.output(o, header=True, sep=args.delim, counts=args.counts, taxid=args.taxid)
		else:
			fu.output(o, header=False, sep=args.delim, counts=args.counts, taxid=args.taxid)

	msg('Done. Thank you for using frakka. Please cite https://github.com/stroehleina/frakka')

if __name__ == "__main__":
    main()
