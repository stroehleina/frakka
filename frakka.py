import frakka.utils
import argparse
import sys
import re

VERSION = '0.0.1'

def set_parsers():

	parser = argparse.ArgumentParser(description='frakka - a tool to filter Kraken output files and calculate read-level\
		and summary confidence score metrics per classified species',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)
	parser.add_argument('--kreport', '-k', help='Kraken2 report file(s) (comma-separated if multiple, required)')
	parser.add_argument('--kout', '-o', help='Kraken2 output file(s) (comma-separated if multiple, required, must be in same order as --kreport / -k)')
	parser.add_argument('--species', '-sp', help='List of corresponding true / known species names or taxids\
		(comma-separated if multiple, optional, must be in same order as -k and -o files,\
		single species provided assumes all files are same true / known species)')
	parser.add_argument('--fof', '-f', help='Tab-separated file of one (-k, -o, -sp)-tuple per line')
	# TODO parser.add_argument('--score', '-s', help='Confidence score threshold, only reads / counts higher than this score are reported', default=0)
	parser.add_argument('--counts', '-c', action='store_true', default=False, help='Report total counts per species with score > --score / -s instead of per-read reporting')
	parser.add_argument('--taxid', '-t', action='store_true', default=False, help='Species input and output are NCBI taxIDs instead of species names')
	# TODO parser.add_argument('--plot', '-p', action='store_true', default=False, help='Plot distribution of score per species')
	# TODO parser.add_argument('--directory', '-d', default='.', help='Specify output directory')
	# TODO parser.add_argument('--prefix', '-x', default='.', help='Specify output file prefix')
	parser.add_argument('--delim', '-t', default='\t', help='Specify output file delimiter')


	args = parser.parse_args()

	return args

def main():

	args = set_parsers()
	file_s = []

	if len({a: args[a] for a in ('kreport', 'kout')}) < 2 or not args.fof:
		parser.print_help(sys.stderr)
		err('You need to specify either both --kreport and --kout or alteratively, --fof. Exiting.')
		
	if args.fof:
		file_s = fileOfFiles(args.fof)
		if args.species and len(file_s[0]) == 3:
			err('Species provided in column 3 of file specified via --fof but also via --species. Provide one or the other. Exiting')
	else:
		f_kreport = args.kreport.split(',')
		f_kout = args.species.split(',')
		f_species = ['N/A']

		if args.species:
			f_species = args.species.split(',')

		same_spec = False

		if len(f_species) == 1 and len(f_kreport) > 1:
			msg(f'Only one species provided but {len(f_kreport)} report files. Assuming true / known species is {f_species[0]} for all reports.')
			same_spec = True

		if len(f_kreport) != len(f_kout) or ((len(f_kreport) != len(f_species) or len(f_kout) != len(f_species)) and same_spec):
			err('The number of provided comma-separated files for --kreport and --kout do not match')

		for i,f in f_kreport:
			try:
				file_lst = [f, f_kout[i], f_species[i]]
			except IndexError:
				file_lst = [f, f_kout[i], f_species[0]]

			file_s.append(file_lst)

	# List of output records (dicts)
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
			counts = getCounts(report, out)

			# One line per species
			for taxid in counts:

				if args.taxid:
					kspec = taxid
				else:
					kspec = counts[taxid]['name']

				rec = [report, truespec, kspec, counts[taxid]['read_count'], counts[taxid]['median_score']]
				rec_lst.append(rec)
		else:
			report = readKReport(report)
			# report[taxid] = {'count' : count, 'name' : name}

			krak = readKraken(out)
			# krak[read]={'taxid' : taxid, 'conf' : conf}

			for read in krak:
				kspec = report[krak[read]['taxid']]['name']
				conf = krak[read]['conf']

				if args.taxid:
					kspec = krak[read]['taxid']


				rec = [f[0], truespec, kspec, read, conf]
				rec_lst.append(rec)

		outrec = {'name' : f[0].split('/')[-1:][0] 'records' : rec_lst}

		outlst.append(outrec)


	for c,o in enumerate(outlst):
		if c == 0:
			output(o, header=True, sep=args.delim, counts=args.counts, taxid=args.taxid)
		else:
			output(o, header=False, sep=args.delim, counts=args.counts, taxid=args.taxid)


if __name__ == "__main__":
    main()
