import frakka.utils
import argparse
import sys

VERSION = '0.0.1'

def set_parsers():

	parser = argparse.ArgumentParser(description='frakka - a tool to filter Kraken output files and calculate read-level\
		and summary confidence score metrics per classified species',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)
	parser.add_argument('--kreport', '-k', help='Kraken2 report file(s) (comma-separated if multiple, required)')
	parser.add_argument('--kout', '-o', help='Kraken2 output file(s) (comma-separated if multiple, required, must be in same order as --kreport / -k)')
	parser.add_argument('--species', '-sp', help='List of corresponding true / known species names or taxids\
		(comma-separated if multiple, optional, must be in same order as -k and -o files)')
	parser.add_argument('--fof', '-f', help='Tab-separated file of one (-k, -o, -sp)-tuple per line')
	parser.add_argument('--score', '-s', help='Confidence score threshold, only reads / counts higher than this score are reported', default=0)
	parser.add_argument('--counts', '-c', action='store_true', default=False, help='Report total counts per species with score > --score / -s instead of per-read reporting')
	parser.add_argument('--taxid', '-t', action='store_true', default=False, help='Species input and output are NCBI taxIDs instead of species names')
	# TODO parser.add_argument('--plot', '-p', action='store_true', default=False, help='Plot distribution of score per species'

	args = parser.parse_args()
	print(f'arguments: {args}')

	if vars(args) == {}:
		# FIXME it will never be empty because defaults are set
		parser.print_help(sys.stderr)
	else:
		return args

def main():
	args = set_parsers()

	# TODO Steps that need to be run
	# TODO depend on args that are set
	# TODO catch cases of args not set that are required
	# Use args.full_name

	# msg(*args, **kwargs):
	# err(*args, **kwargs):
	# checkPath(path):
	# readTabSep(file, type):
	# readKraken(file):
	# readKReport(file):
	# fileOfFiles(file):
	# getConfidence(taxid, kmerstr):
	# getCounts(krak_f, report_f):
	# output(record, name, header=True, sep='\t', counts=False, taxid=False):


if __name__ == "__main__":
    main()
