from utils import Logger
from utils import CountRecord
# from utils import ReadRecord # FIXME
# from utils import Output # FIXME is this needed if all we do is pass around Output objects but do not create new ones? 
import operator
import matplotlib.pyplot as plt
import sys
import os
from math import log

# Importing messaging and error methods from FrakkaUtils
msg=Logger.msg
err=Logger.err

class Plotter:
	'''A class for creating plots from frakka output objects'''
	def __init__(self, outdir):
		self.outdir = outdir

class CountPlotter(Plotter):
	'''A Plotter to plot per-species confidence summaries'''
	def __init__(self, outdir, file, counts, other_co=0, drop=0, score=0):
		super().__init__(outdir)
		self.file = file
		self.counts = counts
		self.other_co = other_co
		self.drop = drop
		self.score = score

	def plot(self):
		'''Reads a list of CountsRecord objects and creates a plot from each'''

		# Applying minimum read count (drop species with less than drop reads)
		msg(f'Applying minium read count cut-off. Of {len(self.counts)} species, those with less than {self.drop} reads will not be plotted.')
		self.counts = [c for c in self.counts if c.read_count >= self.drop]
		msg(f'Applied minium read count cut-off. {len(self.counts)} species remaining.')

		# Merging counts for low-abundance species read_counts
		msg(f'Merging species with read count below {self.other_co} into category "Other".')
		other = [c.read_count for c in self.counts if c.read_count < self.other_co]
		other_sum = sum(other)
		msg(f'Merged {len(other)} species (total reads: {other_sum}) into category "Other".')
		self.counts = [c for c in self.counts if c.read_count >= self.other_co]

		if self.other_co != 0:
			cr = CountRecord(file=self.file, truespec='Other', kspec='Other', read_count=other_sum, median_score='')
			# INFO cannot provide median of the merged "Other" category
			# INFO There is not a way to calculate the median based on medians of subsets of the whole and still be statistically accurate.
			# INFO We also cannot use the means of the subsets, given that they are not of equal size.
			self.counts = [self.counts, cr]

		plot_lst = []
		# Sort the list of counts objects by read_count
		for c in self.counts:
			color='darkgrey'

			if c.kspec == '9606' or c.kspec == 'Homo sapiens':
				color='darkred'

			if c.kspec == c.truespec:
				color='darkgreen'

			plot_lst.append({'name' : c.kspec, 'read_count' : int(c.read_count), 'color' : color})

		# Sort by reverse numerical 'read_count' then by 'name'
		plot_lst.sort(key=lambda e: (-e['read_count'], e['name']))

		# Adjust upper x_bound based on entry with highest read count
		max_rc = int(sorted(plot_lst, key=lambda e: e['read_count'], reverse=True)[0]['read_count'])

		# Adjust image size based on longest read name
		max_name = len(sorted(plot_lst, key=lambda e: len(e['name']), reverse=True)[0]['name'])

		plt.rcdefaults()
		ax = plt.subplots(figsize=(5 + 0.1 * max_name, 0.25 * len(plot_lst)))[1]
		ax.set_xbound(lower=1, upper=max_rc)
		plt.xscale('log', base=10)

		species = [d['name'] for d in plot_lst]
		y_pos = [i+1 for i,v in enumerate(species)]
		plt_counts = [d['read_count'] for d in plot_lst]
		plt_colors = [d['color'] for d in plot_lst]

		ax.barh(y_pos, plt_counts, color=plt_colors)
		ax.set_yticks(y_pos, labels=species, fontstyle='italic')
		ax.invert_yaxis()
		ax.set_xlabel('Read counts')
		ax.set_title(f'Distribution of reads per species\nwith confidence score \u2265 {self.score}')
		ax.set_xmargin(0.3)
		
		for i, v in enumerate(plt_counts):
			ax.text(v, i+1, str(v), ha='left', va='center')

		plt.subplots_adjust(left=0.4)
		
		# chop file ending, replace all "." by "__" and add .pdf
		outfile = '__'.join(self.file.split('/')[-1].split('.')[:-1]) + f'_c-o_{self.score}.pdf'

		figpath = '/'.join([self.outdir, outfile])
		msg(f'Saving filtered read counts-per-species plot for file {self.file} to {figpath}')
		plt.savefig(figpath)

class ReadPlotter(Plotter):
	'''A Plotter to plot per-file read confidence distributions'''
	def __init__(self, scores):
		self.scores = scores
		# TODO a list of scores / ReadCount objects

	def plot(self):
		'''Reads a list of ReadsOut objects (which have a confidence score and a species and a taxid)'''

		# FIXME splitting this by species (i.e. all, species_1, species_2, species_3 etc.)
		# TODO allow to filter for input species only

		pass
		# TODO density plot, one for each species, confidence score on x axis density on y axis
		# TODO Median as vertical line
		# TODO Cut-off in title
		# TODO file name in title
		# TODO Species name in title
		# TODO name file according to species + taxid + cut-off
