from utils import Logger
from utils import CountRecord
# from utils import ReadRecord # FIXME
# from utils import Output # FIXME is this needed if all we do is pass around Output objects but do not create new ones? 
from operator import itemgetter
import matplotlib.pyplot as plt
import sys
import os

# Importing messaging and error methods from FrakkaUtils
msg=Logger.msg
err=Logger.err

class Plotter:
	'''A class for creating plots from frakka output objects'''
	def __init__(self, dest):
		self.dest = dest

class CountPlotter(Plotter):
	'''A Plotter to plot per-species confidence summaries'''
		def __init__(self, outdir, file, counts, other_co=0, drop=0, cutoff=0):
			super().__init__(outdir)
			self.counts = counts
			self.file = file
			self.counts = counts
			self.other_co = other_co
			self.drop = drop
			self.cutoff = cutoff

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

		if other_co != 0:
			o_obj = CountRecord(file=self.file, truespec='Other', kspec='Other', read_count=other_sum, median_score='')
			# INFO cannot provide median of the merged "Other" category
			# INFO No, unfortunately there is not a way to calculate the median based on medians of subsets of the whole and still be statistically accurate.
			# INFO If you wanted to calculate the mean, however, you could use the means of subsets, given that they are of equal size.
			self.counts = [self.counts, o_obj]

		plot_lst = []
		# Sort the list of counts objects by read_count
		for c in self.counts:
			plot_lst.append({'name' : c.kspec, 'read_count' : int(c.read_count)})

		# Sort by 'read_count' then by 'name'
		plot_lst.sort(key=operator.itemgetter('read_count'), reverse=True)
		plot_lst.sort(key=operator.itemgetter('name'))
		
		# plt.rcdefaults() # INFO Do we need this?
		ax = plt.subplots()[1]

		# Example data
		# ordered and unchangeable
		species = [d['name'] for d in plot_lst]
		y_pos = [i+1 for i,v in enumerate(species)]
		plt_counts = [d['read_count'] for d in plot_lst]

		ax.barh(y_pos, plt_counts)
		ax.set_yticks(y_pos, labels=species)
		ax.invert_yaxis()
		for i, v in enumerate(plt_counts):
    		ax.text(v + 100, i+1, str(v), ha='left', va='center')
		ax.set_xlabel('Read counts')
		ax.set_title(f'Distribution of reads per Species with confidence score > {self.cutoff}')

		# chop file ending, replace all "." by "__" and add .pdf
		outfile = '__'.join(self.file.split('\.').pop()) + f'_c-o_{self.cutoff}.pdf'
		figpath = '/'.join(self.outdir, outfile)
		msg(f'Saving filtered read counts-per-species plot for file {self.file} to {figpath}')
		plt.savefig(figpath)

class ReadPlotter(Plotter):
	'''A Plotter to plot per-file read confidence distributions'''
	def __init__(self):
		pass

	def plotConfDist(self, confs):
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
