from utils import Logger
from utils import CountRecord
from math import log, ceil
from statistics import median
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import os

# Importing messaging and error methods from FrakkaUtils
msg=Logger.msg
err=Logger.err

class Plotter:
	'''A class for creating plots from frakka output objects'''
	def __init__(self, outdir, file, other_co, drop, score, prefix):
		self.outdir = outdir
		self.file = file
		self.other_co = other_co
		self.drop = drop
		self.score = score
		self.prefix = prefix

class CountPlotter(Plotter):
	'''A Plotter to plot per-species confidence summaries'''
	def __init__(self, outdir, file, counts, other_co, drop, score, prefix):
		super().__init__(outdir, file, other_co, drop, score, prefix)
		self.counts = counts

	def plot(self):
		'''Reads a list of CountsRecord objects and creates a plot from each'''

		# Applying minimum read count (drop species with less than drop reads)
		msg(f'Applying minimum read count cut-off. Of {len(self.counts)} species, those with less than {self.drop} reads will not be plotted.')
		self.counts = [c for c in self.counts if c.read_count >= self.drop]
		msg(f'Applied minimum read count cut-off. {len(self.counts)} species remaining.')

		# Merging counts for low-abundance species read_counts
		msg(f'Merging species with read count below {self.other_co} into category "Other".')
		other = [c.read_count for c in self.counts if c.read_count < self.other_co]
		other_sum = sum(other)
		msg(f'Merged {len(other)} species (total reads: {other_sum}) into category "Other".')
		self.counts = [c for c in self.counts if c.read_count >= self.other_co]

		if self.other_co != 0:
			cr = CountRecord(file=self.file, truespec='Other', kspec='Other', species='Other', read_count=other_sum)
			# INFO cannot provide median of the merged "Other" category
			# INFO There is not a way to calculate the median based on medians of subsets of the whole and still be statistically accurate.
			# INFO We also cannot use the means of the subsets, given that they are not of equal size.
			self.counts += [cr]

		plot_lst = []
		# Sort the list of counts objects by read_count
		for c in self.counts:
			color='darkgrey'

			if c.kspec == '9606' or c.kspec == 'Homo sapiens':
				color='darkred'

			if c.kspec == c.truespec and c.kspec != 'Other':
				color='darkgreen'

			plot_lst.append({'name' : c.kspec, 'read_count' : int(c.read_count), 'color' : color})

		# Sort by reverse numerical 'read_count' then by 'name'
		plot_lst.sort(key=lambda e: (-e['read_count'], e['name']))

		# Adjust upper x_bound based on entry with highest read count
		max_rc = int(sorted(plot_lst, key=lambda e: e['read_count'], reverse=True)[0]['read_count'])

		# Adjust image size based on longest read name
		max_name = len(sorted(plot_lst, key=lambda e: len(e['name']), reverse=True)[0]['name'])

		plt.rcdefaults()
		# msg(f'590 pixels corresponds to:')
		# msg(f'5 + 0.1 * max_name = {5 + 0.1 * max_name}')
		# msg(f'80275 pixels corresponds to:')
		# msg(f'0.25 * len(plot_lst) = {0.25 * len(plot_lst)}')

		# maximum is 65536 pixels
		# NOTE reduced dpi?

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
		ax.set_title(f'Distribution of reads per species\nwith confidence score \u2265 {self.score}\n(file {self.file.split("/")[-1]})')
		ax.set_xmargin(0.3)
		
		for i, v in enumerate(plt_counts):
			ax.text(v, i+1, str(v), ha='left', va='center')

		plt.subplots_adjust(left=0.4)
		
		# chop file ending, replace all "." by "__" and add .pdf
		outfile = self.prefix + 'counts_' + '__'.join(self.file.split('/')[-1].split('.')[:-1]) + f'_c-o_{self.score}.pdf'

		figpath = '/'.join([self.outdir, outfile])
		msg(f'Saving filtered read counts-per-species plot for file {self.file} to {figpath}')
		plt.savefig(figpath)

class ReadPlotter(Plotter):
	'''A Plotter to plot per-file (handled in main) read confidence distributions'''
	def __init__(self, outdir, file, rcl, other_co, drop, score, prefix):
		super().__init__(outdir, file, other_co, drop, score, prefix)
		# a list of ReadRecord objects
		self.rcl = rcl
		
	def plot(self):
		'''Reads a list of ReadRecord objects (which have a confidence score and a species and a taxid)'''

		# Non-redundant set of all species present in the list of read records
		# species = set([rc.kspec for rc in rcl]) # don't need this we will just fill dynamically by iterating over all ReadRecords

		species = {}

		for rc in self.rcl:
			# file, truespec, kspec, read_id, score
			# TODO allow to filter for truespec only
			# TODO color addition for truespec
			try:
				species[rc.kspec]['scores'].append(rc.score)
			except KeyError:
				species[rc.kspec] = {'scores' : [rc.score]}

		for s in list(species):

			# remove species for which the length of the list of values is shorter than self.drop (--minreads)
			if len(species[s]['scores']) <= self.drop:
				# msg(f'Length of the score list for species {s} is {len(species[s]["scores"])}. Deleting record.')
				del species[s]
			else:
				try:
					species['All']['scores'] += species[s]['scores']
				except KeyError:
					species['All'] = {'scores' : species[s]['scores']}

			# concatenate all lists of values from species for which the length of the list is shorter than other_co=0
			if len(species[s]['scores']) <= self.other_co:
				# species['Other']['scores'] += species[s]['scores'] # INFO removed this after KeyError: 'Other'
				try:
					species['Other']['scores'] += species[s]['scores']
				except KeyError:
					species['Other'] = {'scores' : species[s]['scores']}

			species[s]['median'] = round(median(species[s]['scores']), 2)

		if 'All' in species.keys():
			species['All']['median'] = round(median(species[s]['scores']), 2)
		else:
			err(f'No reads were added for plotting, try reducing the minimal required read count {self.drop} per species (set via --minreads).')

		if 'Other' in species.keys():
			species['Other']['median'] = round(median(species[s]['scores']), 2)

		# FIXME 'All' should always be the first plot
		# FIXME 'All' median does not seem to be correct (too many added? or not enough added?). Could be correct though but double check.

		# sort species by length of scores list
		species = dict(sorted(species.items(), key=lambda e: len(e[1]['scores']), reverse=True))

		# Create long format pandas object and feed that to displot (should automatically create the facetting we want)
		longlst = []
		for s in species:
			for c in species[s]['scores']:
				longlst.append({'species' : s, 'score' : c, 'median' : species[s]['median']})

		plot_df = pd.DataFrame(longlst)

		plot = sns.displot(plot_df, x='score', col='species', kind='kde', col_wrap=3, color='black', facet_kws={'sharey': False, 'sharex' : False}, warn_singular=False)

		title_f = self.file.split('/')[-1]
		plot.fig.subplots_adjust(top=0.95)
		plot.fig.suptitle(f'Confidence score distribution for file {title_f} (score cut-off: {self.score})')
		plot.set(xlim=(0, 1))
		plot.set_titles(col_template = '{col_name}', fontstyle='italic')
		plot.set_axis_labels('Confidence score')
		spec_median = plot_df[["species", "median"]].drop_duplicates()
		axes = plot.fig.axes

		for i,ax in enumerate(axes):
			ax.axvline(spec_median.iloc[i]['median'], color='red')

		# chop file ending, replace all "." by "__" and add .pdf
		outfile = self.prefix + 'reads_' + '__'.join(self.file.split('/')[-1].split('.')[:-1]) + f'_species_distr_{self.score}.pdf'

		figpath = '/'.join([self.outdir, outfile])
		msg(f'Saving filtered distribution plots for {len(axes)} species (including "All" and "Other" categories, if specified) for file {self.file} to {figpath}')
		plot.figure.savefig(figpath)
