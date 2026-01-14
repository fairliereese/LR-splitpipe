import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
import numpy as np
# from utils import *

def plot_hist(x, **kwargs):
	ax = sns.histplot(x, binwidth=1)

def plot_mm_lines(x, **kwargs):
	col = x.unique().tolist()[0]
	# linker 1
	if 'l1' in col:
		mismatch_lines = [22, 21, 20, 19]
	# linker 2
	elif 'l2' in col:
		mismatch_lines = [30, 29, 28, 27]
	for l in mismatch_lines:
		plt.axvline(l, color='k', linestyle='-', linewidth=1)

def plot_linker_scores(df, oprefix):

	val_vars = ['l1_score', 'l2_score', 'l1_rc_score', 'l2_rc_score']
	cols = ['index'] + val_vars
	temp = df[cols].melt(id_vars='index', value_vars=val_vars)

	g = sns.FacetGrid(temp, col='variable')
	g.map(plot_hist, 'value')
	g.map(plot_mm_lines, 'variable')

	fname = oprefix+'_linker_score_dists.png'
	plt.savefig(fname)

	plt.clf()

def plot_linker_dists(df, oprefix, xlim=None):
	"""
	Plot a histogram of the distance from the end of each read that
	each linker is
	"""
	for l in ['l1_dist', 'l2_dist']:
		ax = sns.displot(df, x=l, kind='hist', binwidth=5)
		if xlim:
			_ = ax.set(xlim=(0,xlim))
		fname = '{}_{}.png'.format(oprefix, l)
		plt.savefig(fname)

# plot heatmap of number of reads recovered with different linker
# mismatch allowances
def plot_linker_heatmap(df, oprefix, how='integer'):

	if how == 'integer':
		m = [0, 1, 2, 3, 4, 5]
		data = [[0 for i in range(len(m))] for j in range(len(m))]
		m_df = pd.DataFrame(data, index=m, columns=m)

		for i in m: # l1
			for j in m: # l2
				l1_min = 22-i
				l2_min = 30-j
				fwd_df = df.loc[(df.l1_score>=l1_min)&(df.l2_score>=l2_min)]
				rev_df = df.loc[(df.l1_rc_score>=l1_min)&(df.l2_rc_score>=l2_min)]
				one_dir_reads = df.loc[list(set(fwd_df.index)^set(rev_df.index))]
				m_df.at[i, j] = len(one_dir_reads.index)
		ax = sns.heatmap(m_df, annot=True)
		_ = ax.set(xlabel='Mismatches/indels allowed in l2',
				   ylabel='Mismatches/indels allowed in l1',
				   title='Reads recovered')

	elif how == 'proportion':
		# proportional allowances
		p = [0, 4, 7, 10, 13, 16]
		data = [[0 for i in range(len(p))] for j in range(len(p))]
		p_df = pd.DataFrame(data, index=p, columns=p)

		for i in p: # l1
			for j in p: # l2
				l1_min = 22-math.ceil((22/100)*i)
				l2_min = 30-math.ceil((30/100)*j)
				fwd_df = df.loc[(df.l1_score>=l1_min)&(df.l2_score>=l2_min)]
				rev_df = df.loc[(df.l1_rc_score>=l1_min)&(df.l2_rc_score>=l2_min)]
				one_dir_reads = df.loc[list(set(fwd_df.index)^set(rev_df.index))]
				p_df.at[i, j] = len(one_dir_reads.index)

		ax = sns.heatmap(p_df, annot=True)
		_ = ax.set(xlabel='Percent mismatches/indels allowed in l2',
				   ylabel='Percent mismatches/indels allowed in l1',
				   title='Reads recovered')

	fname = '{}_{}_linker_mismatch_heatmap.png'.format(oprefix, how)
	fig = plt.gcf()
	fig.set_size_inches(8, 8)
	plt.savefig(fname, bbox_inches='tight')
	plt.clf()

def plot_read_length(df, oprefix, xlim=None):
	"""
	Plot read length distributions after trimming bc / linker
	construct from reads.

	Parameters:
		df (pandas DataFrame): DataFrame from file with suffix
			"_seq_umi_len.tsv" with cols for read length,
			umi length, bc, and umi
		oprefix (str): Output file prefix
	"""
	ax = sns.displot(data=df, x='read_len', color='#CF91A3')
	ax.set(xlabel='Read length', ylabel='Number of reads')
	fname = '{}_read_length_dist.png'.format(oprefix)
	# plt.xlim((0,10000))
	if xlim:
		plt.xlim((0,xlim))
	plt.savefig(fname)
	plt.clf()

def plot_knee_plot(df, oprefix, kind):
	bc_cols = ['bc']
	# only want unique bc/umi combos
	temp = df[bc_cols+['umi']].drop_duplicates()
	# get the number of unique bc/umi combos
	temp = temp[bc_cols+['umi']].groupby(bc_cols).count()
	temp.reset_index(inplace=True)
	temp.rename({'umi':'counts'}, axis=1, inplace=True)
	temp.sort_values(by='counts', ascending=False, inplace=True)
	# plot
	counts = temp['counts'].tolist()
	plt.plot(range(len(counts)),
			counts,
			color='lightgray',
			linewidth=2)
	ax = plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Ranked cells by # UMIs')
	ax.set_ylabel('# UMIs (logscale)')
	ax.set_title(kind)
	if kind == 'Pre-correction':
		kind = 'pre_correction'
	elif kind == 'Post-correction':
		kind = 'post_correction'
	elif kind == 'Illumina':
		kind = 'illumina'
	plt.tight_layout()
	fname = '{}_{}_knee_plot.png'.format(oprefix, kind)
	plt.savefig(fname)
	plt.clf()
	return temp

def plot_sequencing_sat(df, oprefix, kind):
	bc_cols = ['bc']
	# get the number of reads per barcode
	temp1 = df[bc_cols+['umi']].groupby(bc_cols).count()
	temp1.reset_index(inplace=True)
	temp1.rename({'umi':'reads'}, axis=1, inplace=True)
	temp1.sort_values(by='reads', ascending=False, inplace=True)
	# get the number of unique umis per barcode
	temp2 = df[bc_cols+['umi']].drop_duplicates()
	temp2 = temp2[bc_cols+['umi']].groupby(bc_cols).count()
	temp2.reset_index(inplace=True)
	temp2.rename({'umi':'umis'}, axis=1, inplace=True)
	temp2.sort_values(by='umis', ascending=False, inplace=True)
	# merge on barcode
	temp = temp1.merge(temp2, on=bc_cols)
	bins = [i for i in range(0, temp.reads.max(), 1000)]
	temp_reads = temp.reads.values.tolist()
	temp_bins = np.digitize(temp_reads, bins)
	temp['bin'] = temp_bins
	bins.append(temp.reads.max())
	bin_dict = dict([(bin, bin_num) for bin, bin_num in zip([i for i in range(1,len(bins))], bins)])
	temp['bin_total'] = temp['bin'].map(temp['bin'].value_counts())
	temp['bin_reads'] = temp.bin.map(bin_dict)
	# groupby the bin and get the median of number of umis
	temp = temp[['bin_total', 'bin_reads', 'umis']].groupby(['bin_total', 'bin_reads']).median()
	temp.reset_index(inplace=True)
	temp.rename({'umis': 'median_umis'}, axis=1, inplace=True)
	temp.sort_values(by='bin_reads', inplace=True)
	# plot de plot
	ax = sns.lineplot(x='bin_reads', y='median_umis', marker='o', data=temp)
	ax.set_ylabel('Median UMIs per Cell')
	ax.set_xlabel('Sequencing Reads per Cell')
	ax.set_title(kind)
	plt.draw()
	plt.tight_layout()
	if kind == 'Pre-correction':
		kind = 'pre_correction'
	elif kind == 'Post-correction':
		kind = 'post_correction'
	fname = '{}_{}_sequencing_saturation.png'.format(oprefix, kind)
	plt.savefig(fname)
	plt.clf()
