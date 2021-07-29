import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from pandarallel import pandarallel
from Bio import pairwise2
import time 
import math
import sys
import argparse
import os

def get_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('-f', dest='fastq',
		help='FASTQ file output from Lima with LR-Split-seq reads.')
	parser.add_argument('-steps', dest='steps', default='all',
		help='Comma separated list of steps to perform. Default is all. Options '\
			'include: score_linkers, align_linkers, correct_bcs, trim, i_filt, '\
            'write_fastq')
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser.add_argument('-t', dest='threads', 
		help='Number of threads to run on (multithreading is recommended)')
	parser.add_argument('-i_file', dest='i_file', default=None,
		help='Barcodes from Illumina to filter PacBio barcodes on. '+\
			 'Default: None')
	parser.add_argument('-rc', dest='rc', default=500,
		help='Number of reads/bc to require for filtering')
	# parser.add_argument('-v', dest='verbose', default=False,
	# 	help='Display ')

	args = parser.parse_args()
	return args

def check_steps(steps):
	steps = steps.split(',')
	ctrl = ['score_linkers', 'align_linkers', 'correct_bcs', \
		'trim', 'i_filt', 'write_fastq']
	if steps == ['all']:
		steps = ctrl 
		
	# order based on ctrl order
	ord_steps = [step for step in ctrl if step in steps]
	invalid_steps = [step for step in steps if step not in ctrl]

	# step doesn't exist
	if invalid_steps:
		raise ValueError('Steps {} are not valid choices.'.format(invalid_steps))

	# if a number of steps are given are they all in a row?
	if len(ord_steps) >= 1:
		gaps = [1 if step in ord_steps else 0 for step in ctrl]
		gaps = ''.join([str(x) for x in gaps])
		if '01' in gaps:
			raise ValueError('Multiple steps given without all intermediate steps.')

	# set up dictionary with steps
	s = dict()
	if steps == 'all':
		for step in ctrl:
			s[step] = True
	else:
		for step in ord_steps:
			s[step] = True
	return s

def get_linkers():
	# TODO
	# might need to reverse rc and norm
	# linker sequence between barcodes 1 and 2
	l1_rc = 'CCACAGTCTCAAGCACGTGGAT'
	l1 = rev_comp(l1_rc)
	# linker sequence between barcodes 2 and 3
	l2_rc = 'AGTCGTACGCCGATGCGAAACATCGGCCAC'
	l2 = rev_comp(l2_rc)
	return l1, l1_rc, l2, l2_rc

# From the Parse biosciences pipeline
def load_barcodes():
	pkg_path = os.path.dirname(__file__)
	with open(pkg_path + '/barcodes/bc_dict_v1.pkl', 'rb') as f:
		edit_dict_v1 = pickle.load(f)
	with open(pkg_path + '/barcodes/bc_dict_v2.pkl', 'rb') as f:
		edit_dict_v2 = pickle.load(f)
	
	bc1_edit_dict = edit_dict_v1
	bc2_edit_dict = edit_dict_v1
	bc3_edit_dict = edit_dict_v2
	
	return bc1_edit_dict, bc2_edit_dict, bc3_edit_dict

# From the Parse biosciences pipeline
def load_barcodes_set():
	pkg_path = os.path.dirname(__file__)
	with open(pkg_path + '/barcodes/bc_dict_v1.pkl', 'rb') as f:
		edit_dict_v1 = pickle.load(f)
	with open(pkg_path + '/barcodes/bc_dict_v2.pkl', 'rb') as f:
		edit_dict_v2 = pickle.load(f)
		
	# read in barcode sequences
	bc_8nt_v1 = pd.read_csv(pkg_path + '/barcodes/bc_8nt_v1.csv',names=['barcode'],index_col=0).barcode.values
	bc_8nt_v2 = pd.read_csv(pkg_path + '/barcodes/bc_8nt_v2.csv',names=['barcode'],index_col=0).barcode.values
	bc1_edit_dict = edit_dict_v1
	bc2_edit_dict = edit_dict_v1
	bc3_edit_dict = edit_dict_v2
	bc_8nt_set_dict = {}
	bc_8nt_set_dict['bc1'] = set(bc_8nt_v1)
	bc_8nt_set_dict['bc2'] = set(bc_8nt_v1)
	bc_8nt_set_dict['bc3'] = set(bc_8nt_v2)
	bc_8nt_bc1 = bc_8nt_set_dict['bc1']
	bc_8nt_bc2 = bc_8nt_set_dict['bc2']
	bc_8nt_bc3 = bc_8nt_set_dict['bc3']
	return list(bc_8nt_bc1), list(bc_8nt_bc2), list(bc_8nt_bc3)

# From the Parse biosciences pipeline
def get_bc1_matches():
	# from spclass.py - barcodes and their well/primer type identity
	pkg_path = os.path.dirname(__file__)
	bc_file = pkg_path+'/barcodes/bc_8nt_v2.csv'
	bc_df = pd.read_csv(bc_file, index_col=0, names=['bc'])
	bc_df['well'] = [i for i in range(0, 48)]+[i for i in range(0, 48)]
	bc_df['primer_type'] = ['dt' for i in range(0, 48)]+['randhex' for i in range(0, 48)]

	# pivot on well to get df that matches bcs with one another from the same well
	bc_df = bc_df.pivot(index='well', columns='primer_type', values='bc')
	bc_df = bc_df.rename_axis(None, axis=1).reset_index()
	bc_df.rename({'dt': 'bc1_dt', 'randhex': 'bc1_randhex'}, axis=1, inplace=True)

	return bc_df

def rev_comp(s):
	rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 
		   'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 
		   'n': 'n', '*': '*'}
	rev_comp = [rc_map[i] for i in s]
	rev_comp.reverse()
	rev_comp = ''.join(rev_comp)
	return rev_comp

def align_linker_score(x, l_seqs, l_prefs):
	
	entry = {}
	
	# +1 for match
	# -1 for mismatch
	# -1 for gap open
	# -1 for gap extend
	# this scoring schema allows us to know exactly 
	# how many errors are in each best alignment
	for l_seq, pref in zip(l_seqs, l_prefs):
		a = pairwise2.align.localms(x.seq, l_seq,
				1, -1, -1, -1,
				one_alignment_only=True, 
				score_only=True)
		score = a
		entry['{}_score'.format(pref)] = score
	return entry

def find_linkers(df, t=1):

	l1, l1_rc, l2, l2_rc = get_linkers()
	l_seqs = [l1, l1_rc, l2, l2_rc]
	l_prefs = ['l1', 'l1_rc', 'l2', 'l2_rc']
	
	if t == 1:
		l_df = df.apply(align_linker_score, args=(l_seqs, l_prefs),
			axis=1, result_type='expand')
	else:
		pandarallel.initialize(nb_workers=t)	
		l_df = df.parallel_apply(align_linker_score, args=(l_seqs, l_prefs),
			axis=1, result_type='expand')
	
	df = pd.concat([df, l_df], axis=1)
	
	return df  


def fastq_to_df(fastq):
	# TODO remove '@'s from read names first next time!
	# get the sequences from each read
	seqs = []
	read_names = []
	strands = []
	i = 0
	with open(fastq, 'r') as f:
		while True:
			read_name = f.readline().strip()
			read_name = read_name[1:]
			if len(read_name)==0:
				break
			seq = f.readline().strip()
			strand = f.readline().strip()
			qual = f.readline()
			seqs.append(seq)
			strands.append(strand)
			read_names.append(read_name)
			i += 1
			if i % 1000000==0:
				print('Processed {} reads'.format(i))

	# pack everything into a dataframe
	df = pd.DataFrame(seqs)
	df.columns = ['seq']
	df['read_name'] = read_names
	df['strand'] = strands
	
	return df

def plot_hist(x, **kwargs):
	
	# linker 1
	if x.max() == 22:
		mismatch_lines = [22, 21, 20, 19]
	# linker 2
	elif x.max() == 30:
		mismatch_lines = [30, 29, 28, 27]
	ax = sns.histplot(x, binwidth=1)
	for l in mismatch_lines: 
		plt.axvline(l, color='k', linestyle='-', linewidth=1)


def plot_linker_scores(df, oprefix):
	val_vars = ['l1_score', 'l2_score', 'l1_rc_score', 'l2_rc_score']
	cols = ['read_name'] + val_vars
	temp = df[cols].melt(id_vars='read_name', value_vars=val_vars)
	
	g = sns.FacetGrid(temp, col='variable')
	g.map(plot_hist, 'value')

	fname = oprefix+'_linker_score_dists.png'
	plt.savefig(fname)

	plt.clf()

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
	plt.savefig(fname)

	plt.clf()

# get the linker alignments with the seq and the start and end of each alignment
def align_linker_seqs(x, l1, l2, l1_rc, l2_rc):
	
	entry = {}
	
	# forward or reverse search?
	if x.l_dir == '-':
		l1 = l1_rc
		l2 = l2_rc
		
	# compute alignments for both linkers in the correct orientation
	# use the default alignment
	l1_a = pairwise2.align.localms(x.seq, l1, 
				1, -1, -1, -1,
				one_alignment_only=True)
	l2_a = pairwise2.align.localms(x.seq, l2,
				1, -1, -1, -1,
				one_alignment_only=True)

	# grab start and end of each linker
	try:
		l1_start = l1_a[0].start
	except:
		print('REEEEE')
		print(x.read_name)
		raise ValueError('REEEE')
	l1_stop = l1_a[0].end
	l2_start = l2_a[0].start
	l2_stop = l2_a[0].end
	
	# calculate some metrics
	l1_len = l1_stop-l1_start
	l2_len = l2_stop-l2_start
	if x.l_dir == '+':
		bc2_len = l1_start-l2_stop
	elif x.l_dir == '-':
		bc2_len = l2_start-l1_stop
	else: 
		bc2_len = np.nan
		

	# construct an entry
	entry['l1_start'] = l1_start
	entry['l1_stop'] = l1_stop
	entry['l2_start'] = l2_start
	entry['l2_stop'] = l2_stop
	entry['l1_len'] = l1_len
	entry['l2_len'] = l2_len
	entry['bc2_len'] = bc2_len
	 
	return entry

# df: dataframe of each read with its linker sequences scored 
# choose l#_m args or l#_p args
# l1_m (optional): number of mismatches permitted in linker 1
# l2_m (optional): number of mismatches permitted in linker 2
# l1_p (optional): percentage of mismatches permitted in linker 1
# l2_p (optional): percentage of mismatches permitted in linker 2
def get_linker_alignments(df, t=1, l1_m=None, l2_m=None, l1_p=None, l2_p=None,
		keep_dupes=False, verbose=False):
	
	# calculate lowest scores viable for l1 and l2
	if l1_m and l2_m:
		if float(l1_m).is_integer() and float(l2_m).is_integer():
			l1_min = 22-l1_m 
			l2_min = 30-l2_m
		else:
			raise TypeError('l1_m and l2_m must be integers.')
	elif l1_p and l2_p:
		if isnumeric(l1_p) and isnumeric(l2_p):
			l1_min = 22-math.ceil((22/100)*i)
			l2_min = 30-math.ceil((30/100)*j)
		else:
			raise TypeError('l1_p and l2_p must be numeric.')
	else:
		raise Exception('Please provide input for l1_m and l2_m '
			'or l1_p and l2_p')

	# init some dfs
	fwd_df = pd.DataFrame()
	rev_df = pd.DataFrame()
	both_df = pd.DataFrame()
	fwd_tie_df = pd.DataFrame()
	rev_tie_df = pd.DataFrame()
#	 cand_dfs = [fwd_df, rev_df, both_df, fwd_tie_df, rev_tie_df]
#	 for d in cand_dfs:
#		 d = pd.DataFrame()

	# find viable reads and format into a new dataframe
	fwd_df = df.loc[(df.l1_score>=l1_min)&(df.l2_score>=l2_min)]
	rev_df = df.loc[(df.l1_rc_score>=l1_min)&(df.l2_rc_score>=l2_min)]
	both_df = df.loc[list(set(fwd_df.index)&set(rev_df.index))]
	
	if verbose:
		print('Found {} reads with linkers found in the fwd direction'.format(len(fwd_df.index)))
		print('Found {} reads with linkers found in the rev direction'.format(len(rev_df.index)))
		print('Found {} reads with linkers found in both directions'.format(len(both_df.index)))	
	
	# first remove those that are duplicated across fwd and rev dfs
	fwd_df = fwd_df.loc[~fwd_df.index.isin(both_df.index)]
	rev_df = rev_df.loc[~rev_df.index.isin(both_df.index)]
	
	if verbose:
		print('Number of fwd reads after removing reads from both dirs: {}'.format(len(fwd_df.index)))
		print('Number of rev reads after removing reads from both dirs: {}'.format(len(rev_df.index)))
	
	# then assign orientation of found linkers
	fwd_df['l_dir'] = '+'
	rev_df['l_dir'] = '-'
	
	# tie break those that we can that found forward and reverse linkers
	# with scores over the threshold using the sum of scores
	if keep_dupes == True:

		if not both_df.empty:
			both_df['fwd_score'] = both_df.l1_score + both_df.l2_score
			both_df['rev_score'] = both_df.l1_rc_score + both_df.l2_rc_score

			# find real ties and remove them from the df
			tie_df = both_df.loc[both_df.fwd_score == both_df.rev_score]
			both_df = both_df.loc[~both_df.index.isin(tie_df.index)]

			if verbose:
				print('Found {} reads with the same score in fwd and rev directions'.format(len(tie_df.index)))
				print('Number of reads with both fwd and rev linkers that were tie-broken: {}'.format(len(both_df.index)))

			# choose which linkers to use from the non-ties
			both_df['temp'] = both_df[['fwd_score', 'rev_score']].idxmax(axis=1)
			dir_map = {'fwd_score':'+','rev_score':'-'}
			both_df['l_dir'] = both_df.temp.map(dir_map)
			both_df.drop(['fwd_score', 'rev_score', 'temp'], axis=1, inplace=True)

			# we'll examine both directions for the true ties
			tie_df.drop(['fwd_score', 'rev_score'], axis=1, inplace=True)
			fwd_tie_df = tie_df.copy(deep=True)
			rev_tie_df = tie_df.copy(deep=True)
			fwd_tie_df['l_dir'] = '+'
			rev_tie_df['l_dir'] = '-'

			if verbose:
				print('Number of true tied reads {}'.format(len(fwd_tie_df.index)))
				print('Number of true tied reads {}'.format(len(rev_tie_df.index)))

	# concatenate all the dfs 
	dfs = []
	cand_dfs = [fwd_df, rev_df, both_df, fwd_tie_df, rev_tie_df]
	for d in cand_dfs:
		if not d.empty:
			dfs.append(d)
			
	df = pd.concat(dfs)

	if verbose:
		print('Aligning linkers for {} reads'.format(len(df.index)))
	
	# find the start and end of each alignment for the valid reads
	l1, l1_rc, l2, l2_rc = get_linkers()
	if t == 1:
		l_df = df.apply(align_linker_seqs, args=(l1, l2, l1_rc, l2_rc),
			axis=1, result_type='expand')
	else:
		pandarallel.initialize(nb_workers=t)
		l_df = df.parallel_apply(align_linker_seqs, args=(l1, l2, l1_rc, l2_rc),
			axis=1, result_type='expand')
	df = pd.concat([df, l_df], axis=1)
	
	return df

# TODO fix which barcode is which
# TODO also use pandas slice instead of this slow shit maybe
def get_seq_bcs_umis(x):
#	 print(x.read_name)
#	 print(x.seq)
#	 print('l1 loc: {}:{}'.format(x.l1_start, x.l1_stop))
#	 print('l2 loc: {}:{}'.format(x.l2_start, x.l2_stop))
#	 bc1 = x.seq[x.l2_start-8:x.l2_start]
#	 bc2 = x.seq[x.l2_stop:x.l2_stop+8]
#	 bc3 = x.seq[x.l1_stop:x.l1_stop+8]
	# TODO 10/22/20 - this checks out (but rev/fwd nomenclature is WRONG)
	# try:
	if x.l_dir == '+':
		bc1 = x.seq[x.l1_stop:x.l1_stop+8]
		bc2 = x.seq[x.l2_stop:x.l2_stop+8]
		bc3 = x.seq[x.l2_start-8:x.l2_start]
		umi = x.seq[x.l2_start-8-10:x.l2_start-8]
	elif x.l_dir == '-':
		bc2 = rev_comp(x.seq[x.l2_start-8:x.l2_start])
		bc3 = rev_comp(x.seq[x.l2_stop:x.l2_stop+8])
		bc1 = rev_comp(x.seq[x.l1_start-8:x.l1_start])
		umi = rev_comp(x.seq[x.l2_stop+8:x.l2_stop+8+10])
	# except Exception:
	# 	print(x.read_name)
	# 	# print(E)
	# 	sys.exit(1)
	# 	return
	# construct an entry
	entry = {}
	entry['bc1'] = bc1
	entry['bc2'] = bc2
	entry['bc3'] = bc3 
	entry['umi'] = umi
	return entry

def get_bcs_umis(df, t=1):

	# remove fwd/rev ones
	df = df.loc[(df.l_dir == '+')|(df.l_dir == '-')]

	# make sure indices we're gonna use are integers
	df.l1_start = df.l1_start.astype('int')
	df.l1_stop = df.l1_stop.astype('int')
	df.l2_start = df.l2_start.astype('int')
	df.l2_stop = df.l2_stop.astype('int')
		
	# naive method, just look 8 bp up or downstream of the various linkers
	if t == 1:
		temp = df.apply(lambda x: get_seq_bcs_umis(x), axis=1, result_type='expand')
	else: 
		pandarallel.initialize(nb_workers=t)   
		temp = df.parallel_apply(get_seq_bcs_umis, axis=1, result_type='expand')

	df = pd.concat([df, temp], axis=1)

	# # loop through each boi 
	# df['bc1'] = [seq[l1_stop:l1_stop+8] if d == '+' \
	# 	else rev_comp(seq[l2_start-8:l2_start]) \
	# 	for seq,l1_stop,l2_start,d in zip(df.seq,df.l1_stop,df.l2_start,df.l_dir)]
	# df['bc2'] = [seq[l2_stop:l2_stop+8] if d == '+' \
	# 	else rev_comp(seq[l2_start-8:l2_start]) \
	# 	for seq,l2_stop,l2_start,d in zip(df.seq,df.l2_stop,df.l2_start,df.l_dir)]
	# df['bc3'] = [seq[l2_start-8:l2_start] if d == '+' \
	# 	else rev_comp(seq[l2_stop:l2_stop+8]) \
	# 	for seq,l2_start,l2_stop,d in zip(df.seq,df.l2_start,df.l2_stop,df.l_dir)]
	# df['umi'] = [seq[l2_start-8-10:l2_start-8] if d == '+' \
	# 	else rev_comp(seq[l2_stop+8:l2_stop+8+10]) \
	# 	for seq,l2_start,l2_stop,d in zip(df.seq,df.l2_start,df.l2_stop,df.l_dir)]
	
	return df

# From the Parse Biosciences pipeline
def get_min_edit_dists(bc,edit_dict,max_d=3):
	"""Returns a list of nearest edit dist seqs
	Input 8nt barcode, edit_dist_dictionary
	Output <list of nearest edit distance seqs>, <edit dist>"""
	bc_matches = edit_dict[0][bc]
	edit_dist = 0
	if (len(bc_matches)==0) and (max_d>=1):
		edit_dist+=1
		bc_matches = edit_dict[1][bc]
	if (len(bc_matches)==0) and (max_d>=2):
		edit_dist+=1
		bc_matches = edit_dict[2][bc]
	if (len(bc_matches)==0) and (max_d>=3):
		edit_dist+=1
		bc_matches = edit_dict[3][bc]
	return bc_matches,edit_dist

# From the Parse Biosciences pipeline
def get_perfect_bc_counts(df, verbose=False):
	reads_in_cells_thresh = 0.92
	
	# TODO
	bc_8nt_bc3, bc_8nt_bc2, bc_8nt_bc1 = load_barcodes_set()

	df['bc1_valid'] = df.bc1.isin(list(bc_8nt_bc1))
	df['bc2_valid'] = df.bc2.isin(list(bc_8nt_bc2))
	df['bc3_valid'] = df.bc3.isin(list(bc_8nt_bc3))

	n_bc1_valid = len(df.loc[df.bc1_valid == True].index)
	n_bc2_valid = len(df.loc[df.bc2_valid == True].index)
	n_bc3_valid = len(df.loc[df.bc3_valid == True].index)
	n_all_valid = len(df.loc[(df[['bc1_valid', 'bc2_valid', 'bc3_valid']].all(axis=1))])

	if verbose:
		print('Valid bc1 counts: {}'.format(n_bc1_valid))
		print('Valid bc2 counts: {}'.format(n_bc2_valid))
		print('Valid bc3 counts: {}'.format(n_bc3_valid))
		print('All 3 valid barcode counts: {}'.format(n_all_valid))

	# calculate thresholds based on what we see in the data
	# my way
	counts = df.loc[(df.bc1_valid)&(df.bc2_valid)&(df.bc3_valid)]
	counts = counts[['bc1', 'bc2', 'bc3']].value_counts()
	count_threshold = max(2, counts.iloc[abs(counts.cumsum()/counts.sum()-reads_in_cells_thresh).values.argmin()])
	
	return df, counts, count_threshold

# From the Parse Biosciences pipeline
def correct_seq_barcodes(x,
		counts, count_thresh,
		bc_edit_dist,
		bc1_dict,
		bc2_dict,
		bc3_dict):
	
	bc1 = x.bc1
	bc2 = x.bc2
	bc3 = x.bc3
	
	debug = False
	if debug:
		print('Attempting to correct barcodes : {}, {}, {}'.format(bc1, bc2, bc3))
	
	bc1_matches,edit_dist1 = get_min_edit_dists(bc1,bc1_dict,max_d=bc_edit_dist)
	bc2_matches,edit_dist2  = get_min_edit_dists(bc2,bc2_dict,max_d=bc_edit_dist)
	bc3_matches,edit_dist3  = get_min_edit_dists(bc3,bc3_dict,max_d=bc_edit_dist)
	
	# Check if any barcode matches have a counts above the threshold
	if 0 == edit_dist1 == edit_dist2 == edit_dist3:
		bc1 = bc1_matches[0]
		bc2 = bc2_matches[0]
		bc3 = bc3_matches[0]
	else:
		matches = 0
		for bc1_m in bc1_matches:
			for bc2_m in bc2_matches:
				for bc3_m in bc3_matches:

					if debug:
						print(bc1_m)
						print(bc2_m)
						print(bc3_m)
					try:
						cur_counts = counts.loc[(bc1_m,bc2_m,bc3_m)]
					except:
						cur_counts = 0
					if cur_counts>count_thresh:
						bc1_fixed = bc1_m
						bc2_fixed = bc2_m
						bc3_fixed = bc3_m
						matches += 1
		if matches==1:
			bc1 = bc1_fixed
			bc2 = bc2_fixed
			bc3 = bc3_fixed
		else:
			bc1 = bc2 = bc3 = np.nan
		
	entry = {}
	entry['bc1'] = bc1
	entry['bc2'] = bc2
	entry['bc3'] = bc3

	if debug:
		print()
		
	return entry 

def correct_barcodes(df, counts, count_thresh, bc_edit_dist, t=1):
	
	bc3_dict, bc2_dict, bc1_dict = load_barcodes()
	
	if t == 1:
		temp = df.parallel_apply(correct_seq_barcodes,
				args=(counts, count_thresh, bc_edit_dist, bc1_dict, bc2_dict, bc3_dict),
				axis=1, result_type='expand')
	else: 
		temp = df.apply(lambda x: correct_seq_barcodes(x, counts, count_thresh,
					bc_edit_dist, bc1_dict, bc2_dict, bc3_dict),
					axis=1, result_type='expand')
		pandarallel.initialize(nb_workers=t)
	
	df.dropna(axis=0, subset=['bc1', 'bc2', 'bc3'], inplace=True)
	df.drop(['bc1', 'bc2', 'bc3'], axis=1, inplace=True)
	df = pd.concat([df, temp], axis=1)
	
	return df

def plot_umis_v_barcodes(df, oprefix, kind):	
	bc_cols = ['bc1', 'bc2', 'bc3']

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
	ax.set_xlabel('Ranked cells by # UMIs (logscale)')
	ax.set_ylabel('# UMIs (logscale)')
	ax.set_title(kind)

	if kind == 'Pre-correction':
		kind = 'pre_correction'
	elif kind == 'Post-correction':
		kind = 'post_correction'
	elif kind == 'Illumina':
		kinda = 'illumina'

	plt.tight_layout()

	fname = '{}_{}_umis_v_barcodes.png'.format(oprefix, kind)
	plt.savefig(fname)
	plt.clf()

def plot_umis_per_cell(df, oprefix, kind):	   
	bc_cols = ['bc1', 'bc2', 'bc3']

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
	fname = '{}_{}_median_umis_v_reads.png'.format(oprefix, kind)
	plt.savefig(fname)

	plt.clf()

def trim_bcs_x(x):
#	 if x.dir == '+':
#		 trim_start = x.l1_stop+8
#	 elif x.dir == '-':
#		 trim_start = x.l1_start-8
	if x.l_dir == '+':
		trim_seq = x.seq[x.trim_start:]
	elif x.l_dir == '-':
		trim_seq = x.seq[:x.trim_start]
	return trim_seq  

# trim off the barcodes from the start of bc1 until the end of the read
def trim_bcs(df, t=1, verbose=False):
	
	# add placeholders
	df.trim_start = np.nan
	
	# separate into fwd/rev 
	fwd = df.loc[df.l_dir == '+'].copy(deep=True)
	rev = df.loc[df.l_dir == '-'].copy(deep=True)
		
	# where should we start trimming the sequence?
	# TODO remember that these directions are wrong
	fwd['trim_start'] = fwd.l1_stop+8
	rev['trim_start'] = rev.l1_start-8
	fwd.trim_start = fwd.trim_start.astype('int')
	rev.trim_start = rev.trim_start.astype('int')
	
	# put the dfs back together
	dfs = []
	cand_dfs = [fwd, rev]
	for d in cand_dfs:
		if not d.empty:
			dfs.append(d)	
	df = pd.concat(dfs)
	
	if t == 1:
		df['trim_seq'] = df.apply(trim_bcs_x, axis=1)
	else:
		pandarallel.initialize(nb_workers=t)
		df['trim_seq'] =  df.parallel_apply(trim_bcs_x, axis=1)	

	# remove sequences that are empty now because of weird linker things
	df['trim_seq_len'] = df.apply(lambda x: len(x.trim_seq), axis=1)
	df = df.loc[df.trim_seq_len != 0]

	if verbose:
		n = len(df.index)
		print('Number of reads after removing truncated sequences: {}'.format(n))
		
	return df

def flip_reads_x(x):
	# TODO remember that this is the wrong direction (need to verify 11/12/20)
	if x.l_dir == '+':
		s = rev_comp(x.trim_seq)
	else:
		s = x.trim_seq
	return s
	
def flip_reads(df, t=1):
	if t == 1:
		df['trim_seq'] = df.apply(flip_reads_x, axis=1)
	else:
		pandarallel.initialize(nb_workers=t)
		df['trim_seq'] = df.parallel_apply(flip_reads_x, axis=1)

	df.seq = df.trim_seq
	df.drop('trim_seq', axis=1, inplace=True)
	
	return df

# plot read lengths after removing the barcode construct
def plot_read_length(df, oprefix):
	ax = sns.displot(data=df, x='trim_seq_len', color='#CF91A3')
	ax.set(xlabel='Read length', ylabel='Number of reads')

	fname = '{}_read_length_dist.png'.format(oprefix)
	plt.savefig(fname)
	plt.clf()

# remove bcs that aren't in the corresponding Illumina set of barcodes
def filter_on_illumina(df, i_df):

    # get polydT and randhex barcodes
    dt_bcs = i_df.bc3+i_df.bc2+i_df.bc1_dt
    randhex_bcs = i_df.bc3+i_df.bc2+i_df.bc1_randhex
    i_bcs = dt_bcs.tolist()+randhex_bcs.tolist()

    # subset long-read barcodes on those that are present in illumina data
    df = df[['read_name', 'seq', 'bc1', 'bc2', 'bc3', 'umi']]
    df['bc'] = df.bc3+df.bc2+df.bc1

    df = df.loc[df.bc.isin(i_bcs)]

    return df

# remove all combinations of barcodes that don't have enough reads
# if we're also filtering on Illumina bcs, take all cells that pass the 
# threshold in either dt or randhex
def filter_on_read_count(df, read_thresh):
    
    bc_df = get_bc1_matches()
    
    df['bc'] = df.bc3+df.bc2+df.bc1

    # if we have a number of reads threshold, filter on that
    temp = df.copy(deep=True)
    temp = temp.value_counts(['bc', 'bc1', 'bc2', 'bc3']).reset_index(name='bc_counts')
    # temp = temp.loc[temp.bc_counts>read_thresh]
    temp['bc1_partner'] = temp.apply(lambda x: bc_df.loc[bc_df.bc1_dt == x.bc1, 'bc1_randhex'].values[0] \
        if x.bc1 in bc_df.bc1_dt.tolist() else bc_df.loc[bc_df.bc1_randhex == x.bc1, 'bc1_dt'].values[0], \
        axis=1)        
    temp['bc_partner'] = temp.bc3+temp.bc2+temp.bc1_partner
    temp['bc_partner_counts'] = temp.apply(lambda x: temp.loc[temp.bc == x.bc_partner, 'bc_counts'].values[0] \
        if x.bc_partner in temp.bc.tolist() else 0, axis=1)
    temp['total_counts'] = temp.bc_counts + temp.bc_partner_counts
    temp = temp.loc[temp.total_counts>read_thresh]
    valid_bcs = temp.bc.tolist()+temp.bc_partner.tolist()
    
    df = df.loc[df.bc.isin(valid_bcs)]

    return df

# process illumina barcodes to add the randhex bc1s to 
# the list of valid barcodes as well
# Returns a df with barcode combinations that are possible
def process_illumina_bcs(ifile):

	# read in illumina barcodes
	i_df = pd.read_csv(ifile, header=None)
	i_df.columns = ['bc']
	i_df['bc3'] = i_df.bc.str.slice(start=0, stop=8)
	i_df['bc2'] = i_df.bc.str.slice(start=8, stop=16)
	i_df['bc1'] = i_df.bc.str.slice(start=16, stop=24)
	
	bc_df = get_bc1_matches()

	# then merge on dt bc1 with illumina barcodes
	i_df = i_df.merge(bc_df, how='left', left_on='bc1', right_on='bc1_dt')
	
	return i_df


# format/write fastq
def write_fastq(df, oprefix):

	# create the read name header with the bc and umi information
	df.fillna(value='_', inplace=True)
	df['header'] = '@'+df.read_name+':'+df.bc+'_'+df.umi
	df = df[['header', 'seq']]

	# write the fastq
	fname = oprefix+'.fastq'
	ofile = open(fname, 'w')
	for ind, entry in df.iterrows():
		try:
			ofile.write(entry.header+'\n')
		except:
			print(entry.header)
		ofile.write(entry.seq+'\n')
		ofile.write('+\n')
		ofile.write(''.join(['5' for i in range(len(entry.seq))])+'\n')

	ofile.close()

def main():

	args = get_args()
	fastq = args.fastq
	oprefix = args.oprefix
	t = int(args.threads)
	i_file = args.i_file
	rc = int(args.rc)
	steps = check_steps(args.steps)

	# verbose = args.verbose

	sns.set_context("paper", font_scale=1.8)

	# score linkers in each read
	if steps['score_linkers']:
		df = fastq_to_df(fastq)	
		df = find_linkers(df, t=t)

		fname = oprefix+'_seq_linker_alignment_scores.tsv'
		df.to_csv(fname, sep='\t', index=False)

	    # some qc plots
		plot_linker_scores(df, oprefix)
		plot_linker_heatmap(df, oprefix, how='integer')
		plot_linker_heatmap(df, oprefix, how='proportion')


	# TODO - look for dist of occurrences of each linker within each read
	# we can maybe decrease the search space for each read by truncating 
	# the read
	# align linkers that scored high enough to get position of 
	# each linker within each read - this step can take a while!
	if steps['align_linkers'] and not steps['score_linkers']:
		fname = oprefix+'_seq_linker_alignment_scores.tsv'
		df = pd.read_csv(fname, sep='\t')
	if steps['align_linkers']: 
		df = get_linker_alignments(df, t=t, l1_m=3, l2_m=3, verbose=True)
		fname = oprefix+'_seq_linker_alignments.tsv'
		df.to_csv(fname, sep='\t', index=False)

	# correct barcodes 
	if steps['correct_bcs'] and not steps['align_linkers']:
		fname = oprefix+'_seq_linker_alignment_scores.tsv'
		df = pd.read_csv(fname, sep='\t')
	if steps['correct_bcs']:
		# get barcode information from each read
		df = get_bcs_umis(df, t=t)
		df, counts, count_thresh = get_perfect_bc_counts(df)

		fname = oprefix+'_seq_bcs.tsv'
		df.to_csv(fname, sep='\t', index=False)

		# some more qc plots
		plot_umis_v_barcodes(df, oprefix, 'Pre-correction')
		plot_umis_per_cell(df, oprefix, 'Pre-correction')

		edit_dist = 3
		df = correct_barcodes(df, counts, count_thresh, edit_dist, t=t)
		fname = oprefix+'_seq_corrected_bcs.tsv'

		# ***TODO probably want to drop nans here.... not sure what's going on***
		df.to_csv(fname, sep='\t', index=False)

		# some more qc plots
		plot_umis_v_barcodes(df, oprefix, 'Post-correction')
		plot_umis_per_cell(df, oprefix, 'Post-correction')

	# trim reads of their linker + barcode construct
	if steps['trim'] and not steps['correct_bcs']:
		fname = oprefix+'_seq_corrected_bcs.tsv'
		df = pd.read_csv(fname, sep='\t')
	if steps['trim']:
		df = trim_bcs(df, t=t)
		df = flip_reads(df, t=t)

		fname = oprefix+'_trimmed_flipped.tsv'
		df.to_csv(fname, sep='\t', index=False)

		# what do the read lengths look like after this?
		plot_read_length(df, oprefix)

	# filter based on which bc combinations were also seen in Illumina
	if steps['i_filt'] and not steps['correct_bcs']:
		fname = oprefix+'_trimmed_flipped.tsv'
		df = pd.read_csv(fname, sep='\t')
	if steps['i_filt']:
		if i_file:
			i_bcs = process_illumina_bcs(i_file)
			df = filter_on_illumina(df, i_bcs)
			plot_umis_v_barcodes(df, oprefix, 'Illumina')
		# df = filter_on_read_count(df, rc)
		fname = oprefix+'_filtered.tsv'
		df.to_csv(fname, sep='\t', index=False)

	# write the fastq with the cell ID and UMI in the read name
	if steps['write_fastq'] and not steps['i_filt']:
		fname = oprefix+'_filtered.tsv'
		df = pd.read_csv(fname, sep='\t')
	if steps['write_fastq']:
		df = write_fastq(df, oprefix)

if __name__ == '__main__':
	main()
