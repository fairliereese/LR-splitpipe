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

###################################################################################
############################### Helper functions ##################################
###################################################################################

	# https://stackoverflow.com/questions/25962114/how-do-i-read-a-large-csv-file-with-pandas

def load_barcodes():
	"""
	Load the barcodes. Adapted from the Parse biosciences pipeline.

	Returns:
		bc#_edit_dict (dict): Dict for barcode<1,2,3> with
			key: query bc
			item: corrected bc
	"""
	pkg_path = os.path.dirname(__file__)
	pkg_path = '/'.join(pkg_path.split('/')[:-1])
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
	"""
	Load the barcodes. Adapted from the Parse biosciences pipeline.
	"""
	pkg_path = os.path.dirname(__file__)
	pkg_path = '/'.join(pkg_path.split('/')[:-1])
	# pkg_path = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/211206_lr_splitpipe_speedup/'

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

def rev_comp(s):
	"""
	Compute the reverse complement of a given sequence

	Parameters:
		s (str): Input sequence
	Returns:
		rev_comp (str): Reverse complement of input sequence
	"""
	rc_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
			'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
			'n': 'n', '*': '*'}
	rev_comp = [rc_map[i] for i in s]
	rev_comp.reverse()
	rev_comp = ''.join(rev_comp)
	return rev_comp

def get_linkers():
	"""
	Get the Parse biosciences linker sequences
	"""
	# linker sequence between barcodes 1 and 2
	l1_rc = 'CCACAGTCTCAAGCACGTGGAT'
	l1 = rev_comp(l1_rc)
	# linker sequence between barcodes 2 and 3
	l2_rc = 'AGTCGTACGCCGATGCGAAACATCGGCCAC'
	l2 = rev_comp(l2_rc)
	return l1, l1_rc, l2, l2_rc

def get_min_edit_dists(bc, edit_dict, max_d=3):
	"""
	Returns a list of nearest edit dist seqs.
	Adapted from Parse biosciences

	Parameters:
		bc (str): 8nt barcode
		edit_dict (dict): Dict of bcs within
			edit distance of bc
		max_d (int): Edit distance

	Returns:
		bc_matches
		edit_dist
	"""
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
	return bc_matches, edit_dist

###################################################################################
############################### Processing steps ##################################
###################################################################################

def fastq_to_df(fname, oprefix, verbose=1):
	"""
	Save the input fastq file into a table format into a file called
	oprefix_table.tsv.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress

	Returns:
		ofile (str): Name of output file
	"""
	# get the sequences from each read
	seqs = []
	read_names = []
	strands = []
	i = 1
	ofile = '{}_table.tsv'.format(oprefix)
	with open(fname, 'r') as f:
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

			# print out notification and write to file
			chunksize = 100000
			if i % chunksize == 0 and i != 1 :
				if verbose == 2:
					print('Processed {} reads'.format(i))

				# pack everything into a dataframe
				df = pd.DataFrame(seqs)
				df.columns = ['seq']
				df['read_name'] = read_names
				df['strand'] = strands

				# first write
				if i == chunksize:
					df.to_csv(ofile, sep='\t', index=False)
					read_names = []
					strands = []
					seqs = []
				else:
					df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')
					read_names = []
					strands = []
					seqs = []
			i += 1

	# cleanup
	if len(seqs) > 0:
		df = pd.DataFrame(seqs)
		df.columns = ['seq']
		df['read_name'] = read_names
		df['strand'] = strands
		df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

	return ofile

def score_linkers(fname, oprefix, t=1,
				 verbose=1, chunksize=10**5,
				 delete_input=False):
	"""
	Find and report highest linker scores in each read.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	l1, l1_rc, l2, l2_rc = get_linkers()
	l_seqs = [l1, l1_rc, l2, l2_rc]
	l_prefs = ['l1', 'l1_rc', 'l2', 'l2_rc']

	# loop through chunks of the file
	i = 0
	ofile = '{}_seq_linker_alignment_scores.tsv'.format(oprefix)

	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):
		if t == 1:
			l_df = df.apply(score_linkers_x, args=(l_seqs, l_prefs),
				axis=1, result_type='expand')
		else:
			pandarallel.initialize(nb_workers=t, verbose=1)
			l_df = df.parallel_apply(score_linkers_x, args=(l_seqs, l_prefs),
				axis=1, result_type='expand')
		df = pd.concat([df, l_df], axis=1)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize
		if verbose == 2:
			print('Found linker scores for {} reads'.format(i))

	# delete input file
	if delete_input:
		os.remove(fname)

	return ofile

def align_linkers(fname, oprefix,
						  t=1, chunksize=10**5,
						  l1_m=None, l2_m=None,
						  l1_p=None, l2_p=None,
						  keep_dupes=False, verbose=1,
						  delete_input=False):
	"""
	Find indices of highest-scoring linker in each read.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		t (int): Number of threads to run on
		l1_m (int): Number of allowable mismatches in linker 1
		l2_m (int): Number of allowable mismatches in linker 2
		l1_p (float): Proportion of allowable mismatches in linker 1
		l2_p (float): Proportion of allowable mismatches in linker 2
		keep_dupes (bool):
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

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

	# loop through chunks of the file
	i = 0
	n_alignments = 0
	ofile = '{}_seq_linker_alignments.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		# init some dfs
		fwd_df = pd.DataFrame()
		rev_df = pd.DataFrame()
		both_df = pd.DataFrame()
		fwd_tie_df = pd.DataFrame()
		rev_tie_df = pd.DataFrame()

		# find viable reads and format into a new dataframe
		fwd_df = df.loc[(df.l1_score>=l1_min)&(df.l2_score>=l2_min)]
		rev_df = df.loc[(df.l1_rc_score>=l1_min)&(df.l2_rc_score>=l2_min)]
		both_df = df.loc[list(set(fwd_df.index)&set(rev_df.index))]

		if verbose == 2:
			print('Found {} reads with linkers found in the fwd direction'.format(len(fwd_df.index)))
			print('Found {} reads with linkers found in the rev direction'.format(len(rev_df.index)))
			print('Found {} reads with linkers found in both directions'.format(len(both_df.index)))

		# first remove those that are duplicated across fwd and rev dfs
		fwd_df = fwd_df.loc[~fwd_df.index.isin(both_df.index)]
		rev_df = rev_df.loc[~rev_df.index.isin(both_df.index)]

		if verbose == 2:
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

				if verbose == 2:
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

				if verbose == 2:
					print('Number of true tied reads {}'.format(len(fwd_tie_df.index)))
					print('Number of true tied reads {}'.format(len(rev_tie_df.index)))

		# concatenate all the dfs
		dfs = []
		cand_dfs = [fwd_df, rev_df, both_df, fwd_tie_df, rev_tie_df]
		for d in cand_dfs:
			if not d.empty:
				dfs.append(d)

		df = pd.concat(dfs)

		n_alignments += len(df.index)

		# find the start and end of each alignment for the valid reads
		l1, l1_rc, l2, l2_rc = get_linkers()
		if t == 1:
			l_df = df.apply(align_linkers_x, args=(l1, l2, l1_rc, l2_rc),
				axis=1, result_type='expand')
		else:
			pandarallel.initialize(nb_workers=t, verbose=1)
			l_df = df.parallel_apply(align_linkers_x, args=(l1, l2, l1_rc, l2_rc),
				axis=1, result_type='expand')
		df = pd.concat([df, l_df], axis=1)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize
		if verbose == 2:
			print('Found linkers for {} reads'.format(i))

	if verbose < 0:
		print('Found {} reads with valid linker combinations'.format(n_alignments))

	# delete input file
	if delete_input:
		os.remove(fname)

	return ofile

def get_bcs_umis(fname, oprefix, t=1,
				 verbose=1, chunksize=10**5,
				 delete_input=False):
	"""
	Extract the barcodes and UMIs based on linker location
	in each read

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		t (int): Number of threads to run on
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	i = 0
	ofile = '{}_bcs.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		# remove fwd/rev ones
		df = df.loc[(df.l_dir == '+')|(df.l_dir == '-')]

		# make sure indices we're gonna use are integers
		df.l1_start = df.l1_start.astype('int')
		df.l1_stop = df.l1_stop.astype('int')
		df.l2_start = df.l2_start.astype('int')
		df.l2_stop = df.l2_stop.astype('int')

		# naive method, just look 8 bp up or downstream of the various linkers
		if t == 1:
			temp = df.apply(lambda x: get_bcs_umis_x(x), axis=1, result_type='expand')
		else:
			pandarallel.initialize(nb_workers=t, verbose=1)
			temp = df.parallel_apply(get_bcs_umis_x, axis=1, result_type='expand')

		df = pd.concat([df, temp], axis=1)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize
		if verbose == 2:
			print('Found bcs and umis for {} reads'.format(i))

	return ofile

def get_perfect_bc_counts(fnames, verbose=1):
	"""
	Count how many reads with valid barcodes within edit
	distance 3 there are. Adapted from Parse biosciences.

	Parameters:
		fnames (list of str): Files to process
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress

	Returns:
		df (pandas DataFrame): DataFrame detailing the validity
			of each barcode for each read
		counts (pandas DataFrame): The number of reads with
			each barcode
		count_threshold (int): Count threshold to use when
			correcting barcodes
	"""
	reads_in_cells_thresh = 0.92

	bc_8nt_bc3, bc_8nt_bc2, bc_8nt_bc1 = load_barcodes_set()

	# load in just bcs from everything
	if type(fnames) == list:
		df = pd.DataFrame()
		for i, fname in fnames:
			   temp = pd.read_csv(fname, sep='\t', usecols=[15,16,17,18])
			   df = pd.concat([df, temp])
	else:
		fname = fnames
		df = pd.read_csv(fname, sep='\t', usecols=[15,16,17,18])

	df['bc1_valid'] = df.bc1.isin(list(bc_8nt_bc1))
	df['bc2_valid'] = df.bc2.isin(list(bc_8nt_bc2))
	df['bc3_valid'] = df.bc3.isin(list(bc_8nt_bc3))

	n_bc1_valid = len(df.loc[df.bc1_valid == True].index)
	n_bc2_valid = len(df.loc[df.bc2_valid == True].index)
	n_bc3_valid = len(df.loc[df.bc3_valid == True].index)
	n_all_valid = len(df.loc[(df[['bc1_valid', 'bc2_valid', 'bc3_valid']].all(axis=1))])

	if verbose > 0:
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

def correct_barcodes(fname, oprefix,
					 counts, count_thresh,
					 bc_edit_dist=3, t=1,
					 chunksize=10**5,
					 verbose=1,
					 delete_input=False):
	"""
	Correct barcodes based on abundance of those already
	present in the dataset and on a maximum edit distance.

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		counts (pandas DataFrame): Output from get_perfect_bc_counts
			with number of reads observed for each barcode
		count_thresh (int): Minimum number of reads
		bc_edit_dist (int): Maximum edit distance of bc to correct
		t (int): Number of threads to run on
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	bc3_dict, bc2_dict, bc1_dict = load_barcodes()

	i = 0
	n_corrected_bcs = 0
	ofile = '{}_seq_corrected_bcs.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		if t == 1:
			temp = df.parallel_apply(correct_bcs_x,
					args=(counts, count_thresh, bc_edit_dist, bc1_dict, bc2_dict, bc3_dict),
					axis=1, result_type='expand')
		else:
			temp = df.apply(lambda x: correct_bcs_x(x, counts, count_thresh,
						bc_edit_dist, bc1_dict, bc2_dict, bc3_dict),
						axis=1, result_type='expand')
			pandarallel.initialize(nb_workers=t, verbose=1)

		# replace old bcs with corrected bcs and remove bcs
		# that were not corrected
		df.drop(['bc1', 'bc2', 'bc3'], axis=1, inplace=True)
		df = pd.concat([df, temp], axis=1)
		df.dropna(subset=['bc1', 'bc2', 'bc3'],
				  how='any', axis=0, inplace=True)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize
		n_corrected_bcs += len(df.index)
		if verbose == 2:
			print('Corrected bcs for {} reads'.format(i))

	if verbose > 0:
		print('Corrected barcodes for {} reads'.format(n_corrected_bcs))

	# delete input file
	if delete_input:
		os.remove(fname)

	return ofile

# trim off the barcodes from the start of bc1 until the end of the read
def trim_bcs(fname, oprefix, t=1,
			 chunksize=10**5, verbose=1,
			 delete_input=True):
	"""
	Trim BC + UMI construct from each read (from start of bc1 to end)

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		t (int): Number of threads to run on
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	i = 0
	ofile = '{}_trimmed.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		# add placeholders
		df.trim_start = np.nan

		# separate into fwd/rev
		fwd = df.loc[df.l_dir == '+'].copy(deep=True)
		rev = df.loc[df.l_dir == '-'].copy(deep=True)

		# where should we start trimming the sequence?
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
			pandarallel.initialize(nb_workers=t, verbose=1)
			df['trim_seq'] =  df.parallel_apply(trim_bcs_x, axis=1)

		# remove sequences that are empty now because of weird linker things
		df['trim_seq_len'] = df.apply(lambda x: len(x.trim_seq), axis=1)
		df = df.loc[df.trim_seq_len != 0]

		if verbose == 2:
			n = len(df.index)
			print('Number of reads after removing truncated sequences: {}'.format(n))

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize
		if verbose  == 2:
			print('Trimmed {} reads'.format(i))

	# delete input file
	if delete_input:
		os.remove(fname)

	return ofile

def flip_reads(fname, oprefix, t=1,
			   verbose=1, chunksize=10**5,
			   delete_input=False):
	"""
	Trim BC + UMI construct from each read (from start of bc1 to end)

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		t (int): Number of threads to run on
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	i = 0
	ofile = '{}_trimmed_flipped.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):
		if t == 1:
			df['trim_seq'] = df.apply(flip_reads_x, axis=1)
		else:
			pandarallel.initialize(nb_workers=t, verbose=1)
			df['trim_seq'] = df.parallel_apply(flip_reads_x, axis=1)

		df.seq = df.trim_seq
		df.drop('trim_seq', axis=1, inplace=True)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize
		if verbose == 2:
			print('Flipped {} reads'.format(i))

	# delete input file
	if delete_input:
		os.remove(fname)

	return ofile

def filter_dupe_umis(fname, oprefix, verbose=True,
					 chunksize=10**5,
					 delete_input=True):
	"""
	Filters duplicate UMIs such that the longest read is kept

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		t (int): Number of threads to run on
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	i = 0
	ofile = '{}_seq_umi_len.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		# get combined bc
		df['bc'] = df.bc3+df.bc2+df.bc1

		# calc read length and sort reads from longest to smallest
		df['read_len']  = df.seq.str.len()
		df = df.sort_values(by='read_len', ascending=False)

		# get umi len so we can leave incomplete umis alone
		df['umi_len'] = df.umi.str.len()

		df = df[['read_name', 'seq', 'read_len',
				 'umi_len', 'bc',
				 'bc1', 'bc2', 'bc3', 'umi']]
		df.reset_index(inplace=True)
		df.rename({'index': 'row'}, axis=1, inplace=True)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize

	# read whole thing in
	df = pd.read_csv(ofile, sep='\t')
	if verbose > 0:
		print('Number of reads before filtering dupe UMIs: {}'.format(len(df.index)))

	# drop dupe UMI + bc combos for full len umis
	inds = df.loc[df.umi_len == 10][['umi', 'bc']].drop_duplicates(keep='first').index.tolist()
	inds += df.loc[df.umi_len != 10].index.tolist()
	df = df.loc[inds]

	# remove unnecessary columns
	df.drop(['umi_len', 'read_len'], axis=1, inplace=True)

	# now loop through the just saved file in chunks and remove
	# those that do not belong to the deduplicated indices
	fname = ofile
	i = 0
	n_dedup = 0
	ofile = '{}_dedup_umi.tsv'.format(oprefix)
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		# unique UMI + BC combos
		df = df.loc[df.index.isin(inds)]

		n_dedup += len(df.index)

		# first write
		if i == 0:
			df.to_csv(ofile, sep='\t', index=False)
		else:
			df.to_csv(ofile, sep='\t', header=None, index=False, mode='a')

		i += chunksize

	if verbose > 0:
		print('Number of reads after filtering dupe UMIs: {}'.format(n_dedup))

	# remove len tracking info
	os.remove(fname)

	return ofile

# format/write fastq
def write_fastq(fname, oprefix,
				chunksize=10**5,
				delete_input=False):
	"""
	Write a new fastq with cell barcode + umi info in each
	read name

	Parameters:
		fname (str): File to process
		oprefix (str): Where to save output
		verbose (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		ofile (str): Name of output file
	"""

	# output fastq
	ofile = oprefix+'_demux.fastq'
	ofile = open(ofile, 'w')

	i = 0
	for df in pd.read_csv(fname, chunksize=chunksize, sep='\t'):

		# create the read name header with the bc and umi information
		df.fillna(value='_', inplace=True)
		ex_read_name = df.read_name.tolist()[0]
		if ' ' in ex_read_name:
			df.read_name = df.read_name.str.split(' ', n=1, expand=True)[0]
		df['header'] = '@'+df.read_name+':'+df.bc1+'_'+df.bc2+'_'+df.bc3+'_'+df.umi
		df = df[['header', 'seq']]

		# write the fastq
		for ind, entry in df.iterrows():
			try:
				ofile.write(entry.header+'\n')
			except:
				print(entry.header)

			ofile.write(entry.seq+'\n')
			ofile.write('+\n')
			ofile.write(''.join(['5' for i in range(len(entry.seq))])+'\n')

	ofile.close()
	return ofile

###################################################################################
############################### Lambda functions ##################################
###################################################################################

def score_linkers_x(x, l_seqs, l_prefs):
	"""
	Function to apply across rows of a dataframe

	Parameters:
		x (pandas Series): Row from parent dataframe
		l_seqs (list of str): Linker sequences
		l_prefs (list of str): Names of linkers (l1, l2)

	Returns:
		entry (pandas Series): Row with linker alignment scores
	"""

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

def align_linkers_x(x, l1, l2, l1_rc, l2_rc):
	"""
	Function to find inds of linkers across rows of a dataframe

	Parameters:
		x (pandas Series): Row from parent dataframe
		l1 (str): Linker 1
		l2 (str): Linker 2
		l1_rc (str): Linker 1 reverse complement
		l2_rc (str): Linker 2 reverse complement

	Returns:
		entry (pandas Series): Row with linker alignment indices
	"""
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

def get_bcs_umis_x(x):
	"""
	Function to find bcs / umis of reads across rows of a dataframe

	Parameters:
		x (pandas Series): Row from parent dataframe

	Returns:
		entry (pandas Series): Row with extracted barcodes / umis
	"""
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

	# construct an entry
	entry = {}
	entry['bc1'] = bc1
	entry['bc2'] = bc2
	entry['bc3'] = bc3
	entry['umi'] = umi
	return entry

# From the Parse Biosciences pipeline
def correct_bcs_x(x,
		counts, count_thresh,
		bc_edit_dist,
		bc1_dict,
		bc2_dict,
		bc3_dict):
	"""
	Function to correct bcs across rows of a dataframe.
	Adapted from Parse biosciences

	Parameters:
		x (pandas Series): Row from parent dataframe
		counts (pandas DataFrame): Output from get_perfect_bc_counts
			with number of reads observed for each barcode
		count_thresh (int): Minimum number of reads
		bc_edit_dist (int): Maximum edit distance of bc to correct
		bc#_edit_dict (dict): Dict for barcode<1,2,3> with
			key: query bc
			item: corrected bc

	Returns:
		entry (pandas Series): Row with corrected bcs
	"""

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

def trim_bcs_x(x):
	"""
	Function to trim reads across each read in dataframe

	Parameters:
		x (pandas Series): Row from parent dataframe

	Returns:
		entry (pandas Series): Row with trimmed reads
	"""
	if x.l_dir == '+':
		trim_seq = x.seq[x.trim_start:]
	elif x.l_dir == '-':
		trim_seq = x.seq[:x.trim_start]
	return trim_seq

def flip_reads_x(x):
	"""
	Function to flip reads across each read in dataframe

	Parameters:
		x (pandas Series): Row from parent dataframe

	Returns:
		entry (pandas Series): Row with flipped reads
	"""
	if x.l_dir == '+':
		s = rev_comp(x.trim_seq)
	else:
		s = x.trim_seq
	return s
