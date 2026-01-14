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
import pdb

from utils import *
from plotting import *

###################################################################################
############################ Actual function calls ################################
###################################################################################
def score(fastq, oprefix, t,
					  chunksize, verbosity,
					  delete_input):
	"""
	Runs the steps of the pipeline that can be done in parallel including
		* processing the fastq into a dataframe
		* linker scoring and locating
		* bc calling

	Parameters:
		fastq (str): File to process
		oprefix (str): Where to save output
		verbosity (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		t (int): Number of threads to run on
		keep_dupes (bool):
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		fname (str): Last file that is created
	"""

	# portion of the pipeline that can be run in parallel
	fname = fastq_to_df(fastq,
				oprefix,
				verbose=verbosity)

	fname = score_linkers(fname, oprefix,
					 t=t,
					 verbose=verbosity,
					 chunksize=chunksize,
					 delete_input=delete_input)

	# make some plots
	df = pd.read_csv(fname, sep='\t', usecols=[3,4,5,6,7])
	df.reset_index(inplace=True)
	plot_post_score_plots(df, oprefix)

def find_bcs(fastq, oprefix, t,
					  l1_mm, l2_mm,
					  max_dist,
					  max_len, min_len,
					  chunksize, verbosity,
					  delete_input):
	"""
	Runs the steps of the pipeline that can be done in parallel including
		* processing the fastq into a dataframe
		* linker scoring and locating
		* bc calling

	Parameters:
		fastq (str): File to process
		oprefix (str): Where to save output
		verbosity (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		t (int): Number of threads to run on
		l1_m (int): Number of allowable mismatches in linker 1
		l2_m (int): Number of allowable mismatches in linker 2
		max_dist (int): Max. distance that a linker can be from the end of the read
		max_len (int): Max. length of a read to be considered
		min_len (int): Min. length of a read to be considered
		l1_p (float): Proportion of allowable mismatches in linker 1
		l2_p (float): Proportion of allowable mismatches in linker 2
		keep_dupes (bool):
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution

	Returns:
		fname (str): Last file that is created
	"""

	# portion of the pipeline that can be run in parallel
	fname = fastq_to_df(fastq,
				oprefix,
				verbose=verbosity)

	fname = score_linkers(fname, oprefix,
					 t=t,
					 verbose=verbosity,
					 chunksize=chunksize,
					 delete_input=delete_input)

	# make some plots
	df = pd.read_csv(fname, sep='\t', usecols=[3,4,5,6,7])
	df.reset_index(inplace=True)
	plot_post_score_plots(df, oprefix)

	fname = align_linkers(fname, oprefix,
					l1_m=l1_mm, l2_m=l2_mm,
					max_dist=max_dist,
					max_len=max_len,
					min_len=min_len,
					t=t,
					verbose=verbosity,
					chunksize=chunksize,
					delete_input=delete_input)

	# make linker dist plots
	df = pd.read_csv(fname, sep='\t', usecols=[3,4,5,6,9,10])
	fwd, rev = get_fwd_rev(df)
	fwd = get_linker_dists(fwd, '+')
	rev = get_linker_dists(rev, '-')
	df = pd.concat([rev, fwd], axis=0)
	plot_linker_dists(df, oprefix)

	fname = get_bcs_umis(fname, oprefix,
				 t=t,
				 verbose=verbosity,
				 chunksize=chunksize,
				 delete_input=delete_input)

	return fname

def process_bcs(fnames, oprefix,
				  kit, chemistry, t,
				  chunksize, verbosity,
				  delete_input):

	"""
	Runs the steps of the pipeline that should be done all at once including
	* barcode correction
	* read trimming + flipping
	* filtering duplicate umis
	* writing fastq

	Parameters:
		fnames (list of str): Files to process
		oprefix (str): Where to save output
		kit (str): Which kit was used
		chemistry (str): Which chemistry was used
		verbosity (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		t (int): Number of threads to run on
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution
	"""
	_, counts, count_thresh = get_perfect_bc_counts(fnames, kit,
													chemistry,
													verbose=verbosity)

	fname = correct_barcodes(fnames, oprefix,
						 kit, chemistry,
						 counts, count_thresh,
						 bc_edit_dist=3,
						 t=t,
						 verbose=verbosity,
						 chunksize=chunksize,
						 delete_input=False)

	fname = trim_bcs(fname, oprefix,
				t=t,
				verbose=verbosity,
				chunksize=chunksize,
				delete_input=delete_input)

	fname = flip_reads(fname, oprefix,
				t=t,
				verbose=verbosity,
				chunksize=chunksize,
				delete_input=delete_input)

	fname = filter_dupe_umis(fname, oprefix,
				verbose=verbosity,
				chunksize=chunksize,
				delete_input=delete_input)

	# plot read lengths as they are after trimming
	df = pd.read_csv('{}_seq_umi_len.tsv'.format(oprefix), sep='\t',
		usecols=[3,4,5,9])
	plot_read_length(df, oprefix+'_post_bc')

	# for UMI plots, only consider reads with FL UMIs
	kind = 'Post-correction'
	df = df.loc[df.umi_len == 10]
	temp = plot_knee_plot(df, oprefix, kind)
	plot_sequencing_sat(df, oprefix, kind)

	fname = write_fastq(fname, oprefix,
				chunksize=chunksize,
				delete_input=delete_input)

###################################################################################
############################# Argparsing functions ################################
###################################################################################
def get_args():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers(dest='mode')

	# all steps
	parser_all = subparsers.add_parser('all', help='Run all steps')
	parser_all.add_argument('-f', dest='fastq',
		help='FASTQ file output from Lima with LR-Split-seq reads.')
	parser_all.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser_all.add_argument('-k', dest='kit',
		help='Kit used, {WT, WT_mini, WT_mega}'),
	parser_all.add_argument('-c', dest='chemistry', default='v1',
		help='Chemistry used, {v1, v2}'),
	parser_all.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser_all.add_argument('--l1_mm', dest='l1_mm', default=3,
		help='Number of allowable mismatches in linker1')
	parser_all.add_argument('--l2_mm', dest='l2_mm', default=3,
		help='Number of allowable mismatches in linker2')
	parser_all.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in / process at a time')
	parser_all.add_argument('--max_linker_dist', dest='max_dist', default=None,
		help='Maximum distance that a linker can be from the end of a read')
	parser_all.add_argument('--max_read_len', dest='max_len', default=None)
	parser_all.add_argument('--min_read_len', dest='min_len', default=None)
	parser_all.add_argument('--verbosity', dest='verbosity', default=1,
		help="""Verbosity setting.
			    0: No output
				1: QC statistics
				2: QC statistics + progress""")
	parser_all.add_argument('--delete_input', dest='delete_input',
		action='store_true',
		help='Delete temporary files (recommended!)', default=False)
	# parser_all.add_argument('--filt_umi', dest='filt_umi', default=False,
	# 	help='Filter out duplicate UMIs using longest read heuristic')

	parser_score_linkers = subparsers.add_parser('score_linkers', help='Run score_linkers step only')
	parser_score_linkers.add_argument('-f', dest='fastq',
		help='FASTQ file output from Lima with LR-Split-seq reads.')
	parser_score_linkers.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser_score_linkers.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser_score_linkers.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in at any given time')
	parser_score_linkers.add_argument('--verbosity', dest='verbosity', default=1,
		help='Verbosity setting. Higher number = more messages')
	parser_score_linkers.add_argument('--delete_input', dest='delete_input',
		action='store_true', help='Delete temporary files', default=False)
	# parser_find_bcs.add_argument('--filt_umi', dest='filt_umi', default=False,
	# 	help='Filter out duplicate UMIs using longest read heuristic')

	# find bcs
	parser_find_bcs = subparsers.add_parser('find_bcs', help='Run find_bcs step')
	parser_find_bcs.add_argument('-f', dest='fastq',
		help='FASTQ file output from Lima with LR-Split-seq reads.')
	parser_find_bcs.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser_find_bcs.add_argument('-k', dest='kit',
		help='Kit used, {custom_1, WT, WT_mini, WT_mega}', default='WT')
	parser_find_bcs.add_argument('-c', dest='chemistry', default='v1',
		help='Chemistry used, {v1, v2}')
	parser_find_bcs.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser_find_bcs.add_argument('--l1_mm', dest='l1_mm', default=3,
		help='Number of allowable mismatches in linker1')
	parser_find_bcs.add_argument('--l2_mm', dest='l2_mm', default=3,
		help='Number of allowable mismatches in linker2')
	parser_find_bcs.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in at any given time')
	parser_find_bcs.add_argument('--max_linker_dist', dest='max_dist', default=None,
		help='Maximum distance that a linker can be from the end of a read')
	parser_find_bcs.add_argument('--max_read_len', dest='max_len', default=None)
	parser_find_bcs.add_argument('--min_read_len', dest='min_len', default=None)
	parser_find_bcs.add_argument('--verbosity', dest='verbosity', default=1,
		help='Verbosity setting. Higher number = more messages')
	parser_find_bcs.add_argument('--delete_input', dest='delete_input',
		action='store_true',
		help='Delete temporary files', default=False)
	# parser_find_bcs.add_argument('--filt_umi', dest='filt_umi', default=False,
	# 	help='Filter out duplicate UMIs using longest read heuristic')

	# process bcs
	parser_process_bcs = subparsers.add_parser('process_bcs', help='Run process_bcs step')
	parser_process_bcs.add_argument('-f', dest='fnames',
		help='Comma-separated list of files from "find_bcs" step with suffix "_bcs.tsv".')
	parser_process_bcs.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser_process_bcs.add_argument('-k', dest='kit',
		help='Kit used, {custom_1, WT, WT_mini, WT_mega}')
	parser_process_bcs.add_argument('-c', dest='chemistry', default='v1',
		help='Chemistry used, {v1, v2}')
	parser_process_bcs.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser_process_bcs.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in at any given time')
	parser_process_bcs.add_argument('--verbosity', dest='verbosity', default=1,
		help='Verbosity setting. Higher number = more messages')
	parser_process_bcs.add_argument('--delete_input', dest='delete_input',
		action='store_true',
		help='Delete temporary files', default=False)
	# parser_process_bcs.add_argument('--filt_umi', dest='filt_umi', default=False,
	# 	help='Filter out duplicate UMIs using longest read heuristic')

	args = parser.parse_args()
	return args

def main():
	args = get_args()
	mode = args.mode
	oprefix = args.oprefix
	t = int(args.threads)
	v = int(args.verbosity)
	delete_input = args.delete_input

	def format_chunksize(c):
		if isinstance(c, str) and '**' in c:
			i, j = c.split('**')
			i = float(i)
			j = float(j)
			c = i**j
		c = int(c)
		return c

	chunksize = format_chunksize(args.chunksize)

	if mode == 'all' or mode == 'find_bcs' or mode == 'score_linkers':
		fastq = args.fastq
		kit = args.kit
		chemistry = args.chemistry
		if args.max_dist:
			max_dist = int(args.max_dist)
		else:
			max_dist = None
		if args.max_len:
			max_len = int(args.max_len)
		else:
			max_len = None
		if args.min_len:
			min_len = int(args.min_len)
		else:
			min_len = None

		if mode == 'all' or mode == 'find_bcs':
			l1_mm = int(args.l1_mm)
			l2_mm = int(args.l2_mm)
	elif mode == 'process_bcs':
		kit = args.kit
		chemistry = args.chemistry
		fnames = args.fnames
		fnames = fnames.split(',')

	if mode == 'all' or mode == 'find_bcs':
		fname = find_bcs(fastq, oprefix, t,
								  l1_mm, l2_mm,
								  max_dist,
								  max_len, min_len,
								  chunksize, v,
							  	  delete_input)

		if mode == 'all':
			fname = process_bcs(fname, oprefix,
										 kit, chemistry, t,
										 chunksize, v,
										 delete_input)

	elif mode == 'process_bcs':
		fname = process_bcs(fnames, oprefix,
									 kit, chemistry, t,
									 chunksize, v,
									 delete_input)

	elif mode == 'score_linkers':
		score(fastq, oprefix, t,
						  chunksize, v,
						  delete_input)

if __name__ == '__main__': main()
