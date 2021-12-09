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

from utils import *

###################################################################################
############################ Actual function calls ################################
###################################################################################
def find_bcs_function(fastq, oprefix, t,
					  l1_mm, l2_mm,
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

	fname = align_linkers(fname, oprefix,
					l1_m=l1_mm, l2_m=l2_mm,
					t=t,
					verbose=verbosity,
					chunksize=chunksize,
					delete_input=delete_input)

	fname = get_bcs_umis(fname, oprefix,
				 t=t,
				 verbose=verbosity,
				 chunksize=chunksize,
				 delete_input=delete_input)

	return fname

def process_bcs_function(fnames, oprefix, t,
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
		verbosity (int): How much output to show
			0: none
			1: only QC statistics
			2: QC stats + progress
		t (int): Number of threads to run on
		chunksize (int): Number of lines to process at a time
		delete_input (bool): Whether or not to delete input file
			after execution
	"""

	_, counts, count_thresh = get_perfect_bc_counts(fnames, verbose=verbosity)

	fname = correct_barcodes(fname, oprefix,
						 counts, count_thresh,
						 bc_edit_dist=3,
						 t=t,
						 verbose=verbosity,
						 chunksize=chunksize,
						 delete_input=delete_input)

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

	fname = write_fastq(fname, oprefix,
				chunksize=chunksize,
				delete_input=delete_input)

###################################################################################
############################# Argparsing functions ################################
###################################################################################
def find_bcs():
	parser = argparse.ArgumentParser()
	parser.add_argument('find_bcs')
	parser.add_argument('-f', dest='fastq',
		help='FASTQ file output from Lima with LR-Split-seq reads.')
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser.add_argument('--l1_mm', dest='l1_mm', default=3,
		help='Number of allowable mismatches in linker1')
	parser.add_argument('--l2_mm', dest='l2_mm', default=3,
		help='Number of allowable mismatches in linker2')
	parser.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in at any given time')
	parser.add_argument('--verbosity', dest='verbosity', default=1,
		help='Verbosity setting. Higher number = more messages')
	parser.add_argument('--delete_input', dest='delete_input',
		help='Delete temporary files', default=False)

	args = parser.parse_args()

	fastq = args.fastq
	oprefix = args.oprefix
	t = int(args.threads)
	l1_mm = int(args.l1_mm)
	l2_mm = int(args.l2_mm)
	v = int(args.verbosity)
	delete_input = args.delete_input
	chunk = int(args.chunksize)

	# run just the bc finding
	find_bcs_function(fastq, oprefix, t,
						  l1_mm, l2_mm,
						  chunksize, verbosity,
						  delete_input)

def process_bcs():
	parser = argparse.ArgumentParser()
	parser.add_argument('process_bcs')

	parser.add_argument('-f', dest='fnames',
		help='Comma-separated list of files from "find_bcs" step with suffix "seq_linker_alignments.tsv".')
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in at any given time')
	parser.add_argument('--verbosity', dest='verbosity', default=1,
		help='Verbosity setting. Higher number = more messages')
	parser.add_argument('--delete_input', dest='delete_input',
		help='Delete temporary files', default=False)

	args = parser.parse_args()

	fnames = args.fnames
	fnames = fnames.split(',')

	oprefix = args.oprefix
	t = int(args.threads)
	v = int(args.verbosity)
	delete_input = args.delete_input
	chunk = int(args.chunksize)

	fname = process_bcs_function(fname, oprefix, t,
								 chunksize, v,
								 delete_input)

def get_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('-f', dest='fastq',
		help='FASTQ file output from Lima with LR-Split-seq reads.')
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')
	parser.add_argument('-t', dest='threads',
		help='Number of threads to run on (multithreading is recommended)')
	parser.add_argument('--l1_mm', dest='l1_mm', default=3,
		help='Number of allowable mismatches in linker1')
	parser.add_argument('--l2_mm', dest='l2_mm', default=3,
		help='Number of allowable mismatches in linker2')
	parser.add_argument('--chunksize', dest='chunksize', default=10**5,
		help='Number of lines to read in at any given time')
	parser.add_argument('--verbosity', dest='verbosity', default=1,
		help='Verbosity setting. Higher number = more messages')
	parser.add_argument('--delete_input', dest='delete_input',
		help='Delete temporary files', default=False)
	# parser.add_argument('--filt_umi', dest='filt_umi', default=False,
	# 	help='Filter out duplicate UMIs using longest read heuristic')

	args = parser.parse_args()
	return args

def main():
	args = get_args()
	fastq = args.fastq
	oprefix = args.oprefix
	t = int(args.threads)
	l1_mm = int(args.l1_mm)
	l2_mm = int(args.l2_mm)
	v = int(args.verbosity)
	delete_input = args.delete_input
	chunksize = int(args.chunksize)

	fname = find_bcs_function(fastq, oprefix, t,
							  l1_mm, l2_mm,
							  chunksize, v,
						  	  delete_input)

	fname = process_bcs_function(fname, oprefix, t,
								 chunksize, v,
								 delete_input)

	# # portion of the pipeline that can be run in parallel
	# fname = fastq_to_df(fastq,
	# 			oprefix,
	# 			verbose=v)
	#
	# fname = score_linkers(fname, oprefix,
	# 				 t=t,
	# 				 verbose=v,
	# 				 chunksize=chunk,
	# 				 delete_input=delete_input)
	#
	# fname = align_linkers(fname, oprefix,
	# 				l1_m=l1_mm, l2_m=l2_mm,
	# 				t=t,
	# 				verbose=v,
	# 				chunksize=chunk,
	# 				delete_input=delete_input)
	#
	# fname = get_bcs_umis(fname, oprefix,
	# 			 t=t,
	# 			 verbose=v,
	# 			 chunksize=chunk,
	# 			 delete_input=delete_input)

	# _, counts, count_thresh = get_perfect_bc_counts(fname, verbose=v)
	#
	# fname = correct_barcodes(fname, oprefix,
	# 					 counts, count_thresh,
	# 					 bc_edit_dist=3,
	# 					 t=t,
	# 					 verbose=v,
	# 					 chunksize=chunk,
	# 					 delete_input=delete_input)
	#
	# fname = trim_bcs(fname, oprefix,
	# 			t=t,
	# 			verbose=v,
	# 			chunksize=chunk,
	# 			delete_input=delete_input)
	#
	# fname = flip_reads(fname, oprefix,
	# 			t=t,
	# 			verbose=v,
	# 			chunksize=chunk,
	# 			delete_input=delete_input)
	#
	# fname = filter_dupe_umis(fname, oprefix,
	# 			verbose=v,
	# 			chunksize=chunk,
	# 			delete_input=delete_input)
	#
	# fname = write_fastq(fname, oprefix,
	# 			chunksize=chunk,
	# 			delete_input=delete_input)

if __name__ == '__main__': main()
