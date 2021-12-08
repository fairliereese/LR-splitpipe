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
	parser.add_argument('--delete_tmp', dest='delete_tmp',
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
	delete_tmp = args.delete_tmp
	chunk = int(args.chunksize)

	fname = fastq_to_df(fastq,
				oprefix,
				verbose=v)

	fname = score_linkers(fname, oprefix,
					 t=t,
					 verbose=v,
					 chunksize=chunk,
					 delete_input=delete_tmp)

	fname = align_linkers(fname, oprefix,
					l1_m=l1_mm, l2_m=l2_mm,
					t=t,
					verbose=v,
					chunksize=chunk,
					delete_input=delete_tmp)

	fname = get_bcs_umis(fname, oprefix,
				 t=t,
				 verbose=v,
				 chunksize=chunk,
				 delete_input=delete_tmp)

	_, counts, count_thresh = get_perfect_bc_counts(fname, verbose=v)

	fname = correct_barcodes(fname, oprefix,
						 counts, count_thresh,
						 bc_edit_dist=3,
						 t=t,
						 verbose=v,
						 chunksize=chunk,
						 delete_input=delete_tmp)

	fname = trim_bcs(fname, oprefix,
				t=t,
				verbose=v,
				chunksize=chunk,
				delete_input=delete_tmp)

	fname = flip_reads(fname, oprefix,
				t=t,
				verbose=v,
				chunksize=chunk,
				delete_input=delete_tmp)

	fname = filter_dupe_umis(fname, oprefix,
				verbose=v,
				chunksize=chunk,
				delete_input=delete_tmp)

	fname = write_fastq(fname, oprefix,
				chunksize=chunk,
				delete_input=delete_tmp)

if __name__ == '__main__': main()
