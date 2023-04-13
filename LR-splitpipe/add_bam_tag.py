import pandas as pd
import argparse
import os
from utils import *

def get_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('-s', dest='samfile',
		help='SAM file with Split-seq barcode+UMI information in the read name')
	parser.add_argument('-k', dest='kit',
		help='Kit used {WT, WT_mini, WT_mega}')
	parser.add_argument('--merge_primers', dest='merge_primers',
		default=False, action='store_true',
		help='Merge reads that come from the same cell from different priming strategies')
	parser.add_argument('--suffix', dest='suff', default=None,
		help='Suffix to add to cell barcodes. Useful if merging separate LR-Split-seq experiments '+
		'that might have overlapping barcodes otherwise.')
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')

	args = parser.parse_args()
	return args

def get_bc1_matches(kit):
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    # pkg_path = '/'.join(pkg_path.split('/')[:-1])
    bc_round_set = get_bc_round_set('WT')

    # determine file to use
    for entry in bc_round_set:
        if entry[0] == 'bc1':
            ver = entry[1]

    # read in and restructure such that each dt bc is
    # matched with its randhex partner from the same well
    fname = pkg_path+'/barcodes/bc_data_{}.csv'.format(ver)
    df = pd.read_csv(fname)
    df.loc[df.well.duplicated(keep=False)].sort_values(by='well')
    drop_cols = ['bci', 'uid', 'type']
    bc1_dt = df.loc[df['type'] == 'T'].drop(drop_cols, axis=1)
    bc1_dt.rename({'sequence': 'bc1_dt'}, axis=1, inplace=True)
    bc1_randhex = df.loc[df['type'] == 'R'].drop(drop_cols, axis=1)
    bc1_randhex.rename({'sequence': 'bc1_randhex'}, axis=1, inplace=True)
    bc_df = bc1_dt.merge(bc1_randhex, on='well')

    return bc_df

def get_read_info(line):
	''' From a line in a sam file, returns the read name,
	    barcode, and UMI as formatted by demultiplex.py
	'''
	read_bc = line[0] # sam line
	read_bc = read_bc.split(':')
	read_name = read_bc[0]
	bc_umi = read_bc[1]
	bc_umi = bc_umi.split('_')
	bc = ''.join(bc_umi[0:-1][::-1])
	umi =  bc_umi[-1]

	return read_name, bc, umi

def main():
	args = get_args()
	samfile = args.samfile
	kit = args.kit
	oprefix = args.oprefix
	suff = args.suff

	merge_primers = args.merge_primers

	if merge_primers:
		bc_df = get_bc1_matches(kit)
		fname = '{}_merged_primers.sam'.format(oprefix)
		print(bc_df.head())
	else:
		fname = '{}.sam'.format(oprefix)

	ofile = open(fname, 'w')
	ifile = open(samfile, 'r')

	for line in ifile:
		if line.startswith('@'):
			ofile.write(line)
		else:
			line = line.strip().split('\t')
			read_name, bc, umi = get_read_info(line)

			# replace bc with merged bc if necessary
			if merge_primers:
				bc3 = bc[:8]
				bc2 = bc[8:16]
				bc1 = bc[16:]
				if bc1 in bc_df.bc1_randhex.tolist():
					bc1 = bc_df.loc[bc_df.bc1_randhex==bc1, 'bc1_dt'].values[0]
					bc = bc3+bc2+bc1

			line[0] = read_name

			if suff:
				bc = '{}-{}'.format(bc, suff)

			cell_tag = 'CB:Z:{}'.format(bc)
			umi_tag = 'MI:Z:{}'.format(umi)
			line.append(cell_tag)
			line.append(umi_tag)
			line = '\t'.join(line)+'\n'
			ofile.write(line)
	ofile.close()
	ifile.close()

if __name__ == '__main__': main()
