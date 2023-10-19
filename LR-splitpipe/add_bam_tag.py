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
	parser.add_argument('-c', dest='chemistry', default='v1',
		help='Chemistry used, {v1, v2}')
	parser.add_argument('--merge_primers', dest='merge_primers',
		default=False, action='store_true',
		help='Merge reads that come from the same cell from different priming strategies')
	parser.add_argument('--min_umi', dest='min_umi',
		default=0, help='Minimum number of UMIs/cell to retain reads')
	parser.add_argument('--suffix', dest='suff', default=None,
		help='Suffix to add to cell barcodes. Useful if merging separate LR-Split-seq experiments '+
		'that might have overlapping barcodes otherwise.')
	parser.add_argument('-t', dest='threads',
		help='Number of threads to run on.', default=1)
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')

	args = parser.parse_args()
	return args

def get_bc1_matches(kit, chemistry):
    pkg_path = os.path.dirname(os.path.dirname(__file__))
    # pkg_path = '/'.join(pkg_path.split('/')[:-1])
    bc_round_set = get_bc_round_set(kit, chemistry)

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

	# def get_read_info(raw_read_name, merge_primers, bc_df, suff):
	# 	''' From a line in a sam file, returns the read name,
	# 	    barcode, and UMI as formatted by demultiplex.py
	# 		Merge primers if requested. Add suffix as requested.
	# 	'''
	# 	read_bc = raw_read_name.split(':')
	# 	read_name = read_bc[0]
	# 	bc_umi = read_bc[1]
	# 	bc_umi = bc_umi.split('_')
	# 	bc = ''.join(bc_umi[0:-1][::-1])
	# 	umi =  bc_umi[-1]
	#
	# 	# replace bc with merged bc if necessary
	# 	if merge_primers:
	# 		bc3 = bc[:8]
	# 		bc2 = bc[8:16]
	# 		bc1 = bc[16:]
	# 		if bc1 in bc_df.bc1_randhex.tolist():
	# 			bc1 = bc_df.loc[bc_df.bc1_randhex==bc1, 'bc1_dt'].values[0]
	# 			bc = bc3+bc2+bc1
	#
	# 	if suff:
	# 		bc = '{}-{}'.format(bc, suff)
	#
	# 	return raw_read_name, read_name, bc, umi

def main():
	args = get_args()
	samfile = args.samfile
	kit = args.kit
	chemistry = args.chemistry
	oprefix = args.oprefix
	suff = args.suff

	merge_primers = args.merge_primers
	min_umi = int(args.min_umi)
	threads = int(args.threads)

	if merge_primers:
		fname = '{}_merged_primers.sam'.format(oprefix)
	else:
		fname = '{}.sam'.format(oprefix)

	bc_df = get_bc1_matches(kit, chemistry)

	# get the updated info for each read in the sam file
	ifile = open(samfile, 'r')
	read_names = []
	for line in ifile:
		if line.startswith('@'):
			continue
		read_names.append(line.split('\t')[0])

	df = pd.DataFrame(data=read_names, columns=['raw_read_name'])
	# df[['read_name', 'bc_umi']] = df['raw_read_name'].str.split(':', expand=True)
	df[['read_name', 'bc_umi']] = df['raw_read_name'].str.rsplit(':', n=1, expand=True)
	df[['bc1', 'bc2', 'bc3', 'umi']] = df['bc_umi'].str.split('_', n=3, expand=True)
	df['bc'] = df['bc3']+df['bc2']+df['bc1']

	# get merged versions of bcs
	if merge_primers:
	    df = df.merge(bc_df, how='left', left_on='bc1', right_on='bc1_randhex')
	    rand_inds = df.loc[~df.bc1_randhex.isnull()].index
	    df.loc[rand_inds, 'bc1_merge'] = df.loc[rand_inds, 'bc1_dt']
	    dt_inds = df.loc[df.bc1_randhex.isnull()].index
	    df.loc[dt_inds, 'bc1_merge'] = df.loc[dt_inds, 'bc1']
	    df['bc'] = df.bc3+df.bc2+df.bc1_merge

	# add suffix
	if suff:
	    df['bc'] = df.bc+'-'+suff

	# if we're imposing min. umi / cell or nuc figure that out
	if min_umi != 0:
		df2 = df[['bc', 'umi']].groupby('bc').nunique().reset_index().rename({'umi':'n_umi'}, axis=1)
		df2 = df2.loc[df2.n_umi >= min_umi]
		# read_info = df.loc[df.bc.isin(df2.bc.tolist())]
		df['pass'] = df.bc.isin(df2.bc.tolist())
		n_pass = len(df.loc[df['pass']==True].index)
		print(f'{n_pass} reads that belong to cells w/ >= {min_umi} umi')
	else:
		df['pass'] = True
	# else:
	# 	read_info = df

	# table of updated read name etc
	df = df[['raw_read_name', 'read_name', 'bc', 'umi', 'pass']]
	ifile.close()

	# now loop through reads again and replace / add info as needed
	# ifile = pysam.AlignmentFile(samfile, 'r', threads=threads)
	# ofile = pysam.AlignmentFile(fname, 'w', threads=threads, template=ifile)
	ifile = open(samfile, 'r')
	ofile = open(fname, 'w')

	df_ind = 0
	for line in ifile:
		if line.startswith('@'):
			ofile.write(line)
		else:

			# df contains one entry for each read, read_info contains only
			# the reads we want to output
			read_name = line.split('\t')[0]
			temp = df.iloc[df_ind]

			if temp['pass'] == True:

				raw_read_name = temp.raw_read_name
				assert read_name == raw_read_name

				new_read_name = temp.read_name
				bc = temp.bc
				umi = temp.umi

				# make new sam line
				line = line.strip().split('\t')
				line[0] = new_read_name
				cell_tag = 'CB:Z:{}'.format(bc)
				umi_tag = 'MI:Z:{}'.format(umi)
				line.append(cell_tag)
				line.append(umi_tag)

				# write out
				line = '\t'.join(line)+'\n'
				ofile.write(line)

			df_ind += 1

	ofile.close()
	ifile.close()

if __name__ == '__main__': main()
