import pandas as pd
import argparse
import os

def get_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('-s', dest='samfile',
		help='SAM file output from Minimap2/TranscriptClean with splitseq '+\
			 'barcode+UMI information in the read name')
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

def get_bc1_matches():
	pkg_path = os.path.dirname(__file__)
	pkg_path = '/'.join(pkg_path.split('/')[:-1])
	bc_file = pkg_path+'/barcodes/bc_8nt_v2.csv'
	bc_df = pd.read_csv(bc_file, index_col=0, names=['bc'])
	bc_df['well'] = [i for i in range(0, 48)]+[i for i in range(0, 48)]
	bc_df['primer_type'] = ['dt' for i in range(0, 48)]+['randhex' for i in range(0, 48)]

	# pivot on well to get df that matches bcs with one another from the same well
	bc_df = bc_df.pivot(index='well', columns='primer_type', values='bc')
	bc_df = bc_df.rename_axis(None, axis=1).reset_index()
	bc_df.rename({'dt': 'bc1_dt', 'randhex': 'bc1_randhex'}, axis=1, inplace=True)

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
	oprefix = args.oprefix
	suff = args.suff

	merge_primers = args.merge_primers

	if merge_primers:
		bc_df = get_bc1_matches()
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
