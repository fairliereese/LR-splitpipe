import pandas as pd
import argparse

def get_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('-s', dest='samfile',
		help='SAM file output from Minimap2/TranscriptClean with splitseq '+\
			 'barcode+UMI information in the read name')
	parser.add_argument('-o', dest='oprefix',
		help='Output file path/prefix')

	args = parser.parse_args()
	return args

def get_read_info(line):
''' From a line in a sam file, returns the read name, 
    barcode, and UMI as formatted by demultiplex.py '''
	read_bc = line[0]
	read_bc = read_bc.split(':')
	read_name = read_bc[0]
	bc_umi = read_bc[1]
	bc_umi = bc_umi.split('_')
	bc = bc_umi[0]
	umi =  bc_umi[1]

	return read_name, bc, umi


def main():
	args = get_args()
	samfile = args.samfile
	oprefix = args.oprefix

	fname = '{}.sam'
	ofile = open(fname, 'w')
	ifile = open(samfile, 'r')

	for line in ifile:
		if line.startswith('@'):
			ofile.write(line)
		else:
			line = line.strip().split('\t')
			read_name, bc, umi = get_read_info(line)
			line[0] = read_name
			cell_tag = 'RG:Z:{}'.format(bc)
			umi_tag = 'MI:Z:{}'.format(umi)
			line.append(cell_tag)
			line.append(umi)
			line = '\t'.join(line)+'\n'
			ofile.write(line)
	ofile.close()
	ifile.close()

if __name__ == '__main__': main()