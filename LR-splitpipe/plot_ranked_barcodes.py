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
