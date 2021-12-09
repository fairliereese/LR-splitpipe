import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math


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
    cols = ['index'] + val_vars
    temp = df[cols].melt(id_vars='index', value_vars=val_vars)

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
	fig = plt.gcf()
	fig.set_size_inches(8, 8)
	plt.savefig(fname, bbox_inches='tight')
	plt.clf()
