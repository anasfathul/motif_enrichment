#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import sys

args = sys.argv

# Import metadata (pass by args perhaps?)
motifs_metadata = pd.read_table('https://resources.altius.org/~jvierstra/projects/motif-clustering-v2.1beta/metadata.tsv')

# Import all coefficients dataset (Change to by args)
coeffs = pd.read_csv('/net/seq/data2/projects/afathul/motif_enhancement/output/all.coeff.tsv')

# Cluster name by Serj
cluster_names = {}
with open('/home/jvierstra/proj/code/motif-clustering/results/consensus_pwms.uniprobe') as f:
    for line in f:
        if line.startswith('AC'):
            cluster_names[line.split(':')[0]] = line.split(':')[1]
motifs_metadata['cluster_name'] = motifs_metadata.apply(lambda row: cluster_names[row['cluster']], axis=1)
resolved = motifs_metadata[['cluster', 'cluster_name']].drop_duplicates()
resolved['dup_index'] = resolved.groupby('cluster_name').cumcount()
resolved['unique_cluster_name'] = resolved.apply(lambda row: f"{row['cluster_name']}  #{row['dup_index'] + 1}", axis=1)
motifs_meta = motifs_metadata.reset_index().merge(resolved, on=['cluster', 'cluster_name'])

# Color by Sasha, might not relevant since we're running for >16 components
import matplotlib.colors as mcolors

def get_colors_order(n_components):
    # code from Wouter's library, to have the similar colors as in barplot
    colors_order = ['#FFE500', '#FE8102', '#FF0000', '#07AF00', '#4C7D14', '#414613', '#05C1D9', '#0467FD', '#009588', '#BB2DD4', '#7A00FF', '#4A6876', '#08245B', '#B9461D', '#692108', '#C3C3C3']
    neworder = np.array([16,10,7,11,2,12,1,8,4,15,14,5,9,6,3,13]).astype(int) - 1
    colors_order = list(np.array(colors_order)[neworder])

    if (n_components>len(colors_order)):
        colornames = np.sort(list(mcolors.CSS4_COLORS.keys()))
        count = len(colors_order)
        np.random.seed(10)
        myrandint = np.random.randint(len(colornames))
        while (count < n_components):
            myrandint =    np.random.randint(len(colornames))
            newcolor = colornames[myrandint]
            trialcount = 0
            while ((newcolor in colors_order) and (trialcount < 100)):
                newcolor = colornames[np.random.randint(0,len(colornames))]
                trialcount+=1
            #print('new color ',count,newcolor)
            colors_order.append(newcolor)
            count+=1
    return colors_order

colors_order = get_colors_order(16) # Colors by order

# Function to aggregate p-value
def cauchy_combination(p_val):
    return stats.cauchy.sf(np.sum(np.tan(np.multiply(np.subtract(0.5, p_val), np.pi))) / len(p_val))

# Function to generate plot
def generate_top(modal, cluster, component, **kwargs):
    X_components = "X" + str(component)
    vis = modal[modal["Components"] == X_components].head(10).sort_values(by="Estimation")
    # vis.plot(x="tf_name", y="Estimation", kind="bar")
    plt.barh(data=vis, y=cluster, width='Estimation', color = colors_order[component - 1] ,  **kwargs)
    plt.title(X_components)
    # plt.savefig('plot.png', dpi=300, bbox_inches='tight')
    plt.show()

# combine coeffs with some columns with motifs_metadata
new_df = pd.merge(coeffs,
                  motifs_metadata[['motif_id','tf_name']],
                  on='motif_id', how='left')


# new_df["q_val"] = multipletests(coeffs['Pr(>|z|)'], alpha=0.05, method='fdr_by', maxiter=1)[1]

# Dropping Intercept before cauchy?
# new_df.drop(columns=['(Intercept)'])

# Cauchy aggregation
qval_df = new_df.pivot_table(index='tf_name', columns='Unnamed: 0',
                             values='Pr(>|z|)', aggfunc=cauchy_combination)

# Do fdr corrected on p-val, fdr_by
qval_df_fdr = qval_df.applymap(lambda x: multipletests(x, alpha=0.05, method='fdr_by', maxiter=1)[1][0])

# Average estimate based on tf_name
tf_name_est_df = new_df.pivot_table(index='tf_name', columns='Unnamed: 0',
                             values='estimate', aggfunc=np.mean)

# Plotting
pvaltest_df = pd.DataFrame({'index': np.argmax(tf_name_est_df.to_numpy(), axis=1),
                            'value': tf_name_est_df.max(axis=1)})

max_df = pvaltest_df.sort_values(by=['index', 'value'], ascending=[True, False])

est = tf_name_est_df.loc[max_df.index, :]
qvals = qval_df_fdr.loc[max_df.index, :]

# fdr_by on P-value, then aggregate q-val by tf_name, plotting estimation by tf_name
filtered_est = est.where(~((est<0) | (qvals>0.05)))
filtering = filtered_est.reset_index()

tf_result_cluster_name = filtering.melt('tf_name',
                                     var_name='Components',
                                     value_name='Estimation').sort_values('Estimation', ascending=False)

# generate_top(modal, cluster, component)
for i in range(1,17):
    generate_top(tf_result_cluster_name, "tf_name", i)