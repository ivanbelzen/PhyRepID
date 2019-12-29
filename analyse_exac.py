# IvB 
# Edit 03-09-2018 
# Input: exac data file 'fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt' 
# documentation: ftp://ftp.broadinstitute.org/pub/ExAC_release/current/functional_gene_constraint/README_fordist_cleaned_exac_r03_z_data_pLI_2016_01_13.txt
# Input: human_gene_symbol.tsv mapping of human ensembl genes to gene symbol 
# Input: pfam_summary_df.json meme_summary_df.json previously made by analyse_evo_events.py

# Exac df is defined per human gene, hence per OG

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np
import seaborn as sns
from scipy.stats.mstats import mode, gmean, hmean
import matplotlib.ticker as ticker

#analysis_path = pre.root+'analysis/exac/'
#analysis_path_exac = pre.root+'analysis/exac_10102018/'
analysis_path_exac = pre.root+'analysis/exac_01112018/'

#File names for loading/saving dataframes
exac_data_file = pre.root+'analysis/exac/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt'
exac_df_file = 'exac_df.json'
human_gene_symbol_file = pre.root+'human_gene_symbol.tsv'

human_df_file = pre.root+'analysis/human_evo_events/full_lineage_trace/human_summary_df.json'
combined_summary_df_file = pre.root+'analysis/evo_events_06102018/pre_dataset_df.json'

##Make EXAC dataframe
'''
exac_cols = ['gene_symbol', 'syn_z','mis_z','pLI']
exac_rows = []

with open(human_gene_symbol_file, 'r') as gene_data:
	mapping_gene_symbol = {}
	for line in gene_data:
		cols = line.strip().split('\t')
		mapping_gene_symbol[cols[1]]=cols[0]

with open(exac_data_file, 'r') as exac_data:
	exac = {}
	for line in exac_data:
		cols = line.strip().split('\t')
		if cols[1] == 'gene': continue
		exac_rows.append( [ cols[1], float(cols[16]), float(cols[17]), float(cols[19]) ])
		
exac_df = pd.DataFrame(exac_rows, columns=exac_cols)
exac_df['gene_id'] = exac_df['gene_symbol'].map(mapping_gene_symbol)

with open(analysis_path_exac+exac_df_file, 'w') as output:
	output.write(json.dumps(exac_df.to_dict()))
'''	
## Analysis
	
#read exac df
with open(pre.root+'analysis/exac_10102018/'+exac_df_file, 'r') as output:
	exac = json.load(output)
	exac_df = pd.DataFrame.from_dict(exac)
#read datasets
with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)	
with open(human_df_file, 'r') as output:
	human = json.load(output)
	human_df = pd.DataFrame.from_dict(human)

#Full dataset
'''
sys.stdout = open(analysis_path_exac+'full_dataset.txt', 'w')


combined_summary_df['gene_id'] = combined_summary_df['identifier'].str.split('_').str[1]
exac_summary_df = pd.merge(combined_summary_df, exac_df, on=['gene_id'], how='inner')

print('#OGs in total',len(combined_summary_df.index), 'exac entries in total',len(exac_df.index), '#OGs in Exac', len(exac_summary_df.index))
print('genetrees',combined_summary_df.genetree_id.nunique(),exac_summary_df.genetree_id.nunique())

positive_set = combined_summary_df.loc[combined_summary_df['duplication_score']>0]
fast_set = combined_summary_df.loc[combined_summary_df['duplication_score']>combined_summary_df['duplication_score'].std()]
conserved_set = combined_summary_df.loc[combined_summary_df['duplication_score']<0]

exac_positive_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(positive_set['identifier'])]
exac_fast_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(fast_set['identifier'])]
exac_conserved_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(conserved_set['identifier'])]

print('#OGs in positive set and in Exac',len(exac_positive_set.index),exac_positive_set.genetree_id.nunique())
print('#OGs in fast set and in Exac',len(exac_fast_set.index),exac_fast_set.genetree_id.nunique())

#Interpretation of EXAC scores: from their FAQ
## For synonymous and missense, we created a signed Z score for the deviation of observed counts from the expected number. 
### Positive Z scores indicate increased constraint (intolerance to variation). Negative Z scores are given to genes that had a more variants than expected.
## For LoF, we assume that there are three classes of genes with respect to tolerance to LoF variation: null (tolerated), recessive (where heterozygous LoFs are tolerated), 
### and haploinsufficient (where heterozygous LoFs are not tolerated). Probability that a falls into the third category. 
### The closer pLI is to one, the more LoF intolerant the gene appears to be. We consider pLI >= 0.9 as an extremely LoF intolerant set of genes.


print('Enrichment analyses')
print("Increased constraint to missense mis_z > 0:",len(exac_summary_df.loc[exac_summary_df['mis_z']>0].index),'of total dataset',float(len(exac_summary_df.loc[exac_summary_df['mis_z']>0].index))/len(exac_summary_df.index))
print("positive set",len(exac_positive_set.loc[exac_summary_df['mis_z']>0].index),float(len(exac_positive_set.loc[exac_summary_df['mis_z']>0].index))/len(exac_positive_set.index))
print("fast set",len(exac_fast_set.loc[exac_summary_df['mis_z']>0].index),float(len(exac_fast_set.loc[exac_summary_df['mis_z']>0].index))/len(exac_fast_set.index))
print("conserved set",len(exac_conserved_set.loc[exac_summary_df['mis_z']>0].index),float(len(exac_conserved_set.loc[exac_summary_df['mis_z']>0].index))/len(exac_conserved_set.index))

print("LoF intolerant pLI > 0.9:",len(exac_summary_df.loc[exac_summary_df['pLI']>0.9].index),'total OGs','of total dataset',float(len(exac_summary_df.loc[exac_summary_df['pLI']>0.9].index))/len(exac_summary_df.index))
print("positive set",len(exac_positive_set.loc[exac_summary_df['pLI']>0.9].index),float(len(exac_positive_set.loc[exac_summary_df['pLI']>0.9].index))/len(exac_positive_set.index))
print("fast set",len(exac_fast_set.loc[exac_summary_df['pLI']>0.9].index),float(len(exac_fast_set.loc[exac_summary_df['pLI']>0.9].index))/len(exac_fast_set.index))
print("conserved set",len(exac_conserved_set.loc[exac_summary_df['pLI']>0.9].index),float(len(exac_conserved_set.loc[exac_summary_df['pLI']>0.9].index))/len(exac_conserved_set.index))



print('Enrichment analyses,but with gene trees')
print("Increased constraint to missense mis_z > 0:",exac_summary_df.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),'of total dataset',float(exac_summary_df.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_summary_df.genetree_id.nunique())
print("positive set",exac_positive_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),float(exac_positive_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_positive_set.genetree_id.nunique())
print("fast set",exac_fast_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),float(exac_fast_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_fast_set.genetree_id.nunique())
print("conserved set",exac_conserved_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),float(exac_conserved_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_conserved_set.genetree_id.nunique())

print("LoF intolerant pLI > 0.9:",exac_summary_df.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),'total OGs','of total dataset',float(exac_summary_df.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_summary_df.genetree_id.nunique())
print("positive set",exac_positive_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),float(exac_positive_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_positive_set.genetree_id.nunique())
print("fast set",exac_fast_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),float(exac_fast_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_fast_set.genetree_id.nunique())
print("conserved set",exac_conserved_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),float(exac_conserved_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_conserved_set.genetree_id.nunique())

sns.scatterplot(data=exac_summary_df, x='duplication_score',y='mis_z')
plt.savefig(analysis_path_exac+'mis_z_vs_dupscore.png',format='png',bbox_inches='tight')
plt.xscale('symlog')
#plt.show()
plt.clf()


sns.scatterplot(data=exac_summary_df, x='duplication_score',y='pLI')
plt.savefig(analysis_path_exac+'pLI_vs_dupscore.png',format='png',bbox_inches='tight')
plt.xscale('symlog')
#plt.show()
plt.clf()


#Plots

sns.distplot(exac_summary_df['syn_z'], color='blue',label='negative',norm_hist=True)
sns.distplot(exac_positive_set['syn_z'], color='green',label='positive',norm_hist=True)
plt.legend()
plt.savefig(analysis_path_exac+'syn_z.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

sns.distplot(exac_summary_df['mis_z'], color='blue',label='negative',norm_hist=True)
sns.distplot(exac_positive_set['mis_z'], color='green',label='positive',norm_hist=True)
plt.legend()
plt.savefig(analysis_path_exac+'mis_z.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

weights = np.ones_like(np.array(exac_summary_df['pLI']))/float(len(np.array(exac_summary_df['pLI'])))
plt.hist(exac_summary_df['pLI'],weights=weights,color='blue',label='negative',bins=20,alpha=0.6)
weights = np.ones_like(np.array(exac_positive_set['pLI']))/float(len(np.array(exac_positive_set['pLI'])))
plt.hist(exac_positive_set['pLI'],weights=weights,color='green',label='positive',bins=20,alpha=0.6)
#sns.rugplot(exac_summary_df['pLI'], color='blue')
#sns.rugplot(exac_positive_set['pLI'], color='green')
plt.legend()
plt.savefig(analysis_path_exac+'pli.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

'''

#Human

sys.stdout = open(analysis_path_exac+'human_dataset.txt', 'w')

#necessary for overlap identifiers
positive_set = combined_summary_df.loc[combined_summary_df['duplication_score']>0]
fast_set = combined_summary_df.loc[combined_summary_df['duplication_score']>combined_summary_df['duplication_score'].std()]


human_df['gene_id'] = human_df['identifier'].str.split('_').str[1]
human_df['genetree_id'] = human_df['identifier'].str.split('_').str[0]
exac_summary_df = pd.merge(human_df, exac_df, on=['gene_id'], how='inner')

print('#OGs in total',len(human_df.index), 'exac entries in total',len(exac_df.index), '#OGs in Exac', len(exac_summary_df.index))
print('genetrees',human_df.genetree_id.nunique(),exac_summary_df.genetree_id.nunique())


'''
positive_set = human_df.loc[human_df['dev_from_mean']>0]
fast_set = human_df.loc[human_df['dev_from_mean']>human_df['dev_from_mean'].std()]
conserved_set = human_df.loc[human_df['dev_from_mean']<0]

exac_positive_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(positive_set['identifier'])]
exac_fast_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(fast_set['identifier'])]
exac_conserved_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(conserved_set['identifier'])]
'''

human_positive_set = human_df.loc[(human_df['dev_from_mean']>0)&(human_df['identifier'].isin(positive_set['identifier']))]
human_fast_set = human_df.loc[(human_df['dev_from_mean']>0)&(human_df['identifier'].isin(fast_set['identifier']))]
human_conserved_set = human_df.loc[(human_df['dev_from_mean']<0)&(~human_df['identifier'].isin(positive_set['identifier']))]


exac_positive_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(human_positive_set['identifier'])]
exac_fast_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(human_fast_set['identifier'])]
exac_conserved_set = exac_summary_df.loc[exac_summary_df['identifier'].isin(human_conserved_set['identifier'])]

print('#OGs in positive set and in Exac',len(exac_positive_set.index),exac_positive_set.genetree_id.nunique())
print('#OGs in fast set and in Exac',len(exac_fast_set.index),exac_fast_set.genetree_id.nunique())

#Interpretation of EXAC scores: from their FAQ
## For synonymous and missense, we created a signed Z score for the deviation of observed counts from the expected number. 
### Positive Z scores indicate increased constraint (intolerance to variation). Negative Z scores are given to genes that had a more variants than expected.
## For LoF, we assume that there are three classes of genes with respect to tolerance to LoF variation: null (tolerated), recessive (where heterozygous LoFs are tolerated), 
### and haploinsufficient (where heterozygous LoFs are not tolerated). Probability that a falls into the third category. 
### The closer pLI is to one, the more LoF intolerant the gene appears to be. We consider pLI >= 0.9 as an extremely LoF intolerant set of genes.


print('Enrichment analyses')
print("Increased constraint to missense mis_z > 0:",len(exac_summary_df.loc[exac_summary_df['mis_z']>0].index),'of total dataset',float(len(exac_summary_df.loc[exac_summary_df['mis_z']>0].index))/len(exac_summary_df.index))
print("positive set",len(exac_positive_set.loc[exac_summary_df['mis_z']>0].index),float(len(exac_positive_set.loc[exac_summary_df['mis_z']>0].index))/len(exac_positive_set.index))
print("fast set",len(exac_fast_set.loc[exac_summary_df['mis_z']>0].index),float(len(exac_fast_set.loc[exac_summary_df['mis_z']>0].index))/len(exac_fast_set.index))
print("conserved set",len(exac_conserved_set.loc[exac_summary_df['mis_z']>0].index),float(len(exac_conserved_set.loc[exac_summary_df['mis_z']>0].index))/len(exac_conserved_set.index))

print("LoF intolerant pLI > 0.9:",len(exac_summary_df.loc[exac_summary_df['pLI']>0.9].index),'total OGs','of total dataset',float(len(exac_summary_df.loc[exac_summary_df['pLI']>0.9].index))/len(exac_summary_df.index))
print("positive set",len(exac_positive_set.loc[exac_summary_df['pLI']>0.9].index),float(len(exac_positive_set.loc[exac_summary_df['pLI']>0.9].index))/len(exac_positive_set.index))
print("fast set",len(exac_fast_set.loc[exac_summary_df['pLI']>0.9].index),float(len(exac_fast_set.loc[exac_summary_df['pLI']>0.9].index))/len(exac_fast_set.index))
print("conserved set",len(exac_conserved_set.loc[exac_summary_df['pLI']>0.9].index),float(len(exac_conserved_set.loc[exac_summary_df['pLI']>0.9].index))/len(exac_conserved_set.index))



print('Enrichment analyses,but with gene trees')
print("Increased constraint to missense mis_z > 0:",exac_summary_df.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),'of total dataset',float(exac_summary_df.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_summary_df.genetree_id.nunique())
print("positive set",exac_positive_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),float(exac_positive_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_positive_set.genetree_id.nunique())
print("fast set",exac_fast_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),float(exac_fast_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_fast_set.genetree_id.nunique())
print("conserved set",exac_conserved_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique(),float(exac_conserved_set.loc[exac_summary_df['mis_z']>0].genetree_id.nunique())/exac_conserved_set.genetree_id.nunique())

print("LoF intolerant pLI > 0.9:",exac_summary_df.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),'total OGs','of total dataset',float(exac_summary_df.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_summary_df.genetree_id.nunique())
print("positive set",exac_positive_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),float(exac_positive_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_positive_set.genetree_id.nunique())
print("fast set",exac_fast_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),float(exac_fast_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_fast_set.genetree_id.nunique())
print("conserved set",exac_conserved_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique(),float(exac_conserved_set.loc[exac_summary_df['pLI']>0.9].genetree_id.nunique())/exac_conserved_set.genetree_id.nunique())

sns.scatterplot(data=exac_summary_df, x='dev_from_mean',y='mis_z')
plt.savefig(analysis_path_exac+'human_mis_z_vs_dupscore.png',format='png',bbox_inches='tight')
plt.xscale('symlog')
#plt.show()
plt.clf()


sns.scatterplot(data=exac_summary_df, x='dev_from_mean',y='pLI')
plt.savefig(analysis_path_exac+'human_pLI_vs_dupscore.png',format='png',bbox_inches='tight')
plt.xscale('symlog')
#plt.show()
plt.clf()


#Plots
sns.distplot(exac_summary_df['syn_z'], color='blue',label='full',norm_hist=True)
sns.distplot(exac_positive_set['syn_z'], color='green',label='dynamic',norm_hist=True)
plt.legend()
plt.savefig(analysis_path_exac+'human_syn_z.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

sns.distplot(exac_summary_df['mis_z'], color='blue',label='full',norm_hist=True)
sns.distplot(exac_positive_set['mis_z'], color='green',label='dynamic',norm_hist=True)
plt.legend()
plt.savefig(analysis_path_exac+'human_mis_z.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

weights = np.ones_like(np.array(exac_summary_df['pLI']))/float(len(np.array(exac_summary_df['pLI'])))
plt.hist(exac_summary_df['pLI'],weights=weights,color='blue',label='full',bins=20,alpha=0.6)
weights = np.ones_like(np.array(exac_positive_set['pLI']))/float(len(np.array(exac_positive_set['pLI'])))
plt.hist(exac_positive_set['pLI'],weights=weights,color='green',label='dynamic',bins=20,alpha=0.6)
#sns.rugplot(exac_summary_df['pLI'], color='blue')
#sns.rugplot(exac_positive_set['pLI'], color='green')
plt.legend()
plt.savefig(analysis_path_exac+'human_pli.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
