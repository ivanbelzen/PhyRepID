# IvB 

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np
import seaborn as sns
from scipy.stats.mstats import mode, gmean, hmean
import matplotlib.ticker as ticker

#Combined summary df made in analyse_merge_export.py
analysis_path = pre.root+'analysis/evo_events_14102018/'
sys.stdout = open(analysis_path+'genetrees.txt', 'w')

#File names for loading/saving dataframes
combined_summary_df_file = pre.root+'analysis/evo_events_06102018/'+'pre_dataset_df.json'
genetree_df_file = pre.root+'genetree_df_reroot.json'

#Make gene tree dataframe for pfam and meme
#'''
genetree_df_cols = ['og_id', 'gt_dup', 'gt_loss','gt_root']
genetree_df_rows = []

pfam_gt_file = pre.root+'pfam_evo_events_reroot_gt.json'
meme_gt_file = pre.root+'meme_evo_events_reroot_gt.json'

with open(pfam_gt_file, 'r') as results:
	pfam_gt = json.load(results)

with open(meme_gt_file, 'r') as results:
	meme_gt = json.load(results)

for og_id,data in pfam_gt['genetree'].items(): 
	netto_dup = data['netto_dup']
	loss = data['loss']
	root = data['root']
	genetree_df_rows.append([og_id,netto_dup,loss,root])
for og_id,data in meme_gt['genetree'].items(): 
	netto_dup = data['netto_dup']
	loss = data['loss']
	root = data['root']
	genetree_df_rows.append([og_id,netto_dup,loss,root])
		
genetree_df = pd.DataFrame(genetree_df_rows, columns=genetree_df_cols)

with open(genetree_df_file, 'w') as output:
	output.write(json.dumps(genetree_df.to_dict()))	
#'''
## Analyse genetree dataframe
with open(genetree_df_file,'r') as output:
	gt = json.load(output)
	genetree_df = pd.DataFrame.from_dict(gt)

with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)

stats = combined_summary_df.describe()
std=stats['duplication_score']['std']

combined_summary_df['genetree_id'] = combined_summary_df['identifier'].str.split('_').str[0]	
combined_summary_df['gene_id'] = combined_summary_df['identifier'].str.split('_').str[1]	
combined_summary_df['og_id'] = combined_summary_df['genetree_id'] +'_'+ combined_summary_df['gene_id']
print(combined_summary_df['og_id'].head())
merged_df = combined_summary_df.merge(genetree_df, on='og_id')

print("Analysis gene trees")
print("Entries",len(genetree_df.index), "unique OGs",genetree_df.og_id.nunique())
print("OGs with gene duplication",genetree_df.loc[genetree_df['gt_dup']>0].og_id.nunique())
print("Entries overlap with combined summary",len(merged_df.index), "unique OGs",merged_df.og_id.nunique())
print("OGs with gene duplication",merged_df.loc[merged_df['gt_dup']>0].og_id.nunique())

positive_set = merged_df.loc[merged_df['duplication_score']>0]
fast_set =merged_df.loc[merged_df['duplication_score']>std]

print("OGs with gene duplication in positive set",positive_set.loc[positive_set['gt_dup']>0].og_id.nunique())
print("OGs with gene duplication in fast set",fast_set.loc[fast_set['gt_dup']>0].og_id.nunique())

print('Fraction OGs with gene dup total',float(len(merged_df.loc[merged_df['gt_dup']>0].index))/len(merged_df.index))
print('Fraction OGs with gene dup positive',float(len(positive_set.loc[positive_set['gt_dup']>0].index))/len(positive_set.index))
print('Fraction OGs with gene dup fast',float(len(fast_set.loc[fast_set['gt_dup']>0].index))/len(fast_set.index))

print("OGs with non-euteleostomi root",len(merged_df.loc[merged_df['gt_root']!='t117571'].index))
print(float(len(merged_df.loc[merged_df['gt_root']!='t117571'].index))/len(merged_df.index))
print("OGs with non-euteleostomi root in fast set",len(fast_set.loc[fast_set['gt_root']!='t117571'].index))
print(float(len(fast_set.loc[fast_set['gt_root']!='t117571'].index))/len(fast_set.index))
print("Average gene duplications total vs fast",merged_df['gt_dup'].mean(),fast_set['gt_dup'].mean())
print("Average gene duplications non-euteleostomi root total vs fast",merged_df.loc[merged_df['gt_root']!='t117571']['gt_dup'].mean(),fast_set.loc[fast_set['gt_root']!='t117571']['gt_dup'].mean())
