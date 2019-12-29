# IvB 
# September 2018
# Input: pfam_evo_events.json meme_evo_events.json
# PRE evo events summary contains dup/loss per node/ortholog in an OG and in total 
# Need to parse the taxons in human lineage

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np
import seaborn as sns
from scipy.stats.mstats import mode, gmean, hmean


## Settings
#taken from analyse_maptotree.py
node_projection = { 't32525': ['t9347', 't1437010', 't314147'], 't8287': ['t32523'], 't186625':['t186626'], 't41665': ['t1489872'], 't117571':['t7742', 't7711', 't33213']}
root_nodes = ['t117571','t7742', 't7711', 't33213'] #bony vertebrates and before

#File names for loading/saving dataframes
#input
human_protein_file = 'og_human_mapping.json'
pfam_results_file = 'pfam_evo_events_reroot.json'
meme_results_file = 'meme_evo_events_reroot.json'
clans_dict = pre.get_pfam_clans(pre.pfam_clans_file)
#clans_dict = pre.get_pfam_clans('Pfam-A.clans.tsv')

#output
human_summary_df_file = 'human_summary_df.json'

#File names for comparison purposes
selectome_df_file = 'selectome_df.json'
combined_summary_df_file = 'combined_summary_df.json'


## Analysis specific settings
#Full lineage
analysis_path_full_lineage = pre.root+'analysis/human_evo_events/full_lineage_trace/'
taxa_full_lineage = ['ENSP','t9606','t314146','t32525','t40674','t32524','t8287']
taxa_full_lineage.extend(node_projection['t32525'])
taxa_full_lineage.extend(node_projection['t8287'])

#human_only
analysis_path_human_only = pre.root+'analysis/human_evo_events/human_only/'
taxa_human_only = ['ENSP','t9606']

## Figure definitions
def make_hist(data,output_path='',rug=False,log_scale=False,bin_cnt=100):
	weights = np.ones_like(np.array(data))/float(len(np.array(data)))
	plt.hist(data,weights=weights,color=pre.blue_colour,bins=bin_cnt)
	if rug: sns.rugplot(data,color=pre.green_colour, height=0.01)
	if log_scale: plt.xscale(log_scale)
	if output_path != '': 
		plt.savefig(output_path,format='png',bbox_inches='tight')
		plt.clf()
	else:
		return(plt)
	

def make_human_summary_df(data,taxa,analysis_output_path=''): #no output files if this is left empty
	skipped_ogs = [] #some not in human og mapping?? 
	skipped_ogs_file = 'excluded_ogs_nomapping.txt'
	
	human_summary_df_cols = ['identifier', 'human_protein_cnt', 'domain', 'clan', 'netto_dup','root_dup']
	human_summary_df_rows = []
	human_dataframe_file = 'human_summary_df.json'
	
	for og_id_domain,data in data.items(): 
		og_id = '_'.join(og_id_domain.split('_')[0:2])
		pre_model = og_id_domain[len(og_id)+1:]
		if pre_model != '':
			clan = clans_dict[pre_model] if pre_model in clans_dict else ''
			data_type = 'pfam'
		else: 
			clan = 'motif'
			data_type = 'meme'
		
		if og_id_domain in og_human_mapping[data_type]: 
			human_list = og_human_mapping[data_type][og_id_domain]
		else: 
			skipped_ogs.append(og_id_domain)
			continue
		
		root_dup = 0
		netto_dup = 0
		
		for node,cnt in data['events']['dup'].items():
			node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))
				
			if node in root_nodes: root_dup += cnt
			elif node in taxa: netto_dup += cnt
		
		if len(human_list) > 0:
			human_summary_df_rows.append( [og_id_domain, len(human_list), pre_model, clan, netto_dup,root_dup])
	
	human_summary_df = pd.DataFrame(human_summary_df_rows, columns=human_summary_df_cols)
	
	
	if analysis_output_path != '': 
		with open(analysis_output_path + skipped_ogs_file ,'w') as output:
			for x in skipped_ogs:
				output.write(x+'\n')
	
		with open(analysis_output_path + human_dataframe_file, 'w') as output:
			output.write(json.dumps(human_summary_df.to_dict()))	
		
	return(human_summary_df)


## Load data if needed in multiple analyses
with open(pfam_results_file, 'r') as results:
	pfam_results = json.load(results)
with open(meme_results_file, 'r') as results:
	meme_results = json.load(results)
with open(human_protein_file, 'r') as human_proteins:
	og_human_mapping = json.load(human_proteins)
	
repeat_data = pfam_results['repeat']
repeat_data.update(meme_results['repeat'])


with open(combined_summary_df_file, 'r') as output:
	combined_summary_data = json.load(output)
	

## Analysis: full lineage

sys.stdout = open(analysis_path_full_lineage+'output.txt', 'w')

#Make full lineage dataframe
if pre.file_notempty(analysis_path_full_lineage+human_summary_df_file):
	with open(analysis_path_full_lineage+human_summary_df_file, 'r') as output:
		human_summary = json.load(output)
		human_summary_df = pd.DataFrame.from_dict(human_summary)
else:
	print('nothing')
	#human_summary_df = make_human_summary_df(repeat_data,taxa_full_lineage,analysis_path_full_lineage)

#Print some stats

print('#OGs in total',human_summary_df.identifier.nunique(),'(same as index?',len(human_summary_df.index),') human proteins',human_summary_df.human_protein_cnt.sum())

#calculate dev from mean as netto dup / human cnt - mean
mean_netto_dup = (human_summary_df['netto_dup']/human_summary_df['human_protein_cnt']).mean()
human_summary_df['dev_from_mean'] =(human_summary_df['netto_dup']/human_summary_df['human_protein_cnt'])-mean_netto_dup
print(human_summary_df.describe().to_string())

std=human_summary_df.describe()['dev_from_mean']['std']

print('#OGs >> netto dup / human cnt - mean ', len(human_summary_df.loc[human_summary_df['dev_from_mean']>0].index))
print('#OGs >> netto dup > 0 ', len(human_summary_df.loc[human_summary_df['netto_dup']>0].index))
print('#OGs >> netto dup > 1 ', len(human_summary_df.loc[human_summary_df['netto_dup']>1].index))
print('#OGs >> dev_from_mean > std',len(human_summary_df.loc[human_summary_df['dev_from_mean']>std].index))

mean_netto_dup = human_summary_df['netto_dup'].mean()
human_summary_df['dev_from_mean'] =(human_summary_df['netto_dup']-mean_netto_dup)/human_summary_df['human_protein_cnt']
print(human_summary_df.describe().to_string())

print('#OGs >> netto dup - mean / human cnt ', len(human_summary_df.loc[human_summary_df['dev_from_mean']>0].index))
print('#OGs >> netto dup > 0 ', len(human_summary_df.loc[human_summary_df['netto_dup']>0].index))
print('#OGs >> netto dup > 1 ', len(human_summary_df.loc[human_summary_df['netto_dup']>1].index))
print('#OGs >> dev_from_mean > std',len(human_summary_df.loc[human_summary_df['dev_from_mean']>std].index))

human_summary_df['human_dup'] = human_summary_df['netto_dup']
human_summary_df['human_duplication_score'] = human_summary_df['dev_from_mean']


quit()

#Histogram of netto duplications
make_hist(human_summary_df['netto_dup'],output_path=analysis_path_full_lineage+'netto_dup.png',rug=True,log_scale=False,bin_cnt=100)
#Histogram of duplication score
make_hist(human_summary_df['dev_from_mean'],output_path=analysis_path_full_lineage+'dup_score.png',rug=True,log_scale=False,bin_cnt=200)

#Comparison with full dataset 
human_summary_df['human_dup'] = human_summary_df['netto_dup']
combined_summary_df = pd.DataFrame.from_dict(combined_summary_data)
combined_summary_df = combined_summary_df.merge(human_summary_df[['identifier','human_dup','dev_from_mean','human_protein_cnt']], on='identifier')
combined_summary_df['human_frac']=np.where(combined_summary_df['netto_dup'] == 0, 0, combined_summary_df['human_dup']/combined_summary_df['netto_dup'])

#Save with human frac
new_human_summary_df = human_summary_df.merge(combined_summary_df[['identifier','human_frac']], on='identifier')
with open(analysis_path_full_lineage + 'human_summary_df.json', 'w') as output:
	output.write(json.dumps(new_human_summary_df.to_dict()))	
quit()

#Scatterplot of full dataset duplication score vs human dev from mean
sns.scatterplot(data=combined_summary_df,x='duplication_score',y='dev_from_mean')
plt.savefig(analysis_path_full_lineage+'comparision_full_dataset.png',format='png',bbox_inches='tight')
plt.clf()



combined_summary_df[['gene_symbol','identifier','clan','duplication_score','netto_dup','orthologs_cnt','human_dup','human_frac','dev_from_mean','human_protein_cnt']].to_csv(analysis_path_full_lineage+'human_full_lineage.csv')

## Analysis: human only
sys.stdout = open(analysis_path_human_only+'output.txt', 'w')


#Make full lineage dataframe
if pre.file_notempty(analysis_path_human_only+human_summary_df_file):
	with open(analysis_path_human_only+human_summary_df_file, 'r') as output:
		human_summary = json.load(output)
		human_summary_df = pd.DataFrame.from_dict(human_summary)
else:
	human_summary_df = make_human_summary_df(repeat_data,taxa_human_only,analysis_path_human_only)

#Print some stats

print('#OGs in total',human_summary_df.identifier.nunique(),'(same as index?',len(human_summary_df.index),') human proteins',human_summary_df.human_protein_cnt.sum())

#calculate dev from mean as netto dup / human cnt - mean
mean_netto_dup = (human_summary_df['netto_dup']/human_summary_df['human_protein_cnt']).mean()
human_summary_df['dev_from_mean'] =(human_summary_df['netto_dup']/human_summary_df['human_protein_cnt'])-mean_netto_dup
print(human_summary_df.describe())

std=human_summary_df.describe()['dev_from_mean']['std']

print('#OGs >> netto dup / human cnt - mean ', len(human_summary_df.loc[human_summary_df['dev_from_mean']>0].index))
print('#OGs >> netto dup > 0 ', len(human_summary_df.loc[human_summary_df['netto_dup']>0].index))
print('#OGs >> netto dup > 1 ', len(human_summary_df.loc[human_summary_df['netto_dup']>1].index))
print('#OGs >> dev_from_mean > std',len(human_summary_df.loc[human_summary_df['dev_from_mean']>std].index))


#Histogram of netto duplications
make_hist(human_summary_df['netto_dup'],output_path=analysis_path_human_only+'netto_dup.png',rug=True,log_scale=False,bin_cnt=20)

#Histogram of duplication score
sns.distplot(human_summary_df['dev_from_mean'], kde_kws={'cut':0,'bw':0.1})
plt.savefig(analysis_path_human_only+'dup_score.png',format='png',bbox_inches='tight')
plt.clf()

#Comparison with full dataset 
human_summary_df['human_dup'] = human_summary_df['netto_dup']
combined_summary_df = pd.DataFrame.from_dict(combined_summary_data)
combined_summary_df = combined_summary_df.merge(human_summary_df[['identifier','human_dup','dev_from_mean','human_protein_cnt']], on='identifier')
combined_summary_df['human_frac']=np.where(combined_summary_df['netto_dup'] == 0, 0, combined_summary_df['human_dup']/combined_summary_df['netto_dup'])

#Scatterplot of full dataset duplication score vs human dev from mean
sns.scatterplot(data=combined_summary_df,x='duplication_score',y='dev_from_mean')
plt.savefig(analysis_path_human_only+'comparision_full_dataset.png',format='png',bbox_inches='tight')
plt.clf()

combined_summary_df[['gene_symbol','identifier','clan','duplication_score','netto_dup','orthologs_cnt','human_dup','human_frac','dev_from_mean','human_protein_cnt']].to_csv(analysis_path_human_only+'human_only.csv')

