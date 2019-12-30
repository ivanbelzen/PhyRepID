# IvB 
# Export results table of PhyRepID pipeline
# Calculate PRD score

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np

#Input
clans_dict = pre.get_pfam_clans(pre.pfam_clans_file)

#repeat stats
pfam_repeat_stats_file = pre.root+'pfam_repeat_stats.json'
meme_repeat_stats_file = pre.root+'meme_repeat_stats.json'

#evo events
pfam_results_file = pre.root+'pfam_evo_events.json'
meme_results_file = pre.root+'meme_evo_events.json'
pfam_gt_file = pre.root+"pfam_evo_events_genetrees.json"
meme_gt_file = pre.root+"meme_evo_events_genetrees.json"

selectome_file = pre.root+'selectome.tsv'
human_mapping_file = pre.root+'og_human_mapping.json' #Generate OG to human gene mapping from generate_human_protein_mapping.py
human_gene_symbol_file = pre.root+'human_gene_symbol.tsv'


#schaper_summary_df_file = pre.root+'schaper_summary_df_reroot.json'
#exac_df_file = pre.root+'analysis/exac_10102018/exac_df.json'



## Make dataframe evo events summary (repeat tree reconcilition)

with open(pfam_repeat_stats_file, 'r') as pfam_output:
	pfam_repeat_stats = json.load(pfam_output)
	
with open(meme_repeat_stats_file, 'r') as meme_output:
	meme_repeat_stats = json.load(meme_output)

with open(pfam_results_file, 'r') as results:
	pfam_results = json.load(results)

with open(meme_results_file, 'r') as results:
	meme_results = json.load(results)

pfam_summary_df_cols = ['og_id_domain', 'og_id', 'pre_model', 'clan', 'netto_dup', 'loss', 'orthologs_cnt']
pfam_summary_df_rows = []

meme_summary_df_cols = ['og_id', 'netto_dup', 'loss', 'orthologs_cnt']
meme_summary_df_rows = []


for og_id_domain,data in pfam_results['repeat'].items(): 
	og_id = '_'.join(og_id_domain.split('_')[0:2])
	pre_model = '_'.join(og_id_domain.split('_')[2:])
	clan = clans_dict[pre_model] if pre_model in clans_dict else ''
	netto_dup = data['netto_dup']
	loss = data['loss']

	if og_id_domain not in pfam_repeat_stats: continue
	orthologs_list = pfam_repeat_stats[og_id_domain][pre_model]['orthologs_dict'].keys()
	pfam_summary_df_rows.append( [og_id_domain, og_id, pre_model, clan, netto_dup, loss, len(orthologs_list)] )
	
pfam_summary_df = pd.DataFrame(pfam_summary_df_rows, columns=pfam_summary_df_cols)

for og_id,data in meme_results['repeat'].items(): 
	netto_dup = data['netto_dup']
	loss = data['loss']
	
	if og_id not in meme_repeat_stats: continue
	orthologs_list = meme_repeat_stats[og_id]['motif']['orthologs_dict'].keys()
	meme_summary_df_rows.append( [og_id, netto_dup, loss, len(orthologs_list)] )
	
meme_summary_df = pd.DataFrame(meme_summary_df_rows, columns=meme_summary_df_cols)

##


## Make combined summary dataframe 

combined_summary_df = pd.concat([pfam_summary_df,meme_summary_df], ignore_index=True, sort=False)
#Add clan for meme
combined_summary_df['clan'] = np.where(combined_summary_df['clan'].isna(), 'motif', combined_summary_df['clan'])

#Add identifier with has both meme and pfam annotation
combined_summary_df['identifier'] = np.where(combined_summary_df['og_id_domain'].isna(),combined_summary_df['og_id'],combined_summary_df['og_id_domain'])

combined_summary_df['genetree_id'] = combined_summary_df['identifier'].str.split('_').str[0]
combined_summary_df['gene_id'] = combined_summary_df['identifier'].str.split('_').str[1]	
combined_summary_df['og_id'] = combined_summary_df['genetree_id'] +'_'+ combined_summary_df['gene_id']

#Add gene symbol
with open(human_gene_symbol_file, 'r') as gene_data:
	mapping_gene_symbol = {}
	for line in gene_data:
		cols = line.strip().split('\t')
		mapping_gene_symbol[cols[0]]=cols[1]
		
combined_summary_df['gene_symbol'] = combined_summary_df['gene_id'].map(mapping_gene_symbol)

# Calculate PRD score 
mean_netto_dup = combined_summary_df['netto_dup'].mean()
combined_summary_df['PRD_score']= (combined_summary_df['netto_dup']-mean_netto_dup)/combined_summary_df['orthologs_cnt']

##Export

#with open(combined_summary_df_file, 'w') as output:
#	output.write(json.dumps(combined_summary_df.to_dict()))	

combined_summary_df[['gene_symbol','PRD_score','identifier','netto_dup','loss','orthologs_cnt','clan']]\
.to_csv(pre.root+'phyrepid_results_simple.csv')#, separator='\t')	

#Still need merging with other files


## Selectome

selectome = []
with open(selectome_file, 'r') as human_prot_selected:
	for line in human_prot_selected:
		cols = line.split('\t')
		selectome.append(cols[0])

with open(human_mapping_file , 'r') as human_proteins:
	og_human_mapping = json.load(human_proteins)

selectome_df_cols = ['identifier', 'protein_id', 'positive_selection']
selectome_df_rows = []

for og_id,human_list in og_human_mapping['meme'].items(): 
	for human_protein,orthologs_list in human_list.items():	 
		if human_protein in selectome: selectome_df_rows.append([og_id, human_protein,1])
		else: selectome_df_rows.append([og_id, human_protein,0])
		
for og_id,human_list in og_human_mapping['pfam'].items(): 
	for human_protein,orthologs_list in human_list.items():	 
		if human_protein in selectome: selectome_df_rows.append([og_id, human_protein,1])
		else: selectome_df_rows.append([og_id, human_protein,0])
		
selectome_df = pd.DataFrame(selectome_df_rows, columns=selectome_df_cols)

##

## Make dataframe evo events of gene tree reconciliation. 
with open(pfam_gt_file, 'r') as results:
	pfam_gt = json.load(results)

with open(meme_gt_file, 'r') as results:
	meme_gt = json.load(results)

genetree_df_cols = ['og_id', 'gt_dup', 'gt_loss','gt_root']
genetree_df_rows = []

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

##

quit()

#Make combined summary df with updated Duplication Score, unit CV and repeat unit duplication score

#Merge with genetree_df
combined_summary_df = combined_summary_df.merge(genetree_df, on='og_id')

#unit count 
with open(pre.root+'analysis/dataset_04102018/'+'unit_cnt_df.json', 'r') as output:
	unit_cnt = json.load(output)
	unit_cnt_df = pd.DataFrame.from_dict(unit_cnt)

combined_summary_df = combined_summary_df.merge(unit_cnt_df, on='identifier', how='left')
combined_summary_df['unit_cv'] = combined_summary_df['cv']
combined_summary_df['unit_sum'] = combined_summary_df['sum']

combined_summary_df['unit_scaled_dup']  = combined_summary_df['netto_dup'] / combined_summary_df['unit_sum']
unit_dup_mean = combined_summary_df['unit_scaled_dup'].mean()
combined_summary_df['unit_duplication_score']= combined_summary_df['unit_scaled_dup']-unit_dup_mean

combined_summary_df[['genetree_id','gene_symbol','duplication_score','normalized_score','expectancy_score','unit_duplication_score','identifier','netto_dup','loss','orthologs_cnt','clan','unit_cv','unit_sum','gt_dup','gt_root']]\
.to_csv(analysis_evo_events_path+'evo_events.csv')


##Combine summary df with other datasets

with open(schaper_summary_df_file, 'r') as output:
	schaper_summ = json.load(output)
	schaper_summary_df = pd.DataFrame.from_dict(schaper_summ)

with open(exac_df_file, 'r') as output:
	exac = json.load(output)
	exac_df = pd.DataFrame.from_dict(exac)
			
#add schaper results
combined_summary_df = combined_summary_df.merge(schaper_summary_df[['identifier','schaper_pos']], on='identifier', how='left')

# dnds selectome
#selectome_df['identifier'] = np.where(selectome_df['og_id_domain'].empty,selectome_df['og_id'],selectome_df['og_id_domain'])
#add dnds selectome

selectome_positive_set = selectome_df.loc[selectome_df['positive_selection'] == 1]
combined_summary_df['selectome'] = np.where(combined_summary_df['identifier'].isin(selectome_positive_set['identifier']),1,0)

#exac
combined_summary_df = combined_summary_df.merge(exac_df[['gene_id','mis_z','pLI']], on='gene_id', how='left')

#export
combined_summary_df[['genetree_id','gene_symbol','duplication_score','identifier','netto_dup','loss','orthologs_cnt','clan','unit_cv','selectome','schaper_pos','mis_z','pLI','gt_dup','gt_root']]\
.to_csv(analysis_evo_events_path+'evo_events_merged.csv')


## Text
std=combined_summary_df.describe()['duplication_score']['std']

print('#OGs in total',len(combined_summary_df.index))
print('#OGs >> netto dup - mean / orth cnt ', \
len(combined_summary_df.loc[combined_summary_df['duplication_score']>0].index), \
'>1stdev above mean',\
len(combined_summary_df.loc[combined_summary_df['duplication_score']>std].index))
