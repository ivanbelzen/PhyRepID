# IvB 
## Export results table of PhyRepID pipeline
# Calculate PRD score from evo events reconciliation of repeat tree with gene tree
# Incorporates additional info from 
# - Reconciliation of gene tree with species tree
# - Selectome
# - ExAC
# - Comparison to Schaper et al 
## Output files: 
# - phyrepid_results_simple
# - phyrepid_results_full

#Unstage 
#What outputs?
#(take out unit_cv only additional analysis but not required for table)
# Needs output of analyse schaper comparison if available.


import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np

#Output files
phyrepid_results_simple=pre.phyrepid_results_simple
phyrepid_results_full=pre.phyrepid_results_full

#Input
clans_dict = pre.get_pfam_clans(pre.pfam_clans_file)

#repeat stats
pfam_repeat_stats_file = pre.pfam_repeat_stats_file
meme_repeat_stats_file = pre.meme_repeat_stats_file

#evo events
pfam_results_file = pre.pfam_evo_events_file
meme_results_file = pre.meme_evo_events_file
pfam_gt_file = pre.pfam_evo_gt_file
meme_gt_file = pre.meme_evo_gt_file

# Additional datasets - fils not required 
selectome_file = pre.selectome_file #needs human mapping 
human_mapping_file = pre.human_mapping_file #OG to human gene mapping from generate_human_protein_mapping.py
exac_data_file = pre.exac_file #needs gene symbol
human_gene_symbol_file = pre.human_gene_symbol_file
schaper_comparison_file = pre.schaper_comparison_file #output from analyse_schaper_comparison.py

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

#Add gene symbol if available (necessary for ExAC)
combined_summary_df['gene_symbol'] = ''
if pre.file_notempty(human_gene_symbol_file):
	with open(human_gene_symbol_file, 'r') as gene_data:
		mapping_gene_symbol = {}
		for line in gene_data:
			cols = line.strip().split('\t')
			mapping_gene_symbol[cols[0]]=cols[1]

	combined_summary_df['gene_symbol'] = combined_summary_df['gene_id'].map(mapping_gene_symbol)

# Calculate PRD score 
mean_netto_dup = combined_summary_df['netto_dup'].mean()
combined_summary_df['PRD_score']= (combined_summary_df['netto_dup']-mean_netto_dup)/combined_summary_df['orthologs_cnt']

# Export

combined_summary_df.to_csv(phyrepid_results_simple, sep='\t')	

print('#OGs in total ',len(combined_summary_df.index))
print('#OGs >0 PRD score ', len(combined_summary_df.loc[combined_summary_df['PRD_score']>0].index))


### Merging with additional data if available

## Evo events of gene tree reconciliation

if pre.file_notempty(pfam_gt_file) and pre.file_notempty(meme_gt_file):		
	
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
	combined_summary_df = combined_summary_df.merge(genetree_df, on='og_id')

##

## Selectome
#list identifiers under positive selection
selectome_identifiers = [] 

if not pre.file_notempty(human_mapping_file):
	print("Cannot merge without human mapping file")
	quit()
	
if pre.file_notempty(selectome_file):
	selectome = []
	with open(selectome_file, 'r') as human_prot_selected:
		for line in human_prot_selected:
			cols = line.split('\t')
			selectome.append(cols[0])

	with open(human_mapping_file , 'r') as human_proteins:
		og_human_mapping = json.load(human_proteins)

	for og_id,human_list in og_human_mapping['meme'].items(): 
		for human_protein,orthologs_list in human_list.items():	 
			if human_protein in selectome: selectome_identifiers.append(og_id)
				
	for og_id,human_list in og_human_mapping['pfam'].items(): 
		for human_protein,orthologs_list in human_list.items():	 
			if human_protein in selectome: selectome_identifiers.append(og_id)
		
	combined_summary_df['selectome'] = np.where(combined_summary_df['identifier'].isin(selectome_identifiers),1,0)
		
##

## ExAC 
if pre.file_notempty(exac_data_file):
		
	exac_cols = ['gene_symbol', 'syn_z','mis_z','pLI']
	exac_rows = []

	with open(exac_data_file, 'r') as exac_data:
		exac = {}
		for line in exac_data:
			cols = line.strip().split('\t')
			if cols[1] == 'gene': continue
			exac_rows.append( [ cols[1], float(cols[16]), float(cols[17]), float(cols[19]) ])
			
	exac_df = pd.DataFrame(exac_rows, columns=exac_cols)
	exac_df['gene_id'] = exac_df['gene_symbol'].map(mapping_gene_symbol)

	combined_summary_df = combined_summary_df.merge(exac_df[['gene_id','mis_z','pLI']], on='gene_id', how='left')

##

## Comparison to Schaper et al results
if pre.file_notempty(schaper_comparison_file):		
	with open(schaper_comparison_file, 'r') as output:
		schaper_summ = json.load(output)
		schaper_summary_df = pd.DataFrame.from_dict(schaper_summ)
					
	combined_summary_df = combined_summary_df.merge(schaper_summary_df[['identifier','schaper_positive']], on='identifier', how='left')


#Export
combined_summary_df.to_csv(phyrepid_results_full, sep='\t')	


#ToDo: re-instate unit count

'''
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

combined_summary_df.to_csv()
'''
