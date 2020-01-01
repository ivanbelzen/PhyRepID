# IvB 
# September 2018
# Last edit: 01-01-2020
# Input: pfam_evo_events.json meme_evo_events.json
# PRE evo events summary contains dup/loss per node/ortholog in an OG and in total 
# Need to parse the taxons in human lineage
# Full-lineage analysis: consider duplications in amniote to human branch. For comparison to human datasets, use duplications in this lineage trace only (like require 1+)
# Human only analysis: consider duplications in the direct human branch - after split from MRCA with mouse. 

import sys,json
import pandas as pd
import pipeline_methods as pre
import numpy as np
from scipy.stats.mstats import mode, gmean, hmean

#Input files
human_mapping_file = pre.root+'og_human_mapping.json'
pfam_results_file = pre.root+'pfam_evo_events.json'
meme_results_file = pre.root+'meme_evo_events.json'
clans_dict = pre.get_pfam_clans(pre.pfam_clans_file)

#Output
output_full_lineage = pre.root+'phyrepid_results_human_full_lineage.tsv'
output_human_only = pre.root+'phyrepid_results_human_only.tsv'


## Taxon id settings

## Mappings
node_projection = { 't32525': ['t9347', 't1437010', 't314147'], 't8287': ['t32523'], 't186625':['t186626'], 't41665': ['t1489872'], 't117571':['t7742', 't7711', 't33213']}
root_nodes = ['t117571','t7742', 't7711', 't33213'] #bony vertebrates and before

#Full lineage
taxa_full_lineage = ['ENSP','t9606','t314146','t32525','t40674','t32524','t8287']
taxa_full_lineage.extend(node_projection['t32525'])
taxa_full_lineage.extend(node_projection['t8287'])

#Human only
taxa_human_only = ['ENSP','t9606']


## Load data 
with open(pfam_results_file, 'r') as results:
	pfam_results = json.load(results)
with open(meme_results_file, 'r') as results:
	meme_results = json.load(results)
with open(human_mapping_file, 'r') as human_proteins:
	og_human_mapping = json.load(human_proteins)
	
evo_events_data = pfam_results['repeat']
evo_events_data.update(meme_results['repeat'])



def make_human_summary_df(data,taxa): 

	human_summary_df_cols = ['identifier', 'human_protein_cnt', 'domain', 'clan', 'human_dup','root_dup','netto_dup','human_frac']
	human_summary_df_rows = []
	
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
			continue
		
		root_dup = 0
		human_dup = 0
		netto_dup = 0
		human_frac = 0
		
		for node,cnt in data['events']['dup'].items():
			node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))
				
			if node in root_nodes: root_dup += cnt
			elif node in taxa: human_dup += cnt
			else: netto_dup += cnt
			
		netto_dup += human_dup
		if netto_dup>0: human_frac=human_dup/netto_dup
		
		if len(human_list) > 0:
			human_summary_df_rows.append( [og_id_domain, len(human_list), pre_model, clan, human_dup, root_dup, netto_dup,human_frac])
	
	human_summary_df = pd.DataFrame(human_summary_df_rows, columns=human_summary_df_cols)
			
	return(human_summary_df)



## Analysis: full lineage

human_full_lineage_df = make_human_summary_df(evo_events_data,taxa_full_lineage)
human_full_lineage_df.to_csv(output_full_lineage,sep='\t')


## Analysis: human only
human_only_df = make_human_summary_df(evo_events_data,taxa_human_only)
human_only_df.to_csv(output_human_only,sep='\t')

