# IvB analyse schaper overlap with PRE denovo data
# Edit 17-09-2018
# Input: pfam_evo_events.json meme_evo_events.json schaper_summary_denovo.json  schaper_summary.json (for pfam)

# Schaper summary contains info per pair of orthologs, protein IDs
# PRE evo events summary contains dup/loss per node/ortholog in an OG and in total 

# To compare schaper and PRE collapse info on OG level
# So for every OG, consider the human proteins in it
# Look at what Schaper has to say about the cherries with these human proteins and vertebrates
# Ignore none and not assigned, only look at  slightly/perfectly separated and conserved 

#Both PFAM and MEME on one big pile of data
# TODO search for MEME ogs also in pfam data and vice versa

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np
import seaborn as sns

#analysis_path_schaper = pre.root+'analysis/schaper_25102018/'
analysis_path_schaper = pre.root+'analysis/schaper_01112018/'
human_protein_file = 'og_human_mapping_reroot.json'
#Generate OG to human gene mapping --> file generate_human_protein_mapping.py

schaper_file = 'schaper_summary.json'
schaper_file_denovo = 'schaper_summary_denovo.json'

human_df_file = pre.root+'analysis/human_evo_events/full_lineage_trace/human_summary_df.json'
combined_summary_df_file = pre.root+'analysis/evo_events_06102018/pre_dataset_df.json'

schaper_summary_df_file = 'schaper_summary_df_reroot.json'
schaper_df_file = 'schaper_detailed_df_reroot.json'

#make schaper dataframe
'''

with open(schaper_file, 'r') as schaper:
	schaper_dict = json.load(schaper)
	
with open(schaper_file_denovo, 'r') as schaper:
	schaper_denovo_dict = json.load(schaper)
	
with open(human_protein_file, 'r') as human_proteins:
	og_human_mapping = json.load(human_proteins)

#clans_dict = pre.get_pfam_clans('Pfam-A.clans.tsv')

schaper_df_columns = ['identifier', 'human_protein', 'ortholog_id', 'schaper_model', 'schaper_conclusion']
schaper_df_rows = []

schaper_summary_cols = ['identifier','schaper_model', 'schaper_sli_sep', 'schaper_per_sep', 'schaper_model_cnt'] 
schaper_summary_rows = [] #per OG: model, clan, netto dup, loss, perfectly separated, slightly separated

meme_in_pfam = []
for og_id,human_list in og_human_mapping['meme'].items(): 
	models_list = [];schaper_sli_sep = 0;schaper_per_sep = 0;common_model = '';
	to_add_or_not_to_add = False
	
	for human_protein,orthologs_list in human_list.items():	 
		if human_protein in schaper_denovo_dict:
			to_add_or_not_to_add = True
			schaper_data = { orth:vals for (orth,vals) in schaper_denovo_dict[human_protein].items() if orth in orthologs_list }
			cherries_cnt = len(schaper_data)
			for orth,schaper_conclusion in schaper_data.items():
				if 'perfectly_separated' in schaper_conclusion: schaper_per_sep+=1
				elif 'strongly_separated' in schaper_conclusion: schaper_sli_sep+=1
				schaper_df_rows.append([og_id, human_protein, orth, schaper_conclusion]) 
		
		elif human_protein in schaper_dict: #exists appearantly
			to_add_or_not_to_add = True
			schaper_data = { orth:vals for (orth,vals) in schaper_dict[human_protein].items() if orth in orthologs_list }
			cherries_cnt = len(schaper_data)
			for orth,vals in schaper_data.items():
				for clan,hits in vals.items():
					for model,schaper_conclusion in hits.items():
						if model not in models_list: models_list.append(model)
				
						if 'perfectly_separated' in schaper_conclusion: schaper_per_sep+=1
						elif 'strongly_separated' in schaper_conclusion: schaper_sli_sep+=1	
			if len(models_list) > 0:
				common_model = max(set(models_list), key=models_list.count)	
			
	if to_add_or_not_to_add:		
		schaper_summary_rows.append([og_id,common_model,schaper_sli_sep,schaper_per_sep,cherries_cnt] )


pfam_in_denovo = []		
for og_id_domain,human_list in og_human_mapping['pfam'].items(): 
	models_list = [];schaper_sli_sep = 0;schaper_per_sep = 0;common_model = '';
	to_add_or_not_to_add = False
	
	for human_protein,orthologs_list in human_list.items():	 	 
		if human_protein in schaper_dict:
			to_add_or_not_to_add = True
			schaper_data = { orth:vals for (orth,vals) in schaper_dict[human_protein].items() if orth in orthologs_list }
			cherries_cnt = len(schaper_data)
			for orth,vals in schaper_data.items():
				for clan,hits in vals.items():
					for model,schaper_conclusion in hits.items():
						if model not in models_list: models_list.append(model)
				
						if 'perfectly_separated' in schaper_conclusion: schaper_per_sep+=1
						elif 'strongly_separated' in schaper_conclusion: schaper_sli_sep+=1
						schaper_df_rows.append([og_id_domain, human_protein, orth, model, schaper_conclusion]) 
					
			if len(models_list) > 0:
				common_model = max(set(models_list), key=models_list.count)		
		
		elif human_protein in schaper_denovo_dict:
			to_add_or_not_to_add = True
			schaper_data = { orth:vals for (orth,vals) in schaper_denovo_dict[human_protein].items() if orth in orthologs_list }
			cherries_cnt = len(schaper_data)
			for orth,schaper_conclusion in schaper_data.items():
				if 'perfectly_separated' in schaper_conclusion: schaper_per_sep+=1
				elif 'strongly_separated' in schaper_conclusion: schaper_sli_sep+=1
				#schaper_df_rows.append([og_id, human_protein, orth, schaper_conclusion]) 
			
	if to_add_or_not_to_add:		
		schaper_summary_rows.append([og_id_domain,common_model,schaper_sli_sep,schaper_per_sep,cherries_cnt] )			

#write detailed df
schaper_df = pd.DataFrame(schaper_df_rows, columns=schaper_df_columns)
#write schaper df
with open(schaper_df_file, 'w') as output:
	output.write(json.dumps(schaper_df.to_dict()))


#summary df
schaper_summary_df = pd.DataFrame(schaper_summary_rows, columns=schaper_summary_cols)
schaper_summary_df['schaper_pos'] = np.where(((schaper_summary_df['schaper_sli_sep'] > 0) |\
 (schaper_summary_df['schaper_per_sep'] > 0) ), 'True', 'False')

#write schaper df
with open(schaper_summary_df_file, 'w') as output:
	output.write(json.dumps(schaper_summary_df.to_dict()))

'''

## Analysis
	
#read schaper summary df
with open(schaper_summary_df_file, 'r') as output:
	schaper_summ = json.load(output)
	schaper_summary_df = pd.DataFrame.from_dict(schaper_summ)	
#read schaper df
with open(schaper_df_file, 'r') as output:
	schaper = json.load(output)
	schaper_detailed_df = pd.DataFrame.from_dict(schaper)	
#read datasets
with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)	
with open(human_df_file, 'r') as output:
	human = json.load(output)
	human_df = pd.DataFrame.from_dict(human)
with open(pre.root+'analysis/dataset_04102018/unit_cnt_difference_df.json','r') as output:
	unit_cnt_diff = json.load(output)
	unit_cnt_diff_df = pd.DataFrame.from_dict(unit_cnt_diff)
	unit_cnt_diff_df['identifier'] = unit_cnt_diff_df['og_id_domain']
with open(pre.root+'analysis/dataset_04102018/'+'unit_cnt_df.json', 'r') as output:
	unit_cnt = json.load(output)
	unit_cnt_df = pd.DataFrame.from_dict(unit_cnt)	
	
sys.stdout = open(analysis_path_schaper+'output.txt', 'w')

schaper_summary_df = schaper_summary_df.merge(unit_cnt_df[['identifier','cv']], on='identifier', how='left')
schaper_summary_df = schaper_summary_df.merge(unit_cnt_diff_df[['identifier','diff']], on='identifier', how='left')
	
#full dataset comparison
print("Comparison with full dataset")
print('#OGs in combined summary',len(combined_summary_df.index), 'human',len(human_df.index), 'schaper entries',len(schaper_summary_df.index))

positive_set = combined_summary_df.loc[combined_summary_df['duplication_score']>0]
fast_set = combined_summary_df.loc[combined_summary_df['duplication_score']>combined_summary_df['duplication_score'].std()]
schaper_positive_set = schaper_summary_df.loc[schaper_summary_df['schaper_pos'] == 'True']
print('#OGs positive in Full', len(positive_set.index))
print('#OGs positive in Schaper', len(schaper_positive_set.index))
print('#OGs overlap:',len(schaper_positive_set.loc[schaper_positive_set['identifier'].isin(positive_set['identifier'])]))

non_overlap_set = schaper_positive_set.loc[~schaper_positive_set['identifier'].isin(positive_set['identifier'])]
non_overlap_set = non_overlap_set.merge(combined_summary_df, on='identifier', how='left')
print("nonoverlap:",len(non_overlap_set.index))
non_overlap_set.to_csv(analysis_path_schaper+'non_overlap_full_schaper.csv')
#[['identifier','duplication_score','schaper_model','schaper_model_cnt','clan']].to_string())

#human dataset comparison
print("Comparison with human dataset")
#human_df['human_duplication_score']= (human_df['netto_dup']-human_df['netto_dup'].mean()) / human_df['human_protein_cnt']	

#human_positive_set = human_df.loc[human_df['human_duplication_score']>0]
human_positive_set = human_df.loc[human_df['dev_from_mean']>0]
schaper_positive_set = schaper_summary_df.loc[schaper_summary_df['schaper_pos'] == 'True']
human_fast_set = human_df.loc[human_df['dev_from_mean']>human_df['dev_from_mean'].std()]

print('#OGs positive in Human lineage', len(human_positive_set.index))
print('#OGs fast in Human lineage', len(human_fast_set.index))
print('#OGs overlap:',len(schaper_positive_set.loc[schaper_positive_set['identifier'].isin(human_positive_set['identifier'])]))

print("Redefine human positive set as only the overlapping part with the dynamic set")
human_positive_set = human_df.loc[(human_df['dev_from_mean']>0)&(human_df['identifier'].isin(positive_set['identifier']))]
human_fast_set = human_df.loc[(human_df['dev_from_mean']>0)&(human_df['identifier'].isin(fast_set['identifier']))]
print('#OGs positive in Human lineage', len(human_positive_set.index))
print('#OGs overlap:',len(schaper_positive_set.loc[schaper_positive_set['identifier'].isin(human_positive_set['identifier'])]))

print('#OGs fast in Human lineage', len(human_fast_set.index))
print('#OGs overlap:',len(schaper_positive_set.loc[schaper_positive_set['identifier'].isin(human_fast_set['identifier'])]))


'''
non_overlap_set = schaper_positive_set.loc[~schaper_positive_set['identifier'].isin(human_positive_set['identifier'])]
non_overlap_set = non_overlap_set.merge(human_df, on='identifier', how='left')
print("nonoverlap:",len(non_overlap_set.index))
#print(non_overlap_set[['identifier','human_duplication_score','schaper_model','schaper_model_cnt','domain','clan']].to_string())
non_overlap_set.to_csv(analysis_path_schaper+'non_overlap_human_schaper.csv')


print("Overlap human and full")
fast_set = combined_summary_df.loc[combined_summary_df['duplication_score']>combined_summary_df['duplication_score'].std()]

print('#OGs overlap positive:',len(positive_set.loc[positive_set['identifier'].isin(human_positive_set['identifier'])]))
print('#OGs overlap fast:',len(human_positive_set.loc[human_positive_set['identifier'].isin(fast_set['identifier'])]))
'''

##Schaper detailed df
print("Schaper detailed df")
print('#human proteins',schaper_detailed_df.human_protein.nunique())
print('#human proteins positive',schaper_detailed_df.loc[(schaper_detailed_df['schaper_conclusion']=='perfectly_separated') | (schaper_detailed_df['schaper_conclusion']=='strongly_separated')].human_protein.nunique())
