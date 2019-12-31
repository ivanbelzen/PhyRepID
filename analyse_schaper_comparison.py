## NOTE PARTIAL

# IvB Analyse schaper overlap with PhyRepID data
# 17-09-2018
# Last edit 31-12-2019
# Input PhyRepID: pfam_evo_events.json meme_evo_events.json 
# Input Schaper: eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick
# Output: schaper_comparison.json

# PhyRepID evo events contains dup/loss per node/ortholog in an OG and in total 

# Schaper summary contains info per pair of orthologs, protein IDs
# Hence to compare it is necessary to collapse info on OG level: so for every OG, consider the human proteins in it
# Uses OG to human gene mapping --> file generate_human_protein_mapping.py
# Look at what Schaper has to say about the cherries with these human proteins and vertebrates we also used
# Ignore none and not assigned, only look at  slightly/perfectly separated and conserved 

import sys,json
import pandas as pd
import pipeline_methods as pre
import numpy as np
import gzip 

schaper_files = [pre.root+'resources/eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick.gz', pre.root+'resources/eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick.gz']
human_mapping_file = pre.root+'og_human_mapping.json'

schaper_data = {} # human ortholog clan pfam_hit relationship

schaper_summary_df_file = 'schaper_comparison.json' #used in phyrepid export
schaper_df_file = 'schaper_detailed_df_reroot.json' #maybe interesting as export table s well

## Build Schaper data dictionary

#Only compare to species in our pipeline
species_selection_lst = []
with open(pre.species_file,'r') as sp_mapping:
	species_dict = json.load(sp_mapping)
	for item in species_dict:
		species_selection_lst.append(item['ensembl_stable_id'])					
 
#Load Pfam data
clans = pre.get_pfam_clans(pre.pfam_clans_file) #pfam hit to clan mapping
pfam_accession = {} #maps pfam accession ID to hit name
with open(pre.pfam_clans_file, 'r') as pfam_file: 
	for rows in pfam_file:
		cols = rows.strip().split('\t')
		pfam_accession[cols[0]] = cols[3]	
			
for filename in schaper_files:
	with gzip.open(filename, 'r') as f:
		for line in f: print(line)
			if line[0] != '>': continue
			param = {}
				
			for col in line.strip().split(' '):
				key,value = col.split(':')
				param[key] = value
				
			#ignore all inconclusive entries		
			if 'Pairwise_phylogeny_type' not in param: continue
			
			phylogeny_type = param['Pairwise_phylogeny_type']
			human_protein_id = param['Ensembl_Protein_ID_Human_Ortholog']
			ortholog = param['Ensembl_Protein_ID_Second_Ortholog']
			identifier = ortholog.split('0') #assumption correct that always starts with 0??
			identifier = identifier[0][:-1] #only need first part (species) and remove trailing P as indication for protein

			#Only compare species also in our dataset					
			if identifier not in species_selection_lst: continue
				
			if human_protein_id not in schaper_data:
				schaper_data[human_protein_id] = {}
			
			if ortholog not in schaper_data[human_protein_id]:
				schaper_data[human_protein_id][ortholog] = {}
			
			if "PFAM" in filename: 
				pfam_id = param['TR_detection_type']
				pfam_name = pfam_accession[pfam_id] if pfam_id in pfam_accession else 'denovo'
				pfam_clan = clans[pfam_name] if pfam_name in clans else 'denovo'
				
				if pfam_clan not in schaper_data[human_protein_id][ortholog]:
					schaper_data[human_protein_id][ortholog][pfam_clan] = {}
				schaper_data[human_protein_id][ortholog][pfam_clan][pfam_name]=phylogeny_type
			else:
				schaper_data[human_protein_id][ortholog]["denovo"]["denovo"]=phylogeny_type

###


# Make Schaper dataframe

with open(human_mapping_file, 'r') as human_proteins:
	mapping = json.load(human_proteins)

og_human_mapping = mapping["pfam"].update(mapping["meme"])  #make one big dictionary to prevent looping twice

schaper_df_columns = ['identifier', 'human_protein', 'ortholog_id', 'schaper_model', 'schaper_conclusion']
schaper_df_rows = []

schaper_summary_cols = ['identifier','schaper_model', 'schaper_sli_sep', 'schaper_per_sep', 'schaper_model_cnt'] 
schaper_summary_rows = [] #per OG: model, clan, netto dup, loss, perfectly separated, slightly separated

for og_id_domain,human_list in og_human_mapping.items(): 
	models_list = [];schaper_sli_sep = 0;schaper_per_sep = 0;common_model = '';
	
	for human_protein,orthologs_list in human_list.items():	 	 
		if human_protein in schaper_data:
			schaper_entry = { orth:vals for (orth,vals) in schaper_data[human_protein].items() if orth in orthologs_list }
			cherries_cnt = len(schaper_entry)
			for orth,vals in schaper_entry.items():
				for clan,hits in vals.items():
					for model,schaper_conclusion in hits.items():
						if model not in models_list: models_list.append(model)
						if 'perfectly_separated' in schaper_conclusion: schaper_per_sep+=1
						elif 'strongly_separated' in schaper_conclusion: schaper_sli_sep+=1
						schaper_df_rows.append([og_id_domain, human_protein, orth, model, schaper_conclusion]) 
					
			if len(models_list) > 0:
				common_model = max(set(models_list), key=models_list.count)		
			schaper_summary_rows.append([og_id_domain,common_model,schaper_sli_sep,schaper_per_sep,cherries_cnt] )
		
quit()

#write detailed df
schaper_df = pd.DataFrame(schaper_df_rows, columns=schaper_df_columns)
#write schaper df
with open(schaper_df_file, 'w') as output:
	output.write(json.dumps(schaper_df.to_dict()))


#summary df
schaper_summary_df = pd.DataFrame(schaper_summary_rows, columns=schaper_summary_cols)
schaper_summary_df['schaper_positive'] = np.where(((schaper_summary_df['schaper_sli_sep'] > 0) | (schaper_summary_df['schaper_per_sep'] > 0) ), 'True', 'False')

#write schaper df
with open(schaper_summary_df_file, 'w') as output:
	output.write(json.dumps(schaper_summary_df.to_dict()))


####
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
