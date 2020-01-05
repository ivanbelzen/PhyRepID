# IvB Analyse schaper overlap with PhyRepID data
# 17-09-2018
# Last edit 01-01-2020
# Input PhyRepID: og_human_mapping.json 
# Input Schaper: eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick
# Output: schaper_comparison.json - summary used in generate_phyrepid_results.py 
# Output: phyrepid_schaper_comparison.tsv - detailed df with info per protein (schaper_detailed_comparison.json)

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

schaper_files = [pre.schaper_input_file_1, pre.schaper_input_file_2]
human_mapping_file = pre.human_mapping_file

schaper_data = {} # human ortholog clan pfam_hit relationship

schaper_summary_df_file = pre.schaper_comparison_file #used in phyrepid export
schaper_detailed_df_file = pre.root+"schaper_detailed_comparison.json"
schaper_detailed_export = pre.root+"phyrepid_schaper_comparison.tsv"

## Build Schaper data dictionary

#Only compare to species in our pipeline
species_selection_lst = []
with open(pre.species_mapping_file,'r') as sp_mapping:
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


## Make Schaper dataframe

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
		

# Export detailed, protein-wise comparision

schaper_df = pd.DataFrame(schaper_df_rows, columns=schaper_df_columns)
schaper_df.to_csv(schaper_detailed_export, separator="\t")

with open(schaper_detailed_df_file, 'w') as output:
	output.write(json.dumps(schaper_df.to_dict()))


# Export summary of comparision per OG
schaper_summary_df = pd.DataFrame(schaper_summary_rows, columns=schaper_summary_cols)
schaper_summary_df['schaper_positive'] = np.where(((schaper_summary_df['schaper_sli_sep'] > 0) | (schaper_summary_df['schaper_per_sep'] > 0) ), 'True', 'False')

with open(schaper_summary_df_file, 'w') as output:
	output.write(json.dumps(schaper_summary_df.to_dict()))
