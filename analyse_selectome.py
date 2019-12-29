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

#analysis_path_selectome = pre.root+'analysis/selectome_07102018/'
analysis_path_selectome = pre.root+'analysis/selectome_01112018/'
human_protein_file = 'og_human_mapping_reroot.json'
#Generate OG to human gene mapping --> file generate_human_protein_mapping.py
selectome_df_file = 'selectome_df_reroot.json'

human_df_file = pre.root+'analysis/human_evo_events/full_lineage_trace/human_summary_df.json'
combined_summary_df_file = pre.root+'analysis/evo_events_06102018/pre_dataset_df.json'

human_prot_selected_file = 'human_prot_selected.tsv'
		
#make selectome df
'''
selectome = []
with open(human_prot_selected_file, 'r') as human_prot_selected:
	selectome = []
	for line in human_prot_selected:
		cols = line.split('\t')
		selectome.append(cols[0])

with open(human_protein_file, 'r') as human_proteins:
	og_human_mapping = json.load(human_proteins)

selectome_df_cols = ['identifier', 'protein_id', 'positive_selection']
selectome_df_rows = []

#og_human_mapping['meme'].update(og_human_mapping['pfam'])

for og_id,human_list in og_human_mapping['meme'].items(): 
	for human_protein,orthologs_list in human_list.items():	 
		if human_protein in selectome: selectome_df_rows.append([og_id, human_protein,1])
		else: selectome_df_rows.append([og_id, human_protein,0])
for og_id,human_list in og_human_mapping['pfam'].items(): 
	for human_protein,orthologs_list in human_list.items():	 
		if human_protein in selectome: selectome_df_rows.append([og_id, human_protein,1])
		else: selectome_df_rows.append([og_id, human_protein,0])
		

selectome_df = pd.DataFrame(selectome_df_rows, columns=selectome_df_cols)
		
with open(selectome_df_file, 'w') as output:
	output.write(json.dumps(selectome_df.to_dict()))
'''
## Analysis
	
#read selectome df
with open(selectome_df_file, 'r') as output:
	selectome = json.load(output)
	selectome_df = pd.DataFrame.from_dict(selectome)
#read datasets
with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)	
with open(human_df_file, 'r') as output:
	human = json.load(output)
	human_df = pd.DataFrame.from_dict(human)

sys.stdout = open(analysis_path_selectome+'output.txt', 'w')

selectome_positive_set = selectome_df.loc[selectome_df['positive_selection'] == 1]

#full dataset comparison
print("Comparison with full dataset")
print('#OGs in combined summary',len(combined_summary_df.index), 'selectome',selectome_df.identifier.nunique())

positive_set = combined_summary_df.loc[combined_summary_df['duplication_score']>0]
fast_set = combined_summary_df.loc[combined_summary_df['duplication_score']>combined_summary_df['duplication_score'].std()]
print('#OGs positive in Full', len(positive_set.index))
print('#OGs positive in Selectome', selectome_positive_set.identifier.nunique())
print('#OGs overlap:',positive_set.loc[positive_set['identifier'].isin(selectome_positive_set['identifier'])].identifier.nunique())
print('#OGs overlap fast set:',fast_set.loc[fast_set['identifier'].isin(selectome_positive_set['identifier'])].identifier.nunique())

print(fast_set.loc[fast_set['identifier'].isin(selectome_positive_set['identifier'])].gene_symbol.unique())

#human dataset comparison
print("Comparison with human dataset")
#human_df['human_duplication_score']= (human_df['netto_dup']-human_df['netto_dup'].mean()) / human_df['human_protein_cnt']	

#human_positive_set = human_df.loc[human_df['human_duplication_score']>0]
human_positive_set = human_df.loc[human_df['dev_from_mean']>0]
print('#OGs positive in Human lineage', len(human_positive_set.index))
print('#OGs overlap:',selectome_positive_set.loc[selectome_positive_set['identifier'].isin(human_positive_set['identifier'])].identifier.nunique())
human_positive_set = human_df.loc[human_df['netto_dup']>0]
print('#OGs overlap netto dup > 0:',selectome_positive_set.loc[selectome_positive_set['identifier'].isin(human_positive_set['identifier'])].identifier.nunique())

print("Redefine human positive set as only the overlapping part with the dynamic set")
human_positive_set = human_df.loc[(human_df['dev_from_mean']>0)&(human_df['identifier'].isin(positive_set['identifier']))]
human_fast_set = human_df.loc[(human_df['dev_from_mean']>0)&(human_df['identifier'].isin(fast_set['identifier']))]
print('#OGs positive in Human lineage', len(human_positive_set.index))
print('#OGs fast in Human lineage', len(human_fast_set.index))

print('#OGs overlap:',selectome_positive_set.loc[selectome_positive_set['identifier'].isin(human_positive_set['identifier'])].identifier.nunique())
print('#OGs overlap fast:',selectome_positive_set.loc[selectome_positive_set['identifier'].isin(human_fast_set['identifier'])].identifier.nunique())


print("Check controls")
controls = {'PRDM9': 'ENSGT00530000063157_ENSG00000164256_zf-H2C2_2',\
'KNL1': 'ENSGT00410000025918_ENSG00000137812',\
'VCAM1': 'ENSGT00830000128299_ENSG00000162692_Ig_3',\
 'CDC20': 'ENSGT00870000136444_ENSG00000117399_WD40',
 'AHNAK': 'ENSGT00530000063716_ENSG00000124942_DUF945',\
 'PRX': 'ENSGT00530000063716_ENSG00000105227_Cornifin',\
 'SPATA31A3':'ENSGT00530000063191_ENSG00000275969',\
 'SPATA31A6':'ENSGT00530000063191_ENSG00000185775' }

for gene,og_id in controls.items():
	in_dataset = selectome_positive_set.loc[selectome_positive_set['identifier']==og_id]
	print(gene)
	if not in_dataset.empty: print("yes")
	else: print("no")


quit()
human_fast_set = human_df.loc[human_df['human_frac']>0.2]
print(selectome_positive_set.loc[selectome_positive_set['identifier'].isin(human_fast_set['identifier'])].identifier.unique())
