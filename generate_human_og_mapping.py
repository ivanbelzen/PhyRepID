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

human_protein_file = pre.root+'og_human_mapping_reroot.json'

##Generate OG to human gene mapping
#make dict of lists of og_id -> [human protein ids] for pfam and meme. easier to search schaper and then merge with evo events dataframe
# also need to include orthologs list...

og_id_human_dict = {'pfam':{},'meme':{}}

pfam_results_file = pre.root+'pfam_evo_events_reroot.json'
meme_results_file = pre.root+'meme_evo_events_reroot.json'
pfam_summary_file = pre.root+'pfam_summary.json'
meme_summary_file = pre.root+'meme_summary.json'

with open(pfam_results_file, 'r') as results:
	pfam_results = json.load(results)

with open(meme_results_file, 'r') as results:
	meme_results = json.load(results)

with open(pfam_summary_file, 'r') as pfam_output:
	pfam_summary = json.load(pfam_output)
	
with open(meme_summary_file, 'r') as meme_output:
	meme_summary = json.load(meme_output)

pfam_repeat_data = pfam_results['repeat']
meme_repeat_data = meme_results['repeat']

print(len(meme_summary), len(pfam_summary))

list_of_all_proteins = []

for og_id_domain,data in pfam_summary.items(): 
	pre_model = '_'.join(og_id_domain.split('_')[2:])
	orthologs_list = pfam_summary[og_id_domain][pre_model]['orthologs_dict'].keys()
	human_proteins = [x for x in orthologs_list if x.startswith('ENSP0')]
	og_id_human_dict['pfam'][og_id_domain] = {}
	for hprot in human_proteins:
		og_id_human_dict['pfam'][og_id_domain][hprot]=orthologs_list
		list_of_all_proteins.append(hprot)
		list_of_all_proteins.extend(orthologs_list)
		
for og_id,data in meme_summary.items(): 
	orthologs_list = meme_summary[og_id]['motif']['orthologs_dict'].keys()
	human_proteins = [x for x in orthologs_list if x.startswith('ENSP0')]
	og_id_human_dict['meme'][og_id]={}
	for hprot in human_proteins:
		og_id_human_dict['meme'][og_id][hprot]=orthologs_list
		list_of_all_proteins.append(hprot)
		list_of_all_proteins.extend(orthologs_list)

print(len(np.unique(np.array(list_of_all_proteins))))
print(len(list_of_all_proteins))

'''		
with open(human_protein_file,'w') as output:
	output.write(json.dumps(og_id_human_dict))
'''	
