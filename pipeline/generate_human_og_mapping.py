# IvB 
# 17-09-2018
##Generate OG to human gene mapping
# Mapping for Schaper comparison and Human-lineage overview 
# make dict of lists of og_id -> [human protein ids] for Pfam and MEME.
# Last edit 30-12-2019

# To compare schaper and PRE collapse info on OG level, so consider the human proteins in every OG

import sys,json
import pandas as pd
import pipeline_methods as pre
import numpy as np

human_mapping_file = pre.human_mapping_file

og_id_human_dict = {'pfam':{},'meme':{}}

pfam_repeat_stats_file = pre.pfam_repeat_stats_file
meme_repeat_stats_file = pre.meme_repeat_stats_file

with open(pfam_repeat_stats_file, 'r') as pfam_output:
	pfam_repeat_stats = json.load(pfam_output)
	
with open(meme_repeat_stats_file, 'r') as meme_output:
	meme_repeat_stats = json.load(meme_output)

for og_id_domain,data in pfam_repeat_stats.items(): 
	pre_model = '_'.join(og_id_domain.split('_')[2:])
	orthologs_list = pfam_repeat_stats[og_id_domain][pre_model]['orthologs_dict'].keys()
	human_proteins = [x for x in orthologs_list if x.startswith('ENSP0')]
	og_id_human_dict['pfam'][og_id_domain] = {}
	for hprot in human_proteins:
		og_id_human_dict['pfam'][og_id_domain][hprot]=orthologs_list
		
		
for og_id,data in meme_repeat_stats.items(): 
	orthologs_list = meme_repeat_stats[og_id]['motif']['orthologs_dict'].keys()
	human_proteins = [x for x in orthologs_list if x.startswith('ENSP0')]
	og_id_human_dict['meme'][og_id]={}
	for hprot in human_proteins:
		og_id_human_dict['meme'][og_id][hprot]=orthologs_list

with open(human_mapping_file,'w') as output:
	output.write(json.dumps(og_id_human_dict))
