## IvB August 2018
# Last edit: 30 december 2019
## Generate repeat stats from tblout for both Pfam and Meme
# Parse only finished OGs that have Treefix annotate output files: glob treefix path *.treefix.mpr.recon  
# output:  pfam_repeat_stats.json,  meme_repeat_stats.json
# output optional: pfam_repeat_stats_initial.json pfam_repeat_stats_initial_filtered.json
# output optional (false default): report OGs that have repeats in less than 4 proteins 
# possibly exclude later because of less trustworthy trees: exclude_ogs_lt4_pfam.lst exclude_ogs_lt4_meme.lst

import sys, os, glob, json
import numpy as np
import pipeline_methods as pre
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

species_mapping = pre.get_species_mapping(pre.species_mapping_file)

pfam_output_file="pfam_repeat_stats.json"
meme_output_file="meme_repeat_stats.json"

flag_generate_initial=True
pfam_initial_path = pre.pfam_hmm_path
pfam_initial_output_file="pfam_repeat_stats_initial.json"
pfam_initial_filtered_output_file="pfam_repeat_stats_initial_filtered.json"

#Report OGs that have repeats in less than 4 proteins, possibly exclude if desired.
flag_report_lt4=False
pfam_exclude_ogs_file="exclude_ogs_lt4_pfam.lst"
meme_exclude_ogs_file="exclude_ogs_lt4_meme.lst"

## Generate PFAM summary
pfam_summary = {}
pfam_file_list = glob.glob(pre.pfam_treefix_path + '*.treefix.mpr.recon' ) 

excluded_ogs = []
for pfam_file in pfam_file_list:	
	og_hit_id = pfam_file[len(pre.pfam_treefix_path):-len('.treefix.mpr.recon')]
	pfam_tblout = pre.pfam_tblout_path+og_hit_id+'.tblout'
	if not pre.file_notempty(pfam_tblout):continue
	repeats = pre.parse_domtblout_stats(pfam_tblout)
	
	domain = '_'.join(og_hit_id.split('_')[2:])
	
	if len(repeats[domain]['orthologs_dict'].items()) < 4:
		excluded_ogs.append(og_hit_id)
		continue
			
	pfam_summary[og_hit_id] = repeats	

print("Completed Pfam output: ",len(pfam_file_list), " entries in stats: ",len(pfam_summary))

with open(pre.root+pfam_output_file, 'w') as output:
	output.write(json.dumps(pfam_summary))

if flag_report_lt4:
	print("Pfam OGs that have less than 4 proteins with repeats:")
	with open(pre.root+pfam_exclude_ogs_file,'w') as output:
		for x in excluded_ogs:
			output.write(x+'\n')

##

## Generate PFAM summary before repeat detection optimization
if flag_generate_initial: 
	pfam_file_list = glob.glob(pfam_initial_path + '*.tblout' ) 
		
	pfam_summary_initial = {}
	for pfam_tblout in pfam_file_list:	
		og_hit_id = pfam_tblout[len(pfam_initial_path):-len('.tblout')]
		repeats = pre.parse_domtblout_stats(pfam_tblout, True)
		pfam_summary_initial[og_hit_id] = repeats
		
	with open(pre.root+pfam_initial_output_file, 'w') as output:
		output.write(json.dumps(pfam_summary_initial))

	#Filter: match up OG ids from initial to the OG-domain ids from later to compare improvement in repeat detection
	filtered_pfam_summary_initial = {}
	for og_id_domain in pfam_summary.keys():
		og_id = '_'.join(og_id_domain.split('_')[0:2])
		domain = og_id_domain[len(og_id)+1:]
		if og_id_domain not in filtered_pfam_summary_initial: filtered_pfam_summary_initial[og_id_domain] = {}
		filtered_pfam_summary_initial[og_id_domain][domain] = pfam_summary_initial[og_id][domain]

	with open(pre.root+pfam_initial_filtered_output_file,'w') as output:
		output.write(json.dumps(filtered_pfam_summary_initial))
		
	print("Pfam filtered and initial: ",len(filtered_pfam_summary_initial),len(pfam_summary_initial))

## Generate MEME summary
meme_summary = {}
meme_file_list = glob.glob(pre.denovo_treefix_path + '*.treefix.mpr.recon' ) 

excluded_ogs = []		
for meme_repeats in meme_file_list:	
	og_id = meme_repeats[len(pre.denovo_treefix_path):-len('.treefix.mpr.recon')]
	meme_tblout = pre.denovo_tblout_path+og_id+'.tblout'
	
	repeats = pre.parse_domtblout_stats(meme_tblout)
	meme_summary[og_id] = repeats
	meme_summary[og_id]['consensus'] = ''
	
	if len(repeats[domain]['orthologs_dict'].items()) < 4:
		excluded_ogs.append(og_id)
		continue
	
with open(pre.root+meme_output_file, 'w') as output:
	output.write(json.dumps(meme_summary))

print("Completed MEME output: ",len(meme_file_list), " entries in stats: ",len(meme_summary))


if flag_report_lt4:
	print("MEME OGs that have less than 4 proteins with repeats:")
	with open(pre.root+meme_exclude_ogs_file,'w') as output:
		for x in excluded_ogs:
			output.write(x+'\n')
			
##

	
