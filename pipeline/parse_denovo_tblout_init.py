## IvB March 2018
# Parses HMMscan domain table output
# Filter OGs based on minimum #repeats in #species
# Prepare set of repeats for iterative domain annotation
# Writes a fasta file for each clan with all of the hits 

# Input: HMMscan domain tblout, fasta file for OGs
# Outputs: Fasta file with repeats for each OG-hit combination 
# Logs: meme_hmm_results.json for statistics, dict of repeat proteins because of species threshold.

# ENV +-5 
import sys, os
import csv
import operator
import re
import json
import pipeline_methods as pre

hmmr_tbl_input_path = pre.denovo_tblout_path
fasta_output_path = pre.denovo_repeats_path
input_ext = '.tblout'
output_ext = '.fa' 
outfile_path = "" #defined based on hmm input path, fasta output path and extension

species_mapping_file = pre.species_mapping_file
species_mapping = {}
results_file = pre.meme_hmm_results_file
results_dict = {}
excluded_ogs_file = pre.meme_excluded_ogs_file
excluded_ogs={}

repeat_threshold = 3 # >=
species_threshold = 1 # >=
protein_threshold = 4 # >=

#strict slicing of repeats
spacing = 'ali'
padding = 0

hmmr_tbl = str(sys.argv[1]) #hammer output
repeats = {}
pfam_hits = {}
fasta_file = str(sys.argv[2]) #orthologs file
fasta = {}

if hmmr_tbl_input_path in hmmr_tbl:
	outfile_path = hmmr_tbl.replace(hmmr_tbl_input_path, fasta_output_path) 
	outfile_path = outfile_path.replace(input_ext, output_ext) 
	gene_name = hmmr_tbl[ len(hmmr_tbl_input_path):-(len(input_ext)) ]
else: 
	quit('Need valid .tblout from full pfam scan')

if pre.file_notempty(results_file):
	with open(results_file,'r') as log:
		results_dict = json.load(log)
else:
	results_dict = {}	
if gene_name not in results_dict: results_dict[gene_name]={}

#init log
#init excluded file
if pre.file_notempty(excluded_ogs_file):
	with open(excluded_ogs_file,'r') as log:
		excluded_ogs = json.load(log)
else: 
	excluded_ogs['settings']={'repeat_threshold':repeat_threshold, 'species_threshold':species_threshold }
	excluded_ogs['og_id'] = []		


#read clans
pfam_clans = {}

#read fasta
fasta = pre.read_fasta(fasta_file)

#hmm data
repeats, pfam_hits = pre.parse_domtblout(hmmr_tbl, 'iterative')

#filter repeats based on thresholds
#save excldued in seperate file
#what if it would be pfam_hits with hit -> {protein uri: number}
for hit,protein_uri in pfam_hits.items():
	clan = 'clan'	
	
	#logging for statistics/analytics
	if clan not in results_dict[gene_name]: results_dict[gene_name][clan] = {}
	results_dict[gene_name][clan][hit]={'protein_uri':protein_uri,'score':repeats[clan][hit]['score']}
	
	#Repeat must be present in human and mouse, and comply to set thresholds
	
	if len([x for x in pfam_hits[hit] if pfam_hits[hit][x] >= repeat_threshold and 'ENSP0' in x]) < 1 or \
	   len([x for x in pfam_hits[hit] if 'ENSMUSP' in x]) < 1 or \
	   len([x for x in pfam_hits[hit]]) < protein_threshold or \
	   (len([i_total for i_total in protein_uri.values() if i_total >= repeat_threshold]) < species_threshold):
		   
		del(repeats[clan][hit])
		continue

#save log
with open(results_file,'w') as log:
	log.write(json.dumps(results_dict))

best_hit_repeats = pre.get_best_hit(repeats)
if best_hit_repeats is None or len(best_hit_repeats) == 0: 
	excluded_ogs['og_id'].append(gene_name)
	
#store excluded 
with open(excluded_ogs_file,'w') as excl:
	excl.write(json.dumps(excluded_ogs))
if best_hit_repeats is None or len(best_hit_repeats) == 0: quit()	

pre.write_fasta_file(outfile_path, best_hit_repeats, fasta, padding)
