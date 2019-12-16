## IvB March 2018
# Parses HMMscan domain table output
# Use gathing cutoff for sequences
# Select best Pfam model for each clan
# Based on sequence bitscore maximalisation
# Filter OGs based on minimum #repeats in #species
# Prepare set of repeats for iterative domain annotation
# Writes a fasta file for each clan with all of the hits 

# Input: HMMscan domain tblout, fasta file for OGs
# Needs: Pfam gathering cutoff, pfam clans, 
# Outputs: Fasta file with repeats for each OG-hit combination 
# Logs: HMM_results.json for statistics, dict of repeat proteins because of species threshold.

# ENV +-5 
import sys, os
import csv
import operator
import re
import json
import pipeline_methods as pre

root = '/home/ianthe/protein-repeat-evolution/'
hmmr_tbl_input_path = root+'pfam/hmm/'
fasta_output_path = root+'pfam/repeats/'
input_ext = '.tblout'
output_ext = '_hit.fa' #hit replaced later by best pfam model of clan
outfile_path = "" #defined based on hmm input path, fasta output path and extension

root_old = '/home/ianthe/pipeline/'
species_mapping_file = root+'ensembl_stable_id_species.json'
species_mapping = {}
pfam_gacutoff_file = root+'ga_cutoff.tsv'
pfam_gacutoff = {}
pfam_clans_file = root+'Pfam-A.clans.tsv' 
pfam_clans = {}
hmm_results_file = root+'hmm_results.json'
hmm_results_dict = {}
excluded_repeats_file = root+'excluded_repeats.json'
excluded_repeats={}

repeat_threshold = 3 # >=
species_threshold = 1 # >=

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


def file_notempty(filepath):
	if(os.path.exists(filepath) and os.path.getsize(filepath) > 0):
		return True
	else:
		return False
if file_notempty(hmm_results_file):
	with open(hmm_results_file,'r') as log:
		hmm_results_dict = json.load(log)
else:
	hmm_results_dict = {}	
if gene_name not in hmm_results_dict: hmm_results_dict[gene_name]={}

#init log
#init excluded repeats
if file_notempty(excluded_repeats_file):
	with open(excluded_repeats_file,'r') as log:
		excluded_repeats = json.load(log)
else: 
	excluded_repeats['settings']={'repeat_threshold':repeat_threshold, 'species_threshold':species_threshold }
 
		
		

#read clans
pfam_clans = pre.get_pfam_clans(pfam_clans_file)

#read fasta
fasta = pre.read_fasta(fasta_file)

#hmm data
repeats, pfam_hits = pre.parse_domtblout(hmmr_tbl)

#filter repeats based on thresholds
#save excldued in seperate file
#what if it would be pfam_hits with hit -> {protein uri: number}
for hit,protein_uri in pfam_hits.items():
	clan = pfam_clans[hit]	
	
	#logging for statistics/analytics
	if clan not in hmm_results_dict[gene_name]: hmm_results_dict[gene_name][clan] = {}
	hmm_results_dict[gene_name][clan][hit]={'protein_uri':protein_uri,'score':repeats[clan][hit]['score']}
	
	#Repeat must be present in human and mouse, and comply to set thresholds
	
	if len([x for x in pfam_hits[hit] if pfam_hits[hit][x] >= repeat_threshold and 'ENSP0' in x]) < 1 or \
	   len([x for x in pfam_hits[hit] if 'ENSMUSP' in x]) < 1 or \
	   (len([i_total for i_total in protein_uri.values() if i_total >= repeat_threshold]) < species_threshold):
		   
		if clan not in excluded_repeats: excluded_repeats[clan] = {}
		excluded_repeats[clan][hit] = repeats[clan][hit]
		del(repeats[clan][hit])
		continue

#save log
with open(hmm_results_file,'w') as log:
	log.write(json.dumps(hmm_results_dict))

#store excluded repeats
with open(excluded_repeats_file,'w') as excl:
	excl.write(json.dumps(excluded_repeats))


best_hit_repeats = pre.get_best_hit(repeats)
if best_hit_repeats is None or len(best_hit_repeats) == 0: quit()

hit_name = '_'.join(best_hit_repeats.keys()[0].split('_')[1:-1]) 
print(hit_name)

outfile = outfile_path.replace('hit', hit_name) #for each clan, best hit
pre.write_fasta_file(outfile, best_hit_repeats, fasta, padding)
