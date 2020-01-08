## IvB March 2018
# Parses HMMscan domain table output
# Use to finalize the iterative domain detection
# Writes a fasta file with repeats, using env and +- 5 padding

# Input: HMMscan domain tblout, fasta file for OGs
# Outputs: Fasta file with repeats for each OG-hit combination 
# Logs: meme_results_final.json for statistics

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
results_file = pre.meme_hmm_results_final_file
results_dict = {}

#strict slicing of repeats
spacing = 'env'
padding = 5

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

#read clans
pfam_clans = {}

#read fasta
fasta = pre.read_fasta(fasta_file)

#hmm data
repeats, pfam_hits = pre.parse_domtblout(hmmr_tbl, 'iterative', spacing)

#should only be one clan and one hit

best_hit_repeats = pre.get_best_hit(repeats)
pre.write_fasta_file(outfile_path, best_hit_repeats, fasta, padding)

#print(repeats)
#logging for statistics
for hit,protein_uri in pfam_hits.items():
	clan = 'clan'	
	
	#logging for statistics/analytics
	if clan not in results_dict[gene_name]: results_dict[gene_name][clan] = {}
	results_dict[gene_name][clan][hit]={'protein_uri':protein_uri}
	
#save log
with open(results_file,'w') as log:
	log.write(json.dumps(results_dict))

