## IvB May 2018
# Mask Fasta sequences for PFAM hits
# Derived form parse_tblout_init.py
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

input_ext = '.tblout'
output_ext = '.fa' #hit replaced later by best pfam model of clan
outfile_path = "" #defined based on hmm input path, fasta output path and extension

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

if pre.pfam_hmm_path in hmmr_tbl:
	outfile_path = hmmr_tbl.replace(pre.pfam_hmm_path, pre.fasta_chopped_path) 
	outfile_path = outfile_path.replace(input_ext, output_ext) 
	og_id = hmmr_tbl[ len(pre.pfam_hmm_path):-(len(input_ext)) ]
else: 
	quit('Need valid .tblout from full pfam scan')

#read clans
pfam_clans = pre.get_pfam_clans(pre.pfam_clans_file)

#read fasta
fasta = pre.read_fasta(fasta_file)

#hmm data
repeats, pfam_hits = pre.parse_domtblout(hmmr_tbl)

best_hit_repeats = pre.get_best_hit(repeats)

#write entire sequences if no best hits available
if best_hit_repeats is None or len(best_hit_repeats) == 0: 
	with open(outfile_path, 'w') as outfile:	
		for protein_uri in fasta:
			outfile.write('>'+protein_uri+'\n'+fasta[protein_uri]+'\n')
	quit()
			
hit_name = '_'.join(best_hit_repeats.keys()[0].split('_')[1:-1]) 
print(og_id,hit_name)

#otherwise write masked output file
pre.write_chopped_file(outfile_path, best_hit_repeats, fasta)

