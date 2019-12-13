## IvB March 2018
# Parses HMMscan domain table output of one hit only
# Use in iterative process to refine profile HMM model
# Writes a fasta file with repeats
# No filtering possible on model level, only domain > 0

# Input: HMMscan domain tblout of one hit only, fasta file for OGs
# Needs: Pfam clans, 
# Outputs: Fasta file with repeats for the input OG-hit combination 
# Logs: none?

import sys, os
import csv
import operator
import re
import json
import pipeline_methods as pre

#TODO INCLUDE new best hits function?

root = '/home/ianthe/protein-repeat-evolution/'
hmmr_tbl_input_path = root+'pfam/tblout/'
fasta_output_path = root+'pfam/repeats/'
input_ext = '.tblout'
output_ext = '.fa' 
outfile_path = "" #defined based on hmm input path, fasta output path and extension

root_old = '/home/ianthe/pipeline/'
species_mapping_file = root+'ensembl_stable_id_species.json'
species_mapping = {}
pfam_clans_file = root_old+'Pfam-A.clans.tsv' 
pfam_clans = {}

repeat_threshold = 3 # >=
species_threshold = 1 # >=

#strict slicing of repeats
spacing = 'ali'
padding = 0

hmmr_tbl = str(sys.argv[1]) #hammer output, of one hit only
repeats = {}
pfam_hits = {}
fasta_file = str(sys.argv[2]) #orthologs file
fasta = {}

if hmmr_tbl_input_path in hmmr_tbl:
	outfile_path = hmmr_tbl.replace(hmmr_tbl_input_path, fasta_output_path) 
	outfile_path = outfile_path.replace(input_ext, output_ext) 
else: 
	quit('Need valid .tblout from profile pfam scan of one hit only')	
		

#read clans
pfam_clans = pre.get_pfam_clans(pfam_clans_file)

#read fasta
fasta = pre.read_fasta(fasta_file)

#hmm data
repeats, pfam_hits = pre.parse_domtblout(hmmr_tbl, 'iterative')

#should only be one clan and one hit
best_hit_repeats = pre.get_best_hit(repeats)
pre.write_fasta_file(outfile_path, best_hit_repeats, fasta, padding)
