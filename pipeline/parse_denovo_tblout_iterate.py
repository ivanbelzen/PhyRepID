## IvB March 2018
# Parses HMMscan domain table output

# Input: HMMscan domain tblout, fasta file for OGs
# Needs: Pfam gathering cutoff, pfam clans, 
# Outputs: Fasta file with repeats for each OG-hit combination 

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
	quit('Need valid .tblout from full domain scan')

#read fasta
fasta = pre.read_fasta(fasta_file)

#hmm data
repeats, pfam_hits = pre.parse_domtblout(hmmr_tbl, 'iterative')

best_hit_repeats = pre.get_best_hit(repeats)
if best_hit_repeats is None or len(best_hit_repeats) == 0: quit()

hit_name = '_'.join(best_hit_repeats.keys()[0].split('_')[1:-1]) 
print(hit_name)

pre.write_fasta_file(outfile_path, best_hit_repeats, fasta, padding)
