## Ivb August 2018
## Make list of OGs on which MEME should be run

import sys, json, glob
import operator
import numpy as np

import pipeline_methods as pre

meme_og_set = []

hmm_results = {} #raw pfam domain detection
hmm_final_results = {} #output of iterative domain detection
og_hit_list = []
ogs_set = [] #All OGs with human protein coding gene and ... 
pfam_ogs_set = [] #OGs which have tree as output / hence PFAM hit with >3 repeats in human and ...
disjunct_set = [] #OGs which do not have PFAM tree as output
filtered_fasta_files = [] # OG >= sequences 

analysis_img_path = pre.root+'analysis/'

repeat_threshold = 3 # >=
seq_threshold = 6 # >=



with open(pre.hmm_results_file,'r') as log:
	hmm_results = json.load(log)

def get_og_id_from_og_hit(og_hit_list):
	og_list = []
	for og_hit in og_hit_list:
		og_id = '_'.join(og_hit.split('_')[0:2] )
		if og_id not in og_list:
			og_list.append(og_id)
	return(og_list)	
	
fasta_files = glob.glob(pre.fasta_path+"*.fa")
ogs_set = [ og_id[len(pre.fasta_path):-len(".fa")] for og_id in fasta_files ]
print('OG total: ',len(ogs_set))
print('OG ids missing from HMM results.json',len( [og_id for og_id in ogs_set if og_id not in hmm_results ]))
#should be 0

pfam_trees_files = glob.glob(pre.pfam_trees_path+"*.treefile")
pfam_trees_set = [  og_id[len(pre.pfam_trees_path):-len(".treefile")] for og_id in  pfam_trees_files]
pfam_ogs_set = get_og_id_from_og_hit(pfam_trees_set)

disjunct_set = [og_id for og_id in ogs_set if og_id not in pfam_ogs_set ]	
print('OG that do not have a treefile as output: ',len(disjunct_set))
print('OG that do have a treefile as output: ',len(pfam_ogs_set))

#Filter fasta files on OGs >= sequences
for f in fasta_files:
	with open(f,'r') as fasta_file:
		fasta = fasta_file.read() 
		seq_count = fasta.count('>ENS')
		if seq_count >= seq_threshold:
			filtered_fasta_files.append(f)

ogs_set = [ og_id[len(pre.fasta_path):-len(".fa")] for og_id in filtered_fasta_files ]	
print('OG total & seq count >= ',seq_threshold, ': ', len(ogs_set))

disjunct_set = [og_id for og_id in ogs_set if og_id not in pfam_ogs_set ]	
print('OG that do not have a treefile as output & seq count >= ',seq_threshold, ': ', len(disjunct_set))

if 'ENSGT00410000025918_ENSG00000137812' in disjunct_set: print('KNL1')

motif_files = glob.glob(pre.denovo_trees_path+"*.treefile")
motif_done_set = [ og_id[len(pre.denovo_trees_path):-len(".treefile")] for og_id in motif_files ]

'''
#For Snakefork meme run
with open(pre.root+'meme_dataset.json', 'w') as log:
	log.write(json.dumps(disjunct_set))
'''

todo_set = [og_id for og_id in disjunct_set if og_id not in motif_done_set ]	
print('OGs that were not yet done in previous run', len(todo_set))

#For Snakefork meme run
with open(pre.root+'meme_todo.json', 'w') as log:
	log.write(json.dumps(todo_set))

