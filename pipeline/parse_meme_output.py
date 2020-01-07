## IvB 
## August 2018
# Parse meme output and derive repeated motifs
# Log in meme_log_chop.json

## Input
# MEME html file in pre.meme_output_path
## Output
# Fasta file of hits in pre.denovo_meme_repeats_path
	
import sys, os
import csv
import operator
import re
import json
import pipeline_methods as pre

meme_html = str(sys.argv[1]) #meme output
meme_hits = {}
input_ext = '/meme.html'
output_ext = '.fa'
padding = 5
#log_file = pre.root+'meme_log_nox.json'
log_file = pre.root+'meme_log_chop.json'

if pre.meme_output_path in meme_html:
	outfile_path = meme_html.replace(pre.meme_output_path, pre.denovo_meme_repeats_path) 
	outfile_path = outfile_path.replace(input_ext, output_ext) 
	
	gene_name = meme_html[ len(pre.meme_output_path):-(len(input_ext)) ]
	print(gene_name)
else: 
	quit('Need valid MEME output file')

with open(meme_html, 'r') as meme: 
	meme_html = meme.read()
	json_text = re.compile('var data = ({.*?});', re.DOTALL)
	matches = json_text.search(meme_html)
	meme_json = json.loads(matches.group(1))

if pre.file_notempty(log_file):
	with open(log_file,'r') as log:
		meme_results = json.load(log)
else:
	meme_results = {}	
if gene_name not in meme_results: meme_results[gene_name]={'consensus':meme_json['motifs'][0]['id'], 'nsites':meme_json['motifs'][0]['nsites']}

with open(log_file,'w') as log:
	log.write(json.dumps(meme_results))

#minimum motifs to make a tree? 10?
if meme_json['motifs'][0]['nsites'] < 10 : quit()


#TODO: 
# Check if motifs contain 'X' sequences, count number of X'es in string
# Cutoff, max 10-20% of sequence. 
# ['motifs'][0]['len'] and ['motifs'][0]['id']
# consensus = ['motifs'][0]['id'].count('X')

#should not included X'es in chopped
if meme_json['motifs'][0]['id'].count('X') > (int(meme_json['motifs'][0]['len'])*0.1): quit()


motif = meme_json['motifs'][0]['sites']
#only one motif

with open(outfile_path, 'w') as outfile:
	for hits in motif:
		outfile.write('>seq'+str(hits['seq'])+'_pos'+str(hits['pos'])+'\n'+ \
		hits['lflank'][padding:]+hits['match']+hits['rflank'][0:padding]+'\n') #last 5 of first and first 5 of last

quit()

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


