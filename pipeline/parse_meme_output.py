## IvB 
## August 2018
# Parse meme output and derive repeated motifs
# Log in meme_results_parse_log.json

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
log_file = pre.meme_results_parse_log_file

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
