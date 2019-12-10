## IvB April 2018
## Retrieve gene trees from API for selected orthologs
## Query returns subtree of only indicated taxons from a larger gene tree
## The gene trees can be seen as orthologuous groups, containing both orthologs and paralogs
## Gene tree json is saved, and to be parsed by parse_genetree.py

## Gene tree json is saved as NHX tree, and to a mapping file
## Gene tree identifier: [ (gene_id, protein_id)  of all from subtree }
## log file: human gene_uri -> gene tree identifier
## Also writes fasta file for orthoguous groups

import requests, sys, json, os.path
from SPARQLWrapper import SPARQLWrapper, JSON
import pipeline_methods as pre
		
root = '/home/ianthe/protein-repeat-evolution/'
genetree_path = root+'genetrees_nhx/'
ensembl_path = root+'ensembl_api/'     
filtered_ogs_file = root+'orthologs_filtered.json'
log_file = 'log_retrieve_genetrees.txt'

#species = ['9606','10090','9258','13616','9305','8839','9031','13735','8364','7918','7955','28377','31033','69293']

species_list = pre.get_taxon_list('', [])
species_str = '&prune_taxon='.join([str(x) for x in species_list])

gene_id_uri = 'http://rdf.ebi.ac.uk/resource/ensembl/' 
ensembl_server = "http://rest.ensembl.org"
sparql=SPARQLWrapper('https://www.ebi.ac.uk/rdf/services/sparql')
sparql.setReturnFormat(JSON)

genes_to_genetrees_file = root+'genes_to_genetrees.json'
genes_to_genetrees = {}
redo_path = root+'redo_orth_parse.json'
redo = {}
genetrees_mapping_file = root+'genetrees_ogs.json'
genetrees_mapping = {}

if(pre.file_notempty(filtered_ogs_file)):
	with open(filtered_ogs_file,'r') as filtered_ogs:
		mapping = json.load(filtered_ogs)
else: 	
	print("Requires input file ortholog_filtered.json")
	quit()	

#reset log
with open(log_file,'w') as log:
	log.write('')
	

#Retrieve gene trees for OGs
# extract subtree containing all orthologs
for gene_uri in mapping:	
	
	#json contains raw input of REST
	gene_id = gene_uri[len(gene_id_uri):]
	orthologs_list = [gene_id]
	for orthologs in mapping[gene_uri].values():
		orthologs_list.extend([x[len(gene_id_uri):] for x in orthologs])

	api_file_path = ensembl_path+gene_id+'.json'
	
	if(pre.file_notempty(api_file_path)):
		with open(api_file_path, 'r') as gene_api_file: 
			genetree_json = json.load(gene_api_file)
	else:
		print("redo request")
		
		continue
		
		ext = "/genetree/member/id/"+gene_id+"?sequence=protein&prune_taxon="+species_str
		r = requests.get(ensembl_server+ext, headers={  "Content-Type" : "application/json"})

		if not r.ok:
			with open(log_file, mode='ab') as log:
				log.write(r.content)
			
			redo[gene_uri]=mapping[gene_uri]
			continue
		
		genetree_json = r.json() 
		with open(api_file_path, 'w') as gene_api_file: 
			gene_api_file.write(json.dumps(genetree_json))
		
				
	genetree_id = genetree_json['id']
	og_id = genetree_id+'_'+gene_id
	print(og_id)
	
	genes_to_genetrees[gene_id]=genetree_id
	
	if genetree_id not in genetrees_mapping: genetrees_mapping[og_id] = []
	genetrees_mapping[og_id].append(get_identifiers(genetree_json))
	
with open(redo_path, 'w') as output_redo:		#part of gene tree in NH(x?) format
	output_redo.write(json.dumps(redo))

with open(genes_to_genetrees_file, 'w') as output_g2g:		#part of gene tree in NH(x?) format
	output_g2g.write(json.dumps(genes_to_genetrees))				 

with open(genetrees_mapping_file, 'w') as output_gt:		#part of gene tree in NH(x?) format
	output_gt.write(json.dumps(genetrees_mapping))

		
	
