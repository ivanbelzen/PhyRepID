## IvB April 2018
## Parse gene tree retrieved from API by retrieve_genetrees.py
## The gene trees can be seen as orthologuous groups, containing both orthologs and paralogs
## Gene tree json is saved as NHX tree, and to a mapping file
## Gene tree identifier: [ (gene_id, protein_id)  of all from subtree }
## Also writes fasta file for orthoguous groupswaaro

import requests, sys, json, os.path
import pipeline_methods as pre

ensembl_path = pre.ensembl_path
genetree_path = pre.genetree_path
fasta_path = pre.fasta_path
filtered_ogs_file = pre.orthologs_filtered_file
log_file = pre.log_parse_genetree_file

gene_uri_prefix = 'http://rdf.ebi.ac.uk/resource/ensembl/' 
gene_id = '' #human gene id
genetree_json_file = sys.argv[1] #api output

if ensembl_path in genetree_json_file and pre.file_notempty(genetree_json_file):
	gene_id = genetree_json_file[len(ensembl_path):-len('.json')]
	with open(genetree_json_file, 'r') as gt_file: 
		genetree_json = json.load(gt_file)
else:
	with open(log_file, mode='ab') as log:
		log.write("Could not parse: "+genetree_json_file+"\n")
	quit()
	
if(pre.file_notempty(filtered_ogs_file)):
	with open(filtered_ogs_file,'r') as filtered_ogs:
		mapping = json.load(filtered_ogs)
else: 	
	quit()	


#Extract subtree containing all orthologs
#write to NHX tree and to fasta file
orthologs_list = [gene_id]
for orthologs in mapping[gene_uri_prefix+gene_id].values():
	orthologs_list.extend([x[len(gene_uri_prefix):] for x in orthologs])

genetree_id = genetree_json['id']	
og_id = genetree_id+'_'+gene_id
nhx_file_path = genetree_path+og_id+'.nhx'
fasta_file_path = fasta_path+og_id+'.fa'	
	
print(og_id)
	
genetree_json = pre.get_subtrees(genetree_json, orthologs_list)
if genetree_json is None:
	with open(log_file, mode='ab') as log:
		log.write("No subtree: "+og_id+'\n')
	quit()
	
nhx_tree = pre.make_valid_nhx(pre.get_nhxtree(genetree_json))

if nhx_tree is None:
	with open(log_file, mode='ab') as log:
		log.write("No NHX tree: "+og_id+'\n')
	quit()
			
with open(nhx_file_path,'w') as genetree_file:
	genetree_file.write(nhx_tree)


with open(fasta_file_path, 'w') as output_fa:	 #fasta file with all sequences
	output_fa.write(pre.get_proteinseq(genetree_json))

