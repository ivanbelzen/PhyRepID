## IvB April 2018
## Input: 
# pre.orthologs_file - orthologs.json   generated by retrieve_orthologs.py
# All human protein coding genes with their orthology relationships with species in pipeline
## Output: 
# pre.orthologs_filtered_file - orthologs_filtered.json 
# minimum of 10 proteins in group and max 3 orthologs of same species
# should have a mouse ortholog and comply to pairs of species, these groups are defined in ensembl_stable_id_species.json
## Also outputs: excluded_orth_rel.json   excluded_pairs.json
# for human protein coding gene ids which do not comply to these restrictions

import requests, sys, json, os.path
import numpy as np
import pipeline_methods as pre

species_shared_cutoff = 10
orth_rel_cutoff = 3
number_orthologs_cutoff = 1

#Load orthologs
if(pre.file_notempty(pre.orthologs_file)):
	with open(pre.orthologs_file,'r') as ogs:
		mapping = json.load(ogs)
else: 
	print ("Requires ",pre.orthologs_file," file")
	quit()


print("Filter OGs, number of OGs at start: ",len(mapping))

species_groups = pre.get_species_group_list(pre.species_mapping_file).values()

#should have a mouse ortholog
for gene_uri in list(mapping):
	if pre.taxon_uri+'10090' not in mapping[gene_uri]:
		mapping.pop(gene_uri)

print("#OGs with mouse orth: ",len(mapping))

(mapping,excluded_pairs) = pre.filter_close_pairs(mapping, species_groups, True)
print("#OGs with close pair restriction: ",len(mapping), "excluded :", len(excluded_pairs))

mapping = pre.filter_min_proteins(mapping,number_orthologs_cutoff)
print("#OGs with at least "+str(number_orthologs_cutoff)+" sequences: ",len(mapping))

(mapping,excluded_orth) = pre.filter_orth_rel(mapping, orth_rel_cutoff, True)
mapping = pre.filter_close_pairs(mapping, species_groups)
mapping = pre.filter_min_proteins(mapping, number_orthologs_cutoff)		
print("#OGs with max "+str(orth_rel_cutoff)+" orth per species: " ,len(mapping), "excluded: ",len(excluded_orth))

with open(pre.excluded_orth_rel_file,'w') as excluded:
	excluded.write(json.dumps(excluded_orth))

with open(pre.excluded_pairs_file,'w') as excluded:
	excluded.write(json.dumps(excluded_pairs))
					
with open(pre.filtered_ogs_file,'w') as filtered_ogs:
	filtered_ogs.write(json.dumps(mapping))	

