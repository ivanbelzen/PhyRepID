## IvB April 2018
## Input: orthologs.json
# All human protein coding genes with their orthology relationships with species in pipeline
## Output: orthologs_filtered.json 
# minimum of 10 proteins in group and max 3 orthologs of same species
# should have a mouse ortholog and comply to pairs of species, these groups are defined in ensembl_stable_id_species.json
## Also outputs: excluded_orth_rel.json   excluded_pairs.json
# for human protein coding gene ids which do not comply to these restrictions

import requests, sys, json, os.path
import matplotlib.pyplot as plt
import numpy as np
import pipeline_methods as pre
		
root = '/home/ianthe/protein-repeat-evolution/'
ogs_file = root+'orthologs.json'
species_file = root+'ensembl_stable_id_species.json'

filtered_ogs_file = root+'orthologs_filtered.json'
excluded_orth_rel_file = root + 'excluded_orth_rel.json'
excluded_pairs_file = root + 'excluded_pairs.json'

taxon_uri =  'http://identifiers.org/taxonomy/'
#species_groups = [['13616','9305'],['8839','9031'],['31033','69293']]

species_shared_cutoff = 10
orth_rel_cutoff = 3
number_orthologs_cutoff = 1

#Load orthologs
if(pre.file_notempty(ogs_file)):
	with open(ogs_file,'r') as ogs:
		mapping = json.load(ogs)
else: 
	print ("Requires orthologs.json file")
	quit()


print("Filter OGs, number of OGs at start: ",len(mapping))

species_groups = pre.get_species_group_list(species_file).values()

#should have a mouse ortholog
for gene_uri in list(mapping):
	if taxon_uri+'10090' not in mapping[gene_uri]:
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

if 'http://rdf.ebi.ac.uk/resource/ensembl/ENSG00000164256' in mapping: print('prdm9')

if 'http://rdf.ebi.ac.uk/resource/ensembl/ENSG00000137812' in mapping: print('knl1') 

quit()
with open(excluded_orth_rel_file,'w') as excluded:
	excluded.write(json.dumps(excluded_orth))

with open(excluded_pairs_file,'w') as excluded:
	excluded.write(json.dumps(excluded_pairs))
					
with open(filtered_ogs_file,'w') as filtered_ogs:
	filtered_ogs.write(json.dumps(mapping))	

