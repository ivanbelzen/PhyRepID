# IvB extract conserved/separated TRs 
# Edit 27-05-2018 
# 

import sys
import pipeline_methods as pre
import json

#pairwise_flag = True if (len(sys.argv) > 1 and sys.argv[1] == 'pairwise') else False
#summary_flag = True  if (len(sys.argv) > 1 and sys.argv[1] == 'summary') else False

#files = ['eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick', 'eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick']
#schaper_files = [pre.root+'schaper/eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick']
schaper_files = [pre.root+'schaper/eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick']

#summary_file = 'schaper_summary.json'
summary_file = 'schaper_summary_denovo.json'
#pfam_summary_file = 'schaper_pfam_summary.json'

species_dict = {} #pre.species_file full dict
ensembl_stable_id_list = []

pfam_accession = {} #maps pfam accession ID to hit name
clans = {} #pre function mapping of hit to clan
 
schaper_summary = {} # human ortholog clan pfam_hit relationship
#pfam_summary = {} #clan hit relationship #count

'''
all_stable_ids = [] #in Schaper results
protein_not_found = [] #unassociated with gene id atm
ensembl_stable_id_list = [] #ids from ensembl, determines what gets included/recognised
ensembl_mapping = {}
pfam_mapping = {}
pfam_not_found = [] #unassociated with pfam accession
summary = {}
'''

with open(pre.species_file,'r') as sp_mapping:
	species_dict = json.load(sp_mapping)
	for item in species_dict:
		ensembl_stable_id_list.append(item['ensembl_stable_id'])					

'''
clans = pre.get_pfam_clans(pre.pfam_clans_file) #pfam hit to clan mapping

#pfam accession -> hit mapping
with open(pre.pfam_clans_file, 'r') as pfam_file: 
	for rows in pfam_file:
		cols = rows.strip().split('\t')
		pfam_accession[cols[0]] = cols[3]	
'''
#if pairwise_flag: pairwise_output = open(pairwise_file,'w')

for filename in schaper_files:
	with open(filename, 'r') as f:
		for line in f:
			if line[0] != '>': continue
			param = {}
				
			for col in line.strip().split(' '):
				key,value = col.split(':')
				param[key] = value
				
			#ignore all inconclusive entries		
			if 'Pairwise_phylogeny_type' not in param: continue
			
			phylogeny_type = param['Pairwise_phylogeny_type']
			human_protein_id = param['Ensembl_Protein_ID_Human_Ortholog']
			ortholog = param['Ensembl_Protein_ID_Second_Ortholog']
			identifier = ortholog.split('0') #assumption correct that always starts with 0??
			identifier = identifier[0][:-1] #only need first part (species) and remove trailing P as indication for protein

			'''
			pfam_id = param['TR_detection_type']
			pfam_name = pfam_accession[pfam_id] if pfam_id in pfam_accession else 'denovo'
			pfam_clan = clans[pfam_name] if pfam_name in clans else 'denovo'
			'''
			#species should be in our list					
			if identifier not in ensembl_stable_id_list: continue
			
			#alternative is None which is less informative	
			#if not ('conserved' in phylogeny_type or 'separated' in phylogeny_type): continue
			
			if human_protein_id not in schaper_summary:
				schaper_summary[human_protein_id] = {}
			if ortholog not in schaper_summary[human_protein_id]:
				schaper_summary[human_protein_id][ortholog] = {}
			schaper_summary[human_protein_id][ortholog]=phylogeny_type
			'''
			if pfam_clan not in schaper_summary[human_protein_id][ortholog]:
				schaper_summary[human_protein_id][ortholog][pfam_clan] = {}
			schaper_summary[human_protein_id][ortholog][pfam_clan][pfam_name]=phylogeny_type
			'''				
			
			
			'''
			if pfam_clan not in pfam_summary:
				pfam_summary[pfam_clan]={}
			if pfam_name not in pfam_summary[pfam_clan]:
				pfam_summary[pfam_clan][pfam_name] = {}
			if phylogeny_type not in pfam_summary[pfam_clan][pfam_name]:
				pfam_summary[pfam_clan][pfam_name][phylogeny_type] = 0
			pfam_summary[pfam_clan][pfam_name][phylogeny_type]+=1
			'''
			

with open(summary_file, 'w') as summary_file: 	
	summary_file.write(json.dumps(schaper_summary))

'''	
with open(pfam_summary_file, 'w') as pfam_summary_file: 	
	pfam_summary_file.write(json.dumps(pfam_summary))
'''	

			
