## Ivb September 2018
## Analyse output of evo events in dataset
## {og_hit: {nD, nL, root_dup, events = {node = {d,l}}, species = {d,l}

# pfam_evo_events.json / meme_evo_events.json reroots
# repeat { og_id_domain : {netto_dup, loss, events = {dup = {node:cnt, ...} , loss = {node:cnt, ...} } } }, genetree : {}

import sys, json
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import StandardScaler, Normalizer, RobustScaler, Imputer
import pandas as pd
from sklearn.decomposition import PCA

from ete3 import Tree, TreeStyle, NodeStyle, TextFace, PhyloTree, AttrFace, CircleFace, faces

import pipeline_methods as pre

analysis_img_path = pre.root+'analysis/evo_events/'
output_prefix = ''
output_postfix = '_mya'
identifier = 'ENS' #general prefix of OG ids
duplication_absolute = False #default is per mya
genetree = False #species tree by default

species_mapping = pre.get_species_mapping_full(pre.species_mapping_file)

if len(sys.argv) > 1: 
	identifier = sys.argv[1]
	analysis_img_path = pre.root+'analysis/maptotree/'
	output_prefix = identifier+'_'
	duplication_absolute = True
	genetree = True
	output_postfix = '_abs'
else:
	quit()

pfam_results_file = 'pfam_evo_events_reroot.json'
meme_results_file = 'meme_evo_events_reroot.json'
iqtree_results_file = 'iqtree_evo_events.json'

with open(pfam_results_file, 'r') as results: pfam_results = json.load(results)
with open(meme_results_file, 'r') as results: meme_results = json.load(results)
#with open(iqtree_results_file, 'r') as results: iqtree_results = json.load(results)

pfam_repeat_data = pfam_results['repeat']
meme_repeat_data = meme_results['repeat']
#iqtree_repeat_data = iqtree_results

node_projection = { 't32525': ['t9347', 't1437010', 't314147'], 't8287': ['t32523'], 't186625':['t186626'], 't41665': ['t1489872'], 't117571':['t7742', 't7711', 't33213']}
root_nodes = ['t117571','t7742', 't7711', 't33213'] #bony vertebrates and before

map_to_species_tree = {}

#PFAM
for og_id_domain,data in pfam_repeat_data.items():
	if identifier not in og_id_domain: continue 
	if not 'dup' in data['events']: continue
	
	for node,cnt in data['events']['dup'].items():
		if genetree is False:
			node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))[:-1]
			if node == 'ENS': node = 'ENSP'
	
		if node not in map_to_species_tree:
			map_to_species_tree[node] = {'d':0,'l':0}
				
		map_to_species_tree[node]['d'] += cnt
	
#MEME
for og_id,data in meme_repeat_data.items(): 
	if identifier not in og_id: continue 
	if not 'dup' in data['events']: continue
	
	for node,cnt in data['events']['dup'].items():
		if genetree is False:	
			node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))[:-1]
			if node == 'ENS': node = 'ENSP'
	
		if node not in map_to_species_tree:
			map_to_species_tree[node] = {'d':0,'l':0}
				
		map_to_species_tree[node]['d'] += cnt

total_dups = sum([int(val['d']) for x,val in map_to_species_tree.items()])
netto_dup = total_dups-sum([int(val['d']) for x,val in map_to_species_tree.items() if x in root_nodes])


print('Total',total_dups,'netto',netto_dup)

'''
#IQTREE
print(len( iqtree_repeat_data))
map_to_species_tree = {}

for og_id,data in iqtree_repeat_data.items(): 
	#if 'ENSGT00410000025918_ENSG00000137812' not in og_id: continue 
	if not 'events' in data: continue
	if not 'dup' in data['events']: continue
	
	for node,cnt in data['events']['dup'].items():
		node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))[:-1]
		if node == 'ENS': node = 'ENSP'
	
		if node not in map_to_species_tree:
			map_to_species_tree[node] = {'d':0,'l':0}
				
		map_to_species_tree[node]['d'] += cnt
		
'''

rates_dict = {}

#print(total_dups)
total_dups = sum([int(val['d']) for x,val in map_to_species_tree.items()])
netto_dup = total_dups-sum([int(val['d']) for x,val in map_to_species_tree.items() if x in root_nodes])


print('Total',total_dups,'netto',netto_dup)

#quit()

if genetree is False:
	species_tree_file = pre.root+'ensembl_taxon_species_tree_mya.phy'
	with open(species_tree_file, 'r') as species_tree_string:
		species_tree = PhyloTree(species_tree_string.read(), format=1)
	
	for x in root_nodes: 
		if x in map_to_species_tree: del map_to_species_tree[x];
		
if genetree is True:
	og_id = '_'.join(identifier.split('_')[0:2])
	def parse_sp_name(node_name): return filter(str.isalpha, node_name.split('_')[0])[:-1] # return node_name.split('_')[0]
	genetree_file = pre.genetree_path+og_id+'.nhx'
	#genetree_file = pre.pfam_treefix_path+og_id+'.stree'
	with open(genetree_file, 'r') as gene_tree_string:
		gene_tree_string = gene_tree_string.read()
		gene_tree_string = pre.make_valid_nhx(gene_tree_string)
	species_tree = PhyloTree(gene_tree_string, format=1, sp_naming_function=parse_sp_name)
	


for node in species_tree.traverse():
	duplications = 0
	if node.name in map_to_species_tree:
		duplications += map_to_species_tree[node.name]['d']
		del(map_to_species_tree[node.name])
	
	if node.name in node_projection:
		for taxa in node_projection[node.name]:
			if taxa in map_to_species_tree:
				duplications += map_to_species_tree[taxa]['d']
				del(map_to_species_tree[taxa])
		
	if any(node.name.startswith(ensembl_id[:-1]) for ensembl_id in species_mapping.keys())  and genetree is False:
		if node.name == 'ENSP': node.name = 'ENS'
		taxon_uri = 't'+str(species_mapping[node.name+'P']['taxon_uri'])
		ns = NodeStyle()
		ns['fgcolor'] = species_mapping[node.name+'P']['colour']
		ns['size'] = 12
		node.name = species_mapping[node.name+'P']['common_name']
		node.set_style(ns)
		
		
		if taxon_uri in map_to_species_tree:
			duplications += map_to_species_tree[taxon_uri]['d']
			del(map_to_species_tree[taxon_uri])
	
	if genetree is True and hasattr(node,'species') and node.species in map_to_species_tree:
		duplications += map_to_species_tree[node.species]['d']
		del(map_to_species_tree[node.species])
	
	if duplications > 0:
		relative_dup = round(float(duplications)/netto_dup,2)
		if node.dist > 0:
			timed_dup = round(float(duplications)/(float(node.dist)/2),2)
		else: timed_dup = 0
		rates_dict[node.name] = {}
		rates_dict[node.name]['fraction'] = relative_dup
		rates_dict[node.name]['timed'] = timed_dup
		rates_dict[node.name]['absolute'] = duplications
		
		if duplication_absolute is False:
			node.add_feature('dup', timed_dup)	
		elif duplication_absolute is True:
			node.add_feature('dup', duplications)
						
def duploss(node):
	if hasattr(node, 'dup'):
		node.add_face(AttrFace('dup', fsize=5, fgcolor='blue'), column=0, position="branch-bottom")
	if hasattr(node, 'loss'):
		node.add_face(AttrFace('loss', fsize=5, fgcolor='red'), column=0, position="branch-bottom")
 
print('unmapped nodes',map_to_species_tree)


ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False
ts.layout_fn = [duploss]

#species_tree.show(tree_style=ts)

img_name = analysis_img_path+output_prefix+"maptotree"+output_postfix+".pdf"
species_tree.render(img_name, tree_style = ts)

outfile_name = analysis_img_path+output_prefix+"maptotree"+output_postfix
species_tree.write(features=['dup'], outfile=outfile_name+".nhx", format_root_node=True)

for node in species_tree.traverse():
	if hasattr(node, 'dup'):
		node.dup = np.log(node.dup)
		
species_tree.write(features=['dup'], outfile=outfile_name+"_log.nhx", format_root_node=True)

