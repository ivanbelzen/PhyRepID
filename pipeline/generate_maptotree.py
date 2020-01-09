# IvB 
# September 2018
# Last edit 02-01-2020
# Project duplications on gene tree, from inferred evo events resulting from reconciliation of repeat tree with gene tree
# Input: 
# - $1 identifier (og_id_hit) with gene tree 
# - pfam_evo_events.json / meme_evo_events.json 
# Output: 
# - /maptotree/{identifier}_maptotree.pdf
# - /maptotree/{identifier}_maptotree.nhx

import sys, json
import numpy as np
import pandas as pd

from ete3 import Tree, TreeStyle, NodeStyle, TextFace, PhyloTree, AttrFace, CircleFace, faces

import pipeline_methods as pre

if len(sys.argv) < 1:
	print("Need OG-id as input")
	quit()

identifier = sys.argv[1]

img_path = pre.maptotree_path
species_mapping = pre.get_species_mapping_full(pre.species_mapping_file)
root_nodes = pre.root_nodes 
node_projection = pre.node_projection
map_to_species_tree = {}

# Load data
pfam_results_file = pre.pfam_evo_events_file
meme_results_file = pre.meme_evo_events_file

with open(pfam_results_file, 'r') as results: pfam_results = json.load(results)
with open(meme_results_file, 'r') as results: meme_results = json.load(results)


# Load gene tree
og_id = '_'.join(identifier.split('_')[0:2])
def parse_sp_name(node_name): return filter(str.isalpha, node_name.split('_')[0])[:-1] # return node_name.split('_')[0]
genetree_file = pre.genetree_path+og_id+'.nhx'

with open(genetree_file, 'r') as gene_tree_string:
	gene_tree_string = gene_tree_string.read()
	gene_tree_string = pre.make_valid_nhx(gene_tree_string)
gene_tree = PhyloTree(gene_tree_string, format=1, sp_naming_function=parse_sp_name)
	
## Load evo events
pfam_repeat_data = pfam_results['repeat']
meme_repeat_data = meme_results['repeat']
	
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

## Project on gene tree

for node in gene_tree.traverse():
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
 
## Export 

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False
ts.layout_fn = [duploss]

gene_tree.render(img_path+identifier+"_maptotree.pdf", tree_style = ts)
gene_tree.write(features=['dup'], outfile=img_path+identifier+"_maptotree.nhx", format_root_node=True)
