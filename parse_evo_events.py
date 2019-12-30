## IvB August 2018
# Last edit 30-12-2019
# Parse Treefix by removing inconsistent duplications and counting dup/loss 
# PFAM and MEME version => differences in paths and og ids

# input: .treefix.nhx.tree (output of Treefix Annotate) and .nhx gene tree and ensembl species tree 
# output: adjusted repeat tree and image - paths defined in pipeline methods.py
# pre.pfam_treefix_adjusted_path and pre.pfam_treefix_adjusted_images_path
# pre.denovo_treefix_adjusted_path and pre.denovo_treefix_adjusted_images_path

# output dict rearranged: log_pfam_rearranged_nodes.json log_meme_rearranged_nodes.json
# output dict dup/loss repeat tree vs gene tree: pfam_evo_events.json meme_evo_events.json
# output dict dup/loss gene tree vs species tree: pfam_evo_events_genetrees.json meme_evo_events_genetrees.json

## Inconsistent duplication correction
# Remove spurious duplications from repeat tree that have a consistency score below the threshold
# Default threshold = 0
# Duplication consistency depends on number of orthologs (count from gene tree)
# Prune/graft gene tree on repeat-tree parts with inconsistent duplications

## Basic algorithm
# Count duplications/losses using tree reconciliation with gene tree
# Count dup/loss in gene tree compared to species tree => 
# Calculate number of netto duplications (after vertebrate ancestor as defined in gene tree)

## Notes:
# Dup threshold >1 does not work well it seems (27-08-18)
# Note layout species coloring changes name so should be left on??

from ete3 import Tree, TreeStyle, NodeStyle, TextFace, PhyloTree, AttrFace, CircleFace, faces
import sys, glob, json 
import os.path
import pipeline_methods as pre
import numpy as np


## Settings
dupconsistency_threshold = 0
root_nodes = ['t117571','t7742', 't7711', 't33213']
node_projection = { 't32525': ['t9347', 't1437010', 't314147'], 't8287': ['t32523'], 't186625':['t186626'], 't41665': ['t1489872'], 't117571':['t7742', 't7711', 't33213']}
colour_list = ['lightgrey','dimgrey','darkgrey','rosybrown','sienna','powderblue','olive'] #used for coloring nodes
species_mapping = pre.get_species_mapping_full(pre.species_mapping_file)
protein_dict = {} #for colours and naming


## Define functions for parsing ETE trees

#parsing species names for different input trees
def parse_sp_name(node_name):
	return node_name.split('_')[0]

def parse_genetree_sp(node_name):
	return filter(str.isalpha, node_name.split('_')[0])

def parse_ensembl_sp(node_name):
	if node_name == 'ENSP': return node_name
	return node_name+'P'    

#functions for parsing events and reconciled trees
def find_taxon_id(node_list):
	for node in node_list:
		if node.name != '': 
			return(node.name.split('_')[0]) #only protein id or taxon id, not lead number and location
		else: 
			return(find_taxon_id(node.get_children()))

def loss_leaf(node):
	#if you are loss node and only have children with loss
	if getattr(node,"evoltype", None) != "L":
		return(False) #no loss node
	elif not node.is_leaf(): #need to check if both children are only loss nodes as well
		[child1, child2] = node.get_children()
		return(loss_leaf(child1) and loss_leaf(child2))
	else: #is loss node and leaf node
		return(True)

def get_evo_events(recon_tree,root=''):			
	duplications=[];losses=[]
	if root == '': root = 't117571'
	
	for node in recon_tree.traverse():
		if getattr(node,"evoltype", None) == "D":
			duplications.append(find_taxon_id([node]))
		
	for leaf in recon_tree.iter_leaves(is_leaf_fn=loss_leaf):
		losses.append(find_taxon_id([leaf]))
	
	values,counts = np.unique(duplications, return_counts=True)
	duplication_events = dict(zip(values,counts))
	values,counts = np.unique(losses, return_counts=True)
	loss_events = dict(zip(values,counts))
	
	total_dup = sum(duplication_events.values())
	
	root_list = [root]
	if root in node_projection: root_list.extend(node_projection[root])
	
	netto_dup = total_dup 
	for node in root_list:
		if node in duplication_events:
			netto_dup -= duplication_events[node]
		
	loss = sum(loss_events.values())

	return({'total_dup':total_dup, 'netto_dup': netto_dup, 'loss':loss, 'events':{'dup':duplication_events, 'loss':loss_events}})

#layout function to colour and format nodes
def species_colouring(node):
	
	if node.is_leaf() and not node.name.isdigit() and not node.name is '' and not getattr(node,"evoltype", None) == "L":
		#node.name [Protein_id]_[hit (can contain underscores) ]_[number]_[position]
		frag = node.name.split('_')
		protein_id = frag[0]
		species_id = ''.join(filter(str.isalpha, protein_id))
		ns = NodeStyle()
		#if protein_id in protein_dict:
		#	protein_face = CircleFace(8, protein_dict[protein_id], style='circle', label=None)
		#	faces.add_face_to_node(protein_face, node, column=1)
		if species_id in species_mapping:
			ns['fgcolor'] = species_mapping[species_id]['colour']
			ns['size'] = 12
		node.name = protein_id+' '+frag[-2]+' ('+frag[-1]+')' #last and second to last elements.
		name_face = AttrFace("name", fsize=8)
		faces.add_face_to_node(name_face, node, column=0, position="branch-right")
		node.set_style(ns)
			
def duplication_colouring(node):
	if getattr(node,"evoltype", None) == "D":
		ns = NodeStyle()	
		ns["shape"] = "square"
		ns['size'] = 8
		if find_taxon_id([node]) in root_nodes:
			ns['fgcolor'] = 'grey' 
		else:
			ns['fgcolor'] = 'blue' 
			#name_face = TextFace(find_taxon_id([node]), fsize=5)
			#faces.add_face_to_node(name_face, node, column=0, position="branch-right")
		node.set_style(ns)	
		

## Reconcile gene tree with species tree 
# do not merge with other reconciliation loop to prevent duplications or overwriting in dict
# need to reload and parse genetree with other species names to reconcile with Ensembl species tree

def reconcile_gene_species_tree():
	for treefix_file in treefix_file_list:
	og_id_domain = treefix_file[len(treefix_path):-len(treefix_ext)]
	og_id = '_'.join(og_id_domain.split('_')[0:2])
	
	genetree_file = pre.genetree_path+og_id+genetree_ext
	if not pre.file_notempty(genetree_file):
		continue
	
	with open(genetree_file, 'r') as gene_tree_string:
		gene_tree_string = gene_tree_string.read()
		gene_tree_string = pre.make_valid_nhx(gene_tree_string)
	gene_tree = PhyloTree(gene_tree_string, format=1, sp_naming_function=parse_genetree_sp)
	
	evo_events_log['genetree'][og_id]={}
		
	with open(pre.species_tree_file, 'r') as species_tree_string:
		species_tree_string = species_tree_string.read()
	species_tree = PhyloTree(species_tree_string, format=1, sp_naming_function=parse_ensembl_sp)
		
	for node in species_tree.traverse():
		node.name = parse_ensembl_sp(node.name)
		
	recon_genetree, events = gene_tree.reconcile(species_tree)
	events = get_evo_events(recon_genetree) #define myself instead of using ETE events
	evo_events_log['genetree'][og_id] = events
	evo_events_log['genetree'][og_id]['root']=gene_tree.name

	with open(evo_events_genetrees_log_file, 'w') as log:
		log.write(json.dumps(evo_events_log))

## Reconcile repeat tree with gene tree

def reconcile_repeat_gene_tree():
	for treefix_file in treefix_file_list:
	og_id_domain = treefix_file[len(treefix_path):-len(treefix_ext)]
	og_id = '_'.join(og_id_domain.split('_')[0:2])
		
	genetree_file = pre.genetree_path+og_id+genetree_ext
	if not pre.file_notempty(genetree_file):
		continue
	
	print(og_id_domain)
	
	rearranged_log[og_id_domain] = {'affected_nodes':[]}
	evo_events_log['repeat'][og_id_domain]={}

	# open Files
	with open(treefix_file, 'r') as repeat_tree_string:
		repeat_tree_string = repeat_tree_string.read()
	tree = PhyloTree(repeat_tree_string, format=1, sp_naming_function=parse_sp_name)
	tree.ladderize()

	with open(genetree_file, 'r') as gene_tree_string:
		gene_tree_string = gene_tree_string.read()
		gene_tree_string = pre.make_valid_nhx(gene_tree_string)
	gene_tree = PhyloTree(gene_tree_string, format=1, sp_naming_function=parse_sp_name)
	
	
	i=0
	for l in gene_tree.get_leaves():
		protein_dict[l.name] = colour_list[i% len(colour_list)]
		i+=1
	
	## Adjustment algorithm for inconsistent duplications in repeat tree
	for node in tree.traverse():
		node.dist=1
		if hasattr(node, 'T'): print(node.T)
		if not node.is_leaf() and hasattr(node, 'D') and node.D == 'Y':
			children = node.get_children()
			c1 = children[0]
			c2 = children[1]
		
			consistency = None
			intersection = None
			total_species = []
				
			if c1.is_leaf() and c2.is_leaf():
				if c1.S == c2.S: 
					consistency = 1; intersection = c1.S; total_species.append(c1.S)
				else: 
					consistency = 0; intersection = None; total_species.append(c1.S); total_species.append(c2.S)
			
			elif c2.is_leaf(): 
				print('yeah, c2 can be leaf')
				#does this ever happen???
				if hasattr(c1, 'intersection'):
					if c1.intersection == c2.S: 
						consistency = 1; intersection = c1.intersection; #totalspecies already has been appended
					else:
						consistency = 0; intersection = None; total_species.append(c2.S)
				else: 
					print( c1+'no intersection attr')
				
			else: 
				species_c1 = []
				for l in c1.get_leaves():
					sp = l.name.split('_')[0]
					if sp not in species_c1: species_c1.append(sp)
				
				species_c2 = []
				for l in c2.get_leaves():
					sp = l.name.split('_')[0]
					if sp not in species_c2: species_c2.append(sp)
							
				intersection =  [v for v in species_c1 if v in species_c2]	
				consistency = round(float(len(intersection))/max(len(species_c1),len(species_c2)),3)
				total_species += species_c1 
				total_species += species_c2
		
			if intersection is not None:
				node.add_feature('intersection', intersection)
			if consistency is not None:
				node.add_feature('consistency', consistency)	
			
			if node.consistency <= dupconsistency_threshold and len(children) < 8: #only rearrange terminal dups
				#print( c1, c2, node )
				
				rearranged = gene_tree.copy()
				rearranged.ladderize()
				all_nodes = []
				
				for n in rearranged.traverse():
					n.dist = 1
					if n.name in total_species:
						overlap_node = node.search_nodes(species=n.name)[0]
						n.name = overlap_node.name
						all_nodes.append(n.name)
							
				rearranged.prune(all_nodes)
				
				rearranged_log[og_id_domain]['affected_nodes'].append(all_nodes)
				
				ns = NodeStyle()	
				ns["shape"] = "square"
				ns['size'] = 5	
				ns['fgcolor'] = 'blue'
				rearranged.add_face(TextFace('rd', fsize=5), column=0, position="branch-bottom")
				rearranged.set_style(ns)
				
				node.add_sister(rearranged)
				node.detach()
				node.delete()

	##Parsing evolutionary events
	#Repeat tree events
	recon_tree, events = tree.reconcile(gene_tree) 
	events = get_evo_events(recon_tree,gene_tree.name) #define myself instead of using ETE events
	evo_events_log['repeat'][og_id_domain] = events
	
	#output adjusted repeat tree
	if output_adjusted_repeat_tree is True:	
		recon_tree.write(format=1, features=['evoltype','T'], outfile=output_adjusted_path+og_id_domain+'.ete')
	
	if output_image is True:
		ts = TreeStyle()
		ts.show_leaf_name = False
		ts.show_branch_length = False
		ts.show_branch_support = False
		ts.complete_branch_lines_when_necessary = False
		ts.title.add_face(TextFace(og_id_domain, fsize=18), column=0)
		ts.layout_fn = [species_colouring, duplication_colouring]
		#recon_tree.show( tree_style = ts)
		img_name = output_images_path+og_id_domain+".pdf"
		recon_tree.render(img_name, tree_style = ts)
	
	#Write logs
	with open(rearranged_nodes_log_file, 'w') as log:
		log.write(json.dumps(rearranged_log))
	
	with open(evo_events_log_file, 'w') as log:
		log.write(json.dumps(evo_events_log))


##Settings for PFAM and MEME
output_images_path = pre.pfam_treefix_adjusted_images_path
output_adjusted_path = pre.pfam_treefix_adjusted_path
treefix_path = pre.pfam_treefix_path

output_image = True #Outputs tree as PDF, True for debugging and visualisation, False for production purposes
output_adjusted_repeat_tree = True

rearranged_nodes_log_file = pre.root+'log_pfam_rearranged_nodes.json'
rearranged_log = {'settings':{'threshold':dupconsistency_threshold}}

evo_events_log_file = pre.root+'pfam_evo_events.json'
evo_events_genetrees_log_file = pre.root+'pfam_evo_events_genetrees.json'
evo_events_log = {'repeat':{}, 'genetree':{}}

treefix_ext = '.treefix.nhx.tree'
genetree_ext = '.nhx'
genetree_file_list = glob.glob(pre.genetree_path + '*' + genetree_ext )
treefix_file_list = glob.glob(treefix_path + '*' + treefix_ext ) 		

reconcile_gene_species_tree()
reconcile_repeat_gene_tree()

## MEME

rearranged_nodes_log_file = pre.root+'log_meme_rearranged_nodes.json'
rearranged_log = {'settings':{'threshold':dupconsistency_threshold}}

evo_events_log_file = pre.root+'meme_evo_events.json'
evo_events_genetrees_log_file = pre.root+'meme_evo_events_genetrees.json'
evo_events_log = {'repeat':{}, 'genetree':{}}

output_images_path = pre.denovo_treefix_adjusted_images_path
output_adjusted_path = pre.denovo_treefix_adjusted_path
treefix_path = pre.denovo_treefix_path

reconcile_gene_species_tree()
reconcile_repeat_gene_tree()

## Finished
