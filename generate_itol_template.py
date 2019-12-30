## IvB September 2018
# Generate domain templates for easy visualisation in ITOL
# using pfam summary.json and meme_summary.json

import sys, os, glob, json
import numpy as np
import pipeline_methods as pre

if len(sys.argv) < 2: 
	print 'Please provide an OG ID'
	quit()
	
og_id = sys.argv[1] 

itol_template_path = pre.root+'analysis/itol/'
	
shape = 'RE'
colour = '#00aedb'
padding = 5

with open(pre.root+"pfam_repeat_stats.json", 'r') as output: pfam_summary = json.load(output)
with open(pre.root+"meme_repeat_stats.json", 'r') as output: meme_summary = json.load(output)
	
pfam_initial_filtered_output_file=pre.root+"pfam_repeat_stats_initial_filtered.json"

def generate_domain_string(orthologs_dict,label,colour='#0099ff'):
	template_string = ''
	for protein_id,orth_val in orthologs_dict.items():
		template_string += protein_id+','+str(orth_val['length'])+','
			
		units_list = [int(x) for x in orth_val['units_dict'].keys()]
		units_list.sort()
		for units_id in units_list:
			units_val=orth_val['units_dict'][str(units_id)]
			b,e = units_val['seq_coordinates']
			b = int(b)-padding
			e = int(e)+padding
			template_string += shape+'|'+str(b)+'|'+str(e)+'|'+colour+'|'+label+','
			
		template_string = template_string[:-1]+'\n'	
	return(template_string)


def generate_label_string_repeats(orthologs_dict,label):
	species_mapping = pre.get_species_mapping_full(pre.species_mapping_file)
	
	template_string = ''
	for protein_id,orth_val in orthologs_dict.items():
				
		units_list = [int(x) for x in orth_val['units_dict'].keys()]
		units_list.sort()
		for units_id in units_list:
			units_val=orth_val['units_dict'][str(units_id)]
			b,e = units_val['seq_coordinates']
			b = int(b)-padding
			e = min(int(e)+padding,orth_val['length'])
			ensembl_id = filter(str.isalpha,str(protein_id))
			if ensembl_id in species_mapping: 
				name = species_mapping[ensembl_id]['common_name']
			else: name = ensembl_id
			
			template_string += protein_id+'_'+label+'_'+str(units_id)+"_"+str(b)+'-'+str(e)+','+name+'_'+str(units_id)+'\n'	
			
	return(template_string)

def generate_label_string_genetree(orthologs_dict):
	species_mapping = pre.get_species_mapping_full(pre.species_mapping_file)
	
	template_string = ''
	for protein_id,orth_val in orthologs_dict.items():			
		ensembl_id = filter(str.isalpha,str(protein_id))
		if ensembl_id in species_mapping: 
			name = species_mapping[ensembl_id]['common_name']
		else: name = ensembl_id
			
		template_string += protein_id+','+name+'\n'	
			
	return(template_string)


if not pre.file_notempty(itol_template_path+og_id+'_domains.txt'):
	template_string = 'DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL,domains\nCOLOR,#ff0000\nDATA\n'

	if og_id in pfam_summary:
		for domain, dom_val in pfam_summary[og_id].items():
			label = domain
			template_string += generate_domain_string(dom_val['orthologs_dict'], label)

	if og_id in meme_summary:
		label = ''
		template_string += generate_domain_string(meme_summary[og_id]['motif']['orthologs_dict'], label)			

	if og_id not in pfam_summary and og_id not in meme_summary: print("Unable to find ",og_id)
	
	with open(itol_template_path+og_id+'_domains.txt','w') as output:
		output.write(template_string)

else:
	print('File already exists for',og_id)

if pre.file_notempty(itol_template_path+og_id+'_domains_initial.txt'):
	print('File initial already exists for',og_id)


elif len(sys.argv)>2 and 'initial' in sys.argv[2]:
	with open(pfam_initial_filtered_output_file, 'r') as output:
		pfam_summary_initial = json.load(output)

	real_og_id = og_id #need the og_id_domain for the filtered mapping and otherwise real og id (without domain)
	colour_list = ['#4cc88a']
	template_string = 'DATASET_DOMAINS\nSEPARATOR COMMA\nDATASET_LABEL,initial\nCOLOR,#4cc88a\nDATA\n'
	if real_og_id in pfam_summary_initial:
		print(pfam_summary_initial[real_og_id])
		
		i=0
		for domain, dom_val in pfam_summary_initial[real_og_id].items():
			colour_temp = colour_list[i% len(colour_list)]
			label = domain
			template_string += generate_domain_string(dom_val['orthologs_dict'], label, colour_temp)
			i+=1
	with open(itol_template_path+og_id+'_domains_initial.txt','w') as output:
		output.write(template_string)

if not pre.file_notempty(itol_template_path+og_id+'_repeat_labels.txt'):
	template_string = 'LABELS\nSEPARATOR COMMA\nDATA\n'
	
	if og_id in pfam_summary:
		for domain, dom_val in pfam_summary[og_id].items():
			label = domain#.replace('_',' ')
			template_string += generate_label_string_repeats(dom_val['orthologs_dict'], label)

	if og_id in meme_summary:
		label = ''
		template_string += generate_label_string_repeats(meme_summary[og_id]['motif']['orthologs_dict'], label)			

	if og_id not in pfam_summary and og_id not in meme_summary: print("Unable to find ",og_id)
	
	with open(itol_template_path+og_id+'_repeat_labels.txt','w') as output:
		output.write(template_string)
else:
	print('Repeat labels file already exists for',og_id)

if not pre.file_notempty(itol_template_path+og_id+'_genetree_labels.txt'):
	template_string = 'LABELS\nSEPARATOR COMMA\nDATA\n'
	

	if og_id in pfam_summary:
		for domain, dom_val in pfam_summary[og_id].items():
			template_string += generate_label_string_genetree(dom_val['orthologs_dict'])

	if og_id in meme_summary:
		template_string += generate_label_string_genetree(meme_summary[og_id]['motif']['orthologs_dict'])

	if og_id not in pfam_summary and og_id not in meme_summary: print("Unable to find ",og_id)
	
	with open(itol_template_path+og_id+'_genetree_labels.txt','w') as output:
		output.write(template_string)

else:
	print('Genetree labels file already exists for',og_id)
