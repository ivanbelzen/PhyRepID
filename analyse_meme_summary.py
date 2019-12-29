## IvB September 2018
# making plots of pfam meme.json

import sys, os, glob, json
import numpy as np
import pipeline_methods as pre
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

species_mapping = pre.get_species_mapping(pre.species_mapping_file)
denovo_repeats_path = pre.root+'denovo/repeats/'
denovo_tblout_path = pre.root+'denovo/tblout/'
analysis_path = pre.root+'analysis/meme/'

meme_summary = {}

fig, axes = plt.subplots()

## Generate MEME summary
'''
meme_file_list = glob.glob(pre.denovo_treefix_path + '*.treefix.mpr.recon' ) 
print(len(meme_file_list))
		
for meme_repeats in meme_file_list:	
	og_id = meme_repeats[len(pre.denovo_treefix_path):-len('.treefix.mpr.recon')]
	meme_tblout = pre.denovo_tblout_path+og_id+'.tblout'
	
	repeats = pre.parse_domtblout_stats(meme_tblout)
	meme_summary[og_id] = repeats
	
	meme_summary[og_id]['consensus'] = ''
print(len(meme_summary))
with open(pre.root+'meme_summary.json', 'w') as meme_output:
	meme_output.write(json.dumps(meme_summary))

'''	

## Analyse meme summary
def calculate_repeat_cov(units_dict):
	repeat_aa = 0;begin=0;end=0;end_prev=0;max_dist=0;overlap=False
	units_list = [int(x) for x in units_dict.keys()]
	units_list.sort()
	
	for units_id in units_list:
		units_val=units_dict[str(units_id)]
		b,e = units_val['seq_coordinates']
		b,e = int(b), int(e)
		unit_length = e-b
		repeat_aa += unit_length
		
		if units_id == 1:
			begin = b
		elif units_id == len(units_list):
			end = e
		
		if units_id > 1:
			if (end_prev-b) >5: overlap=True
			max_dist = max((b-end_prev),max_dist)
		
		end_prev = e	
		
			
	repeat_cov = float(repeat_aa)/(float(end)-float(begin)) if float(end)>float(begin) else 0
	return({'repeat_cov':repeat_cov, 'repeat_aa':repeat_aa, 'max_dist':max_dist, 'overlap':overlap})

def create_dataframe_meme_summary(meme_summary):
	detailed_df_columns = ['og_id', 'consensus', 'model_length', 'ortholog_id', 'species', 'unit_cnt']
	detailed_df_rows = []
	exclude_ogs_list = []
	for og_id,repeats in meme_summary.items():
		consensus = repeats['consensus']
		dom_val = repeats['motif']
	
		if len(dom_val['orthologs_dict'].items())<4: continue
		for orth,orth_val in dom_val['orthologs_dict'].items():
			#if orth_val['unit_count'] < 3: continue
			species = species_mapping[filter(str.isalpha, str(orth))]
			#unit cnt is number of elements, also seq_bitscore and length (protein length)
			
			'''	
			#calculate how much of the protein is repeat
			o = calculate_repeat_cov(orth_val['units_dict'])
			repeat_cov = float(o['repeat_cov'])
			repeat_aa = float(o['repeat_aa'])
			max_dist = float(o['max_dist'])
			protein_cov = repeat_aa/float(orth_val['length'])
			if o['overlap'] == True: 
				exclude_ogs_list.append(og_id_domain)
								
			if repeat_cov >1 or protein_cov > 1: 
				exclude_ogs_list.append(og_id_domain)
			'''	
			detailed_df_rows.append([og_id, consensus, dom_val['length'], orth, species, orth_val['unit_count']])	
	
	#print('excluded ogs cnt:',len(np.unique(exclude_ogs_list)))
	#detailed_df_rows = [row for row in detailed_df_rows if row[0] not in exclude_ogs_list]
	meme_detailed_df = pd.DataFrame(detailed_df_rows, columns=detailed_df_columns)
	return(meme_detailed_df)


#Open files to make dataframes
#'''
with open(pre.root+'meme_summary.json', 'r') as meme_output:
	meme_summary = json.load(meme_output)
	
#make MEME df

meme_detailed_df = create_dataframe_meme_summary(meme_summary)

with open('meme_detailed_df.json','w') as output:
	output.write(json.dumps(meme_detailed_df.to_dict()))
#'''
#Read meme df
#'''
with open('meme_detailed_df.json','r') as output:
	meme_detailed_df = pd.DataFrame.from_dict(json.load(output))
final_og_cnt = len(meme_detailed_df.og_id.unique())
print(final_og_cnt)
#'''


#Distribution of MEME motif length
plt.title('MEME motif length ')
weights = np.ones_like(np.array(meme_detailed_df['model_length']))/float(len(np.array(meme_detailed_df['model_length'])))
plt.hist(meme_detailed_df['model_length'],weights=weights,color=pre.blue_colour,label='meme',bins=30)#,alpha=0.6)
plt.savefig(analysis_path+'plots/motif_length.png',format='png',bbox_inches='tight')
plt.show()
