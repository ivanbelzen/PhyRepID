# IvB 

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np

analysis_path = pre.root+'analysis/evo-events/'

meme_summary_df_file = pre.root+'meme_summary_df_reroot.json'
pfam_summary_df_file = pre.root+'pfam_summary_df_reroot.json'
combined_summary_df_file = 'combined_summary_df_reroot.json'

genetree_df_file = pre.root+'genetree_df_reroot.json'

schaper_summary_df_file = pre.root+'schaper_summary_df_reroot.json'
selectome_df_file = pre.root+'selectome_df_reroot.json'
exac_df_file = pre.root+'analysis/exac_10102018/exac_df.json'



'''
human_gene_symbol_file = pre.root+'human_gene_symbol.tsv'

with open(human_gene_symbol_file, 'r') as gene_data:
	mapping_gene_symbol = {}
	for line in gene_data:
		cols = line.strip().split('\t')
		mapping_gene_symbol[cols[0]]=cols[1]


'''

#Make combined summary dataframe with Duplication score  combined_summary_df_reroot.json
'''
combined_summary_df = pd.concat([pfam_summary_df,meme_summary_df], ignore_index=True, sort=False)

mean_combined_dup = combined_summary_df['netto_dup'].mean()
combined_summary_df['duplication_score']= (combined_summary_df['netto_dup']-mean_combined_dup)/combined_summary_df['orthologs_cnt']
std=combined_summary_df.describe()['duplication_score']['std']

print('#OGs in total',len(combined_summary_df.index))
print('#OGs >> netto dup - mean / orth cnt ', \
len(combined_summary_df.loc[combined_summary_df['duplication_score']>0].index), \
'>1stdev above mean',\
len(combined_summary_df.loc[combined_summary_df['duplication_score']>std].index))

#add gene symbol
combined_summary_df['gene_id'] = combined_summary_df['og_id'].str.split('_').str[1]
combined_summary_df['gene_symbol'] = combined_summary_df['gene_id'].map(mapping_gene_symbol)

#add clan for meme
combined_summary_df['clan'] = np.where(combined_summary_df['clan'].isna(), 'motif', combined_summary_df['clan'])
#add identifier with has both meme and pfam annotation
combined_summary_df['identifier'] = np.where(combined_summary_df['og_id_domain'].isna(),combined_summary_df['og_id'],combined_summary_df['og_id_domain'])

with open(combined_summary_df_file, 'w') as output:
	output.write(json.dumps(combined_summary_df.to_dict()))	

combined_summary_df[['gene_symbol','duplication_score','identifier','netto_dup','loss','orthologs_cnt','clan']]\
.to_csv(pre.root+'evo_events_summary_simple_reroot.csv')#, separator='\t')	
'''


'''
#For gene ontology, gene list 

analysis_path_datachar = pre.root+'analysis/dataset_04102018/'

pre_ensembl_gene_file = analysis_path_datachar+'pre_ensembl_gene_list.txt'
positive_gene_file = analysis_path_datachar+'pre_positive_gene_list.txt'

with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)	
	

full_dataset_gene_list = combined_summary_df['gene_id'].tolist()
positive_gene_list = combined_summary_df.loc[combined_summary_df['duplication_score']>0]['gene_id'].tolist()

print(len(combined_summary_df.loc[combined_summary_df['netto_dup']==0].index))
print(len(combined_summary_df.loc[combined_summary_df['duplication_score']>0].index))
quit()
with open(pre_ensembl_gene_file,'w') as output:
	for identifier in full_dataset_gene_list:
		output.write(identifier+"\n")

with open(positive_gene_file,'w') as output:
	for identifier in positive_gene_list:
		output.write(identifier+"\n")
		
'''

#Make combined summary df with updated Duplication Score, unit CV and repeat unit duplication score
with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)

combined_summary_df['genetree_id'] = combined_summary_df['identifier'].str.split('_').str[0]
combined_summary_df['gene_id'] = combined_summary_df['identifier'].str.split('_').str[1]	
combined_summary_df['og_id'] = combined_summary_df['genetree_id'] +'_'+ combined_summary_df['gene_id']

#add genetree
with open(genetree_df_file,'r') as output:
	gt = json.load(output)
	genetree_df = pd.DataFrame.from_dict(gt)

combined_summary_df = combined_summary_df.merge(genetree_df, on='og_id')

analysis_evo_events_path = pre.root+'analysis/evo_events_06102018/'

combined_df_file = 'pre_dataset_df.json'

mean_netto_dup = combined_summary_df['netto_dup'].mean()
combined_summary_df['duplication_score']= (combined_summary_df['netto_dup']-mean_netto_dup)/combined_summary_df['orthologs_cnt']

#alternative scores
combined_summary_df['scaled_dup']= combined_summary_df['netto_dup']/combined_summary_df['orthologs_cnt']
mean_scaled_dup = combined_summary_df['scaled_dup'].mean()
combined_summary_df['normalized_score']= combined_summary_df['scaled_dup']-mean_scaled_dup
combined_summary_df['expectancy_score']= combined_summary_df['netto_dup']-(mean_scaled_dup*combined_summary_df['orthologs_cnt'])

#unit count 
with open(pre.root+'analysis/dataset_04102018/'+'unit_cnt_df.json', 'r') as output:
	unit_cnt = json.load(output)
	unit_cnt_df = pd.DataFrame.from_dict(unit_cnt)

combined_summary_df = combined_summary_df.merge(unit_cnt_df, on='identifier', how='left')
combined_summary_df['unit_cv'] = combined_summary_df['cv']
combined_summary_df['unit_sum'] = combined_summary_df['sum']

combined_summary_df['unit_scaled_dup']  = combined_summary_df['netto_dup'] / combined_summary_df['unit_sum']
unit_dup_mean = combined_summary_df['unit_scaled_dup'].mean()
combined_summary_df['unit_duplication_score']= combined_summary_df['unit_scaled_dup']-unit_dup_mean

combined_summary_df[['genetree_id','gene_symbol','duplication_score','normalized_score','expectancy_score','unit_duplication_score','identifier','netto_dup','loss','orthologs_cnt','clan','unit_cv','unit_sum','gt_dup','gt_root']]\
.to_csv(analysis_evo_events_path+'evo_events.csv')

with open(analysis_evo_events_path+combined_df_file, 'w') as output:
	output.write(json.dumps(combined_summary_df.to_dict()))	

##Combine summary df with other datasets
with open(analysis_evo_events_path+combined_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)

with open(schaper_summary_df_file, 'r') as output:
	schaper_summ = json.load(output)
	schaper_summary_df = pd.DataFrame.from_dict(schaper_summ)
with open(selectome_df_file, 'r') as output:
	selectome = json.load(output)
	selectome_df = pd.DataFrame.from_dict(selectome)	
with open(exac_df_file, 'r') as output:
	exac = json.load(output)
	exac_df = pd.DataFrame.from_dict(exac)
			
#add schaper results
combined_summary_df = combined_summary_df.merge(schaper_summary_df[['identifier','schaper_pos']], on='identifier', how='left')

# dnds selectome
#selectome_df['identifier'] = np.where(selectome_df['og_id_domain'].empty,selectome_df['og_id'],selectome_df['og_id_domain'])
#add dnds selectome

selectome_positive_set = selectome_df.loc[selectome_df['positive_selection'] == 1]
combined_summary_df['selectome'] = np.where(combined_summary_df['identifier'].isin(selectome_positive_set['identifier']),1,0)

#exac
combined_summary_df = combined_summary_df.merge(exac_df[['gene_id','mis_z','pLI']], on='gene_id', how='left')

#export
combined_summary_df[['genetree_id','gene_symbol','duplication_score','identifier','netto_dup','loss','orthologs_cnt','clan','unit_cv','selectome','schaper_pos','mis_z','pLI','gt_dup','gt_root']]\
.to_csv(analysis_evo_events_path+'evo_events_merged.csv')


## Genetree analysis

analysis_evo_events_path = pre.root+'analysis/evo_events_18102018/'

sys.stdout = open(analysis_evo_events_path+'output.txt', 'w')

print('#OGs in total',len(combined_summary_df.index))
print('#genetrees in total',combined_summary_df['genetree_id'].nunique())

genetree_df = combined_summary_df.groupby(by=['genetree_id']).max()['netto_dup'].reset_index().sort_values('netto_dup', ascending=False)
mean_dup = genetree_df['netto_dup'].mean()

genetree_df = genetree_df.merge(combined_summary_df, on='genetree_id',how='inner')

genetree_df = genetree_df.drop_duplicates(subset='genetree_id')
genetree_df['dev_from_mean']= (genetree_df['netto_dup_x']-mean_dup)/genetree_df['orthologs_cnt']


std=genetree_df.describe()['dev_from_mean']['std']

print('#OGs in total',len(genetree_df.index))
print("mean netto dup",mean_dup)
print('#OGs >> netto dup - mean / orth cnt ', \
len(genetree_df.loc[genetree_df['dev_from_mean']>0].index), \
'>1stdev above mean',\
len(genetree_df.loc[genetree_df['dev_from_mean']>std].index))

genetree_df[['genetree_id','gene_symbol','duplication_score','identifier','netto_dup_x','loss','orthologs_cnt','clan','unit_cv','selectome','schaper_pos','mis_z','pLI','gt_dup','gt_root']]\
.to_csv(analysis_evo_events_path+'evo_events_no_genetree_duplicates.csv')


with open(analysis_evo_events_path+combined_df_file, 'w') as output:
	output.write(json.dumps(genetree_df.to_dict()))	
