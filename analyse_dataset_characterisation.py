## IvB August 2018
#Needs refactoring still!
# Note old meme_summary.json  pfam summary.json etc
 
import sys, os, glob, json
import numpy as np
import pipeline_methods as pre
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

species_mapping = pre.get_species_mapping(pre.species_mapping_file)
pfam_tblout_path = pre.root+'pfam/tblout/'
pfam_initial_path = pre.root+'pfam/hmm/'
analysis_path = pre.root+'analysis/pfam/'
output_file_prefix = 'pfam_data_' #for images and output files

pfam_summary = {}

species_mapping = pre.get_species_mapping(pre.species_mapping_file)

fig, axes = plt.subplots()

#See analyse_pfam_summary.py for
#	# Generate PFAM summary 
#	# Generate PFAM summary before iteration
#See analyse meme summary for generate MEME summary
# But make the dataframes here

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

def create_dataframe_pfam_summary(pfam_summary,orth_filter=True):
	detailed_df_columns = ['og_id_domain', 'model', 'clan', 'model_length', 'ortholog_id','species', 'protein_length', 'unit_cnt', 'repeat_aa', 'repeat_cov', 'protein_cov','max_dist']
	detailed_df_rows = []
	exclude_ogs_list = []
	for og_id_domain,repeats in pfam_summary.items():
		for domain,dom_val in repeats.items():
			if orth_filter == True:
				if len(dom_val['orthologs_dict'].items())<4: continue
			for orth,orth_val in dom_val['orthologs_dict'].items():
				#if orth_val['unit_count'] < 3: continue
				species = species_mapping[filter(str.isalpha, str(orth))]
				#unit cnt is number of elements, also seq_bitscore and length (protein length)
				
				#calculate how much of the protein is repeat
				o = calculate_repeat_cov(orth_val['units_dict'])
				repeat_cov = float(o['repeat_cov'])
				repeat_aa = float(o['repeat_aa'])
				max_dist = float(o['max_dist'])
				protein_cov = repeat_aa/float(orth_val['length'])
				if o['overlap'] == True: 
					exclude_ogs_list.append(og_id_domain)
				'''
				#if og_id_domain == 'ENSGT00890000139348_ENSG00000162631_Laminin_EGF': 
					print(og_id_domain, max_dist)
					units_list = [int(x) for x in orth_val['units_dict'].keys()]
					units_list.sort()
					print(units_list)
					for units_id in units_list:
						print(orth_val['units_dict'][str(units_id)]['seq_coordinates'])
					quit()
				#else: continue
				#if max_dist <20: print(og_id_domain, domain, max_dist, orth_val['unit_count'], repeat_cov)
				'''
				
				
				if repeat_cov >1 or protein_cov > 1: 
					exclude_ogs_list.append(og_id_domain)
				
				detailed_df_rows.append([og_id_domain, domain, dom_val['clan'], dom_val['length'], orth, species, orth_val['length'], orth_val['unit_count'], repeat_aa, repeat_cov, protein_cov, max_dist])	
	
	#print('excluded ogs cnt:',len(np.unique(exclude_ogs_list)))
	#detailed_df_rows = [row for row in detailed_df_rows if row[0] not in exclude_ogs_list]
	pfam_detailed_df = pd.DataFrame(detailed_df_rows, columns=detailed_df_columns)
	return(pfam_detailed_df)

def create_dataframe_meme_summary(meme_summary,orth_filter=True):
	detailed_df_columns = ['og_id', 'consensus', 'model_length', 'ortholog_id', 'protein_length', 'species', 'unit_cnt', 'repeat_aa', 'repeat_cov', 'protein_cov','max_dist']
	detailed_df_rows = []
	exclude_ogs_list = []
	for og_id,repeats in meme_summary.items():
		consensus = repeats['consensus']
		dom_val = repeats['motif']
	
		if orth_filter == True:
			if len(dom_val['orthologs_dict'].items())<4: continue
		for orth,orth_val in dom_val['orthologs_dict'].items():
			#if orth_val['unit_count'] < 3: continue
			species = species_mapping[filter(str.isalpha, str(orth))]
			#unit cnt is number of elements, also seq_bitscore and length (protein length)
			
			#calculate how much of the protein is repeat
			o = calculate_repeat_cov(orth_val['units_dict'])
			repeat_cov = float(o['repeat_cov'])
			repeat_aa = float(o['repeat_aa'])
			max_dist = float(o['max_dist'])
			protein_cov = repeat_aa/float(orth_val['length'])
			if o['overlap'] == True: 
				exclude_ogs_list.append(og_id)
								
			if repeat_cov >1 or protein_cov > 1: 
				exclude_ogs_list.append(og_id)
			
			detailed_df_rows.append([og_id, consensus, dom_val['length'], orth, orth_val['length'],species, orth_val['unit_count'],repeat_aa, repeat_cov, protein_cov, max_dist])	
	
	#print('excluded ogs cnt:',len(np.unique(exclude_ogs_list)))
	#detailed_df_rows = [row for row in detailed_df_rows if row[0] not in exclude_ogs_list]
	meme_detailed_df = pd.DataFrame(detailed_df_rows, columns=detailed_df_columns)
	return(meme_detailed_df)

#### DATAFRAMES
'''
## Make dataframes of PFAM
with open(pre.root+'pfam_summary.json', 'r') as pfam_output:
	pfam_summary = json.load(pfam_output)

with open(pre.root+'pfam_repeat_stats_initial.json', 'r') as pfam_output:
	pfam_summary_initial = json.load(pfam_output)

#make general Pfam dataframe
pfam_detailed_df = create_dataframe_pfam_summary(pfam_summary)
with open('pfam_detailed_df.json','w') as output:
	output.write(json.dumps(pfam_detailed_df.to_dict()))

#Make Pfam initial dataframes
#match up OG ids from initial to the OG-domain ids from later to compare improvement in repeat detection
filtered_pfam_summary_initial = {}
for og_id_domain in pfam_summary.keys():
	og_id = '_'.join(og_id_domain.split('_')[0:2])
	domain = og_id_domain[len(og_id)+1:]
	#print(og_id, domain)
	if og_id_domain not in filtered_pfam_summary_initial: filtered_pfam_summary_initial[og_id_domain] = {}
	filtered_pfam_summary_initial[og_id_domain][domain] = pfam_summary_initial[og_id][domain]

with open('pfam_repeat_stats_initial_filtered.json','w') as output:
	output.write(json.dumps(filtered_pfam_summary_initial))

initial_detailed_df = create_dataframe_pfam_summary(filtered_pfam_summary_initial)	
with open('initial_detailed_df.json','w') as output:
	output.write(json.dumps(initial_detailed_df.to_dict()))

#Make dataframe for unfiltered
unfiltered_initial_detailed_df = create_dataframe_pfam_summary(pfam_summary_initial)	
with open('unfiltered_initial_detailed_df.json','w') as output:
	output.write(json.dumps(unfiltered_initial_detailed_df.to_dict()))

## Make MEME dataframe
with open(pre.root+'meme_summary.json', 'r') as meme_output:
	meme_summary = json.load(meme_output)

meme_detailed_df = create_dataframe_meme_summary(meme_summary)
with open('meme_detailed_df.json','w') as output:
	output.write(json.dumps(meme_detailed_df.to_dict()))
'''

## Fraction of human repeat proteins in full dataset, pfam and meme
analysis_path_datachar = pre.root+'analysis/dataset_24102018/'
sys.stdout = open(analysis_path_datachar+'regular_dataset.txt', 'w')

with open(pre.root+'pfam_detailed_df.json','r') as output:
	pfam_df = pd.DataFrame.from_dict(json.load(output))
pfam_og_cnt = pfam_df.og_id_domain.nunique()

with open(pre.root+'meme_detailed_df.json','r') as output:
	meme_df = pd.DataFrame.from_dict(json.load(output))
meme_og_cnt = meme_df.og_id.nunique()

with open(pre.root+'unfiltered_initial_detailed_df.json','r') as output:
	initial_df = pd.DataFrame.from_dict(json.load(output))
initial_og_cnt = initial_df.og_id_domain.nunique()

print("#OGs:")
print("initial ", initial_og_cnt,"Pfam ",pfam_og_cnt,"MEME ", meme_og_cnt)

print("#human proteins:")
print("initial ", len(initial_df.loc[initial_df['species']=='HSAP'].index),\
"pfam ", len(pfam_df.loc[pfam_df['species']=='HSAP'].index),\
"meme ", len(meme_df.loc[meme_df['species']=='HSAP'].index),\
"total ", len(pfam_df.loc[pfam_df['species']=='HSAP'].index)+len(meme_df.loc[meme_df['species']=='HSAP'].index))

print("#human proteins >=3 unit cnt:")
print("initial ", len(initial_df.loc[(initial_df['species']=='HSAP') & (initial_df['unit_cnt']>=3)].index),\
"pfam ", len(pfam_df.loc[(pfam_df['species']=='HSAP') & (pfam_df['unit_cnt']>=3)].index),\
"meme ", len(meme_df.loc[(meme_df['species']=='HSAP') & (meme_df['unit_cnt']>=3)].index),\
"total ", len(pfam_df.loc[(pfam_df['species']=='HSAP') & (pfam_df['unit_cnt']>=3)].index)+len(meme_df.loc[(meme_df['species']=='HSAP') & (meme_df['unit_cnt']>=3)].index))

print("Count the actual number of unique protein ids instead of index")
#Now count the actual number of unique protein ids
print("#human proteins:")
print("initial ", initial_df.loc[initial_df['species']=='HSAP'].ortholog_id.nunique(),\
"pfam ", pfam_df.loc[pfam_df['species']=='HSAP'].ortholog_id.nunique(),\
"meme ", meme_df.loc[meme_df['species']=='HSAP'].ortholog_id.nunique(),\
"total ", pfam_df.loc[pfam_df['species']=='HSAP'].ortholog_id.nunique()+meme_df.loc[meme_df['species']=='HSAP'].ortholog_id.nunique())

print("#human proteins >=3 unit cnt:")
print("initial ", initial_df.loc[(initial_df['species']=='HSAP') & (initial_df['unit_cnt']>=3)].ortholog_id.nunique(),\
"pfam ", pfam_df.loc[(pfam_df['species']=='HSAP') & (pfam_df['unit_cnt']>=3)].ortholog_id.nunique(),\
"meme ", meme_df.loc[(meme_df['species']=='HSAP') & (meme_df['unit_cnt']>=3)].ortholog_id.nunique(),\
"total ", pfam_df.loc[(pfam_df['species']=='HSAP') & (pfam_df['unit_cnt']>=3)].ortholog_id.nunique()+meme_df.loc[(meme_df['species']=='HSAP') & (meme_df['unit_cnt']>=3)].ortholog_id.nunique())


#OK... and now without orthology constraints?
## Conclusion: I don't have that data on final pfam/meme where it has iterated but not orth constraints ofcourse
'''
#Make no filter dataframes

with open(pre.root+'meme_summary.json', 'r') as meme_output:
	meme_summary = json.load(meme_output)
meme_df = create_dataframe_meme_summary(meme_summary,orth_filter=False)
with open(analysis_path_datachar+'meme_df_nofilter.json','w') as output:
	output.write(json.dumps(meme_df.to_dict()))

with open(pre.root+'pfam_summary.json', 'r') as pfam_output:
	pfam_summary = json.load(pfam_output)
pfam_df = create_dataframe_pfam_summary(pfam_summary,orth_filter=False)

with open(analysis_path_datachar+'pfam_df_nofilter.json','w') as output:
	output.write(json.dumps(pfam_df.to_dict()))

with open(pre.root+'pfam_repeat_stats_initial.json', 'r') as pfam_output:
	pfam_summary_initial = json.load(pfam_output)
initial_df = create_dataframe_pfam_summary(pfam_summary_initial,orth_filter=False)
with open(analysis_path_datachar+'initial_df_nofilter.json','w') as output:
	output.write(json.dumps(initial_df.to_dict()))	

#Actual analysis
sys.stdout = open(analysis_path_datachar+'nofilter_dataset.txt', 'w')

with open(analysis_path_datachar+'pfam_df_nofilter.json','r') as output:
	pfam_df = pd.DataFrame.from_dict(json.load(output))
pfam_og_cnt = pfam_df.og_id_domain.nunique()

with open(analysis_path_datachar+'meme_df_nofilter.json','r') as output:
	meme_df = pd.DataFrame.from_dict(json.load(output))
meme_og_cnt = meme_df.og_id.nunique()

with open(analysis_path_datachar+'initial_df_nofilter.json','r') as output:
	initial_df = pd.DataFrame.from_dict(json.load(output))
initial_og_cnt = initial_df.og_id_domain.nunique()

print("#OGs:")
print("initial ", initial_og_cnt,"Pfam ",pfam_og_cnt,"MEME ", meme_og_cnt)

print("#human proteins:")
print("initial ", len(initial_df.loc[initial_df['species']=='HSAP'].index),\
"pfam ", len(pfam_df.loc[pfam_df['species']=='HSAP'].index),\
"meme ", len(meme_df.loc[meme_df['species']=='HSAP'].index),\
"total ", len(pfam_df.loc[pfam_df['species']=='HSAP'].index)+len(meme_df.loc[meme_df['species']=='HSAP'].index))

print("#human proteins >=3 unit cnt:")
print("initial ", len(initial_df.loc[(initial_df['species']=='HSAP') & (initial_df['unit_cnt']>=3)].index),\
"pfam ", len(pfam_df.loc[(pfam_df['species']=='HSAP') & (pfam_df['unit_cnt']>=3)].index),\
"meme ", len(meme_df.loc[(meme_df['species']=='HSAP') & (meme_df['unit_cnt']>=3)].index),\
"total ", len(pfam_df.loc[(pfam_df['species']=='HSAP') & (pfam_df['unit_cnt']>=3)].index)+len(meme_df.loc[(meme_df['species']=='HSAP') & (meme_df['unit_cnt']>=3)].index))

'''
	
## Read dataframes used in multiple analyses
'''
#General pfam dataframe
with open(pre.root+'pfam_detailed_df.json','r') as output:
	pfam_detailed_df = pd.DataFrame.from_dict(json.load(output))
final_og_cnt = pfam_detailed_df.og_id_domain.nunique()

#Get top 10 clans
og_cnt_clan = pfam_detailed_df.groupby(by='clan', as_index=False).agg({'og_id_domain': pd.Series.nunique}).sort_values('og_id_domain', ascending=False)
og_cnt_clan['fraction'] = og_cnt_clan['og_id_domain'].divide(final_og_cnt)
og_cnt_clan['group'] = 'final'
top_clans = og_cnt_clan.head(n=10)

#Read meme df
with open('meme_detailed_df.json','r') as output:
	meme_detailed_df = pd.DataFrame.from_dict(json.load(output))
meme_og_cnt = meme_detailed_df.og_id.nunique()

#Read initial df,filtered on OGs in final dataset
with open(pre.root+'initial_detailed_df.json','r') as output:
	initial_detailed_df = pd.DataFrame.from_dict(json.load(output))
initial_filtered_og_cnt = pfam_detailed_df.og_id_domain.nunique()

#Evo events dataframe
with open(pre.root+'combined_summary_df_reroot.json', 'r') as output:
	evo_events = json.load(output)
	evo_events_df = pd.DataFrame.from_dict(evo_events)

pfam_evo_events_df = evo_events_df.loc[evo_events_df['clan']!='motif']
pfam_evo_events_cnt = len(pfam_evo_events_df.index)

meme_evo_events_df = evo_events_df.loc[evo_events_df['clan']=='motif']
meme_evo_events_cnt = len(meme_evo_events_df.index)

overlap_evo_events_df = pfam_evo_events_df.loc[pfam_evo_events_df['og_id'].isin(meme_evo_events_df['og_id'])]
overlap_evo_events_cnt = len(overlap_evo_events_df.index)
'''
## Dataset characterisation analysis 

'''
analysis_path_datachar = pre.root+'analysis/dataset_04102018/'
sys.stdout = open(analysis_path_datachar+'output.txt', 'w')
print(len(evo_events_df.index))
print(overlap_evo_events_df.to_string())


#Write dataframes to csv files
meme_detailed_df.to_csv(analysis_path_datachar+'meme_detailed_df.csv')
pfam_detailed_df.to_csv(analysis_path_datachar+'pfam_detailed_df.csv')
'''

#PFAM dataset
#Stats
#print("PFAM dataset")

#Clans, plot top 10 initial vs final pfam dataset
'''
#Read initial unfiltered df
with open(pre.root+'unfiltered_initial_detailed_df.json','r') as output:
	unfiltered_initial_detailed_df = pd.DataFrame.from_dict(json.load(output))
initial_og_cnt = unfiltered_initial_detailed_df.og_id_domain.nunique()

print("Final dataset #OGs",final_og_cnt,"Initial set filtered on final ogs",initial_filtered_og_cnt,"Initial set full", initial_og_cnt,"Pfam evo events",pfam_evo_events_cnt)

og_cnt_clan_initial = unfiltered_initial_detailed_df.groupby(by='clan', as_index=False).agg({'og_id_domain': pd.Series.nunique}).sort_values('og_id_domain', ascending=False)
og_cnt_clan_initial['fraction'] = og_cnt_clan_initial['og_id_domain'].divide(initial_og_cnt)
og_cnt_clan_initial['group'] = 'initial'
#print("Top 10 clans initial dataset:", og_cnt_clan_initial.head(n=10)) 

#Plot of #OGs per clan, show the top 10

view_clan_df = pd.concat([top_clans,og_cnt_clan_initial.loc[og_cnt_clan_initial['clan'].isin(top_clans['clan'])]])
sns.barplot(data=view_clan_df, y='clan', x='fraction', hue='group', palette=[pre.green_colour,pre.blue_colour], hue_order=['final','initial'])
axes.set_xlabel('%OGs')
plt.title('Top 10 clans in dataset')
plt.savefig(analysis_path_datachar+'top10_clans_initial_vs_final.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()

enrichment_df = og_cnt_clan.merge(og_cnt_clan_initial, on='clan')
enrichment_df['enrichment'] = enrichment_df['fraction_x'] / enrichment_df['fraction_y']
print(enrichment_df.loc[enrichment_df['clan'].isin(top_clans['clan'])].to_string())

#for report
print(enrichment_df.loc[enrichment_df['clan'].isin(top_clans['clan'])][['clan','enrichment']].to_string())
'''

#Compare to positives in analyse_evo_events_summary.py
#'''


#MEME dataset
# Stats
'''
print("MEME #OGs",meme_og_cnt,"MEME evo events",meme_evo_events_cnt,"overlap evo event",overlap_evo_events_cnt)

#Distribution of MEME motif length
plt.title('MEME motif length ')
weights = np.ones_like(np.array(meme_detailed_df['model_length']))/float(len(np.array(meme_detailed_df['model_length'])))
plt.hist(meme_detailed_df['model_length'],weights=weights,color=pre.blue_colour,label='meme',bins=30)
sns.rugplot(meme_detailed_df['model_length'], color=pre.green_colour,height=0.02)
plt.savefig(analysis_path_datachar+'motif_length.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
'''

# Visualize differences between top 10 clans and MEME
'''
top_clans_list = top_clans['clan'].tolist()
meme_detailed_df['clan'] = 'motif'
meme_detailed_df['identifier'] = meme_detailed_df['og_id']
pfam_detailed_df['identifier'] = pfam_detailed_df['og_id_domain']
top_clans_meme = pd.concat([pfam_detailed_df.loc[pfam_detailed_df['clan'].isin(top_clans_list)],meme_detailed_df], sort=False)
top_clans_list = top_clans_list.append(['meme'])
#top_clans_meme = pfam_detailed_df.loc[pfam_detailed_df['clan'].isin(top_clans_list)]

#Unit count distribution for clans in top 10
sns.violinplot(data=top_clans_meme, x='unit_cnt', y='clan', cut=0, scale='count', order=top_clans_list)
axes.set_xlabel('#units in repeat')
plt.title('Unit count per clan')
plt.xscale('log')
plt.savefig(analysis_path_datachar+'unit_cnt_per_clan.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Model length distribution for clans in top 10
sns.violinplot(data=top_clans_meme, x='model_length', y='clan', cut=0, scale='count', order=top_clans_list)
axes.set_xlabel('model length in aa')
plt.title('Model length per clan')
plt.savefig(analysis_path_datachar+'model_length_per_clan.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Protein length distribution for clans in top 10
sns.violinplot(data=top_clans_meme, x='protein_length', y='clan', cut=0, scale='count', order=top_clans_list)
axes.set_xlabel('protein length in aa')
plt.title('Protein length per clan')
plt.xscale('log')
plt.savefig(analysis_path_datachar+'protein_length_per_clan.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

## Data distributions boxplots top 10 only
'''
analysis_path_datachar = pre.root+'analysis/dataset_16102018/'
sys.stdout = open(analysis_path_datachar+'output.txt', 'w')

top_clans_list = top_clans['clan'].tolist()
top_clans = pfam_detailed_df.loc[pfam_detailed_df['clan'].isin(top_clans_list)]

print('#OGs',pfam_detailed_df.og_id_domain.nunique())
print('#OGs in top 10 clans',top_clans.og_id_domain.nunique())
print('Figures are about top 10 clans only')
print('#OGs with less than 3 repeat units in some protein',top_clans.loc[top_clans['unit_cnt']<3].og_id_domain.nunique())
print("Filter is not used for regular per clan figures")
'''
#Plots
'''
#Unit count distribution for clans in top 10
axes = sns.boxplot(data=top_clans, x='unit_cnt', y='clan',   order=top_clans_list)
axes.set_xlabel('#units in repeat')
plt.title('Unit count per clan')
plt.xscale('log')
axes.set_xticks([1, 10, 100])
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'unit_cnt_per_clan_boxplot_log.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Model length distribution for clans in top 10
axes = sns.boxplot(data=top_clans, x='model_length', y='clan',   order=top_clans_list)
axes.set_xlabel('model length in aa')
plt.title('Model length per clan')
plt.xscale('log')
axes.set_xticks([20, 30, 40, 50, 100, 200])
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'model_length_per_clan_boxplot_log.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Protein length distribution for clans in top 10
axes = sns.boxplot(data=top_clans, x='protein_length', y='clan',   order=top_clans_list)
axes.set_xlabel('protein length in aa')
plt.title('Protein length per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'protein_length_per_clan_boxplot_log.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''
## HUMAN only, supplementary figures 1
'''
human_top_clans = pfam_detailed_df.loc[ (pfam_detailed_df['clan'].isin(top_clans_list)) & (pfam_detailed_df['species'] == 'HSAP') ]
print('#OGs with a human protein less than 3 repeat units',human_top_clans.loc[human_top_clans['unit_cnt']<3].og_id_domain.nunique())

print('#OGs with at least 3 repeat units in 1 human protein',human_top_clans.loc[human_top_clans['unit_cnt']>=3].og_id_domain.nunique())

human_top_clans = pfam_detailed_df.loc[ (pfam_detailed_df['clan'].isin(top_clans_list)) & (pfam_detailed_df['species'] == 'HSAP') & (pfam_detailed_df['unit_cnt']>=3)]
print('#OGs in top 10 clans with which the human figures were made ',human_top_clans.og_id_domain.nunique())

#Unit count distribution for clans in top 10
axes = sns.boxplot(data=human_top_clans, x='unit_cnt', y='clan', order=top_clans_list)
axes.set_xlabel('#units in repeat')
plt.title('Unit count per clan')
plt.xscale('log')
axes.set_xticks([3, 7, 10, 30, 50])
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'unit_cnt_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Model length distribution for clans in top 10
axes = sns.boxplot(data=human_top_clans, x='model_length', y='clan',   order=top_clans_list)
axes.set_xlabel('model length in aa')
plt.title('Model length per clan')
plt.xscale('log')
axes.set_xticks([20, 30, 40, 50, 100, 200])
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'model_length_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Protein length distribution for clans in top 10
axes = sns.boxplot(data=human_top_clans, x='protein_length', y='clan',   order=top_clans_list)
axes.set_xlabel('protein length in aa')
plt.title('Protein length per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'protein_length_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Protein coverage for clans in top 10
axes = sns.boxplot(data=human_top_clans, x='protein_cov', y='clan',   order=top_clans_list)
axes.set_xlabel('protein coverage in aa')
plt.title('protein coverage per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'protein_coverage_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Repeat coverage for clans in top 10
axes = sns.boxplot(data=human_top_clans, x='repeat_cov', y='clan',   order=top_clans_list)
axes.set_xlabel('repeat coverage in aa')
plt.title('repeat coverage per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'repeat_coverage_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Max dist distribution for clans in top 10
axes = sns.boxplot(data=human_top_clans, x='max_dist', y='clan',   order=top_clans_list)
axes.set_xlabel('Max dist between repeats in aa')
plt.title('Max dist between repeats per clan')
plt.xscale('log')
axes.set_xticks([0, 1, 10, 100, 1000])
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'max_dist_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

# VIOLIN, useless
'''
#Protein coverage for clans in top 10
axes = sns.violinplot(data=human_top_clans, x='protein_cov', y='clan', cut=0, scale='width',  order=top_clans_list)
axes.set_xlabel('protein coverage in aa')
plt.title('protein coverage per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'protein_coverage_per_clan_violinplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Repeat coverage for clans in top 10
axes = sns.violinplot(data=human_top_clans, x='repeat_cov', y='clan', cut=0, scale='width', order=top_clans_list)
axes.set_xlabel('repeat coverage in aa')
plt.title('repeat coverage per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'repeat_coverage_per_clan_violinplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Max dist distribution for clans in top 10
axes = sns.violinplot(data=human_top_clans, x='max_dist', y='clan', cut=0, scale='width', order=top_clans_list)
axes.set_xlabel('Max dist between repeats in aa')
plt.title('Max dist between repeats per clan')
plt.xscale('log')
axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'max_dist_per_clan_violinplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''


## Comparison distribution figures meme/initial/pfam
'''
#Distributions for full dataset, initial, meme
#Exploratory figures, not used
distributions = ['unit_cnt','model_length','protein_length','protein_cov','repeat_cov']
#distributions = ['protein_cov','repeat_cov']
for dist in distributions:
	sns.distplot(initial_detailed_df[dist], color=pre.blue_colour, label='initial')
	sns.distplot(pfam_detailed_df[dist], color=pre.green_colour, label='pfam')
	sns.distplot(meme_detailed_df[dist], color=pre.purple_colour, label='meme')
	plt.title('Compare distribution: '+dist)
	plt.legend()
	plt.savefig(analysis_path_datachar+'compare_dist_'+dist+'.png',format='png',bbox_inches='tight')
	plt.clf()
	
	plt.savefig(analysis_path_datachar+'compare_dist_'+dist+'_cut.png',format='png',bbox_inches='tight')
'''

## Figure 1 Model length and repeat coverage
'''
##Model length meme vs pfam cut at 350
meme_detailed_df['model_length'] = np.where(meme_detailed_df['model_length']>350,350,meme_detailed_df['model_length'])
pfam_detailed_df['model_length'] = np.where(pfam_detailed_df['model_length']>350,350,pfam_detailed_df['model_length'])

weights = np.ones_like(np.array(pfam_detailed_df['model_length']))/float(len(np.array(pfam_detailed_df['model_length'])))
plt.hist(pfam_detailed_df['model_length'],weights=weights,color=pre.green_colour,label='pfam',bins=40,alpha=0.8)
#sns.rugplot(pfam_detailed_df['model_length'],color=pre.green_colour,height=0.01)

weights = np.ones_like(np.array(meme_detailed_df['model_length']))/float(len(np.array(meme_detailed_df['model_length'])))
plt.hist(meme_detailed_df['model_length'],weights=weights,color=pre.purple_colour,label='meme',bins=20,alpha=0.8)
#sns.rugplot(meme_detailed_df['model_length'],color=pre.purple_colour,height=0.01)

plt.xlabel('model length (aa)')
plt.title('Model length')
plt.legend()
plt.savefig(analysis_path_datachar+'compare_model_length_meme_pfam_cut350.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

##Repeat coverage meme vs pfam
weights = np.ones_like(np.array(pfam_detailed_df['repeat_cov']))/float(len(np.array(pfam_detailed_df['repeat_cov'])))
plt.hist(pfam_detailed_df['repeat_cov'],weights=weights,color=pre.green_colour,label='pfam',bins=40,alpha=0.8)
#sns.rugplot(pfam_detailed_df['repeat_cov'],color=pre.green_colour,height=0.01)

weights = np.ones_like(np.array(meme_detailed_df['repeat_cov']))/float(len(np.array(meme_detailed_df['repeat_cov'])))
plt.hist(meme_detailed_df['repeat_cov'],weights=weights,color=pre.purple_colour,label='meme',bins=35,alpha=0.8)
#sns.rugplot(meme_detailed_df['repeat_cov'],color=pre.purple_colour,height=0.01)

plt.xlabel('repeat coverage')
plt.title('Repeat coverage')
plt.legend()
plt.savefig(analysis_path_datachar+'compare_repeat_cov_meme_pfam.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

#USELESS
#Distribution of unit count stdev
'''
stdev_unit_cnt = {}
og_id_domain_list = pfam_detailed_df.og_id_domain.unique()
pfam_std_list, initial_std_list = [],[]
for og_id_domain in og_id_domain_list:
	if pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.mean() < 3: continue
	if initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.mean() < 3: continue
	std_pfam = pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.std()
	std_initial = initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.std()
	stdev_unit_cnt[og_id_domain]={'final':std_pfam,'initial':std_initial}
	pfam_std_list.append(std_pfam)
	initial_std_list.append(std_initial)

#'' '	
with open(analysis_path_datachar+'stdev_unit_cnt.json','w') as output:
	output.write(json.dumps(stdev_unit_cnt))
#'' '

weights = np.ones_like(np.array(pfam_std_list))/float(len(np.array(pfam_std_list)))
plt.hist(pfam_std_list,weights=weights,color=pre.green_colour,label='pfam',bins=20,alpha=0.8)
sns.rugplot(pfam_std_list,color=pre.green_colour)
weights = np.ones_like(np.array(initial_std_list))/float(len(np.array(initial_std_list)))
plt.hist(initial_std_list,weights=weights,color=pre.blue_colour,label='initial',bins=20,alpha=0.8)
sns.rugplot(initial_std_list,color=pre.blue_colour)
plt.legend()
plt.savefig(analysis_path_datachar+'stdev_rug.png',format='png',bbox_inches='tight')
#plt.show()
'''
#USELESS
#Max distance between repeats
'''
#pfam_detailed_df['max_dist'] = np.where(pfam_detailed_df['max_dist'] > 50 , 50,pfam_detailed_df['max_dist'])
#initial_detailed_df['max_dist'] = np.where(initial_detailed_df['max_dist'] > 50 , 50,initial_detailed_df['max_dist'])

print(pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'] == 'ENSGT00830000128299_ENSG00000162692_Ig_3']['max_dist'])
print(initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == 'ENSGT00830000128299_ENSG00000162692_Ig_3']['max_dist'])
print(pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'] == 'ENSGT00830000128299_ENSG00000162692_Ig_3']['max_dist'].mean())
print(pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'] == 'ENSGT00830000128299_ENSG00000162692_Ig_3']['max_dist'].std())
print(initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == 'ENSGT00830000128299_ENSG00000162692_Ig_3']['max_dist'].mean())
print(initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == 'ENSGT00830000128299_ENSG00000162692_Ig_3']['max_dist'].std())

initial_maxdist_stdev = initial_detailed_df.groupby(by='og_id_domain')['max_dist'].std()
pfam_maxdist_stdev = pfam_detailed_df.groupby(by='og_id_domain')['max_dist'].std()
initial_maxdist_stdev.to_csv(analysis_path_datachar+'maxdist_initial.csv')
pfam_maxdist_stdev.to_csv(analysis_path_datachar+'maxdist_pfam.csv')

weights = np.ones_like(np.array(initial_maxdist_stdev))/float(len(np.array(initial_maxdist_stdev)))
plt.hist(initial_maxdist_stdev,weights=weights,color=pre.blue_colour,label='initial',bins=20,alpha=0.6)
weights = np.ones_like(np.array(pfam_maxdist_stdev))/float(len(np.array(pfam_maxdist_stdev)))
plt.hist(pfam_maxdist_stdev, weights=weights,color=pre.green_colour,label='final',bins=20,alpha=0.6)
plt.legend()
#plt.yscale('symlog')
plt.xlabel('maxdist')
plt.savefig(analysis_path_datachar+'maxdist_std.png',format='png',bbox_inches='tight')
#plt.savefig(analysis_path+'plots/compare_repeat_maxdist_cut50.png',format='png',bbox_inches='tight')
'''

# THIS SHOWS AN IMPROVEMENT!!
#Repeat unit count coefficient of variation
'''
sys.stdout = open(analysis_path_datachar+'coefficient_variation.txt', 'w')

print('Repeat unit count coefficient of variation')

initial_unit_cnt = initial_detailed_df.groupby(by='og_id_domain')['unit_cnt'].agg([np.mean,np.std])
initial_unit_cnt.reset_index(inplace=True)
initial_unit_cnt['cv'] = initial_unit_cnt['std']/initial_unit_cnt['mean']

pfam_unit_cnt = pfam_detailed_df.groupby(by='og_id_domain')['unit_cnt'].agg([np.mean,np.std,sum])
pfam_unit_cnt.reset_index(inplace=True)
pfam_unit_cnt['cv'] = pfam_unit_cnt['std']/pfam_unit_cnt['mean']

merged_df = initial_unit_cnt.merge(pfam_unit_cnt, on='og_id_domain')
merged_df['diff'] = merged_df['cv_y']-merged_df['cv_x']


with open(analysis_path_datachar+'unit_cnt_difference_df.json', 'w') as output:
	output.write(json.dumps(merged_df.to_dict()))

print('#OGs decrease in CV',len(merged_df.loc[merged_df['diff']<0].index),float(len(merged_df.loc[merged_df['diff']<0].index))/len(merged_df.index))
print('#OGs become 0', len(merged_df.loc[(merged_df['cv_x']>0) & (merged_df['cv_y']==0)].index),\
float(len(merged_df.loc[(merged_df['cv_x']>0) & (merged_df['cv_y']==0)].index))/len(merged_df.index))
print('#OGs decrease or stay same in CV',len(merged_df.loc[merged_df['diff']<=0].index),float(len(merged_df.loc[merged_df['diff']<=0].index))/len(merged_df.index))
'''

#Export unit cnt dataframe
'''
pfam_unit_cnt = pfam_detailed_df.groupby(by='og_id_domain')['unit_cnt'].agg([np.mean,np.std,sum])
pfam_unit_cnt['cv'] = pfam_unit_cnt['std']/pfam_unit_cnt['mean']
pfam_unit_cnt.reset_index(inplace=True)

meme_unit_cnt = meme_detailed_df.groupby(by='og_id')['unit_cnt'].agg([np.mean,np.std,sum])
meme_unit_cnt['cv'] = meme_unit_cnt['std']/meme_unit_cnt['mean']
meme_unit_cnt.reset_index(inplace=True)

unit_cnt_df = pd.concat([pfam_unit_cnt, meme_unit_cnt], ignore_index=True, sort=False)
unit_cnt_df['identifier'] = np.where(unit_cnt_df['og_id_domain'].isna(),unit_cnt_df['og_id'],unit_cnt_df['og_id_domain'])

with open(analysis_path_datachar+'unit_cnt_df.json', 'w') as output:
	output.write(json.dumps(unit_cnt_df.to_dict()))
'''

#Plots
'''
weights = np.ones_like(np.array(merged_df['diff']))/float(len(np.array(merged_df['diff'])))
plt.hist(merged_df['diff'],weights=weights,color=pre.blue_colour,bins=50)
plt.title('Pfam repeat unit coefficient of variation')
plt.xlabel('CV final - CV initial')
plt.axvline(x=0,color='black')
plt.savefig(analysis_path_datachar+'unit_cnt_cv_diff.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

weights = np.ones_like(np.array(initial_unit_cnt['cv']))/float(len(np.array(initial_unit_cnt['cv'])))
plt.hist(initial_unit_cnt['cv'],weights=weights,color=pre.blue_colour,label='initial',bins=30,alpha=0.6)
sns.rugplot(initial_unit_cnt['cv'],color=pre.blue_colour,height=0.02)
weights = np.ones_like(np.array(pfam_unit_cnt['cv']))/float(len(np.array(pfam_unit_cnt['cv'])))
plt.hist(pfam_unit_cnt['cv'], weights=weights,color=pre.green_colour,label='final',bins=30,alpha=0.6)
sns.rugplot(pfam_unit_cnt['cv'],color=pre.green_colour,height=0.02)
plt.legend()
#plt.yscale('symlog')
plt.xlabel('repeat unit cnt coefficient of variation')
plt.savefig(analysis_path_datachar+'unit_cnt_cv.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

# Plot CV per clan
'''
human_top_clans = human_top_clans.merge(merged_df, on='og_id_domain', how='left')

axes = sns.boxplot(data=human_top_clans, x='diff', y='clan', order=top_clans_list)
axes.set_xlabel('CV final - CV initial')
plt.title('Unit CV difference per clan')
plt.xscale('symlog')
#axes.set_xticks([3, 7, 10, 100])
plt.axvline(x=0, color='grey')
#axes.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
plt.savefig(analysis_path_datachar+'unit_cv_diff_per_clan_boxplot_log_human.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

