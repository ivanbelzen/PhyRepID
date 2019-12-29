## IvB August 2018
# making plots of pfam summary.json

import sys, os, glob, json
import numpy as np
import pipeline_methods as pre
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

species_mapping = pre.get_species_mapping(pre.species_mapping_file)
pfam_tblout_path = pre.root+'pfam/tblout/'
pfam_initial_path = pre.root+'pfam/hmm/'
analysis_path = pre.root+'analysis/pfam/'
output_file_prefix = 'pfam_data_' #for images and output files

pfam_summary = {}

species_mapping = pre.get_species_mapping(pre.species_mapping_file)

fig, axes = plt.subplots()


## Generate PFAM summary
#'''
pfam_file_list = glob.glob(pre.pfam_treefix_path + '*.treefix.mpr.recon' ) 
print(len(pfam_file_list))

exclude_og_cnt = []	
for pfam_file in pfam_file_list:	
	og_hit_id = pfam_file[len(pre.pfam_treefix_path):-len('.treefix.mpr.recon')]
	pfam_tblout = pre.pfam_tblout_path+og_hit_id+'.tblout'
	if not pre.file_notempty(pfam_tblout):continue
	repeats = pre.parse_domtblout_stats(pfam_tblout)
	
	domain = '_'.join(og_hit_id.split('_')[2:])
	if len(repeats[domain]['orthologs_dict'].items()) < 4:
		exclude_og_cnt.append(og_hit_id)
		continue
			
	pfam_summary[og_hit_id] = repeats	

print(len(pfam_summary))

with open(pre.root+'pfam_exclude_ogs_lt4.txt','w') as output:
	for x in exclude_og_cnt:
		output.write(x+'\n')

with open(pre.root+'pfam_summary.json', 'w') as output:
	output.write(json.dumps(pfam_summary))
#'''

## Generate PFAM summary before iteration
'''
pfam_file_list = glob.glob(pfam_initial_path + '*.tblout' ) 
print(len(pfam_file_list))
	
for pfam_tblout in pfam_file_list:	
	og_hit_id = pfam_tblout[len(pfam_initial_path):-len('.tblout')]
	repeats = pre.parse_domtblout_stats(pfam_tblout, True)
	pfam_summary[og_hit_id] = repeats
	
with open(pre.root+'pfam_summary_initial_gacutoff.json', 'w') as output:
	output.write(json.dumps(pfam_summary))
'''


## Analyse pfam summary
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
		
		'''
		unit_length = e-b
		repeat_aa += unit_length
		begin = b if begin == 0 else min(b,begin)
		end = e if end == 0 else max(e,end) 
		if end_prev == 0: end_prev = b  
		max_dist = b-end_prev if max_dist == 0 else max(b-end_prev,max_dist)
		end_prev = e
		'''
		
	repeat_cov = float(repeat_aa)/(float(end)-float(begin)) if float(end)>float(begin) else 0
	return({'repeat_cov':repeat_cov, 'repeat_aa':repeat_aa, 'max_dist':max_dist, 'overlap':overlap})

def create_dataframe_pfam_summary(pfam_summary):
	detailed_df_columns = ['og_id_domain', 'model', 'clan', 'model_length', 'ortholog_id','species', 'protein_length', 'unit_cnt', 'repeat_aa', 'repeat_cov', 'protein_cov','max_dist']
	detailed_df_rows = []
	exclude_ogs_list = []
	for og_id_domain,repeats in pfam_summary.items():
		for domain,dom_val in repeats.items():
			if len(dom_val['orthologs_dict'].items())<4: continue
			for orth,orth_val in dom_val['orthologs_dict'].items():
				if orth_val['unit_count'] < 3: continue
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

'''
#Open files to make dataframes

with open(pre.root+'pfam_summary.json', 'r') as pfam_output:
	pfam_summary = json.load(pfam_output)

with open(pre.root+'pfam_summary_initial.json', 'r') as pfam_output:
	pfam_summary_initial = json.load(pfam_output)

#make Pfam df
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

with open('initial_pfam_filtered.json','w') as output:
	output.write(json.dumps(filtered_pfam_summary_initial))
	
#print(len(filtered_pfam_summary_initial),len(pfam_summary_initial))
initial_detailed_df = create_dataframe_pfam_summary(filtered_pfam_summary_initial)	

with open('initial_detailed_df.json','w') as output:
	output.write(json.dumps(initial_detailed_df.to_dict()))

#Make dataframe for unfiltered
unfiltered_initial_detailed_df = create_dataframe_pfam_summary(pfam_summary_initial)	
#print(len(unfiltered_initial_detailed_df.index))

with open('unfiltered_initial_detailed_df.json','w') as output:
	output.write(json.dumps(unfiltered_initial_detailed_df.to_dict()))

#'''

#Read pfam df
#'''
with open('pfam_detailed_df.json','r') as output:
	pfam_detailed_df = pd.DataFrame.from_dict(json.load(output))
final_og_cnt = len(pfam_detailed_df.og_id_domain.unique())
print(final_og_cnt)

#Read initial df
with open('initial_detailed_df.json','r') as output:
	initial_detailed_df = pd.DataFrame.from_dict(json.load(output))

#Read initial unfiltered df
with open('unfiltered_initial_detailed_df.json','r') as output:
	unfiltered_initial_detailed_df = pd.DataFrame.from_dict(json.load(output))
initial_og_cnt = len(unfiltered_initial_detailed_df.og_id_domain.unique())
print(initial_og_cnt)
#'''

#read pfam summary df from analyse evo events
with open('pfam_summary_df.json', 'r') as output:
	evo_events = json.load(output)
	evo_events_df = pd.DataFrame.from_dict(evo_events)
positive_set = evo_events_df.loc[evo_events_df['dev_from_mean']>0]
	
#Clans, plot top 10 initial vs final pfam dataset, and vs positive set
#'''
og_cnt_clan = pfam_detailed_df.groupby(by='clan', as_index=False).agg({'og_id_domain': pd.Series.nunique}).sort_values('og_id_domain', ascending=False)
og_cnt_clan['fraction'] = og_cnt_clan['og_id_domain'].divide(final_og_cnt)
og_cnt_clan['group'] = 'final'
top_clans = og_cnt_clan.head(n=10)
#print("Top 10 clans final dataset:",top_clans) 

og_cnt_clan_initial = unfiltered_initial_detailed_df.groupby(by='clan', as_index=False).agg({'og_id_domain': pd.Series.nunique}).sort_values('og_id_domain', ascending=False)
og_cnt_clan_initial['fraction'] = og_cnt_clan_initial['og_id_domain'].divide(initial_og_cnt)
og_cnt_clan_initial['group'] = 'initial'
#print("Top 10 clans initial dataset:", og_cnt_clan_initial.head(n=10)) 

view_clan_df = pd.concat([top_clans,og_cnt_clan_initial.loc[og_cnt_clan_initial['clan'].isin(top_clans['clan'])]])

#Plot of #OGs per clan, show the top 10
sns.barplot(data=view_clan_df, y='clan', x='fraction', hue='group')
axes.set_xlabel('%OGs')
plt.title('Top 10 clans in dataset')
plt.savefig(analysis_path+'plots/clans_initial_vs_final.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()

#Compare to positives
#positives, different or not?
pfam_detailed_df_pos = pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'].isin(positive_set['og_id_domain'])]
og_cnt_clan_positive = pfam_detailed_df_pos.groupby(by='clan', as_index=False).agg({'og_id_domain': pd.Series.nunique}).sort_values('og_id_domain', ascending=False)
og_cnt_clan_positive['fraction'] = og_cnt_clan_positive['og_id_domain'].divide(len(positive_set.index))
og_cnt_clan_positive['group'] = 'positive'
top_clans=og_cnt_clan_positive.head(n=10) #top from positive

view_clan_df = pd.concat([top_clans,og_cnt_clan.loc[og_cnt_clan['clan'].isin(top_clans['clan'])]])
sns.barplot(data=view_clan_df, y='clan', x='fraction', hue='group', hue_order=['final','positive'])
axes.set_xlabel('%OGs')
plt.title('Top 10 clans in positive dataset')
plt.savefig(analysis_path+'plots/clans_final_vs_positive.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()

#to prevent going from 1 to 1 and being very enriched, minimum of fraction of 1%
#mean_og_cnt = og_cnt_clan_positive.og_id_domain.mean()
#sns.distplot(og_cnt_clan_positive.fraction)

enrichment_df = og_cnt_clan_positive.loc[og_cnt_clan_positive['fraction']>0.01].merge(og_cnt_clan, on='clan')
enrichment_df['enrichment'] = enrichment_df['fraction_x'] / enrichment_df['fraction_y']
enrichment_df = enrichment_df.sort_values(by='enrichment', ascending=False)
enrichment_df = enrichment_df.reset_index(drop=True)
top_enrichment = enrichment_df.head(n=5).clan
top_depleted = enrichment_df.tail(n=7).clan

view_clan_df = pd.concat( [ og_cnt_clan_positive.loc[og_cnt_clan_positive['clan'].isin(top_enrichment)],\
og_cnt_clan.loc[og_cnt_clan['clan'].isin(top_enrichment)],\
og_cnt_clan_positive.loc[og_cnt_clan_positive['clan'].isin(top_depleted)],\
og_cnt_clan.loc[og_cnt_clan['clan'].isin(top_depleted)],\
])

clan_order = top_enrichment.to_dict().values()
sns.barplot(data=view_clan_df, y='clan', x='fraction', hue='group', order=clan_order, hue_order=['final','positive'])
axes.set_xlabel('%OGs')
plt.title('Enriched in positive dataset')
#plt.savefig(analysis_path+'plots/clan_enrichment_positive.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()

clan_order = top_depleted.to_dict().values()[::-1]
sns.barplot(data=view_clan_df, y='clan', x='fraction', hue='group', order=clan_order, hue_order=['final','positive'])
axes.set_xlabel('%OGs')
plt.title('Depleted in positive dataset')
#plt.savefig(analysis_path+'plots/clan_depletion_positive.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
#'''


'''
#Unit count distribution for clans in top 10
sns.violinplot(data=pfam_detailed_df.loc[pfam_detailed_df['clan'].isin(top_clans)], x='unit_cnt', y='clan', cut=0, scale='count', order=top_clans)
axes.set_xlabel('#units in repeat')
plt.title('Unit count per clan')
plt.xscale('log')
plt.savefig(analysis_path+'plots/unit_cnt_per_clan.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()

#Model length distribution for clans in top 10
sns.violinplot(data=pfam_detailed_df.loc[pfam_detailed_df['clan'].isin(top_clans)], x='model_length', y='clan', cut=0, scale='count', order=top_clans)
axes.set_xlabel('model length in aa')
plt.title('Model length per clan')
plt.savefig(analysis_path+'plots/model_length_per_clan.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
#Protein length distribution for clans in top 10
sns.violinplot(data=pfam_detailed_df.loc[pfam_detailed_df['clan'].isin(top_clans)], x='protein_length', y='clan', cut=0, scale='count', order=top_clans)
axes.set_xlabel('protein length in aa')
plt.title('Protein length per clan')
plt.xscale('log')
plt.savefig(analysis_path+'plots/protein_length_per_clan.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
'''

'''
#repeat cov per species
sns.violinplot(data=pfam_detailed_df, x='repeat_cov', y='species', cut=0, scale='count')
plt.title('Distribution: repeat cov per species')
#plt.xscale('symlog')
plt.savefig(analysis_path+'plots/repeat_cov_per_species.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''	

'''
#Distributions for full dataset initial
distributions = ['unit_cnt','model_length','protein_length','protein_cov','repeat_cov']
#distributions = ['protein_cov','repeat_cov']
for dist in distributions:
	#sns.violinplot(data=initial_detailed_df, x=dist, cut=0, scale='count')
	ax1 = sns.distplot(initial_detailed_df[dist], color='blue', label='initial')
	ax2 = sns.distplot(pfam_detailed_df[dist], color='green', label='final')
	plt.title('Compare distribution: '+dist)
	#if initial_detailed_df[dist].max() > 100*max(1,initial_detailed_df[dist].min()): 
		#plt.xscale('symlog')
	plt.legend()
	plt.savefig(analysis_path+'plots/compare_dist_'+dist+'.png',format='png',bbox_inches='tight')
	
	plt.show()
	plt.clf()
'''

'''
#Distribution of unit count stdev

pfam_std_list = []
initial_std_list = []
for og_id_domain in pfam_summary.keys():
	if initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.mean() < 3: continue
	std_pfam = pfam_detailed_df.loc[pfam_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.std()
	if std_pfam > 0:
		pfam_std_list.append(float(std_pfam))
	else: continue
	std_initial = initial_detailed_df.loc[initial_detailed_df['og_id_domain'] == og_id_domain].unit_cnt.std()
	initial_std_list.append(float(std_initial))

#sns.distplot(initial_std_list, color='blue')
sns.distplot([20 if x>20 else x for x in pfam_std_list if x>0], color='green', label='final',kde=False)
sns.distplot([20 if x>20 else x for x in initial_std_list if x>0], color='blue', label='initial',kde=False)
plt.title('Compare distribution: std unit count')
#if initial_detailed_df[dist].max() > 100*max(1,initial_detailed_df[dist].min()): 
	#plt.xscale('symlog')
#plt.axvline(x=np.mean(initial_std_list), color='blue')
#plt.axvline(x=np.mean(pfam_std_list), color='green')

plt.legend()
plt.savefig(analysis_path+'plots/compare_unit_std.png',format='png',bbox_inches='tight')
plt.show()

plt.clf()
'''

#try max dist
#print(pfam_detailed_df['max_dist'].loc[pfam_detailed_df['max_dist']<20].value_counts())
'''

pfam_detailed_df['max_dist'] = np.where(pfam_detailed_df['max_dist'] > 50 , 50,pfam_detailed_df['max_dist'])
initial_detailed_df['max_dist'] = np.where(initial_detailed_df['max_dist'] > 50 , 50,initial_detailed_df['max_dist'])

weights = np.ones_like(np.array(initial_detailed_df['max_dist']))/float(len(np.array(initial_detailed_df['max_dist'])))
plt.hist(initial_detailed_df['max_dist'],weights=weights,color='blue',label='initial',bins=20,alpha=0.6)
weights = np.ones_like(np.array(pfam_detailed_df['max_dist']))/float(len(np.array(pfam_detailed_df['max_dist'])))
plt.hist(pfam_detailed_df['max_dist'],weights=weights,color='green',label='final',bins=20,alpha=0.6)
plt.legend()
#plt.yscale('symlog')
#plt.savefig(analysis_path+'plots/compare_repeat_maxdist.png',format='png',bbox_inches='tight')
plt.savefig(analysis_path+'plots/compare_repeat_maxdist_cut50.png',format='png',bbox_inches='tight')
plt.xlabel('maxdist')
plt.show()

'''

#'''
#Model length distribution for zinc fingers only 
for dist in distributions:
	ax1 = sns.distplot(initial_detailed_df['model_length'], color='blue', label='initial')
	ax2 = sns.distplot(pfam_detailed_df['model_length'], color='green', label='final')
	plt.title('Zinc fingers model length distribution')
	plt.legend()
	plt.savefig(analysis_path+'plots/model_length_znfingers.png',format='png',bbox_inches='tight')
	plt.show()
	plt.clf()
#'''
#Plot model length difference initial vs final
merged = pfam_detailed_df.merge(initial_detailed_df, on='og_id_domain')
merged['model_difference'] = merged['model_length_x']-merged['model_length_y']
print(merged.groupby(by='og_id_domain',as_index=False)['model_difference'])
sns.distplot(merged.groupby(by='og_id_domain')['model_difference'])
plt.show()


