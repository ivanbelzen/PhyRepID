# IvB 
# Edit 03-09-2018 
# Input: pfam_evo_events.json meme_evo_events.json
# PRE evo events summary contains dup/loss per node/ortholog in an OG and in total 

import sys,json
import pandas as pd
import matplotlib.pyplot as plt
import pipeline_methods as pre
import numpy as np
import seaborn as sns
from scipy.stats.mstats import mode, gmean, hmean
import matplotlib.ticker as ticker

#analysis for duplication score 1
#analysis_path = pre.root+'analysis/evo_events/reroot/'
#pfam_summary_df_file = 'pfam_summary_df_reroot.json'
#meme_summary_df_file = 'meme_summary_df_reroot.json'
#combined_summary_df_file = 'combined_summary_df_reroot.json'

#Combined summary df now made in analyse_merge_export.py
analysis_path = pre.root+'analysis/evo_events_06102018/'
sys.stdout = open(analysis_path+'output.txt', 'w')

#File names for loading/saving dataframes
combined_summary_df_file = analysis_path+'pre_dataset_df.json'
selectome_df_file = 'selectome_df.json'
#genetree_df_file = 'genetree_df.json' #now just pfam but should be both
genetree_df_file = 'genetree_df_reroot.json'

'''
#Make dataframes
pfam_results_file = 'pfam_evo_events_reroot.json'
meme_results_file = 'meme_evo_events_reroot.json'
pfam_summary_file = 'pfam_summary.json'
meme_summary_file = 'meme_summary.json'
human_prot_selected_file = 'human_prot_selected.tsv'

pfam_gt_results_file = 'pfam_evo_events_gt.json'

with open(pfam_gt_results_file, 'r') as results:
	genetree_results = json.load(results)

with open(pfam_results_file, 'r') as results:
	pfam_results = json.load(results)

with open(meme_results_file, 'r') as results:
	meme_results = json.load(results)

with open(pfam_summary_file, 'r') as pfam_output:
	pfam_summary = json.load(pfam_output)
	
with open(meme_summary_file, 'r') as meme_output:
	meme_summary = json.load(meme_output)

with open(human_prot_selected_file, 'r') as human_prot_selected:
	selectome = []
	for line in human_prot_selected:
		cols = line.split('\t')
		selectome.append(cols[0])
'''

'''		
clans_dict = pre.get_pfam_clans(pre.pfam_clans_file)
#clans_dict = pre.get_pfam_clans('Pfam-A.clans.tsv')

genetree_df_cols = ['og_id', 'gt_dup', 'gt_loss']
genetree_df_rows = []

pfam_summary_df_cols = ['og_id_domain', 'og_id', 'pre_model', 'clan', 'netto_dup', 'loss', 'orthologs_cnt']
pfam_summary_df_rows = []

meme_summary_df_cols = ['og_id', 'netto_dup', 'loss', 'orthologs_cnt']
meme_summary_df_rows = []

duplications_df_cols = ['og_id_domain', 'pre_model', 'clan', 'node', 'dup_cnt']
duplications_df_rows = []

selectome_df_cols = ['og_id_domain', 'og_id', 'protein_id', 'positive_selection']
selectome_df_rows = []

pfam_repeat_data = pfam_results['repeat']
meme_repeat_data = meme_results['repeat']
genetree_results = genetree_results['genetree']


for og_id_domain,data in pfam_repeat_data.items(): 
	og_id = '_'.join(og_id_domain.split('_')[0:2])
	pre_model = '_'.join(og_id_domain.split('_')[2:])
	clan = clans_dict[pre_model] if pre_model in clans_dict else ''
	netto_dup = data['netto_dup']
	loss = data['loss']

	if og_id_domain not in pfam_summary: continue
	orthologs_list = pfam_summary[og_id_domain][pre_model]['orthologs_dict'].keys()
	#if len(orthologs_list)<5: continue
	human_proteins = [x for x in orthologs_list if x.startswith('ENSP0')]
	for human_protein in human_proteins:
		if human_protein in selectome: selectome_df_rows.append([og_id_domain, og_id, human_protein,1])
		else: selectome_df_rows.append([og_id_domain, og_id, human_protein,0])
		
	pfam_summary_df_rows.append( [og_id_domain, og_id, pre_model, clan, netto_dup, loss, len(orthologs_list)] )
	
	if 'dup' in data['events']:
		for node,cnt in data['events']['dup'].items():
			node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))
			duplications_df_rows.append([og_id_domain,pre_model,clan,node,cnt])		

pfam_summary_df = pd.DataFrame(pfam_summary_df_rows, columns=pfam_summary_df_cols)

for og_id,data in meme_repeat_data.items(): 
	netto_dup = data['netto_dup']
	loss = data['loss']
	
	if og_id not in meme_summary: continue
	orthologs_list = meme_summary[og_id]['motif']['orthologs_dict'].keys()
	#if len(orthologs_list)<4: continue
	human_proteins = [x for x in orthologs_list if x.startswith('ENSP0')]
	for human_protein in human_proteins:
		if human_protein in selectome: selectome_df_rows.append(['', og_id, human_protein,1])
		else: selectome_df_rows.append(['', og_id, human_protein,0])
	
	meme_summary_df_rows.append( [og_id, netto_dup, loss, len(orthologs_list)] )
	
	if 'dup' in data['events']:
		for node,cnt in data['events']['dup'].items():
			node = node if filter(str.isalpha,str(node)) == 't' else filter(str.isalpha,str(node))
			duplications_df_rows.append([og_id,'','motif',node,cnt])

meme_summary_df = pd.DataFrame(meme_summary_df_rows, columns=meme_summary_df_cols)

#duplications_df = pd.DataFrame(duplications_df_rows, columns=duplications_df_cols)
selectome_df = pd.DataFrame(selectome_df_rows, columns=selectome_df_cols)

for og_id,data in genetree_results.items(): 
	netto_dup = data['netto_dup']
	loss = data['loss']
	genetree_df_rows.append([og_id,netto_dup,loss])
	
genetree_df = pd.DataFrame(genetree_df_rows, columns=genetree_df_cols)

#Write dataframes to file for faster running of script

with open(meme_summary_df_file, 'w') as output:
	output.write(json.dumps(meme_summary_df.to_dict()))	

with open(pfam_summary_df_file, 'w') as output:
	output.write(json.dumps(pfam_summary_df.to_dict()))	
		
with open(selectome_df_file, 'w') as output:
	output.write(json.dumps(selectome_df.to_dict()))

with open(duplications_df_file, 'w') as output:
	output.write(json.dumps(duplications_df.to_dict()))	
with open(genetree_df_file, 'w') as output:
	output.write(json.dumps(genetree_df.to_dict()))	
'''


#Deprecated, make combined summary df from pfam and meme
'''
#Read dataframes from file
with open(selectome_df_file, 'r') as output:
	selectome = json.load(output)
	selectome_df = pd.DataFrame.from_dict(selectome)
with open(pfam_summary_df_file, 'r') as output:
	pfam_summary = json.load(output)
	pfam_summary_df = pd.DataFrame.from_dict(pfam_summary)
with open(meme_summary_df_file, 'r') as output:
	meme_summary = json.load(output)
	meme_summary_df = pd.DataFrame.from_dict(meme_summary)
with open(genetree_df_file, 'r') as output:
	genetree = json.load(output)
	genetree_df = pd.DataFrame.from_dict(genetree)
with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)

combined_summary_df = pd.concat([pfam_summary_df,meme_summary_df], ignore_index=True, sort=False)
combined_summary_df['clan'] = np.where(combined_summary_df['clan'].isna(), 'motif', combined_summary_df['clan'])
mean_combined_dup = combined_summary_df['netto_dup'].mean()
combined_summary_df['duplication_score']= (combined_summary_df['netto_dup']-mean_combined_dup)/combined_summary_df['orthologs_cnt']

'''	


#Use combined dataframe from analyse merge export .py 
#'''


with open(combined_summary_df_file, 'r') as output:
	combined_summary = json.load(output)
	combined_summary_df = pd.DataFrame.from_dict(combined_summary)

stats = combined_summary_df.describe()
std=stats['duplication_score']['std']

pfam_subset = combined_summary_df.loc[combined_summary_df['clan']!='motif']
meme_subset = combined_summary_df.loc[combined_summary_df['clan']=='motif']

sys.stdout = open(analysis_path+'dataframe_stats.txt', 'w')
print(combined_summary_df_file)
print('#OGs in total',len(combined_summary_df.index))
print('#OGs in pfam',len(pfam_subset.index))
print('#OGs in meme',len(meme_subset.index))
print('#gene trees in total',combined_summary_df.genetree_id.nunique())
print('#gene trees in pfam',pfam_subset.genetree_id.nunique())
print('#gene trees in meme',meme_subset.genetree_id.nunique())

#quit()
## Compare different scoring measures
'''
sys.stdout = open(analysis_path+'comparing_scores.txt', 'w')

std2=stats['normalized_score']['std']
std3=stats['expectancy_score']['std']

print("Compare different scoring measures")
print("Duplication score")
print("Netto dup mean",stats['netto_dup']['mean'])
print('#OGs >> netto dup - mean / orth cnt ', \
len(combined_summary_df.loc[combined_summary_df['duplication_score']>0].index), \
'>1stdev above mean',\
len(combined_summary_df.loc[combined_summary_df['duplication_score']>std].index))

print("Normalized score")
print("Scaled dup mean",stats['scaled_dup']['mean'])
print('#OGs >> netto dup / orth cnt - scaled mean ', \
len(combined_summary_df.loc[combined_summary_df['normalized_score']>0].index), \
'>1stdev above mean',\
len(combined_summary_df.loc[combined_summary_df['normalized_score']>std2].index))

print("Expectancy score")
print("Netto dup per protein",stats['scaled_dup']['mean'])
print('#OGs >> netto dup - (scaled mean * orth cnt) ', \
len(combined_summary_df.loc[combined_summary_df['expectancy_score']>0].index), \
'>1stdev above mean',\
len(combined_summary_df.loc[combined_summary_df['expectancy_score']>std3].index))
'''

## Histograms of distribution of duplications
'''
#Pfam and meme separately plot scores
pfam_subset = combined_summary_df.loc[combined_summary_df['clan']!='motif']
meme_subset = combined_summary_df.loc[combined_summary_df['clan']=='motif']

weights = np.ones_like(np.array(pfam_subset['duplication_score']))/float(len(np.array(pfam_subset['duplication_score'])))
plt.hist(pfam_subset['duplication_score'],weights=weights,color=pre.blue_colour,bins=300,label='pfam')
weights = np.ones_like(np.array(meme_subset['duplication_score']))/float(len(np.array(meme_subset['duplication_score'])))
plt.hist(meme_subset['duplication_score'],weights=weights,color=pre.green_colour,bins=300,label='meme')
plt.xscale('symlog')
plt.legend()
plt.savefig(analysis_path+'duplication_score_pfam_meme.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()


weights = np.ones_like(np.array(pfam_subset['normalized_score']))/float(len(np.array(pfam_subset['normalized_score'])))
plt.hist(pfam_subset['normalized_score'],weights=weights,color=pre.blue_colour,bins=300,label='pfam')
weights = np.ones_like(np.array(meme_subset['normalized_score']))/float(len(np.array(meme_subset['normalized_score'])))
plt.hist(meme_subset['normalized_score'],weights=weights,color=pre.green_colour,bins=300,label='meme')
plt.xscale('symlog')
plt.legend()
plt.savefig(analysis_path+'normalized_score_pfam_meme.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()


weights = np.ones_like(np.array(pfam_subset['expectancy_score']))/float(len(np.array(pfam_subset['expectancy_score'])))
plt.hist(pfam_subset['expectancy_score'],weights=weights,color=pre.blue_colour,bins=300,label='pfam')
weights = np.ones_like(np.array(meme_subset['expectancy_score']))/float(len(np.array(meme_subset['expectancy_score'])))
plt.hist(meme_subset['expectancy_score'],weights=weights,color=pre.green_colour,bins=300,label='meme')
plt.xscale('symlog')
plt.legend()
plt.savefig(analysis_path+'expectancy_score_pfam_meme.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

#Full dataset
weights = np.ones_like(np.array(combined_summary_df['duplication_score']))/float(len(np.array(combined_summary_df['duplication_score'])))
plt.hist(combined_summary_df['duplication_score'],weights=weights,color=pre.blue_colour,bins=420)
plt.xscale('symlog')
plt.savefig(analysis_path+'duplication_score.png',format='png',bbox_inches='tight')
#Combined with RUG
sns.rugplot(combined_summary_df['duplication_score'], color=pre.green_colour, height=0.01)
plt.savefig(analysis_path+'duplication_score_rug.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

weights = np.ones_like(np.array(combined_summary_df['normalized_score']))/float(len(np.array(combined_summary_df['normalized_score'])))
plt.hist(combined_summary_df['normalized_score'],weights=weights,color=pre.blue_colour,bins=420)
plt.xscale('symlog')
plt.savefig(analysis_path+'normalized_score.png',format='png',bbox_inches='tight')
#Combined with RUG
sns.rugplot(combined_summary_df['normalized_score'], color=pre.green_colour, height=0.01)
plt.savefig(analysis_path+'normalized_score_rug.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()

weights = np.ones_like(np.array(combined_summary_df['expectancy_score']))/float(len(np.array(combined_summary_df['expectancy_score'])))
plt.hist(combined_summary_df['expectancy_score'],weights=weights,color=pre.blue_colour,bins=420)
plt.xscale('symlog')
plt.savefig(analysis_path+'expectancy_score.png',format='png',bbox_inches='tight')
#Combined with RUG
sns.rugplot(combined_summary_df['expectancy_score'], color=pre.green_colour, height=0.01)
plt.savefig(analysis_path+'expectancy_score_rug.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

#Netto dup distribution
'''
weights = np.ones_like(np.array(pfam_subset['netto_dup']))/float(len(np.array(pfam_subset['netto_dup'])))
plt.hist(pfam_subset['netto_dup'],weights=weights,color=pre.blue_colour,bins=50,label='pfam',alpha=0.6)
weights = np.ones_like(np.array(meme_subset['netto_dup']))/float(len(np.array(meme_subset['netto_dup'])))
plt.hist(meme_subset['netto_dup'],weights=weights,color=pre.green_colour,bins=15,label='meme',alpha=0.6)
plt.yscale('log')
plt.title('Net duplication distribution')
plt.xlabel('# net duplications')
plt.ylabel('relative frequency')
plt.legend()
plt.savefig(analysis_path+'netto_dup_pfam_meme.png',format='png',bbox_inches='tight')
plt.clf()
'''
# Distribution with examples

controls = {'PRDM9': 'ENSGT00530000063157_ENSG00000164256_zf-H2C2_2',\
'KNL1': 'ENSGT00410000025918_ENSG00000137812',\
'VCAM1': 'ENSGT00830000128299_ENSG00000162692_Ig_3',\
 'CDC20': 'ENSGT00870000136444_ENSG00000117399_WD40', \
 'AHNAK': 'ENSGT00530000063716_ENSG00000124942_DUF945', \
 'PRX': 'ENSGT00530000063716_ENSG00000105227_Cornifin', \
 'BRCA2': 'ENSGT00390000003602_ENSG00000139618_BRCA2' }

colour_list = ['#00b159','#bf4dea','#f37735','#ce1b2e']

print("Duplication score")

control_values = {}
for gene,og_id in controls.items():
	print(gene)
	in_dataset = combined_summary_df.loc[combined_summary_df['identifier']==og_id]
	if (len(in_dataset.index)>1): in_dataset = in_dataset.iloc[0]# = pfam_subset.loc[pfam_subset['identifier']==og_id]
	if not in_dataset.empty: control_values[gene]= float(in_dataset['duplication_score'])

print(control_values)
plt.figure(figsize=(8,4))
i=0
weights = np.ones_like(np.array(combined_summary_df['duplication_score']))/float(len(np.array(combined_summary_df['duplication_score'])))
plt.hist(combined_summary_df['duplication_score'],weights=weights,color=pre.blue_colour,bins=420)
plt.xscale('symlog')
for label, value in control_values.items():	plt.axvline(x=value, color=colour_list[i% len(colour_list)], label=label); i+=1
plt.legend()
plt.ylabel('relative frequency')
plt.title("Duplication score landscape")
plt.xlabel("duplication score")
plt.savefig(analysis_path+'duplication_score_examples.png',format='png',bbox_inches='tight')
plt.savefig(analysis_path+'duplication_score_examples.svg',format='svg',bbox_inches='tight')
plt.clf()

quit()

'''
print("Normalized score")

control_values = {}
for gene,og_id in controls.items():
	in_dataset = combined_summary_df.loc[combined_summary_df['identifier']==og_id]
	if not in_dataset.empty: control_values[gene]= float(in_dataset['normalized_score'])

print(control_values)

i=0
weights = np.ones_like(np.array(combined_summary_df['normalized_score']))/float(len(np.array(combined_summary_df['normalized_score'])))
plt.hist(combined_summary_df['normalized_score'],weights=weights,color=pre.blue_colour,bins=420)
plt.xscale('symlog')
for label, value in control_values.items():	plt.axvline(x=value, color=colour_list[i% len(colour_list)], label=label); i+=1
plt.legend()
plt.savefig(analysis_path+'normalized_score_examples.png',format='png',bbox_inches='tight')
plt.savefig(analysis_path+'normalized_score_examples.svg',format='svg',bbox_inches='tight')
plt.clf()

print("Expectancy score")

control_values = {}
for gene,og_id in controls.items():
	in_dataset = combined_summary_df.loc[combined_summary_df['identifier']==og_id]
	if not in_dataset.empty: control_values[gene]= float(in_dataset['expectancy_score'])

print(control_values)

i=0
weights = np.ones_like(np.array(combined_summary_df['expectancy_score']))/float(len(np.array(combined_summary_df['expectancy_score'])))
plt.hist(combined_summary_df['expectancy_score'],weights=weights,color=pre.blue_colour,bins=420)
plt.xscale('symlog')
for label, value in control_values.items():	plt.axvline(x=value, color=colour_list[i% len(colour_list)], label=label); i+=1
plt.legend()
plt.savefig(analysis_path+'expectancy_score_examples.png',format='png',bbox_inches='tight')
plt.savefig(analysis_path+'expectancy_score_examples.svg',format='svg',bbox_inches='tight')
plt.clf()
'''

## Correlate duplication score with unit CV
'''
sns.jointplot(data=combined_summary_df, x ='duplication_score', y='unit_cv',kind='reg')#,marginal_kws={'cut':0,'bw':0.1})
#plt.xscale('symlog')
plt.savefig(analysis_path+'dupl_score_vs_unit_cv.png',format='png',bbox_inches='tight')
plt.xscale('symlog')
plt.savefig(analysis_path+'dupl_score_vs_unit_cv_log.png',format='png',bbox_inches='tight')
#plt.show()
plt.clf()
'''

#Characterize positive set

#Clans, plot top 10 initial vs final pfam dataset, and vs positive set
'''
sys.stdout = open(analysis_path+'characterize_positive.txt', 'w')

positive_set = combined_summary_df.loc[combined_summary_df['duplication_score']>0]
fast_set = combined_summary_df.loc[combined_summary_df['duplication_score']>std]

og_cnt_clan = combined_summary_df.groupby(by='clan', as_index=False).agg({'identifier': pd.Series.nunique}).sort_values('identifier', ascending=False)
og_cnt_clan['fraction'] = og_cnt_clan['identifier'].divide(len(combined_summary_df.index))

og_cnt_clan_pos = positive_set.groupby(by='clan', as_index=False).agg({'identifier': pd.Series.nunique}).sort_values('identifier', ascending=False)
og_cnt_clan_pos['fraction'] = og_cnt_clan_pos['identifier'].divide(len(positive_set.index))
og_cnt_clan_pos['group'] = 'positive'
top_clans = og_cnt_clan_pos.head(n=5)
print("Top 10 clans positive:",og_cnt_clan_pos.head(n=10)) 

og_cnt_clan_pos_fast = fast_set.groupby(by='clan', as_index=False).agg({'identifier': pd.Series.nunique}).sort_values('identifier', ascending=False)
og_cnt_clan_pos_fast['fraction'] = og_cnt_clan_pos_fast['identifier'].divide(len(fast_set.index))
og_cnt_clan_pos_fast['group'] = 'fast'
print("Top 10 clans fast dataset:", og_cnt_clan_pos_fast.head(n=10)) 

view_clan_df = pd.concat([top_clans,og_cnt_clan_pos_fast.loc[og_cnt_clan_pos_fast['clan'].isin(top_clans['clan'])]])

#Plot of #OGs per clan, show the top 10

' ''
fig,axes = plt.subplots()
sns.barplot(data=view_clan_df, y='clan', x='fraction', hue='group')
axes.set_xlabel('%OGs')
plt.title('Top 5 clans in positive dataset')
plt.savefig(analysis_path+'clans_positive_vs_fast.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
'''

## Calculate unit duplication score per clan
'''
sys.stdout = open(analysis_path+'unit_duplication_score.txt', 'w')
clan_selection = og_cnt_clan.loc[og_cnt_clan['fraction']>0.01]['clan']
unit_dupscore_df = combined_summary_df.groupby(by='clan')['unit_duplication_score'].agg([np.mean,np.std])
unit_dupscore_df.reset_index(inplace=True)
unit_dupscore_df = unit_dupscore_df.sort_values('mean', ascending=False)
print(unit_dupscore_df.loc[unit_dupscore_df['clan'].isin(clan_selection)].to_string())
'''
'''
sns.distplot(combined_summary_df.loc[combined_summary_df['clan']=='C2H2-zf']['unit_duplication_score'], color=pre.blue_colour,kde=True,kde_kws={'cut':0})
sns.distplot(positive_set.loc[positive_set['clan']=='C2H2-zf']['unit_duplication_score'], color=pre.green_colour,kde=True,kde_kws={'cut':0})
sns.distplot(fast_set.loc[fast_set['clan']=='C2H2-zf']['unit_duplication_score'], color=pre.purple_colour,kde=True,kde_kws={'cut':0})
plt.legend()
plt.savefig(analysis_path+'unit_duprate_zinc_finger.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
'''
'''
weights = np.ones_like(np.array(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score']))/float(len(np.array(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score'])))
plt.hist(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score'],weights=weights,color=pre.blue_colour,bins=50,label='top 10 clans',alpha=0.6)
weights = np.ones_like(np.array(combined_summary_df.loc[combined_summary_df['clan']=='C2H2-zf']['unit_duplication_score']))/float(len(np.array(combined_summary_df.loc[combined_summary_df['clan']=='C2H2-zf']['unit_duplication_score'])))
plt.hist(combined_summary_df.loc[combined_summary_df['clan']=='C2H2-zf']['unit_duplication_score'],weights=weights,color=pre.green_colour,bins=50,label='zinc fingers',alpha=0.6)
sns.rugplot(fast_set.loc[fast_set['clan']=='C2H2-zf']['unit_duplication_score'], color=pre.purple_colour,label='fast zf')
plt.legend()
plt.savefig(analysis_path+'unit_duprate_comparison_zf.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()

weights = np.ones_like(np.array(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score']))/float(len(np.array(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score'])))
plt.hist(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score'],weights=weights,color=pre.blue_colour,bins=50,label='top 10 clans',alpha=0.6)
weights = np.ones_like(np.array(combined_summary_df.loc[combined_summary_df['clan']=='motif']['unit_duplication_score']))/float(len(np.array(combined_summary_df.loc[combined_summary_df['clan']=='motif']['unit_duplication_score'])))
plt.hist(combined_summary_df.loc[combined_summary_df['clan']=='motif']['unit_duplication_score'],weights=weights,color=pre.green_colour,bins=50,label='motifs',alpha=0.6)
sns.rugplot(fast_set.loc[fast_set['clan']=='motif']['unit_duplication_score'], color=pre.purple_colour,label='fast motifs')
plt.legend()
plt.savefig(analysis_path+'unit_duprate_comparison_motif.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()


weights = np.ones_like(np.array(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score']))/float(len(np.array(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score'])))
plt.hist(combined_summary_df.loc[combined_summary_df['clan'].isin(clan_selection)]['unit_duplication_score'],weights=weights,color=pre.blue_colour,bins=50,label='top 10 clans',alpha=0.6)
weights = np.ones_like(np.array(combined_summary_df.loc[combined_summary_df['clan']=='Beta_propeller']['unit_duplication_score']))/float(len(np.array(combined_summary_df.loc[combined_summary_df['clan']=='Beta_propeller']['unit_duplication_score'])))
plt.hist(combined_summary_df.loc[combined_summary_df['clan']=='Beta_propeller']['unit_duplication_score'],weights=weights,color=pre.green_colour,bins=10,label='Beta propellers',alpha=0.6)
sns.rugplot(fast_set.loc[fast_set['clan']=='Beta_propeller']['unit_duplication_score'], color=pre.purple_colour,label='fast Beta propellers')
plt.legend()
plt.savefig(analysis_path+'unit_duprate_comparison_beta_propeller.png',format='png',bbox_inches='tight')
plt.show()
plt.clf()
'''

## Analyse genetree dataframe

with open(genetree_df_file,'r') as output:
	gt = json.load(output)
	genetree_df = pd.DataFrame.from_dict(gt)
	
#genetree_df = genetree_df.loc[genetree_df['og_id'].isin(combined_summary_df['og_id'])]	
sys.stdout = open(analysis_path+'genetree.txt', 'w')
print("Analysis gene trees")
print("Entries",len(genetree_df.index), "unique OGs",genetree_df.og_id.nunique())
print("OGs with gene duplication",genetree_df.loc[genetree_df['gt_dup']>0].og_id.nunique())
