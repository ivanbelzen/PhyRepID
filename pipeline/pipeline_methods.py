import csv,json, os
import numpy as np
root = '/home/ianthe/protein-repeat-evolution/'

species_mapping_file = root+'ensembl_stable_id_species.json'
species_tree_file = root+'ensembl_taxon_species_tree.phy'

pfam_gacutoff_file = root+'ga_cutoff.tsv'
pfam_clans_file = root+'Pfam-A.clans.tsv' 
pfam_file=root+'/pfam/Pfam-A.hmm'

taxon_uri =  'http://identifiers.org/taxonomy/'
gene_uri_prefix = 'http://rdf.ebi.ac.uk/resource/ensembl/' 

#Folder paths
genetree_path = root+'genetrees_nhx/'
ensembl_path = root+'ensembl_api/'     
fasta_path = root+'fasta_ogs/'
fasta_masked_path=root+'fasta_ogs_masked/'
fasta_chopped_path=root+'fasta_ogs_chopped/'
full_alignment_path = root+'full_alignment/'
n_alignment_path = root+'n_alignment/'
c_alignment_path = root+'c_alignment/'
maptotree_path = root+'maptotree/'

pfam_hmm_path=root+'pfam/hmm/'
pfam_aligned_path=root+'pfam/aligned/'
pfam_tblout_path=root+'pfam/tblout/'
pfam_profiles_path=root+'pfam/profiles/'
pfam_trees_path=root+'pfam/trees/'
pfam_treefix_path = root+'pfam/treefix/'
pfam_treefix_adjusted_path = root+'pfam/treefix/adjusted/'
pfam_treefix_adjusted_images_path = root+'pfam/treefix/adjusted/images/'
#pfam_alignments_path=root+'pfam/alignments/'
progress_path=root+'progress_files/'

meme_output_path = root+'denovo/meme/'
denovo_meme_repeats_path = root+'denovo/meme/repeats/'
#denovo_repeats_path = root+'denovo/repeats/'
denovo_tblout_path=root+'denovo/tblout/'
denovo_profiles_path=root+'denovo/profiles/'
denovo_aligned_path=root+'denovo/aligned/'
#denovo_final_alignments_path=root+'denovo/alignments/'
denovo_trees_path=root+'denovo/trees/'
denovo_treefix_path = root+'denovo/treefix/'
denovo_treefix_adjusted_path = root+'denovo/treefix/adjusted/'
denovo_treefix_adjusted_images_path = root+'denovo/treefix/adjusted/images/'

#Mapping files
orthologs_file = root+'orthologs.json'
orthologs_filtered_file = root+'orthologs_filtered.json' 					#output of filter_orthologs.py
genes_to_genetrees_file = root+'genes_to_genetrees.json'	#output of retrieve_genetrees.py
ogs_to_hits_file = root+'hmm_results_final.json'			#output of iterative domain detection
hmm_results_file = root+'hmm_results.json'

### Output files ###

#Pipeline result files
phyrepid_results_simple=root+'phyrepid_results_simple.tsv' 
phyrepid_results_full=root+'phyrepid_results_full.tsv'
phyrepid_results_human_full_lineage = root+'phyrepid_results_human_full_lineage.tsv'
phyrepid_results_human_only = root+'phyrepid_results_human_only.tsv'

#repeat stats
pfam_repeat_stats_file = root+'pfam_repeat_stats.json'
meme_repeat_stats_file = root+'meme_repeat_stats.json'
pfam_repeat_stats_initial_file = root+'pfam_repeat_stats_initial.json'
pfam_repeat_stats_initial_filtered_file = root+'pfam_repeat_stats_initial_filtered.json'
		
#evo events
pfam_evo_events_file = root+'pfam_evo_events.json'
meme_evo_events_file = root+'meme_evo_events.json'
pfam_evo_gt_file = root+"pfam_evo_events_genetrees.json"
meme_evo_gt_file = root+"meme_evo_events_genetrees.json"

# Additional datasets - files not required 
selectome_file = root+'resources/selectome.tsv' #needs human mapping 
human_mapping_file = root+'og_human_mapping.json' #OG to human gene mapping from generate_human_protein_mapping.py
exac_file = root+'resources/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt' #needs gene symbol
human_gene_symbol_file = root+'resources/human_gene_symbol.tsv'
schaper_input_file_1 = root+'resources/eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick.gz',
schaper_input_file_2 = root+'resources/eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick.gz'
schaper_comparison_file = root+'schaper_comparison.json' #output from analyse_schaper_comparison.py
		

### ###

root_nodes = ['t117571','t7742', 't7711', 't33213']
node_projection = { 't32525': ['t9347', 't1437010', 't314147'], 't8287': ['t32523'], 't186625':['t186626'], 't41665': ['t1489872'], 't117571':['t7742', 't7711', 't33213']}


#Constants
blue_colour = '#00aedb'
green_colour = '#00b159'
purple_colour = '#bf4dea'

def file_notempty(filepath):
	if(os.path.exists(filepath) and os.path.getsize(filepath) > 0):
		return True
	else:
		return False

def get_taxon_list(species_file_path = '', exclude_list = []):
	if species_file_path == '': species_file_path = species_file
	taxon_list = []		
	
	with open(species_file_path,'r') as sp_mapping:
		species_dict = json.load(sp_mapping)
		for item in species_dict:
			if item["taxon_uri"] in exclude_list: continue
			taxon_list.append( item["taxon_uri"] ) #it is actually an id...
	return(taxon_list)

def get_species_group_list(species_file):
	group_list = {}
	with open(species_file,'r') as sp_mapping:
		species_dict = json.load(sp_mapping)
		for item in species_dict:
			if "group_id" in item:
				if item["group_id"] in group_list:
					group_list[item["group_id"]].append(item["taxon_uri"])
				else: 
					group_list[item["group_id"]] = [item["taxon_uri"]]				
	return(group_list)

def filter_species_shared(mapping, threshold):
	for gene_uri in list(mapping):
		if len(mapping[gene_uri]) < threshold:
			mapping.pop(gene_uri)
			continue
	return(mapping)

def filter_min_proteins(mapping, threshold):
	for gene_uri in list(mapping):
		if sum([len(mapping[gene_uri][x]) for x in mapping[gene_uri]]) < threshold:
			mapping.pop(gene_uri)
			continue
	return(mapping)


def filter_orth_rel(mapping, threshold, return_disjunct = False ):
	disjunct_dict = {} #not just removed taxons but entire gene
	for gene_uri in list(mapping):	
		for taxon in list(mapping[gene_uri]):
			if gene_uri in mapping and len(mapping[gene_uri][taxon]) > threshold: 
				disjunct_dict[gene_uri]=mapping[gene_uri]
									
				#mapping[gene_uri].pop(taxon)
				mapping.pop(gene_uri)
				continue
	if return_disjunct: return(mapping,disjunct_dict)
	return(mapping)
	
	
def filter_close_pairs(mapping, pairs, return_disjunct = False ):
	disjunct_dict = {} 
	#Pairs restriction means that either gene in both or in none
	for gene_uri in list(mapping):
		for species in pairs:
			if gene_uri in mapping and any(taxon_uri+str(x) in mapping[gene_uri] for x in species) and not all(taxon_uri+str(x) in mapping[gene_uri] for x in species):
				disjunct_dict[gene_uri]= mapping[gene_uri]
				mapping.pop(gene_uri)
			continue
	if return_disjunct: return(mapping,disjunct_dict)	
	return(mapping)


def check_list_str(li, st):
	return(all(e in st for e in li))

	
def get_nhxtree(gene_tree):
	if 'sequence' in gene_tree: #leaf node
		return( gene_tree['sequence']['id'][0]['accession']+':'+str(gene_tree['branch_length'])+'[&&NHX:T='+str(gene_tree['taxonomy']['id'])+']' )
		
	elif 'children' in gene_tree:
		nhx_tag = '[&&NHX:T='+str(gene_tree['taxonomy']['id'])+']'
		internal_node = 't'+str(gene_tree['taxonomy']['id'])
		
		if gene_tree['events']['type'] == 'duplication':
			nhx_tag = nhx_tag[:-1]+':D=Y]'
		#elif gene_tree['events']['type'] == 'dubious':
		#	nhx_tag = '[&&NHX:D=D]'
		elif gene_tree['events']['type'] == 'speciation':
			nhx_tag = nhx_tag[:-1]+':D=N]'
			
		return('('+get_nhxtree(gene_tree['children'][0])+','+get_nhxtree(gene_tree['children'][1])+')'+internal_node+':'+str(gene_tree['branch_length'])+nhx_tag)
		
	elif 'tree' in gene_tree:
		return(get_nhxtree(gene_tree['tree']))
	
	
	return(None)


def make_valid_nhx(tree):
	if tree[0] =='\"':
		tree = tree[1:]
	if tree[-1] != ');':
		split_tree = tree.split(')')
		partial_tree = ')'.join(split_tree[:-1])+')' #split on last )
		tree = partial_tree+split_tree[-1:][0].split(':')[0]+';' #last element, split on first occurance of : to get the internal node name but not length and NHX tag
		#tree = ':'.join(tree.split(':')[:-3]) + ';' 	#remove trailing branch length	and nhx tag
		
	return(tree)
			
def get_subtrees(gene_tree,orthologs_list):
	event_time = 440 #look for first speciation node with more recent timing
	
	if not check_list_str(orthologs_list, json.dumps(gene_tree)):
		return(None)
		
	if 'tree' in gene_tree:
		return(get_subtrees(gene_tree['tree'],orthologs_list))

	if 'children' in gene_tree:
		if 'timetree_mya' not in gene_tree['taxonomy']: return(None)
		
		if gene_tree['taxonomy']['timetree_mya'] > event_time:
			#can only be in one of the two children branches.
			if check_list_str(orthologs_list, json.dumps(gene_tree['children'][0])): return(get_subtrees(gene_tree['children'][0], orthologs_list))
			elif check_list_str(orthologs_list, json.dumps(gene_tree['children'][1])): return(get_subtrees(gene_tree['children'][1], orthologs_list))
		
		if gene_tree['taxonomy']['id'] == 117571 and gene_tree['events']['type'] == 'duplication': #bony vertebrates node
			if gene_tree['children'][0]['taxonomy']['id'] == 117571 and check_list_str(orthologs_list, json.dumps(gene_tree['children'][0])): return(get_subtrees(gene_tree['children'][0], orthologs_list))
			elif gene_tree['children'][1]['taxonomy']['id'] == 117571 and check_list_str(orthologs_list, json.dumps(gene_tree['children'][1])): return(get_subtrees(gene_tree['children'][1], orthologs_list))
			else: return (gene_tree)
		
		return(gene_tree)
			
	return(None)

	
def get_proteinseq(gene_tree):
	if 'sequence' not in json.dumps(gene_tree): return(None)
		
	if 'sequence' in gene_tree: #leaf node
		return('>'+gene_tree['sequence']['id'][0]['accession']+'\n'+gene_tree['sequence']['mol_seq']['seq']+'\n')
	elif 'children' in gene_tree:
		return(get_proteinseq(gene_tree['children'][0])+'\n'+get_proteinseq(gene_tree['children'][1]))
	elif 'tree' in gene_tree:
		return(get_proteinseq(gene_tree['tree']))
		
	return(None)

	
def get_identifiers(gene_tree):
	if 'sequence' not in json.dumps(gene_tree): return(None)
		
	if 'sequence' in gene_tree: #leaf node
		return('('+gene_tree['id']['accession']+','+gene_tree['sequence']['id'][0]['accession']+')')
	elif 'children' in gene_tree:
		return(get_identifiers(gene_tree['children'][0])+','+get_identifiers(gene_tree['children'][1]))
	elif 'tree' in gene_tree:
		return(get_identifiers(gene_tree['tree']))
		
	return(None)


def read_fasta(fasta_file):
	fasta = {}
	buffer = ""
	protein_uri = None 
	with open(fasta_file, 'r') as seq:	
		for row in seq:
			if row[0] == '>':
				if buffer != "": 
					fasta[protein_uri] = buffer
					buffer = ""
				protein_uri = row[1:].strip()
				
			else: 
				buffer += row.strip()
		fasta[protein_uri] = buffer
	return(fasta)

def get_pfam_clans(pfam_clans_file):
	pfam_clans = {}
	#read pfam clans
	with open(pfam_clans_file, 'r') as pfam_clans_file:
		pfam_clan_data = csv.reader(pfam_clans_file, delimiter='\t')
		for row in pfam_clan_data:
			if row[2] is "": row[2] = row[3]
			pfam_clans[row[3]] = row[2] #domain -> clan
	return(pfam_clans)

def get_pfam_gacutoff(pfam_gacutoff_file):	
	pfam_gacutoff = {}
	#read gathering cutoffs
	with open(pfam_gacutoff_file, 'r') as pfam_gacutoff_file:
		pfam_gacutoff_data = csv.reader(pfam_gacutoff_file, delimiter='\t')
		for row in pfam_gacutoff_data:
			pfam_gacutoff[row[0]] = float(row[1]) #acc -> cutoff
	return(pfam_gacutoff)

def get_species_mapping(species_mapping_file):
	species_mapping = {}		
	#read species mapping
	with open(species_mapping_file,'r') as sp_mapping:
		species_dict = json.load(sp_mapping)
		for item in species_dict:
			species_mapping[item["ensembl_stable_id"]+'P'] = item["abbr"]
	return(species_mapping)

def get_species_mapping_full(species_mapping_file):
	species_mapping = {}		
	#read species mapping
	with open(species_mapping_file,'r') as sp_mapping:
		species_dict = json.load(sp_mapping)
		for item in species_dict:
			species_mapping[item["ensembl_stable_id"]+'P'] = {}
			for attr,val in item.items():
				if 'ensembl_stable_id' in attr: continue
				species_mapping[item["ensembl_stable_id"]+'P'][attr] = val
	return(species_mapping)
		
def parse_domtblout(hmm_file, state = 'init', spacing = 'ali'):
	pfam_clans = get_pfam_clans(pfam_clans_file)
	pfam_gacutoff = get_pfam_gacutoff(pfam_gacutoff_file)
	species_mapping = get_species_mapping(species_mapping_file)
	
	repeats = {} #dict with detected repeats for pfam accessions per clan, with summerized sequence score over proteins and coordinates per repeat
	pfam_hits = {} #unique pfam accessions
	
	with open(hmm_file,'r') as hmmr_tbl:
		for row in hmmr_tbl:
			if row[:1] == '#': continue; #skip header rows
			cols = row.split() #0 gene uri, 3 HMM match 'hit', 9 number, 12 i-e
			
			pfam_accession = cols[1]
			sequence_score = float(cols[7])
			dom_score = float(cols[13])
			
			if state is 'init' and pfam_accession is not '-':
				if pfam_gacutoff[pfam_accession] > sequence_score: 
					#print(pfam_gacutoff[pfam_accession], sequence_score)
					continue; #skip sequences which are not according to threshold
			if dom_score < 0: continue #skip domains < 0
			
			hit, protein_uri = cols[0], cols[3] #vice versa in HMMsearch..
			identifier = protein_uri+'_'+hit	#followed by lead number
			ensembl_stable_id = filter(str.isalpha, protein_uri)
			species = species_mapping[ensembl_stable_id]
			if hit in pfam_clans:
				clan = pfam_clans[hit]
			else: clan = 'clan'
			
			if spacing is 'env': begin, end = cols[19], cols[20] #env begin and end
			elif spacing is 'ali': begin, end = cols[17], cols[18] #uses alignment instead of env 
			i = int(cols[9])
			i_total = int(cols[10])
			
			
			if clan not in repeats:
				repeats[clan] = {}
			if hit not in repeats[clan]:
				repeats[clan][hit] = {} 
				repeats[clan][hit]['score'] = 0
			
			repeats[clan][hit][identifier+'_'+str(i)] = [begin, end] #save detected repeat to get fasta later
			
			#not necessary in initial parsing
			if state is 'iterative' and i > 1:
				#compare to previous domain for continuity
				#sometimes does not exist because of too low domain score or le
				#resolving overlaps
				if identifier+'_'+str(i-1) in repeats[clan][hit]:
					begin_prev, end_prev = repeats[clan][hit][identifier+'_'+str(i-1)]
					if int(end_prev) > int(begin): 
						avg_coord = (int(end_prev)+int(begin))/2
						repeats[clan][hit][identifier+'_'+str(i)] = [avg_coord, end]
						repeats[clan][hit][identifier+'_'+str(i-1)] = [begin_prev, (avg_coord-1)]
	
	
			#keep track of unique hits for filtering
			if hit not in pfam_hits: 
				pfam_hits[hit] = {}
				
			#species can occur multiple times... use protein_uri?	
			if protein_uri not in pfam_hits[hit]:
				pfam_hits[hit][protein_uri] = i_total #uses tblout number
				repeats[clan][hit]['score'] += sequence_score #add once for each species
			
			
	return(repeats, pfam_hits)

def parse_domtblout_stats(hmm_file, use_gacutoff = False):
	#returns per domain: clan, hmm length, query length, units per protein id with bitscore, model coordinates (to derive length/coverage), and ali sequence coordinates
	#only one domain per clan expected
	repeats = {}
	pfam_clans = get_pfam_clans(pfam_clans_file)
	pfam_gacutoff = get_pfam_gacutoff(pfam_gacutoff_file)
	
	with open(hmm_file,'r') as hmmr_tbl:
		for row in hmmr_tbl:
			if row[:1] == '#': continue; #skip header rows
			cols = row.split() #0 gene uri, 3 HMM match 'hit', 9 number, 12 i-e
			
			domain = cols[0] #name of domain, pfam accession is #1 but - after iteration
			protein_uri = cols[3]
			hmm_length = int(cols[2]) #tlen
			prot_length = int(cols[5]) #qlen
			seq_bitscore = float(cols[7])
			
			if use_gacutoff and cols[1] is not '-':
				if pfam_gacutoff[cols[1]] > seq_bitscore: continue
			unit_number = int(cols[9])
			unit_count = int(cols[10])
			dom_bitscore = float(cols[13])			
			model_coordinates = (cols[15], cols[16]) #alignment coordinates
			seq_coordinates = (cols[17], cols[18]) #alignment coordinates
			
			clan = pfam_clans[domain] if domain in pfam_clans else 'clan'
			
			#identifier = protein_uri+'_'+hit	#followed by lead number
			
			#ensembl_stable_id = filter(str.isalpha, protein_uri)
			#species = species_mapping[ensembl_stable_id]
					
			if domain not in repeats:
				repeats[domain] = {'clan':clan, 'length':hmm_length, 'orthologs_dict': {}} 
		
			if protein_uri not in repeats[domain]['orthologs_dict']:
				repeats[domain]['orthologs_dict'][protein_uri] = {'length':prot_length, 'unit_count':unit_count,'seq_bitscore':seq_bitscore,'units_dict':{}}
				
			repeats[domain]['orthologs_dict'][protein_uri]['units_dict'][unit_number] = {'dom_bitscore': dom_bitscore,'model_coordinates':model_coordinates, 'seq_coordinates': seq_coordinates}
			
	return(repeats)


def get_best_hit(repeats):
	#repeats[clan][hits]['score']
	
	best = (0, None) #score, hit values
	for clan in repeats:
		
		if len(repeats[clan]) == 0: continue
		
		for hit in repeats[clan]:
			
			if best[1] is None or repeats[clan][hit]['score'] > best[0]:
				
				best = (repeats[clan][hit]['score'] , repeats[clan][hit])
	
	return(best[1])

def write_fasta_file(outfile, best_hit_repeats, fasta, padding = 0):		
	best_hit_repeats.pop('score')
	with open(outfile, 'w') as outfile:	
		for identifier,value in best_hit_repeats.items():
			begin, end = value
			protein_uri = identifier.split('_')[0] #protein uri _ hit _ lead number
			sequence = fasta[protein_uri]
					
			#option to make selection larger
			begin = max((int(begin)-padding), 1)
			end = min((int(end)+padding), len(sequence)) #length -1?
					
			outfile.write('>'+identifier+'_'+str(begin)+'-'+str(end)+'\n')
			outfile.write(sequence[begin-1:end]+'\n')
			


def write_masked_file(outfile, best_hit_repeats, fasta):
	#replace pfam domain with XXX'es in protein sequences
			
	best_hit_repeats.pop('score')
	with open(outfile, 'w') as outfile:	
		for protein_uri in fasta:
			sequence = fasta[protein_uri]
			outfile.write('>'+protein_uri+'\n')
			protein_uri_hits = [best_hit_repeats[identifier] for identifier in best_hit_repeats if protein_uri in identifier]
			for value in protein_uri_hits: #protein uri _ hit _ lead number
				begin = int(value[0]); end = int(value[1])
				#print(len(sequence[begin-1:end]), len('x'*(end-begin+1)))	#equal length
				sequence = sequence.replace(sequence[begin-1:end], 'x'*(end-begin+1))
				
			outfile.write(sequence+'\n')
			
def write_chopped_file(outfile, best_hit_repeats, fasta):
	minimum_length_meme = 30 
	species_mapping = get_species_mapping(species_mapping_file)	
			
	#chop up protein sequences if there is a pfam domain 		
	best_hit_repeats.pop('score')
	buffer=''
	human = False
	mouse = False
	third_prot = False
	
	for protein_uri in fasta:
		sequence = fasta[protein_uri]
		species = species_mapping[filter(str.isalpha, str(protein_uri))]
		
		protein_uri_hits = [best_hit_repeats[identifier] for identifier in best_hit_repeats if protein_uri in identifier]
		#print(protein_uri,protein_uri_hits, len(sequence))
		if len(protein_uri_hits) == 0: 
			#write full length sequence
			buffer += '>'+protein_uri+'\n'+sequence+'\n'
			if species == 'HSAP': human = True
			elif species == 'MMUS': mouse = True
			else: third_prot = True
			
		else:
			i=0
			offset=0
			
			for value in protein_uri_hits: # identifier = protein uri _domain_ lead number
				i+=1
				begin = int(value[0])-1; end = int(value[1])-1 #1-based coordinates from HMM, python needs zero-based
					
				#minimum length of sequence necessary for MEME?
				if (begin-offset) > minimum_length_meme:
					#write n-terminal sequence	
					buffer += '>'+protein_uri+'_'+str(i)+'\n'
					buffer += sequence[offset:begin]+'\n' #up to but not including begin
					
					if species == 'HSAP': human = True
					elif species == 'MMUS': mouse = True
					else: third_prot = True
					
				if i == len(protein_uri_hits) and (len(sequence)-(end+1)) > minimum_length_meme:
					if buffer == '':
						buffer += '>'+protein_uri+'\n'
					else: 
						buffer += '>'+protein_uri+'_'+str(i+1)+'\n'		
					buffer += sequence[(end+1):]+'\n'
					
					if species == 'HSAP': human = True
					elif species == 'MMUS': mouse = True
					else: third_prot = True
							
				offset = end+1 #new begin, including slice so need to offset end+1
				
	print(human,mouse,third_prot)			
	if len(buffer) > 0 and human and mouse and third_prot:
		with open(outfile, 'w') as outfile:	
			outfile.write(buffer)
						


def parse_pfam_to_columns(pfam_summary):
	pfam_stats = {} #to be made into pandas dataframe
	species_dist, unit_count_dist, model_coverage_dist = {},{},{}
	species_mapping = get_species_mapping(species_mapping_file)
	pfam_clans = get_pfam_clans(pfam_clans_file)		
	for og_id in pfam_summary:
		#expect only one domain
		domain_name = [x for x in pfam_summary[og_id].keys()][0]
		domain = pfam_summary[og_id][domain_name]
		
		for species in species_mapping.values():
			species_dist[species] = np.nan	#orthologs per species
			unit_count_dist[species] = [np.nan]	#averaged over species, unit count
			model_coverage_dist[species] = [np.nan]	#averaged over species, coverage of units		
			
		for protein_uri in domain['orthologs_dict']:
			species = species_mapping[filter(str.isalpha, str(protein_uri))]
			if species_dist[species] is np.nan: species_dist[species] = 1
			else: species_dist[species]+=1
						
			protein_domain = domain['orthologs_dict'][protein_uri]
			#protein_domain['seq_bitscore']
			 
			unit_count_dist[species].append(protein_domain['unit_count'])
					
			model_coverage_list = []
			for unit_number,unit_data in protein_domain['units_dict'].items():
				#unit_data['dom_bitscore']
				#unit_data['model_coordinates']
				#unit_data['seq_coordinates']
				
				#use model coordinates and length to calculate percentage covered, model coordinates max (1, length)
				model_coverage_list.append( (float(unit_data['model_coordinates'][1])-float(unit_data['model_coordinates'][0])+1)/float(domain['length']))
			model_coverage = np.mean(model_coverage_list)
			model_coverage_dist[species].append(model_coverage)
		
		if 'clan' not in domain: 
			domain['clan'] = '-' if domain_name not in pfam_clans else pfam_clans[domain_name]
			
		pfam_stats[og_id] = { 'length':domain['length'], 'clan':domain['clan'], 'domain':domain_name }
		
		for species, number in species_dist.items():
			pfam_stats[og_id]['orthologs_'+species] = number
		for species, number_list in unit_count_dist.items():
			pfam_stats[og_id]['unit_count_'+species] = np.nanmean(number_list) if number_list != [np.nan] else 0
		for species, number_list in model_coverage_dist.items():
			pfam_stats[og_id]['model_coverage_'+species] = np.nanmean(number_list) if number_list != [np.nan] else 0
		
		pfam_stats[og_id]['avg_unit_count'] = np.nanmean([np.nanmean(unit_count_dist[sp]) for sp in unit_count_dist if unit_count_dist[sp] != [np.nan]])
		pfam_stats[og_id]['avg_model_coverage'] = np.nanmean([np.nanmean(model_coverage_dist[sp]) for sp in model_coverage_dist if unit_count_dist[sp] != [np.nan]])
		pfam_stats[og_id]['total_orthologs'] = sum([x for x in species_dist.values() if x is not np.nan])	

	#pfam_df = pd.DataFrame.from_dict(pfam_stats, orient='index')
	return(pfam_stats)
