# PhyRepID pipeline documentation

 
## Overview of the workflow

The workflow is defined in five Snakefiles which link together tools, algorithms, Python and Bash scripts. In addition, comparisons with other datasets are done with stand-alone Python and R scripts [(documentation of analysis scripts)](docs/analysis_documentation.md).

 1. Data collection
 2. Domain detection (Pfam) & optimization
 3. Motif detection (MEME) & optimization
 4. Phylogenetic analysis
 5. Post-processing

## Setup and running

 - Define root in *pipeline_methods.py* and *config.sh*, if desired, (path) variables can also be changed. 
 - Install the required tools and dependencies (see [README.md](README.md))
 - Run the Snakefiles in the order defined above.

##   Data collection 

**setup_pipeline**

Makes the required folder structure and downloads the Pfam-A hidden markov models (.hmm files).
Note: fist do the manual steps of installing dependencies (see [README.md](README.md)) and defining the root directory.

**retrieve\_orthologs(.py)**

Retrieve orthologs of human protein coding genes from selected species using a SPARQL query on the ENSEMBL endpont. 

- Input: 
   - ensembl\_stable\_id\_species.json
- Output: 
  - logs/orthologs.json  
     - mapping of human\_gene\_id to orthologs gene ids

Default 14 species, if possible, in pairs of closely related species: human-mouse, opossum-tasmanian devil, duck-chicken, stickleback-takifugu, spotted gar, zebrafish, xenopus, platypus, anole lizard, turtle.

**retrieve_orthologs_filter(.py)**

Filter criteria for the OG that is formed by a human protein coding gene that:
- should have an ortholog in mouse
-   should either have orthologs in all members of a group of closely related species or none of them (pairs of closely related species defined in ensembl\_stable\_id\_species.json) (removed are in “excluded pairs”)
-   maximum of 3 orthologs in a species (removed are in “excluded orth rel”)


- Input: 
   - logs/orthologs.json
 - Output:
   - logs/orthologs\_filtered.json
   -   logs/excluded\_pairs.json
   -   logs/excluded\_orth\_rel.json

**retrieve\_genetrees(.py)**

Queries ensembl API for genetrees of the filtered orthologs. Downloading and saving all raw data for reproducibility.

- Input: 
  - logs/orthologs\_filtered.json
- Output:
   -  ensembl\_api/{gene\_id}.json
  -  logs/genes\_to\_genetrees.json 
  -  logs/genetrees\_ogs.json 
  -  logs/redo\_orth\_parse.json

Completed if: ensembl\_api/{gene\_id}.json file is created for all gene\_ids in *orthologs\_filtered.json*.

## Domain detection (Pfam)
Check if next steps are done for all protein families that are retrieved  during data collection by using *genes\_to\_genetrees.json* as mapping. For each *og\_id* in genes\_to\_genetrees.json where *{og\_id} = {genetree\_id}\_{gene\_id}*.

**parse\_genetree(.py)** 

Extract smallest gene tree containing all orthologs by parsing the unfiltered, large gene tree retrieved in the previous step. Also acquire fasta sequences of the orthologs. 

- Input:
  -    ensembl\_api/{gene\_id}.json
  -  logs/orthologs\_filtered.json
- Output:
  -   genetrees\_nhx/{genetree\_id}\_{gene\_id}.nhx
  -    fasta\_ogs/{genetree\_id}\_{gene\_id}.fa

**hmmpress, hmmscan**

Detect Pfam domains using HMMER3 using liberal thresholds: hmmscan using a sequence bit score of 12.5 and domain score 0. 

- Input:
  -  fasta\_ogs/{og\_id}.fa
  -  resources/Pfam-A.hmm, Pfam-A.hmm.h3i 
- Output:
   -   pfam/hmm/{og\_id}.tblout

**repeat_detection_iteration**

Improvement on default HMMscan for more accurate repeat detection by making OG-specific hidden markov models (HMMs) for a repeat unit by using an iterative process. 

run\_pfam\_iteration.sh *\$1: pfam/hmm/{og\_id}.tblout*

Scripts:
-   parse\_tblout\_init.py *\$1: pfam/hmm/{og\_id}.tblout \$2: fasta\_ogs/{og\_id}.fa*

    -   Use sequences gathering cutoff, select best Pfam model for each clan based on sequence bitscore maximalisation (over-all).

    -   Filter Pfam hits based on minimum \#repeats (&gt;=3) in human and presence in mouse.
    -   Writes a fasta file for each best hit of clan with all of the  repeats
    -   logs/hmm\_results.json - log before any filtering
-   parse\_tblout\_iterate.py *\$1: pfam/tblout{og\_hit\_id}.tblout \$2: fasta\_ogs/{og\_id}.fa*
    -   extract repeats from tblout, same as happens in parse\_tblout\_init
-   parse\_tblout\_final.py *\$1: pfam/tblout{og\_hit\_id}.tblout \$2: fasta\_ogs/{og\_id}.fa*
    -   Use more liberal slicing for making tree, envelope spacing and padding >5
    -   logs/hmm\_results\_final.json - log after filtering

Input:
-   pfam/hmm/{og\_id}.tblout
-   fasta\_ogs/{og\_id}.fa
-   Pfam clans and gathering cutoff:  resources/Pfam-A.clans.tsv, resources/ga\_cutoff.tsv

Output:
-   parse\_tblout\_init.py:
    -   pfam/repeats/{og\_id\_hit}.fa *(serves as input for iteration)*
    -  logs/hmm\_results.json
    -  logs/pfam_excluded_ogs.json'
-   pfam/aligned/{og\_id\_hit}.linsi.fa *(with last iteration a copy to pfam/alignments/{og\_id\_hit}.linsi.fa)*
-   pfam/tblout/{og\_id\_hit}.tblout
-   parse\_tblout\_final.py
    - logs/hmm\_results\_final.json 
    - progress\_files/{og\_id}.pfamhmm.done *empty file to track progress*

Alignment with MAFFT linsi: *mafft --quiet --localpair --maxiterate 1000*

Build profile HMMs based on alignment: *hmmbuild --fast --symfrac 0.6 -n*

Check for length of profile, prevent profile drift: enforce minimum length of 20

Scan with these profiles: *hmmscan --domT 0 --domtblout *
Iterates max 20 times or until convergence, comparing the new output sequences to previous round.

## Motif detection (MEME)

Run motif detection on all OGs minus those with Pfam repeats  (>= units 3 in human, presence in mouse).
		
**mask_sequences:**

Prepare for motif detection by masking Pfam sequences. Select best matching Pfam domain and chop it from the protein sequence from the FASTA files to prevent MEME detecting known domains. 

- Input:
   - pfam/hmm/{og_id}.tblout 
   - fasta_ogs/{og_id}.fa
- Output:
    - fasta_ogs_chopped/{og\_id}.fa
    - progress_files/{og_id}.chopped_fasta.done 
	
**run_meme_detection:**

Run MEME to detect denovo repeats.

- Input:
    - fasta_ogs_chopped/{og_id}.fa
- Output:
    - denovo/meme/{og_id}/meme.html
  
*meme {input} -oc {output} -protein -mod anr -minsites 4 -minw 12 -maxw 30 -time 300 -maxiter 10 -nostatus -maxsize 10000000*

**parse_meme:**

Parse meme output and derive repeated motifs

- Input:
    - denovo/meme/{og_id}/meme.html
- Output:
    - denovo/meme/repeats/{og_id}.fa
    - progress_files/{og_id}.meme_parse_chop.done
    - logs/meme_results_parse_log.json

**select_repeats:**

As a first step to repeat motif detection, build HMM profile of detected motif by MEME and scan sequence.
Note: only one motif-repeat is currently analysed per orthologuous group.

- Input:
    - denovo/meme/repeats/{og_id}.fa
- Output:
    - denovo/aligned/{og_id}.linsi.fa
    - denovo/hmm/{og_id}.hmm
    - denovo/tblout/{og_id_hit}.tblout 
    - progress_files/{og_id}.meme_filter_chop.done
    - parse_denovo_tblout_init.py 
        - denovo/repeats/{og_id}.linsi.fa
        - logs/meme_hmm_results.json - log before any filtering
        - logs/meme_excluded_ogs.json
      
Analoguous to **repeat_detection_iteration**  for Pfam domains:  build profile HMMs based on alignment with MAFFT linsi. Use HMMscan to scan sequences in OG with this model and generate first set of repeat sequences to start iteration.  

**repeat_detection_iteration**

Similar to *run_pfam_iteration.sh* using scripts *parse_denovo_tblout_iterate.py* *parse_denovo_tblout_final.py*. However, repeats are already filtered and parsed by *parse_denovo_tblout_init.py* in the previous step.

- Input
    - denovo/repeats/{og_id}.fa	
- Output:
   -    denovo/aligned/{og_id}.linsi.fa *(with last iteration a copy to denovo/alignments/{og_id}.linsi.fa)*
   -   denovo/tblout/{og_id_hit}.tblout
    - logs/meme_hmm_results\_final.json  - log after filtering
    - progress_files/{og_id}.meme_hmm.done *empty file to track progress*
	

## Phylogenetic analysis

Finish both Pfam and MEME repeat detection first.
This snakefile determines the OG ids and OG-hit ids to run pipeline on all files in */pfam/aligned/{og_id_hit}.linsi.fa* and  */denovo/aligned/{og_id}.linsi.fa* 

Hence the *"pfam"* can be replaced by *"denovo"* in the paths outlined below unless otherwise indicated.

**run_iqtree** and **mv_iqtree**
Generate phylogenetic tree from repeats using IQtree. Note that IQTree generates more files in the *pfam/aligned/* directory we do not use.

- Input:
   - pfam/aligned/{og\_id\_hit}.linsi.fa
- Output:
   - pfam/aligned/{og\_id\_hit}.linsi.fa.treefile 
	   - copied to pfam/trees/{og\_id\_hit}.treefile
   - pfam/aligned/{og\_id\_hit}.iqtree 
  
*iqtree-omp -s {input} -bb 1000 -nt AUTO -mset LG,WAG,VT,Dayhoff,JTT
-redo*


**generate_treefix_input(.py)**

- Input:
  - pfam/aligned/{og\_id\_hit}.linsi.fa.treefile
  - genetrees\_nhx/{og\_id}.nhx
- Output:
   - pfam/aligned/{og\_id\_hit}.linsi.fa.treefile.rooted
   - pfam/treefix/{og\_id}.stree 
   - pfam/treefix/{og\_id}.smap
 
**run_treefix** and **mv_treefix**

Note that Treefix generates more files in the *pfam/aligned/* directory we do not use.

- Input:
  - pfam/aligned/{og\_id\_hit}.linsi.fa.treefile.rooted
  - pfam/treefix/{og\_id}.stree 
  - pfam/treefix/{og\_id}.smap
  - pfam/aligned/{og\_id\_hit}.linsi.fa.iqtree
     - to derive best-fit model according to BIC
  - *(also uses alignment but not specified as input)*
- Output:
  - pfam/treefix/{og\_id\_hit}.treefix.tree
     - copied from pfam/aligned/{og\_id\_hit}.treefix.tree   

*treefix -s {input.stree} -S {input.smap} -A .linsi.fa -o .linsi.fa.treefile.rooted -n .treefix.tree -V 0 -m treefix.models.iqtreemodel.CoarseModel -e \\'-t AA -m  \\'$model\\'\\' {input.treefile}*
(Uses the best-fit model determined before. )

**run\_treefix\_annotate**

Note that Treefix Annotate generates more files in the *pfam/treefix/* directory we do not use.

- Input:
    -  pfam/treefix/{og\_id\_hit}.treefix.tree
    -  pfam/treefix/{og\_id}.stree 
    - pfam/treefix/{og\_id}.smap
- Output
    -   pfam/treefix/{og\_id\_hit}.treefix.nhx.tree
      
*tree-annotate -s {input.stree} -S {input.smap} {input.tree}*


## Post-processing

**parse_evo_events(.py)**

Parse TreeFix output by removing inconsistent duplications, inferring events and summarizing duplications/losses of repeat trees  reconciled with gene trees. 
A protein repeat duplication (PRD) score is calculated for each OG as relative measure of repeat evolution. 
This PRD score can be used to rank protein families, compare them based on repeat evolution and find OGs with rapidly-evolving repeats. 
 
<p align="center"> 
<b>PRD score = (x - u) / n;</b><br />
 with x = netto duplications since most recent common ancestor (gene tree root); u = mean duplications in the full dataset, n = number of proteins in the OG.</p>

Count duplications/losses using tree reconciliation of repeat trees with gene trees. Note: the netto duplications are non-ancestral duplications, taking into account the genetree root instead of vertebrate common ancestor

 For the gene duplications, count duplications/losses  in the gene tree compared to species tree
 
- Input (glob):
    -  pfam/treefix/{og\_id\_hit}.treefix.nhx.tree
       - denovo/treefix/{og\_id\_hit}.treefix.nhx.tree
    - genetrees\_nhx/{og\_id}.nhx
  - ensembl species tree 
- Output: 
  - pfam/treefix/adjusted/{og\_id\_hit}.ete
     - denovo/treefix/adjusted/{og\_id\_hit}.ete
  - pfam/treefix/adjusted/images/{og\_id\_hit}.pdf
    - denovo/treefix/adjusted/images/{og\_id\_hit}.pdf
  - dict with what nodes are rearranged:
   log_pfam_rearranged_nodes.json / log_meme_rearranged_nodes.json
  - output dict dup/loss repeat tree vs gene tree: pfam_evo_events.json / meme_evo_events.json
  - output dict dup/loss gene tree vs species tree: pfam_evo_events_genetrees.json / meme_evo_events_genetrees.json

 
*Duplication consistency correction:*
Remove spurious duplications from repeat tree that have a consistency score of 0 ( adjustable threshold)
Duplication consistency depends on number of orthologs (count from gene tree). 
Prune/graft gene tree on repeat-tree parts with inconsistent duplications 

*pfam/meme_evo_events.json*
repeat { og_id_hit : {netto_dup, loss, events = {dup = {node:cnt, ...} , loss = {node:cnt, ...} } } }, 
*pfam/meme_evo_events_genetrees.json*:
genetree {  og_id_hit : {root, netto_dup, loss, events = {dup = {node:cnt, ...} , loss = {node:cnt, ...} } } }

**analyse_repeat_stats(.py)**

Generate repeat stats from tblout for both Pfam and MEME after repeat detection/optimization. Parse only finished OGs that have Treefix annotate output files.

Initial: Pfam summary before repeat detection optimization. 
Initial filtered: match up OG ids from initial to the OG-domain ids from later to compare improvement in repeat detection

- Input (glob):
    - pfam/treefix/{og\_id\_hit}.treefix.mpr.recon /denovo/treefix/{og\_id\_hit}.treefix.mpr.recon 
    - pfam/tblout/{og\_id\_hit}.tblout / denovo/tblout/{og\_id\_hit}.tblout
    - (glob) pfam/hmm/*.tblout
- Output:
    - pfam/meme_repeat_stats.json
   -   pfam_repeat_stats_initial.json
   - pfam_repeat_stats_initial_filtered.json

Output optional (not by default): report OGs that have repeats in less than 4 proteins - possibly exclude later because of less trustworthy trees: exclude_ogs_lt4_pfam.lst  / exclude_ogs_lt4_meme.lst

**_repeat_stats.json*
{ og\_id\_domain {domain = { clan, length, orthologs\_dict = { protein\_uri : {unit\_count, seq bitscore, units\_dict = { unit\_number: { dom bitscore, model coordinates, seq coordinates }, ... } }, ...} }

**generate_human_og_mapping.py**

Generate OG to human gene mapping
for the human-lineage overview and comparison to external datasets.

- Input:
  - pfam/meme_repeat_stats.json
- Output: 
  - og_human_mapping.json

**generate_phyrepid_results(.py)**

Generates PRD score ranking and annotate with comparisons to external datasets if available [(documentation)](docs/analysis_documentation.md).
Makes a table with all columns: ExAC, Selectome, Schaper, gene tree annotation.

- Input
  - pfam/meme_evo_events.json
  - pfam/meme_evo_events_genetrees.json  
  - pfam/meme_repeat_stats.json
- Optional input:
  -  resources/selectome.tsv (for Selectome) 
     - og_human_mapping.json (from generate_human_og_mapping.py)
    -  resources/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt (for ExAC) 
       - resources/human_gene_symbol.tsv
  - schaper_comparison.json (from analyse_schaper_comparison.py)
- Output
   - phyrepid_results_simple.tsv
   - phyrepid_results_full.tsv
  
  
**generate_phyrepid_results_human**

Analyses only evens in human lineage (full lineage) or branch (human only) and outputs table with duplications on those branches only.

- Input: 
  - og_human_mapping.json (from generate_human_og_mapping.py)
  - pfam/meme_evo_events.json
 - Output
   -    phyrepid_results_human_full_lineage.tsv
   - phyrepid_results_human_only.tsv
