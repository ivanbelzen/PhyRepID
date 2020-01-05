# PhyRepID pipeline documentation


##   Data collection (Snakefile)

**retrieve\_orthologs(.py)**

Retrieve orthologs of human protein coding genes from selected species using a SPARQL query on the ENSEMBL endpont. 

Input: 

 - ensembl\_stable\_id\_species.json

Json dictionary with list of selected species. Default 14 species if possible in pairs of closely related species: human-mouse, opossum-tasmanian devil, duck-chicken, stickleback-takifugu, spotted gar, zebrafish, xenopus, platypus, anole lizard, turtle.

Output: 

 - orthologs.json  

Json dictionary with human\_gene\_id as key and orthologs gene ids as values


**filter\_orthologs(.py)**

Filter criteria for the OG which is formed by a human protein coding gene that:
- should have an ortholog in mouse
-   should have orthologs for all species in group or none (pairs defined in ensembl\_stable\_id\_species.json) “excluded pairs”
-   minimum of 10 orthologs (removed if &lt;10)
-   maximum of 3 orthologs in a species (removed if &gt;3 “excluded orth rel”)
Input: 
 - orthologs.json

Output:
-   orthologs\_filtered.json
-   excluded\_pairs.json
-   excluded\_orth\_rel.json

**retrieve\_genetrees(.py)**

Queries ensembl API for genetrees of the filtered orthologs. Downloading and saving all raw data for reproducibility.

Input: 
- orthologs\_filtered.json

Output:
-  ensembl\_api/{gene\_id}.json
-   genes\_to\_genetrees.json *(mapping)*
-   genetrees\_ogs.json *(mapping)*
-   redo\_orth\_parse.json

Completed if: ensembl\_api/{gene\_id}.json file is created for all gene\_ids in orthologs\_filtered.json

## Repeat detection (Snakefile)
Files  to check if next steps are done for all proteins are generated during data collection.

**parse\_genetree(.py)** 

Subtract smallest gene tree  containing all orthologs; and acquire fasta sequences for proteins in gene tree. 
For each *og\_id* in genes\_to\_genetrees.json where *{og\_id} = {genetree\_id}\_{gene\_id}*.

Input:
-    ensembl\_api/{gene\_id}.json
-  orthologs\_filtered.json

Output:
-   genetrees\_nhx/{genetree\_id}\_{gene\_id}.nhx
-    fasta\_ogs/{genetree\_id}\_{gene\_id}.fa

**hmmpress, hmmscan**
Detect Pfam domains using HMMER3 using liberal thresholds: sequence bit score of 12.5 and domain score 0. 
For each *og\_id* in genes\_to\_genetrees.json

Input:
-  fasta\_ogs/{og\_id}.fa
-  Pfam-A.hmm, Pfam-A.hmm.h3i 

Output:
-   pfam/hmm/{og\_id}.tblout

**repeat_detection_iteration**
Iterative pfam domain annotation. Improvement on regular hmmscan for more accurate repeat detection.

run\_pfam\_iteration\_single.bash *\$1: pfam/hmm/{og\_id}.tblout*

Completed if: dummy file is made for all {og\_id} from
genes\_to\_genetrees.json


Scripts:
-   parse\_tblout\_init.py *\$1: pfam/hmm/{og\_id}.tblout \$2: fasta\_ogs/{og\_id}.fa*

    -   Use sequences gathering cutoff, select best Pfam model for each clan based on sequence bitscore maximalisation (over-all).

    -   Filter Pfam hits based on minimum \#repeats (&gt;=3) in human, presence in mouse.
    -   Writes a fasta file for each best hit of clan with all of the  repeats
    -   hmm\_results.json before any filtering
-   parse\_tblout\_iterate.py *\$1: pfam/tblout{og\_hit\_id}.tblout \$2: fasta\_ogs/{og\_id}.fa*
    -   extract repeats from tblout, same as happens in parse\_tblout\_init
-   parse\_tblout\_final.py *\$1: pfam/tblout{og\_hit\_id}.tblout \$2: fasta\_ogs/{og\_id}.fa*
    -   Use more liberal slicing for making tree, envelope spacing and padding >5
    -   hmm\_results\_final.json

Input:
-   pfam/hmm/{og\_id}.tblout
-   fasta\_ogs/{og\_id}.fa
-   Pfam clans and gathering cutoff:  Pfam-A.clans.tsv ga\_cutoff.tsv

Output:
-   parse\_tblout\_init.py:
    -   pfam/repeats/{og\_id\_hit}.fa *(serves as input for iteration)*
    -  hmm\_results.json
    -  excluded\_repeats.json
-   pfam/aligned/{og\_id\_hit}.linsi.fa *(with last iteration a copy to pfam/alignments/{og\_id\_hit}.linsi.fa)*
-   pfam/tblout/{og\_id\_hit}.tblout
-   parse\_tblout\_final.py
    - hmm\_results\_final.json 
    - progress\_files/{og\_id}.pfamhmm.done


Alignment with mafft linsi: *mafft --quiet --localpair --maxiterate 1000*

Build profile HMMs based on alignment: *hmmbuild --fast --symfrac 0.6 -n*

Check for length of profile, prevent profile drift: enforce minimum length of 20

Iterates max 20 times or until convergence 
.fa.old as temporary file in /pfam/repeats/


**run_iqtree** and **mv_iqtree**

Input:

-   pfam/aligned/{og\_id\_hit}.linsi.fa

Output:

- pfam/aligned/{og\_id\_hit}.linsi.fa.treefile 
- pfam/aligned/{og\_id\_hit}.iqtree 
- ....
-   pfam/trees/{og\_id\_hit}.treefile

*iqtree-omp -s {input} -bb 1000 -nt AUTO -mset LG,WAG,VT,Dayhoff,JTT
-redo*

Completed if: .treefile is made for all {og\_id\_hit} from
hmm\_results\_final.json

(copy to pfam/trees is separate rule)

**reparse_genetrees**
parse_genetree.py

Input:
-    ensembl\_api/{gene\_id}.json

Output:
-   genetrees\_nhx/{genetree\_id}\_{gene\_id}.nhx
-    fasta\_ogs/{genetree\_id}\_{gene\_id}.fa


**prepare_treefix**
Script: generate_treefix_input.py {input.treefile} {input.genetree}'

Input:
- pfam/aligned/{og\_id\_hit}.linsi.fa.treefile
- genetrees\_nhx/{og\_id}.nhx

Output:
- pfam/aligned/{og\_id\_hit}.linsi.fa.treefile.rooted
- pfam/treefix/{og\_id}.stree 
- pfam/treefi/{og\_id}.smap
 
**run_treefix** and **mv_treefix**
Input:
- pfam/aligned/{og\_id\_hit}.linsi.fa.treefile.rooted
- pfam/treefix/{og\_id}.stree 
-  pfam/treefi/{og\_id}.smap
   - pfam/aligned/{og\_id\_hit}.linsi.fa.iqtree
  - model, extract best-fit model according to BIC
  - *(also uses alignment but not specified as input)*

Output:
- pfam/treefix/{og\_id\_hit}.treefix.tree
  - copied from pfam/aligned/{og\_id\_hit}.treefix.tree
- ..
   

*treefix -s {input.stree} -S {input.smap} -A .linsi.fa -o
.linsi.fa.treefile.rooted -n .treefix.tree -V 0 -m
treefix.models.iqtreemodel.CoarseModel -e \\'-t AA -m \\'\$model\\'\\'
{input.treefile}*

**run\_treefix\_annotate**
- input:
    -  pfam/treefix/{og\_id\_hit}.treefix.tree
    -  pfam/treefix/{og\_id}.stree pfam/treefi/{og\_id}.smap
- output
    -  pfam/treefix/{og\_id\_hit}.treefix.mpr.recon
    -   pfam/treefix/{og\_id\_hit}.treefix.nhx.tree
    -  ....txt
    
*tree-annotate -s {input.stree} -S {input.smap} {input.tree}*


## Meme dataset construction:
*DRAFT*

Snakefork\_chop

input todo .json which is a subset of meme\_dataset.json that had not be
run/completed before

Snakefork\_xed

run meme on the fasta\_ogs\_chopped

which are the ones from todo that have &gt;=3 seq of which one mouse and
human


## Post-processing
**parse_evo_events(.py)**
Parse TreeFix output by removing inconsistent duplications, inferring events and summarizing duplications/losses of repeat tree  reconciled with gene tree. 
Netto duplications are non-ancestral, taking into account the genetree root instead of vertebrate common ancestor

Count duplications/losses using tree reconciliation of repeat tree with gene tree
For the gene duplications, count dup/loss in gene tree compared to species tree
 
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
for Schaper comparison and human-lineage overview 

- Input:
  - pfam/meme_repeat_stats.json
- Output: 
  - og_human_mapping.json

**generate_phyrepid_results(.py)**
 Makes table with all columns: ExAC, Selectome, Schaper, gene tree annotation.

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
