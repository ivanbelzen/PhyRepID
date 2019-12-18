# PhyRepID pipeline documentation
[https://docs.google.com/document/d/1SfTu8sbzZAmX9UJLnoqzd-pvcVvzN4v4PypTHRK2vnY/edit#](https://docs.google.com/document/d/1SfTu8sbzZAmX9UJLnoqzd-pvcVvzN4v4PypTHRK2vnY/edit#)

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
** ** 

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
- pfam/aligned/{og}\_{hit}.linsi.fa.treefile
- genetrees\_nhx/{og\_id}.nhx

Output:
- pfam/aligned/{og}\_{hit}.linsi.fa.treefile.rooted
- pfam/treefix/{og}.stree 
- pfam/treefi/{og}.smap
 
**run_treefix** and **mv_treefix**
Input:
- pfam/aligned/{og}\_{hit}.linsi.fa.treefile.rooted
- pfam/treefix/{og}.stree 
-  pfam/treefi/{og}.smap
   - pfam/aligned/{og}\_{hit}.linsi.fa.iqtree
  - model, extract best-fit model according to BIC
  - *(also uses alignment but not specified as input)*

Output:
- pfam/treefix/{og}\_{hit}.treefix.tree
  - copied from pfam/aligned/{og}\_{hit}.treefix.tree
- ..
   

*treefix -s {input.stree} -S {input.smap} -A .linsi.fa -o
.linsi.fa.treefile.rooted -n .treefix.tree -V 0 -m
treefix.models.iqtreemodel.CoarseModel -e \\'-t AA -m \\'\$model\\'\\'
{input.treefile}*

**run\_treefix\_annotate**
- input:
    -  pfam/treefix/{og}\_{hit}.treefix.tree
    -  pfam/treefix/{og}.stree pfam/treefi/{og}.smap
- output
    -  pfam/treefix/{og}\_{hit}.treefix.mpr.recon
    -  ...
    
*tree-annotate -s {input.stree} -S {input.smap} {input.tree}*



****


**parse\_treefix.bash**

tryin gto make notung work

first with example

cdc20 ENSGT00870000136444\_ENSG00000117399

had to get rid of the dubious annotation in gene trees

( NHX duplication = dubious) to make notung work

parsing the parsable.txt and reconciled to get root dups

don’t include the empty internal nodes, but do include all tips/proteins

also include events per species

**notung\_output.json**

pfam\_results\_final.json

made by parse\_tblout\_final\_stats.py

I want to add length of hmm models, and number of orthologs per species

analyse pfam .py

new dict in

og\_id (genetree+gene\_id) = { domain, clan, hmm length, orthologs per
species \[\], units per species \[\] }

made separate python script for generation of this json file

**analyse\_pfam\_tblout.py** generate pfam\_summary.json and hists?

analyse\_meme

meme\_summary.json 2677

dict with og\_id as key, consensus and repeats dict from
parse\_tblout\_stats

This repeats dict contains:

-   domain = { clan, length, orthologs\_dict = {

    -   protein\_uri : {unit\_count, seq bitscore, units\_dict = {

        -   unit\_number: { dom bitscore, model coordinates, seq
            > coordinates }

        -   ... } }

    -   ...}

-   ...

analyse sequence coverage .py



## Meme dataset construction:

Snakefork\_chop

input todo .json which is a subset of meme\_dataset.json that had not be
run/completed before

Snakefork\_xed

run meme on the fasta\_ogs\_chopped

which are the ones from todo that have &gt;=3 seq of which one mouse and
human

**analysis pfam dataset**

generate\_pfam\_dataframe.py

makes analysis/pfam/pfam\_dataframe.csv

**Figures:**

Projection on species tree:

analyse pfam adjusted to map to species tree

repeat - integration is done automatically because you map multiple
genes to the same species if they have paralogs

('Unknown vs total', {'d': 7742, 'l': 19183}, 30042, 107535)

**Dictionaries**

-   pfam\_summary.json

    -   key: og\_id\_domain

    -   repeats dict containing: domain = { length,

        -   orthologs\_dict = { protein\_uri : {unit\_count, seq
            > bitscore,

            -   units\_dict = { unit\_number: { dom bitscore, model
                > coordinates, seq coordinates }, ... } }, ...} }

-   notung\_ete\_summary.json

    -   contains seperate dicts for reconciled repeat / gene tree\
        > identifier respectively og\_id\_domain and og\_id

    -   identifier = {nD, nL, root\_dup,

        -   events = { node\_name: {d,l} , … },

        -   species = { species\_abbrev : {d,l } ... }, }

-   pfam\_evo\_events\_summary.json

    -   contains seperate dicts for reconciled repeat / gene tree\
        > identifier respectively og\_id\_domain and og\_id


