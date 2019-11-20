# Pipeline documentation

## Introduction

The protein repeat evolution (PRE/PhyRepD) pipeline quantifies repeat evolution using comparative (phylo)genomics of orthologous groups (OGs). 

The pipeline consists of ~~three~~ four components: 

 1. **Repeat detection:** Detection of a broad spectrum of protein repeats consisting of either structural domains or short linear motif sequences.
 2. **Repeat optimization:** Improving the sensitivity and precision of detection by making OG-specific repeat models.
 3. **Tree reconciliation:** Inferring evolutionary events in the repeat region through phylogenetic comparison of repeat trees to gene trees.
 4. **Post-processing / duplication analysis:** Finally, to quantify repeat evolution and compare protein families, a duplication score was derived from the post-ancestral duplications and the number of proteins in the OG.

## Overview of the workflow

The workflow for components 1-3 is defined in two Snakemake files which link together Python and Bash scripts. Post-processing and analysis of the inferred repeat duplications, as well as comparisons with other datasets are done with various Python and R scripts. 

 - Data collection (Snakemake)
 - Snakemake: repeat detection
 - Scripts collection: repeat duplication analysis 


## Prerequisites

### Python
Version 2.7 with pip (9.0.1)
 *pipeline\_methods.py*
contains configuration like paths and environmental variables, as well as methods used by other Python scripts

**Data retrieval and processing**
SPARQLWrapper (1.8.1)
requests (2.18.4)
pandas (0.23.1)
numpy (1.13.1)

**Tree reconciliation** 
anaconda-client (1.2.2)
conda (4.3.27)
treefix (1.1.10)
ete3 (3.0.0b36)

**Analysis:**
matplotlib (2.2.2)
matplotlib-venn (0.11.5)
seaborn (0.9.0)
scipy (1.1.0)
scikit-learn (0.19.1)
goatools (0.8.4)

### Software
HMMer
TreeFix


##   Data collection (Snakefile)

### Retrieve orthologs
Retrieve orthologs of human protein coding genes from selected species using a SPARQL query on the ENSEMBL endpont. 

Script: retrieve\_orthologs.py

Input: ensembl\_stable\_id\_species.json with list of selected species. Default 14 species if possible in pairs of closely related species: human-mouse, opossum-tasmanian devil, duck-chicken, stickleback-takifugu, spotted gar, zebrafish, xenopus, platypus, anole lizard, turtle (14 species)

Output: orthologs.json  dictionary with human\_gene\_id as key and orthologs gene ids as values


**filter\_orthologs.py**

Input: orthologs.json

Filter criteria for human protein coding gene:

-   should have an ortholog in mouse

-   should have orthologs for all species in group or none (pairs
    > defined in ensembl\_stable\_id\_species.json) “excluded pairs”

-   minimum of 10 orthologs (removed if &lt;10)

-   maximum of 3 orthologs in a species (removed if &gt;3 “excluded orth
    > rel”

Output:

-   orthologs\_filtered.json

-   excluded\_pairs.json

-   excluded\_orth\_rel.json

How to check whether is completed: Not. Filtered set &lt; full set


**retrieve\_genetrees.py**

Queries ensembl API for genetrees of the filtered orthologs.

Input: orthologs\_filtered.json

Output:

-   ensembl\_api/{gene\_id}.json

-   genes\_to\_genetrees.json (mapping)

-   genetrees\_ogs.json (mapping)

-   redo\_orth\_parse.json

Completed if: ensembl\_api/{gene\_id}.json file is created for all
gene\_ids in orthologs\_filtered.json

Downloading and saving all raw data is necessary for reproducibility

--- Snakefile 1 ended, all files necessary to check if next steps are
done for all proteins.

## Repeat detection (Snakefile)

**parse\_genetree.py** *\$1: ensembl\_api/{gene\_id}.json*

Subtract smallest gene tree still containing all orthologs. Acquire
fasta sequences for proteins in gene tree.

Advantage: can be called for a single gene tree.

Disadvantage: need to load mapping each time.

Input:

-   orthologs\_filtered.json

-   ensembl\_api/{gene\_id}.json

Output:

-   genetrees\_nhx/{genetree\_id}\_{gene\_id}.nhx

-   fasta\_ogs/{genetree\_id}\_{gene\_id}.fa

Completed if: genetrees\_nhx/{og\_id}.nhx and fasta\_ogs/{og\_id}.fa for
every og\_id in

genes\_to\_genetrees.json where {og\_id} = {genetree\_id}\_{gene\_id}

**(DEP) hmmscan\_pfam.bash **

fully replaced by snakemake rule:

Input:

-   Pfam-A.hmm (+ pressed files)

-   fasta\_ogs/{og\_id}.fa

Output:

-   pfam/hmm/{og\_id}.tblout

'hmmscan -T 12.5 --domT 0 --domtblout {output} '+pfam\_file+' {input}
&&gt;/dev/null'

Completed if: pfam/hmm/{og\_id}.tblout for og\_id in
genes\_to\_genetrees.json

**run\_pfam\_iteration\_single.bash** *\$1: pfam/hmm/{og\_id}.tblout*

Iterative pfam domain annotation. Improvement on regular hmmscan for
more accurate repeat detection.

Scripts:

-   parse\_tblout\_init.py *\$1: pfam/hmm/{og\_id}.tblout \$2:
    > fasta\_ogs/{og\_id}.fa*

    -   Use sequences gathering cutoff, select best Pfam model for each
        > clan based on sequence bitscore maximalisation (over-all).

    -   Filter Pfam hits based on minimum \#repeats (&gt;=3) in human,
        > presence in mouse and general filtering rule: &gt;\# repeats
        > in &gt; \#species (3 and 1, respectively so not currently
        > used.

    -   Writes a fasta file for each best hit of clan with all of the
        > repeats

    -   hmm\_results.json before any filtering

<!-- -->

-   parse\_tblout\_iterate.py *\$1: pfam/tblout{og\_hit\_id}.tblout \$2:
    > fasta\_ogs/{og\_id}.fa*

    -   extract repeats from tblout, same as happens in
        > parse\_tblout\_init

-   parse\_tblout\_final.py *\$1: pfam/tblout{og\_hit\_id}.tblout \$2:
    > fasta\_ogs/{og\_id}.fa*

    -   Use more liberal slicing for making tree, envelope spacing and
        > padding 5

    -   hmm\_results\_final.json

Input:

-   pfam/hmm/{og\_id}.tblout

-   fasta\_ogs/{og\_id}.fa

-   Pfam gathering cutoff: ga\_cutoff.tsv

-   Pfam clans: Pfam-A.clans.tsv

Output:

-   progress\_files/{og\_id}.pfamhmm.done

-   pfam/repeats/{og\_id\_hit}.fa parse\_tblout\_init.py

    -   serves as input for iteration

-   pfam/aligned/{og\_id\_hit}.linsi.fa

-   pfam/alignments/{og\_id\_hit}.linsi.fa last copy of pfam/aligned

-   pfam/tblout/{og\_id\_hit}.tblout

-   hmm\_results.json parse\_tblout\_init.py

-   excluded\_repeats.json parse\_tblout\_init.py

-   hmm\_results\_final.json parse\_tblout\_final.py

Alignment with mafft linsi: mafft --quiet --localpair --maxiterate 1000

.fa.old as temporary file in /pfam/repeats/

Iterates max 20 times or until convergence (identical output file as
.fa.old)

Build profile HMMs based on alignment hmmbuild --fast --symfrac 0.6 -n

Check for length of profile, prevent profile drift: enforce minimum
length of 20

Completed if: dummy file is made for all {og\_id} from
genes\_to\_genetrees.json

**run\_iqtree.bash**

Can be replaced by snakemake rule

Input:

-   pfam/aligned/{og\_id\_hit}.linsi.fa

Output:

-   pfam/aligned/{og\_id\_hit}.linsi.fa.treefile and .iqtree and lots of
    > others

-   pfam/trees/{og\_id\_hit}.treefile

iqtree-omp -s {input} -bb 1000 -nt AUTO -mset LG,WAG,VT,Dayhoff,JTT
-redo

copy to pfam/trees is separate rule

Completed if: .treefile is made for all {og\_id\_hit} from
hmm\_results\_final.json

**run\_treefix.bash**

split into three snakemake rules:

-   prepare\_treefix:

    -   input:

        -   pfam/aligned/{og\_hit}.linsi.fa.treefile

        -   genetrees\_nhx/{og}.nhx

    -   output:

        -   pfam/aligned/{og}\_{hit}.linsi.fa.treefile.rooted

        -   pfam/treefix/{og}.stree pfam/treefi/{og}.smap (*not in
            > snakemake)*

    -   script: python generate\_treefix\_input.py {input.treefile}
        > {input.genetree}

<!-- -->

-   run\_treefix:

    -   input:

        -   pfam/aligned/{og}\_{hit}.linsi.fa.treefile.rooted

        -   pfam/treefix/{og}.stree pfam/treefi/{og}.smap

        -   pfam/aligned/{og}\_{hit}.linsi.fa.iqtree

        -   *(also uses alignment but not specified as input)*

    -   output:

        -   pfam/treefix/{og}\_{hit}.treefix.tree

        -   *+ other files?*

    -   shell:

model=\`cat {input.iqtree} | grep -oP "Best-fit model according to BIC:
\\K.\*"\`; treefix -s {input.stree} -S {input.smap} -A .linsi.fa -o
.linsi.fa.treefile.rooted -n .treefix.tree -V 0 -m
treefix.models.iqtreemodel.CoarseModel -e \\'-t AA -m \\'\$model\\'\\'
{input.treefile}

move with seperate rule

-   run\_treefix\_annotate:

    -   input:

        -   pfam/treefix/{og}\_{hit}.treefix.tree

        -   pfam/treefix/{og}.stree pfam/treefi/{og}.smap

    -   output

        -   pfam/treefix/{og}\_{hit}.treefix.mpr.recon

        -   *+ others*

    -   shell:

tree-annotate -s {input.stree} -S {input.smap} {input.tree}

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

**Meme dataset construction:**

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


