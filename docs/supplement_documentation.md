# Supplement documentation

## Sup. table 1: PRD score ranking

Ranking of all OGs according to PRD score, annotated with duplications and losses as inferred by tree reconciliation of the repeat tree with the gene tree. 

PRD score is the relative measure of protein domain evolution for the repeat region the OG.

<p align="center"> 
<b>PRD score = (x - u) / n;</b><br />
 with x = netto duplications since most recent common ancestor (gene tree root); u = mean duplications in the full dataset, n = number of proteins in the OG.</p>
 
File made during **Post-processing - generate_phyrepid_results**

Columns:
-   genetree_id: \<identifier> 
Persistent identifier (PID) in ENSEMBL database: prefix.
-   gene_symbol: \<text>
HGNC approved gene symbol, used for display purposes only
-   prd_score: \<float>
Quantification of protein repeat evolution according in this repeat region the OG
-   identifier: \<identifier>
 Internal identifier used in the pipeline \<genetree_id>\_<gene_id>\_\<hit>
 -   netto_dup: \<integer>
Duplications since the most recent common ancestor in the gene tree
-   loss: \<integer>
Number of losses as inferred by tree reconciliation of the repeat tree with the gene tree.
-   orthologs_cnt: \<integer>
Number of proteins in the orthologous group.
-   clan: \<string> 
Pfam clan of the Pfam model which makes up the repeat.
-   unit_cv: \<integer>
Co-variance of the number of repeat units (occurrences) there are in the OG
-   selectome: [0,1]
Boolean for positive selection according to Selectome database (Moretti et al. 2014)
-   schaper_pos: [0,1]
Boolean for >=1 human protein “perfectly or strongly separated” in (Schaper et al. 2014) from another protein in the OG.
-   mis_z: \<float>
Missense mutation z-score in ExAC (Karczewski et al. 2016)
-   pLi: \<float>
Sensitivity to loss of function heterozygosity in ExAC (Karczewski et al. 2016).
-   gt_dup: \<integer>
 Number of duplications after reconciliation of gene tree with the species tree.
-   gt_root: \<identifier>
Most recent common ancestor in the gene tree

## Sup. table 2: Human full lineage  

File made during **Post-processing - generate_phyrepid_results_human**

First columns are the same as in Sup. table 1.

 - gene_symbol
 - identifier
 - clan
 - prd_score
 - netto_dup
 - orthologs_cnt
-   human_dup:  \<integer> Duplications on the full human lineage.
-   human_frac:  \<float> Fraction of netto_dup from human_dup.
 - human_protein_cnt: \<integer>
Number of human proteins in the orthologous group.

## Sup. table 3: Human branch only 

File made during **Post-processing - generate_phyrepid_results_human**

First columns are the same as in Sup. table 1.

Ranking: filtered on duplication_score > 0 and sorted by human_frac, descending.

 - gene_symbol
 - identifier
 - clan
 - prd_score
 - netto_dup
 - orthologs_cnt
-   human_dup:  \<integer> Duplications on human branch.
-   human_frac:  \<float> Fraction of netto_dup from human_dup.
 - human_protein_cnt: \<integer>
Number of human proteins in the orthologous group.

## Sup. tables 5 and 6 Repeat analysis Pfam and MEME

File made during **Post-processing - analyse_repeat_stats**

First columns are the same as in Sup. table 1.
 
- row_id:  \<identifier>. Unique row number. 
- identifier: \<identifier>. Internal identifier used in the pipeline \<genetree_id>\_<gene_id>\_\<hit>
- ortholog_id: \<identifier> 
Persistent identifier (PID) of gene in ENSEMBL database.
- species: \<identifier> Species abbreviation.
- clan: \<string> 
*(Pfam only)* Pfam clan of the Pfam model which makes up the repeat.
- model: \<string> 
*(Pfam only)* Pfam model which makes up the repeat.
- unit_cnt: \<integer> Number of repeat units.
- model_length: \<integer> Length of repeat unit hidden markov model in  amino acids.
- max_dist: \<integer> Number of amino acids between two repeat units, maximum.
- protein_length:  \<integer> Protein length in amino acids.
- protein_cov: \<float> Fraction of protein length covered by repeat units.
- repeat_aa :  \<integer> Total number of amino acids covered by repeat units.
- repeat_cov: \<float> Fraction of total repeat length covered by repeat units.
