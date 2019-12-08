### Sup. Table 1 PRD score ranking

SupTable_1_PRD_score_ranking.tsv

Ranking of all OGs according to PRD score, annotated with duplications and losses as inferred by tree reconciliation of the repeat tree with the gene tree. 

PRD score is the relative measure of protein domain evolution for the repeat region the OG.

PRD_score = (netto_dup - mean_dup) / othologs_cnt

where mean_dup (arithmetic average of netto duplications in the dataset) is a constant

  

File made by script: …

  
  

Columns

-   genetree_id: <identifier>
    

Persistent identifier (PID) in ENSEMBL database: prefix

  

-   gene_symbol: <text>
    

HGNC approved gene symbol, used for display purposes only

  

-   prd_score: <float>
    

Quantification of protein repeat evolution according in this repeat region the OG

prd_score = (netto_dup - mean_dup) / othologs_cnt

where mean_dup (mean duplications in the dataset) is a constant

  

-   identifier: <identifier>
    

Internal identifier used in the pipeline <genetree_id>_<gene_id>_<pfam_model_id>

  

-   netto_dup: <integer>
    

Duplications since the most recent common ancestor in the gene tree

  

-   loss: <integer>
    

Number of losses as inferred by tree reconciliation of the repeat tree with the gene tree

  

-   orthologs_cnt: <integer>
    

Number of proteins in the orthologous group

  

-   clan
    

PFAM clan that the Pfam model is associated which makes op the repeat

  

-   unit_cv: <integer>
    

Co-variance of the number of repeat units (occurrences) there are in the OG

  

-   selectome: [0,1]
    

Boolean for positive selection according to Selectome database (Moretti et al. 2014)

  

-   schaper_pos: [0,1]
    

Boolean for >=1 human protein “perfectly or strongly separated” in (Schaper et al. 2014) from another protein in the OG.

  

-   mis_z: <float>
    

Missense mutation z-score in ExAC (Karczewski et al. 2016)

  

-   pLi: <float>
    

Sensitivity to loss of function heterozygosity in ExAC (Karczewski et al. 2016)

  

-   gt_dup: <integer>
    

Number of duplications after reconciliation of gene tree with the species tree

  

-   gt_root: <identifier>
    

Most recent common ancestor in the gene tree
