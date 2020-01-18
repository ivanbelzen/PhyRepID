# Analysis scripts PhyRepID


## generate_maptotree.py
Generate gene trees with reconciliation annotatation. 
Note: not always correct placement if multiple branches of gene tree have the same branch name, e.g. in cases of gene duplications.

- Input:
  - Look up: $1 {og\_id\_hit} 
  - genetrees\_nhx/{og\_id}.nhx
  - pfam_evo_events.json
  - meme_evo_events.json

- Outputs in /maptotree/
  - {og\_id\_hit}_maptotree.nhx
  - {og\_id\_hit}_maptotree.pdf


## analyse_schaper_comparison.py

Analyse Schaper *et. al.* output for overlap with PhyloRepID output.
To compare: collapse on OG level so for every OG, consider the human proteins and what Schaper's annotation was on the cherries with these human proteins and vertebrates. Ignore none and not assigned, only look at  slightly/perfectly separated and conserved.  

Schaper summary contains info per pair of orthologs by protein ID.
PhyloRepID evo events summary contains dup/loss per node/ortholog in an OG and in total. 
Note: "positive" in Schaper if any of species in OG is separated from human for any domain according to their dataset.

- Input:
  - Schaper files: 
    - resources/eukaryotic_pairwise_repeat_unit_phylogenies_PFAM.newick.gz 
    - resources/eukaryotic_pairwise_repeat_unit_phylogenies_denovo.newick.gz
  - og_human_mapping_reroot.json
(made by generate_human_protein_mapping.py)
- Output:
  - schaper_comparison.json (used in phyrepid export) 
  - schaper_detailed_comparison.json
  - phyrepid_schaper_comparison.tsv"

Schaper E, Gascuel O, Anisimova M. Deep conservation of human protein tandem repeats within the eukaryotes. Mol Biol Evol. 2014;31: 1132–1148.

## generate_itol_domain_templates.py
Generate domain templates for easy visualisation in ITOL
- Input
    - $1: identifier {og_id_hit} 
    - pfam/meme_repeat_stats.json
    - pfam_repeat_stats_initial_filtered.json (optional)
 - Output path /itol/
    - {og_id_hit}_domains.txt
    - {og_id_hit}_domains_initial.txt
    - {og_id_hit}_repeat_labels.txt
    - {og_id_hit}_genetree_labels.txt
 


# Draft ##
## analyse_dataset_characteristion.py
Makes dataframes and extensive analysis 
in detailed_df. Compares initial dataset to final dataset (figures are deprecated)

*Makes unit_cv dataframe used in export*
Experiments with lot of different analyses and metrics.

- Input (from analyse_pfam/meme_summary.py) 
  - pfam_summary.json
  - pfam_summary_initial.json
  - meme_summary.json
- Output:
  - pfam_detailed_df.json and pfam_detailed_df.csv 
  - initial_pfam_filtered.json
(filtered_pfam_summary - matched on OG ids from initial to the OG-domain ids from later to compare improvement in repeat detection - can only be fore pfam )
  - initial_detailed_df.json (filtered dataframe - pfam only)
  - unfiltered_initial_detailed_df.json
(UNfiltered dataframe -pfam only )
  - meme_detailed_df.json and meme_detailed_df.csv 


Files with numbers and stats like analysis/dataset_24102018/regular_dataset.txt

- Files used for follow-up:
  - unit_cv dataframe  *used in export*
 
 Detailed dataframes only used in the file itself. 

 


