# Analysis scripts PhyRepID


## analyse_maptotree.py

Generate gene trees with reconciliation annotatation. 
Note: not always correct placement if multiple branches of gene tree have the same branch name, e.g. in cases of gene duplications.

- Input:
  - Look up: $1 {og}\_{hit} 
  - genetrees\_nhx/{og\_id}.nhx
  - pfam_evo_events_reroot.json
  - meme_evo_events_reroot.json
 (Uses events and not summary stats so root node does not matter)

- Outputs in analysis/maptotree/
  - ${og_id}_maptotree_abs.nhx
  - ${og_id}_maptotree_abs.pdf


## analyse_merge_export.py
 Makes table with all columns: ExAC, Selectome, Schaper, 

- Input
  - pfam/meme_summary_df_reroot.json
 (made by analyse_evo_events_summary.py)
   - schaper_summary_df_reroot.json
(made by analyse_schaper_comparision.py)
  - selectome_df_reroot.json (made analyse_selectome.py)
  - analysis/exac_10102018/exac_df.json
(made by analyse_exac.py)
  - genetree_df_reroot.json
(made by analyse_genetrees.py from  reconciled gene trees) 

- Output
   - /analysis/evo_events/combined_summary_df_reroot.json
   - /analysis/evo_events_06102018/pre_dataset_df.json
   - evo_events_merged.csv (combined_summary_df export)
    


## analyse_schaper_comparison.py
 
Analyse Schaper et al output for overlap with PhyloRepID output.
To compare: collapse on OG level so for every OG, consider the human proteins and what Schaper's annotation was on the cherries with these human proteins and vertebrates. Ignore none and not assigned, only look at  slightly/perfectly separated and conserved 

Schaper summary contains info per pair of orthologs by protein ID.
PhyloRepID evo events summary contains dup/loss per node/ortholog in an OG and in total 

- Input:
  - analysis/human_evo_events/full_lineage_trace/human_summary_df.json
  - og_human_mapping_reroot.json'
(OG to human gene mapping made by generate_human_protein_mapping.py)
  - schaper_summary(_denovo).json
(made by extract_schaper.py)

combined_summary_df_file = pre.root+'analysis/evo_events_06102018/pre_dataset_df.json'


- Output:
  - schaper_detailed_df_reroot.json
  - schaper_summary_df_reroot.json
(used in merge_export)



## get_human_evo_events.py
Analyses only evens in human lineage (full lineage) or branch (human only) and calculates adjusted PRD score. 

- Input: 
  - og_human_mapping_reroot.json (from generate_human_og_mapping.py)
  - pfam_evo_events_reroot.json
  - meme_evo_events_reroot.json
- Output: 
  - human_summary_df_human_only.json
  - human_summary_df_full_lineage.json
  - human_full_lineage.csv
  - human_only.csv

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
 
## analyse_pfam(/meme)_summary.py 

Glob tblout files to derive characteristics of repeats and put in summary.json file.
Only process tblout from OGs that have Treefix annotate output.

Note: dataframes made in analyse_dataset_characteristion.py

- Input:
  - pre.pfam_tblout_path+og_hit_id+'.tblout'
- Output:
   - pfam_summary.json
  - pfam_exclude_ogs_lt4.txt (less than 4 repeats)

For meme same

## analyse_evo_events_summary.py 

Makes pfam/meme_summary_df_reroot.json
from
 pfam/meme_evo_events_reroot.json 

And additional analyses/figures which are deprecated and replaced by R script

Note: Confusing naming but *_summary_df_reroot.json* is not *_summary.json*. Former contains evo_event summary and latter repeat unit detection statistics
 
## generate_itol_domain_templates.py
Generate domain templates for easy visualisation in ITOL
using pfam summary.json and meme_summary.json
=> iqtree is un-reconciled and only for comparison

---

**Dictionaries**

-   pfam\_summary.json

    -   key: og\_id\_domain
    -   repeats dict containing: domain = { clan, length,
        -   orthologs\_dict = { protein\_uri : {unit\_count, seq bitscore,
            -   units\_dict = { unit\_number: { dom bitscore, model coordinates, seq coordinates }, ... } }, ...} }


