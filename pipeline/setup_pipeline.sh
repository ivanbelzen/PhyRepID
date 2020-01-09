#!/bin/bash

# IvB 
# Setup pipeline  => also adjust pipeline_methods.py root 
# make directories
# download Pfam hmm profiles

#configuration file
source config.sh

mkdir -p ${root}resources/
mkdir -p ${root}docs/
mkdir -p ${root}logs/

mkdir -p ${progress_path}
mkdir -p ${orthologs_path}
mkdir -p ${root}fasta_ogs_chopped/
mkdir -p ${root}genetrees_nhx/
mkdir -p ${root}ensembl_api/  
mkdir -p ${root}maptotree/

mkdir -p ${pfam_path}
mkdir -p ${pfam_tblout_path}
mkdir -p ${pfam_hmm_path}
mkdir -p ${pfam_repeats_path}
mkdir -p ${pfam_aligned_path}
mkdir -p ${pfam_final_alignments_path}
mkdir -p ${root}pfam/trees/
mkdir -p ${root}pfam/treefix/adjusted/images/

mkdir -p ${denovo_tblout_path}
mkdir -p ${denovo_hmm_path}
mkdir -p ${meme_repeats_path}
mkdir -p ${denovo_repeats_path}
mkdir -p ${denovo_aligned_path}
mkdir -p ${denovo_final_alignments_path}
mkdir -p ${root}denovo/trees/
mkdir -p ${root}denovo/treefix/adjusted/images/

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
mv Pfam-A.hmm.gz ${root}resources/.
gunzip ${root}resources/Pfam-A.hmm.gz 
