#!/bin/bash

# IvB 
# Setup pipeline  => also adjust pipeline_methods.py root 
# configuration file

root=

orthologs_path=${root}fasta_ogs/
progress_path=${root}progress_files/

pfam_path=${root}pfam/hmm/
pfam_tblout_path=${root}pfam/tblout/
pfam_hmm_path=${root}pfam/profiles/
pfam_repeats_path=${root}pfam/repeats/
pfam_aligned_path=${root}pfam/aligned/
pfam_final_alignments_path=${root}pfam/alignments/

denovo_tblout_path=${root}denovo/tblout/
denovo_hmm_path=${root}denovo/profiles/
meme_repeats_path=${root}denovo/meme/repeats/
denovo_repeats_path=${root}denovo/repeats/
denovo_aligned_path=${root}denovo/aligned/
denovo_final_alignments_path=${root}denovo/alignments/
