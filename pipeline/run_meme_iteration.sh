#!/bin/bash

# IvB May 2018
# Adjustment from Pfam iteration but for meme repeats in denovo folder

# Iterative domain annotation
# parse_tblout scripts output to /denovo/repeats/*.fa
# repeats are already filtered and parsed by parse_denovo_tblout_init.py
# iterate/final scripts take tblouts /denovo/tblout/ from iteration
# to detect iteration is done, final mafft is copied to /denovo/alignments and progress_files

# Alignment with mafft linsi
# Needs: motif hmm .tblout  ogs fasta file
# scripts: parse_denovo_tblout_iterate.py parse_denovo_tblout_final.py 
# Outputs: repeats fasta, profile hmm, aligned repeats linsi.fa
# Logs: meme_hmm_results.json before filter/iteration and meme_hmm_results_final.json after iteration and filtering.


shopt -s nullglob

# configuration file
source config.sh

threshold_length=75	#minimum percentage of start hmm profile length 

file=$1			#denovo/repeats/{genetree}_{gene}.fa
echo ${file}
ogs_id=${file//${denovo_repeats_path}/}
ogs_id=${ogs_id//.fa/} 	#genetree_gene

i=20
while : # while files are different 
do
	if [[ i -eq 0 ]]; then break; fi

	cp ${denovo_repeats_path}${ogs_id}.fa ${denovo_repeats_path}${ogs_id}.fa.old
	#compare new .fa with old .fa in loop
	
	mafft --quiet --localpair --maxiterate 1000 ${denovo_repeats_path}${ogs_id}.fa > ${denovo_aligned_path}${ogs_id}.linsi.fa 			
				
	if [[ -f ${denovo_hmm_path}${ogs_id}.hmm ]]; then
		rm ${denovo_hmm_path}${ogs_id}.hmm
	fi		
			
	hmmbuild --fast --symfrac 0.6 -n 'motif' ${denovo_hmm_path}${ogs_id}.hmm ${denovo_aligned_path}${ogs_id}.linsi.fa 		
			
	#check for length of profile, prevent profile drift
	leng_str=`grep -e "LENG [0-9]*" ${denovo_hmm_path}${ogs_id}.hmm | head -1`
	length=${leng_str/LENG /}
			
	#enforce minimum length of 15
	if [[ ${length} -lt 15 ]]; then	
		break
	fi
					
	hmmpress -f ${denovo_hmm_path}${ogs_id}.hmm 
	hmmscan --domT 0 --domtblout ${denovo_tblout_path}${ogs_id}.tblout ${denovo_hmm_path}${ogs_id}.hmm ${orthologs_path}${ogs_id}.fa &>/dev/null 
			
	python ${root}parse_denovo_tblout_iterate.py ${denovo_tblout_path}${ogs_id}.tblout ${orthologs_path}${ogs_id}.fa
			
	diff ${denovo_repeats_path}${ogs_id}.fa ${denovo_repeats_path}${ogs_id}.fa.old  &>/dev/null

	if [ $? -eq 0 ]; then
		echo "Convergence reached"
		break				
	fi
	((i--))
		
done
		
rm ${denovo_repeats_path}${ogs_id}.fa.old
	
python ${root}parse_denovo_tblout_final.py ${denovo_tblout_path}${ogs_id}.tblout ${orthologs_path}${ogs_id}.fa
	
mafft --quiet --localpair --maxiterate 1000 ${denovo_repeats_path}${ogs_id}.fa > ${denovo_aligned_path}${ogs_id}.linsi.fa 
cp ${denovo_aligned_path}${ogs_id}.linsi.fa ${denovo_final_alignments_path}${ogs_id}.linsi.fa  		

touch ${progress_path}${ogs_id}.meme_hmm.done	
