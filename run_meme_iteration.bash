#!/bin/bash

# IvB May 2018
# Adjustment from Pfam iteration but for meme repeats in denovo folder

# Iterative pfam domain annotation
# Starts from domain tblout from hmmscan, parse it with python and orthologs fasta file.
# Alignment with mafft linsi
# Needs: motif hmm .tblout  ogs fasta file
# scripts: parse_denovo_tblout_iterate.py parse_denovo_tblout_final.py 
# Outputs: repeats fasta, profile hmm, aligned repeats linsi.fa
# Logs: meme_hmm_results.json before filter/iteration and hmm_meme_results_final.json after iteration and filtering.

# parse_tblout scripts output to /denovo/repeats/*.fa
# repeats are already filtered and parsed by parse_denovo_tblout_init.py
# iterate/final scripts take tblouts /denovo/tblout/ from iteration
# to detect iteration is done, final mafft is copied to /denovo/alignments and progress_files

shopt -s nullglob

# configuration file

root=/home/ianthe/protein-repeat-evolution/
orthologs_path=${root}fasta_ogs/
tblout_path=${root}denovo/tblout/
hmm_path=${root}denovo/profiles/
meme_repeats_path=${root}denovo/meme/repeats/
repeats_path=${root}denovo/repeats/
aligned_path=${root}denovo/aligned/
final_alignments_path=${root}denovo/alignments/
progress_path=${root}/progress_files/

threshold_length=75	#minimum percentage of start hmm profile length 

file=$1			#denovo/repeats/{genetree}_{gene}.fa
echo ${file}
ogs_id=${file//${repeats_path}/}
ogs_id=${ogs_id//.fa/} 	#genetree_gene

i=20
while : # while files are different 
do
	if [[ i -eq 0 ]]; then break; fi

	cp ${repeats_path}${ogs_id}.fa ${repeats_path}${ogs_id}.fa.old
	#compare new .fa with old .fa in loop
	
	mafft --quiet --localpair --maxiterate 1000 ${repeats_path}${ogs_id}.fa > ${aligned_path}${ogs_id}.linsi.fa 			
				
	if [[ -f ${hmm_path}${ogs_id}.hmm ]]; then
		rm ${hmm_path}${ogs_id}.hmm
	fi		
			
	hmmbuild --fast --symfrac 0.6 -n 'motif' ${hmm_path}${ogs_id}.hmm ${aligned_path}${ogs_id}.linsi.fa 		
			
	#check for length of profile, prevent profile drift
	leng_str=`grep -e "LENG [0-9]*" ${hmm_path}${ogs_id}.hmm | head -1`
	length=${leng_str/LENG /}
			
	#enforce minimum length of 15
	if [[ ${length} -lt 15 ]]; then	
		break
	fi
					
	hmmpress -f ${hmm_path}${ogs_id}.hmm 
	hmmscan --domT 0 --domtblout ${tblout_path}${ogs_id}.tblout ${hmm_path}${ogs_id}.hmm ${orthologs_path}${ogs_id}.fa &>/dev/null 
			
	python ${root}parse_denovo_tblout_iterate.py ${tblout_path}${ogs_id}.tblout ${orthologs_path}${ogs_id}.fa
			
	diff ${repeats_path}${ogs_id}.fa ${repeats_path}${ogs_id}.fa.old  &>/dev/null

	if [ $? -eq 0 ]; then
		echo "Convergence reached"
		break				
	fi
	((i--))
		
done
		
rm ${repeats_path}${ogs_id}.fa.old
	
python ${root}parse_denovo_tblout_final.py ${tblout_path}${ogs_id}.tblout ${orthologs_path}${ogs_id}.fa
	
mafft --quiet --localpair --maxiterate 1000 ${repeats_path}${ogs_id}.fa > ${aligned_path}${ogs_id}.linsi.fa 
cp ${aligned_path}${ogs_id}.linsi.fa ${final_alignments_path}${ogs_id}.linsi.fa  		

touch ${progress_path}${ogs_id}.meme_hmm.done	
