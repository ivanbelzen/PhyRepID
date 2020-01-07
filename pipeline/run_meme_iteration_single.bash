#!/bin/bash

# IvB April 2018
# Iterative pfam domain annotation
# Starts from domain tblout from hmmscan, parse it with python and orthologs fasta file.
# Alignment with mafft linsi
# Needs: full pfam hmm .tblout  ogs fasta file
# scripts: parse_tblout_init.py parse_tblout_iterate.py parse_tblout_final.py 
# Outputs: repeats fasta, profile hmm, aligned repeats linsi.fa
# Logs: hmm_results.json before iteration and hmm_results_final.json after iteration and filtering.

# parse_tblout scripts output to /pfam/repeats/*.fa
# init takes full tblout as input /pfam/hmm/ from Pfam-A scan and iterate/final scripts take partial tblouts /pfam/tblout/ from iteration

# to detect iteration is done, final mafft is copied to /pfam/alignments
shopt -s nullglob

# configuration file

root=/home/ianthe/protein-repeat-evolution/

orthologs_path=${root}fasta_ogs/

pfam_path=${root}denovo/meme/hmm/
tblout_path=${root}denovo/meme/tblout/
hmm_path=${root}denovo/meme/profiles/
repeats_path=${root}denovo/meme/repeats/
aligned_path=${root}denovo/meme/aligned/

final_alignments_path=${root}denovo/meme/alignments/


threshold_length=75	#minimum percentage of start hmm profile length 

file=$1			#denovo/meme/repeats/{ogs_id}.fa

ogs_id=${file//${repeats_path}/}
ogs_id=${ogs_id//.fa/} 	#genetree_gene
	
echo ${ogs_id} #progress indicator

i=20
while : # while files are different 
do
	if [[ i -eq 0 ]]; then break; fi

	cp ${file} ${file}.old
	#compare new .fa with old .fa in loop
	
	mafft --quiet --localpair --maxiterate 1000 ${file} > ${aligned_path}${ogs_id}.linsi.fa 			
				
	if [[ -f ${hmm_path}${ogs_id}.hmm ]]; then
		rm ${hmm_path}${ogs_id}.hmm
	else 
		start_length=0 	#if first time iteration, reset start_length to save later
	fi		
			
	hmmbuild --fast --symfrac 0.6 -n motif ${hmm_path}${ogs_id}.hmm ${aligned_path}${ogs_id}.linsi.fa 		
			
	#check for length of profile, prevent profile drift
	leng_str=`grep -e "LENG [0-9]*" ${hmm_path}${ogs_id}.hmm | head -1`
	length=${leng_str/LENG /}
			
	#enforce minimum length of 15
	if [[ ${length} -lt 15 ]]; then	
		break
	fi
					
	hmmpress -f ${hmm_path}${ogs_id}.hmm 
	hmmscan --domT 0 --domtblout ${tblout_path}${ogs_id}.tblout ${hmm_path}${ogs_id}.hmm ${orthologs_path}${ogs_id}.fa &>/dev/null 
			
	python ${root}parse_tblout_iterate.py ${tblout_path}${ogs_id}.tblout ${orthologs_path}${ogs_id}.fa
			
	diff ${file} ${file}.old  &>/dev/null

	if [ $? -eq 0 ]; then
		echo "Convergence reached"
		break				
	fi
	((i--))
done
		
rm ${file}.old
		
python ${root}parse_tblout_final.py ${tblout_path}${ogs_id}.tblout ${orthologs_path}${ogs_id}.fa
	
mafft --quiet --localpair --maxiterate 1000 ${file} > ${aligned_path}${ogs_id}.linsi.fa 
cp ${aligned_path}${ogs_id}.linsi.fa ${final_alignments_path}${ogs_id}.linsi.fa  		
			

