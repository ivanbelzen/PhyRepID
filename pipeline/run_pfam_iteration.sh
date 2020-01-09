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
source config.sh

threshold_length=75	#minimum percentage of start hmm profile length 

file=$1			#pfam/hmm/{genetree}{gene}.tblout
echo ${file}
ogs_id=${file//${pfam_path}/}
ogs_id=${ogs_id//.tblout/} 	#genetree_gene
	
echo ${ogs_id} #progress indicator
	
python ${root}parse_tblout_init.py ${file} ${orthologs_path}${ogs_id}.fa 
#full tblout as input
#outputs to pfam/repeats/{ogs_id}.fa

for file in ${pfam_repeats_path}${ogs_id}*	#should be only one, but you don't know name
do
	og_hit=${file//${pfam_repeats_path}/} # genetree _ human gene _ hit pfam domain 
	og_hit=${og_hit//.fa/}
	hit_only=${og_hit//${ogs_id}_/}
	echo ${hit_only}
		

	i=20
	while : # while files are different 
	do
		if [[ i -eq 0 ]]; then break; fi

		cp ${file} ${file}.old
		#compare new .fa with old .fa in loop
	
		mafft --quiet --localpair --maxiterate 1000 ${file} > ${pfam_aligned_path}${og_hit}.linsi.fa 			
				
		if [[ -f ${pfam_hmm_path}${og_hit}.hmm ]]; then
			rm ${pfam_hmm_path}${og_hit}.hmm
		else 
			start_length=0 	#if first time iteration, reset start_length to save later
		fi		
			
		hmmbuild --fast --symfrac 0.6 -n ${hit_only} ${pfam_hmm_path}${og_hit}.hmm ${pfam_aligned_path}${og_hit}.linsi.fa 		
			
		#check for length of profile, prevent profile drift
		leng_str=`grep -e "LENG [0-9]*" ${pfam_hmm_path}${og_hit}.hmm | head -1`
		length=${leng_str/LENG /}
			
		#enforce minimum length of 20
		if [[ ${length} -lt 20 ]]; then	
			break
		fi
					
		hmmpress -f ${pfam_hmm_path}${og_hit}.hmm 
		hmmscan --domT 0 --domtblout ${pfam_tblout_path}${og_hit}.tblout ${pfam_hmm_path}${og_hit}.hmm ${orthologs_path}${ogs_id}.fa &>/dev/null 
			
		python ${root}parse_tblout_iterate.py ${pfam_tblout_path}${og_hit}.tblout ${orthologs_path}${ogs_id}.fa
			
		diff ${file} ${file}.old  &>/dev/null

		if [ $? -eq 0 ]; then
			echo "Convergence reached"
			break				
		fi
		((i--))
		
	done
		
	rm ${file}.old
		
	python ${root}parse_tblout_final.py ${pfam_tblout_path}${og_hit}.tblout ${orthologs_path}${ogs_id}.fa
	
	mafft --quiet --localpair --maxiterate 1000 ${file} > ${pfam_aligned_path}${og_hit}.linsi.fa 
	cp ${pfam_aligned_path}${og_hit}.linsi.fa ${pfam_final_alignments_path}${og_hit}.linsi.fa  		
			
done

touch ${progress_path}${ogs_id}.pfamhmm.done	
