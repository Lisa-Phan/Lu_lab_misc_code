#!/bin/bash

#flow control
#Work on assumption that analysis proceed regardless of fastqc results
#1. Get fastqc job id and store in log file
#2. Fill out template file using previous job script and id depencencies

#Usage: job_submission <step number> <name of job> <name of script> <name of previous script> 
job_submission () {
	sleep 5   
	echo step $1: submitting job $2 
	id=`sbatch $3 --dependency=afterok:$id $4 | tail -n 1 | cut -d' ' -f4` 
	echo job id is $id
}

#path starting from main directory, which is one level above 'scripts'
fastqc=scripts/mod_fastqc.sh
assembly=scripts/mod_assembly.sh
cutadapt=scripts/mod_cutadapt.sh
fastapta_clust=scripts/mod_fastaptamer_cluster.sh

echo step 1: submitting job fastqc...
id=`sbatch $fastqc | tail -n 1 | cut -d' ' -f4`
echo job id is $id

job_submission 2 'pair end assembly' $assembly $fastqc
job_submission 3 'cutadapt' $cutadapt $assembly
job_submission 4 'fastaptamer cluster' $fastapta_clust $cutadapt
