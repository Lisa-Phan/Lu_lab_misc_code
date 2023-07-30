#!/bin/bash

#SBATCH -J fastqc
#SBATCH -o ../io_log/%x-%j.out
#SBATCH -e ../io_log/%x-%j.err
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:30:00


##Comment block
##last mod 4/27/2023
##NOTES:
##out directory has to already exists or script will not work
##file name must not contain spaces
##assuming fastqc is already installed and can be called from anywhere


echo Step 0: Define variables, tools, and paths at echo `date`.

## get variable name from hidden log file
## symlink is sym link to the raw data, equivalent of input_dir, path is the working directory foldername
JOB_NAME=`tail -n 1 ../.userlog.txt | cut -d" " -f4`
REP_OPT=`tail -n 1 ../.userlog.txt | cut -d" " -f6`
OUT_DIR=../$JOB_NAME/step1_fastqc

echo job name is $JOB_NAME, with repeat option $REP_OPT

## Two options for input_ext, if rep_opt is 1 then input_ext is FASTQ
## if rep_opt is 2 then input_ext is _001.fastq
input_ext=.fastq

module_type=fastqc

module load $module_type

## PEAR script; Note, no extensions set for output because it auto generates them as .assembled.fastq
echo retrieve files

## $file_list is a file of all sequences to be processed
## If $file_list is uncommented above and set to a value, updates the file_list file which the script uses to get the list of filenames to process

if [ ! -z "$OUT_DIR/file_list" ]; then
   	 ls ../$JOB_NAME/raw_data/*$input_ext > $OUT_DIR/file_list  ## Updates the file_list file
   	 echo Created an updated file_list in $OUT_DIR at `date`.
    else
    sleep 0.1                                       ## Sleep for 0.1 second if not the first task
fi


COUNTER=1

echo `wc -l < $OUT_DIR/file_list`
line_count=`wc -l < $OUT_DIR/file_list`
echo Line_number $line_count

while [ $COUNTER -le $line_count ]
do
echo sample $COUNTER


	LineA=$((2*$COUNTER-1))
	LineB=$((2*$COUNTER))
	A=`echo "$LineA"p`
	B=`echo "$LineB"p`

	echo extracting from line number $A and $B
	R1=`sed -n "$A" $OUT_DIR/file_list`   ## Gets the filename from the $file_list file corresponding to the current $task_id
  	R2=`sed -n "$B" $OUT_DIR/file_list`  

	echo Grabbing filenames from the external file $file_list at `date`.

input_file1=`basename $R1 $input_ext`          ## Removes the directory from the filename to get the basename
input_file2=`basename $R2 $input_ext`
sample_name=${input_file1%%$input_ext}  ## Remove the input_ext from the filename

## Output of variables to help with troubleshooting
echo -e "Variable idx is ${idx}."
echo -e "Variable R1 is $R1."
echo -e "Variable R2 is $R2."
echo -e "Variable sample_name is $sample_name."

echo STEP2: Analyze the quality of the files

mkdir -p $OUT_DIR/out
fastqc -o $OUT_DIR/out -f fastq $R1 $R2

printf "The value of COUNTER=%d\n" $COUNTER
let COUNTER=COUNTER+1
printf "The new value of COUNTER=%d\n" $COUNTER


## Stops the program if the input_file is not found.
if [ ! -s $OUT_DIR/$input_file ]; then
	echo -e "Error: Could not find the \'input_file\'\: $input_file. The \'OUT_DIR\' is\: $OUT_DIR. Check for typos!"
  if [ ! -z "$file_list" ] && [ ! -s $OUT_DIR/$file_list ]; then
    echo The input file list, $file_list, could not be created at $OUT_DIR
  fi
	exit 1
fi

done

echo Finished processing $samplename at `date`.
echo FASTQC $SLURM_JOB_ID >> ../io_log/job_ID.txt

module purge

