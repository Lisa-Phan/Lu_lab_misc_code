#!/bin/bash
#SBATCH -p normal
#SBATCH -J PC
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 48:00:00

echo Step 0: Define variables, tools, and paths at echo `date`.

## Home Directory
project_dir=/scratch/08332/mona163/TH_SELEX
## Input Directory for sequences to be processed
input_dir=/scratch/08332/mona163/WL_LK_SELEX_Raw_data/TH
## Directory for sequences that have been cut
output_dir=$project_dir/PEAR

input_ext=.fastq
fastxmod=fastx_toolkit
all_dir=($output_dir/fastqc $project_dir/merged $project_dir/merged/fastqc)

module load biocontainers

echo Step 1: Check and create directories at `date`.
## Only runs if $all_dir is uncommented above and has something added to it. Otherwise, will create $count_dir and $clust_dir separately later in the script when the "count" or "clust" programs are run.
if [ ! -z $all_dir ]; then
  for d in ${all_dir[@]}; do  ## Loops through all directories listed in array $all_dir
  ## If the directory doesn't exist, creates the directory (and all parent directories).
    if [ ! -d $d ]; then
      mkdir -p $d
      echo Creating directory $d and any parent directories at `date`.
    fi
  done
fi

## PEAR script; Note, no extensions set for output because it auto generates them as .assembled.fastq
echo retrieve files

## If $file_list is uncommented above and set to a value, updates the file_list file which the script uses to get the list of filenames to process

if [ ! -z "$input_dir/file_list" ]; then
   	 ls $input_dir/*$input_ext > $input_dir/file_list  ## Updates the file_list file
   	 echo Created an updated file_list in $input_dir at `date`.
    else
    sleep 0.1                                          ## Sleep for 0.1 second if not the first task
fi

COUNTER=1

echo `wc -l < $input_dir/file_list`
line_count=`wc -l < $input_dir/file_list`
echo Line_number $line_count
max_file_number=$(($line_count / 2))
echo max sample number $max_file_number

while [ $COUNTER -le $max_file_number ]
do
echo sample $COUNTER

	LineA=$((2*$COUNTER-1))
	LineB=$((2*$COUNTER))
	A=`echo "$LineA"p`
	B=`echo "$LineB"p`

	echo extracting from line number $A and $B
	R1=`sed -n "$A" $input_dir/file_list`   ## Gets the filename from the $file_list file corresponding to the current $task_id
  	R2=`sed -n "$B" $input_dir/file_list`  

	echo Grabbing filenames from the external file $file_list at `date`.

input_file1=`basename $R1 $input_ext`          ## Removes the directory from the filename to get the basename
input_file2=`basename $R2 $input_ext`
sample_name=${input_file1%%$input_ext}  ## Remove the input_ext from the filename

## Output of variables to help with troubleshooting
echo -e "Variable R1 is $R1."
echo -e "Variable R2 is $R2."
echo -e "Variable sample_name is $sample_name."

##echo STEP2: Pear/pair the Read1 and Read2 files into a combined 'peared' file

##/work2/08332/mona163/stampede2/pear/bin/pear -f $R1 -r $R2 -o $output_dir/$sample_name

echo STEP3:Duplicate sequence from the second direction
module load $fastxmod

fastx_reverse_complement -i $output_dir/$sample_name.assembled.fastq -o $output_dir/$sample_name.rc.assembled.fastq

echo STEP4:Combine sequences from different directions

cat $output_dir/$sample_name.assembled.fastq $output_dir/$sample_name.rc.assembled.fastq > $project_dir/merged/$sample_name.fastq

echo STEP5:FastQC of the paired files and merged files
module load fastqc
fastqc -t 80 -o $project_dir/merged/fastqc $PEAR_dir/$project_dir/merged/$sample_name.fastq
fastqc -t 80 -o $output_dir/fastqc $output_dir/$sample_name.assembled.fastq

printf "The value of COUNTER=%d\n" $COUNTER
let COUNTER=COUNTER+1
printf "The new value of COUNTER=%d\n" $COUNTER


## Stops the program if the input_file is not found.
if [ ! -s $input_dir/$input_file ]; then
	echo -e "Error: Could not find the \'input_file\'\: $input_file. The \'input_dir\' is\: $input_dir. Check for typos!"
  if [ ! -z "$file_list" ] && [ ! -s $input_dir/$file_list ]; then
    echo The input file list, $file_list, could not be created at $input_dir
  fi
	exit 1
fi


done

echo Finished processing $samplename at `date`.

module purge
