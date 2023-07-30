#!/bin/bash

#SBATCH -p normal
#SBATCH -J assembly
#SBATCH -e %x.err
#SBATCH -o %x.out
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 04:00:00

##Comment block
## Last mod 5/30/23
##Dependent modules: fastx_toolkit, biocontainers, pear 
##Note for future iterations, all paths should be supplied as a command line argument 
## Error ###: R1 option not working
## Error ##: redundant code to handle R1 and R2 options

echo Step 0: Define variables, tools, and paths at echo `date`.

input_ext=.fastq

#input sequence is the same as prior step
JOB_NAME=`tail -n 1 .userlog.txt | cut -d" " -f4`
REP_OPT=`tail -n 1 .userlog.txt | cut -d" " -f6`
OUT_DIR=$JOB_NAME/step2_assembly
INPUT=$JOB_NAME/raw_data/file_list

module load biocontainers

echo variables $JOB_NAME, $REP_OPT, $INPUT, $OUT_DIR

## Unclear if this is a good step?
## Ensure that directory structure is good should be enforced at the init.sh step
## all_dir=(../$JOB_NAME/step1_fastqc $/merged $project_dir/merged/fastqc)


# echo Step 1: Check and create directories at `date`.
# ## Only runs if $all_dir is uncommented above and has something added to it. Otherwise, will create $count_dir and $clust_dir separately later in the script when the "count" or "clust" programs are run.
# if [ ! -z $all_dir ]; then
#   for d in ${all_dir[@]}; do  ## Loops through all directories listed in array $all_dir
#   ## If the directory doesn't exist, creates the directory (and all parent directories).
#     if [ ! -d $d ]; then
#       mkdir -p $d
#       echo Creating directory $d and any parent directories at `date`.
#     fi
#   done
# fi

## PEAR script; Note, no extensions set for output because it auto generates them as .assembled.fastq
echo retrieve files

## If $file_list is uncommented above and set to a value, updates the file_list file which the script uses to get the list of filenames to process

# if [ ! -z "$input_dir/file_list" ]; then
#    	 ls $input_dir/*$input_ext > $input_dir/file_list  ## Updates the file_list file
#    	 echo Created an updated file_list in $input_dir at `date`.
#     else
#     sleep 0.1                                          ## Sleep for 0.1 second if not the first task
# fi

COUNTER=1

#Iterate over file list
#if repeat option is 1, then max file number is the same as file number
echo `wc -l < $INPUT`
line_count=`wc -l < $INPUT`
echo Line_number $line_count

## If there are two directions to read per seq,
## sample count is half the listed count 
if [ "$REPT_OPT" == 2 ]; then 
  max_file_number=$(($line_count / 2))
  while [ $COUNTER -le $max_file_number ]
  do
    echo sample $COUNTER

    LineA=$((2*$COUNTER-1))
    LineB=$((2*$COUNTER))
    A=`echo "$LineA"p`
    B=`echo "$LineB"p`

    echo extracting from line number $A and $B
    R1=`sed -n "$A" $INPUT`   ## Gets the filename from the $file_list file corresponding to the current $task_id
    R2=`sed -n "$B" $INPUT`  

    echo Grabbing filenames from the external file $file_list at `date`.

    input_file1=`basename $R1 $input_ext`          ## Removes the directory from the filename to get the basename
    input_file2=`basename $R2 $input_ext`
    sample_name=${input_file1%%$input_ext}         ## Remove the input_ext from the filename

    ## Output of variables to help with troubleshooting
    echo -e "Variable R1 is $R1."
    echo -e "Variable R2 is $R2."
    echo -e "Variable sample_name is $sample_name."

    echo STEP2: Pear/pair the Read1 and Read2 files into a combined 'peared' file
    dependencies/pear/pear-0.9.11-linux-x86_64/bin/pear -f $R1 -r $R2 -o $OUT_DIR/$sample_name

    echo STEP3:Duplicate sequence from the second direction
    module load fastx_toolkit
    fastx_reverse_complement -i $OUT_DIR/$sample_name.assembled.fastq -o $OUT_DIR/$sample_name.rc.assembled.fastq
    
    echo STEP4:Combine sequences from different directions
    mkdir -p $OUT_DIR/merged
    cat $OUT_DIR/$sample_name.assembled.fastq $OUT_DIR/$sample_name.rc.assembled.fastq > $OUT_DIR/merged/$sample_name.merged.fastq

    echo STEP5:FastQC of the paired files and merged files
    module load fastqc

    #extra_fastqc step
    EXTRAQC=$OUT_DIR/fastqc_extra
    mkdir -p $EXTRAQC/{merged,paired}/

    fastqc -t 80 -o $EXTRAQC/merged $OUT_DIR/merged/$sample_name.merged.fastq
    fastqc -t 80 -o $EXTRAQC/paired $OUT_DIR/$sample_name.assembled.fastq

    printf "The value of COUNTER=%d\n" $COUNTER
    let COUNTER=COUNTER+1
    printf "The new value of COUNTER=%d\n" $COUNTER
  done

else 
  #If there are only one read per lane, then max file number is the same as line count
  max_file_number=$line_count
  while [ $COUNTER -le $max_file_number ]; do
    echo sample $COUNTER
    LineA=$((1*$COUNTER))
    A=`echo "$LineA"p`

    echo extracting from line number $A
    R1=`sed -n "$A" $INPUT`   ## Gets the filename from the $file_list file corresponding to the current $task_id

    echo Grabbing filenames from the external file $file_list at `date`.
    input_file1=`basename $R1 $input_ext`          ## Removes the directory from the filename to get the basename
    sample_name=${input_file1%%$input_ext}  ## Remove the input_ext from the filename
    echo -e "Variable R1 is $R1."
    echo -e "Variable sample_name is $sample_name."
    
    echo STEP2: generate reverse compliment
    module load fastx_toolkit
    fastx_reverse_complement -i $R1 -o $OUT_DIR/$sample_name.rc.fastq

    echo STEP4:Combine sequences from different directions
    mkdir -p $OUT_DIR/merged
    cat $R1 $OUT_DIR/$sample_name.rc.fastq > $OUT_DIR/merged/$sample_name.merged.fastq
    
    echo STEP5:FastQC of the paired files and merged files
    module load fastqc
    
    #extra_fastqc step
    EXTRAQC=$OUT_DIR/fastqc_extra
    mkdir -p $EXTRAQC/{merged}/

    fastqc -t 80 -o $EXTRAQC/merged $OUT_DIR/merged/$sample_name.merged.fastq

    printf "The value of COUNTER=%d\n" $COUNTER
    let COUNTER=COUNTER+1  
    printf "The new value of COUNTER=%d\n" $COUNTER
    echo max sample number $max_file_number
  done
fi

#Underlying assumption: file is already ordered

# Stops the program if the input_file is not found
if [ ! -s $input_dir/$input_file ]; then
	echo -e "Error: Could not find the \'input_file\'\: $input_file. The \'input_dir\' is\: $input_dir. Check for typos!"
  if [ ! -z "$file_list" ] && [ ! -s $input_dir/$file_list ]; then
    echo The input file list, $file_list, could not be created at $input_dir
  fi
	exit 1
fi

echo Finished processing $samplename at `date`.
module purge