#!/bin/bash
#SBATCH -p normal
#SBATCH -J Cutadapt
#SBATCH -D /home1/08332/mona163/src
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -t 20:00:00



echo Step 0: Define variables, tools, primer_sequences and paths at echo `date`.

## Home Directory
project_dir=/scratch/08332/mona163/TH_SELEX
## Input Directory for sequences to be processed
input_dir=$project_dir/merged
## Directory for sequences that have been cut
output_dir=$project_dir/Cutadapt_TH
input_ext=_R1.fastq
filearray=filenames.txt

trimmod=cutadapt
fastqcmod=fastqc
fastxmod=fastx_toolkit

all_dir=($output_dir $output_dir/fastqc)

## Forward primer of selection (or 5' adapter of sequencing)
fwdprimer=GCTGACTAGTACATGACCAGG
## Reverse primer of selection (or 3' adapter of sequencing) in "forward" direction, 5'-3' as you would order from IDT
revprimer=CTGCGATACGTCTACTATGGC

## Maximum error rate for adapters (e.g. 0.1 for fivepadapt of 10nts = 1nt mismatch allowed, it should be percentage. Details could be found in the website of the pipline)
error=0.1
## Minimum length of sequences (sequences below this # of nt are discarded)
minlength=84
## Maximum length of sequences (sequences above this # of nt are discarded, for the synthetic pools, it will be the size of the random regions. Such as 48 min and 52 max for the N50 pool)
maxlength=86

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

echo STEP2: load module, define input files.

module load biocontainers
module load $trimmod

fivepadapt1=$fwdprimer
fivepadapt2=$revprimer
threepadapt1=$( echo $fivepadapt2 | tr ACGTacgt TGCAtgca | rev )

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
max_file_number=$line_count
echo max sample number $max_file_number

while [ $COUNTER -le $max_file_number ]
do
echo sample $COUNTER

	A=`echo "$COUNTER"p`

	echo extracting from line number $A
	R=`sed -n "$A" $input_dir/file_list`   ## Gets the filename from the $file_list file corresponding to the current $task_id

	echo Grabbing filenames from the external file $file_list at `date`.

input_file=`basename $R $input_ext`          ## Removes the directory from the filename to get the basename
sample_name=${input_file%%$input_ext}        ## Remove the input_ext from the filename

## Output of variables to help with troubleshooting
echo -e "Variable R is $R"
echo -e "Variable sample_name is $sample_name."
echo -e "Variable input_file is $input_file."

echo STEP2: Remove Illumina universal adapter sequences at 'date'.

cutadapt -j 80 -a AGATCGGAAGAG -e $error --discard-trimmed -o $output_dir/$sample_name.T1.fastq $R
cutadapt -j 80 -b file:/home1/08332/mona163/src/LK_CoA_SELEX/shorter_complex.fasta -n 4 --discard-trimmed -e 0.2 -o $output_dir/$sample_name.T2.fastq $output_dir/$sample_name.T1.fastq
echo STEP3: Remove adapter sequences and remain sequence with correct size at 'date'.
cutadapt -j 80 -a $fivepadapt1...$threepadapt1 -e $error -m $minlength -M $maxlength --discard-untrimmed -o $output_dir/$sample_name.Trim.fastq $output_dir/$sample_name.T2.fastq
exitcode=$?
 if [ $exitcode -ne 0 ]
 then
    echo The cutadapt program stopped at sample=${samplename}, with Error Code=$exitcode
    exit $exitcode
 fi
echo Finished trimming $samplename in $outputdir at `date`

echo STEP4: FastQC of the merged files to assess quality of selection pool sequences at `date`.
module load $fastqcmod

fastqc -t 80 -o $output_dir/fastqc $output_dir/$sample_name.Trim.fastq

echo Finished processing $R at `date`.

printf "The value of COUNTER=%d\n" $COUNTER
let COUNTER=COUNTER+1
printf "The new value of COUNTER=%d\n" $COUNTER

done

echo Finished processing at `date`.

module purge
