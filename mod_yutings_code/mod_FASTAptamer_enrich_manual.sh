#!/bin/bash
## FASTAptamer Enrich Script created by Ryan Lake, Version 2020-02-26.
## mod by Lisa Phan 2023-05-22

## See bottom of file for Script Description.
## ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 00:30:00
#SBATCH -J fastapt
#SBATCH -o %x.out
#SBATCH -e %x.err
## The below starts this job only after a previous job finishes (enter prev job id & uncomment below)

#get var from param script function
get_var() { grep "$1" params.txt | cut -d '=' -f2 ; }

## -------------Script Options-----------------
## Specify if you want to do enrichment of "count" or "clust" files.
enrich_input=clust

## Specify if you want to automatically generate all possible combinations of the input files 
## specified with "auto" or if you want to manually choose which combinations to create with "manual".
script_mode=manual

echo -e "Starting FASTAptamer Enrich in \'$script_mode mode\' with \'$enrich_input inputs\' at `date`."

## -------------Define Variables-----------------
echo -e "\nPrep Step: Define paths, tools, and variables at `date`. \n"

## Generates date in YYYYMMDD format. Useful to generate file and folder names for outputs.
day=$( date +%Y%m%d )

## Your Home Directory
JOB_NAME=`tail -n 1 .userlog.txt | cut -d" " -f4`

## Input Directory for sequences to be processed; 
## Automatically adjusts for "count" or "clust" based on $enrich_input, or you can change to set sub-directory manually yourself.
input_dir=/scratch/08332/mona163/TH_SELEX/Counster_TH/clust-8001


## Directory for results of FASTAptamer enrich
output_dir=/scratch/08332/mona163/TH_SELEX/Counster_TH/enrich-clust

## Directory where the perl script files of FASTAptamer are stored
fastapt_dir=/work2/08332/mona163/stampede2/FASTAptamer

## Type all output directories that you want to have created into this array variable, all_dir, separated by spaces.
## Note: It will automatically create all parent directories, so inputting "/parent/child" creates both /parent and /parent/child, therefore you don't need to include any parent directories if its child is already listed.
all_dir=($output_dir)

## Filename extensions for input files; Output file extensions are defined within below script.
input_ext=.$enrich_input.fa

## For Manual Mode Only: Set the three filenames you want to compare here (don't put the input_ext at the end).
## Note: If using Automatic Mode, it's fine to leave this as is because they will get overwritten.
rX=TH-18_S52-10
rY=TH-21_S55-6
rZ=TH-24_S58-4

## Module defined for Perl to run FASTAptamer Perl scripts
perl_mod=biocontainers

if [ ! -z $all_dir ]; then
  for d in ${all_dir[@]}; do  ## Loops through all directories listed in array $all_dir
  ## If the directory doesn't exist, creates the directory (and all parent directories).
    if [ ! -d $d ]; then
      mkdir -p $d
      echo Creating directory $d and any parent directories at `date`.
    fi
  done
fi

output_file=$rX-$rY-$rZ.$enrich_input.tsv

module load $perl_mod

perl $fastapt_dir/fastaptamer_enrich -x $input_dir/$rX$input_ext -y $input_dir/$rY$input_ext -z $input_dir/$rZ$input_ext -o $output_dir/$output_file
