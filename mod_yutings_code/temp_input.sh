#!/usr/bin/env bash

#print usage function
usage() { echo "$0 usage:" && grep " .)\ #" $0; exit 0; }

#print usage if no arguments supplied
[ $# -eq 0 ] && usage

# -h for help 
# -d for directory
# -n for name

#default name
NAME=fastqc_data

while getopts ":hn:d:" arg; do
  case $arg in
    d) # Specify directory of sequencing files.
      DIR=${OPTARG}
      ;;
    n) # Specify optional job_name.
      NAME=${OPTARG}
      ;;
    h | *) # Display help.
      usage
      exit 0
      ;;
  esac
done

#Check if directory is provided
if [ -z "$DIR" ]; then
    echo "Usage: $0 -d <dir> Please provide directory"
    exit 1
fi

#Check if the directory exists
if [ ! -d $DIR ]; then
    echo "Directory does not exist :("
    exit 1
fi

#ititialize directories
mkdir -p $NAME/{step1_fastqc_initial_library,step2,step3}/

ABS_PATH=`readlink -f "$DIR"`
ln -s $ABS_PATH $NAME/raw_data