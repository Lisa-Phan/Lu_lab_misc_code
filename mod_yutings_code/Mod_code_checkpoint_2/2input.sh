#!/bin/bash
#usage function
usage() { echo 'Usage: input.sh -d|--dir directory for sequencing file';
          echo 'Optional arguments: -n|--name for job name, default to job_name';
          echo '                    -h|--help for help';
          echo '                    -r|--rept_opt for paired ends read option, default to 2';
          exit 0; }

if [[ $# -eq 0 ]] ; then
  usage
  exit 1
fi

#default variables
NAME='job_name'
REP=2

while [[ $# -gt 0 ]]; do
  opt="$1"
  shift;
  current_arg="$1"
  if [[ "$current_arg" =~ ^-{1,2}.* ]]; then
    echo "WARNING: You may have left an argument blank. Double check your command." 
  fi
  case "$opt" in
    "-d"|"--dir"      ) DIR="$1"; shift;;
    "-n"|"--name"     ) NAME="$1"; shift;;
    "-r"|"--rep_opt"  ) REP="$1"; shift;;
    "-h"|"--help"     ) usage; shift;;
    *                 ) echo "ERROR: Invalid option: \""$opt"\"" >&2
                      exit 1;;
  esac
done

if [[ "$DIR" == "" ]]; then
  echo "ERROR: Options -d is a required argument." >&2
  exit 1
fi


#Check if directory is provided
if [ -z "$DIR" ]; then
usage    
exit 1
fi

#Check if the directory exists
if [ ! -d $DIR ]; then
    echo "Directory does not exist :("
    exit 1
fi

#ititialize directories
mkdir -p $NAME/{step1_fastqc,step2_assembly,step3_cutadapt}/

ABS_PATH=`readlink -f "$DIR"`
ln -s $ABS_PATH $NAME/raw_data

#record variable for later use
echo SYM_LINK $ABS_PATH JOB_NAME $NAME REPEAT_OPTION $REPEAT `date +"%Y_%m_%d_%M:%S"` >> .userlog.txt

#run the flow control script
scripts/control.sh
