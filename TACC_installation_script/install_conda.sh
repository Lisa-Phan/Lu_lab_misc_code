#!/bin/bash -e
  
###### INSTALL CONDA IN WORK DIRECTORY ######

cd $WORK
#get conda from source
#latest version
echo 'current director'
echo $(pwd)

wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh

#install conda
bash Miniconda3-py310_23.3.1-0-Linux-x86_64.sh

#when prompted, just ensure that the installation is in work directory