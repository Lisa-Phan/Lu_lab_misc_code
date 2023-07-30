#! /bin/tcsh -f

#Last mod: 2023-06-22

echo $PWD
#load module gcc 9.4.0
module load gcc/9.4.0

#load module cuda 11.3
module load cuda/11.3

#bash equivalent is export CUDA_HOME=$TACC_CUDA_DIR
setenv CUDA_HOME $TACC_CUDA_DIR
echo $CUDA_HOME

#create a new conda environment in ls6 work directory

#yaml file path
setevn yaml_file '/scratch/09069/dhp563/ESM_fold_setup_test/esm/environment.yml'

#create conda environment
conda env create -n esmfold_20230622 -f $yaml_file
conda init tcsh
conda activate esmfold_20230622

#install pip
conda install -c anaconda pip

#install pytorch
pip install torch==1.12.1+cu113 torchvision==0.13.1+cu113 torchaudio==0.12.1 --extra-index-url 'https://download.pytorch.org/whl/cu113'

#install things listed on github
pip install fair-esm #latest release
pip install 'git+https://github.com/facebookresearch/esm.git' #bleeding edge,current repo main branch

pip install "fair-esm[esmfold]"
# OpenFold and its remaining dependency

pip install 'dllogger @ git+https://github.com/NVIDIA/dllogger.git'
pip install 'openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307'
pip install biotite
pip install rich

