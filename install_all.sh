#!/bin/sh
##################################################
# Setup recipe for RaSE and all its dependencies
# This script assumes current directory as the target installation path 
# Please run this script from the desired installation path
##################################################

# Download miniconda installer 
wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh

#Run miniconda installer
bash miniconda.sh -b -p ./conda
MYCONDAPATH=$(readlink -e ./conda/bin/)
export PATH=$MYCONDAPATH:$PATH

#Update conda

conda update conda
export PYTHONNOUSERSITE=1 

#Activate virtual env for rase
source activate condarase
#Install Rase Dependencies
conda install dill scikit-learn networkx=1.10 docopt toolz joblib pydot graphviz requests matplotlib


# Clone Rase and EdEN repositories, todo: Specify release tag
git clone https://github.com/fabriziocosta/RaSE.git
git clone https://github.com/fabriziocosta/EDeN.git

MYEDENPATH=$(readlink -e ./EDeN/)
MYRASEPATH=$(readlink -e ./RaSE/)
export PYTHONPATH=$MYEDENPATH:$PYTHONPATH

# This is needed for EDeN invoking the vienna package tools
export PATH=/scratch/rna/bisge001/Software/ViennaRNA/2.2.10/bin/:$PATH 
python2 $MYRASEPATH/code/RaSE.py -i 'ACGUGGCUG'
