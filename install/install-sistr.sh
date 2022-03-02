#!/bin/bash
set -e

pip install numpy pandas==1.0.5

conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install sistr_cmd