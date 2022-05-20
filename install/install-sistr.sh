#!/bin/bash
set -e

pip install numpy

conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install sistr_cmd

pip install pandas==1.0.5
