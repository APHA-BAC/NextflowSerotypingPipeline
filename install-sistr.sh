#!/bin/bash

pip install --upgrade pip
pip install wheel
pip install numpy pandas==1.1.5

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install sistr_cmd
