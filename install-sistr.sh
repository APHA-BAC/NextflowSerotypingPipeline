#!/bin/bash

pip3 install --upgrade pip
pip3 install wheel
pip3 install numpy pandas

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels r
conda config --add channels bioconda
conda install sistr_cmd
